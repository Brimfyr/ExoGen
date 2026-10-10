"""Phase 0 spike: energy-balance climate model with Held-Hou cell width and Williams-Kasting transport.

Checks the relations in PLAN.md section 4.6 and docs/parameters.md (F21, F22, F32, F33).
The model is one-dimensional in the angle that the forcing depends on: colatitude for a
rotating body with annual-mean insolation, and the angle from the substellar point for a
tidally locked body. Both are exact reductions of the 2D sphere model when the surface is
uniform, which is all this spike needs. The 2D version lives on the cubesphere in Rust.

Run:  research/.venv/Scripts/python research/ebm_climate.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.linalg import solve_banded  # noqa: E402

OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(exist_ok=True)

SIGMA = 5.670374e-8
OMEGA_EARTH = 7.2921e-5
D_EARTH = 0.58  # W m^-2 K^-1, Williams and Kasting 1997 / North 1975 family


# ---------------------------------------------------------------------------
# Insolation
# ---------------------------------------------------------------------------

def annual_mean_insolation(lat_rad, s0, obliquity_deg, ecc=0.0, n_orbit=720):
    """Annual-mean daily-mean insolation [W/m^2] at each latitude, integrated over the orbit.

    Perihelion is placed at the northern spring equinox; for annual means the choice only
    matters at second order in eccentricity.
    """
    eps = math.radians(obliquity_deg)
    nu = np.linspace(0.0, 2.0 * math.pi, n_orbit, endpoint=False)
    r_over_a = (1.0 - ecc**2) / (1.0 + ecc * np.cos(nu))
    # time weighting: dt proportional to r^2 dnu (Kepler's second law)
    w = r_over_a**2
    w /= w.sum()
    dec = np.arcsin(np.sin(eps) * np.sin(nu))  # solar declination along the orbit
    lat = lat_rad[:, None]
    cos_h0 = np.clip(-np.tan(lat) * np.tan(dec)[None, :], -1.0, 1.0)
    h0 = np.arccos(cos_h0)
    daily = (s0 / math.pi) / r_over_a[None, :] ** 2 * (
        h0 * np.sin(lat) * np.sin(dec)[None, :] + np.cos(lat) * np.cos(dec)[None, :] * np.sin(h0)
    )
    return (daily * w[None, :]).sum(axis=1)


def locked_insolation(theta_rad, s0):
    """Permanent insolation [W/m^2] as a function of the angle from the substellar point."""
    return s0 * np.maximum(0.0, np.cos(theta_rad))


# ---------------------------------------------------------------------------
# Greenhouse and transport parametrisations
# ---------------------------------------------------------------------------

def grey_tau_for_surface_temperature(t_eq, t_s):
    """Grey optical depth that lifts T_eq to T_s in the zero-dimensional estimate T_s^4 = T_eq^4 (1 + 3 tau / 4)."""
    return (4.0 / 3.0) * ((t_s / t_eq) ** 4 - 1.0)


def olr_grey(t, tau):
    return SIGMA * t**4 / (1.0 + 0.75 * tau)


def transport_coefficient(p_bar, omega, cp_ratio=1.0, mu_ratio=1.0, cap=100.0):
    """Williams and Kasting 1997 scaling as used in exoplanet energy-balance models (F33)."""
    d = D_EARTH * p_bar * cp_ratio * mu_ratio * (OMEGA_EARTH / max(omega, 1e-12)) ** 2
    return min(d, cap * D_EARTH)


def held_hou_edge_deg(g, h_trop, delta_h, omega, radius):
    """Held and Hou 1980 small-angle estimate of the Hadley cell edge latitude (F32)."""
    phi2 = 5.0 * g * h_trop * delta_h / (3.0 * omega**2 * radius**2)
    return math.degrees(min(math.pi / 2.0, math.sqrt(max(phi2, 0.0))))


# ---------------------------------------------------------------------------
# One-dimensional energy balance on the sphere in a single polar angle theta in [0, pi]
#   0 = S(theta) (1 - alpha(T)) - OLR(T) + D (1/sin theta) d/dtheta (sin theta dT/dtheta)
# Solved by pseudo-time stepping with implicit diffusion and explicit forcing.
# ---------------------------------------------------------------------------

@dataclass
class EBMResult:
    theta: np.ndarray
    temperature: np.ndarray
    insolation: np.ndarray
    albedo: np.ndarray


def solve_ebm(theta, insolation, tau, d_coeff, albedo_warm, albedo_cold=0.6, t_ice=263.0,
              t_floor=30.0, ice_feedback=True, n_iter=40000, dt=0.05, heat_capacity=1.0,
              t_init=310.0):
    """Relax to equilibrium from a warm uniform start (t_init), which selects the warm branch
    where the ice-albedo feedback admits two solutions. A cold start lands on the snowball branch."""
    n = len(theta)
    dth = theta[1] - theta[0]
    sin_c = np.sin(theta)
    # face values of sin(theta) for the conservative Laplacian
    sin_f = np.sin(np.concatenate(([theta[0] - 0.5 * dth], 0.5 * (theta[1:] + theta[:-1]), [theta[-1] + 0.5 * dth])))
    sin_f[0] = 0.0
    sin_f[-1] = 0.0
    sin_c = np.maximum(sin_c, 1e-6)
    coeff = d_coeff * dt / heat_capacity / (dth * dth)
    lower = -coeff * sin_f[:-1] / sin_c
    upper = -coeff * sin_f[1:] / sin_c
    diag = 1.0 + coeff * (sin_f[:-1] + sin_f[1:]) / sin_c
    ab = np.zeros((3, n))
    ab[0, 1:] = upper[:-1]
    ab[1, :] = diag
    ab[2, :-1] = lower[1:]

    t = np.full(n, float(t_init))
    for _ in range(n_iter):
        alpha = np.where(t < t_ice, albedo_cold, albedo_warm) if ice_feedback else np.full(n, albedo_warm)
        forcing = insolation * (1.0 - alpha) - olr_grey(t, tau)
        rhs = t + dt / heat_capacity * forcing
        t_new = solve_banded((1, 1), ab, rhs)
        t_new = np.maximum(t_new, t_floor)
        if np.max(np.abs(t_new - t)) < 1e-6:
            t = t_new
            break
        t = t_new
    alpha = np.where(t < t_ice, albedo_cold, albedo_warm) if ice_feedback else np.full(n, albedo_warm)
    return EBMResult(theta, t, insolation, alpha)


def area_mean(theta, values):
    w = np.sin(theta)
    return float((values * w).sum() / w.sum())


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------

@dataclass
class Case:
    name: str
    s0: float            # W/m^2 at the body
    albedo: float        # planetary albedo of ice-free surface plus clouds
    t_surface_obs: float  # observed global mean surface temperature K
    p_bar: float
    omega: float
    obliquity_deg: float
    ecc: float
    g: float
    radius: float
    h_trop: float        # tropopause height m
    delta_h: float       # fractional equator-pole radiative contrast for Held-Hou
    locked: bool = False
    ice_feedback: bool = True
    cp_ratio: float = 1.0
    mu_ratio: float = 1.0


CASES = [
    Case("Earth", 1361.0, 0.30, 288.0, 1.013, OMEGA_EARTH, 23.44, 0.0167, 9.81, 6.371e6, 1.5e4, 1.0 / 4.0),
    Case("Mars", 589.0, 0.25, 213.0, 0.006, 7.088e-5, 25.2, 0.0934, 3.71, 3.39e6, 4.0e4, 1.0 / 2.0, ice_feedback=False),
    Case("Venus", 2601.0, 0.77, 737.0, 92.0, 2.99e-7, 2.6, 0.0068, 8.87, 6.052e6, 6.0e4, 1.0 / 10.0, ice_feedback=False, cp_ratio=0.85, mu_ratio=0.66),
    Case("Titan", 15.1, 0.22, 94.0, 1.47, 4.56e-6, 26.7, 0.0288, 1.35, 2.575e6, 4.0e4, 1.0 / 5.0, ice_feedback=False),
    Case("Locked Earth-twin, 10 d orbit, 1 bar", 1361.0, 0.30, float("nan"), 1.013, 2.0 * math.pi / (10 * 86400), 0.0, 0.0, 9.81, 6.371e6, 1.5e4, 1.0, locked=True),
    Case("Locked Earth-twin, airless", 1361.0, 0.12, float("nan"), 0.0, 2.0 * math.pi / (10 * 86400), 0.0, 0.0, 9.81, 6.371e6, 1.5e4, 1.0, locked=True, ice_feedback=False),
]


def main() -> None:
    lines: list[str] = []

    def say(s: str = "") -> None:
        print(s)
        lines.append(s)

    say("ExoGen Phase 0 spike: energy-balance climate")
    say("=" * 60)

    theta = np.linspace(0.0, math.pi, 181)  # 1 degree
    theta = 0.5 * (theta[1:] + theta[:-1])  # cell centres
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    axes = axes.ravel()

    for i, c in enumerate(CASES):
        t_eq = (c.s0 * (1.0 - c.albedo) / (4.0 * SIGMA)) ** 0.25
        say(f"\n{c.name}")
        if c.locked:
            insol = locked_insolation(theta, c.s0)
            tau = grey_tau_for_surface_temperature(t_eq, 288.0) if c.p_bar > 0 else 0.0
            # The Williams-Kasting scaling grows without bound as rotation slows; show the
            # sensitivity to the cap, which Phase 5 calibrates against published GCM contrasts.
            caps = (100.0, 10.0) if c.p_bar > 0 else (0.0,)
            for cap in caps:
                d_coeff = transport_coefficient(c.p_bar, c.omega, cap=cap) if c.p_bar > 0 else 0.0
                res = solve_ebm(theta, insol, tau, d_coeff, c.albedo, ice_feedback=c.ice_feedback, t_ice=263.0)
                t_mean = area_mean(theta, res.temperature)
                t_ss = res.temperature[0]
                t_as = res.temperature[-1]
                t_term = float(np.interp(math.pi / 2.0, theta, res.temperature))
                say(f"  cap {cap:5.0f} x D_E: tau = {tau:6.3f}  D = {d_coeff:7.3f}  global mean {t_mean:6.1f} K  "
                    f"substellar {t_ss:6.1f}  terminator {t_term:6.1f}  antistellar {t_as:6.1f} K")
                frozen = theta[res.temperature < 273.15]
                if len(frozen):
                    say(f"    water freezes beyond {math.degrees(frozen.min()):5.1f} deg from the substellar point")
                else:
                    say("    no freezing anywhere")
        else:
            lat = math.pi / 2.0 - theta
            insol = annual_mean_insolation(lat, c.s0, c.obliquity_deg, c.ecc)
            tau = grey_tau_for_surface_temperature(t_eq, c.t_surface_obs)
            d_coeff = transport_coefficient(c.p_bar, c.omega, c.cp_ratio, c.mu_ratio)
            res = solve_ebm(theta, insol, tau, d_coeff, c.albedo, ice_feedback=c.ice_feedback,
                            t_ice=263.0 if c.name.startswith("Earth") else 150.0)
            t_mean = area_mean(theta, res.temperature)
            say(f"  T_eq = {t_eq:6.1f} K   grey tau = {tau:8.3f}   D = {d_coeff:7.3f} W/m^2/K")
            say(f"  global mean T = {t_mean:6.1f} K  (observed {c.t_surface_obs:.0f})")
            lat_deg = np.degrees(math.pi / 2.0 - theta)
            t_at = lambda phi: float(np.interp(-phi, -lat_deg, res.temperature))  # noqa: E731
            say(f"  equator {t_at(0):6.1f} K   30N {t_at(30):6.1f} K   60N {t_at(60):6.1f} K   pole {t_at(89.5):6.1f} K")
            ice = lat_deg[res.temperature < 263.0]
            if len(ice) and c.ice_feedback:
                say(f"  ice line at {np.abs(ice).min():5.1f} deg latitude")
            edge = held_hou_edge_deg(c.g, c.h_trop, c.delta_h, c.omega, c.radius)
            say(f"  Held-Hou Hadley edge = {edge:5.1f} deg")

        ax = axes[i]
        x = np.degrees(theta) if c.locked else np.degrees(math.pi / 2.0 - theta)
        ax.plot(x, res.temperature, "k-")
        ax.set_title(c.name, fontsize=9)
        ax.set_xlabel("angle from substellar point [deg]" if c.locked else "latitude [deg]")
        ax.set_ylabel("T [K]")
        ax.grid(True, alpha=0.3)

    say("\nObserved for comparison: Earth zonal annual mean about 299 K equator, 293 K 30N, 273 K 60N, 255 K north pole;")
    say("Mars about 220 K equator and 150-160 K poles; Venus uniform 737 K; Titan 94 K with a 3 K equator-pole contrast.")
    say("Locked 1 bar Earth-twin GCMs (Joshi et al. 1997 and later) give a nightside 30-60 K colder than the dayside, not frozen out.")

    fig.tight_layout()
    fig.savefig(OUT / "ebm_profiles.png", dpi=130)
    (OUT / "ebm_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    say(f"\nPlot and report written to {OUT}")


if __name__ == "__main__":
    main()
