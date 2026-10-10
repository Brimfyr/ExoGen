"""Phase 0 spike: crater production, scaling and morphology laws.

Checks the relations in PLAN.md section 4.5 and docs/parameters.md (F24, F30, F31)
against published anchor values, prints a report and writes plots to research/output/.

Run:  research/.venv/Scripts/python research/crater_scaling.py
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

OUT = Path(__file__).resolve().parent / "output"
OUT.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Production function and chronology (Neukum, Ivanov and Hartmann 2001)
# ---------------------------------------------------------------------------

# Cumulative N(>=D) per km^2 for a 1 Gyr surface, log10 polynomial in log10(D [km]),
# valid 0.01 km <= D <= 300 km.
NPF_A = [
    -3.0876, -3.557528, 0.781027, 1.021521, -0.156012, -0.444058,
    0.019977, 0.086850, -0.005874, -0.006809, 8.25e-4, 5.54e-5,
]


def npf_cumulative_1gyr(d_km):
    """Cumulative crater density N(>=D) [per km^2] for a 1 Gyr lunar surface."""
    logd = np.log10(np.asarray(d_km, dtype=float))
    acc = np.zeros_like(logd)
    for j, a in enumerate(NPF_A):
        acc = acc + a * logd**j
    return 10.0**acc


def chronology_n1(t_gyr):
    """N(>=1 km) per km^2 as a function of surface age in Gyr (lunar chronology function)."""
    t = np.asarray(t_gyr, dtype=float)
    return 5.44e-14 * (np.exp(6.93 * t) - 1.0) + 8.38e-4 * t


def isochron(d_km, t_gyr):
    """Cumulative density at diameter D for a surface of age T."""
    return npf_cumulative_1gyr(d_km) * chronology_n1(t_gyr) / chronology_n1(1.0)


# ---------------------------------------------------------------------------
# Impactor to crater: pi-group scaling in the form of Collins, Melosh and Marcus 2005
# ---------------------------------------------------------------------------

def transient_diameter(L, v, g, rho_i=2700.0, rho_t=2700.0, theta_deg=45.0):
    """Transient crater diameter [m] from impactor diameter L [m], speed v [m/s], gravity g [m/s^2]."""
    return (1.161 * (rho_i / rho_t) ** (1.0 / 3.0) * L**0.78 * v**0.44 * g**-0.22
            * math.sin(math.radians(theta_deg)) ** (1.0 / 3.0))


def final_diameter(d_tc, d_sc):
    """Final rim-to-rim diameter [m] from the transient diameter and the simple-complex transition d_sc [m]."""
    d_simple = 1.25 * d_tc
    if d_simple < d_sc:
        return d_simple
    # Collins et al. 2005 eq. 27 (Croft collapse), quantities in km inside the power law.
    return 1.17e3 * (d_tc / 1e3) ** 1.13 / (d_sc / 1e3) ** 0.13


def impactor_for_crater(d_final, v, g, d_sc, **kw):
    """Invert final_diameter(transient_diameter(L)) for L by bisection."""
    lo, hi = 0.01, 1.0e6
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if final_diameter(transient_diameter(mid, v, g, **kw), d_sc) < d_final:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


# ---------------------------------------------------------------------------
# Morphometry (Pike 1977, 1980; lunar). D in km, results in km.
# ---------------------------------------------------------------------------

def depth_simple(d_km):
    return 0.196 * d_km**1.010


def depth_complex(d_km):
    return 1.044 * d_km**0.301


def rim_height_simple(d_km):
    return 0.036 * d_km**1.014


def rim_height_complex(d_km):
    return 0.236 * d_km**0.399


def central_peak_height(d_km):
    return 0.032 * d_km**0.900


def depth_mars_garvin(d_km):
    """Garvin et al. 2003 Mars fits (verify): simple d = 0.21 D^0.81, complex d = 0.36 D^0.49."""
    d = np.asarray(d_km, dtype=float)
    return np.where(d < 7.0, 0.21 * d**0.81, 0.36 * d**0.49)


def ejecta_thickness(r_m, R_m):
    """McGetchin et al. 1973: thickness [m] at range r from centre for a crater of rim radius R (both m), r >= R."""
    return 0.14 * R_m**0.74 * (r_m / R_m) ** -3.0


# ---------------------------------------------------------------------------
# Scaling with gravity: simple-to-complex transition (F30) and atmospheric cut-off (F24)
# ---------------------------------------------------------------------------

def d_sc_predicted_km(g, icy=False):
    return 16.5 * (1.62 / g) * (0.15 if icy else 1.0)


def crater_cutoff_km(p_pa, g, c=None):
    col = p_pa / g
    if c is None:
        c = 0.15 / (101325.0 / 9.81) ** 0.6  # Earth anchored at 150 m
    return c * col**0.6


# ---------------------------------------------------------------------------
# Diffusional degradation (Fassett and Thomson 2014): linear diffusion of a crater profile
# ---------------------------------------------------------------------------

def fresh_profile(r, D):
    """Axisymmetric fresh simple crater: parabolic bowl of depth 0.2 D, rim 0.036 D, ejecta decay outside."""
    R = D / 2.0
    d = 0.2 * D
    h_rim = 0.036 * D
    h = np.empty_like(r)
    inside = r <= R
    h[inside] = -d + (d + h_rim) * (r[inside] / R) ** 2
    out = ~inside
    h[out] = h_rim * (r[out] / R) ** -3.0
    return h


def diffuse_radial(h, r, kt, steps=4000):
    """Explicit radial diffusion of profile h(r) for a total kappa*t = kt [m^2]."""
    dr = r[1] - r[0]
    dt = 0.2 * dr * dr
    n = max(1, int(math.ceil(kt / dt)))
    dt = kt / n
    h = h.copy()
    for _ in range(n):
        dh = np.zeros_like(h)
        # (1/r) d/dr (r dh/dr) with a symmetric centre
        flux = np.zeros(len(h) + 1)
        rf = np.concatenate(([0.0], 0.5 * (r[1:] + r[:-1]), [r[-1] + 0.5 * dr]))
        flux[1:-1] = rf[1:-1] * (h[1:] - h[:-1]) / dr
        dh[:] = (flux[1:] - flux[:-1]) / (np.maximum(r, 0.5 * dr) * dr)
        dh[0] = 4.0 * (h[1] - h[0]) / (dr * dr)  # axis
        dh[-1] = 0.0
        h += dt * dh
    return h


# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------

BODIES = {
    # name: (g m/s^2, icy, observed D_sc km, impact velocity km/s, surface pressure Pa)
    "Moon": (1.62, False, 16.5, 17.0, 0.0),
    "Mercury": (3.70, False, 10.3, 40.0, 0.0),
    "Mars": (3.71, False, 6.5, 10.0, 610.0),
    "Earth": (9.81, False, 3.2, 20.0, 101325.0),
    "Venus": (8.87, False, 3.0, 25.0, 9.2e6),
    "Ganymede": (1.43, True, 2.5, 20.0, 0.0),
    "Callisto": (1.24, True, 2.5, 15.0, 0.0),
    "Titan": (1.35, True, 2.5, 10.5, 1.47e5),
}

# Lunar craters with measured depth and diameter (km) for the morphometry check.
LUNAR_CRATERS = {
    "Linne": (2.2, 0.5), "Bessel": (16.0, 1.7), "Aristarchus": (40.0, 3.7),
    "Tycho": (85.0, 4.8), "Copernicus": (93.0, 3.8), "Theophilus": (99.0, 4.2),
    "Langrenus": (132.0, 4.5), "Petavius": (177.0, 3.3),
}


def main() -> None:
    lines: list[str] = []

    def say(s: str = "") -> None:
        print(s)
        lines.append(s)

    say("ExoGen Phase 0 spike: crater scaling")
    say("=" * 60)

    # Chronology anchors
    say("\nChronology: N(>=1 km) per km^2 by surface age")
    for t in (0.1, 1.0, 3.0, 3.3, 3.5, 3.8, 4.0, 4.2):
        say(f"  T = {t:4.1f} Gyr  N(1) = {chronology_n1(t):.3e}")
    say("  Expect: maria (3.2-3.5 Gyr) a few e-3, highlands (4.0-4.2 Gyr) 5e-2 to 1e-1 (near saturation)")
    say(f"  Internal consistency: NPF(1 km) = {npf_cumulative_1gyr(1.0):.3e} vs chronology N1(1 Gyr) = {chronology_n1(1.0):.3e}")

    # Local slope of the production function
    say("\nProduction function cumulative slope d log N / d log D")
    for d in (0.1, 1.0, 10.0, 100.0):
        s = (np.log10(npf_cumulative_1gyr(d * 1.1)) - np.log10(npf_cumulative_1gyr(d / 1.1))) / (2 * np.log10(1.1))
        say(f"  D = {d:6.1f} km  slope = {s:6.2f}")

    # Transition diameter scaling
    say("\nSimple-to-complex transition: predicted (F30) vs observed [km]")
    for name, (g, icy, dsc_obs, _, _) in BODIES.items():
        say(f"  {name:9s} g={g:4.2f}  predicted {d_sc_predicted_km(g, icy):6.2f}  observed {dsc_obs:5.1f}")

    # Atmospheric cut-off
    say("\nAtmospheric crater cut-off (F24), anchored at Earth 150 m [km]")
    for name in ("Earth", "Venus", "Mars", "Titan"):
        g, _, _, _, p = BODIES[name]
        say(f"  {name:6s} P/g = {p / g:9.3e} kg/m^2  D_min = {crater_cutoff_km(p, g):8.4f}")
    say("  Expect: Venus 1.5-3 km, Mars metres, Titan several km (atmosphere plus erosion)")

    # Pi scaling worked examples
    say("\nPi-scaling worked examples (Collins et al. 2005)")
    d_tc = transient_diameter(50.0, 12800.0, 9.81, rho_i=7800.0, rho_t=2500.0)
    say(f"  Meteor Crater: 50 m iron at 12.8 km/s -> transient {d_tc:.0f} m, final {final_diameter(d_tc, 3200.0):.0f} m (observed 1186 m)")
    d_tc = transient_diameter(10000.0, 20000.0, 9.81)
    say(f"  Chicxulub-class: 10 km rock at 20 km/s -> transient {d_tc / 1e3:.1f} km, final {final_diameter(d_tc, 3200.0) / 1e3:.1f} km (observed ~150-180 km)")
    say("  Impactor size needed for a 100 km crater, by body:")
    for name, (g, icy, dsc_obs, v, _) in BODIES.items():
        rho_t = 920.0 if icy else 2700.0
        rho_i = 600.0 if icy else 2700.0
        L = impactor_for_crater(100e3, v * 1e3, g, dsc_obs * 1e3, rho_i=rho_i, rho_t=rho_t)
        say(f"    {name:9s} v={v:4.1f} km/s  L = {L / 1e3:5.2f} km")

    # Morphometry check
    say("\nLunar depth-diameter check (Pike): predicted vs measured [km]")
    for name, (D, d_obs) in LUNAR_CRATERS.items():
        d_pred = depth_simple(D) if D < 16.5 else depth_complex(D)
        say(f"  {name:11s} D={D:6.1f}  predicted d={d_pred:5.2f}  measured d={d_obs:4.1f}")
    say(f"  Copernicus rim height predicted {rim_height_complex(93.0):.2f} km (measured about 1.2-1.5), central peak {central_peak_height(93.0):.2f} km (measured about 1.2)")

    # Ejecta
    say("\nEjecta thickness (McGetchin) for a 100 km crater [m]")
    for mult in (1.0, 1.5, 2.0, 3.0, 5.0):
        say(f"  r = {mult:3.1f} R  thickness = {ejecta_thickness(mult * 50e3, 50e3):8.1f}")

    # Degradation
    say("\nDiffusional degradation of a 1 km simple crater (lunar kappa = 5.5 m^2/Myr)")
    r = np.linspace(0.0, 3000.0, 601)
    h0 = fresh_profile(r, 1000.0)
    prof = {0.0: h0}
    for kt in (1000.0, 5000.0, 20000.0):
        prof[kt] = diffuse_radial(h0, r, kt)
    for kt, h in prof.items():
        depth = h.max() - h.min()
        say(f"  kappa*t = {kt:7.0f} m^2 (t = {kt / 5.5 / 1000:4.2f} Gyr)  rim-to-floor depth = {depth:6.1f} m  d/D = {depth / 1000:5.3f}")
    say("  Expect: fresh d/D 0.2 falling toward 0.05 and below over 3-4 Gyr for 1 km craters")

    # ------------------------------------------------------------------ plots
    d = np.logspace(-2, np.log10(300), 200)
    fig, ax = plt.subplots(figsize=(6, 5))
    for t in (0.1, 1.0, 3.0, 3.5, 4.0):
        ax.loglog(d, isochron(d, t), label=f"{t:g} Gyr")
    ax.set_xlabel("crater diameter D [km]")
    ax.set_ylabel("cumulative N(>=D) [km$^{-2}$]")
    ax.set_title("Neukum production function isochrons")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "crater_isochrons.png", dpi=130)

    fig, ax = plt.subplots(figsize=(6, 5))
    dd = np.logspace(-0.5, 2.5, 200)
    ax.loglog(dd[dd < 16.5], depth_simple(dd[dd < 16.5]), "k-", label="Moon simple (Pike)")
    ax.loglog(dd[dd >= 16.5], depth_complex(dd[dd >= 16.5]), "k--", label="Moon complex (Pike)")
    ax.loglog(dd, depth_mars_garvin(dd), "r-", alpha=0.7, label="Mars (Garvin, verify)")
    for name, (D, d_obs) in LUNAR_CRATERS.items():
        ax.plot(D, d_obs, "ko", ms=4)
        ax.annotate(name, (D, d_obs), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("D [km]")
    ax.set_ylabel("depth [km]")
    ax.set_title("Depth versus diameter")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "crater_depth_diameter.png", dpi=130)

    fig, ax = plt.subplots(figsize=(6, 5))
    for name, (g, icy, dsc_obs, _, _) in BODIES.items():
        ax.plot(g, dsc_obs, "s" if icy else "o", color="tab:blue" if icy else "tab:red")
        ax.annotate(name, (g, dsc_obs), fontsize=8, xytext=(4, 2), textcoords="offset points")
    gg = np.linspace(1.0, 10.0, 50)
    ax.plot(gg, d_sc_predicted_km(gg), "r-", label="rock: 16.5 km x 1.62/g")
    ax.plot(gg, d_sc_predicted_km(gg, icy=True), "b-", label="ice: x 0.15")
    ax.set_xlabel("surface gravity [m/s$^2$]")
    ax.set_ylabel("simple-complex transition [km]")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "crater_transition_vs_gravity.png", dpi=130)

    fig, ax = plt.subplots(figsize=(6, 4))
    for kt, h in prof.items():
        ax.plot(r, h, label=f"kappa t = {kt:g} m$^2$")
    ax.set_xlabel("radius [m]")
    ax.set_ylabel("height [m]")
    ax.set_title("Diffusional degradation of a 1 km crater")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "crater_degradation.png", dpi=130)

    (OUT / "crater_scaling_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    say(f"\nPlots and report written to {OUT}")


if __name__ == "__main__":
    main()
