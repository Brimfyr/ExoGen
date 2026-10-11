# Phase 4 report: ice worlds

Status: in progress, started 2026-10-11.

Phase 4 gives the icy bodies their own physics. Phase 3 left them as cratered balls of ice with a regime name: Europa a blank sheet, Ganymede one terrain, and both carrying twenty kilometres of basin relief that warm ice cannot hold. The plan's list (PLAN.md, Phase 4 and section 4.3) is composition and shell thickness from density, viscous relaxation of craters into palimpsests, the tidal stress field with lineae, cycloids and bands, grooved terrain, focused polar heating, chaos, cryovolcanism with plume fallout and salt deposits, sublimation lag, hemispheric asymmetry for locked moons, and the ice and salt end-members. The gate is that Ganymede shows two terrains with a crater density contrast, that Europa's crack orientations follow the tidal stress field, and that relaxation removes large craters first.

The principle set in the Phase 3 reviews holds throughout: counts come from budgets and distributions, kinds from probabilities, colours from materials and processes.

## Milestones

1. **Ice foundations.** The crater scaling's ice factor keyed on the surface material rather than the bulk composition (the Phase 2 note); sublimation lag as the maturity law of icy surfaces, from the vapour pressure of ice at the latitude's mean temperature and the regolith's non-ice share, so warm old ice darkens and fresh craters expose clean ice; cryolava floods on ice regimes from the ice magma budget.
2. **Viscous relaxation.** Topography on an ice shell relaxes by wavelength with the viscosity of the ice at the depth the wavelength probes, through the shell's thermal history, in the loop on the coarse grid; large craters flatten into palimpsests while small ones keep their depth; replayed craters relax by their age.
3. **Tidal stress and the cracked shell.** The degree-2 membrane stress of a thin shell over an ocean from the Love numbers, its diurnal and non-synchronous parts; lineae as double ridges along the most tensile direction at their time, cycloids where the diurnal stress turns the crack through an orbit, bands where extension dominates; salts in young cracks fading under frost.
4. **Partial resurfacing.** Grooved terrain on a shell that extended, swaths of parallel grooves in clean ice cutting old dark terrain, with the crater density contrast; the focused polar terrain of a small shell over a deep ocean, with its fissures and fallout.
5. **Chaos, cryovolcanism and asymmetry.** Melt-through chaos where the shell is thinnest; brine domes and bright salt deposits in young crater floors on a shell over rock; plume fallout; leading and trailing hemisphere tints on locked moons.
6. **Documentation and review.** Parameter rows and formulas, this report, the presets' descriptions.

## Gate results

- **Relaxation removes large craters first** (`stages::tests::ganymede_relaxes_its_large_craters_first`). On the Ganymede preset the median depth-to-diameter of craters 30 to 70 km across stays above 0.02 while that of craters over 150 km falls below a quarter of it. At standard quality: 20 to 40 km 0.072, 40 to 80 km 0.044, 80 to 160 km 0.001, over 160 km 0.000; Callisto, whose dirty ice is stiffer, keeps 0.060 at 40 to 80 km. Ganymede's relief fell from 30 km of basin to 11 km.
- The Ganymede two-terrain gate and the Europa crack-orientation gate wait for milestones 3 and 4.

## Measurements

Standard quality, after milestone 2.

| Body | Loop | Finish | What shows |
|---|---|---|---|
| Europa | 3.9 s | 3.9 s | a bright sheet of 0.83 reflectance, relief within 3 km, five craters |
| Ganymede | 6.6 s | 21 s | relief 11 km, large craters flat, frost caps past 40°, mean reflectance 0.65 against 0.83 at the poles (sRGB) |
| Callisto | 4.5 s | 11 s | dark lag at the equator (0.40) grading to frost caps (0.77), bright young craters |
| Ceres | 10 s | 59 s | dark mixture with cryolava on 2% of the surface, large basins relaxed |

## Findings during the phase

- **The lag without a brake went opaque in a few million years at Callisto's equator**, and its caps began at a cliff of latitude, because the free sublimation rate at 120 K times the age is thousands of tonnes of ice per square metre. The lag chokes its own growth once it is thick and grows as the square root of time, which is what Moore et al. 1999 describe, and the transition spreads over twenty degrees.
- **Every ancient cell was excavated clean ice.** The deep bedrock that Phase 2 gave icy bodies made Ganymede white once thousands of overlapping craters had spread their ejecta over everything. Gardening now mixes clean ice back into the regolith over half a gigayear, and a shell over rock has no clean ice below it at all.
- **A ductile shell alone flattened everything.** With the viscous law alone, even 40 km craters relaxed on Ganymede, because the shell was thin and warm early and a step of fifty million years is an age to warm ice. The cold elastic lid holds loads shorter than a few times its thickness without flow, and the rock in a dirty crust stiffens it by orders of magnitude (Fu et al. 2017 on Ceres), which keeps the small craters and the dirty-ice bodies' craters where they are observed.
- **Ceres is buried in its own basin ejecta.** The production function in the belt gives Ceres nine basins over 400 km, as the models in Marchi et al. 2016 also predict; the real body has none, having relaxed them. Ours relaxes their relief but keeps their ejecta blankets, which bury most of the older small craters. The crispness of Ceres's small craters waits for the Phase 5 look at ejecta burial on small bodies.
- **Enceladus keeps four kilometres of primordial relief** at scales the elastic lid holds; the real body has a kilometre. The primordial crust's amplitude on small bodies is a Phase 1 setting to revisit with the polar terrain in milestone 4.
