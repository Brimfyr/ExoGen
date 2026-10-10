# ExoGen: process-based planetary map generator

Project plan, revision 0.1, 2026-10-10. Working title taken from the folder name.

This document is the design baseline. It records what the software is, what it is not, the physical models each stage uses, how the stages fit together, how correctness will be checked, and the order in which it will be built. Numbers quoted in the model sections are starting values recalled from the literature. Every one of them is re-checked against the cited source during Phase 0 before it is coded.

---

## 0. Purpose and scope

ExoGen is a desktop application and command-line tool that generates full-planet texture sets for solid-surface bodies. Surfaces are produced by simulating the processes that shape real planets, scaled to the body's size, composition, thermal state, atmosphere and orbital environment. The outputs are equirectangular (plate carrée) maps at up to 16k: albedo, height, roughness, specular, plus masks and auxiliary layers, with a globe preview.

The defining design choice: terrain and albedo come from a sequence of physical processes acting over simulated time on a layered global field. Noise is used only to seed initial conditions and to add fine detail that is conditioned on the simulated macro-structure. This is what separates the result from a noise field.

### Target output set

| Map | Encoding | Notes |
|---|---|---|
| Albedo | Linear reflectance, EXR float16; sRGB PNG 8/16-bit convenience copy | Equal-energy illuminant. Star colour belongs in the renderer. |
| Height | EXR float32, metres relative to a datum (mean radius or sea level) | Also 16-bit PNG normalised, with min/max in metadata |
| Roughness | 8/16-bit grey, GGX perceptual roughness | |
| Specular | 8-bit grey F0 | Water 0.02, ice 0.018, rock 0.04 |
| Normal | 8/16-bit tangent space | Gradients computed on the sphere, not on the image |
| Masks | 8-bit | Ocean, ice, snow, cloud, lava, dust, vegetation |
| Cloud | One map set per cloud class: optical depth, base and top height; cloud type ID; optional composite, shadow bake, sequences | From the weather and cloud stage. Classes in section 4.7 |
| Atmosphere | Sidecar parameters only | Derived scattering and aerosol parameters for the renderer's own atmosphere. No maps |
| Unit ID | 16-bit integer | Geological unit per pixel, with a legend in the sidecar |
| Auxiliary | float | Temperature, precipitation, surface age, crust age, slope, wind |
| Sequence | Any of the above at a list of epochs | Formation sequence for animation; CLI batch |
| Sidecar | JSON | All parameters, derived values, seeds, software version, legend, epoch |

Single-file export to 16384 x 8192. Tiled export to 32768 x 16384. Cube-map export as an option.

### Non-goals

- 3D terrain meshes, flyover rendering or an in-app renderer beyond the preview globe.
- Gas giants, ice giants, rings.
- Irregular small bodies. Below roughly 400 km radius the spheroid assumption weakens. Smaller radii are allowed with a warning.
- A general circulation model or a mantle convection model. Reduced models are used wherever they produce the right morphology.
- Civilisation features.
- Atmospheric maps. Haze and background aerosol are properties of the atmosphere, not layers that vary across the surface. They are described by derived parameters in the sidecar for the renderer's own scattering term. The scattering toggle on the globe is a preview, and nothing in the export depends on it.
- Browser deployment.

---

## 1. Assessment of the concept

### What works

1. **Process pipeline instead of a noise stack.** Real surfaces are the product of a small number of dominant processes with known scaling laws: impact cratering, tectonics, volcanism, erosion, volatile transport and weathering. Each is well documented. Simulating them in order, over time, is the only route to the "ancient highlands next to young plains" look that noise cannot produce.

2. **Derived parameters with unlock.** A physical derivation graph gives coherence: star and orbit give temperature, temperature and inventory give phase states, heat budget gives tectonic regime, atmosphere gives crater cut-off and weathering. Letting the user override any node gives control without losing the chain downstream. This is the right model and section 3.3 specifies it.

3. **Scale laws as a first-class requirement.** Crater morphology transitions, maximum relief, ejecta range and rift widths all scale with gravity, target strength and body size in known ways. Building the validation suite around published scaling data (section 6) is what will make a 1,500 km moon look different from a 6,000 km planet rather than just smaller.

4. **Composition-based colour.** Assigning each cell a mixture of real materials with measured reflectance spectra, then integrating to RGB, gives correct colours for free and makes "atmosphere changes surface colour" a weathering process rather than a tint.

5. **Signature features as outcomes, not decorations.** The Tharsis to Valles Marineris chain, lunar maria filling the thin-crust basins, Hellas with a bright dusty floor, Titan's dark equatorial dune belt and Pluto's dark equatorial tholin band are each the consequence of one process acting on the output of another. Section 5 lists the chain for each reference body and the pipeline is designed so that those chains run naturally.

### What needs reframing

1. **"Photorealistic on any computer."** Fidelity costs memory and time. The design sets a hardware floor (4 cores, 8 GB RAM, any Vulkan/DX12/Metal GPU optional) and separates simulation resolution from export resolution. Simulation runs at 2k to 8k equivalent; export to 16k is a conditioned detail pass, tiled. Everything has a CPU path. See section 9.

2. **Tectonics and climate fidelity.** Full mantle convection and a GCM are out of reach and unnecessary. A kinematic plate model (plates as rigid caps with Euler poles, boundaries classified by relative motion, crust tracked with age and type) reproduces mountain belts, arcs, rifts, seafloor age-depth and hotspot chains. A two-dimensional energy-balance climate model with parametrised circulation reproduces Hadley cell width, deserts at the right latitudes, rain shadows, polar and hemispherical ice, and gives a time-averaged wind field for dust. A two-layer shallow-water model on the sphere supplies the weather snapshot that clouds are built from: cyclones, fronts, jets, superrotation and the substellar cloud of locked worlds emerge from it rather than being drawn. All three are stated plainly in the UI as reduced models.

3. **Vegetation.** Plausible at the biome level: temperature and precipitation give a Köppen-type class, each class has a coverage and a pigment reflectance. That is what Earth looks like from orbit. Anything finer is out of scope. It stays a late phase because it depends on the climate stage being good.

4. **Parameter space consistency.** Much of the parameter space is physically inconsistent (nitrogen seas at 300 K, a 90 bar atmosphere on a 0.1 g body). The plausibility checker warns and explains. It never blocks. Users will deliberately make inconsistent worlds and that is their call.

5. **Exotic chemistry colours.** Only a handful of surface colouring mechanisms are understood well enough to model: iron oxidation, tholin deposition, sulfur allotropes and SO2 frost, salts and evaporites, space weathering, volatile frosts. For atmospheres outside those regimes the user supplies colour rules through the material system rather than the software inventing chemistry.

6. **"Nothing like this exists."** Mostly true for the combination, but there is prior art to learn from and it is listed so we do not reinvent it: tectonics.js and PlaTec (kinematic plates), WorldEngine (plates plus simple climate), Songs of the Eons (plates and climate for worldbuilding), Cortial et al. 2019 "Procedural Tectonic Planets" (the closest academic precedent for the plate stage), FastScape (fast fluvial erosion), SpaceEngine (noise-based, closed). None of them combine process-based cratering with scaling validation, multiple body classes, volatile chemistry and studio outputs. That gap is the product.

7. **Terminology.** The maps you called "equatorial" are equirectangular. The document uses that word throughout. Height, roughness and specular are the inputs to a physically based renderer, and the photorealism of the final image depends on that renderer as well. We document the intended shading model so the maps are used as designed.

---

## 2. Design principles

1. Every surface feature has a cause in the pipeline. If a feature cannot be explained by a process, it does not go in.
2. Every physical parameter has a derivation, a unit and a reference. The UI shows all three.
3. Derived by default, overridable per node, with artistic multipliers layered on top so that creative control never breaks the derivation chain.
4. Simulate on a sphere. No equirectangular artefacts in simulation. Projection happens at export.
5. Resolution-independent processes. The same seed and parameters produce the same world at any resolution, with more detail at higher resolution.
6. Deterministic. Seed plus parameters plus version reproduces the output bit for bit on the CPU path.
7. Each process is a stage with a declared input and output set. Stages cache. Changing a late parameter does not re-run early stages.
8. Validation against real bodies is part of the test suite, not a one-time exercise.
9. UI copy is technical, unit-labelled and free of adjectives.
10. The simulated history is an output. Any epoch can be previewed, rendered or exported as a sequence.

---

## 3. Architecture

### 3.1 Modules

```
exogen-core      Simulation engine library. No UI dependency. CPU path mandatory, GPU path optional.
  params/        Parameter graph: definitions, derivations, lock state, plausibility rules
  grid/          Cubesphere at any resolution; Fibonacci point set and spherical Voronoi for plates; neighbour tables; resampling; projection
  fields/        Layered field storage, tiling, caching, hashing
  stages/        One module per process (section 4)
  materials/     Spectral end-members, mixing, colour integration
  export/        EXR, PNG, TIFF, cube map, sidecar
  validate/      Scaling tests, statistical tests, reference presets
exogen-cli       Headless: render, render a sequence of epochs, derive, validate, preset management, batch
exogen-app       Desktop UI (parameter panels, preview, globe, pipeline control)
exogen-py        Optional thin bindings for analysis notebooks. Not a production path.
```

The engine is a library first. The GUI and CLI are thin clients. Studios will want the CLI for batch and farm use regardless.

### 3.2 Data model

**Grids.** Two grids, each used where it fits.

- *Cubesphere* for everything that touches the surface at resolution: cratering, erosion, deposition, materials, detail synthesis. Six faces with the equiangular mapping (cell area ratio about 1.3 centre to corner). Face size F gives an equirectangular-equivalent width of about 4F. Neighbour lookup across edges and corners via precomputed tables. GPU friendly: six 2D textures.
- *Coarse spherical point set* (Fibonacci sphere, 10k to 100k points) with its spherical Voronoi tessellation, for plate tectonics only: plates are regions and crust is advected as points, which a point set handles without remeshing. Rasterised to the cubesphere by spherical interpolation at each loop step. The climate, weather and thermal models run on coarse cubespheres (F=64 to 256) with the same grid code as the surface (ADR 0001).

Equirectangular is an export projection only. Poles are oversampled on export, which is harmless. Gradients for normals and slopes are computed on the cubesphere.

**Layered fields** (per cubesphere cell). This is what lets processes interact correctly.

| Layer | Type | Purpose |
|---|---|---|
| bedrock | f32 m | Solid surface elevation |
| regolith | f32 m | Loose material thickness: impact gardening, weathering products, sediment |
| dust | f32 m | Fine mobile deposit (aeolian) |
| ice | f32 m | Surface volatile ice: caps, glaciers, frost |
| liquid | f32 m | Standing liquid level above bedrock (seas, lakes) |
| crust_thickness | f32 km | Drives isostatic elevation |
| crust_type | u8 | Continental, oceanic, primary, resurfaced |
| crust_age | f32 Myr | Formation age |
| surface_age | f32 Myr | Time since last resurfacing: controls crater density, maturity |
| composition | 8 x f16 | Fractions of bedrock and deposit end-members |
| maturity | f16 | Space weathering or chemical weathering state |
| unit_id | u16 | Geological unit label |
| wind, temp, precip | coarse | From climate, interpolated |

Output height = bedrock + regolith + dust + ice, with the liquid surface optionally included. Masks derive from layer thickness. Memory at F=2048 (8k equivalent, 25M cells) is about 100 MB per f32 layer, so the full set is 1 to 2 GB. F=1024 is the default simulation size.

### 3.3 Parameter graph

Every parameter is a node with: id, group, unit, type, range, documentation string with reference, a derivation function over upstream nodes (or none for root inputs), and a lock state.

Lock states:

- **Derived.** Value computed from upstream. Shown dimmed with a lock glyph and the derivation inputs on hover. Default for every non-root node.
- **Set.** User value. Shown in full contrast with an "override" marker and a "reset to derived" action. Downstream nodes continue to derive from the set value.

On top of the physical nodes sit **multipliers**: artistic controls that scale a derived quantity without unlocking it (crater density x, erosion intensity x, volcanic output x, dust cover x, relief x). They default to 1 and are always available. This two-tier scheme gives artists "more craters" without asking them to pick an impact flux and a retention age.

Expert mode reveals internal nodes that are otherwise hidden: production function slope, diffusivity, lithosphere yield stress, stream power constants, plate count, and so on.

Plausibility rules are predicates over nodes with a severity and an explanation string. Examples: liquid nitrogen requested but surface temperature exceeds 77 K at 1 bar; atmosphere unlikely to be retained (Jeans escape parameter below 6 for the main gas); tidal locking time below body age but rotation set unlocked; requested relief exceeds the crustal strength limit. Warnings appear in a panel and against the parameter.

Seeds follow the same rule. One master seed derives a sub-seed per stage. Any sub-seed can be set, so the user can re-roll craters without changing continents.

**History curves.** A few root inputs can vary over simulated time: stellar luminosity (derived from stellar mass and age by default), orbital eccentricity (and so tidal heating), atmospheric surface pressure, volatile inventory, and the impact flux. Each is a piecewise-linear curve over time with a constant or derived default. The coarse loop samples them at every step. This is how an early thick, wet atmosphere that thins to the present one produces valley networks that are later cratered and dust-covered, and how a moon whose eccentricity was pumped and then damped shows a tectonic episode frozen under younger craters.

### 3.4 Pipeline and time model

The pipeline has a coarse forward loop and a high-resolution finishing pass.

**Coarse forward loop** at F=256 to 512 over the body's age (default 4.5 Gyr) in N steps (default 90 steps of 50 Myr). Each step, in order:

1. Sample the history curves for this step. Interior and heat budget update (radiogenic decay, tidal heating from the eccentricity history, secular cooling).
2. Tectonics step (regime-dependent: plates, contraction, loading, ice-shell stress).
3. Volcanism step (events sampled from the magma budget, flows filled by topography).
4. Large impacts for this interval (craters above a size threshold, flux declining with time).
5. Climate update every k steps, if there is an atmosphere.
6. Coarse erosion and deposition step.
7. Bookkeeping: ages, unit labels, composition. Keyframe written.

**Finishing pass** at full simulation resolution (F=1024 default):

1. Resample coarse state.
2. Emplace the small-crater population statistically: ages sampled from the flux history and the per-cell surface age, each crater pre-degraded for its age (section 4.5).
3. Erosion refinement: a bounded number of iterations with intensity from the recent erosion rate field.
4. Aeolian deposit equilibrium, weathering, materials, colour. Weather model and cloud layers at the output epoch (section 4.7).
5. Detail synthesis to export resolution, tiled (section 4.11).
6. Export.

The split keeps the expensive processes at coarse resolution where their interactions matter, and the resolution-dependent processes at full resolution where their statistics matter.

**Where the split falls.** A crater, flow or fault is simulated inside the loop when it spans at least about eight loop cells; below that it is emplaced in the finishing pass. The loop resolution is therefore a quality dial, not an accuracy ceiling. At F=256 on a Mars-sized body the loop cell is about 21 km, so craters above roughly 170 km are simulated and everything smaller is replayed. At F=1024 the cell is about 5 km and craters above about 40 km are simulated. The dial trades time for the size of the truly simulated population. The replay uses the same rate fields the loop produced, so the two populations degrade consistently.

**The finishing pass is an age-ordered replay.** Every sub-threshold feature (small craters, scarps, flows, cones, dune fields) carries a formation age. The pass emplaces them oldest first, so superposition and cross-cutting come out right: a younger crater's ejecta buries an older rim, a scarp cuts the craters older than it and is buried by flows younger than it. Each feature is degraded by the process history of its cell since its own age (section 4.5) and warped by the deformation accumulated since then. It is a replay of the same history at higher resolution, with integrated rates in place of iterated time steps.

**Timeline and epochs.** The loop keeps a keyframe of its state at every step. The preview has a time scrubber: each keyframe gets a cheap finishing pass at preview resolution on demand, a few seconds, cached. This gives a scrubbable evolution of the body at 1k including the small-crater population and colour. An **output epoch** parameter sets which keyframe the full-resolution finishing pass starts from, so any past state can be exported at full quality for the cost of one normal render. The CLI exports sequences of epochs. What is not offered is live scrubbing at export resolution; the finishing pass is minutes, not frames.

Each stage declares its inputs. Results are cached keyed by the hash of the upstream parameter set and seeds, so tweaking the material palette does not re-run tectonics.

### 3.5 Project files and determinism

Project file: TOML. Sections for parameters with lock states, multipliers, seeds, output settings, preset lineage and the software version that wrote it. Presets are project files.

Determinism: the CPU path is bit-exact for a given version. The GPU path is visually identical but not guaranteed bit-exact across vendors. The sidecar records which path produced each map.

---

## 4. Physical models by stage

### 4.0 Body and environment derivation

Root inputs: star (mass or spectral type, luminosity, age), orbit (semi-major axis, eccentricity), parent planet if the body is a moon (mass, distance), body (radius and mass, or radius and density), rotation (period, obliquity, or locked), atmosphere (surface pressure, composition by volume), volatile inventory (mass of H2O, CO2, CH4, N2, NH3, SO2 available at the surface), age.

Derived nodes and their basis:

| Node | Basis |
|---|---|
| Surface gravity, escape velocity | g = GM/R², v_esc = sqrt(2GM/R) |
| Bulk composition (core, mantle, ice mass fractions) | Fit of (M, R) to three-layer mass-radius relations (Zeng et al. 2016, Fortney et al. 2007). Mercury 5.43, Earth 5.51, Mars 3.93, Io 3.53, Europa 3.01, Ganymede 1.94, Titan 1.88, Pluto 1.85 g/cm³ as anchors |
| Composition class | Rock-iron, rock, rock-ice, ice-dominated, from ice fraction and surface temperature |
| Equilibrium temperature | T_eq⁴ = L(1 - A_B) / (16 π σ a²) full redistribution; substellar T = √2 T_eq without redistribution |
| Greenhouse offset | Grey radiative equilibrium T_s⁴ = T_eq⁴ (1 + 3τ/4) with τ from column mass and per-gas opacity (CO2, H2O, CH4, N2 collision-induced). Calibrated on Earth (+33 K), Mars (+5 K), Venus (+500 K), Titan (net +3 K with anti-greenhouse) |
| Tidal heating | Ė = (21/2)(k₂/Q) G M_p² R⁵ n e² / a⁶ (Peale, Cassen and Reynolds 1979); k₂/Q from composition and the presence of an internal ocean |
| Tidal locking time vs age | t_lock ≈ ω a⁶ I Q / (3 G M_p² k₂ R⁵) (Gladman et al. 1996); locked if t_lock < age |
| Radiogenic heat | Rock mass fraction x chondritic heat production decaying with age |
| Surface heat flow | Radiogenic + tidal + secular cooling, per unit area. Earth 0.087 W/m², Io about 2 W/m² |
| Lithosphere and elastic thickness | Conductive thermal thickness from heat flow; elastic thickness from the 600 to 750 K isotherm for rock, 150 K-scaled for ice |
| Tectonic regime | Mobile lid, episodic lid, stagnant lid, ice shell over ocean, ice shell over rock, from heat flow, lithosphere strength, surface water and composition |
| Maximum sustainable relief | Crustal strength over (ρ g), temperature-weakened; Earth about 9 km, Mars about 22 km, Venus 11 km. For ice, relaxation-limited rather than strength-limited |
| Volatile phase states | Vapour pressure curves and triple points per species at (T_s, P_s): H2O 273 K / 612 Pa, CO2 216.6 K / 518 kPa, CH4 90.7 K / 11.7 kPa, N2 63.1 K / 12.5 kPa, NH3 195.4 K / 6.1 kPa. Liquid exists only above the triple-point pressure. Ammonia-water eutectic at 176 K for cryolava |
| Atmospheric retention | Jeans parameter per gas; warning below about 6 |
| Impact velocity | v² = v_∞² + v_esc², with v_∞ from orbital location and parent focusing |
| Crater cut-off diameter | From atmospheric column mass P/g (section 4.5) |
| Space weathering rate | Scaled by 1/a² for solar wind and micrometeorites; zero under a dense atmosphere |
| Scale height | H = kT / (μ g); used by the climate and scattering preview |
| Stellar luminosity history | Main-sequence evolution from stellar mass and age (a Sun-like star was about 70% as bright 4.5 Gyr ago); sampled by the loop's climate steps |

### 4.1 Primordial crust

Initial crustal thickness field: low-degree spherical harmonics (l ≤ 8) with an explicit l=1 term for a hemispheric dichotomy (Mars strong, Moon moderate, Earth weak) plus noise. Option: a giant impact as the dichotomy source (a Borealis-type basin, Andrews-Hanna et al. 2008). Elevation follows from Airy isostasy: a 35 km continental crust at Earth densities gives about +0.8 km, matching Earth's mean land height; a 30 km difference on Mars gives the 5 km dichotomy step. For mobile-lid worlds, continental nuclei (cratons) are seeded here.

### 4.2 Interior and heat budget

Tracked through the coarse loop: radiogenic production, tidal heating (constant unless eccentricity is scripted), mantle temperature with secular cooling, lithosphere thickness, cumulative radial contraction ΔR, magma production rate as the excess of heat flow over the conductive capacity of the lithosphere. Outputs feed tectonics (regime, contraction strain, plume count), volcanism (magma budget), and cratering (relaxation viscosity for ice).

### 4.3 Tectonics

The regime node selects which of these run. Regimes are not exclusive: a stagnant-lid world can have both contraction and plume loading.

**Mobile lid (Earth-type).** Kinematic plates on the coarse point set. Plate count derived from area and heat flow, overridable. Each plate has an Euler pole. Steps of 1 to 5 Myr over the simulated interval, with periodic re-organisation events when plates are consumed or grow too large.

- Divergent boundary: new oceanic crust with age 0; ridge at 2.6 km depth; seafloor deepens as d = 2600 + 365 √t m (t in Myr) to 70 Myr then approaches 5.65 km (Parsons and Sclater 1977), scaled by mantle temperature. This alone produces realistic bathymetry.
- Convergent, oceanic under continental: subduction. Trench at the margin (flexural profile to -8 to -11 km on Earth), volcanic arc 100 to 300 km inboard, crustal thickening and uplift (Andean type).
- Convergent, oceanic under oceanic: island arc and trench, older plate subducts.
- Convergent, continental against continental: collision. Crust thickens toward double, plateau forms behind the range (Himalaya-Tibet type), thickening propagates inboard with time.
- Transform: strike-slip offset of existing features.
- Continental rifting: a divergent boundary through continental crust produces a graben with flank uplift and volcanism; continued spreading produces a new ocean with thinned, tapering passive margins. These margins are the continental shelves.
- Hotspots: fixed points in the mantle frame. A volcano forms over the plume; earlier volcanoes are carried away and subside with crust age, producing seamount chains.
- Orogen age is tracked. Young belts are sharp and high; old belts are low and rounded after the erosion pass.

Study PlaTec and Cortial et al. 2019 before implementing; both solve the advection and re-rasterisation problem.

**Stagnant lid (Mars, Mercury, Moon type).**

- Global contraction: cumulative ΔR from cooling gives strain ε = 2ΔR/R (Mercury: 5 to 7 km, about 0.5%). Expressed as lobate thrust scarps with a power-law length distribution and total shortening equal to the strain, asymmetric profile (steep face, gentle back slope), height scaling with length (Enterprise Rupes: 1,000 km long, 3 km high). Scarps cut and deform existing craters.
- Plume loading: one or more volcanic rises (Tharsis: 5,000 km across, 7 km high) as broad domes. The thin-shell flexural stress field around the load gives radial extension on the flanks and concentric compression at the periphery. Faults form where stress exceeds yield: radial graben sets, concentric wrinkle ridges. One or a few graben can be selected (derived from the load size, overridable) to grow into a canyon system by wall collapse and erosion. Valles Marineris: 4,000 km long, 200 km wide, 7 km deep. Shield volcanoes then form on the rise, sized by the relief limit and magma budget, with nested calderas and landslide aureoles.
- Flood basalts: see volcanism.

**Episodic lid (Venus type).** A global resurfacing event at a chosen time resets the surface age everywhere except a retained fraction of high-standing deformed terrain (tessera, about 8% on Venus). Plains carry flow fields, coronae (annular ridges 100 to 1,000 km across with a central dome or depression), pancake domes, long lava channels. The crater population after the event is sparse, uniformly distributed and un-degraded. Albedo is nearly uniform basalt; the visual character lives in height and roughness.

**Ice shell (Europa, Ganymede, Enceladus type).**

- Tidal stress field from diurnal and non-synchronous rotation terms (Greenberg et al. 1998) gives crack orientation everywhere, including cycloids.
- Lineae along cracks: double ridges with a central trough, 1 to 2 km wide, hundreds of km long. Age sequence: new lineae cross old ones. Composition: hydrated salts, brown-red when young, fading as frost accumulates.
- Bands: dilational spreading zones 10 to 30 km wide where extension dominates.
- Chaos: melt-through regions where local heat flux is highest; rafted polygonal blocks in a hummocky matrix.
- Grooved terrain (Ganymede): extensional resurfacing bands of parallel grooves at 5 to 10 km spacing, bright fresh ice, cutting dark cratered terrain. Resurfaced fraction is a parameter (Ganymede about two thirds).
- Focused heating (Enceladus): a polar terrain with a few parallel fissures (130 km long, 35 km apart), plume fallout brightening the surrounding region, young terrain there and ancient cratered terrain elsewhere.
- Viscous relaxation of craters: relaxation time τ = 4π η / (ρ g λ) for wavelength λ, with η strongly temperature dependent. Large craters relax first, producing palimpsests and shallow, flat-floored forms.

### 4.4 Volcanism

Rate from the magma budget. Distribution from the tectonic stage: arcs, rifts, hotspots, plume rises, or (stagnant lid without plumes) scattered vents and fissures.

- Shield volcanoes: low slopes (4 to 6°), diameter set by magma supply and relief limit. On low-g bodies they are enormous (Olympus Mons: 600 km across, 22 km high). Calderas as nested collapse pits.
- Stratovolcanoes at arcs: steeper cones, 20 to 40 km across, 2 to 4 km high.
- Flood basalts and plains: eruption volume filled by topographic flooding from the vent (priority-flood fill), so lava follows terrain. Basin floors flood first. Older craters partly flooded keep their rims where the rim stands above the fill (ghost craters). Each flooding episode resets the surface age and composition of the filled cells. Lunar maria occur where crust is thin and basins exist, so the nearside-farside asymmetry emerges from the dichotomy.
- Wrinkle ridges on plains: compressional, concentric to basin loads.
- Lava channels and rilles.
- Io regime (heat flow above about 1 W/m² on a rocky body): resurfacing outruns impacts, so no craters survive. Hundreds of paterae (10 to 200 km, Loki 200 km) concentrated at low latitudes. Fresh silicate flows are black, then whiten under SO2 frost over the resurfacing timescale. Plume deposits form rings: red short-chain sulfur near the vent, white SO2 frost further out. Mountains are tectonic blocks (steep isolated massifs up to 17 km) from crustal subsidence under burial, not volcanic edifices. Between the Io extreme and a dead world the model is continuous: the ratio of resurfacing rate to crater production sets how many craters survive.
- Cryovolcanism (icy bodies with heat): ammonia-water or brine flows and domes, smooth resurfaced plains, plume fallout brightening, bright salt deposits in crater floors (Ceres-type). Relaxation and tectonics interact with it.

### 4.5 Impact cratering

**Production function.** Cumulative size-frequency distribution per unit area per unit time. Inner solar system: Neukum production function (Neukum, Ivanov and Hartmann 2001), with chronology N(≥1 km) = 5.44e-14 (e^(6.93 T) - 1) + 8.38e-4 T per km² for surface age T in Gyr. Outer solar system: Zahnle et al. 2003 flux. The top-level "crater density" multiplier scales the whole function. Flux history (early heavy bombardment declining to a constant) is a parameter; a late spike is an option.

**Impactor to crater.** Pi-group scaling in the gravity regime (Holsapple 1993) in the form used by Collins, Melosh and Marcus 2005: transient diameter D_t = 1.161 (ρ_i/ρ_t)^(1/3) L^0.78 v^0.44 g^(-0.22) sin^(1/3)θ. Final rim diameter: 1.25 D_t for simple craters; 1.17 D_t^1.13 / D_sc^0.13 for complex (Croft collapse). Impact velocity from section 4.0.

**Morphology scaling.** Simple-to-complex transition diameter D_sc scales as 1/g for silicate targets (Moon 15 to 18 km, Mercury about 10 km, Mars 6 to 7 km, Earth 3 to 4 km). Icy targets transition at several times smaller diameter for the same g (Ganymede 2 to 3 km). Peak-ring onset about 140 km on the Moon, multi-ring onset about 300 to 400 km, both scaling similarly. Ring spacing ratio about √2.

Dimensions (lunar baseline, D in km, from Pike 1977 and 1980; other bodies get fitted coefficients in Phase 2):

- Simple: depth 0.196 D^1.01, rim height 0.036 D^1.01, bowl profile, depth-diameter about 0.2.
- Complex: depth 1.044 D^0.301, rim height 0.236 D^0.399, flat floor, terraced walls, central peak height 0.032 D^0.9. Volatile-rich targets (Mars) are shallower.
- Peak ring and multi-ring basins: inner ring, rim massifs, radial sculpture, antipodal disrupted terrain for the largest.
- Ejecta: thickness δ(r) = 0.14 R^0.74 (r/R)^-3 (McGetchin et al. 1973), continuous to about one radius beyond the rim, rays to 10 to 15 radii for fresh craters. Rampart (fluidised) ejecta where the target holds volatiles. Secondary chains and clusters outside the continuous ejecta. Oblique impacts give elongated forms and butterfly ejecta at low angles (rare, sampled from the impact angle distribution).

**Time integration.** Craters above the loop threshold are emplaced inside the forward loop at the step they form and are then acted on by everything that follows: diffusion, fluvial breaching, lava flooding (ghost craters), dust infill, scarp cross-cutting, plate transport and deformation, subduction, later impacts. Their position, diameter, type and age are also recorded as features so the finishing pass can add their resolved morphology (terraces, peak ring, rim crest detail, ejecta texture) on top of the evolved low-frequency form carried by the coarse grid. Sub-threshold craters are emplaced in the finishing pass in age order, with ages sampled from the flux history and capped by the surface age of the cell they land on. This is what produces saturated ancient terrain beside sparsely cratered young plains. On a mobile-lid world the age cap alone removes nearly every crater from oceanic crust, and the erosion integral removes most of the rest from wet continents, which is why an Earth-like comes out almost crater-free, as it should.

**Degradation.** For simulated craters, degradation is whatever the subsequent loop steps did to them. For replayed craters it is computed from the cell's process history between the crater's age and the output epoch, using the rate fields the loop produced:

- Diffusion: pre-blur the crater profile with σ² = 2∫κ dt over the crater's lifetime, which is exact for linear diffusion. Lunar κ is about 5.5 m² per Myr (Fassett and Thomson 2014), scaled by regolith production and gravity.
- Infill: floor raised by the deposition integral since the crater's age (dust, later ejecta, lava). Lava that arrives after the crater partly buries it and leaves a ghost rim where the rim stands above the fill. Lava before it leaves the crater fresh on top.
- Fluvial breaching and floor sedimentation proportional to the integrated fluvial intensity since the crater's age, followed by the full-resolution refinement iterations, which carve real channels through rims where drainage demands.
- Deformation: the crater is stamped in the coordinates of its formation time and pushed through the displacement accumulated since then. Rigid plate motion moves it whole. A rift that opened after it tears it in two. A transform offsets it. A contraction scarp younger than it cuts it. Shortening in an orogen flattens it. The loop stores per-step plate assignments and a coarse non-rigid displacement field so this composition is cheap.
- Overprinting by younger craters comes from the age ordering.
- Sublimation degradation on ice (knobby terrain, Callisto type) and viscous relaxation (section 4.3), from the integrated temperature history.
- Space weathering darkens and reddens exposed surfaces with maturity m = 1 - exp(-t/τ_sw), τ_sw of order 1 Gyr at 1 AU. Ray contrast fades as (1 - m).

The result is a continuous degradation sequence by age: fresh rayed craters on top, subdued rimless rings beneath, and basins visible only as circular lows and arcuate massifs. That sequence is what crater-counting geologists read, and it is the main visual signature of an old surface.

**Atmospheric filtering.** Minimum crater diameter grows with atmospheric column mass P/g. Anchors: Venus 1.5 to 3 km, Earth 100 to 200 m for stony impactors, Mars a few metres, Titan several km. A power-law fit on P/g with exponent near 0.6 is the starting point.

**Giant basins.** Sampled from the production function or placed explicitly by the user (count, diameter, location, age). Multi-ring structure, rim massifs, ejecta sculpture, antipodal terrain, subsequent lava fill and later cratering all follow automatically.

### 4.6 Climate

Runs on a coarse cubesphere (F=64 to 128) whenever surface pressure is non-negligible. Also runs in a degenerate mode for airless bodies to give surface temperature for volatile stability.

**Insolation.** Computed for arbitrary rotation state: period, obliquity, eccentricity, and locked. Annual and diurnal averaging by rotation regime. Obliquity above 54° puts the annual insolation maximum at the poles, which moves ice to the equator. Locked bodies get a substellar-centred field.

**Energy balance.** Two-dimensional diffusive energy balance (North 1975 form): heat capacity by surface type, insolation times co-albedo, outgoing longwave A + BT generalised through the greenhouse optical depth, horizontal transport D∇²T with D scaled by pressure and rotation rate following Williams and Kasting 1997. Ice-albedo feedback included. Ocean heat transport added where seas exist.

**Downscaling.** Climate fields are produced on a coarse cubesphere but are never applied to the surface at that resolution. Temperature is downscaled to the surface grid with the lapse rate on the fine elevation plus aspect and slope terms. Precipitation is downscaled with orographic lift from the fine slope and the wind. All thresholds are evaluated on the downscaled fields (section 4.12).

**Circulation.** Hadley cell edge latitude from Held and Hou 1980, scaling as (g H Δ / (Ω² a²))^(1/2). Fast rotators get multiple cells and jets; slow rotators and locked worlds get a single cell from substellar to antistellar or from equator to pole. Surface wind field: zonal bands with the right sign per cell, meridional components, deflection by rotation, and a thermal-gradient component for locked worlds. This wind field drives moisture transport, dust and clouds.

**Moisture and precipitation.** Saturation vapour pressure by Clausius-Clapeyron for the working volatile (water, methane, ammonia). Evaporation from liquid surfaces, advection along the wind field, condensation where converging or lifted orographically (vertical velocity from wind dotted with slope), depletion downwind. This gives rain shadows, subtropical deserts under descending branches, equatorial maxima, and wet windward coasts.

**Clouds.** Produced by the weather and cloud stage (section 4.7), which runs a dynamical model at the output epoch. The climatological circulation above serves the time loop; the weather stage produces the snapshot winds and the cloud layers.

**Biomes (Phase 8).** Köppen-type classification from temperature and precipitation. Each class carries coverage and pigment reflectance. Pigment is a parameter defaulting to chlorophyll.

### 4.7 Weather and clouds

Clouds are the largest single visual element on many worlds and the first thing that gives a planet away as procedural. The approach has three levels: a dynamical weather model that produces the synoptic structure, a classification that turns its fields into amounts of each cloud type, and detail recipes that give each type its texture at export resolution, conditioned on the fine terrain. Nothing is drawn where the dynamics do not call for it.

The stage runs in the finishing pass at the output epoch. The climatological circulation of section 4.6 serves the time loop; this stage produces the snapshot.

**Weather model.** A two-layer moist shallow-water model on the coarse cubesphere (F=128 to 256, about 1° to 0.5°). Rotation from the body's actual rotation state, including the orbital rate for locked bodies. Bottom topography from the smoothed elevation, so ranges block, deflect and launch stationary waves. Thermal forcing by relaxation toward a radiative-equilibrium thickness set by the surface temperature field: zonal for rotators, substellar-centred for locked bodies, offset by the season parameter, with an optional diurnal cycle keyed to a subsolar longitude parameter. Bottom friction. One moisture tracer per working condensable: evaporation over liquid (weakly over ice), advection, saturation by Clausius-Clapeyron, condensation into cloud water with latent heating, removal as precipitation on a timescale. Spin up for about 60 model days until eddy statistics are steady, then take a snapshot or a sequence.

What emerges without being drawn: baroclinic eddies along the storm tracks for rapid rotators (the reason for two layers rather than one), jets and their meanders, a single thermally direct cell for slow rotators, equatorial superrotation and the eastward shift of the substellar cloud for locked bodies (Showman and Polvani 2011, Penn and Vallis 2017), the regime change from banded rapid rotators through Rhines rotators to slow rotators (Haqq-Misra et al. 2018), and stationary waves downstream of large ranges, including a Venus-type bow over a highland. Resolution limits are stated plainly: 1° resolves cyclones of a few thousand km but not hurricanes, the fine structure of fronts or any mesoscale texture. Those come from the next two levels. For a body with no atmosphere to speak of, the stage runs in a degenerate mode: plumes only.

**Cloud type classification.** Per coarse cell, the amount, height and composition of each type, from the dynamics and thermodynamics:

| Type | Where the fields put it | Height |
|---|---|---|
| Deep convection: clusters, anvils | Low-level convergence with conditional instability and moisture: the convergence zone at the thermal equator or substellar point, the warm sector of cyclones, warm land in the afternoon when the diurnal option is on | High tops |
| Frontal stratiform | Ascent along strong horizontal temperature gradients in the eddies; this draws the comma and its trailing band | Mid to high |
| Stratocumulus and stratus | Subsidence over a cool surface with a moist lower layer: cold eastern ocean margins at subtropical latitudes, sea-ice margins, cold-air outbreaks behind fronts | Low |
| Shallow cumulus | Weak subsidence over warm liquid with moderate humidity: the trade belts | Low |
| Cirrus | Upper-layer moisture, anvil outflow advected by the upper wind, jet streaks | High |
| Polar hood and stratospheric condensate | Upper-layer temperature below the condensation point of any atmospheric species over the winter pole, bounded by the polar vortex: water ice on a Mars-type, ethane, HCN and benzene ice on a Titan-type, CO2 ice where cold enough | High |
| Fog and valley stratus | Near-saturated surface air with light wind under an inversion: basins, canyons, cold land at dawn, air crossing cold water or an ice edge | Surface |
| Dust | Lifting rate from the aeolian stage times wind and availability. Regional storms where lifting is strong and coherent: frontal passages, the retreating cap edge in spring, deep basins. A global storm option raises background optical depth to a chosen value | Low to mid |
| Plumes | Active vents from the volcanism stage emit a column to a height set by vent energy, gravity and pressure, bent downwind into a trail whose length is wind speed times settling time. Composition sets colour: SO2 frost white, Triton-type particulates dark, water ice white | By vent |
| Tropical cyclones | Placed as parametric primitives where the dynamics allow: warm liquid surface above a species-specific threshold, 8 to 30° from the equator, low shear, seeded by the vorticity field. Count derived, overridable | Low to high |
| Planetary-scale wave | For superrotating thick atmospheres, an equatorial wave with tilted mid-latitude arms superimposed on the dynamically produced banding (the Venus Y). Contrast derived near zero in visible light, higher for a UV-like look | Deck |

Tropical cyclones and the Venus wave are the two places where a structure is placed rather than emergent. Both are labelled as such in the UI.

**What is not a map.** Maps carry what varies across the surface; parameters carry what is a property of the atmosphere. Photochemical haze (tholin on methane-rich worlds, sulfuric acid aerosol on SO2-rich ones) and the background dust loading of a Mars-type atmosphere are smooth aerosol layers, best rendered by the renderer's own Mie term with a density profile. ExoGen does not generate maps for them. The sidecar carries the derived parameters: Rayleigh coefficients from composition and pressure, scale height, and for each aerosol its optical depth, scale height, asymmetry parameter, colour and a latitude profile of any seasonal asymmetry. The globe preview's scattering toggle uses the same parameters. The dust class below holds only the spatially structured part: storms and fronts.

**Terrain coupling.** The terrain shapes clouds at two scales. In the dynamics it is the lower boundary: large ranges deflect flow, anchor stationary waves and set where storm tracks run. In the recipes the fine terrain controls the edges:

- Windward slopes carry cloud in proportion to the upslope wind component, with a sharp clear edge at the crest on the lee side. This is the hard edge along a range.
- Lee waves as parallel bands downstream of a crest, wavelength scaling with wind speed over stability, resolvable at 16k on an Earth-sized body.
- Cap clouds over isolated peaks; afternoon convection over high ground.
- Fog fills valleys and basins as a level set below the inversion height on the fine terrain, so its edge follows the terrain exactly. Canyon fog on a Mars-type world comes from this rule.
- Vortex streets downstream of islands and isolated mountains that penetrate a stratocumulus deck: alternating vortices at a spacing of four to five obstacle widths, decaying downstream.
- Convergence lines along coastlines from land-sea breezes, which trace the coast.
- Cloud streets starting at a coast or ice edge in cold-air outbreaks, aligned with the wind, so the edge of the street field is the coastline.
- Plume trails from vents, and dust sources tied to the dust and sand fields of the aeolian stage.

**Detail recipes.** Each type has a texture recipe at export resolution: cellular convection with closed and open cells for stratocumulus, streets and popcorn for shallow cumulus, clusters with anvils spreading along the upper wind for deep convection, banding and a dry slot inside the comma for frontal cloud, curl-noise streaks with fall-streak texture for cirrus, spiral streaks for hoods, lobate fronts with convective texture for dust storms, a narrow umbrella and dispersing trail for plumes, an eye, eyewall and log-spiral bands under a cirrus canopy for tropical cyclones. All textures are advected and oriented by the local wind and conditioned on the coarse amounts, so there is no cloud where the dynamics put none.

**Export: one map set per cloud class.** Cloud types are never merged into a single map by default, because they need different treatment downstream: a dust storm is not lit, animated or scattered like cirrus. Types are grouped into classes by material, phase and altitude band, and each class is exported on its own:

| Class | Contents | Material and phase | Maps |
|---|---|---|---|
| low | Fog, stratus, stratocumulus, shallow cumulus | Droplets or small crystals of the working condensable | Optical depth, base height, top height |
| mid | Frontal stratiform, altostratus, altocumulus | Mixed phase | Optical depth, base height, top height |
| convective | Deep convective towers, eyewalls | Mixed phase, tall | Optical depth, base height, top height |
| high | Cirrus, anvils, polar hoods, stratospheric condensate | Ice crystals | Optical depth, height |
| dust | Regional and global dust storms, dust fronts | Mineral particulates, colour from the surface composition at the source | Optical depth, top height, colour |
| plume | Volcanic and cryovolcanic plumes and trails | Particulates or condensate by vent composition | Optical depth, height, colour |

Optical depth is exported linear in float16 EXR, with an 8-bit alpha convenience copy (1 − e^−τ). Heights are in km above the datum. A cloud type ID map records the dominant type per pixel. The sidecar lists discrete features (tropical cyclones with centre and radius, storms, plumes with their vent) and, per class, the species, phase, particle size regime, colour in linear RGB and typical altitude range, so a renderer can set its own phase function and placement. A single composite RGBA map with a height map remains available as a convenience, built by compositing the classes top down with their colours. A cloud shadow bake for a given sun direction is an option. Because the weather model is time-stepped, a sequence is cheap: N snapshots at a chosen interval, per class, for animated layers.

**Cloud regimes by world type.** What the stage should produce for the cases that matter:

| World | Expected result |
|---|---|
| Earth-like | Global cloud fraction near two thirds; a broken convective band near the thermal equator; trade cumulus fields; stratocumulus decks off cold subtropical west coasts; trains of comma-shaped cyclones with frontal bands along both storm tracks; cirrus along the jets; clear subtropical highs over the deserts; tropical cyclones over warm oceans in season; afternoon convection over warm continents |
| Mars-like | Sparse water-ice cirrus; a winter polar hood bounded by the vortex; an equatorial ice-cloud belt near aphelion; afternoon orographic clouds over the big shields with lee streaks; morning fog in canyons and deep basins; dust storms at the cap edge in spring and in the deep basins, with a global storm option |
| Venus-like | A complete deck; dynamically produced banding under superrotation; a faint planetary-scale Y pattern; a polar vortex; a stationary bow over the main highland. Nearly featureless in visible light, which is correct |
| Titan-like | Haze is a renderer parameter, not a map. Under it: a cluster of convective methane clouds over the summer pole, thin mid-latitude methane streaks, a broad stratiform ethane cloud over the winter pole, and a faint stratospheric condensate hood (HCN and benzene ice) that the generic condensation rule produces but that has little visual effect. In visible light the haze hides most of this, which is correct; the maps exist for renderers that look through it |
| Triton-, Io-, Enceladus-like | No weather. Plume columns from active vents bent into thin downwind trails, deposition streaks on the surface from the same vents, thin condensate wisps near the plume heights |
| Locked ocean world | A convective cloud shield over the substellar point for slow rotators, shifting east and breaking into a banded pattern as the rotation period shortens; a cloud and precipitation belt at the terminator where day and night flows converge; a clear night side |
| Locked desert world | Little cloud. Thin ice cloud near the terminator, nightside frost fog over cold traps, dust lifted by the substellar convection and carried toward the terminator, dust storms along the terminator front |
| Hot dry desert world | Sparse fair-weather cumulus; haboob-type dust storms with lobate fronts in the dust class; a persistent background dust loading as a sidecar aerosol parameter rather than a map; cirrus; what little cloud exists at high latitudes and over any remaining water |
| Cold dry frozen world | A dry atmosphere over a bright surface: thin cirrus, polar stratus along a weak storm track, blowing-snow streaks downwind of ridges; low cloud amount overall, and little visual effect against the ice |

**Limits.** The dynamical model is the second largest piece of numerical work after plate tectonics. The structural realism it buys is the point: cyclones, jets, superrotation and substellar shields emerge from it, and hypothetical worlds extrapolate sensibly because the physics does. What it does not buy is photograph-level realism of individual systems. Mesoscale texture at 16k depends on the recipes, tropical cyclones and the Venus wave are placed, and a map is a snapshot of weather on one day. A fallback exists if the model proves troublesome for a body class: the climatological circulation of section 4.6 with the structured primitives, which gives a weaker but acceptable result.

### 4.8 Volatiles and hydrosphere

- Phase state per species at local temperature and pressure.
- Sea level from inventory: find the level z such that the integrated volume below z equals the liquid inventory (bisection on the hypsometric curve). Connected-ocean assumption by default; priority-flood lakes with a volume budget as an option. This is how land-to-sea ratio follows from inventory and topography rather than from a slider.
- Cold traps: ice accumulates where annual mean temperature is below the frost point at the species' partial pressure. Polar for low obliquity, equatorial for high obliquity, hemispherical for locked bodies. With limited inventory and a locked body, the nightside takes it all.
- Sea ice where the liquid's surface temperature is below its freezing point. The edge sits on that isoline with a marginal zone of broken floes. Shelves form where land ice flows onto frozen or near-freezing water.
- Seasonal caps as an option, exported as a second state.
- Volatile transport between traps and warm regions over time (Pluto type), producing bright ice in cold traps and exposed dark substrate elsewhere.
- Subsurface oceans do not render but set ice-shell thickness, which sets the ice tectonics style.

### 4.9 Surface processes

Each process is on or off and scaled by the derivation. Order within a step: mass wasting, fluvial, glacial, aeolian, sublimation, weathering.

- **Mass wasting.** Slope-dependent diffusion everywhere, rate by temperature and regolith. Landslides where slopes exceed the angle of repose (canyon walls, crater walls, volcano aureoles).
- **Fluvial** (requires a surface liquid). Flow routing on the cubesphere (D8 or D-infinity with cross-face neighbours), depression filling by priority flood (Barnes et al. 2014), stream power incision ∂h/∂t = U - K A^m S^n with m about 0.5 and n = 1 solved implicitly (Braun and Willett 2013, which is O(n) and fast), sediment transport with deposition at slope breaks (fans, deltas), lake filling. K scales with precipitation and rock erodibility. Dendritic valley networks on old wet terrain; catastrophic outflow channels as an optional event type.
- **Glacial.** Ice thickness where mass balance is positive; shallow-ice flow with temperature-dependent softness; erosion proportional to sliding speed; U-shaped broadening and overdeepening of valleys, fjords at coasts, moraines at margins; ice sheets flatten and bury terrain. Margins over land are set by flow against ablation, so they are lobate, follow valleys and carry a debris-darkened zone. Where ice reaches open water it ends in a calving front: a cliff tens to a hundred metres high with icebergs as marginal texture. Where it reaches water below freezing it floats as a shelf, with the grounding line visible in the slope field and a straight seaward edge merging into sea ice. Nitrogen glaciers on Pluto-type bodies use the same code with different constants.
- **Aeolian.** Surface wind stress from the climate wind field and atmospheric density. Threshold friction speed from the Iversen-White relation: low gravity and dense air make sand mobile at low wind (Titan), thin air needs strong wind (Mars). Dust accumulates where stress is below threshold (basin floors, lee sides), is stripped where above (plateaus, wind corridors). Wind streaks downwind of craters with length proportional to diameter. Dune fields where sand supply exists in basins and crater floors: dark basaltic sand on Mars, dark tholin sand on Titan; dune orientation from the wind field feeds roughness and normals. This stage produces most of Mars's albedo pattern and Titan's equatorial dark belt.
- **Sublimation and periglacial.** Pitting of volatile ice (Swiss-cheese and scalloped terrain), penitente and bladed terrain as a texture class at high elevation on cold methane-ice bodies, convective cell texture on thick nitrogen ice plains, sublimation lag darkening.
- **Space weathering** (airless). Maturity accrues on exposed cells, resets on fresh ejecta and new lava. Maturity controls albedo and spectral slope.
- **Chemical weathering** (atmosphere-surface interaction). Oxidation of iron-bearing surfaces at a rate set by oxidant availability (O2 and water vapour), temperature and exposure time, producing nanophase oxide and hematite fractions (red to orange, Mars). Tholin deposition from CH4 and N2 photochemistry at a rate set by UV flux and methane mixing ratio, dark red-brown where not buried by fresh frost (Pluto's low latitudes, Titan's plains and dunes). Sulfur and SO2 frost on Io-type bodies. Salts and evaporites where liquid has retreated (bright) and along ice-shell cracks (brown). Carbonaceous dark lag on cold outer bodies. Hemispheric asymmetry for locked moons of giant planets: leading-side darkening from dust, trailing-side reddening from plasma. Thermal segregation runaway between dark and bright terrain (Iapetus type) as an option.

### 4.10 Materials and colour

Each cell holds fractions of end-members. Each end-member has a reflectance spectrum from a measured source (USGS Spectral Library v7, RELAB, published PDS spectra). Albedo RGB = spectrum mixed linearly (Hapke intimate mixing as a later refinement), integrated against the CIE 1931 colour matching functions under an equal-energy illuminant, converted to linear sRGB (ACEScg as an export option). Illuminant and star colour are deliberately left to the renderer.

Initial end-member set (about 20): mature and fresh anorthositic regolith, mature and fresh mare basalt, fresh basaltic lava, oxidised fine dust, dark basaltic sand, hematite, water ice frost (fine) and coarse ice, CO2 frost, N2 ice, CH4 ice, tholin (Titan type and Pluto type), S8 sulfur, SO2 frost, short-chain red sulfur, white salts and carbonates, hydrated brown salts, carbonaceous dark material, quartz sand, soil, snow, sea ice, vegetation pigment. Users can add end-members with a spectrum or an RGB.

Liquid bodies: colour from the absorption spectrum of the liquid and depth, R(λ) = R_bottom e^(-2 a(λ) z) + R_deep(λ)(1 - e^(-2 a(λ) z)). Water goes turquoise over shelves and deep blue over basins. Liquid methane is near black.

Roughness: per-material baseline (fresh lava rough, dust mantle smooth, ice plains smooth, fresh ejecta rough, water near zero) combined with micro-slope statistics from the detail model. Specular F0 by material. Metallic is zero everywhere by default. A wetness option raises specular and darkens albedo near liquid.

### 4.11 Detail synthesis

Simulation resolution is 2k to 8k equivalent. Export resolution is up to 16k single file, 32k tiled. The gap is filled per tile by detail that is conditioned on the simulated state, never by uniform noise:

- Continue the crater production function below the simulated size limit, deterministic per tile from the seed, respecting each cell's surface age and the atmospheric cut-off.
- Extend drainage: tributary synthesis along existing valleys following the flow direction field.
- Ridge detail in mountains aligned with the macro slope and orogen direction.
- Dune textures in dune fields, oriented by the wind field.
- Flow lobes and channels on young lava units; polygonal and hummocky textures on chaos and grooved units; pitted texture on sublimating ice.
- Fault-scarp fine structure along mapped faults.

Each recipe produces height, composition and roughness detail together. Each boundary type in section 4.12 has an edge treatment: crenulation and icebergs on calving fronts, fine lobes and moraines on flow margins, floe texture on sea-ice margins, tapered streak ends on dust boundaries, aspect-dependent patchiness along snow lines, embayments and islands of old terrain along lava contacts.

Each unit ID has a detail recipe. The recipe set is the place where a trained eye adds realism without touching the physics.

### 4.12 Boundaries and variation

The realism of a planetary map lives in its edges. A noise field has one kind of edge everywhere. A real surface has many, and each is the signature of the process that made it. ExoGen does not choose edge styles. They follow from the pair of processes that meet at the boundary, from the age of the feature, and from what has acted on it since. This section makes that explicit so it can be tested.

| Mechanism | Geometry | Height edge | Albedo edge | Examples |
|---|---|---|---|---|
| Level set of a flat liquid surface | Exactly determined by topography; the fine structure is the topography's | Sharp | Sharp from land to water, then graded within the water by depth | Shorelines, lakes, Titan seas, nitrogen ice filling a basin to a level |
| Threshold of a smooth field | Isoline of the downscaled temperature or frost point; follows ridges and valleys through the lapse rate | None or small | Sharp line, soft coverage near it | Snow line, sea-ice edge, frost line, seasonal cap edge |
| Flow front over land | Lobate; ponds in lows, tongues down valleys, stops where it freezes or where ablation wins | Sharp step at the front | Sharp fresh-on-old contact | Lava margins, maria embayments with islands of old terrain, glacier and ice-sheet margins, landslide aureoles, cryolava |
| Flow meeting water | Calving cliff, or a floating shelf with a straight seaward edge | Very sharp | Very sharp | Ice sheet reaching open ocean, ice shelves |
| Fault-bounded | Linear or arcuate traces, one side up | Sharp step | Weak unless fresh material is exposed | Scarps, graben walls, grooved bands, lineae, rift flanks |
| Impact | Circular rim; radial and lobate ejecta; rampart edge | Sharp rim when fresh | Sharp ejecta edge when fresh; rays feathered and radial | Craters, basins, ramparts |
| Transport deposit | Follows wind stress or fallout distance; thins gradually | None to gentle | Graded over tens to hundreds of km; tapered streaks | Dust, sand sheets, plume fallout, tholin mantles |
| Age or exposure controlled | Continuous except across resurfacing contacts, where it steps | None | Graded with age; stepped at contacts | Space weathering, oxidation, ray fading |
| Resurfacing contact | The outline of the younger unit | Inherited from the unit | Sharp step in every age-dependent property | Plains against highlands, grooved against dark terrain |

Two consequences matter most. A sharp threshold applied to a smooth field gives a sharp edge, so phase changes and water levels turn gentle climate gradients into crisp lines. A diffusive transport process gives a graded edge, so dust and fallout never have crisp lines. Those two rules cover most of what makes a map read as real rather than as noise.

**Age softens edges across all feature types.** The age-ordered replay (section 3.4) applies to scarps, flow fronts, lineae and rims alike. Fresh features are sharp. Old ones are rounded, mantled and cut by whatever came later. Sharpness is a function of mechanism, age and later processes, never a style setting.

**Climate fields are never applied at coarse resolution.** Temperature is downscaled to the surface grid with the lapse rate on the fine elevation plus aspect and slope terms. Precipitation is downscaled with fine-scale orographic lift. Every threshold (freezing, frost point, snow line, sea-ice edge) is evaluated on the downscaled field, so its isoline follows real ridges and valleys rather than the smooth coarse field. Flow and transport processes then reshape the thresholded result on the fine grid.

**Worked example: a tidally locked ocean world.** The substellar region is open ocean. Toward the terminator the sea surface crosses its freezing point and sea ice begins: a sharp edge on the ocean-surface freezing isoline, with a marginal zone of broken floes as texture. On nightside land, snow accumulates in the cold trap and the ice sheet flows under its own weight toward the terminator and down into valleys. Over land its margin is lobate, with tongues in valleys and a debris-darkened zone where ablation exposes entrained regolith. Where the ice sheet reaches a coast with open water it ends in a calving cliff: a sharp step in height and albedo, with icebergs as marginal texture. Where it reaches a frozen sea it floats as a shelf that merges into the sea ice. One cap, four different edges, chosen by the land-sea layout under the terminator rather than by a parameter. Convective cloud sits over the substellar point, and the terminator belt is where liquid, ice, land and cloud all meet.

**Variation within units.** Colour and relief inside a unit come from three sources, in this order of weight: the process fields themselves (dust thickness varies continuously with wind stress, maturity with exposure), sub-resolution process statistics from the detail recipes (small-crater gardening, dune patches, flow lobes), and conditioned noise inside the recipes for grain size, mottling and roughness. Every detail recipe produces height, composition and roughness detail together, so colour variation is tied to relief variation rather than painted over it. Plain noise is never applied to albedo.

### 4.13 Outputs

See the table in section 0. Height is exported in metres relative to a chosen datum. Oblate bodies can export relative to the reference ellipsoid. Normals are computed from the spherical gradient so polar pixels are correct. Sidecar JSON carries parameters, derived values, seeds, legend for unit IDs, datum, radius, and units for every map.

Atmosphere: no maps. The sidecar carries the derived scattering parameters (Rayleigh coefficients from composition and pressure, scale height from section 4.0, and per-aerosol optical depth, scale height, asymmetry, colour and latitude profile) so a renderer can build its own atmosphere and haze. The globe preview's scattering toggle uses precomputed single and multiple scattering (Bruneton and Neyret 2008) from the same parameters. It is a preview, not an export, and it is deferrable.

---

## 5. Signature features per reference body

Each line gives the process chain that produces the feature. These are the acceptance targets for the look of each preset, by mechanism rather than by geography.

| Body | Signature features | Process chain |
|---|---|---|
| Mercury | Caloris basin and antipodal terrain, smooth northern plains, lobate scarps cutting craters, bright rays, dark heavily cratered rest | Giant basin + flood plains + contraction scarps + space weathering |
| Venus | Uniform basaltic plains, scattered un-degraded craters, no craters below a few km, coronae, tessera highlands, long channels | Episodic resurfacing + dense-atmosphere filtering + plume features |
| Earth | Continents with shelves, young sharp ranges and old rounded ones, island arcs, hotspot chains, rift valleys, 30° deserts, equatorial rainforest, ice caps and snow lines, turquoise shelves | Mobile lid + seafloor age-depth + climate + fluvial and glacial erosion + sea-level fill |
| Moon | Nearside maria in basins, farside highlands, saturated craters, rays, ghost craters, wrinkle ridges | Crust dichotomy + basins + flooding + space weathering |
| Mars | Dichotomy, Tharsis rise with giant shields, Valles Marineris, Hellas with bright floor, Syrtis-type dark wind-swept highlands, bright dusty regions, polar layered caps, valley networks on old terrain | Dichotomy + plume loading + graben to canyon + aeolian dust + oxidation + cold traps |
| Io | Hundreds of paterae, no craters, black fresh flows whitening under frost, red plume rings, yellow mid-latitudes, steep isolated mountains | Extreme tidal heating + resurfacing rate + sulfur deposition + block tectonics |
| Europa | Global lineae, cycloids, bands, chaos, a few bright-rayed craters | Ice shell over ocean + tidal stress + salts + low surface age |
| Ganymede | Dark ancient cratered terrain beside bright grooved terrain, palimpsests | Ice shell + partial extensional resurfacing + viscous relaxation |
| Callisto | Uniform dense cratering, dark lag with bright crater rims, multi-ring Valhalla | Old ice surface + sublimation lag + basins |
| Enceladus | Bright snow everywhere, south polar fissures, young terrain there, cratered north | Focused tidal heating + plume fallout |
| Titan | Dark equatorial dune seas, bright highlands, polar methane seas skewed to one pole, few craters, fluvial networks | Methane cycle climate + tholin sand + aeolian + dense-atmosphere filtering |
| Triton | Cantaloupe terrain, south polar nitrogen cap with dark wind-aligned streaks, smooth plains | Diapir texture + cold trap + plume streaks along wind |
| Pluto | Bright nitrogen-ice basin with convective cells, dark equatorial tholin band, water-ice mountains, bladed terrain, nitrogen glaciers | Basin + volatile transport + tholin deposition + glacial flow |
| Ceres | Bright salt spots in crater floors, dark cratered surface | Brine cryovolcanism + old surface |
| Iapetus | Two-tone hemispheres, equatorial ridge | Thermal segregation + (ridge as an explicit feature) |
| Locked ocean world | Substellar open ocean, sea-ice belt at the terminator, nightside ice sheet with lobate land margins and calving coasts, convective cloud at the substellar point | Substellar insolation + sea-level fill + cold trap + glacial flow + sea ice |

Intermediate worlds are the point of the parameter model, but the reference set is how we know the model is right.

---

## 6. Scaling and validation suite

Run by `exogen validate` and in CI. Three layers.

**Scaling law tests** (unit tests against tabulated data): simple-complex transition and depth-diameter curves on the Moon, Mercury, Mars and Ganymede; maximum relief limits versus observed maxima (Olympus Mons 22 km, Everest 9 km, Maxwell 11 km, Io mountains 17 km); seafloor age-depth; Hadley edge latitude for Earth (about 30°) and single-cell for Venus and Titan; equilibrium and surface temperature for the reference set; tidal heating for Io and Europa within a factor of two of accepted values; locking times for the known locked moons below their ages.

**Statistical tests** on generated presets: crater size-frequency distributions on 4.0 Gyr highlands and 3.3 Gyr plains within 30% of the Neukum isochrons from 1 to 300 km; hypsometric curve shape (Earth and Mars bimodal, Venus unimodal); topographic power spectral slope between -2 and -3 across two decades (pure fractal noise fails this by being too uniform); albedo histogram comparison against published maps; latitudinal temperature profile against observed means.

**Cloud tests**: the dynamical core passes the standard shallow-water test cases on the sphere, including flow over a mountain (Williamson et al. 1992). Earth-like global cloud fraction between 0.60 and 0.72 with zonal maxima at the thermal equator and both storm tracks and minima in the subtropics (Rossow and Schiffer 1999); stratocumulus decks off cold subtropical west coasts; cyclone scale and spacing of a few thousand km. Mars-like winter polar hood and afternoon orographic cloud over the large shields. Venus-like superrotation under a complete deck. Titan-like winter polar hood. Locked ocean world substellar cloud shield for slow rotation that shifts east as the period shortens, compared against published circulation regimes (Yang, Cowan and Abbot 2013; Haqq-Misra et al. 2018).

**Boundary character tests**: sampled boundaries in each preset are measured for edge width in height and albedo and compared with the table in section 4.12. Shorelines resolve within one simulation cell. Lava and ice flow fronts step within one to two cells. Dust boundaries grade over at least ten cells. Ice margins over land have a sinuosity well above that of the climate isoline they started from. The locked ocean world shows each ice-edge type where its land-sea layout calls for it.

**Regression tests**: perceptual hashes of each preset's output at 1k so a change in one stage that alters another body's look is caught.

Reference inputs:

| Body | R km | ρ g/cm³ | g m/s² | P surface | T mean K |
|---|---|---|---|---|---|
| Mercury | 2440 | 5.43 | 3.70 | ~0 | 440 |
| Venus | 6052 | 5.24 | 8.87 | 92 bar | 737 |
| Earth | 6371 | 5.51 | 9.81 | 1.013 bar | 288 |
| Moon | 1737 | 3.34 | 1.62 | 0 | 250 |
| Mars | 3390 | 3.93 | 3.71 | 6 mbar | 210 |
| Io | 1822 | 3.53 | 1.80 | ~0 | 110 |
| Europa | 1561 | 3.01 | 1.31 | ~0 | 102 |
| Ganymede | 2634 | 1.94 | 1.43 | ~0 | 110 |
| Callisto | 2410 | 1.83 | 1.24 | ~0 | 134 |
| Titan | 2575 | 1.88 | 1.35 | 1.47 bar | 94 |
| Enceladus | 252 | 1.61 | 0.11 | 0 | 75 |
| Triton | 1353 | 2.06 | 0.78 | 14 µbar | 38 |
| Pluto | 1188 | 1.85 | 0.62 | 10 µbar | 44 |
| Ceres | 470 | 2.16 | 0.28 | 0 | 168 |

Public domain reference maps for side-by-side comparison in the developer harness: LOLA and LRO WAC (Moon), MOLA and Viking/MGS colour (Mars), MESSENGER MDIS (Mercury), Magellan (Venus), Galileo and Voyager (Jovian moons), Cassini ISS and RADAR (Saturnian moons), New Horizons LORRI/MVIC (Pluto), Dawn (Ceres).

---

## 7. UI specification

Single window, dockable panels, no splash, no decorative imagery. Visual reference: a scientific tool's property pane, with monospaced numeric fields and units in a column.

**Left: parameters.** Collapsible groups: Star, Orbit, Body, Interior, Rotation, Atmosphere, Volatiles, Tectonics, Volcanism, Cratering, Climate, Weather and clouds, Surface processes, Materials, Detail, Output. Each row: label, value field, unit, lock glyph, override marker. Derived rows are dimmed. Hover shows the derivation, its inputs and the reference. Multipliers sit at the top of each group. A search field filters rows. Expert mode toggle. Preset selector with "save as preset".

**Centre: preview.** Toggle between 2D equirectangular and globe. Channel selector: albedo, height, slope, roughness, specular, normal, composition layers, surface age, crust age, unit ID, temperature, precipitation, wind, clouds, masks. Lighting controls for the globe: sun direction, exposure, atmosphere scattering on or off, clouds on or off, height exaggeration. Preview resolution 1k with progressive refinement. A readout under the cursor gives latitude, longitude, height, unit, composition. A time scrubber under the preview runs from the start of the simulation to the present and shows the age in Gyr. Scrubbing shows the cached preview-resolution state at that keyframe. A "set as output epoch" action makes the full-resolution render start from it.

**Right: pipeline.** Ordered list of stages with enable toggles, status, elapsed time, cache state, "run to here" and "invalidate from here". Below it the plausibility panel: warnings with severity, explanation and the parameter they refer to.

**Bottom: diagnostics.** Log and progress. In-app plots: hypsometric curve, crater size-frequency against isochrons, latitudinal temperature and precipitation, albedo histogram. These serve both the user and validation.

**Copy rules.** Labels name the quantity and unit: "Surface pressure (bar)", "Crater density multiplier", "Surface age (Gyr)". Tooltips state the basis in one or two sentences with the reference. Warnings state the condition and the consequence: "Liquid N2 requested. Surface temperature 142 K exceeds the N2 boiling point at 1.2 bar (79 K). No liquid will form." No adjectives, no exclamation marks.

---

## 8. Technology stack

**Recommendation: Rust.** Engine crate, CLI crate, GUI crate with egui on wgpu.

- Compute: CPU first via rayon. Hot kernels (flow routing, stream power, crater stamping at high resolution, detail synthesis) get wgpu compute shaders later, which run on Vulkan, DX12, Metal and GL. CPU fallback stays mandatory.
- UI: egui with egui_dock for panels and egui_plot for the diagnostic plots. It is dense, immediate-mode and matches the technical look. Less polished than Qt for complex layout, acceptable for this UI.
- I/O: `image` for PNG/TIFF 16-bit, `exr` for OpenEXR, `serde` with TOML for projects, `toml` and `serde_json` for sidecars.
- Packaging: a single binary per platform. Windows first since that is the development machine; macOS and Linux follow at no extra architectural cost.
- Optional: PyO3 bindings for validation notebooks with numpy, scipy and matplotlib. No production code in Python.

Why not the alternatives:

- C++ with Vulkan and Dear ImGui would work and is the conventional studio-tool stack. It costs more in build system, memory safety and cross-platform packaging, with no capability gain for this project. Acceptable if you prefer it; the architecture is identical.
- Python with Qt would be the fastest to prototype and the hardest to ship and to make fast at 8k. It is the wrong production choice, but it is the right place for the Phase 0 equation checks.

---

## 9. Performance and resolution budget

| Cubesphere face F | Equirect equivalent | Cells | One f32 layer | Use |
|---|---|---|---|---|
| 256 | 1k | 0.4M | 1.6 MB | Coarse loop, preview |
| 512 | 2k | 1.6M | 6 MB | Coarse loop high |
| 1024 | 4k | 6.3M | 25 MB | Default finishing resolution |
| 2048 | 8k | 25M | 100 MB | High finishing resolution, 16 GB machines |
| 4096 | 16k | 100M | 400 MB | Detail synthesis only, tiled, never held whole |

Targets on a 4-core CPU without GPU: 1k preview under 10 s after a parameter change in a late stage; full 4k simulation under 5 minutes; 16k bake under 30 minutes. GPU paths aim for a 5 to 10x reduction on the erosion and detail stages. Implicit stream-power solvers keep erosion linear in cell count. Crater stamping cost is bounded because the production function's steepness makes the total crater area per size decade roughly constant.

Memory: default 4k simulation with 20 layers is about 500 MB. 8k is about 2 GB. Export tiles are 2048 x 2048 and streamed.

Keyframes: 90 steps at F=256 with the key layers at 16-bit precision is about 5 MB each, 450 MB for the full history, held in memory. At F=512 it is about 1.7 GB and goes to disk. The preview finishing pass per keyframe takes a few seconds and is cached.

Weather model: two layers at F=128 with about 60 model days of spin-up runs in under a minute on four cores; F=256 takes a few minutes. Each additional frame of a cloud sequence costs seconds.

---

## 10. Roadmap

Sizes are relative (S, M, L, XL). Each phase ends with a shippable generator for a world class and a validation gate that must pass.

### Phase 0: research spikes and decisions (M)

- Confirm the stack with a spike: cubesphere grid with neighbour tables, resampling to equirectangular, one layer rendered in egui with wgpu.
- Parameter graph specification as a document: every node, unit, derivation and reference.
- Re-check every number in section 4 against its source; record the values in a constants file with citations.
- Crater scaling spike in a notebook: implement the production function, Pi-scaling and morphology laws; plot against lunar and Martian crater databases (Head et al. 2010, Robbins and Hynek 2012).
- Climate spike in a notebook: 2D energy balance with Held-Hou circulation; check Earth, Mars, Venus and a locked case.
- Architecture decision records for grid, time model, caching, determinism.

Gate: spikes reproduce the published curves; decisions recorded.

### Phase 1: foundation (L)

- Core crate: parameter graph with derive, set, multipliers, plausibility rules and tests.
- Grid, layered fields, tiling, caching, hashing.
- Stage runner with the coarse-loop and finishing-pass structure, initially with a single placeholder terrain stage (spherical harmonics plus noise).
- Export: height EXR and PNG, albedo PNG, sidecar.
- CLI: render, derive, validate (empty), preset list.
- App: parameter panel, 2D preview, globe without atmosphere, pipeline panel, log.
- Keyframe storage for the coarse loop and a time scrubber over it.
- Determinism test across runs.

Gate: a project file round-trips; the same seed produces identical output twice; derived values for the reference bodies match the Phase 0 tables; scrubbing the timeline shows the placeholder terrain changing step by step.

### Phase 2: airless cratered worlds (L)

- Production functions, chronology, flux history.
- Impactor-to-crater scaling, morphology classes, ejecta, rays, secondaries, oblique forms.
- Time-integrated emplacement in the loop and statistical emplacement in the finishing pass with pre-degradation.
- Surface age field, diffusion, infill, space weathering, regolith roughness.
- Primordial crust with dichotomy and isostasy.
- Materials: regolith end-members, maturity, colour integration.
- Preview-resolution finishing pass per keyframe so the scrubber shows the crater population at any age. Output epoch parameter.
- Presets: Moon, Mercury without scarps, Callisto without relaxation.

Gate: statistical tests of section 6 for crater SFD, depth-diameter and power spectral slope pass on the Moon preset; a 1,500 km and a 6,000 km body produced from the same inputs differ in crater morphology as the scaling laws predict.

### Phase 3: interior, volcanism, stagnant and episodic lid (L)

- Heat budget through time, lithosphere, relief limit, magma budget.
- Contraction scarps. Plume rises with loading stress, radial graben, canyon growth, concentric ridges.
- Shields, calderas, stratovolcanoes, flood basalts by topographic fill, ghost craters, wrinkle ridges, channels.
- Episodic resurfacing with retained tessera, coronae, domes.
- Io regime with sulfur and SO2 materials and block mountains.
- Eccentricity history curve driving tidal heating through the loop.
- Presets: Mercury complete, Moon with maria, Venus, Io, Mars without atmosphere.

Gate: Mars preset produces a rise, a canyon system and giant shields consistent with the relief limit; Venus preset crater population is sparse and uniform; Io preset has no craters and hundreds of paterae.

### Phase 4: ice worlds (L)

- Composition from density, ice shell thickness, relaxation viscosity.
- Tidal stress field, lineae, cycloids, bands, chaos, grooved terrain, focused polar heating.
- Viscous relaxation of craters, palimpsests, sublimation lag.
- Cryovolcanism, plume fallout, salt deposits.
- Hemispheric asymmetry and thermal segregation for locked moons.
- Ice and salt end-members.
- Presets: Europa, Ganymede, Callisto complete, Enceladus, Ceres, Pluto and Triton without a volatile cycle.

Gate: Ganymede preset shows two terrains with a crater density contrast; Europa preset crack orientations follow the tidal stress field; relaxation removes large craters first.

### Phase 5: atmosphere, climate, volatiles and erosion (XL)

- Insolation for arbitrary rotation states including locked. 2D energy balance, circulation, wind field.
- Moisture transport, precipitation, orographic effects.
- Volatile phase states, sea-level fill, cold traps, hemispheric caps, volatile transport.
- Fluvial erosion with implicit stream power, lakes, deltas. Glacial flow and erosion. Aeolian dust transport, streaks, dune fields. Mass wasting. Sublimation textures.
- Chemical weathering: oxidation, tholins. Atmospheric crater filtering.
- Weather model: two-layer moist shallow water on the coarse cubesphere with topography and zonal or day-night forcing, validated on the standard shallow-water test cases.
- Cloud type classification, structured primitives for tropical cyclones, hoods, plumes and dust storms, per-class export.
- History curves for atmospheric pressure, volatile inventory and stellar luminosity.
- Presets: Mars complete, Titan, Pluto complete, Earth-like with simple continents from the primordial crust stage, a tidally locked eyeball world, and an early-wet Mars that dries out along its pressure history.

Gate: Earth-like at 1 AU around a G2V star gives global mean temperature within 5 K of Earth and a Hadley edge within 25 to 35°; deserts and rain shadows appear in the right places; Mars preset dust albedo correlates with topography and wind as on Mars; Titan preset has dark equatorial dunes and polar seas; the locked ocean world has a sharp sea-ice edge on the freezing isoline, lobate ice margins over land and calving fronts at open-water coasts. Earth-like cloud fraction and its zonal profile match the climatology within tolerance; the Mars preset shows a winter polar hood and afternoon orographic clouds over the large shields; the locked ocean world shows a substellar cloud shield that shifts east as the rotation period shortens.

### Phase 6: mobile-lid plate tectonics (XL)

- Kinematic plates, boundary classification, subduction, collision, rifting, transforms, hotspots, re-organisation events.
- Seafloor age-depth, passive margins and shelves, trenches, arcs, orogen ages.
- Integration with erosion so old and young belts differ.
- Preset: Earth-like with convincing continents.

Gate: hypsometry is bimodal with peaks near +0.5 km and -4.5 km; mountain belts sit on convergent boundaries; shelves exist on passive margins; seamount chains align with plate motion.

### Phase 7: detail synthesis and studio export (L)

- Tiled detail pass with per-unit recipes.
- Roughness, specular, normal, ambient occlusion, masks, cloud height, unit ID legend.
- 16k single-file and 32k tiled export, cube maps, colour management options.
- Sequence export of epochs from the CLI and the app.
- Cloud detail recipes at export resolution, per-class export with heights, cloud type ID map, feature list, optional composite, cloud shadow bake, cloud sequences.
- Atmosphere scattering preview on the globe, driven by the derived atmosphere parameters in the sidecar. Optional and deferrable; nothing in the export depends on it.
- GPU kernels for the hot stages.

Gate: a 16k export of each preset completes within the time budget and shows no tile seams or pole artefacts; maps load correctly in Blender and a game engine with the documented shading setup.

### Phase 8: vegetation and biomes (M)

- Köppen-type classification, coverage, pigment parameter, seasonal snow, biome ID map.

Gate: Earth-like preset biome map matches Earth's zonal distribution of forest, grassland, desert and tundra.

### Phase 9: polish and release (M)

- Preset library, documentation with references per parameter, UI pass, test coverage, builds for Windows, macOS and Linux, user guide.

### Risks and mitigations

| Risk | Mitigation |
|---|---|
| Scope | Each phase ships a usable generator. Stage interfaces are fixed in Phase 1 so phases are independent. |
| Model tuning takes longer than coding | Phase 0 notebooks settle the equations before Rust code is written. Validation gates are explicit. |
| High-resolution performance | Coarse loop plus conditioned detail; implicit erosion solver; GPU kernels in Phase 7 not before. |
| Parameter graph state management | The graph lives in the core crate with tests. The UI is a view of it. |
| GPU non-determinism | CPU path is the reference. The sidecar records the path. |
| Colour realism | Measured spectra only; curated set first; user-extensible. |
| Plate tectonics convincingness | Phase 6 is late and isolated. Earth-likes are available earlier through the primordial crust stage. |

---

## 11. Decisions

All resolved on 2026-10-10. Architecture decision records are in `docs/adr/`.

1. **Stack.** Rust with wgpu and egui. Dependencies must be GPL-3.0-compatible; the planned crates are MIT or Apache-2.0 licensed, which is compatible.
2. **Simulation grid.** Cubesphere at every resolution; a Fibonacci point set only for plates (ADR 0001).
3. **Colour management.** Linear sRGB EXR plus sRGB PNG by default, ACEScg optional.
4. **Minimum body radius for the spheroid assumption.** 400 km nominal, warning below.
5. **Platform order.** Windows first, then macOS and Linux.
6. **Licence.** GPL-3.0-or-later.
7. **Name.** ExoGen. Repository: github.com/Brimfyr/ExoGen.

---

## 12. References

Starting points, grouped by stage. Textbooks first because they cover most of the quantitative relations in one place.

**General.** Melosh, Planetary Surface Processes (2011). Melosh, Impact Cratering: A Geologic Process (1989). Turcotte and Schubert, Geodynamics (3rd ed.). Pierrehumbert, Principles of Planetary Climate (2010). Catling and Kasting, Atmospheric Evolution on Inhabited and Lifeless Worlds (2017).

**Body derivation.** Eker et al. 2018 (main-sequence mass–luminosity and mass–radius relations). Zeng, Sasselov and Jacobsen 2016 (mass-radius relations). Fortney, Marley and Barnes 2007. Peale, Cassen and Reynolds 1979 (tidal heating). Gladman et al. 1996 (locking time). Williams and Kasting 1997 (energy balance scaling).

**Cratering.** Neukum, Ivanov and Hartmann 2001 (production function and chronology). Ivanov 2001 (Mars/Moon ratio). Zahnle, Schenk, Levison and Dones 2003 (outer solar system rates). Holsapple 1993 (scaling). Collins, Melosh and Marcus 2005 (impact effects equations). Pike 1977, 1980, 1988 (morphometry: Moon, Mercury). Schenk 2002 (icy satellite morphology transitions). McGetchin, Settle and Head 1973 (ejecta thickness). Fassett and Thomson 2014 (diffusional degradation). Head et al. 2010 (lunar crater database). Robbins and Hynek 2012 (Mars crater database). Melosh and McKinnon 1978 (multi-ring basins).

**Tectonics and volcanism.** Byrne et al. 2014 (Mercury contraction). Watters et al. (lobate scarps). Banerdt, Golombek and Tanaka 1992 (Tharsis stress and faulting). Andrews-Hanna, Zuber and Banerdt 2008 (Borealis basin). Wieczorek and Phillips 1998 (lunar crust and maria). Parsons and Sclater 1977 and Stein and Stein 1992 (seafloor age-depth). Cortial, Peytavie, Galin and Guérin 2019 (Procedural Tectonic Planets). Viitanen 2012 (PlaTec). Greenberg et al. 1998 (Europa tidal cracking). Pappalardo et al. 1998 (Ganymede grooved terrain). Schenk and Bulmer 1998 (Io mountains). Spencer et al. 2006 (Enceladus south polar terrain).

**Climate and surface processes.** North 1975 (energy balance model). Held and Hou 1980 (Hadley circulation). Whipple and Tucker 1999 (stream power). Braun and Willett 2013 (implicit O(n) erosion solver). Barnes, Lehman and Mulla 2014 (priority flood). Cordonnier et al. 2016 (large-scale terrain with uplift and erosion). Iversen and White 1982 (threshold friction speed). Lorenz et al. 2006 (Titan dunes). Hayes 2016 (Titan lakes). Bertrand and Forget 2016 (Pluto volatile cycle). Earle et al. 2017 (Pluto volatile transport). Hynek, Beach and Hoke 2010 (Mars valley networks).

**Weather and clouds.** Vallis, Atmospheric and Oceanic Fluid Dynamics (2nd ed.). Phillips 1954 (two-layer model). Williamson et al. 1992 (shallow-water test cases on the sphere). Ronchi, Iacono and Paolucci 1996 (cubed-sphere shallow water). Showman and Polvani 2011 (superrotation under day-night forcing in shallow water). Penn and Vallis 2017 (shallow-water model of locked and unlocked terrestrial planets). Haqq-Misra et al. 2018 (circulation regimes of locked planets). Yang, Cowan and Abbot 2013 (substellar cloud shield). Rossow and Schiffer 1999 (ISCCP cloud climatology). Clancy et al. 1996 (Mars aphelion cloud belt). Fukuhara et al. 2017 (Venus stationary bow wave). Rodriguez et al. 2009 (Titan cloud distribution). de Kok et al. 2014 (Titan south polar HCN ice cloud). Vinatier et al. 2018 (Titan benzene ice clouds). Soderblom et al. 1990 (Triton plumes). Abbot, Voigt and Koll 2012 (clouds on a snowball planet). Lobo et al. 2023 (locked desert planets).

**Materials.** Kokaly et al. 2017 (USGS Spectral Library v7). RELAB spectral database. Hapke 2012 (Theory of Reflectance and Emittance Spectroscopy). Pope and Fry 1997 (water absorption). Bruneton and Neyret 2008 (precomputed atmospheric scattering).

**Prior art to study.** tectonics.js, PlaTec and WorldEngine, Songs of the Eons, FastScape, SpaceEngine, Gaea and World Machine (for export conventions studios expect).
