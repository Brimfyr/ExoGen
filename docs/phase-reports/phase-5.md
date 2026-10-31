# Phase 5 report: atmosphere, climate, volatiles and erosion

Status: started 2026-10-31, after Phase 4F's close; sub-phase 5A, climate on every body, under way.

Phase 5 gives the bodies their air, their weather's long average, their volatiles and the processes that wear a surface down: sunlight and the temperature it makes on every body, seas, caps and cold traps from the inventory, rivers, glaciers, wind, dust, sand and the chemistry that colours them. The plan (PLAN.md, Phase 5) runs it as seven sub-phases in order, each with its own gate and review: 5A climate on every body, 5B volatiles, 5C water and ice at work, 5D dust, sand and chemistry, 5E surface physics on every body, 5F Io and Thaumasia, 5G weather and clouds. Mars, Titan and Pluto are its test bodies for the look; an Earth-like and a tidally locked ocean world give the climate its numbers.

The order (2026-10-31), the user taking the recommendations: climate first, since every later stage reads temperature and wind, rather than building only what Mars needs and widening later; the weather model as the last sub-phase with its own gate, since nothing on the surface reads it. The materials grow by what processes make, olivine and evaporites added to the plan's list beside quartz sand (PLAN.md 4.10).

The findings continue Phase 4F's numbering from F234.

## The first milestone: the star's light on any body

5A begins with the sunlight itself, since every temperature after it is read from it. The ice stage knew a latitude's mean as S cos(lat) / pi, the equinox's day with no tilt, and nothing else knew more.

- **F234.** The year's sunlight at any point from one integral of the star's path across the sky: the day's mean through the seasons for a body that spins many times an orbit, the time integral over the repeat for one in an exact spin-orbit resonance, a moon's year and seasons its planet's. Checked against a brute-force integration: the Earth's equator 415.9 W/m² and poles 172.4, Mercury's hot longitudes 2.53 times its warm ones, a locked body's substellar hemisphere, the poles passing the equator above 54° of tilt (Pluto, Miranda).

How it went:

- **A moon's tilt.** A moon locked to its planet had no tilt (`rotation.effective_obliquity` was zero for every locked body), so Titan had no seasons; its tilt to the Sun is its planet's where its orbit lies in its planet's equator. `rotation.obliquity` is now that tilt for every body, and the presets set it: the Galilean moons 3.13°, Saturn's 26.73°, Miranda 97.77°, Charon 119.6°, the Moon its own 1.54°. Triton's orbit, retrograde at 157° to Neptune's equator, precesses in some seven centuries, so its tilt to the Sun swings between 5° and 51° and its year's sunlight is the mean over the swing.
- **Pins noticed.** Mercury's and Venus's presets pin `rotation.mode` free at their observed periods (their spins from tides are not modelled, principle 11); the insolation reads Mercury's 3:2 from the pinned period. Their spins from tides, a resonance's capture and Venus's atmospheric tides, go into 5A's list.
- **The renders.** Nothing reads the field yet: Mercury, Titan, Charon, Mars and Triton rendered at 1,024 px are byte-identical to milestone 140's in albedo, height, units and surface age.

## The second milestone: the energy balance

With the sunlight known, every cell's temperature follows from what it absorbs, what it radiates and what its atmosphere carries away.

- **F235.** The year's mean surface temperature from the energy balance of every cell of a coarse cubesphere (North 1975): the sunlight it absorbs and the body's heat against the grey atmosphere's outgoing radiation and the atmosphere's transport. An Earth-like at 1 AU stands at 287.8 K, 301 in its tropics and 257 at its poles; Mars 215.7, Venus 755, Titan 91.1. The temperature is exported as a map.

How it went:

- **The tropics ran away.** A first balance on the zero-dimensional model's radiative scheme found no steady Earth: its water vapour's optical depth, the square root of the column, made the outgoing radiation rise by only 0.31 W/m² per kelvin where Earth's zonal means show 2.1, and the tropics warmed without end. Water's bands saturate as they fill; a logarithm keeps Earth's zero-dimensional greenhouse and gives the measured slope (prototypes in Python before the engine's code, scratchpad ebm_proto*.py).
- **The poles too warm.** With a uniform albedo the Earth-like's poles stood at 270 K, 31 K under its equator where Earth's are some 45 K: the low sun and the clouds brighten the high latitudes (North and Coakley 1979). The albedo now rises as the sun's year-mean height falls, under an atmosphere by its column, and the poles stand at 257 K.
- **A floor that fought the balance.** Io's poles would not converge: a floor at the temperature its own heat flow radiates pushed back against the balance. The body's heat now enters the balance itself, and a cell the star never warms settles at that temperature by the physics.
- **Known limitations.** An airless body's mean is too warm until its days and nights (the next milestone); Titan's and Venus's contrasts wait on the wind field (5C).
- **The renders.** The Moon, Venus and Earth rendered at 1,024 px are byte-identical to milestone 140's in albedo, height, units and surface age, with the temperature map added. Mars moves as its zero-dimensional mean warms by 0.8 K (its albedo by more than 0.2 in log luminance on 0.1 percent of the map, its tallest shield from 25.4 to 26.1 km), and Titan at the level of single pixels, its water vapour's optical depth at 91 K changed in the ninth decimal.
