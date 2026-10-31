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

## The third milestone: the day and the year

The energy balance gave every body the radiative balance with its year's mean sunlight, which an airless body never stands at: the Moon's equator came to 284 K, where Diviner finds noons near 390 K and nights near 95.

- **F236.** Every place runs a column of ground through the star's path, its regolith's thermal inertia from its grains in vacuum and the gas in its pores, the air's heat on its surface, the transport the year's balance found held steady: through the year until its seasons repeat and then through a day at each season, or through the orbits a resonant body repeats in. The Moon's equator now runs from 98 to 381 K, Mercury's hot longitudes reach 696 K, Europa's equator 79 to 124 K, all near the measured ranges.

How it went:

- **One column, two clocks.** A fast spinner's days are thousands to its year (Miranda's 21,700): its year runs on the day's mean sunlight and its days at 24 seasons from the ground the season left; a resonant body's day is its year, and it runs whole. The deep ground is settled at the surface's mean after each year, as it would settle after many.
- **The thermal inertia from the gas.** A regolith conducts through its grains' contacts and the gas in its pores (Presley and Christensen 1997): the Moon's 55 in vacuum, Mars's about 185 at 6 mbar, a bar's dry sand about 470, one process from the airless to the thick.
- **The Earth-like's seasons.** Without its ocean the Earth-like's seasons swing as land's do: its summer poles reach 295 K and its mean falls 2 K under its steady balance, still within the gate. Its ocean and sea ice come with 5B.
- **The renders.** The Moon, Mercury, Mars, Titan and Europa rendered at 1,024 px are byte-identical to the last renders in albedo, height, units and surface age; their temperature maps now carry the days and seasons, Mercury's with its two hot longitudes, and the finishing pass takes some 3 s more.

## The fourth milestone: the climate brought down to the ground

The climate stood on its coarse grid, a summit as warm as the plain under it.

- **F237.** Each place takes the climate's temperature by its own height and slope: the ground's radiated share follows its own sunlight and greenhouse, the rest the air over it. A slope takes the sunlight of the latitude it faces; a height cools by the air's share of the lapse rate and the radiation's share of the thinner greenhouse. An Earth-like's 3 km summit stands some 15 K under its sea level, Venus's 11 km summit near 650 K, and a slope on the airless Moon at the temperature of the latitude it faces.

How it went:

- **The greenhouse alone was too weak.** Thinning the grey greenhouse with height cooled an Earth-like's summit by only 2.3 K/km: water vapour's optical depth saturates (F235) and thins slowly. The cooling a mountain feels comes from the convecting air, which the ground follows by its exchange with it; the radiation's share of the exchange weighs the two.
- **One share for slopes and heights.** The same share, radiation against the air's exchange, decides how far a slope's own sunlight moves it: all the way on the airless Moon, almost not at all under Titan's air.
- **The renders.** Earth, Venus, the Moon, Mercury and Mars rendered at 1,024 px are byte-identical to the last renders in albedo, height, units and surface age; their temperature maps now carry the relief, an Earth-like's continents cooler than its sea, Venus's highlands near 655 K and its rift floors warm, Mars's basins and sunward walls a few kelvin warmer. Noticed: Venus's surface at the tests' preview quality holds one cell whose height is not finite (none at standard quality); the downscaling takes it at the datum, and its source is to be traced.

## The fifth milestone: the spin's history from its present state

The user's direction (2026-10-31): the present rotation is measured and should be the input, its history derived from it, how long a body took to lock and whether it could have; satellites included where they matter, a binary pair above all; anything that looks wrong flagged rather than forbidden.

- **F238.** The present period is the input and `rotation.mode` describes it; the history comes from the tides of the primary and the largest satellite, the satellite's orbit taking up what the spin loses. Earth's day runs back to 22.4 h 620 Myr ago and 4.9 h at its birth, the Moon at its Roche limit near Earth's age; Mercury reached its 3:2 0.81 Gyr after it formed; Pluto has been locked to Charon since its first Myr. No preset pins its rotation any more, and the plausibility rules flag what the history cannot reach.

How it went:

- **The lunar timescale problem.** Run back on Earth's present tidal Q, 12, the Moon reaches Earth 1.5 Gyr ago; Earth's oceans dissipate today some three times their long average (Bills and Ray 1999), and on 34 the Moon stands at its Roche limit near Earth's age, a calibration on Earth.
- **A lock run backward.** Pluto's spin, run back from today against Charon's tide, ran away to a 2 h day: going back in time a lock repels what it attracts. A spin locked to its satellite is held, like a spin locked to its primary, back to its capture.
- **Venus.** Its slow retrograde spin runs back to a 30 h retrograde one, the Sun's tide on it too weak to have locked it in its age; the thermal tides that hold it (Correia and Laskar 2001) want the atmosphere's heating profile, which the weather model will give (5G).
- **The plan.** 5A's last bullet, the spin from tides, became this; 'one temperature for every stage' moved to 5B, where the ice's sublimation is integrated through the days and seasons rather than taken at a mean.
- **The renders.** Mercury, the Moon, Pluto, Earth and Venus rendered at 1,024 px are byte-identical to the last renders in every map, the temperature's included: nothing reads the spin's history yet but the plausibility rules.

## 5A's close

The climate through the loop's history, 5A's last item, moved to 5B (the user, 2026-10-31): with the atmosphere held at today's state until 5B it would add only the faint young Sun, shifting every preset's thermal history a little, where beside the atmosphere's escape and outgassing it carries Mars's thick early atmosphere and an M dwarf's planets losing their water in its bright youth. The user asked whether the history would need every kind of star: the main sequence from 0.08 to 3 solar masses is modelled (F3), a low-mass star's bright contraction before it is added with 5B, and hot stars and stellar remnants stay outside the plan unless chosen. 5A's gate review follows.
