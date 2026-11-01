# Phase 5 report: atmosphere, climate, volatiles and erosion

Status: started 2026-10-31, after Phase 4F's close; 5A, climate on every body, passed the user's review on 2026-10-31 in five milestones; 5B, the volatiles, under way: the frost's seasons (F239), the cold traps (F240), the caps of the inventory's ices (F241) and the frost's own gas's (F242), and the poles' shadows (F243).

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

### 5A passed

The user's verdict on the review page (2026-10-31): "Everything here looks good." The gate as reviewed: the Earth-like's mean 284.8 K, its Hadley edge 32.6°, Mars's mean 213.6 K, the Moon's equator by day and night, Mercury's hot and warm longitudes met; Venus's mean (753 K, 2.2 percent high) and Titan's on its relief (88.9 K against 93.7, its uplands far above its methane seas) near; the poles' permanent shadows left to 5B's cold traps.

## 5B, the volatiles

The order: the frost's seasons first, the atmosphere's main volatile condensing on the cold ground and holding it at its frost point, then the caps the inventory lays in the cold traps with their materials and their albedo, the shadows of the airless poles, sea ice and the locked ocean world, one temperature for every stage, and the climate through the loop's history with the atmosphere's escape and the volatiles' transport.

## The sixth milestone: the frost's seasons

5B's first milestone (the order above): the atmosphere's main condensable gas condensing on the cold ground, holding it at its frost point, and the pressure falling and rising with the caps.

- **F239.** Mars's carbon dioxide holds its winter poles at its frost point, 145 K, where they fell to 80 K; its seasonal caps reach some 50° of latitude in their winters, the northern 1,099 kg/m² at the pole by the end of its winter, about MOLA's metre; its pressure swings from 502 to 696 Pa, 28 percent, as Viking's did. Pluto's and Triton's nitrogen lays thin seasonal caps from their own gas. The frost is drawn over the ground at the snapshot's season, its slopes taking the frost of the places they face, and exported as `frost.exr`.

How it went:

- **The surface's balance.** The first frost took the step's whole shortfall, the column's conduction counted twice, and laid some 10,000 kg/m² at Mars's north pole; the column is now solved under a surface held at the point, and the surface layer's own balance is the frost it takes.
- **One atmosphere.** With every latitude drawing on a gas of its own, held at the body's pressure, Pluto and Triton froze over and their caps grew without end; the latitudes now run in step, the frost the cold ground takes the gas's own.
- **A step of 2.6 years.** Pluto's year is 248 of Earth's, its 96 steps each longer than its air takes to freeze out, and every latitude took the same air in the same step: its frost came to three times its gas. The frost point is now settled each step, at the one point at which the frost taken leaves the gas its pressure.
- **A pressure no gas can reach.** Raising the gas's mass until the year's mean pressure met the body's ran away where the frost's own balance holds the pressure: Pluto's frost reached hundreds of kg/m² and its air still stood under 0.4 Pa of its 1. The gas is now the frost's that comes and goes in the share of the body's pressure the air keeps, all of it on Mars and little on Pluto. Their pressures stand far under the measured, the frost radiating as a black body where nitrogen ice's emissivity is some 0.6: a shortfall reported, not tuned.
- **A locked Mars.** Mars locked to its star at 3 AU laid frost at 0 K on its substellar ground: the steady transport, held through the cycle, took more heat from it at aphelion than the frosted ground absorbed, and frost formed without the gas to form it. Frosted ground now takes the heat the air brings it but gives none, as ground its frost holds below its steady temperature would draw heat rather than shed it; the locked Mars lays its gas on its night side and poles. The diagnostic found the transport's local response, by Williams and Kasting's scaling with rotation, 30 times Pluto's radiation's slope.
- **The renders.** Mars, Pluto and Triton rendered at 2,048 px take their frost at the snapshot's season, the northern spring equinox: Mars's northern cap covers the ground to 50° north, frost on its pole-facing slopes past the edge, and its southern begins at 81° south; Pluto's northern cap reaches 63°, Triton's 61° north and 80° south. The caps' outlines follow their latitudes, broken only by the slopes: the cycle's ground is a latitude's, its albedo and inertia the same along it. The Moon is byte-identical to 5A's renders in every map; the Earth-like's and Titan's albedo and relief are identical and their temperatures within 0.04 K, their latitudes now running in step under gases whose frost points they never reach.

## The seventh milestone: the cold traps

- **F240.** Where each ice of the inventory would last: each runs its own column at every cell of the climate's grid, at its height's pressure and beside the bare ground about it. Mars's water goes to a pole, the north on ground at the datum and a deep southern basin on the model's own Mars; Earth's water to its poles; Pluto's methane to 20° north. A basin traps better than a height. Nothing is drawn yet: the caps that fill the traps come next.

How it went:

- **The Earth-like's equator.** Its water first found its coldest trap at the equator: the year's steady transport, held in the column, took 41 W/m² from the bright ice there as from its bare ground, whatever the ice's temperature. The ice now exchanges heat with the air over the bare ground beside it (F237's), and under the Earth-like's thick air it stands near the air's temperature; the water goes to its poles.
- **Six thousand columns.** Every cell's own column took 8 to 42 s a body, Triton's four tilts the longest; a fast spinner's climate is the same along a latitude, and its traps are run on its latitudes at seven heights, 0.1 to 0.6 s.
- **No residual carbon dioxide cap.** At the seasonal frost's albedo Mars's carbon dioxide ice lasts only in its deepest polar basins; its residual south cap, at 0.8, waits on the ice's own albedo.
- **Mars's water vapour.** The clean ice's albedo holds Mars's air at 2.4e-6 Pa of water vapour, where it holds some 0.1: the residual cap's ice is dusty.
- **The renders.** Mars, the Earth-like, Pluto, Triton and the Moon rendered at 2,048 px are byte-identical to the last milestone's in every map; the traps reach the sidecar only.

## The eighth milestone: the caps of the inventory's ices

- **F241.** Each ice of the inventory but the frost's own gas lays caps in its cold traps, from the coldest, its surface Nye's plastic profile from its margin, until they hold its inventory; a working liquid's only where its ice never melts, the rest left to the sea. Mars's water stands in caps at both poles, up to 4.4 km thick over 2.4 percent of the body, its polar basins filled; Pluto's and Triton's methane in domes near their equators.

How it went:

- **Each cell's own trap.** Read at each fine cell's own height, the traps let the relief of the craters decide: within a few degrees of a pole the latitude hardly changes the trap, and the ice covered the cells a crater's floor made colder one by one, a scatter of thin ice. A cap kilometres thick buries the relief under it: each place's trap is now its region's on the climate's grid, read smoothly between its cells, and the caps stand as domes.
- **Pluto's methane over the whole body.** Laid over every place, ice holds exactly its inventory, and a rounding decided whether the fewest places were sought at all: Pluto's methane lay 4 m thick over the whole body. Ice that may cover every place now spreads over all of them only if no fewer hold it.
- **No ice on the Earth-like.** Its summer poles stand at some 295 K, its seasons land's until its ocean holds them (5A), and its water's ice melts everywhere: it lays no ice sheets, and its sea keeps all its water.
- **The renders.** Mars, Pluto and Triton rendered at 2,048 px take their caps; the Earth-like, Titan and the Moon are byte-identical to the last renders in every map.

## The ninth milestone: the frost's own gas's caps

- **F242.** The frost's own gas lays its caps last, where its ice grows the most: Mars's carbon dioxide in a basin at 77 to 79° north, up to 0.9 km thick over its water's ice; Pluto's nitrogen mostly between 37 and 43° south, Triton's at both poles. Mars's preset holds its carbon dioxide as measured, 0.12 m of its liquid where it held 8 m of ice.

How it went:

- **Mars's carbon dioxide.** The preset held 8 m of carbon dioxide ice over the whole body: laid, it would have stood a hundred times Mars's residual cap and buried deposits. It holds 0.12 m now (Phillips et al. 2011).
- **North where Mars's is south.** The model's Mars's south polar band lies lower than its north, where Mars's stands 3 to 4 km higher, and its basins' air raises the frost's point: its carbon dioxide lies in a northern basin, where Mars's residual cap stands at its south pole at an albedo of some 0.8. 5B's gate wants it in the south.
- **Pluto's nitrogen.** At the pressure the frost holds, Pluto's nitrogen ice loses mass everywhere, and its caps lie where it loses the least, at 40° south; Pluto's gathered in Sputnik Planitia over its history, which comes with the volatiles' transport.
- **The air's reservoir.** The caps of the frost's own gas do not yet hold the air: their exchange with it waits on the volatiles' transport, where the reservoir's history matters.
- **The renders.** Mars, Pluto and Triton rendered at 2,048 px take their caps; the Earth-like and the Moon are byte-identical to the last renders in every map.

## The tenth milestone: the poles' shadows

- **F243.** Ground near an airless pole whose horizon, in every direction, stands above the highest the sun's disc climbs never sees it, and stands at what its walls and the body's heat hold it at: a crater at the Moon's south pole keeps its floor at 41 K; the Moon's polar craters keep shadowed crescents under their equatorward walls, 39 K at the coldest, Mercury's 45 K, Ceres's 29 K. 5A's gate item met only in part, the poles' permanent shadows, is met.

How it went:

- **The sun at a pole.** Near a pole the sun's azimuth is its hour angle, and at the solstice it stands highest toward the equator, at the colatitude and the tilt; a place is shadowed where its horizon in each direction stands over the sun's highest there, the disc's top included, and the equator's direction, where the sun stands highest, is tried first.
- **Too many shadows on the Moon.** Its shadowed share comes to some three times the measured: near the pole, where the sun climbs under 3°, a cell a few hundred metres lower than its neighbours is shadowed, and the finish's cells, some 5 km, do not resolve the smallest craters that hold most of the Moon's shadow.
- **The renders.** The Moon, Mercury, Ceres and Callisto rendered at 2,048 px take their shadows in their temperatures alone; Mars, tilted 25°, is byte-identical.
