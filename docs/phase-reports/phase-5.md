# Phase 5 report: atmosphere, climate, volatiles and erosion

Status: started 2026-10-31, after Phase 4F's close; 5A, climate on every body, passed the user's review on 2026-10-31 in five milestones; 5B, the volatiles, under way: the frost's seasons (F239), the cold traps (F240), the caps of the inventory's ices (F241) and the frost's own gas's (F242), the poles' shadows (F243), the sea's heat and its ice (F244), the seasons' air (F245), which revises 5A's cycle, the ices that come and go (F246), nitrogen frost's emissivity (F247) and the ground without its seasons' ice (F248).

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

## The eleventh milestone: the sea's heat and its ice

- **F244.** A sea of water runs its own year, its 50 m mixed layer holding its seasons and freezing ice that grows and melts by Semtner's thermodynamics; the air over a place mixes its land's and its sea's. The Earth-like's sea at 45° swings 8 K where its land swings 80; its sea ice stands at the end of its northern winter from some 60° to the pole. The ice is exported and reported, not drawn: the climate's ice line near 60° leaves a belt of it there.

How it went:

- **Ice that never melted.** Its top held at the freezing point, as the frost's is, the ice grew some 10 m a winter and never melted; as Semtner's zero layers, its top cooled under its snow, and the air's heat and its melting ice's albedo let the summer take it back.
- **The land's air.** The ice's air at first was its land's, the continent's winter: the sea at 45° froze. The air over a place now mixes its land's and its sea's by their shares, and a sea that is most of its latitude keeps its own.
- **The climate's transport.** Laid on the ice's top, the transport radiated away and left the ice 10 to 13 m thick at 55° south, where the Southern Ocean stands open; through the sea's water, as a slab ocean's heat transport comes, it melts the ice's base. Where the climate's ice line stands, near 60°, the ice still keeps most of the year.
- **A belt of ice near 60°.** Drawn, the Earth-like's sea ice covered 14.9 percent of the body in March, four times Earth's: a belt near 60° in both hemispheres, around an open southern polar sea melted by its polar day, and the map's mean fell 2.1 K under the 5A gate's. The climate's year's balance there stands at 268 K, its ice line (F235); the sea's year follows it. A restoring toward that balance would only hold it frozen: the ice is exported, not drawn, and the map keeps the climate's temperatures, until the ice line stands where Earth's does.
- **A world frozen on average.** The locked ocean world, an Earth-like locked to an M dwarf at 0.0485 AU, takes 0.44 of Earth's sunlight and stands at 209.5 K on average, its water solid at its mean: it has no sea. Its sea waits on the working liquid read from the local climate, the next milestone.
- **The renders.** The Earth-like rendered at 2,048 px is byte-identical to the last renders in every map, its ice exported in sea_ice.exr alone; Mars and the Moon, without seas, likewise.

## The twelfth milestone: the seasons' air

- **F245.** 5A's cycle revised (the user, 2026-11-01): the places run in step under their air, which holds its column's heat, exchanges heat with its land and its sea by their shares and carries the transport between places, solved at every step, so that the transport answers the seasons; a sea's water carries a share of it. The Earth-like's mean through its seasons stands at its steady balance; its air swings as Earth's does, 27 K over 80° north and 10 K over 45°; its sea's ice comes and goes from 55 to 70° and lies all year near the poles; it is drawn, and the sea's cells take the sea's temperatures.

How it went:

- **The sea's ice kept 60°.** With the air answering the seasons, the Earth-like's mean held at 287.0 K and its air swung as Earth's does; but over a sea at every latitude its ice lay all year poleward of 65° and covered 11.2 percent of the body at the snapshot. The sea's water took the transport only through its air, and its ice kept the air's heat off it: so a slab ocean without an ocean's transport freezes, as an aquaplanet's does. Carried 30 percent of the transport through its water across the sea's share of each boundary, the ocean's share at 20 to 30° (Trenberth and Caron 2001), the ice withdrew to 75°, came and went from 55 to 70°, and the mean rose to 288.0 K.
- **The transport through thin air.** Williams and Kasting's coefficient lets a slowly turning thin atmosphere carry as much as Venus's: Mars locked to its star took the cap of ten times Earth's. Capped too at its column carried by the surface wind across the body, it collapses onto its night side, at 29 to 33 Pa, as thin air on a planet locked to its star does (Wordsworth 2015). No preset reaches that cap.
- **A resonant body's air.** Solved at the insolation's 720 steps an orbit, locked Mars's air took 75 of its cycle's 82 s, by conjugate gradients preconditioned along each latitude, some 300 steps a solve. Stepped at a fast spinner's year's steps an orbit, and preconditioned by the system's mean along each latitude, solved exactly by its longitudes' Fourier modes, a solve takes 15 steps and the cycle 23 s; a locked Earth-like's air, its own heat outweighing its transport's spread, takes one or two steps, its cycle under 2 s.
- **The frost's own gas's traps.** Under the cycle's air, which brings the frost its heat by the greenhouse's radiation too, the frost's own gas's ice loses mass at every trap: Mars's carbon dioxide some 315 kg/m² a year at its poles, Pluto's nitrogen everywhere. The caps still stand where it loses the least, Mars's carbon dioxide at its poles, Pluto's nitrogen at 40° south.
- **Ice sheets on a latitude.** With the summers the air now gives it, the Earth-like's polar land kept its ice, and its water laid caps, 17 m over 2.3 percent of the body; but their margins ran along a latitude, the traps' air the datum's at every height. The traps' air cooled at the free air's lapse rate turned Mars's basins into worse traps than its heights, against the frost point the basin's pressure raises: near its ground the air follows the ground, and Mars's barely cools with height. Cooled as the ground is (F237), the caps follow the Earth-like's highlands, 39.8 m over 3.9 percent of the body, and Mars's basins trap the better again.
- **The renders.** The Earth-like rendered at 2,048 px draws its sea's ice, 8.4 percent of the body at the snapshot, and its sea's temperatures, its mean 287.1 K, its water's caps 3.9 percent of the body where none lay before; Mars's mean through its seasons falls from 211.3 to 202.6 K, its equator's from 231.5 to 218.7 K, its pressure 507 to 694 Pa, its caps 2.5 percent of the body in water and 0.01 in carbon dioxide; Pluto's and Triton's pressures 0.005 to 0.026 Pa and 0.0004 to 0.073 Pa, Pluto's nitrogen caps in its southern mid-latitudes; Venus and Titan at their steady balances, 755.3 and 91.1 K.

## The thirteenth milestone: the ices that come and go

- **F246.** The user's dynamic caps (2026-11-01): an RGBA map of the ices that come and go through the year, R the phase each freezes at, G the phase it melts out at, B its thickness at its peak, A the phase of its peak, from which an engine or a film draws the caps at any season. Mars's carbon dioxide frost, the Earth-like's sea ice and Pluto's nitrogen frost come and go as their snapshots show them, the windows at the snapshot's phase matching the snapshot's maps on 98 to 99 percent of their pixels.

How it went:

- **The window.** An ice's year is a run of steps it lies through, round the year's end where a winter spans it; taken as the longest run, a mass past half its colour hidden counting as lying, as the drawn maps count it. Its phases are shares of the year's time from the northern spring equinox, so that an engine's clock reads them, and it lies all year as 0 and 1.
- **Phases between places.** A phase cannot be blended across the year's end: read between the cycle's places on the circle of the year, and nearest the cell when exported, while the thickness blends.
- **A marine pole's thick ice.** The Earth-like's north pole lies wholly under sea, and its ice, its top never melting under its column's radiation, thickens past its cycles, 9.8 m at its peak; over a sea that is 71 percent of every latitude it stands at 4.2 m, the land beside it warming the summer's air.
- **The renders.** The Earth-like and Mars rendered at 2,048 px are byte-identical to the last renders in every map they had, the cryosphere map written beside them, and Pluto with them: Mars's thickest 0.94 m, Pluto's 0.08 m, the Earth-like's 9.8 m.

## The fourteenth milestone: nitrogen frost's emissivity

- **F247.** The user's call (2026-11-01): nitrogen's frost radiates at its published emissivities, 0.3 in its alpha phase and 0.8 in its beta phase above 35.61 K. Pluto's pressure through its year climbs to 0.87 Pa, near the 1 Pa measured, where it stood under 0.03; Triton's to 1.20 Pa, near its 1.4.

How it went:

- **The published numbers.** Nitrogen ice's emissivity is its phase's: its cubic alpha phase's far under its hexagonal beta phase's (Stansberry et al. 1996), 0.3 and 0.8 in Pluto's volatile transport models (Bertrand et al. 2019), which with the beta's 0.8 and an albedo of 0.72 to 0.73 reproduce Pluto's pressure (Meza et al. 2019). Taken as published, with the frost's albedo left its colour's 0.80.
- **A step in the frost's balance.** The emissivity's jump at the transition, taken as a step, would make the frost's radiation leap as its point crosses 35.61 K; blended over half a kelvin, its slope enters the frosted ground's linearised step.
- **Pressures near the measured.** Pluto's frost, radiating 0.8 in its beta phase, holds a frost point warm enough for its pressure to climb to 0.87 Pa; in its alpha phase, at 0.3, it gives back its gas sooner, and the pressure falls only to 0.010 Pa, where it fell to 0.005. Triton's climbs to 1.20 Pa.
- **The renders.** Pluto and Triton rendered at 2,048 px: Pluto's pressure through its year from 0.010 to 0.87 Pa, Triton's from 0.13 to 1.11 Pa, their caps where they stood, 1.35 and 1.78 percent of the body in nitrogen; Mars byte-identical in every map, its carbon dioxide frost still a black body.

## The fifteenth milestone: the ground without its seasons' ice

- **F248.** The user's call (2026-11-01): the albedo and height maps hold no seasonal ice, the cryosphere map's thickness giving the frost's height. The albedo no longer draws the snapshot's frost and sea ice, which stay exported apart; the height never held them, its polar domes the perennial caps' ice (F241), terrain that lies all year. The cryosphere map's sidecar gives each ice's colour.

How it went:

- **Two kinds of cap.** The seasonal ices come and go in the cryosphere map, metres at most; the perennial caps the inventory lays (F241) are kilometres of ice that stay all year, Mars's polar layered deposits and the Earth-like's ice sheets: terrain, which the height and the albedo keep.
- **A line of empty windows.** The user saw a thin line of frost at the end of Mars's southern summer, along its southern frost's margin: the thickness, read between the grid's cells, bled a pixel past the windows, read at the nearest cell, so that the margins held thickness with no window. It is read only where the nearest cell's window lies.
- **Split caps tried.** Carbon dioxide frost brightening under its sun (Guo et al. 2010's fit to the Viking orbiters' albedos of Mars's southern residual cap, 0.532 and 8.72e-4 per W/m²) laid Mars's carbon dioxide all year in the south, capped at the fit's 0.80: but the water followed it, its coldest trap under the frost held at its point, and all of Mars's water went south. The caps lay each ice where it lasts longest, as if time had no end; Mars's water stands at both poles from the obliquity's cycles. Not kept: the caps through the climate's history (5B) will lay it.
- **The renders.** The Earth-like, Mars and Pluto rendered at 2,048 px hold no seasonal ice in their albedo and shaded maps, every other map byte-identical but the cryosphere map, whose margins' thickness now stands within its windows.

## The sixteenth milestone: the locked ocean world

- **F249, F250.** The user's call (2026-11-01): the locked ocean world first of 5B's remaining items, and TRAPPIST-1d and e. The working liquid is read at the surface's warmest place, so that a body locked to its star whose mean stands below freezing keeps its sea; and ice and snow brighten the ground by their star's light, little under a red dwarf's. TRAPPIST-1e, THAI's first habitable case, keeps its sea open over 0.57 of its dayside, its mean and its dayside's within THAI's four models'.

How it went:

- **No sea at the mean.** The Earth-like locked to a red dwarf at two thirds of Earth's sunlight stands at a mean below freezing, and its working liquid, read at the mean, was none: no sea at all. Read at the warmest place, an estimate from the year's sunlight's harmonics held by the radiation and the transport, it is water; the tides keep the liquid at the mean, which the warmest place reads through the transport.
- **A snowball.** With its sea, the world froze over, its ice brightening it as the Earth's does, 0.62 cold and 0.55 melting over a ground of 0.30. With the ice's albedos lowered by hand, 0.50 opened an eye to 40 to 50° from the substellar point, and 0.40 to 60°, as THAI's models do with their ice at 0.25; 0.57, sea ice's snow in two bands under TRAPPIST-1's light, still froze over, 1.6 m of ice under its star.
- **The ice under its star's light.** Two bands split at 0.7 µm, with sea ice's model's albedos (CICE: snow 0.98 and 0.70, bare ice 0.78 and 0.36) weighted by the star's blackbody, dimmed snow's contrast only to 0.84 of the Sun's under TRAPPIST-1: its near band's albedo is the solar near infrared's, mostly below 1.2 µm, where a red dwarf's lies further out, where snow is dark and water vapour absorbs. Shields et al. 2013's planetary albedos over snow and over the ocean, by radiative transfer through Earth's air and clouds, give the contrast itself: under AD Leo half the Sun's. A third band from 1.19 µm, as sea ice's three-band model has (Briegleb and Light 2007), its contrast nil and the near band's weighted 0.7, fits their four stars within 0.07; under TRAPPIST-1, 0.28.
- **Against THAI.** TRAPPIST-1e as THAI's Hab 1, Earth's air and water at 900 W/m² on a circular orbit about a star of 0.0898 solar masses: its sea open over 0.57 of its dayside (THAI: almost half), its mean 245.5 K (232 to 246), its dayside's 265.1 (261 to 268), its planetary albedo 0.31 (0.22 to 0.28), the ground's own albedo Earth's under any star until the clouds (5G). The test runs it in 31 s.
- **The caps under the star.** Rendered with the Earth-like's continents, TRAPPIST-1e came out colder than the aquaplanet, its mean 225.6 K and its eye some 40° across: its caps, laid where water's ice lasts, held every land under up to 5 km of ice, the land under its star too, their ice at its clean colour's albedo, under the Sun's light, in the traps that place it and in the climate over it. Water's ice now brightens the ground as the snow does wherever it lies, so that its caps stand only where it lasts under the star; and snow brightens a place only where it is brighter than what lies there, a cap's ice under the star's light, where under the Sun the snow's planetary albedo against a cap's colour decided it.
- **TRAPPIST-1d.** At 1518 W/m² its sea stays open but at its antistellar point, 352 K under its star. Wolf 2017's model found d too hot for liquid water under any atmosphere, while others keep it temperate under a substellar cloud shield, which the model has not yet (5G), nor the moist greenhouse's limit to the radiation, its grey air radiating more the warmer it stands; and d may never have cooled enough to rain out its water at all (Turbet et al. 2023), which the climate's history will ask (5B).
- **The renders.** TRAPPIST-1e, a project file of the Earth-like at 900 W/m² about a star of 0.0898 solar masses, rendered at 2,048 px: the land under its star bare, up to 309 K; its sea open some 60° about its star, its edge on the sea's freezing isoline in the cycle's 4° steps, the 5B gate's first item; its caps, 273 m of its water, over its land beyond the eye and its night, 17 percent of the body, up to 5.2 km thick; its mean 245.0 K, its sea ice over 57 percent of the body at the snapshot. The Earth-like's mean stays 287.0 K, its sea ice 8.49 percent of the body at the snapshot (8.48 before), its water caps 40.6 m over 3.99 percent (39.8 m, 3.93 before), its star's 5795 K brightening its ice 1.003 times the Sun's 5772 K.

## The seventeenth milestone: the ice's edge between the places

- **F251.** The user's call (2026-11-01): TRAPPIST-1e's sea ice ended in square steps where its caps ended smoothly on its land. Its edge now runs between the cycle's places along the pattern of the year's sunlight, smoothed along itself.

How it went:

- **The steps.** The cycle decides the sea's ice place by place, 4° apart on a body locked to its star; read between the places, the ice covered half the sea between a place under ice and an open one, so that its edge traced the places' outline, stepping by up to 1.2° from one bearing to the next round TRAPPIST-1e's star. The caps were smooth, laid on the map's own cells by their traps read smoothly between the climate's cells (F241).
- **A guide.** The year's steady balance's temperature, smooth between the places, was the first guide: but on the cube's grid it is not round about the star, varying by 1.9 K round it at 64°, some 4° of the edge's distance at its slope. The year's sunlight's own diffusive pattern about its axis, F249's harmonics, is round by construction.
- **Smoothed along the edge.** Between each pair of places the edge stood at the midpoint of their guides, still where the places fell: its distance ranged over 2.8° round the star, as read between the places. Smoothed along the edge over three spacings of the places, it ranges over 1.3°, its steps from one bearing to the next under 0.1°.
- **The renders.** TRAPPIST-1e's eye ends smoothly round its star, its sea ice over 56.4 percent of the body at the snapshot (56.8 before); the Earth-like's edges stay latitudes, its sea ice 8.37 percent (8.49); their other maps byte-identical.

## The eighteenth milestone: the cryosphere's export

- **F252, F253.** The user's call (2026-11-01): another Claude, testing the cryosphere map live, wrote a review of it; its claims checked, its changes made before the rest of 5B, the map's design changed where cleaner while the format is young. The frost's windows at its edge read through the margin and close on their peak, and the export holds three maps of one meaning each with a sidecar that says how to read them (docs/adr/0006).

How it went:

- **The review checked.** An agent checked its seventeen claims against the code and a fresh render of Mars: thirteen held as written, two of them understated (slopes carried the margin's window as far as 68° north, and 11 to 13 percent of the margin's pixels held under 6 mm of frost, down to 14 µm, through a full window); four held in part (the frost's series was not kept, only the pressure's; the size the features take; how the wraps were counted). Its third item, the frost's energy albedo, which makes its edge a cliff in the model and its spring late at 64° N (Ls 76 against Titus 2005's Ls 12), goes with the rest of 5B and a review after it.
- **The snap.** Past the last place the frost lay at, every pixel copied its whole window: the margin of F251, smoothed along the edge, now ends it where the sunlight's pattern crosses between the places, and its windows close on their peak toward the edge. Mars's caps' outlines stand still for under 9 percent of the year, where some 30.
- **The design.** Three maps of one meaning each, no alpha: the window's phases, the peak's thickness, and which ice lies, an index into the sidecar's list of the ices, where the liquid mask had implied it. The phases keep the northern spring equinox for their origin: an origin for each hemisphere would soften the filtering it cannot cure, and baked frames will serve pipelines that filter (Phase 7).
- **The curve.** Each ice's year, rising from its window's opening to its peak and falling to its close, fitted over its places: Mars's frost rises nearly linearly and falls as a parabola, the sea's ice rises fast; the year step by step gives the hemispheres' exact means.
- **The rest.** The parameters in their own units, the temperature's extremes, each species' caps, the features apart, a filter for every map, and readings for ports to check themselves by.
- **The renders.** Mars, the Earth-like, TRAPPIST-1e and Pluto rendered at 2,048 px in the new form; Mars's caps' lying area through the year changes by at most 0.75 percent of the body in 1/192 of a year, where 1.38, its sidecar 1.33 MB with 6.35 MB of features apart.

## The nineteenth milestone: the cycle by longitude

- **F254, F255.** The user's questions (2026-11-01): the Earth-like's sea ice ending at straight latitudes and TRAPPIST-1e's in a near-circle, broken only by land; whether physics or noise would vary them; and why only the Earth-like's sea holds seasonal ice. Physics will: the plan now gives the seas' currents and the ice's drift with the wind to 5C, a locked world's superrotation to 5G, the floes' texture to the noise pass after Phase 6; the land's snow waits for the snowpack the moisture brings (5C). The first step is 5B's own: a fast spinner's cycle by latitude and longitude, its windows opening and closing where its ice crosses the mass it lies at.

How it went:

- **The cost.** By longitude at the resonant grid's 4°, the Earth-like's cycle took 11.2 s where 3.7, Mars's 8.8 where 2.6; at 2° in latitude, as the zonal cycle had, twice that. A render's traps, each climate cell's own, cost Mars about a minute. An airless fast spinner, whose places no air joins, keeps its latitudes.
- **The same climate.** On a uniform ground the cycle by longitude gives every body the zonal cycle's means within 0.1 K: only the ground's variation by longitude changes it.
- **Steps again.** By longitude the edges stepped where neighbouring places' windows, opening and closing at the middle of a step, differed by one step along a latitude: at a phase the edge held one latitude over tens of degrees of longitude, 67.9° north, then jumped a row to 63.9°. Opening and closing where the ice's mass crosses the mass it lies at, between the steps' ends, the windows' phases run on continuously. The last steps traced the boundary between the ice that lies all year and the ice that comes and goes, a pixel lying all year where its nearest place did; toward those places the window now lengthens by their share.
- **The renders.** The Earth-like's sea ice curves with its coasts, from 52° north in its enclosed waters to 84° where its Arctic meets land, along a latitude only over its open ocean; its caps shrink from 3.99 to 1.29 percent of the body, its low polar land warming in summer with its own air. Mars, Pluto and Titan within a percent or two of before, TRAPPIST-1e the same.

## The twentieth milestone: one temperature for every stage

- **F256.** The plan's item from 5A: the ice stage's lag and frost and Io's frost reading the climate's temperatures, their sublimation through the days and seasons rather than at a mean. The stages now read the climate's year-mean balance at their own albedo where their slopes face, and Io's frost its vapour pressure through the cycle's days and seasons.

How it went:

- **The day's sublimation.** Integrated through the day, the gross sublimation of water ice runs 10^4 to 10^5 times the rate at the mean temperature: Callisto's equator sublimes as it would at 154 K, its mean 122 K; Europa's 4 x 10^5 times. Taken so, the lag covered clean old ice over 0.78 of its surface and young ice over 0.83 where none lies, and Callisto's lag reached 0.83 of its half-kilometre layer where 0.02, which Europa's brightness, the young craters and Callisto's bright fresh craters deny. Most of what sublimes at noon settles again nearby and the gardening mixes the lag in: the net sublimation is 5E's, and the lag keeps the year's mean temperature, now the climate's.
- **The spreads.** Each place keeps its surface temperatures through its cycle as 64 quantiles. The Moon's equator's run from its night's to its noon's, their mean its mean.
- **Io.** Its frost's noon under the ground's thermal inertia stands lower than the old balance's and falls faster with latitude, its vapour pressure at 30° 0.15 of the equator's where 0.26. Its frost shows over 29.1 percent of the surface where 27.0, the more between 10° and 45°; the plasma's red-brown begins where it did, and the shield's calibration holds.
- **The renders.** Europa's albedo is unchanged and Ganymede's mean the same to 0.0001; Callisto's darkens by 0.0005 on average, a step of the 8-bit map at most and the more toward its poles, its graph's lag a hair wider, its Bond albedo 0.2075 where 0.2080, and through the ice's temperature its relaxation moves a little of its and Ganymede's terrain.

## The twenty-first milestone: the frost's energy albedo

- **F257.** The cryosphere review's third item: the frost's albedo in the energy balance apart from its colour, from measured broadband albedos, the thresholds its energy, its window and its drawn cover read made one, and its ageing.

How it went:

- **The albedo.** Carbon dioxide's frost reflected 0.69 in the energy balance, its colour's mean. It now reflects 0.5, between the measured broadband albedos of Mars's seasonal caps, the dusty north's typically 0.40 and the south's brighter; not fitted to the TES dates, against which it is checked. Nitrogen's and methane's keep their values until the volatile transport tests them.
- **One cover.** The energy took the frost's albedo at the first gram, where the window counts frost past half its cover and the map draws it by its cover; now all three read the cover. The albedo's share changes nothing on its own, as the review found, a step laying more than a cover's frost at once.
- **The emission, tried and not kept.** Radiating by the cover too, a thin nitrogen frost cooled as bare ground and a thick one at its own low emissivity, so that the frost held at about its cover's mass over a third of Pluto and a fifth of Triton and flickered across its window's threshold, 3.5 percent of Triton turning at once, and Pluto's pressure fell from 0.87 Pa to 0.55 against the 1 Pa or so measured. A frost's transparency in the far infrared is not its visible cover; the frost radiates at its own emissivity, as before.
- **Mars.** At 64° north the frost leaves at Ls 39 where 77, against Ls 12 measured; its edges stand at 61.3° north at Ls 26 and 53.2° south at Ls 180, where 54.8 and 49.8, against about 66 and 59; its pressure swings as before. The north's summer holds no carbon dioxide frost for a quarter of its year, as Mars's from Ls 90 to its autumn.
- **The ageing** goes to 5D with the dust: the frost darkens with the dust it gathers, Mars's northern frost some ten times the dustier, and its grains grow into a translucent slab, the south's cryptic region; the measured spring albedo stands for both meanwhile.

## The twenty-second milestone: the caps that hold the air

- **F258.** The frost's own gas's caps as the reservoir its air draws on through the year, the caps laid again in the traps of the year they hold, and Pluto's own season at perihelion.

How it went:

- **The collapse.** The caps were laid after the year and took no part in it: Pluto's winter froze its air out to 0.011 Pa, its year's mean a sixth of its measured pressure, and a nitrogen ice held at that air's point lost some 1,100 kg/m² a year everywhere, so that no basin could gather it.
- **The reservoir.** Each place's share of the caps now runs on its own column at the frost point, giving the air what it sublimes and taking what the air leaves it, the air's mean the body's. Kept as laid each year, the caps made gas from nothing, which piled up as frost; booked with the frost, the buffered air's frost took fifty of Pluto's years to come to its own; and a step that read how much of the gas added the air keeps set Triton swinging between 0.4 and 37 Pa. The seasons are spun up on the latitudes' rows instead, at a ninetieth of the cost, the places started from them.
- **Pluto.** Its air runs from 0.25 to 2.71 Pa through its year where 0.011 to 0.87; it doubles from Ls 2 to Ls 37 as Pluto's did from 1988 to 2002, but peaks there where Pluto's went on rising to 2015, and its mean is the preset's 2015 value read as its mean. Its nitrogen gathers in its two great ringed basins and in craters near its equator, where its tilt leaves the least sunlight and a basin's air raises its point, as the volatile transport models find, at its measured perihelion, Ls 3.8, where Earth's 283° stood.
- **Triton and Mars.** Triton's air averages its 1.4 Pa where 0.44; its caps stay at its poles. Mars's seasons are as they were; its small cap of carbon dioxide, losing everywhere at its albedo, now lies on the northern water cap, where on the southern pole: the residual cap's place waits for the frost brightening under its sun with the caps' history.
- **Not yet:** the caps' exchange at their own heights in the year and their gathering through the history, the volatile transport. Tried after: the pressure read at the season it was measured, Pluto's July 2015 1.25 Pa at its Ls 63.6; where the frost buffers the air the air at a season is the frost's balance's whatever the gas's mass, and it stood at 1.76 Pa as the year's gas drained. A pressure measured at such a season is a test of the ice's energy balance, not an input; nitrogen's frost albedo is measured, 0.80 (Hofgartner et al. 2023).

## The twenty-third milestone: the ices gathered through the age

- **F259.** Each ice of the inventory that is not the body's working liquid gathered through the body's age from a layer over the whole body, at the rates its traps give at its surface's height, as the seasons turn; its caps laid where it gathered.

How it went:

- **The gathering.** The traps are run on the latitudes' rows at a few heights of the air's pressure, at the present season and at each season at perihelion as the seasons turn; the frost's own gas, moving in years, follows the present, a slow ice the turning's mean. The ices gather from an even layer through the body's age, each exchanging where it lies uppermost, the frost's gas by its growth, any other ice by its saturation pressure against the air's vapour.
- **Mars.** Its water gathers through its 4.5 Gyr all at its north pole, as the caps were laid before. On flat ground the turning's mean holds its two poles alike, but its relief makes the north the better trap, and over the age the south drains to it: Mars's southern layered deposits keep under the lag a sublimating ice leaves, which waits for 5E.
- **The obliquity's swing, tried and not kept.** At 35° Mars's summer pole is its warmest place; the swing's mean, ruled by its high epochs, laid half of its water in its southern mid-latitude basins. A polar deposit keeps through a high epoch under its lag: the swing waits for the lag, in 5E.
- **Pluto and Triton.** Their nitrogen, moving in years, lies where the present season puts it, as before.
- **Not yet:** the swing and the lag (5E); the caps' exchange at their own heights in the year.

## The twenty-fourth milestone: the star and the climate through the body's life

- **F260.** The star's light through the body's life by the BHAC15 tracks up to 1.4 solar masses, a low-mass star's bright contraction before its main sequence among them; the mean surface temperature through it under the present air; the time until which the body stood past the runaway greenhouse's flux.

How it went:

- **The star.** The tracks' share of the star's present luminosity at each time. The Sun shone at 0.69 of its present light when its planets formed. A red dwarf of TRAPPIST-1's mass stands on its track at 5.50e-4 of the Sun's light, against its measured 5.66e-4, and stood at 6.4 times that at 50 Myr.
- **The climate.** Under its present air the young Earth stood at 249 K, below freezing until 2.52 Gyr after it formed: the faint young Sun's paradox, which the air's own history answers (F261).
- **The runaway greenhouse.** By the star's effective temperature and the body's mass (Kopparapu et al. 2014), where the present albedo, a guess by pressure class, is not the steam's. TRAPPIST-1e ran away until 0.34 Gyr after it formed, where Bourrier et al. 2017 find TRAPPIST-1's planets' runaways ending after a few tens to a few hundred Myr; the Earth never did.
- **Not yet:** the air's history (F261); the albedo's and the star's temperature's histories.

## The twenty-fifth milestone: the star's XUV by its mass

- **F261.** A star's XUV light against its whole saturated at a level by its mass until a time by its mass, then falling as the age to the -1.23: a red dwarf's saturated for Gyr at 10^-3.05 of its light, the Sun's for 0.1 Gyr at 10^-3.68. The XUV at the body through its life.

How it went:

- **The measures.** The red dwarfs' saturation time measured on a volume-complete sample (Pass et al. 2024), the F, G and K stars' and their saturated levels by colour (Jackson et al. 2012), the M dwarfs' level (Wright et al. 2018); the Sun's history as before.
- **TRAPPIST-1.** A star of its mass at its age sends 8.5e-4 of its light in XUV against the 6e-4 to 9e-4 measured; TRAPPIST-1e takes 165 times the Earth's XUV today and some 2000 times when it formed.
- **Not yet:** the escape it drives (F262, F263); a star's rotation through its life.

## The twenty-sixth milestone: the water lost to space in a runaway

- **F262.** While its star held it past the runaway greenhouse's flux, a body's water's hydrogen escaped at the energy limit, a tenth of the XUV it took, dragging oxygen with it past the crossover; the water lost, the water's history from its present inventory, and the oxygen left.

How it went:

- **TRAPPIST-1e.** It loses 12.5 Earth oceans through its 0.34 Gyr in the runaway, from 36.8 km of water to its present 2.7, and leaves 1805 bar of oxygen: Luger and Barnes's several oceans and up to 2000 bar. Bourrier et al. 2017, with a third of the model's XUV held at a constant share of its light, find 2.1 to 2.5 oceans by its entry into the habitable zone.
- **The Earth** keeps its water, never past the runaway's flux.
- **Not yet:** the photolysis's and the diffusion's limits; the oxygen's sinks; Venus's water.

## The twenty-seventh milestone: the air lost to space through the body's life

- **F263.** Every body loses its air as Mars does today, measured by MAVEN: its photochemical oxygen as the EUV to the 2.6, by the share of hot oxygen its gravity lets go, and its oxygen ions as the EUV; at most the energy limit. The air's escape today and through the life, and the pressure history it gives.

How it went:

- **The choice** (the user): Mars's measured rates applied to every body, scaled by its gravity and XUV, over thermal escape alone, which misses Mars's main loss.
- **Mars.** 1.73 kg/s today; 0.36 bar 3.5 Ga, against Lillis et al.'s extrapolated 0.48; 0.44 3.6 Ga, under Kite et al.'s 0.9 bar; 0.69 3.8 Ga and 1.28 4.0 Ga. Before some 4.4 Ga the energy limit holds and the air the present would need climbs to 52 bar: an upper bound the outgassing lowers.
- **Earth and Venus** keep their hot oxygen and lose their ions, at some 4 and 30 times their measured rates, a small share of their air; **Pluto and Triton** keep their pressure, their nitrogen ice paying.
- **Outgassing, staged.** Applied with the Earth's measured juvenile yields, Mars's model outgasses some 0.6 to 0.9 bar after 3.5 Ga, as much as it loses: the budget needs its stores of carbonate, and the Earth's its weathering, so the outgassing comes with them (F264).
- **Not yet:** the outgassing and its stores; light gases' thermal escape; sputtering; hydrogen's escape outside a runaway.

## The twenty-eighth milestone: the air's outgassing by its present balance

- **F264.** Each body's melt gives its air at the one yield at which its melt today gives what it loses to space today, its volcanism's history carrying it back: the air's history from its escape and its outgassing.

How it went:

- **The Earth's yields** on Mars's model give some 0.9 bar after 3.5 Ga, more than it lost and holds, a budget only stores of carbonate would close. **The user** asked for something simpler that fits every body roughly.
- **The present balance.** Mars's melt gives 101 ppm of carbon dioxide, within what its meteorites' redox allows; it held 0.29 bar 3.5 Ga, 0.61 3.8 Ga and 1.20 4.0 Ga. The Earth's and Venus's yields are tiny, their air hardly touched; Pluto's and Triton's ice buffers both ways.
- **Tried and set aside:** an air built from nothing by the outgassing, which leaves Mars without air until the present; the climate under the air's history, whose grey greenhouse gives early Mars 289 K under 0.61 bar where 1-D and 3-D models find 210 to 225 K.
- **Not yet:** a melt's gases by its mantle's redox and the carbonate stores, with 5D; the greenhouse of a thick cold carbon dioxide atmosphere.

## Two items moved, the gate next

- **Moved** (the user, 2026-11-02): Io's frost to its coldest ground to 5F, with Io's surface on the general pipeline that replaces the pipeline it would be built in; the caps' exchange with the air at their own heights to 5E, a small effect (some 0.03 K a kilometre for Pluto's nitrogen, 0.7 for Mars's carbon dioxide). 5B's gate check and its review page next.

## The 5B gate check and its review page

Rendered at 2,048 px at 7a94bd4 (Mars, Pluto, Triton, TRAPPIST-1e and the Earth) and measured on the maps:

- **The locked ocean world's sea-ice edge on the freezing isoline.** TRAPPIST-1e's open sea's last cells beside the ice stand at 271.5 K (10 to 90 percent of them 271.4 to 271.7 K), seawater freezing at 271.35 K; its sea open over 50 percent of its dayside, where THAI's models leave almost half. Meets.
- **Mars's water ice at its north pole.** 1.90 percent of the body, all beyond 60°N, on ground lower than 82 percent of the body's. Meets.
- **Mars's carbon dioxide wherever its frost's balance holds it.** 0.015 percent of the body, on the northern water cap at 83°N; the physics able to give either hemisphere (the user, 2 November). Meets.
- **Pluto's nitrogen gathered in its basin.** 1.34 percent of the body, 66 percent of it within 30° of the equator, on ground lower than 97 percent of the body's: its two great ringed basins and their craters. Meets.

**Found by the check:** TRAPPIST-1e at its measured 7.6 Gyr (Burgasser and Mamajek 2017) cools into a stagnant lid (34 mW/m², a 224 km lithosphere) and its shields and floods pile up to 153 km on an Earth's gravity, where its relief limit is 9.5 km, the stagnant lid's shields and floods not held to it; at the Sun's 4.57 Gyr it keeps a mobile lid and relief of -5.5 to 6.5 km. No preset is old enough to show it. The review page draws TRAPPIST-1e's maps at 4.57 Gyr and its history at 7.6.

**The review page:** https://claude.ai/artifact/MPH4ptNVN4GruAwFLmv2Se, awaiting the user's review.

## The twenty-ninth milestone: a shield's flanks shed no more than they hold

- **F265.** Found by the gate check (the user: fix it before 5C): a failing flank sheds no more than the toe beyond its innermost scarp, the rest of the supply erupting as the region's floods, and no aureole lobe stands thicker than a fifth of its fall, a thicker one spreading wider and farther.

How it went:

- **The old Earth-like.** At TRAPPIST-1's measured 7.6 Gyr its lid stalls and its great shields shed all the supply their crust could not hold as landslides 30 to 280 km thick. Its relief now runs -9.1 to 12.6 km where it ran to 153, its temperatures 184 to 310 K where they fell to -550.
- **Olympus Mons.** Its aureole's lobes stand 0.9 to 2.0 km where 2.4 to 6.5, and run up to 600 km, nearer the measured 1 to 2 km and 700 km; its scarp as before.
- **A shield beyond its foot.** The wider aureoles took a shield's lobed outline past its foot, where its form turned back and drew no number at 2,137 of Mars's cells: past the foot the share now keeps rising.
- **Not yet:** a great shield's own slope where its crust holds little.

## The thirtieth milestone: the ice's edge at the coasts

- **F266.** The user saw a hard, nearly right-angled turn in TRAPPIST-1e's seasonal ice and temperature where its open sea met its northern continent's tip. A place now stands at the ice's edge by the share of its ground it and a place of the other kind beside it both hold, and a point reads the guided edge by how much of the places about it stand there, the places' own ice for the rest, on its own ground only.

How it went:

- **The cause.** The sea's margin read its guided edge only where a place stood beside one of the other kind; a coastal place under ice stood beside open sea by a sliver of 0.01, so the points reading it opened their sea and their neighbours kept their ice, the switch along the places' column and row.
- **The coast now.** The sea at the continent's tip keeps the ice the cold land holds it under and thins into the open sea over some ten degrees.
- **Elsewhere.** The Earth's sea ice 12.30 percent of its sea where 12.26, its changes along its coasts.
- **Tried first:** the edges weighted alone, and the reading on its own ground alone, which kept the switch.

## The thirty-first milestone: the climatological wind

5C begins (the user, 2026-11-02, "Start with 5C"), its order in PLAN.md: the wind first, which the moisture, the currents, the sea ice's drift and 5D's sand read; then the moisture, the snowpack, the currents and the ice's drift, the rivers, the glaciers, the slopes' wasting and Titan's methane cycle.

- **F267.** The surface wind at each of the cycle's places at each step of its year. Its pressure follows the air: along each latitude the air's departure from the latitude's mean, a warm place's air lower (Lindzen and Nigam 1987); and the latitudes' belts, rising from the thermal equator through the cell that overturns from it, its edges at the Held and Hou reach from the equator (Lindzen and Hou 1988), the eddies' belts past them at the Rhines scale, as strong as brings the surface's torque on the air to nothing; where the cell reaches the pole, its own westerlies a force where its descending branch brings its momentum down. The ground's drag and the spin turn it across the isobars, and the work of the atmosphere's heat engine sets its strength. The year's mean is exported as `wind.exr` and `wind.png`.

How it went:

- **The belts.** The Earth-like's trades end at its cells' edges at 32.6°, its westerlies strongest at 42° and its polar easterlies at 60°; a faster spinner's belts are the narrower and the more.
- **A locked world.** TRAPPIST-1e's air flows from its night side and its poles to the point under its star. A first scheme put a ridge of pressure where a slowly spinning body's cell has its own westerlies and sent its high latitudes' air to the poles; the cell's westerlies are now a force on the air, and its surface branch flows toward the thermal equator all the way.
- **A cell heated far from its equator.** Mars's air is warmest at its pole in its southern summer, and the cells followed it there. A cell heated off its equator keeps its edges while its rising branch moves; Mars's thermal equator now stands within 42° of the equator, its solstices' cells crossing the equator each way in turn, and Titan's, whose cell reaches its poles, runs from pole to pole.
- **The speeds.** Over the body in the year's mean: the Earth-like 4.0 m/s, the weather's gusts apart; TRAPPIST-1e 8.2; Mars 3.9 (its landers some 2 to 7); Venus 0.52 (Venera some 0.3 to 1); Pluto 1.6; Titan 0.11, where Huygens measured some 0.3: its air all but isothermal under its transport's cap (F33), the contrast 5A left to the wind field.
- **Known limitations.** The Earth-like's westerlies narrower than Earth's and some 5 to 8° nearer the equator, the Rhines scale read at the surface's mean wind; no easterlies at the equator, where the air's momentum would carry the trades'; the exchange of heat and vapour with the ground still at a stated 5 m/s; no winds down cold slopes; the jets aloft wait for the weather model (5G).

## The thirty-second milestone: the moisture

- **F268.** The working liquid's vapour through the cycle's year at its places: evaporated from the sea and the wet land at each place's wind and the weather's gusts about it, carried by the wind and mixed by the eddies, rained out as its column nears saturation as Bretherton et al. (2004) measured over the tropical oceans, the land holding its rain in a bucket. The year's precipitation and aridity are exported.

How it went:

- **The Earth-like.** 1,182 kg/m² a year over the body, near Earth's some 950 to 1,000; its land's runoff 352 kg/m² a year, Earth's some 270; its deserts on its continents' subtropical latitudes, its rain at its thermal equator and under its westerlies.
- **The locked world.** TRAPPIST-1e's rain falls under its star and none on its night side.
- **Three fixes.** The evaporation and the rain solved together on each step's end, where the storm tracks' seas had taken vapour from columns left supersaturated within a step; the gusts about the mean wind, without which the subtropics' calm ridges stopped evaporating; the rain bracketed and read from the column's balance, where a column over cold high ground, a hundred times saturated, had rained 10^20 kg/m² a year.
- **Known limitations.** The rain at the thermal equator too narrow and too strong; no rain shadows yet, the relief only at the places' 4°; Titan's methane some hundred times too much and wettest at its equator, its evaporation drawing on no energy and its air mixed too hard: the next milestone puts the latent heat into the energy balance at each place's wind, before Titan's transport (the user, 2026-11-03).

## The thirty-third milestone: the year at its own wind

- **F269.** The exchange of heat between the ground and its air at each place's own wind, where it read a stated 5 m/s, and the latent heat of what wet ground gives the air tying it to its air, for any working liquid: the body's first year at the stated wind, its wind and moisture from it, and every year after at each place's own.

How it went:

- **The Earth-like.** Its mean 285.5 K, within 5A's gate; its sea ice and its water much as before, its land's evaporation 485 kg/m² a year, Earth's some 490.
- **The locked world.** Its warmest place 12 K cooler under its strong winds; THAI's measures within or near their ranges, its open sea 0.49 of its dayside, nearer THAI's almost half.
- **Mars, Pluto, Venus.** As before: their air's pressures and their caps within a percent or so.
- **Titan.** Its ground loosened from its air at its slow wind, its contrast 0.5 K where 0.1, short of the measured 3 K; its methane still far too much.
- **Tried first:** the latent heat beyond a saturated air's given to the air over the ground, which warmed the air near the ground and choked its own evaporation.
- **Known limitations.** The grey balance has no budget of the ground apart from its column's, nor its air any radiation of its own: the latent heat of the ground's evaporation is neither limited by the sunlight the ground receives nor radiated where it falls. With Titan's transport (the user, 2026-11-03).

## The thirty-fourth milestone: the rain over the relief

- **F270.** The year's precipitation on the finish's grid over its relief, where it was the places' read between them: air carried up rising ground condenses its vapour while it is saturated, the more often the higher the ground stands above its place's (Smith and Barstad 2004, read through Bretherton et al. 2004), carried on before it falls; descending air dries what falls; and the air carried on is short of what it rained, and rains the less downwind until it has made it good: the rain shadows.

How it went:

- **The Earth-like.** Its land rains from 0.57 to 1.45 of the places' rain, its windward slopes wetter and its lee drier; its wettest cell 10,300 kg/m² a year, near Earth's wettest. Its relief is plateaus some 1.2 km over its sea; its ranges come with Phase 6's plates.
- **The locked world.** Its ice sheets wetter where the wind climbs them, drier where it descends toward the star; the plateau behind its continent's windward edge in its shadow.
- **Titan.** Unchanged: its slow wind lifts its air too little.
- **Tried first:** the share of the time the air is saturated read at the places' own height for all their ground, too weak to move the rain by more than a tenth; each place keeping its rain, which dried the ground upwind of a wet edge.
- **Known limitations.** No flow round the relief nor waves aloft; the places' own rain does not read their relief's, the shadow carried on the finer grid alone; one lifting wind for the year.

## The thirty-fifth milestone: the snowpack

- **F271.** The snow the moisture brings lies on the land through its year where water is the working liquid: it falls by the air's temperature, holds the ground at its melting point while it melts, its cover brightening the ground, and its mass in the cryosphere's maps beside the frost and the sea's ice.

How it went:

- **The Earth-like.** Its northern snow peaks just past its spring equinox at some 14.5 kg/m² over its hemisphere, Earth's some 12, and lies some part of the year on nine tenths of its land past 40°; its climate as before, 285.8 K.
- **The locked world.** Snow on three tenths of its land; its mean as before, 246.1 K.
- **Tried first:** the snow at the full brightness of ice and snow until it melted (a cold Earth-like, 268.6 K); the ground held at the melting point through all of a step on a trace of snow (the tropics 22 K too cold); the melting without the freezing's heat in the air (the locked world 1.4 K too cold).
- **Known limitations.** No rain on snow, no refreezing, no sublimation, no snow on the sea's ice; Titan's methane snow with its methane cycle; the perennial snow kept under 1,000 kg/m² until the glaciers carry it.

## The thirty-sixth milestone: the seas' currents

- **F272.** The wind's stress drives the seas' circulation in Stommel's balance: gyres whose return flow presses against their basins' western coasts, the land holding them by a friction that grows as the sea's share of a place falls; the flow carries the seas' heat, and the diffusion keeps only what the wind does not drive.

How it went:

- **The test basin.** A gyre of 18.7 Sv where Sverdrup's balance gives 22.4, its western current poleward.
- **The Earth-like.** Gyres of some 45 to 49 Sv; its seas warmer toward their basins' east, its ice's edge rising some 3 to 5° there; its climate as before, 286.0 K.
- **The locked world.** THAI's measures as before.
- **The choice:** the user chose the total to emerge, the diffusion keeping the measured share the wind does not drive, the stated third's rest if that failed; it held the gates.
- **Known limitations.** The Earth-like's ice edge bends far less than Earth's 45° to 80° across the North Atlantic: its basins are wide and simple; no overturning by density, no flow round islands, no seasons of the currents; the ice's drift next.

## The thirty-seventh milestone: the climate's layers on the terrain's grid

- **F273.** The climate's places read smoothly on the finer grid: the year's sea surface without the rectangles the ice's margin cut in it, the rain's background and its potential evaporation without the places' cells, and the land's snow run again on the terrain's own heights.

How it went:

- **The year's sea surface.** Read smoothly between the places, no edge drawn through its year's mean: the rectangles along the ice's edge gone.
- **The snow.** Each fine cell's year from the places' precipitation and air lapsed to its height; on the Earth-like from some 40° poleward, deeper and longer on its high ground, perennial on its northern heights; on the locked world on four tenths of its land.
- **The snapshot's sea ice** keeps F251's margin, round about a locked world's star (0.10° steps where the sea's heat read smoothly stepped 0.83°).
- **Tried first:** the sea ice from the sea's heat; the melt over the air's degrees where a trace of snow lay (the snow never melted); the height over the finish's datum, not the sea.
- **Known limitations.** The snapshot's sea ice map still shows its margin's places; the snow's sunny slopes and the orographic snow wait.

## The thirty-eighth milestone: a window's ends within their steps

- **F274.** The cryosphere's windows no longer close at the ends of the climate's steps: where the step after a run ends with none of its ice, the ice goes through it no slower than through the run's last step, and comes so at its start. From the user's demo of the maps and a second session's findings and patch.

How it went:

- **Mars.** Its melts spread through their steps as they would evenly, its polar clearing in no jumps but where its polar rows clear together; its 64° north frost arriving nearer Titus's Ls 225°.
- **Known limitations.** The northern autumn's freezing in one frame, from the water cap's frost subliming all it holds every other step, waits for the frost's step split where it runs out.

## The thirty-ninth milestone: the frost's step where it runs out

- **F275.** Where a step's frost runs out within it, the rest of the step is the bare ground's: it goes free at its own gain, and frosts again where that would cool it past the point. A frost darker than the cap under it no longer leaves it bare every other step.

How it went:

- **The cause, confirmed.** With the frost's gain held at the bare ground's for one render, Mars's autumn frame was gone: the frost darker than the water cap, subliming all it held under the low sun.
- **Mars.** Its northern autumn's largest frame 8,078 pixels where 46,577; its air and its frost through the year as before.
- **Known limitations.** The window still estimates a melt's place in its step (F274) rather than reading the share the step now finds.

## The fortieth milestone: the sea ice's drift

- **F276.** The wind drives the sea ice: its free drift some 2 percent of the wind, turned by the spin; where the wind presses it together the ice resists with Hibler's strength, as a fluid that bears a pressure up to its strength and none drawn apart, so that it piles against the coasts the wind drives it on until its pile holds the wind, and opens freely; the ice the drift brings to open water melts there.

How it went:

- **The tests.** The free drift at 1.8 and 2.1 percent of a 3 and a 10 m/s wind, turned 39° and 12°; a metre's ice blown on a coast piling to 2.25 m where the strength's balance gives 2.0 to 2.25, its pressure rising at the wind's stress.
- **The pressure, converged.** With 50 plain sweeps a step the pressure never settled, each step's carried from the last's, and one of the Earth-like's places by its south pole swung between two years' ice for all 40 of the lockstep's years, the climate stage near twice as long. Over-relaxed sweeps that settle end the swing: the year settles in 17, the stage some 17 percent longer than without the drift.
- **The Earth-like.** Its ice drifts at 7 cm/s; its thickest ice thinner, its extent as before: its polar easterlies carry its edge's ice poleward, and the weak westerlies round its poles turn the ice out from them.
- **The locked world.** Its ice drifts toward its star and melts in its eye: open over 0.17 of its dayside where 0.26, as Yang et al. 2019 found sea-ice drift shrinks a locked world's open sea; its ice piles to 19 m against the coasts that face its night.
- **Known limitations.** The water at rest under the ice; no shear strength; the Earth-like's polar easterlies from 53° (Earth's from some 60° to 65°) and the weak westerlies round its poles, from the wind (F267), which turn its ice out from its poles.

## The forty-first milestone: the rivers and the lakes

- **F277.** The land's water runs down the terrain: the basins found by the priority flood and merged at their passes, a lake in each where its water balances what its surface loses to the air, filled and spilling at its pass, the rivers carrying it to the lakes and the sea. The moisture's bucket takes the snow as it melts, so that snow lying all year sheds nothing.

How it went:

- **The tests.** An island sheds all its runoff to the sea; a bowl's lake stands where its balance puts it, dry and wet; a full lake spills into its neighbour.
- **The snow's water, balanced.** The bucket first took the year's snowpack's melt beside its own rain, the snowpack fed by the first year's snowfall: TRAPPIST-1e's land shed 982 kg/m² a year of 518. Its snow kept in a store of its own, melting by the snowpack's share, the balance closes: 268 of 464.
- **The Earth-like.** Its land sheds 1.04 million m³/s, its largest river 128,000; but a lake stands in each small pit of its terrain, 19 percent of its land under lakes, as no river has shaped it.
- **The locked world.** Its rivers run under its star alone.
- **Known limitations.** The fluvial erosion that would shape the terrain (next); the seas' connection; the lakes not yet in the climate or the colour.

## The forty-second milestone: the rivers' erosion

- **F278.** The rivers cut their beds by stream power, implicitly from their mouths up, through 2 Myr of today's climate in the finish (the user's staged plan: the coarse loop's history next), laying what they cut in the basins their water ends in, by the flood's tree as their lakes settle, the ground rebounding by isostasy over the lithosphere's flexural wavelength. The Earth-like's terrain drains: its pits fill, its lakes from 19 percent of its land to 4.3.

How it went:

- **The tests.** A wet island's pits drain, its land lowering by what did not rebound of its loss; a dry plain's broad basin stays closed while the small pits about it fill.
- **The sediment, kept.** Laid in each lake only to its water's level, a third of what the rivers cut was lost; passed over the basins' passes one by one, it ran back and forth between two full halves for all of its rounds. Laid by the flood's tree as the water's lakes settle, by volume, it is all laid or brought to the sea.
- **The rebound.** Holding the relief over the flexural wavelength whole lifted an eroding peninsula's interior by its coasts' whole loss; the ground now rebounds by the crust's density over the mantle's, the rest lowering it.
- **The Earth-like.** Its rivers branch across its continents; 4.3 percent of its land under lakes, 21 percent closed, Earth's 2 to 4 and some 20.
- **The locked world.** Its eye's rivers end in terminal lakes in its peninsula's hyper-arid south.
- **Known limitations.** The erosion not yet through the history in the loop; no deltas or fans; one erodibility for every ground; Titan's methane rivers too wet.

## The forty-third milestone: the rebound under the sea

- **F279.** Asked for the height, shaded and albedo maps after F278, the height map of TRAPPIST-1e showed a pale rim along its peninsula's west coast: ground that had stood at the sea's level raised by up to 2.6 km, the ground inland of it cut by 300 to 700 m. F278's rebound had averaged the ground's loss over the flexural wavelength over the land alone, so that a coast's dry ground rose by nearly the crust's share of its outlets' whole loss each step. It now averages the loss over the plate, the sea's floor losing none.

How it went:

- **The search.** The finish's height is the erosion's. A run printing each step's most risen cell found it neither cut nor laid, lifted 200 to 350 m in each of the ten steps; the sediment's deepest fills were closed basins filling toward their spills, their volume kept.
- **The flexure.** A wet island falls by near the crust's share of its loss under a plate whose wavelength is short beside it, by more under one whose wavelength spans it.
- **The bodies.** The Earth-like's land rising more than 500 m from 0.46 percent of it to 0.21; TRAPPIST-1e's coast's ridge from up to 2.6 km to 1.1.
- **What it left.** Ground whose rain falls short of its place's evaporation sheds none, F277's runoff clipped at none, and is never cut; beside rivers that cut their cells to grade it rises by their rebound, by up to 1.5 km on TRAPPIST-1e's peninsula.
- **The lesson.** F278's check printed its most risen cell, 2,085 m, and did not follow it: a terrain's change is read at its extremes before a milestone closes.

## The forty-fourth milestone: the caps' true distance

- **F280.** The user saw stretched barcodes on TRAPPIST-1e's southern cap, fanning from its margin and rippling toward its pole. A cap's distance in from its margin was summed along the grid's four steps, running long across them: its domes diamonds, its surface ruled along the grid and stepped where one point of the margin gave way to another. Each place now carries the margin's point it rises from and the true distance from it (the other session's patch).

How it went:

- **The test.** Half way in from the margin of a disc of ice, the dome stood 7.7 percent thicker across the grid than along it; it now stands within 0.4 percent. Its centre, on a grid of an even count of cells a face, lies half a cell's diagonal off the pole, which the summed distance's excess had hidden; it is read at its own centre, not by a finer grid.
- **The bodies.** The caps' surfaces stand smooth where they were ruled: Mars's northern caps in their basins, TRAPPIST-1e's southern cap, its surface falling from its crest to its pole by 75 to 110 m a degree, steepening toward its margin, where it fell in runs of 44 m broken by steps of 120 to 200 m, and the Earth-like's. The domes stand thinner on the ground they cover: where a liquid's ice covers all the ground it may, it holds less, TRAPPIST-1e's caps 248 m of its water's global layer where 280, at most 4,775 m thick where 5,234, the Earth-like's 3.3 m where 4.2; where the inventory sets them they cover more ground, Pluto's nitrogen 1.65 percent of the body where 1.53, at most 6,826 m thick where 7,390, its air 0.25 to 2.73 Pa through the year where 0.24 to 2.72; Mars's water 2.53 percent where 2.52, its air through the year unchanged, 506 to 702 Pa.
- **What it left.** Creases where the margin's ground steps between its points; the blocky edges of the caps whose gathered places rank ahead of their traps, next (F281).

## The forty-fifth milestone: a working liquid's ice where its own ground keeps it frozen

- **F282.** The user saw grid squares, larger than the map's pixels, on the edge of TRAPPIST-1e's southern cap inside its peninsula's western side. Its water, the body's working liquid, covers all the ground where its warmest keeps under its melting, and that warmest was its place's on the climate's grid at the place's mean height. Each cell's warmest now stands from its place's as the year's peak at its own height and slope stands from the peak at its place's height (F237's downscaling): the edge follows the ground.

How it went:

- **First, the gathered caps.** The other session's diagnosis of Mars's blocky edges, the places an ice gathered at ranked ahead of all its traps, came first (F281): each region an ice gathered in, holding its share, grown over the coldest of its traps beside it. Its tests passed; on Pluto, its methane laid 1.1 m of its 5, the regions that had gathered the most hemmed into eight cells each by the floods of regions whose seeds were colder. Kept aside as a patch for a design that grows the regions fairly and passes on what one cannot hold; at preview quality Mars's caps fill their basins, their edges the basins' rims either way.
- **TRAPPIST-1e's edge was not that.** Its water's caps cover all the ground they may, so their edge is the line where its warmest reaches its melting.
- **The test.** On the Earth-like, ground 40 K colder than its places in the tropics holds its water's ice there, none without.
- **The bodies.** TRAPPIST-1e's edge follows its ground, in patches where the unresolved craters' roughness sets each cell's height; the Earth-like's land ice thinner on more ground.
- **The coast's steps.** The steps along TRAPPIST-1e's coast that the user saw are the continental crust's outline on the coarse loop's grid, carried to the finish: the same before the rivers' erosion, the same on the ground at a map twice as wide.

## The forty-sixth milestone: the continents' outline on the finish's cells

- **F283.** The stair steps the user saw along TRAPPIST-1e's coast were the continental crust's outline, its units taken from the loop's nearest cells, 72 km across: a strip of oceanic crust stood above the sea between it and the coast, which follows the bedrock smoothed over the loop's cell. The finish draws the outline again from the crust smoothed as the bedrock is, continental past the two crusts' mean (the user, 2026-11-05: before the plan goes on).

How it went:

- **Telling them apart.** At a map twice as wide the coast's own edge is a pixel's staircase while the continents' unit steps by the loop's cells, the same size on the ground.
- **The test.** A margin across the loop's grid at a shallow angle: its units 0.48 of a loop's cell off it taken from the nearest, 0.11 drawn again.
- **The bodies.** TRAPPIST-1e's coastal strip gone but for single cells; the Earth-like's coasts as they were; heights, rivers and lakes untouched.

## The forty-seventh milestone: the snapshot's phase

- **F284.** The other session found the sidecar's snapshot phase at its season's own, 0, while its maps are the record of the step that holds the season, taken at the step's end: a port reading the cryosphere's windows there disagreed with `frost.png` on 4.9 percent of Mars's pixels where at the step's end on 3.2. The sidecar now gives the step's end, the windows' own convention.

How it went:

- **The test.** The snapshot's phase is its step's end, after its season's by less than a step, on a cycle of one orbit and of three.
- **The stash.** Checking the formatting by setting the tree aside, on a tree with nothing to set aside, brought back the oldest experiment kept aside instead (milestone 39's cracks); the user restored the two files it touched. The formatting is checked against the file's own count of differences now.

## The forty-eighth milestone: the eddies' belts at their own wind

- **F285.** The user chose the wind before the erosion through the history (2026-11-05). The wind's belts past its cells' edges were spaced at the Rhines scale of the body's mean wind, some 4 m/s on the Earth-like: its westerlies from 32 to 50 degrees, polar easterlies from 52, westerlies again about its poles, its trades weak. The belts are now spaced by their own wind, the flow the spin's beta arrests: the Earth-like's westerlies to 56 degrees, its polar easterlies from 56 and weaker, none westerly about its poles, its trades a quarter stronger.

How it went:

- **The aquaplanet.** Its westerlies reach past 55 degrees, its easterlies blowing by 65; read at the body's mean speed they turned by 55.
- **The bodies.** Mars's winter westerlies as before; TRAPPIST-1e, Titan and Venus unchanged, their cells reaching past where belts would stand.
- **Downstream.** The Earth-like's gyres and sea ice slower with its polar easterlies, its rain a tenth wetter with its trades.
- **Known.** The westerlies' greatest speed at the middle of their belt, 42 degrees, where Earth's lean poleward; the Earth-like's rain some 15 percent past Earth's.

## The forty-ninth milestone: the rain between its places and the caps by the sea

- **F286, F287, F288, F289.** Asked about the maps (2026-11-05), the user saw TRAPPIST-1e's rain end in steps, its land show through the ice on its frozen side, and its caps and snow in different shades. The steps were Catmull-Rom's cubic ringing below none beside its wet eye, held at none in rectangles: the places' fields are read now by a monotone cubic. The bare land was the shelves the sea bared as it refilled with what the caps had left it, no cap laid on them: the caps and the sea now settle together, a cap's margin stands on its own edge, not on the sea's floor or a cliff's foot, and water's caps wear their snow's colour.

How it went:

- **The rain.** Its terms along a row through a rectangle showed the places' own rain, read between them, at none, the cubic's lobe below none beside a place raining 10,000 kg/m² a year.
- **The frozen side.** Every cell of it was under ice when the caps were laid; the sea's refill after them bared what it had covered, and settled, the sea stands 70 m lower still.
- **The high coasts.** From the sea's surface, a plateau ending in cliffs at its coast stood bare its first cells in: the margin stands on its own edge.
- **The shades.** Water's caps wear their snow where it lies; the Earth-like's, the brighter, cover a little more of its poles.
- **The colours.** The sea's colour shades from its deep to its shallow by its depth over 60 m and an ice hides what lies under it by its mass's cover; the shallow water does not show its bed, and a thin sea ice is as bright as a thick one.
- **Known.** The sea's ice's squares on the places at its edge where the guide draws none, the Earth-like's edge along a latitude east of its southern continent: next.

## The fiftieth milestone: the snow on its own rain and the ice at the guide's edge

- **F290, F291, F292.** Looking closer at the maps (2026-11-05), the user saw three patches of clear ground with sharp, square edges on TRAPPIST-1e's southern continent, and blocky sea ice along its eye. The patches were the places' own precipitation and warmest read between them, the snow and the caps' edge squared and kinked at the places' cells: each cell's snow now falls by its own rain over the relief, and the caps' edge is read by a monotone cubic. The squares along the eye were the share under ice standing where the guide drew none: a place stands at the edge now by every place it shares a cell with, and the share under ice draws its ice from half way between the places, as the guided edge does.

How it went:

- **The rain.** Read on its own rain, the snow left bare more than the squares had, the lee of the peninsula's ranges, its rain some 18 kg/m² a year; organic in outline, the user said, its edge the rain shadow's.
- **The squares.** A probe at one found the cell's corner under ice at the edge and its far corner, open, at none, the share under ice standing 0.062 into the ice where the guide drew none. Each place by its own standing took the squares away but left a notch where each had stood; standing at the edge by the cells' far corners too, the eye's rim is one curve, and the user's choice stands: the guide draws the edge.
- **The coast.** The eye's ice ended along a row of places where its sea meets the peninsula's coast, a coast's place under ice standing at the edge by its third of sea; from half way between the places it ends between them.
- **Known.** Where the guide and the places disagree, at a coast whose sea its land's air freezes, the guide's edge and the coast's ice meet along the places' column, a channel of open water a pixel or two wide south of the eye; the places' rain in TRAPPIST-1e's dry south steps at the places in a logarithmic stretch, at a few kg/m² a year; the Earth-like's sea ice's edge along a latitude east of its southern continent, the guide's sunlight zonal.

## The fifty-first milestone: the air's radiation

- **F293.** Before the loop's erosion reads the climate through the body's life, its greenhouse had to hold for a thick cold carbon dioxide atmosphere (the user, 2026-11-06, choosing a band-physics grey law over a full spectral scheme, which would need the clouds the weather brings). The grey law, fit to Venus and Mars, warmed Mars to some 289 K under its thicker air 3.5 Ga ago, where published models find at most 220 to 230 K. The air's outgoing radiation is now its column's: carbon dioxide's band by its wings, saturating, its molecules' collisions, its gas condensing aloft, a moist adiabat where it holds water, by two streams over wavenumber, the other absorbers grey and set by measured bodies; and its albedo its own scattering of its star's light over the albedo beneath, sulfuric clouds and a haze where its air makes them, where the pressure's classes stood. Venus at its measured 737 K, not the grey law's 756, moved its lid's last overturn on the grid its lids mature on: its friction is set again.

How it went:

- **The prototype.** Carbon dioxide's band alone gave Mars 214 K today and early Mars 201 to 207 K; a Venus-set continuum scaling as collisions do warmed early Mars to 351 K: what holds Venus's emission at its 0.1 to 1 bar is its clouds, which early Mars, at the same pressures and temperatures, has none of. The clouds are its sulfur dioxide's.
- **The scattering.** Earth's air reflects 0.080 of the Sun's light, leaving 0.255 beneath its 0.30: Mars's 0.25, so that one albedo beneath and each air's scattering replace the thin and the Earth-like classes; Titan's haze absorbs what its thick air would scatter, Venus's clouds reflect what its air does not.
- **The levels and the moist air.** Evenly in the logarithm of the pressure, the water vapour lay in the lowest layer, 2.3 W/m² off; crowded toward the ground, 0.3. On the dry adiabat the water thinned too fast with height and the Earth-like's 3 km summit stood 19 K under its sea level: on the moist one, within the 10 to 18 K its check holds.
- **The suite.** Nine checks failed at first: the summit, the history's present, the polar sea ice, and six of Venus's tectonics, all reading its lid's last overturn, which its friction, set while Venus stood at 756 K, put at 0.95 Gyr. Its lids matured on a grid: no friction held its last overturn in its craters' record within ten kelvin. A last lid that nearly broke renewing its weakest share (F294) came by degrees within a step, but left a surface renewed where the mantle took up none of its cold, and four more of Venus's stages that read the lid's state after it failed: dropped. The friction set again, 0.0107, Venus last overturned 0.661 Gyr ago, nine of its ten checks as before; the regional lid, a lid of regions each foundering by its own share, is its own work (the user).
- **Early Mars.** Frozen at 207 K through its thicker air, its greenhouse saturated and its scattering grown; the gate's valley networks will not come from a warm and wet mean climate, as the published models find.
- **TRAPPIST-1e.** Its albedo its air's own, 0.257, within THAI's 0.22 to 0.28 where the class's 0.30 stood: its sea opens over 0.72 of its dayside and its mean in THAI's case 252 K against their 232 to 246, THAI's models' brightness lying in a cloud deck over its star.
- **Mars's warm places.** Its warmest 217.8 K where 235.4: the moist adiabat takes its warm cells' air as saturated, where Mars's holds a small share of that.
- **Known.** The column's humidity saturated wherever the body holds water, where a dry body's air holds what its cold traps allow; carbon dioxide's dimers, other bands and its collisions' temperature not in it, early Mars some 10 to 15 K under the 3D models' cloud-free means; the albedo beneath stated, alike under a locked world's star; Venus's last overturn stepping on its lids' grid and three fifths of its ridges' points in its lowlands where four fifths, waiting for the regional lid.

## The fifty-second milestone: the spectrum's ends

- **F295.** The tests' cheap invariants (the user, 2026-11-07: cut the suite's time without its coverage, and add invariants that check what the slow runs only sample) found the column's spectrum short at both ends: a clear air let out 1.32 times sigma T^4 over a 5 K ground and 0.52 times over a 2,000 K one. Its wavenumbers now run from 0 to 25,000 cm^-1, finer below 61, every bin between as it was: within 0.7 percent from 5 to 2,500 K.

How it went:

- **The invariants.** A clear air lets out a black body's light; conservative scattering layers stack as one of their summed depth; an air's molecules scatter each on its own; the outgoing radiation rises with the ground's warmth and falls under more air; every preset's outgoing radiation rises with its warmth, its year's sunlight falls alike on each hemisphere (the time it spends at each distance taking the distance's dimming back), its balance on a small grid makes no heat, and without a tilt its hemispheres mirror each other; the rivers keep the water on rough ground of nested basins under a wet year, a dry one and a mixed one; the sea ice's drift keeps its ice over partial seas under a turning wind. Each was checked against a bug planted in the code it guards, a closed lake's partial cell short, the drift's flux short, the sunlight table's latitudes off by one: each caught it. The caps' laid ice is their drawn thickness by construction, and was left out.
- **The frost floor.** Over a ground colder than its air's frost point at the surface, the column stands on the frost point, warmer than its ground: a thicker air lets out more there (95 percent carbon dioxide over a 162 K ground at 1 bar), and Venus's whole column under some 300 K lets out the same over any ground, the ground hidden under it. Such a ground would take its air as frost: the invariants hold over the ground the air may stand on, and to rounding.
- **The presets.** Their year means within 0.1 K: the Earth-like's poles 256.0 to 255.9 K, Triton's equator 47.9 to 47.8.
- **The suite.** 373 tests passed with the whole-body tests in 14.6 minutes under nextest, the suite's longest started first and the runs its tests shared read once (the tests' tiers); the fast tier, 338 tests, in 2.8 minutes.

## The fifty-third milestone: the regional lid

- **F296.** Venus's rock lid foundered whole on the times it matured, so that its last overturn, and every stage that read it, stepped by some 0.3 Gyr as its friction or its surface's temperature changed (the user, 2026-11-06: a lid of regions, its own work). The lid is now 256 regions whose friction spreads about the lid's, each foundering on its own and the mantle taking up its cold by its share; the surface follows them. The user chose the spread of friction, its width set by Venus's tesserae (2026-11-08).

How it went:

- **The prototype.** Regions under one mantle, judged one by one within a step, cooled it for the next: with no spread Venus's lid broke 36 times where it had 29. Judged at the mantle as it stands, no spread founders whole, exactly as the lid of one piece did. A spread scatters the share each overturn takes rather than its time: regions that formed together mature together, and each overturn cools the mantle under all of them.
- **The anchor.** Venus's tesserae, the share older than its plains, could not be read off its surface's ages by a cut: the largest age class flipped between none and three fifths with the friction. Their craters could: some 1.4 times as dense as on the rest of the planet (Ivanov and Basilevsky 1993), the strongest eighth 1.44 times as old at a spread of 0.02: a lid that founders nearly at once, as the near-random spread of Venus's craters suggests. The friction keeps 0.0107, the surface 0.62 Gyr old.
- **The surface.** Each overturn since the oldest standing region formed floods its own regions; the regions stand on the old crust's smoothed height, ranked once when that region formed. The resurfacing age is the surface's mean age; the loop's craters, floods and coronae run from the oldest standing region; the thrusts and ridges read their own region's shrinking; what predates its ground is not kept.
- **The checks.** Venus's renewal by its regions' ages, its tesserae older than its plains, its thrusts after their ground; the episodic sweep along its radius at steps of an eightieth, where a fortieth moved its degree by a half where the surface renews fastest, and its surface's mean age by under 0.06 Gyr each thousandth.
- **Known.** Every place reads the lid's mean state: the weakest tenth's renewal thins every place's lithosphere with the mean heat flow by up to a fifth for some 30 Myr, a 65 km crust 0.3 Gyr after the overturn keeping 0.976 of a thousand-kilometre relief where all of it was (held at 0.95); each place under its own region's lid is the next step. Venus's wrinkle ridges stand 64 percent below its mean height, where three fifths since F293, the four fifths Bilotti and Suppe 1999 measured still short (known); its thrusts, 1,011 km of them, four fifths 50 to 100 Myr after the overturn that laid their ground.
- **The suite.** 375 tests passed with the whole-body tests; in one run Mars's, the icy moons' and the Moon's, the Earth-like's and Mercury's bodies, started together, ran the machine out of memory and passed alone.

## The fifty-fourth milestone: each place under its own region's lid

- **F297.** The regional lid's places read the lid's mean state (F296): when Venus's weakest tenth renewed, every place's lithosphere thinned with the body's mean heat flow by up to a fifth for some 30 Myr. The user chose this before anything else (2026-11-08). The history keeps each cohort of regions an overturn renewed together, its lid and its surface's heat flow, the cohorts' heat the body's by their shares; each place's crust relaxes, and the mantle's flow carries it, under its own region's lid.

How it went:

- **The heat.** A cohort's surface loses the mantle's share of the heat made less what the mantle keeps, alike under every region, and what its own lid gives up as it thickens; by their shares the cohorts lose the body's heat at every step, to 1e-4. Read at a step's end, a cohort an overturn had just formed seemed to lose heat the step before it formed: the history's steps are read within themselves.
- **The lid.** Through the weak renewal 0.348 Gyr ago the main overturn's lid stays at 127 to 128 km where the body's mean fell from 128 to 95, the renewed tenth's 30 km; the 65 km crust 0.3 Gyr after the main overturn keeps all its relief again, its check back at all of it.
- **Venus.** Its ridges 66.5 percent below its mean height, where 64 (four fifths measured; known); its thrusts 954 km, 86 percent of them 50 to 100 Myr after the overturn that laid their ground.
- **Still the body's.** The flexure of loads and the stress field, which take one elastic thickness for the body; the coronae's lids, drawn before the regions stand on the ground; every place before the oldest standing region formed, where Venus's lid foundered whole.
- **The suite.** 377 tests passed with the whole-body tests in 13.3 minutes, the bodies' pipelines two at a time (a test group) where three at once had run the machine out of memory.

## The fifty-sixth milestone: the tesserae's own lid

- **F299, F300.** The user found the regional lid's Venus (F296, F297) had lost its tesserae's fabric and drawn hard edges about its tesserae and the regions renewed later, where the lid of one piece's had been smooth, and asked why (2026-11-09); then for both fixed, the tesserae whichever way matches the real Venus more.

How it went:

- **The diagnosis.** Venus rendered under the lid of one piece, F296 and F297, the tesserae's strain sets and fabric printed; then F297 with its drag, and with its relaxation, under the body's lid. The hard edges came with F297 and went with its drag alone under the body's lid; the fabric thinned with F296 and went with F297.
- **The edges (F299).** Each place's lid was read at the step's end, and a lid read before it formed reads as just formed: through the 42 Myr before Venus's main overturn, in its step, the regions it renewed had no stiff lid, and the mantle's flow carried all their crust, 16.9 million cell-substeps where 0.59 million, about the tesserae whose lid stood. Each place's lid now has its formation at the step's start and its end; 0.95 million.
- **The blend not made.** The lid's lateral conduction across a region's edge, proposed to the user for the edges, reaches √(κt), 30 to 55 km, over the time the drag works after a renewal: under one coarse cell. The edges were F299's, and went with it.
- **The fabric (F300).** Three causes: one window for all the tesserae, from the mean of their formations to the plains' mean, which missed the oldest tesserae's first 86 Myr and ran into the plains' new lid; each set drawn with the body's mean lid, the plains' new 11.5 km under the set that strains every tessera, its folds too close for the finish's cells at 2,048 pixels; and, under F297, no strain on the oldest tesserae after the overturn. Their window is now their own, 3.573 to 3.859 Gyr as the lid of one piece's was, each set on their own lid, 19 km under that set.
- **The six.** Six of F296's twenty tesserae renewed with the plains about them, no older than they. The user left the choice to whichever matches the real Venus: its tesserae are its oldest unit, embayed by every plain (Ivanov and Head 2011), so the tesserae are the strongest regions that outlasted every overturn since they formed, their share from the lid's history. At the friction F296 kept they stood on 5.5 percent, where 7.3 to 8 are mapped; the friction, a stated choice set loosely by the surface's mean age, is now set by their share: 0.01074, 7.8 percent, its mean age 0.643 Gyr where 0.629.
- **Venus.** Its finish keeps 7.2 percent tessera, the lid of one piece's 7.8, each tessera's fabric over the whole of it and its margins onlapped as before F297. Venus's wrinkle ridges stand 59.5 percent below its mean height, where 66.5 under F297, whose drag had carried the renewed regions' crust before their overturn (61.3 with F299 at the friction before, about the lid of one piece's three fifths; four fifths measured; known); its thrusts 524 km, 83 percent of them 50 to 100 Myr after the overturn that laid their ground.
- **The suite.** 378 tests passed with the whole-body tests in 15.6 minutes. A broken lid keeps tesserae only where its last overturn left some of it: the sweep's check that tesserae stand exactly where the lid broke now runs one way.

## The fifty-seventh milestone: the plains' own spill point

- **F301.** The user asked what the large landforms in Venus's deepest ground were, with a hard edge into a moat about them (2026-11-09), then for them fixed before the ground's own budget.

How it went:

- **The diagnosis.** The landforms are the plains of Venus's latest overturns, its weakest regions renewed again 0.348 and 0.050 Gyr ago: their outlines are those of the youngest surface age exactly, on the lowest of the old ground, as the regions stand by its height. Profiles across their edges: the old plains falling into a trough up to some 350 m deep and 150 to 250 km wide, then a step of 500 to 730 m onto the young plains, which sag inwards from a lip. F296's Venus had the step, some 600 m; F299 had not touched it.
- **The cause.** Each overturn fills the depressions on its plains towards their spill point, and the spill point came from a flood of the whole body from its lowest cell: the basins the late renewals lay in spilled far up the old plains about them, and the lava stopped at the regions' outlines standing all that above the old lid. Its load bent the lid into the moat (F162).
- **The fix.** The flood takes the ground beyond the plains an overturn lays as outlets, each at its own height: the lava spills onto the old lid wherever that stands lower than its surface, never more than its veneer above it. Its load tapers to the veneer at the edge, and the moat goes with the step.
- **Venus.** The young plains' edges run into the old plains' slopes, the young plains some 500 m lower; its heights spread 1,318 m where 1,279, its lowest hundredth below −2,093 m where −1,634; the tesserae's margins and their share, 7.2 percent, as they were; its ridges 60.4 percent below its mean height, its thrusts 524 km.
- **The suite.** 378 tests passed with the whole-body tests in 16.2 minutes; the depression test checks that a liquid spills onto a lower outlet beside it.

## The fifty-eighth milestone: each plains' own fill level

- **F302.** The user saw a cliff and a moat still about the southernmost of Venus's young floods after F301 (2026-11-09).

How it went:

- **The diagnosis.** Profiles pixel by pixel across the flood's west edge: the older plains falling 300 to 500 m into a moat two pixels wide, then a lip as high two pixels inside the young plains; its south edge smooth. The finish lays each overturn's plains from the level the loop filled towards, read between the coarse cells: within a coarse cell of the edge it blended the young plains' spill point, their edge's own ground since F301, with the older plains' far higher one, the young plains standing thicker there and the older thinner.
- **The fix.** Each plains' level read from its own overturn's coarse cells alone. The lip and the moat went; the contact is a change of texture.
- **Venus.** Its ridges 60.4 percent below its mean height, its thrusts 524 km.
- **The suite.** 378 tests passed with the whole-body tests in 15.7 minutes, the climate's work in the tree set aside for the run. The ground's own budget is now the fifty-ninth milestone.

## The fifty-ninth milestone: the ground's own budget

- **F303.** The ground and its air keep their own radiation, the ground taking its air's downward stream and the air its ground's light outside its window, exchanging heat by turbulence alone: the budget Titan's methane needs (the user, 2026-11-08).
- **F304.** The sunlight split by the sun's height: the haze high, water vapour on its slant path, methane and the clouds grey, set by Titan's ground's 10 percent and Earth's 161 W/m².
- **F305.** Penman and Monteith's drying: a wet surface gives its air the latent heat its deficit under saturation draws, at each place's own wind and saturation.
- **F306.** The frost's own gas snows out of air cooled past its point, its latent heat holding the air there.
- **F307.** A body whose water has no liquid to evaporate holds the vapour of its ice, the column at its equilibrium temperature.
- **F308.** The ground a frost leaves bare within a step settles where its sunlight holds it.
- **F309.** The ice behind its edge read along it, as its position is.

How it went:

- **The split alone.** With the sunlight split and the ground and air apart, four whole-body checks failed: the Earth-like's year 291.2 K against its steady 286.5 and its rain 1,800 kg/m² a year, Titan's contrast, TRAPPIST-1e's eye edge and Mars's water, gathered all at its north pole. The Earth-like's ground ran 4 to 5 K over its air, since its latent heat went only by its warmth over an air held saturated; Penman and Monteith's term for the air's deficit brought it to 288.2 K at its own wind, its rain 1,292.
- **The first pass.** The year's first pass runs before its moisture, its air saturated. The user chose to read the checks of a year's mean, its seasons' air and its water from the pass at its own wind, and to mark Titan's shortfalls known: its ground wet everywhere, its unsaturated air cools it to 90.5 K against its steady 94.0, its contrast none. The Earth-like's seasons' mean stands 1.08 K over its steady balance, its ground over the balance's one temperature by what it sheds to its air: the user widened the bound to 1.5 K.
- **Mars's water.** The vapour over the water's ice in its traps stood four orders of magnitude above what it was, by two causes found in turn. Its polar winter air fell to 61 K, the share of its frost's light it takes read at its own warmth coming to nothing: now its carbon dioxide snows out at its frost point (F306); reading that share at the ground's warmth instead ran its summer air away and was dropped. Its south polar summer air stood at 349 K for a step, over its air's column made opaque by an Earth-like vapour: the user chose the frozen body's vapour (F307), and the spike was the frost's last step bare, linearised about the frost's point (F308), which stood before this milestone too (343 K) without reaching the ice. Mars's water now gathers 0.59 north and 0.41 south.
- **The eye's edge.** The user saw the whole body's map, then the edge's close-up: four notches under places whose ice only just held. With the ice behind the edge read along it the edge steps 0.10 degrees. The eye is round, its wind blowing in from both sides alike; the superrotation that shifts it east waits for the weather model (5G).
- **The invariants.** The sunlight's split keeps what a place absorbs; the air's own radiation rises as it warms where its air stands; a bare step settles under its sunlight; a frozen body's air holds its ice's vapour; a wet surface dries under an unsaturated air; Mars's air never stands under its frost point; a thin place behind an edge does not notch it: each failed under its planted bug.
- **The suite.** 384 tests passed with the whole-body tests in 14.8 minutes; the Earth-like's pole's sea ice 3.8 m thick by the end of its winter, inside the 0.5 to 4.4 measured, where 5.4 under F298: its check's bound restored. FINISH_VERSION 224.

## The sixtieth milestone: the liquid apart and the snapshot

- **F310.** The albedo the ground's own colour under its liquid too; the liquid its own map, its depth, with its optics in the sidecar, an engine composing them at any sea level (the user's choice, 2026-11-09; ADR 0007).
- **F311.** The snapshot: the surface at the snapshot's season, the ground under its liquid and the season's ices over it.
- **F312.** The shaded preview lit over its liquid's flat surface.

How it went:

- **What the user saw.** The full map sets of Earth, Venus, Mars and TRAPPIST-1e: a noise texture over the seas in the shaded maps, the albedo without its season's ice, and planet maps whose beaches turn blue. The shaded preview lit the seabed; the albedo drew a submerged cell in the liquid's colour alone, its depth's tint carrying the seabed's relief into the coastal water. A scan of every map's texture over the liquid against its land found these two alone; the climate's maps lie smooth over the sea.
- **The design.** One meaning per map, as the cryosphere's (ADR 0006): the ground, the liquid's depth, the ices; the snapshot their composite at ExoGen's own sea level and season, for viewing and for an engine that wants it baked.
- **The invariants.** The ground shows through its liquid by the shallow-water model, the red going first; the Earth-like's seabed keeps its own colour, its snapshot the water's deep colour over the open deep sea and the albedo itself on its bare land, and the surface its relief is lit over flat at its sea's level: each failed under its planted bug.
- **Mars.** The user's two spots: a lower triangle with straight edges where the great basin's ejecta, the dichotomy's edge and a curved rim meet, and a staircase just east of it where a highland crater drops to the lowlands, are the next milestone; the jagged patches north of them are lava stopped at the dichotomy's wall under ejecta, as the user read them.
- **Seen on the way.** The snapshot shows snow lying on TRAPPIST-1e's continent under its star, which stands cooler than the sea about it with a warm fringe along its coast; and Mars's ground is near grey, some 88, 84 and 81 of 255 in sRGB's red, green and blue, where Mars's is orange-red. Both are noted to look into.
- **The suite.** 385 tests passed with the whole-body tests in 15.3 minutes. FINISH_VERSION 225.

## The sixty-first milestone: the crust continuous across a basin's rim

- **F313.** The crust a basin throws rises across its wall to its rim and thins beyond it: no step at the rim, the dichotomy's edge no stairs, no triangle where two collars overlap.
- **F314.** The mega-basin buries the ground it leaves, giving way to its levelled floor across its wall.

How it went:

- **The user's spots.** On Mars's height map a lower triangle with straight edges in the highlands near 0°E 25°N, and a staircase just east of it where a highland crater drops to the lowlands. The keyframes showed the triangle born between 0.60 and 0.80 Gyr into the body's life, and the dumps after each impact over 500 km showed no impact drew it: it appeared at a step's boundary after the dichotomy's impact, in the lid's relaxation. The crust's map showed why: a step of some 40 km at every basin's rim, a staircase along the loop's cells, and a wedge 130 km thick where the dichotomy's collar and a 3,860 km basin's overlapped, which relaxed its relief fastest.
- **A collar from the rim outward, tried.** Rising beyond the rim it moved the crust's peak a wall's width out; the crust's thickest stands at a basin's rim, so it rises across the wall instead.
- **A trough, found and removed.** Burying the levelled floor's relief inside the wall cut a trough a kilometre deep and more along the dichotomy's foot, where the burial drew the floor towards the old ground's lows; the burial now takes only the ground the basin leaves. A trough some half a kilometre deep remains at 14 to 16°E, as before this milestone.
- **The invariants.** A mega-basin throws its crust continuously across its rim; Mars's crust steps under 30 km between neighbouring cells of the loop's grid, 48 km under the step planted back: each failed under its planted bug.
- **The suite.** 386 tests passed with the whole-body tests in 14.9 minutes. The impacts' VERSION 39.

## The sixty-second milestone: the ice's edge through the cycle

- **F315.** An ice's edge drawn at 96 phases of the cycle as its yearly edge is, its position smoothed along the guide; a point freezes and melts as the edge passes it.
- **F316.** The places about a point count by their share of the ice's ground.
- **F317.** A window open for under a ten-thousandth of the cycle lies all year.

How it went:

- **What the user saw.** On the snapshot, the first map of the ices at one phase: steps in the Earth-like's southern sea ice east and west of its continent, and a lobe off the tip of TRAPPIST-1e's northern continent; and asked how steps fixed before had come back. The yearly edge had not changed: the windows' lengthening toward the ice that lies all year (F255), by the places' share, had never been drawn at a phase.
- **The tries.** A guided edge between two kinds of places, those the ice lies at all year and the rest, took the steps away; the user saw a kink at the eastern end of the southern cap, where places open 1.8 and 3.9 percent of the year, of the other kind than their neighbours, moved the edge a row. Weighing the edges by the openings' contrast bent the edge elsewhere. The snapshot drawn from the ice's mass between the places smudged the edge and opened a gap at the coast (the user). The user chose the proper fix: the windows read between the places with no kinds, nothing committed until it looked right on the Earth-like, Mars and TRAPPIST-1e.
- **A value smoothed, then a position.** The share of the cycle the ice is away, smoothed along the guide as a value, has no kinds and passed its tests; the Earth-like's render showed slits and tongues a row wide, the long openings' tails carried into the ice that lies all year, where at a summer's end an opening of a day draws the edge. What the yearly edge smooths is a position: the edge is now drawn at each of 96 phases and its position smoothed, and a point's window read from its passing.
- **Found in the renders.** Slivers of ice off the northern edge at its winter's end: no edge stood past the outermost places the ice comes to, so the ground past them froze when the last of their neighbours had; the edge now runs on to half way toward the places the ice never lies at, as its window closed on its peak. A synthetic opening nearer the pole than other ice took a bite of that ice: the guided edge is trusted only by the edges on the side of a point its places' own edge may stand on. That left a blank spot where the ice met the southern continent's western coast, which the user saw: those edges count whole, not by their share of ground.
- **The coast's ice on TRAPPIST-1e.** The user found the short straight edge left at the continent's tip natural, and the like of it on the southern continent's western shore, to be left if physical. The cold land holds that strip of sea frozen, 2.5 to 7 t/m² at the cycle's places by the shore; where it ends is the places' column. Left, and recorded as the layout's resolution.
- **The sea's floor.** The user asked whether the ground under a long-standing liquid is reworked: the Earth-like's and TRAPPIST-1e's floors carry craters and fine relief. Nothing does it today, the sea being filled in the finish alone; an item stands in the plan under 5C with the loop's erosion through the history: sediment, the liquid's shielding of its floor from small impactors, its shores' wear; Phase 6's plates renew the Earth-like's floors besides.
- **The invariants.** A window opens from the ice that lies all year to its neighbours' own; the edge at a phase runs along the guide between its places; a place that all but lies all year does not move it; an ice lies the more briefly toward its farthest reach; an edge colder than a point's places leaves their ice; the open sea's edge runs on to a coast; a coast's ice ends short of half way to its open places; an opening too short to write is none: each failed under its planted bug.
- **The suite.** 393 tests passed with the whole-body tests in 14.8 minutes. FINISH_VERSION 226.

## The sixty-third milestone: the surface layer's stability

- **F318.** The heat and vapour a ground exchanges with its air scaled over its neutral by the surface layer's stability, by Monin and Obukhov's similarity with the measured profiles (Dyer and Hicks where the ground is the warmer, Beljaars and Holtslag 1991 where the air is), with no critical Richardson number: over land, sea and frost and in a wet surface's drying, in the cycle and in the moisture, on every body.
- **F319.** The moisture spun up until the body's rain and evaporation over the year agree.

How it went:

- **Why.** Titan's methane loop ran away: 380 kg/m² a year, some 6 W/m² of latent heat where its ground absorbs some 0.35 W/m² of sunlight. Its ground stands under its warmer air, and the neutral exchange drew the air's heat down into the evaporation as freely as a warm ground gives its own. The loop's erosion on Titan waits until its methane's cycle is credible (the user, 2026-11-08).
- **The sources.** The profiles from Paulson 1970, Dyer 1974, Högström 1988 and Beljaars and Holtslag 1991; ECMWF's IFS documentation for their use in a model's surface layer; Grachev and Fairall 1997 for the bulk Richardson number against the linear forms' critical 0.2.
- **Titan.** Its rain 1.0 kg/m² a year, balanced by its evaporation once the moisture ran until it balanced (F319: two years had rained its column's start out, 10 against 3.8); its contrast 1.7 K under air flat at 96.5 K; its land 46 percent arid, 27 hyper-arid. Its rain stands under the 5 to 45 of its GCMs, its ground 1 to 3 K under its air: the column's troposphere on the adiabat leaves its ground radiating 0.46 W/m² net where 0.38 (McKay et al. 1989). The user asked whether a troposphere would have to be defined for every planet: no, a radiative-convective column is one process whose profile follows from each body's own air; it changes little where convection rules the lower air (Earth, Venus, TRAPPIST-1e) and most where the air is radiative (Titan). The user chose to commit the stability as it stands, Titan's low rain known, that column to follow the sea ice's milestone, then the loop's erosion on Titan (2026-11-10).
- **The Earth-like.** Its land at 45°N swings 39.5 K where 43.6, its summer's afternoons 7 K cooler and its winter's nights 3 K colder: under 40, F244's check failed. Probed: with the exchange held neutral it swung 43.6 again. The user restated the check (2026-11-10): over 30 K and over four times its sea's 7.5. The day's replay holds its air at the day's mean, which overstates both its noon's convection and its night's inversion; known.
- **The pole's sea ice.** The whole-body tier failed one check: the Earth-like's sea ice at 85°N 4.6 m by the end of its winter, over the 4.4 measured. Probed at the pole: the neutral exchange had given the ice 45 W/m² of its air's heat each winter, five times SHEBA's strongest month (Persson et al. 2002); under the stability 25 at the stated wind, 12.5 at its own, its surface 16 to 20 K under its air. That surplus had hidden the missing snow on the sea's ice and the polar air's thin winter infrared; the Earth-like's poles stand some 7 K colder in the year's mean, 14 in winter's. And the year at each place's own wind, which the maps use, leaves the pole's ice still thickening when the run stops, before this milestone too. The user chose to commit, the check held under 5,000 kg/m² as known, and make the sea ice the next milestone: its snow, its year run until it settles, its check at the maps' year with its measured bound back; then the radiative-convective column (2026-11-10).
- **The bodies.** At 2048 against the last milestone's renders: the Earth-like's mean 286.4 K where 286.9, its land's 282.0 where 283.8, its 88°N 243.8 K where 251.2 and its winter's coldest there 219.9 where 234.4, its rain 1,254 kg/m² a year where 1,251; Mars within half a kelvin, its frost-held winters as they were; Venus 735.5 K where 735.8; TRAPPIST-1e 253.8 K where 253.9, its rain 507 kg/m² a year where 499; Titan's ground 92.8 K where 88.5 (Huygens 93.7), its runaway evaporation having cooled it, its rain 0.74 kg/m² a year where 253, its swing 3.0 K where 0.7, its equator the driest, aridity 0.005, its poles the wettest.
- **The invariants.** The exchange falls monotonically with the Richardson number, its table read the other way; near neutral it follows the measured profiles' slopes; it dies away by degrees under an inversion, with no critical number; a ground under its air exchanges less than at neutral over its land, its sea and its frost and dries the less, one over it the more, on Earth, Mars and Titan, each failing under its planted bug (the stability dropped from all the sites, from the land's, from the drying's); Titan's methane evaporates within the sunlight its ground absorbs and its year balances, which the exchange held neutral fails (380 kg/m² a year, 6.1 W/m² of latent heat) and the moisture stopped at its spin-up's years fails (its rain 10.1 against its evaporation's 3.8).
- **Known.** Titan's low rain until the radiative-convective column; the Earth-like's pole's sea ice and its poles' winters until the sea ice's milestone; the day's replay at its day's mean air; the finer cells' downscaling at neutral; the wind's drag not corrected.
- **The suite.** 398 tests ran with the whole-body tests in 20.3 minutes, 397 passing: every test some 40 percent slower than the milestone before, those this milestone does not touch among them, the machine loaded from elsewhere: the Earth-like's year alone took some 245 s then, 70 to 90 s since; the one failure the pole's sea ice, held as known with the user and passing since. FINISH_VERSION 227.

## The sixty-fourth milestone: the sea's ice

- **F320.** Snow lies on the seas' ice: it insulates it, melts before it, floods into it past its freeboard and drifts with it.
- **F321.** The year runs until each place's sea ice settles, not only its temperatures; F244's pole check reads the year the maps use.
- **F322.** A place's ice approaching its balance geometrically is carried ahead to it, then run on until it settles.

How it went:

- **Why.** M63's whole-body tier found the Earth-like's pole's ice 4.6 m, over the measured 4.4, and the probe behind it found two gaps: no snow on the ice, and the year the maps use stopping while its ice still thickened, unseen because the check read another year (the user chose to commit M63 and make this the next milestone).
- **A finding on the way.** Maykut and Untersteiner's model found snow under 70 cm changes the ice's balance little: its bright albedo delays the summer's melt as much as its insulation slows the winter's growth. Snow was never going to bring the pole under 4.4 m alone; the larger lever is the polar winter's downward infrared, the radiative-convective column's (told to the user as the milestone began).
- **Probes.** Run its full 40 years, the year at each place's own wind headed for some 7,650 kg/m² at 84°N with no snow, its year's change shrinking by some a tenth a year; the year at the stated wind settles near 3,540. The drift is not the cause: switched off, the ice ran the same to within 15 kg/m².
- **Found by the checks.** The flooding first measured the excess snow against the ice before it flooded, and left the ice under-floated, 24 kg/m² of snow where it floats 33; the unit check caught it, and the flooding now solves the balance with the flooded snow among the ice.
- **TRAPPIST-1e.** The tier found its eye open over 0.83 of its dayside, past the known four fifths. Probed: no snow reaches its ice, no rain or snow falling on its night; run on, its night's ice thickens some 100 kg/m² a year without slowing, toward a balance of some 37 m set by the 2 W/m² its sea gives the ice's base alone, and blankets its sea, so that its day warms by some 0.1 K a year and its eye opens: 0.52 in its first year, 0.76 where the old run stopped, 0.83 by its 80th. The user held it under nine tenths as known until the glaciers, whose sea glaciers would carry the night's ice to its day.
- **The time.** Run plain until it settles, the Earth-like's year at its own wind took 54 years and the tier 25 minutes. The user chose to hasten it (F322). Carried on every place's change, the edge of the ice that lies all year, wandering by 50 to 100 kg/m² from year to year with its drift, was carried on its wander and the year never settled in 80; carried only where its last two shares agree, three times at most, it settles in 31 over windows of four years.
- **The render's year.** The Earth-like's render took 13 minutes where 6: its finish's two tilts ran their 80 years unsettled, one held by a place at 80°S swinging with a period of four years about its balance, the other by the pole's thick ice thickening under its drift's noise. Settled over windows of four years, the swing settles and the year ends in 38 to 42 years, its finish 355 s.
- **The pole.** In the maps' year, settled by 54 years, its ice 7,014 kg/m², 7.6 m, under 54.5 kg/m² of snow, where the stated wind's year 4,136. The user held the check under 8,000 kg/m² as known until the radiative-convective column, and accepted the minute it adds to the Earth-like's year (2026-11-10).
- **The bodies.** At 2048 against the last milestone's renders: the Earth-like's sea ice at the snapshot over 7.5 percent of the body where 8.7, its mean north of 60° 2,354 kg/m² of the sea where 3,180 and south of 60° 586 where 1,204, the snow thinning it but where the pole's thickest settles thicker; its mean 286.8 K where 286.4. TRAPPIST-1e's night's ice thickens through its 80 years, north of 60° 3,091 kg/m² where 2,124, its mean 254.0 K where 253.8. The Earth-like's finish takes 653 s where 341, its year's tilts running 38 to 42 years where some 16 to 18; TRAPPIST-1e's 271 s.
- **The invariants.** Sea ice swinging about its balance with a period of four years or two settles over windows of years, ice drifting does not, which year-by-year settling fails; snow on a sea's ice insulates it, its top the colder; it melts before its ice; it floods into a thin ice until the ice floats it at a freeboard of none, the ice and snow together kept; none lies on open water; it drifts with its ice, a tenth of each place's staying a tenth: each failed under its planted bug (the snow left out of the conduction, the ice melting first, no flooding, snow on open water, the drift leaving it behind). The maps' year ends with its ice settled, which a year stopping on its temperatures alone fails, after 23 years, its ice unsettled. The hastened ice ends where the plain ends, in fewer years: carried the wrong way it took 37 years where 32, and not carried, as many.
- **Known.** The pole's ice 7.6 m until the column; TRAPPIST-1e's night's ice and its eye until the glaciers; the snow's own albedo; snow-ice without the sea water it takes in.
- **The suite.** 402 tests passed with the whole-body tests in 22.9 minutes, where 14.7 unloaded at M63: the years the sea ice takes to settle in the sea worlds' checks, which load the machine for the rest. FINISH_VERSION 228.

## The sixty-fifth milestone: the radiative-convective column

- **F323.** Each column's profile stands in radiative-convective equilibrium, solved by Newton's steps at nodes of its ground's temperature, kept for every column alike.
- **F324.** Carbon dioxide's, water's and sulfur dioxide's lines from HITRAN, carbon dioxide's collisions and continuum, its heat capacity by its temperature; Venus's measured water and sulfur dioxide.
- **F325.** Each column takes the sunlight its planet absorbs by its absorbers' depths; the haze's heat returned, a stated share, is gone.
- **F326.** The sulfuric deck takes its sunlight in its upper cloud and its infrared depth evenly in the pressure's logarithm, mixed in its mass in a thick air; it stands in the air over a summit; its depth set by Venus's 737 K.
- **F327.** The temperature history solves its column at five of the star's points.
- **F328.** A preset's graph is computed once under its values.
- **F329.** Earth's clouds' law set again on its anchors.
- **F330.** A thick air stands on F293's profile, the adiabat up to its skin, by degrees over its infrared depth.

How it went:

- **Why.** The user chose the radiative-convective column next (2026-11-10): the column stood on the adiabat up to its skin, and the polar winter's downward infrared, Titan's ground and the deep atmospheres of Phase 6G want the air's own profile. The user's plan, all in this milestone: the sunlight laid by its absorbers, the haze's thermal opacity, the haze's stated return removed, the tables from the nodes' profiles, and the anchors checked again; committed only when Titan and Venus stood near their measurements again.
- **Venus's windows.** The first column let Venus's deep air out through its near-infrared windows: 700 K. Carbon dioxide's lines and far wings, its collisions over five spans, Venus's measured water and sulfur dioxide, and the continuum its windows measure closed them (F324).
- **The time.** The first fast tier timed out, a dozen tests past 5 minutes, where it took 2.8 minutes: each column's nodes solved afresh for its mean, its tables and its heights, a graph's 49 history points each a column, the defaults' graph solved before every preset's, and each relaxation step some 2 ms of 328 bins and 47 layers. Shared nodes (F323), five history points (F327), the graph computed once (F328), Planck's function once a level and the line tables read by their logarithms brought it back.
- **The solve.** A probe of Venus's hot nodes found the relaxation stopping some 13 K from their equilibrium (F323). Newton's steps with the analytic slopes cycled where a level stood exactly at its frost point or on its adiabat, held by neither; held by both, and freed where its own radiation would warm it, they converge in 2 to 8 steps a node over most of the bodies. Venus's coldest node, at its equilibrium temperature under 92 bar, lay under its air's frost point and cost each column 3,000 steps; its chain begins at its frost point now, and a solve ends after 400 steps, its balance the same to the hundredth of a kelvin.
- **The deck.** Converged, Venus's ground stood at 693 K. The deck spread by its air's mass held its depth low; by the pressure's logarithm, as Venus's measured clouds stand, 725. Its outgoing radiation then rose and fell with the ground: a deck placed where its acid condenses on the ground's adiabat rises as the ground warms, its top the colder, and Venus's balance stood anywhere from 682 to 725 K with the nodes' spacing. Convection as a smooth flux, to let Newton's steps reach every node, was tried and not kept: its steps took 30 to 200 rounds and its balances moved with its ramp. The user asked whether a coarser model holding the bulk of the airs' effects would do (2026-11-12). The deck's infrared depth spread through the whole air and set by Venus's measured 737 K held its air above its clouds at some 580 K: not kept. The user chose the placed deck, Venus known at 725 K.
- **The sweep.** The whole-body tier then found Venus's ground leaping between some 696 and 726 K from one body to the next along a sweep of its radius (principle 15); a deck held where Venus's measured clouds stand leapt too, 642 to 660 K: the solved profile is no smooth function of the ground's temperature in an air so thick and so near its adiabat. The user chose the blend by the air's infrared depth (F330): Venus on F293's profile, its deck mixed in its mass and set by 737 K, its sweep smooth.
- **The summit.** Venus's 11 km summit stood under more greenhouse than its lowlands, 688 K: its deck was placed from the summit's own ground, and rose over every height. It stands in the air now, 661 K.
- **Titan.** Its ground at 97.3 K where 93.65 (F325): held as known, the user's choice for a coarser model; until methane's moist convection and its sunlight's bands.
- **Earth.** The column left Earth's anchors at 292.2 K; the clouds' law holds them again (F329). The defaults' Earth, with no ocean, stands at 294.8 K, its drier air on a steeper lapse; its test says so.
- **The invariants.** The net flux's slopes meet its differences over a dry carbon dioxide air within a fifth, which their downward stream's sign turned fails; columns alike share their nodes and a column under other sunlight keeps its own, which a key without its absorbed sunlight fails. The outgoing radiation rises with its ground over a ground its air may stand on, at or over its frost point, as the thinner air's check already was; a frozen air's profile is held at its frost. Carbon dioxide's forcing is checked against the band's own law over the column's own stratosphere.
- **The pole.** The Earth-like's pole's sea ice, in the year the maps use, 3,123 kg/m² by 21 years under 60 kg/m² of snow, some 3.4 m, where 7.6 (M64): within the Arctic's measured 0.5 to 4.4 m, as M64 expected of the column; its check is the measured range again. The hastened ice settles as fast plain, both in 33 years, its check now no more years than the plain's (a carry the wrong way took 37 where 32).
- **Known.** Titan's ground 3.7 K warm, and its methane's evaporation 0.455 W/m² where its ground's sunlight 0.35, its warm troposphere sending its ground the more (F325); Venus's air on F293's profile, without its deck's own heating (F330); the stratospheres unwarmed by ozone. Not done: the profile's table in the sidecar and the column with no ground assumed (PLAN.md), the polar winter's checks with the clouds' milestone.
- **The bodies.** At 2048 against the last milestone's renders: the Earth-like's mean 288.3 K where 286.8, its sea ice at the snapshot over 6.2 percent of the body where 7.5, its mean north of 60° 1,121 kg/m² of the sea where 2,354 and south of 60° 114 where 586, its thickest 5,817 where 18,069; its finish 387 s where 653, its ice settling sooner. TRAPPIST-1e's mean 253.4 K where 254.0, its sea ice over 59.6 percent where 59.1; its finish 268 s.
- **The suite.** 404 tests passed with the whole-body tests in 18.0 minutes, where 22.9 at M64; the fast tier in 3.1, where 2.8 before the column. FINISH_VERSION 229.

## The sixty-sixth milestone: the clouds

- **F331.** The water's clouds by place and season from the year's own air: convective cloud by its rain (Slingo 1987, as ECMWF ran it), its towers and anvils; layer cloud by the large-scale rain and the column's saturation; low cloud by the inversion's strength (Wood and Bretherton 2006) where the air under it is moist.
- **F332.** The inversion's strength under a free air that stands at the overturning cell's warmest within its reach, the tropics' weak temperature gradients.
- **F333.** Each kind's light by its optical depth and its drops' asymmetry, its infrared by its emissivity from its top to space and from its base to the ground, overlapping at random.
- **F334.** The clouds a diagnostic, exported as `clouds.exr` and `clouds.png` and the sidecar's `clouds` block; their coupling into the year tried and not kept.

How it went:

- **Why.** The user chose the clouds as M65 began, and split them in two (2026-11-12): M66 the clouds by place and season with the planet's means kept, M67 the means from the clouds themselves; the other condensables' decks later.
- **The scheme.** The first version found no layer cloud, a column's mean saturation never reaching the 0.8 a model layer's humidity reaches under a cloud, and anvils and infrared far past Earth's. The layers now come from the large-scale rain's systems by the convective law and from the saturation over 0.7; the anvils begin at C = 0.3 without ECMWF's step at 0.4. The rain is the convection's within the overturning cell's reach and the large-scale's beyond, by degrees.
- **The free air.** The subtropics' low cloud needs their cool sea under the deep tropics' warm free air: the free air stands at the warmest of the cell's latitudes within its reach, its own latitude's beyond, by degrees over a quarter of the reach.
- **The coupling.** Held to the planet's means, the clouds' departures from their means alone cooled the Earth-like's poles 7.5 K and its mean 3.0: the stated law's greenhouse already falls toward the dry poles, and the clouds' fell there again. Taken beyond the stated law at each place, the year stood 0.9 K colder, its storm tracks 3.3, its sea's ice over 0.136 of it where 0.100; the whole-body tier then found its rain at 1,491 kg/m² a year where its check allows 1,400 and its sea at 85° S open some of the year. With the law's brightness held off the ice, as its ice's albedo replaces its ground's, its poles froze instead: 286.0 K, the ice unsettled after 80 years. A probe of the parts found the infrared to space did it alone: the deep convection, under the near-aquaplanet's 18 mm a day of rain where Earth's is 6 to 7, keeps 67 W/m² beyond the law where Earth's clouds keep some 10, and with the mean held the poles lose some 23. The user chose the clouds as a diagnostic now and their coupling in M67, after the deep convection's high cloud is brought to the amount measured (2026-11-13).
- **The Earth-like's clouds.** On the year the maps use: cover 0.497 (some 0.67 measured), 0.59 over its deep convection, 0.29 over its subtropics, 0.61 in its storm tracks; their albedo 0.189 over its ground by area, 41 W/m² kept from space where CERES measures 28, 31.6 sent down where some 27. Its deserts are not seen: the test Earth's places are all 71 percent sea.
- **The locked world's.** Under its star 0.77, its dayside 0.50, its night 0.03 (Yang et al. 2013).
- **The invariants.** Every kind covers between none and all of its place, which a low cloud unclamped fails; each grows with its cause, which a convective law turned fails; on the Earth-like's sky a high anvil keeps more from space than a low deck and sends less down, which the anvils' tops read at the low deck's fail, and none change nothing; its belts and the locked world's clouds under its star, which both fail with the rain taken from the scheme.
- **Known.** The deep convection's high cloud, some twice Earth's in the infrared; the cover 0.50 where 0.67.
- **TRAPPIST-1e's night.** Its clouds cover 0.15 of it: under its star, and over its night's land, whose columns stand at 0.92 to 0.98 of saturation at some 230 K, its snow keeping them so, where its frozen sea's stand at 0.1: the land under full low cloud and layers. Such a column holds some 0.15 kg/m² of vapour, where a deck of 70 g/m² cannot stand: known until the clouds' water follows their column's (M67). The cold air's own fix in a climate model's scheme (Vavrus and Waliser 2008), cloud by the air's water under 3 g/kg, would clear the Earth-like's polar winter as well, whose cloud stands near the Arctic's measured.
- **The bodies.** At 2048 against the last milestone's renders the climate is the same, the Earth-like's mean 288.3 K, TRAPPIST-1e's 253.4. The Earth-like's clouds lie over 0.50 of it in its belts: the deep convection's band, its subtropical land clear, its storm tracks and its poles cloudy, its tropical land under its convection's towers and anvils. Its finish 428 s where 387; TRAPPIST-1e's 294 where 268.
- **The suite.** 407 tests passed with the whole-body tests in 18.8 minutes, where 18.0 at M65; the fast tier in 3.2. FINISH_VERSION 230.

## The sixty-seventh milestone: the clouds' physics

- **F335.** The clouds' water held to their air: the decks' the adiabatic water of their measured depth and adiabaticity from their base, the towers' a share of what their ascent condenses, liquid by degrees from 241 to 271 K.
- **F336.** The convection's high cloud by its climatological rain, as ISCCP-H and GPCP measure it over the tropical oceans: all of it 0.054 + 0.185 ln P, its towers 0.0095 P, its anvils thickening with the rain.
- **F337.** The anvils from ISCCP's 440 hPa to the tropopause, letting out from the middle of the span.
- **F338.** The clouds' own absorption under the vapour over them; ozone by the air's oxygen and Lacis and Hansen's fits.
- **F339.** The year under its clouds, built and checked: its own ground's albedos, its clouds' deck at each step, the light's path through its absorbers split at the clouds' tops, the clear air's other absorbers set by its measured 73 W/m².
- **F340.** The year's own planetary albedo, in the sidecar.
- **F341.** The coupling waits on the tropics' heat (the user, 2026-11-13).

How it went:

- **The research.** Two searches (notes in `research/research_notes`, untracked). The first computed, from public climatologies (CERES EBAF Ed4.1, GPCP v3.2, ISCCP-H, CALIPSO-GOCCP), the clouds' radiative effects and high cloud against the seasons' mean rain over the tropical oceans, CERES's zonal means and CMIP6's aquaplanets; the second the adiabatic water's rate, the clouds' phase, their absorption, the clear air's measured 73 W/m² and ozone by oxygen. Slingo's law was fitted to a model's rain at its steps; the climate's rain is the seasons' mean, which the measured relations are made on.
- **The first coupling.** With the clouds' stated paths and Slingo's law, the Earth-like's year under its clouds stood at 289.9 K, its rain 1,616 kg/m² a year where its check allows 1,400 and its pole's ice 6,263 kg/m²: the stated sunlight absorber gone, its ground took what its air had.
- **The high cloud.** The measured high cloud and towers took the tropics' infrared from 100 W/m² to some 69, which is right for the near-aquaplanet's 10 mm a day of rain at its tropics, and left their light at half the measured: the anvils thicken with the rain. With them the Earth-like stood at 288.2 K, but its pole's ice did not settle: its clouds over its sea's ice reflected as if no air lay between them, 0.70 where CERES's Arctic is some 0.55.
- **The light's path.** The air's share of what a place absorbs after it reflects gave Earth's clear air 53 W/m² of the 73 measured and a deck less; every absorber over the clouds gave 73 but ran away to 319 K, the vapour darkening the clouds' own light as it warmed; the absorbers beneath the clouds, split at their tops by the vapour's share above them, keep the clouds' light and Earth's 73 with the clear air's other absorbers.
- **The coupling's state.** Under its clouds the Earth-like reflects 0.293 of its light, CERES's 0.29, but stands at 296.2 K: its clouds keep some 37 W/m² from space where CERES measures 26 to 28. The tests' Earth-like, the sea over 0.71 of every place, stands at 308 K at its equator under the stated law and its clouds alike; the finish's Earth at 301.4 K under the stated law and 307.2 under its clouds (a trial render, not kept: its mean 294.7 K, its rain 1,647 kg/m² a year where 1,295, its clouds keeping 37.1 W/m²). Under their clouds the tropics warm, their deep convection keeping the more infrared over the warmer ground; the polar clouds, with no inversion under them, the more over their cold ground. The user chose the clouds' physics now, the year on the stated law, and the tropics' heat under their clouds and the polar inversion next, then the coupling (2026-11-13; asked on the tests' Earth-like's 308 K, before the finish's 301.4 was read).
- **The locked world.** Under its clouds TRAPPIST-1e's aquaplanet reflects 0.265 of its light where its first year's air alone 0.255 (THAI 0.22 to 0.28), its dayside 274.9 K, 0.66 of it open where 0.82 (THAI almost half); its mean 259.0 K, 4 over its first year's, over THAI's 232 to 246 as its night's ice blankets its sea (KNOWN, M64).
- **The diagnostic.** On the year the maps use, the Earth-like's clouds cover 0.465 of it, 0.49 over its deep convection, 0.26 over its subtropics and 0.61 in its storm tracks; they keep 31.6 W/m² from space where M66's kept 41.1 (CERES 26 to 28), send 29.0 down (some 27), and its poles' 8 (CERES 6 to 14). TRAPPIST-1e's night's land, its columns at 230 K, holds some 3 g/m² where 70.
- **The invariants.** A deck reflects none of the light where it covers none, the more the thicker and a thick deck the more as the sun sinks, which the sun's height term turned fails; a place under clouds keeps all it does not reflect and its ground what they send down, which the ground without their radiation fails; the adiabatic water as Albrecht's rate and the notes' table give it, which the pressure's fall dropped fails; ozone's share and column, which its ultraviolet bands dropped fail; Earth's clear air takes its measured 73 W/m², which the clear air's other absorbers dropped fails.
- **The bodies.** At 2048 against the last milestone's renders the climate is the same, the Earth-like's mean 288.3 K, TRAPPIST-1e's 253.4. The Earth-like's clouds cover 0.474 of it (low cloud 0.242, layers 0.141, anvils 0.130, towers 0.024), their albedo 0.154 over its ground under its air, 32.1 W/m² kept from space where M66's 40.9 and 29.6 sent down; its year's planetary albedo 0.308. TRAPPIST-1e's cover 0.144, their albedo 0.032, 9.2 W/m² kept and 6.5 sent down where 13.0 and 9.9: its night's land's decks hold little water now, but the map's cover still counts them (each kind's optical depth with the weather's export, 5G); its planetary albedo 0.272.
- **The suite.** 412 tests passed with the whole-body tests in 20.0 minutes, where 18.8 at M66; the fast tier, 375, in 3.2. FINISH_VERSION 231.

## The sixty-eighth milestone: the clouds' own year and the convection's clouds

- **F342.** Each year under the clouds of the year before it, started from where it ended, until the body's mean settles to 0.3 K.
- **F343.** The anvils' optical depths ISCCP's albedo-equivalent means: thin cirrus 1.0, medium-thick 8.5.
- **F344.** The convection's own mid-level cloud, as much as ISCCP sees topped between 440 and 680 hPa against the rain.
- **F345.** The convection's own cumulus, 0.33 of the place where it rains, one boundary layer's cloud with the inversion's.
- **F346.** The clear sky over the warm tropics lets out too much; its own moisture next (M69), the coupling after (M70).

How it went:

- **The diagnosis.** A temporary check binned the Earth-like's tropics by their rain as CERES and GPCP measure them (the M67 research notes' table) and by latitude against CERES's zonal means. Under the stated law its clouds kept 10 to 16 W/m² more from space at each rain and reflected 5 to 8 less: some 15 W/m² too warming at every rain from 2 to 10 mm a day; its deep tropics had no low or mid-level cloud, where ISCCP sees some 12 and 16 percent. The coupled year had run under its first year's clouds: that year rained 2,016 kg/m² a year and its clouds kept 42 W/m² from space. Run under the clouds of the year before it, the year came in three years from 296.2 K to 292.8 and settled. The air's transport was not at fault: the heat was in the clouds and the clear sky.
- **The research.** Two passes (notes in `research/research_notes`, untracked). The first computed from public climatologies ISCCP-H's nine classes against GPCP's rain over the tropical oceans, their albedo-equivalent optical depths and cloud tops, ERA-Interim's temperatures there, and CERES's clear sky against HadISST's sea temperature and ERA-Interim's vapour, all reproducing the M67 notes' table exactly. The second gathered the clouds' measured inhomogeneity (Shonk et al. 2010, Oreopoulos and Cahalan 2005), the trade cumulus and the mid-level clouds, and each cloud type's effect from published radiative kernels over ISCCP-H; several classic papers could not be read and their numbers stand unverified.
- **The closure.** ISCCP's own clouds, each class at its albedo-equivalent depth over the sea under the model's optics through the day, reflect 1.34 to 1.41 times what CERES measures by the rain; at the measured inhomogeneity's 0.7 on their depth, 1.13 to 1.15. A regression of CERES on the nine classes gave impossible coefficients, the classes moving together.
- **The convection's clouds.** With the anvils at ISCCP's depths and the convection's own mid-level cloud and cumulus, the Earth-like's clouds under the stated law change what the tropics let out by the rain within 3 to 5 W/m² of CERES (at 9 mm a day -24.4 where -19.4), slightly too cooling where M67's were 15 too warming; but they reflect 99 W/m² where 71 and keep 75 where 51. Over the body they cover 0.600 of it, keep 37.3 W/m² from space and reflect 60.4 W/m² (CERES +26 to 28 and -45).
- **The clear sky.** The model's clear sky lets out 294 W/m² at 300 K, CERES's 290 to 292, and 303 at 303 K where 283; held at 300 K's vapour it rises 4.5 W/m² a kelvin, and it lets out 24 W/m² less a unit of the logarithm of its vapour where CERES's 36. Its vapour is one law of the ground's temperature (53 kg/m² at 300 K, Earth's columns 40): Earth's rising columns hold the more, its sinking the less, and over 300 K the more of its warm sea lies under rising air. Every cloud hides the clear sky, and the low and mid-level clouds keep some twice what the kernels give a unit of their cover.
- **The coupling's state.** Under its own clouds the Earth-like reflects 0.297 of its light (CERES 0.291) and its clouds take 23 W/m² from it in all (CERES 18), but it stands at 293.3 K: its clouds keep 37 W/m² from space where 26 to 28, 80 over its deep tropics where 35. TRAPPIST-1e under its own clouds, which settle in four years, reflects 0.247 of its light (THAI 0.22 to 0.28), its dayside 273.7 K (THAI 261 to 268), 0.61 of it open (THAI almost half): M67's check that it reflect more than its first year under the stated law, 0.255, is gone, that law no measure of its clouds. Without F343 to F345 its own clouds had settled at 0.187. The user chose (2026-11-13) to land this, the climate still on the stated law, and the clear sky's own moisture with the clouds' inhomogeneity next.
- **The invariants.** A year under its own clouds is under its own: its own clouds stand within some thousandths of the albedo and a W/m² or two of the infrared of those it ran under, over the body and over its tropics, which the clouds laid once and never again fail (the first year's clouds laid once and never again, as M67 ran the year, 0.029 and 11.7 W/m² apart, 0.061 and 19.0 over the tropics; the stated year's clouds laid once and never again come within 0.06 K of the year under its own and pass, as they should); the convection's own clouds show as ISCCP sees them, its mid-level cloud within some thousandths of ISCCP's fit and its cumulus within some hundredths, and each level's cloud keeps its measured optical depth, which the mid-level cloud not hidden under the high cloud and the stated depths' water read as liquid alone fail.
- **The bodies.** At 2048 the climate is the same, still on the stated law: the Earth-like's mean 288.3 K, TRAPPIST-1e's 253.4. The Earth-like's clouds cover 0.597 of it where M67's 0.474 (low cloud 0.351, layers 0.262, anvils 0.130, towers 0.024), the cumulus over its tropics and subtropics, its subtropical deserts clear; their albedo 0.178 over its ground under its air where 0.154, 37.2 W/m² kept from space where 32.1, 36.4 sent down where 29.6. TRAPPIST-1e's cover 0.222 where 0.144 (low cloud 0.144, layers 0.099), their albedo 0.054, 11.8 W/m² kept and 13.6 sent down.
- **The suite.** 413 tests passed with the whole-body tests in 22.5 minutes, where 20.0 at M67, the locked ocean world's year now the longest (its night's ice runs each of its years under its clouds to their most); the fast tier, 376, in 2.9. FINISH_VERSION 232.

## The sixty-ninth milestone: the column's Earth profile

- **F347.** The column lifts its surface air at 0.8 of saturation, its troposphere ERA-Interim's over the wet tropics.
- **F348.** Ozone in the column by its air's oxygen: its sunlight heats its stratosphere, its 9.6 um band absorbs by Malkmus bands.
- **F349.** Each place's clear sky at its own column's vapour, the moisture's.
- **F350.** The clouds read their air's own profile.
- **F351.** The stated law set again on the new column.
- **F352.** The coupling's state; water's lines from HITRAN tried, not kept.

How it went:

- **The vapour first.** The plan was to give each place's clear sky its own moisture. The moisture's columns stood 10 to 20 percent over ERA-Interim's where the law's 25 to 35, its saturation rising with the sea's warmth as Earth's rising columns' does; but drier columns let out the more, and the warm tropics' clear sky, 303 W/m² at 303 K where CERES's 283, could not come down so.
- **The bands.** Water's main bands from HITRAN's lines' medians from 0 to 2,500 cm^-1, the way its other gases' are read, in place of Williams's fit (whose flat rotation band under-reads its peak): the outgoing radiation moved by under a W/m². The bands were not the cause.
- **The profile.** The column's troposphere stood 5 K warm at 500 hPa and 10 K at 300 over ERA-Interim's wet tropics: it lifted saturated air from its ground. Lifting its surface air at 0.8 of saturation brought it within a kelvin; its warm tropics' clear sky came down 10 to 12 W/m², rising 2.2 W/m² a kelvin where 3.0. Its steam airs then rose dry through most of their mass, read against the half of the pressure their vapour is held to, and held Earth's whole table at 271 W/m²; their saturation is now read without it.
- **The stratosphere.** The column's stratosphere stood at 150 to 170 K, Manabe and Strickler's air without ozone; held at Earth's standard stratosphere its outgoing radiation rose 10.5 W/m² at 288 K. M65's column had matched CERES's clear sky by two errors, a troposphere too warm and a stratosphere too cold. The user chose (2026-11-13) to fix the column whole: ozone by its oxygen, from WACCM6's layers, its sunlight by Lacis and Hansen's fits and its 9.6 um band by Malkmus bands the research fitted to HITRAN's lines (the column's median tables left the stratosphere's ozone nearly clear). Its stratosphere now stands at 235 K at 15 hPa; its upper stratosphere too warm, carbon dioxide's line cores clear at low pressure (KNOWN, Phase 6G).
- **The clouds' air.** The clouds read their air down a stated lapse rate, 5 K cold at 553 hPa over the warm tropics; they now read their column's own adiabat, and keep 59.9 W/m² from space at 9 mm a day where 67.9 (CERES 51.2).
- **The clouds' tops.** On the air's own profile a thin anvil kept 60 W/m² from space against a low deck's 21, under what the check of M66 asks; its thin cirrus now let out from the tropopause, as ISCCP's do, and keep 88.8, the kernels' 88.
- **The coupling's state.** Under the stated law the Earth-like stands at 287.5 K; under its own clouds at 291.9 K where 293.3 at M68, the rest of its warmth its clear sky's albedo, its clear sky's heat and its clouds' (F352), M70's.
- **Two checks set to their measures.** The Earth-like's pole's sea ice under the stated law, 2 K cooler, stands at 4,021 kg/m², 4.39 m: its check now holds it to its measured 4.4 m at the sea ice's 917 kg/m³, where 4,000 kg/m² stood a rounding inside. The hastened sea ice took 26 years where the plain 24, ending where it does: held to half a settling window more, where a carry the wrong way took five.
- **The invariants.** The column lifts its surface air, its troposphere within a few kelvin of ERA-Interim's and its table rising through 303 K, which saturated air from the ground and the steam air's capped vapour fail; ozone sits where WACCM6 puts it, which its centre's power turned fails; its band lets through what its lines do, which its path without the diffusivity fails; it warms the stratosphere and takes its share of the sunlight high, which its sunlight dropped fails; the clouds read their air's own profile, which the stated lapse rate fails; the Earth-like's clear sky follows its sea's warmth as CERES measures it, which the law's vapour in place of each place's own fails (272.7 W/m² over a 294.5 K sea where 286.3).
- **The bodies.** At 2048 under the stated law the finish's Earth stands at 286.1 K where 288.3 at M68, its equator at 299.7 where 301.4 (Earth's some 300), its poles at 247.4 where 250.8, its sea's ice over 0.070 of it where 0.059: its places' own columns drier than the law's, which the stated law, set on the zero-dimensional balance at the law's vapour, does not see (KNOWN, until the coupling). Its clouds cover 0.598 of it (low cloud 0.363, layers 0.255, anvils 0.125, towers 0.023), keep 34.4 W/m² from space where 37.2 and send 38.2 down. TRAPPIST-1e stands at 251.7 K where 253.4; under its own clouds its planetary albedo 0.237 (THAI 0.22 to 0.28), its dayside 271.8 K (261 to 268), 0.57 of it open (almost half); its years' means swing a kelvin from one to the next as its night's ice never settles (M64), and run to CLOUD_PASSES's six.
- **The suite.** 418 tests passed with the whole-body tests in 34.5 minutes, where 22.5 at M68, the locked ocean world's year the longest, its six years under its clouds; the fast tier, 381, in 2.8. FINISH_VERSION 233.

## The fifty-fifth milestone: the air's spectrum

- **F298.** The ground's own budget (the user's choice, 2026-11-08, for Titan's methane cycle) needs what the air sends down to its ground, which the grey column (F293) got wrong: Earth's ground took 207 W/m² where 342 is measured, Titan's saw space through 0.30 of its radiation where some 0.04 is real. The user chose the published bands first. The column now gives its outgoing, window and downward streams, its layers' source linear in their depth; water vapour absorbs by Williams 2026's bands and continuum, the gases' collisions by HITRAN's tables, the clouds and the air's other absorbers stated and grey until 5G, Titan's haze sending a share of the light it absorbs back to space.

How it went:

- **The research.** Two searches: water's bands (Williams 2026, the paper behind the column's carbon dioxide band, gives water's from the same fit; Koll et al. 2023 and Cohen and Pincus 2025 beside it, their continua a factor of ten apart) and Titan's collisions and budget (HITRAN's tables, which on the Huygens profile give the ground a net 0.35 to 0.41 W/m² against McKay et al. 1989's 0.379). Notes in research/research_notes/Air spectrum band models.
- **Earth.** Its clear air lets out 268.5 W/m² at 288 K against CERES's 267. Clouds lying as the vapour lies kept its 291.6 K but sent down 390 W/m² and took its slope to 1.3 W/m² per kelvin. Clouds well mixed and filling as the vapour grows kept 291.6 K and 2.1 and sent down 341.5 W/m² where 342 is measured, but thinned to nothing over the poles: the whole-body tier found the Earth-like at 280.8 K, its Arctic sea iced all year and the water's ice in its tropics. Clouds as a power of the vapour keep both anchors and the 341.5 W/m², and a polar air's kilogram of vapour keeps some: the year's balance runs 300.4 K at the equator to 271.1 at 55 degrees and 262.4 at 65 (Earth's some 299, 275 and 266).
- **Titan.** Its window 3.4 percent, its air sending down 3.92 W/m², its ground losing 0.46 net. Under all the light it absorbs it stood at 111 K: its troposphere on the adiabat (known) makes its greenhouse 29 K where McKay's is 21. Its haze's anti-greenhouse, 0.703 of the light its stratosphere takes sent back to space, inside the half to all a stratosphere sends, keeps it at 93.7 K; its year at 94.0 K, its methane rain 811 kg/m² a year where 1,172.
- **The year.** The Earth-like at 284.3 K where 285.9, its sea at 65 degrees north freezing and melting with the seasons, its pole's ice 5.4 m by the end of its winter where 0.5 to 4.4 are measured (known, held under 6.5 m until the ground's own budget); TRAPPIST-1e at 254 K, its sea open over 0.77 of its dayside (known, held under four fifths until 5G); Mars at 216.5 K where 215.5.
- **The plateau.** The balance's invariant that the outgoing radiation rises with the warmth found water's runaway plateau stalling it within 0.1 percent at some 355 K on Earth (and on Titan, at temperatures it never sees, through its ice's inventory); held at its highest past it.
- **The suite.** 378 tests passed with the whole-body tests in 12.5 minutes.
