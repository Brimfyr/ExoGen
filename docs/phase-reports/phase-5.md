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
