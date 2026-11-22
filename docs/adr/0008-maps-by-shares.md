# ADR 0008: The maps by shares, not cells

Status: accepted, 2026-11-22. Supersedes ADR 0007's decisions 2 (the liquid's depth), 4 (the shaded surface) and 5 (`mask_liquid`) as decisions 2 and 4 below state.

## Context

The user looked at the Earth-like's drafts (2026-11-21) and found its lakes, its coasts, its frost and its snow drawn as speckled single-pixel dots, the north-east continent's most of all, and asked whether the physics was limited by the map's resolution or whether some natural anti-aliasing was missing.

Each of those maps drew a threshold cell by cell, all or none: a cell under the liquid where its ground stood below the level (F26), a lake where its pit's water stood, ice or frost where its cell's own temperature kept it. The ground carries relief at the grid's own scale: on the Earth-like's coastal plains some 117 m at a pixel's scale, 82 m of it without the unresolved craters' roughness (F54), the rest the other processes' relief at that scale. Where the ground hovers near a threshold, that relief decided each cell's side on its own, and the margin speckled: 228 lone cells of sea and 268 of land on the Earth-like's draft, 1,961 lone lakes over a tenth of the deepest on Titan's, 262 lone cells in the Earth-like's caps' ice, the snapshot's snow edge a scatter. The relief at a cell's scale deciding where a margin lies is noise placing, against PLAN principle 16 (physics places, noise realises). The user asked for the maps to be drawn by shares.

## Decision

1. A cell's ground is read as spread about its relief's mean (F416): the mean of the ground over a Gaussian of `RELIEF_SPREAD_CELLS`, two cells, about it, and the spread the root mean square of the ground less that mean over the same Gaussian. An even slope is its own mean and spreads by none. Two cells is the finest relief the grid draws as shape, twice the finest the relief at the cell's scale reaches (F54's octaves to two cells).
2. `liquid` holds the sea's mean depth over each cell, the ground under it spread as a normal distribution by that spread about that mean; `mask_liquid` the share of each pixel under the sea, 0 to 255 for 0 to 1; both read between cells.
3. `lakes` holds the lakes' mean depth over each cell by the same rule, each lake's level carried to the cells within two cells of it, each taking the highest lake beside it; read between cells, where it was nearest-cell.
4. `snapshot` draws the liquid at its mean depth, and the ices at the snapshot's phase by their share of each cell and their colour, each a mean over the same Gaussian of the cells' own; `shaded` lights the sea's flat level over the share of each cell it covers and the ground over the rest.
5. `ice`, each `ice_<species>` and `frost` hold their means over the same Gaussian.
6. The model's own fields are unchanged: the liquid the climate and the rivers read, the caps, the frost, and the cryosphere's window maps (ADR 0006), which stay per cell so that an engine can draw the ices at any phase from them.

## Consequences

- A pipeline that wants a binary coast thresholds `mask_liquid` at a half, which lies where the ground's shape crosses the level, within a cell on a plain the relief at the cell's scale speckles (the shares' test).
- An engine that sets its own sea level and draws its coast from `height` alone gets the cells' edge, all or none; the shares at another level need the ground's mean and spread, which the sidecar does not yet export.
- The maps soften every margin over some two cells, a sharp one (a cliff's coast, a cap's wall) too.
- The model's decisions on its own grid still read each cell all or none: the hypsometric fill's hollows under the level (F26), the frost's and the caps' cells. The maps no longer show that relief deciding them; the model's own decisions there wait for the processes that wear the relief down (the loop's erosion, creep).
