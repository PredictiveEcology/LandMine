Known issues: <https://github.com/PredictiveEcology/LandMine/issues>

# LandMine (development version)

Unreleased changes on the `development` branch since the 1.0.1 release:

* Standardized fire object names to the scfm/fireSense convention: `rstFlammable` renamed to `flammableMap`, and `rstCurrentBurnCumulative` renamed to `burnMap`.
* Re-enabled the `registerOutputs()` calls (previously disabled pending a fix).
* Fixed a `terra::values()` argument typo ("wat" to "mat").
* Fire-spread speedup: `spreadProb` and `spreadProbRel` are now passed to `spread2()` as numeric vectors rather than rasters. `spread2()` re-materialises a raster in full on every iteration, and since `landmine_burn1()` calls `spread2(iterations = 1L)` in a loop, that was one O(ncell) read per spread step. Measured 11-13x faster end-to-end, bit-identical (the same cells burn). Requires SpaDES.tools >= 2.1.2.9000 (PredictiveEcology/SpaDES.tools#106).
* `omitPixels` is now computed once per `Burn()` event rather than on every reburn iteration; `sim$flammableMap` is not modified within the event, so it was loop-invariant.
* The eligible start-cell pool is now built once per `Burn()` event rather than re-filtered and re-`na.omit()`ed on every reburn iteration. On a 4.5 Mpix study area this is ~306 ms to ~24 ms per round; at ~8 rounds/yr that is roughly 2.2 s per simulated year. Verified to produce the same start cells in the same order for the same seed.
* Dropped `magrittr` from `reqdPkgs`: the start-cell pipeline was its only remaining use.
* `LandMine.Rmd`'s calibration now specifies fire sizes in **hectares** (`landmine_optim_fireSizes()`) rather than pixel counts, and uses the repaired `objective = "andison"` with common random numbers. A fixed pixel-count vector meant a different physical fire at every resolution, which put the v3 (240 m) calibration outside the 10-3,000 ha range Andison fitted.
* The fire-spread calibration moved into LandWebUtils (>= 1.0.3.9024) as `landmine_optim_unpack()`, `landmine_optim_params_read()`/`_append()`, `landmine_optim_landscape()` and `landmine_optim_calibrate()`. `LandMine.R` now reads its parameters through those functions instead of re-implementing the CSV schema and the `10^` convention inline.
* `LandMine.Rmd`'s calibration section was rewritten against those functions and its stale chunks removed. Because the vignette is knitted with `eval = FALSE`, that code had never run and had drifted: it used the removed `raster` package API, `raster`'s argument names on a `terra` call, and an invalid `read.csv(file, row.names = FALSE)`. Parameter values shown there were also years out of date (`sizeCutoffs` of 8000/20000 against the 1629/52016 actually in use).

# LandMine 1.0.1 (2026-06-30)

* Fire-spread re-optimization: reran the LandMine DEoptim fire-spread optimization and adopted the latest optimized parameter values for 240 m pixels.
* Corrected the usage of `registerOutputs()` (temporarily disabled at first until it was fixed).
* Reproducible argument rename: `EstimateTruncPareto()` now uses `cachePath` instead of the deprecated `cacheRepo`.
* Cleanup and verification of the terra-based updates.

# LandMine 1.0.0 (2025-10-06)

* Major terra migration: reworked the module to use `terra` and stop using `quickPlot`/`raster`. `expectsInput`/`createsOutput` now declare `SpatRaster` instead of `Raster`/`RasterLayer`; `compareRaster()` to `terra::compareGeom()`, `raster()` to `terra::rast()`, and use of `terra::freq()`, `terra::values()`, `terra::mask()`, and `terra::deepcopy()`.
* Updated `reqdPkgs`: dropped `quickPlot`, `ggspatial`, `gridExtra`, and `fasterize`; added `terra`, `tidyterra`, `cli`, `RColorBrewer`, `VGAM`, and `stats`.
* Rebuilt all plotting on ggplot2 + tidyterra (`geom_spatraster()`, `geom_spatvector()`, `scale_fill_*` palettes) in place of `quickPlot::Plot()` and manual device management; figures are written to `figurePath()` and tracked via `registerOutputs()`.
* Minor Rmd chunk-option formatting.

# LandMine 0.0.6 (2025-08-13)

* "Burny" discontinuous-fuel spread: fixed and improved the `'burny'` `ROStype` fire spread for areas with discontinuous fuels, enabling the previously-disabled `burny` branch keyed on non-flammable (NA or 0) pixels.
* Robustness for study areas adjacent to polygons with no LTHFC, plus Douglas-fir-related fixes.
* Added checks and reporting for cases where no fires are started in a period.
* Reburn logic: allow different `maxReburns` values for each reburn phase, and updated the default `maxReburns` values.
* Additional DEoptim runs for fire-spread parameters.
* Plotting and raster fixes: corrected `rstCurrentBurnCumulative` being made all-zeros in two places; addressed `quickPlot` device proliferation (use the "transparent" option, limit to a single device for LandWeb); ensure `fasterize` is available so `rasterToMatch` has data; handle `sppEquiv` sizing and sf-vs-Spatial `studyArea` classes.
* Fixed `studyAreaName` inconsistencies and suppressed messages when loading simLists.

# LandMine 0.0.3 (2023-09-14)

* Two-phase reburn logic: documented and reworked reburning into two phases. In phase one, fires that did not reach their target size are reignited from new pixels within the FRI zone (less likely to stay stuck in regions with sinuous fires or discontinuous fuels); in phase two, new smaller fires are ignited with target sizes set to the burned-area shortfall.
* Parameter changes: `maxReburns` default changed from 10 to 5; `maxRetriesPerID` upper bound raised from 20 to 299 to allow larger fire sizes with discontinuous fuels.
* Added the `optimParsRowID` parameter to select which row of `LandMine_DEoptim_params.csv` to use (row 1 is the original 2018 values at 100 m pixels; other rows were computed at 250 m pixels).
* Added a `LandWebUtils` (>= 0.1.7) dependency; documentation cleanups.

# LandMine 0.0.2 (2023-09-07)

* Added the `ROStype` parameter ('burny', 'equal', 'log', or 'default') to increase spread in areas with discontinuous fuels (e.g., shield).
* Added an LTHFC/FRI summary plot and table.
* Manual/vignette updates and added the andiso2019 citation.
