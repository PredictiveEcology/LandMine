Known issues: <https://github.com/PredictiveEcology/LandMine/issues>

# LandMine (development version)

Unreleased changes on the `development` branch since the 1.0.1 release:

* Standardized fire object names to the scfm/fireSense convention: `rstFlammable` renamed to `flammableMap`, and `rstCurrentBurnCumulative` renamed to `burnMap`.
* Re-enabled the `registerOutputs()` calls (previously disabled pending a fix).
* Fixed a `terra::values()` argument typo ("wat" to "mat").

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
