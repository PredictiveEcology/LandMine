Known issues: <https://github.com/PredictiveEcology/LandMine/issues>

# LandMine (development version)

Unreleased changes on the `development` branch since the 1.0.1 release:

* Standardized fire object names to the scfm/fireSense convention: `rstFlammable` renamed to `flammableMap`, and `rstCurrentBurnCumulative` renamed to `burnMap`.
* Re-enabled the `registerOutputs()` calls (previously disabled pending a fix).
* Fixed a `terra::values()` argument typo ("wat" to "mat").
* Fire-spread speedup: `spreadProb` and `spreadProbRel` are now passed to `spread2()` as numeric vectors rather than rasters. `spread2()` re-materialises a raster in full on every iteration, and since `landmine_burn1()` calls `spread2(iterations = 1L)` in a loop, that was one O(ncell) read per spread step. Measured 11-13x faster end-to-end, bit-identical (the same cells burn). Requires SpaDES.tools >= 2.1.2.9000 (PredictiveEcology/SpaDES.tools#106).
* `omitPixels` is now computed once per `Burn()` event rather than on every reburn iteration; `sim$flammableMap` is not modified within the event, so it was loop-invariant.
* The eligible start-cell pool is now built once per `Burn()` event rather than re-filtered and re-`na.omit()`ed on every reburn iteration. On a 4.5 Mpix study area this is ~306 ms to ~24 ms per round; at ~8 rounds/yr that is roughly 2.2 s per simulated year. Verified to produce the same start cells in the same order for the same seed.
* **The same extent mismatch was in the model, not only in the metric.** `Init()` now passes `studyArea` to `landmine_ignition_budget()`. The per-zone fire budget is `(zone area / mean fire size) / FRI` computed from a `fireReturnInterval` raster masked only by flammability, so area outside the study area -- which `ROSmap` and `mod$spreadProb` are both masked to, and where a fire can therefore burn only its own start cell -- was inflating every zone's expected fires per year. The start-cell pool is built from that same raster, so those pixels were also drawn as ignition points. On WesternAlbertaUpland the budget was 2,918 fires/yr against a correct 2,015: **31% of every year's ignitions aimed at unburnable ground**, with one zone (99.95% outside the polygon) over-allocated 2,008-fold, 45.6 fires/yr against 0.023. NOTE: this changes simulated fire regimes, so existing `mainSim` outputs are not comparable to future ones.
* **Fixed a metric that made the fire model look broken when it was not.** `SummarizeFRI`/`SummarizeFRImulti` now pass `studyArea` to `landmine_fri_summary()`. `ROSmap` and `mod$spreadProb` are both masked to `sim$studyArea`, so a pixel outside that polygon can be neither ignited nor spread into -- but `fireReturnInterval` is not masked and routinely overhangs it, and every such pixel was entering the achieved-interval denominator while contributing no opportunity to burn. On WesternAlbertaUpland that was 988,437 of 3,361,240 flammable zone pixels (29.4%), burning at a mean rate of 0.51 against 16.59 inside, with one zone down to 73 burnable pixels out of 146,627. Two zones appeared to under-burn by 3.5x and 12.2x; with the correct denominator **every zone falls between 0.90 and 1.07**. `LandMine_FRI_summary.csv`/`_multi.csv` from earlier runs are wrong for any study area whose LTHFC raster overhangs its polygon.
* New per-zone attainment diagnostics, written every run so nobody has to read a log to learn a zone missed its target: `LandMine_FRI_diagnostics.csv` (`_multi.csv`), a two-panel map `LandMine_FRI_zones_diagnostic.png` (target interval beside achieved/target ratio), `LandMine_FRI_drivers.png` (what separates the off-target zones), and a one-line verdict messaged at the end of the summary. The table reports the ratio over the burnable area alongside the uncorrected one and each zone's share inside the study area, so a discrepancy against older outputs is explainable rather than mysterious. New `friDiagnostics` output object; requires LandWebUtils >= 1.0.3.9030.
* `ROStype`'s metadata description claimed four valid values ('burny', 'equal', 'log', 'default'), but only 'burny' was ever branched on and the other two behaved exactly like 'default'. LandWeb v3 uses only 'default' and 'burny', so the description now says so and `landmine_fire_ros()` rejects the rest.
* `fireROS()` moved into LandWebUtils (>= 1.0.3.9029) as `landmine_fire_ros()`, together with the leading-species-to-fuel-type table as `landmine_known_species()`. It was the module's largest untested internal (~147 lines reading eight values off `sim`/`mod`/`P(sim)`), so it could only be exercised by running a simulation; it is now a pure function of rasters and tables, covered by unit tests and by equivalence testing against the previous inline implementation over 200 randomised landscapes.
* Two long-standing behaviours in that code are now **documented and pinned by tests, but deliberately unchanged**, because either would alter the fire regime and both lie in the code path currently suspected in the FRI zone 55/170 under-burning: mixedwood stands never receive the `ROSTable`'s `mixed` rates (they fall through to `ROSother`, i.e. the mature-spruce rate at every age, leaving both `mixed` rows of the table dead), and the young age class's species filter is applied under the opposite condition to the mature and immature ones -- an unfinished copy-paste, per the git history. `landmine_fire_ros()`'s `youngGuard` argument exists so the latter can be measured rather than argued about.
* Dropped `magrittr` from `reqdPkgs`: the start-cell pipeline was its only remaining use.
* Re-calibrated against the repaired objective, at 120 m and 240 m. `LandMine_DEoptim_params.csv` now holds only those two fits: every earlier row used a different objective (a size penalty that never fired, a constant perimeter:area target, and fire sizes in pixels rather than hectares) and is not comparable to them. The old rows remain in git history.
* `optimParsRowID` now defaults to `1L`, the **120 m** fit, at every resolution. The 120 m fit transfers to 240 m (0.70 vs 0.81, 0.33 pooled sd) with a fifth of the variance, while the 240 m fit is much worse at 120 m (1.78 vs 0.61, 4.91 sd). It is also far better determined: repeat runs land within 4.6% of the search interval against 24.8% at 240 m. Both 240 m fits score identically despite differing by 24.8%, i.e. that resolution's optimum is a broad plateau rather than a point.
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
