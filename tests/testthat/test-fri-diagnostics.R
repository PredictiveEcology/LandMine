## Regression test: the per-zone FRI diagnostics reported `pctFlam` = 100 and
## `pctInStudyArea` = 100 for every zone. `landmine_fri_metrics()` measures both on the
## interval raster the summaries read, and that raster had already been masked to flammable
## pixels (in `.inputObjects`) and to flammable-and-inside-`studyArea` pixels (Init() wrote the
## ignition budget's copy back to `sim$fireReturnInterval`). The fire code now keeps the
## masked copy to itself, so the summaries see each zone's whole extent.
##
## The landscape (see helper-landscape.R) has a zone that is 75% flammable and a zone that is
## 50% inside `studyArea`, so both columns have a known value other than 100. One short run
## serves both checks; it takes about half a minute, most of it loading packages.

test_that("FRI diagnostics see whole zones while fires use the masked copy", {
  landscape <- landmineTestLandscape()
  sim <- runLandMineTest(landscape, paths = testPaths)

  ## the diagnostics count non-flammable pixels and pixels outside `studyArea`
  diag <- sim$friDiagnostics
  expect_equal(sort(diag$LTHFC), c(60, 100, 120, 250))

  pctFlam <- stats::setNames(diag$pctFlam, diag$LTHFC)
  pctInSA <- stats::setNames(diag$pctInStudyArea, diag$LTHFC)
  expect_equal(pctFlam[["100"]], 75)       ## 200 of 800 pixels are the lake
  expect_equal(pctInSA[["250"]], 50)       ## 400 of 800 pixels are outside
  expect_equal(unname(pctFlam[c("60", "120", "250")]), c(100, 100, 100))
  expect_equal(unname(pctInSA[c("60", "100", "120")]), c(100, 100, 100))

  ## `sim$fireReturnInterval` is the input, except that FRI 0 ("no interval") becomes NA
  fri <- terra::values(sim$fireReturnInterval, mat = FALSE)
  expect_true(all(fri[landscape$lake] == 100L))
  expect_true(all(fri[landscape$outside] == 250L))
  expect_true(all(is.na(fri[landscape$zero])))

  ## no ignition or spread on non-flammable pixels or outside `studyArea`: the start-cell
  ## pool and the fire landscape must still be the masked copy
  burned <- terra::values(sim$burnMap, mat = FALSE)
  expect_gt(sum(burned), 0)
  expect_true(all(burned[c(landscape$lake, landscape$outside)] == 0))
})
