## Runs of the whole module: on its own defaults, and the plot event schedule.

test_that("the module runs on its own species defaults", {
  ## Only the landscape is supplied; species, cohorts, `sppEquiv`, `sppColorVect` and
  ## `sppEquivCol` come from the module. The default `sppEquiv` used to be
  ## `LandR::sppEquivalencies_CA`, which has no `LandWeb` column, so Init() failed reading it.
  landscape <- landmineTestLandscape(speciesInputs = FALSE)
  sim <- runLandMineTest(landscape, paths = testPaths, end = 3)

  expect_s4_class(sim, "simList")
  expect_false(anyNA(sim$sppEquiv$LandMine))
  cmp <- as.data.frame(SpaDES.core::completed(sim))
  expect_equal(as.numeric(cmp$eventTime[cmp$eventType == "Burn"]), c(1, 2, 3))
})

test_that("the plot event recurs every .plotInterval after .plotInitialTime", {
  ## It was rescheduled at `P(sim)$.plotInterval` instead of `time(sim) + .plotInterval`.
  ## SpaDES.core reads a unitless event time as `eventTime - start(sim) + time(sim)`, so that
  ## worked by accident when `start(sim)` is 0. A nonzero start exposes it.
  landscape <- landmineTestLandscape()
  sim <- runLandMineTest(landscape, paths = testPaths, start = 10, end = 16,
                         params = list(.plotInterval = 2))

  ## `.plotInitialTime` defaults to start + 1; the last one is the end-of-run plot
  cmp <- as.data.frame(SpaDES.core::completed(sim))
  expect_equal(as.numeric(cmp$eventTime[cmp$eventType == "plot"]), c(11, 13, 15, 16))
})
