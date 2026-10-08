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

test_that("species simulated separately run with fuel from their reporting group", {
  ## With one code per species, the species-name fuel typing has no rate of spread for most codes
  ## (e.g. Pinu_con) and stops at the first fire; sppEquivColFuel types fuel from the leading
  ## reporting group instead.
  landscape <- landmineTestLandscape()
  spp <- c("Pice_gla", "Pinu_con")
  sppEquiv <- LandWebUtils::landweb_species_sppEquiv(LandR::sppEquivalencies_CA)
  sppEquiv <- sppEquiv[sppEquiv[["LandWeb"]] %in% spp, ]
  landscape$objects$species <- data.table::data.table(species = spp, speciesCode = seq_along(spp))
  landscape$objects$sppColorVect <- LandR::sppColors(sppEquiv, "LandWeb", newVals = "Mixed", palette = "Accent")

  ## a copy per run: Init() adds its `LandMine` column to the given table by reference
  landscape$objects$sppEquiv <- data.table::copy(sppEquiv)
  expect_error(runLandMineTest(landscape, paths = testPaths, end = 2), "Pinu_con")

  landscape$objects$sppEquiv <- data.table::copy(sppEquiv)
  sim <- runLandMineTest(landscape, paths = testPaths, end = 3,
                         params = list(sppEquivColFuel = "LandWebReport"))
  expect_s4_class(sim, "simList")
  cmp <- as.data.frame(SpaDES.core::completed(sim))
  expect_equal(as.numeric(cmp$eventTime[cmp$eventType == "Burn"]), c(1, 2, 3))
  expect_false("LandMine" %in% names(sim$sppEquiv)) ## sppEquiv is left as given
})

