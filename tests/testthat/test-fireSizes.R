## Regression test: a burn year with no fires stored `data.table(size = 0L, maxSize = 0L)` in
## `sim$fireSizes`, two columns against the five of a year with fires. burnSummaries binds the
## years with `rbindlist(sim$fireSizes, idcol = "year")`, which stopped on the mismatch ("Item 446
## has 2 columns, inconsistent with item 1 which has 5 columns"), and replicates in four
## lower-fire study areas failed. A fire-free year is now a zero-row table with the same columns.
##
## Fire-return intervals five times the landscape's (300-1,250 years) make most years fire-free.
## A year is fire-free when nothing burned in it, which holds for the old placeholder too, so the
## schema checks below run under either. The run must contain both kinds of year, so the test
## cannot pass on a run without one.

test_that("a fire-free year has the columns of a fire year and adds no rows", {
  landscape <- landmineTestLandscape()
  landscape$objects$fireReturnInterval <- landscape$objects$fireReturnInterval * 5L
  sim <- runLandMineTest(landscape, paths = testPaths, end = 10)

  burned <- vapply(sim$fireSizes, function(x) sum(x$size), numeric(1))
  expect_true(any(burned == 0))
  expect_true(any(burned > 0))

  fireYear <- vapply(sim$fireSizes[[which(burned > 0)[1]]], class, character(1))
  for (year in names(burned)[burned == 0]) {
    expect_identical(vapply(sim$fireSizes[[year]], class, character(1)), fireYear)
    expect_identical(nrow(sim$fireSizes[[year]]), 0L)
  }

  fs <- data.table::rbindlist(sim$fireSizes, idcol = "year")
  expect_identical(nrow(fs), sum(vapply(sim$fireSizes, nrow, integer(1))))
  expect_setequal(unique(fs$year), names(burned)[burned > 0])
})
