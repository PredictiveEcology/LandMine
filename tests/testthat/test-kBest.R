## The fire-size shape `kBest` is fitted at init by optimize() over an objective that draws 1e6
## random fire sizes per evaluation, so unseeded replicates each get their own fire-size
## distribution. `useSeed` makes the fit deterministic without moving the simulation's random
## stream, and the `kBest` parameter skips the fit. Each fit takes a second or two.

test_that("the same useSeed gives the same kBest, whatever the simulation's stream", {
  landscape <- landmineTestLandscape()
  a <- initLandMineTest(landscape, testPaths, seed = 1L, params = list(useSeed = 11L))
  b <- initLandMineTest(landscape, testPaths, seed = 2L, params = list(useSeed = 11L))
  ## a session on another generator, as parallel workers often are
  c <- withr::with_seed(
    3L, initLandMineTest(landscape, testPaths, seed = 3L, params = list(useSeed = 11L)),
    .rng_kind = "L'Ecuyer-CMRG"
  )
  d <- initLandMineTest(landscape, testPaths, seed = 1L, params = list(useSeed = 12L))

  expect_identical(b$sim$kBest, a$sim$kBest)
  expect_identical(c$sim$kBest, a$sim$kBest)
  expect_false(identical(d$sim$kBest, a$sim$kBest))

  ## the fit a default R session gives after `set.seed(11L)`, so it can be reproduced outside a run
  expect_identical(
    a$sim$kBest,
    withr::with_seed(
      11L, LandWebUtils::landmine_estimate_kBest(biggestPossibleFireSizeHa = 1e4),
      .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection"
    )
  )
})

test_that("a seeded fit leaves the simulation's random stream as it was", {
  landscape <- landmineTestLandscape()
  seeded <- initLandMineTest(landscape, testPaths, params = list(useSeed = 11L))
  otherSeed <- initLandMineTest(landscape, testPaths, params = list(useSeed = 12L))
  supplied <- initLandMineTest(landscape, testPaths, params = list(kBest = seeded$sim$kBest))
  unseeded <- initLandMineTest(landscape, testPaths)

  ## with a seed, the stream after init is where it is when no fit runs, whatever the seed.
  ## `expect_true(identical())` because testthat's 2nd-edition `expect_identical()` errors while
  ## describing a difference between two `.Random.seed`s.
  expect_true(identical(seeded$rng, supplied$rng))
  expect_true(identical(otherSeed$rng, supplied$rng))
  ## the default still fits from the simulation's stream, as before `useSeed` was read
  expect_false(identical(unseeded$rng, supplied$rng))

  ## the fit runs on R's default generator; the session's own is back afterwards
  lecuyer <- withr::with_seed(
    3L, initLandMineTest(landscape, testPaths, seed = 3L, params = list(useSeed = 11L)),
    .rng_kind = "L'Ecuyer-CMRG"
  )
  lecuyerSupplied <- withr::with_seed(
    3L, initLandMineTest(landscape, testPaths, seed = 3L, params = list(kBest = seeded$sim$kBest)),
    .rng_kind = "L'Ecuyer-CMRG"
  )
  expect_true(identical(lecuyer$rng, lecuyerSupplied$rng))
})

test_that("useSeed burns the same fires as the kBest it fits, supplied", {
  landscape <- landmineTestLandscape()
  seeded <- runLandMineTest(landscape, paths = testPaths, end = 3, params = list(useSeed = 11L))
  supplied <- runLandMineTest(landscape, paths = testPaths, end = 3,
                              params = list(kBest = seeded$kBest))

  expect_identical(supplied$kBest, seeded$kBest)
  expect_gt(sum(terra::values(seeded$burnMap, mat = FALSE), na.rm = TRUE), 0)
  expect_identical(terra::values(supplied$burnMap, mat = FALSE),
                   terra::values(seeded$burnMap, mat = FALSE))
  ## a year without fires is recorded as `size` and `maxSize` only, hence `fill`
  fires <- function(sim) {
    as.data.frame(data.table::rbindlist(sim$fireSizes, idcol = "year", fill = TRUE))
  }
  expect_identical(fires(supplied), fires(seeded))
})

test_that("a supplied kBest is used as given, and the fit does not run", {
  landscape <- landmineTestLandscape()
  testthat::local_mocked_bindings(
    landmine_estimate_kBest = function(...) stop("the kBest fit ran"),
    .package = "LandWebUtils"
  )
  ## the mock is what init calls when it fits
  expect_error(initLandMineTest(landscape, testPaths), "the kBest fit ran")
  expect_error(initLandMineTest(landscape, testPaths, params = list(useSeed = 11L)),
               "the kBest fit ran")

  res <- initLandMineTest(landscape, testPaths, params = list(kBest = 0.6, useSeed = 11L))
  expect_identical(res$sim$kBest, 0.6)
})

test_that("an unusable kBest is refused", {
  landscape <- landmineTestLandscape()
  msg <- "`kBest` must be a single positive number"
  expect_error(initLandMineTest(landscape, testPaths, params = list(kBest = -0.5)), msg)
  expect_error(initLandMineTest(landscape, testPaths, params = list(kBest = c(0.6, 0.7))), msg)
  expect_error(initLandMineTest(landscape, testPaths, params = list(kBest = "0.6")), msg)
})
