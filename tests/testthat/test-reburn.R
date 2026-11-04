## Regression test: in an extra reburn round (any round after `maxReburns[1]`) in which some fires
## stalled, Burn() recorded only the stalled fires. The fires that reached their target in that
## same round were never added to `fireSizes` or to the burned pixels, and their pixels stayed open
## to later rounds: they burned, but the year did not record them. Seen in production runs, e.g. a
## 4,552-pixel remainder lost in one year of a LacSeulUpland replicate.
##
## The landscape is one fire-return-interval zone in two parts: a solid block of flammable pixels
## (columns 1-25), and to its right single flammable pixels ("islands") at odd rows and columns,
## each ringed by non-flammable pixels. A fire ignited on an island with a target above one pixel
## stalls: fires cannot spread into non-flammable pixels, and `spread2()` only jumps on every 10th
## retry, beyond the default `maxRetriesPerID` of 4. A stalled fire is re-ignited from a random
## start cell in the zone, so of several re-ignited fires some land in the block and reach target
## while others land on an island and stall again. With `maxReburns = c(1L, 5L)` the second round
## is already an extra round. A 2-year interval gives enough fires; with this seed rounds 2 and 3
## each have fires of both kinds.
##
## `landmine_burn1()` is wrapped, not replaced: the wrapper returns the real burn unchanged and
## keeps a copy of each round, so the test knows what each round burned without relying on the
## bookkeeping it checks.

test_that("fires that reach target in an extra reburn round are recorded and burned", {
  landscape <- landmineIslandLandscape()

  rounds <- list()
  burn1 <- LandWebUtils::landmine_burn1
  local_mocked_bindings(
    landmine_burn1 = function(...) {
      fires <- burn1(...)
      ## copies: Burn() modifies `fires` by reference
      rounds[[length(rounds) + 1L]] <<- list(
        clusters = data.table::copy(attr(fires, "spreadState")$clusterDT),
        pixels = data.table::copy(fires)
      )
      fires
    },
    .package = "LandWebUtils"
  )

  ## one burn year, so round r is the r-th burn, and every fire it burns carries `attempt == r`
  sim <- runLandMineTest(landscape, paths = testPaths, end = 1,
                         params = list(maxReburns = c(1L, 5L)))
  fs <- sim$fireSizes[["1"]]
  burned <- terra::values(sim$rstCurrentBurn, mat = FALSE)

  ## the module's own test: a fire stalled when `size` differs from its rounded `maxSize`
  reached <- function(cl) cl$size == LandR::asInteger(cl$maxSize)
  mixed <- Filter(function(r) {
    ok <- reached(rounds[[r]]$clusters)
    any(ok) && any(!ok)
  }, seq_along(rounds)[-1L])
  expect_gt(length(mixed), 0) ## without such a round, the checks below test nothing

  for (r in mixed) {
    cl <- rounds[[r]]$clusters
    ## every fire of the round has a row with what it burned: those that stalled, and those that
    ## reached target
    expect_equal(sort(as.numeric(fs$size[fs$attempt == r])), sort(as.numeric(cl$size)),
                 info = paste("round", r))
    ## the pixels of the fires that reached target are burned
    px <- rounds[[r]]$pixels
    px <- px$pixels[px$initialPixels %in% cl$initialPixels[reached(cl)]]
    expect_gt(length(px), 0)
    expect_true(all(burned[px] == 1L), info = paste("round", r))
  }
})

## Regression test: a later burn round of a year could burn a pixel an earlier round had burned and
## recorded, so the year recorded that pixel twice. Burn() masks those pixels by setting their
## `spreadProb` to NA, but the start cells of re-issued fires were drawn from the whole zone, and
## `landmine_burn1()` spread its first step with probability 1, ignoring the mask. With this seed
## the year recorded 867 pixels for 847 burned (this file run on its own; the draws depend on what
## ran before in the session).
test_that("a pixel burned in one round of a year is not burned again in a later round", {
  rounds <- list()
  burn1 <- LandWebUtils::landmine_burn1
  local_mocked_bindings(
    landmine_burn1 = function(...) {
      args <- list(...)
      fires <- burn1(...)
      rounds[[length(rounds) + 1L]] <<- list(
        startCells = args$startCells,
        pixels = fires$pixels,
        ## what earlier rounds burned: flammable pixels the module masked out of `spreadProb`
        masked = which(is.na(args$spreadProb) & !is.na(args$spreadProbRel))
      )
      fires
    },
    .package = "LandWebUtils"
  )

  sim <- runLandMineTest(landmineIslandLandscape(), paths = testPaths, end = 1, seed = 2L,
                         params = list(maxReburns = c(1L, 5L)))
  fs <- sim$fireSizes[["1"]]
  burned <- terra::values(sim$rstCurrentBurn, mat = FALSE)

  expect_gt(length(rounds), 2L) ## without later rounds, the checks below test nothing
  expect_gt(length(rounds[[length(rounds)]]$masked), 0L)
  for (r in seq_along(rounds)[-1L]) {
    expect_false(any(rounds[[r]]$startCells %in% rounds[[r]]$masked), info = paste("round", r))
    expect_false(any(rounds[[r]]$pixels %in% rounds[[r]]$masked), info = paste("round", r))
  }
  ## the area the year recorded is the area it burned
  expect_identical(as.integer(sum(fs$size)), sum(burned == 1L, na.rm = TRUE))
})

## A zone burned completely within a year can take no more fires: with this seed the lake island
## burns out and still needs a fire. Its start cell was drawn on a burned pixel, re-burning it each
## round to the reburn ceiling: 19 pixels recorded for 9 burned (this file run on its own).
test_that("a zone with no unburned pixel left gets no more fires that year, and says so", {
  landscape <- landmineLakeIslandLandscape()
  ## 2 years, though only year 1 is checked: with one year and several zones, the module's
  ## area-burned plot (`geom_area()`) warns
  expect_message(
    sim <- runLandMineTest(landscape, paths = testPaths, end = 2, seed = 4L,
                           params = list(maxReburns = c(1L, 5L))),
    "reburn-no-start-cell year=1 FRI=1 nFires=1"
  )
  fs <- sim$fireSizes[["1"]]
  ## year 1's burn: the cumulative burn map less year 2's
  burned <- terra::values(sim$burnMap, mat = FALSE) - terra::values(sim$rstCurrentBurn, mat = FALSE)
  expect_true(all(burned[landscape$island] == 1L))
  expect_identical(as.integer(sum(fs$size)), sum(burned == 1L, na.rm = TRUE))
})
