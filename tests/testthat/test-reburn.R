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
  landscape <- landmineTestLandscape()
  rtm <- landscape$objects$rasterToMatch
  cells <- seq_len(terra::ncell(rtm))
  cols <- terra::colFromCell(rtm, cells)
  rows <- terra::rowFromCell(rtm, cells)
  flam <- terra::rast(rtm)
  terra::values(flam) <- as.integer(cols <= 25L | (cols > 26L & cols %% 2L == 1L & rows %% 2L == 1L))
  fri <- terra::rast(rtm)
  terra::values(fri) <- 2L
  landscape$objects$flammableMap <- flam
  landscape$objects$fireReturnInterval <- fri

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
