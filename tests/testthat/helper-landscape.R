## A small synthetic landscape for tests that run the module end to end, and the run itself.
##
## Four fire-return-interval zones in vertical bands of 20 columns (40 rows x 80 columns of
## 250 m pixels, 20,000 ha), with the same interval values as the module's default map:
##
##   columns  1-20  FRI  60   fully flammable, inside `studyArea`
##   columns 21-40  FRI 100   a 10 x 20 non-flammable block (a lake): 200 of its 800 pixels
##   columns 41-60  FRI 120   fully flammable; its top two rows carry FRI 0 ("no interval")
##   columns 61-80  FRI 250   columns 71-80 lie outside `studyArea`: 400 of its 800 pixels
##
## `studyArea` is the rectangle over columns 1-70, with edges on cell boundaries so that
## rasterizing it selects exactly those columns.
##
## With `speciesInputs = TRUE` the species inputs are supplied too, with `sppEquivCol =
## "LandWeb"`, so that LandMine 1.0.13 -- whose species defaults did not run -- can run the same
## landscape for comparison. Pixel groups, cohorts, time since fire and the ROS table always
## come from `.inputObjects`.
landmineTestLandscape <- function(res = 250, nrow = 40L, ncol = 80L, speciesInputs = TRUE) {
  crs <- "EPSG:3400"
  rtm <- terra::rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = ncol * res,
                     ymin = 0, ymax = nrow * res, crs = crs)
  terra::values(rtm) <- 1L

  cells <- seq_len(terra::ncell(rtm))
  cols <- terra::colFromCell(rtm, cells)
  rows <- terra::rowFromCell(rtm, cells)

  lake <- cols %in% 26:35 & rows %in% 11:30
  zero <- cols %in% 41:60 & rows %in% 1:2
  outside <- cols > 70L

  fri <- terra::rast(rtm)
  terra::values(fri) <- ifelse(zero, 0L, c(60L, 100L, 120L, 250L)[(cols - 1L) %/% 20L + 1L])

  flam <- terra::rast(rtm)
  terra::values(flam) <- ifelse(lake, 0L, 1L)

  sa <- terra::as.polygons(terra::ext(0, 70 * res, 0, nrow * res), crs = crs)

  objects <- list(
    studyArea = sa,
    studyAreaReporting = sa,
    rasterToMatch = rtm,
    flammableMap = flam,
    fireReturnInterval = fri
  )
  params <- list()
  if (isTRUE(speciesInputs)) {
    spp <- c("Pice_gla", "Pinu_spp")
    sppEquiv <- LandWebUtils::landweb_sppEquiv(LandR::sppEquivalencies_CA)
    sppEquiv <- sppEquiv[sppEquiv[["LandWeb"]] %in% spp, ]
    objects$species <- data.table::data.table(species = spp, speciesCode = seq_along(spp))
    objects$sppEquiv <- sppEquiv
    objects$sppColorVect <- LandR::sppColors(sppEquiv, "LandWeb", newVals = "Mixed",
                                             palette = "Accent")
    params$sppEquivCol <- "LandWeb"
  }

  list(
    objects = objects,
    params = params,
    lake = which(lake),
    zero = which(zero),
    outside = which(outside)
  )
}

## Run LandMine in "single" mode on `landscape` from `start` to `end` with a fixed seed.
## `params` overrides or adds LandMine parameters, after the landscape's own. Plots are off.
runLandMineTest <- function(landscape, paths, start = 0, end = 5, seed = 20260928L,
                            params = list()) {
  withr::local_seed(seed)
  withr::local_options(LandR.verbose = 0)
  SpaDES.core::simInitAndSpades(
    times = list(start = start, end = end),
    params = list(LandMine = utils::modifyList(utils::modifyList(list(
      .plots = NA_character_,
      .useCache = FALSE,
      .useParallel = 1,
      biggestPossibleFireSizeHa = 1e4,
      maxReburns = c(1L, 2L),
      mode = "single"
    ), landscape$params), params)),
    modules = "LandMine",
    objects = landscape$objects,
    paths = paths
  )
}
