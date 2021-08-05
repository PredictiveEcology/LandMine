defineModule(sim, list(
  name = "LandMine",
  description = "Reimplementation of Andison (1999) LandMine fire model",
  keywords = c("Fire", "Landscape", "Percolation", "Pixel-based"),
  authors = c(
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@canada.ca", role = c("aut", "cre")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.3.9009", LandMine = numeric_version("0.0.1")),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "LandMine.Rmd"),
  reqdPkgs = list("assertthat", "data.table", "grDevices", "magrittr", "ggplot2",
                  "raster", "RColorBrewer", "VGAM",
                  "PredictiveEcology/LandR@development",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/SpaDES.tools@development"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
    defineParameter("biggestPossibleFireSizeHa", "numeric", 1e6, 1e4, 2e6,
                    "An upper limit, in hectares, of the truncated Pareto distribution of fire sizes"),
    defineParameter("burnInitialTime", "numeric", start(sim, "year") + 1, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter("fireTimestep", "numeric", 1, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter("flushCachedRandomFRI", "logical", FALSE, NA, NA,
                    "If no Fire Return Interval map is supplied, then a random one will be created and cached. Use this to make a new one."),
    defineParameter("minPropBurn", "numeric", 0.90, 0.00, 1.00,
                    "Minimum proportion burned pixels to use when triggering warnings about simulated fires."),
    defineParameter("mixedType", "numeric", 2,
                    desc = paste("How to define mixed stands: 1 for any species admixture;",
                                 "2 for deciduous > conifer. See ?vegTypeMapGenerator.")),
    defineParameter("maxRetriesPerID", "integer", 10L, 0L, 20L,
                    "Number of attempts that will be made per event ID, before abandoning. See `?SpaDES.tools::spread2`."),
    defineParameter("ROSother", "integer", 30L, NA, NA,
                    paste0("default ROS value for non-forest vegetation classes.",
                           "this is needed when passing a modified ROSTable, e.g. using log-transformed values.")),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    "The column in sim$specieEquivalency data.table to use as a naming convention"),
    defineParameter("useSeed", "integer", NULL, NA, NA,
                    paste("Only used for creating a starting cohortData dataset.",
                          "If NULL, then it will be randomly generated;",
                          "If non-NULL, will pass this value to set.seed and be deterministic and identical each time.",
                          "WARNING: setting the seed to a specific value will cause all simulations to be identical!")),
    defineParameter("vegLeadingProportion", "numeric", 0.8, 0, 1,
                    "a number that define whether a species is leading for a given pixel"),
    defineParameter(".plotInitialTime", "numeric", start(sim, "year") + 1, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should this entire module be run with caching activated?",
                          "This is generally intended for data-type modules,",
                          "where stochasticity and time are not relevant")),
    defineParameter(".unitTest", "logical", getOption("LandR.assertions", TRUE), NA, NA,
                    "Some functions can have internal testing. This will turn those on or off, if any exist"),
    defineParameter(".useParallel", "numeric", 2, NA, NA,
                    paste("Used in burning. Will be passed to data.table::setDTthreads().",
                          "NOTE: should be <= 2 as the additonal RAM overhead too high given marginal speedup."))
  ),
  inputObjects = bindrows(
    expectsInput("cohortData", "data.table",
                 desc = paste("Columns: B, pixelGroup, speciesCode (as a factor of the names), age.",
                              "indicating several features about the current vegetation of stand."),
                 sourceURL = NA),
    expectsInput("fireReturnInterval", "Raster",
                 desc = paste("A raster layer that is a factor raster, with at least 1 column called",
                              "'fireReturnInterval', representing the fire return interval in years."),
                 sourceURL = NA),
    expectsInput("pixelGroupMap", "RasterLayer",
                 desc = "Pixels with identical values share identical stand features",
                 sourceURL = NA),
    expectsInput("rasterToMatch", "RasterLayer",
                 #desc = "this raster contains two pieces of information: Full study area with fire return interval attribute",
                 desc = "DESCRIPTION NEEDED", # TODO: is this correct?
                 sourceURL = NA),
    expectsInput("rasterToMatchReporting", "RasterLayer",
                 desc = paste("Raster layer of study area used for plotting and reporting only.",
                              "Defaults to the kNN biomass map masked with `studyArea`"),
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"),
    expectsInput("ROSTable", "data.table",
                 desc = paste("A data.table with 3 columns, 'age', 'leading', and 'ros'.",
                              "The values under the 'age' column can be 'mature', 'immature',",
                              "'young' and compound versions of these, e.g., 'immature_young'",
                              "which can be used when 2 or more age classes share same 'ros'.",
                              "'leading' should be vegetation type.",
                              "'ros' gives the rate of spread values for each age and type."),
                 sourceURL = NA),
    expectsInput("rstFlammable", "Raster",
                 desc = paste("A raster layer, with 0, 1 and NA, where 1 indicates areas",
                              "that are flammable, 0 not flammable (e.g., lakes)",
                              "and NA not applicable (e.g., masked)")),
    expectsInput("rstTimeSinceFire", "Raster",
                 desc = "a time since fire raster layer",
                 sourceURL = NA),
    expectsInput("species", "data.table",
                 desc = "Columns: species, speciesCode, Indicating several features about species",
                 sourceURL = NA),
    expectsInput("sppColorVect", "character",
                 desc = "named character vector of hex colour codes corresponding to each species",
                 sourceURL = NA),
    expectsInput("sppEquiv", "data.table",
                 desc = paste("Multi-columned data.table indicating species name equivalencies.",
                              "Default taken from LandR sppEquivalencies_CA which has names for",
                              "species of trees in Canada"),
                 sourceURL = NA),
    expectsInput("studyAreaReporting", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (typically smaller/unbuffered than studyArea) to use for plotting/reporting.",
                              "Defaults to an area in Southwestern Alberta, Canada."),
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput("fireInitialTime", "numeric", paste(
      "The initial event time of the burn event.",
      "This is simply a reassignment from P(sim)$burnInitialTime.")
    ),
    createsOutput("fireSizes", "list", paste(
      "A list of data.tables, one per burn event, each with two columns, size and maxSize.",
      " These indicate the actual sizes and expected sizes burned, respectively.",
      "These can be put into a single data.table with rbindlist(sim$fireSizes, idcol = 'year')")
    ),
    createsOutput("fireReturnInterval", "RasterLayer", paste(
      "A Raster map showing the fire return interval. This is created from the rstCurrentBurn.")
    ),
    createsOutput("fireReturnIntervalsByPolygonNumeric", "numeric", paste(
      "A vector of the fire return intervals, ordered by the numeric representation of polygon ID")
    ),
    createsOutput("fireTimestep", "numeric", paste(
      "The number of time units between successive fire events in a fire module.")
    ),
    createsOutput("kBest", "numeric", paste(
      "A numeric scalar that is the optimal value of K in the Truncated Pareto distribution (rtruncpareto)")),
    createsOutput("numFiresPerYear", "numeric", paste(
      "The average number of fires per year, by fire return interval level on rstCurrentBurn.")
    ),
    createsOutput("rstCurrentBurn", "RasterLayer", paste(
      "A raster layer, produced at each timestep, where each",
      "pixel is either 1 or 0 indicating burned or not burned.")
    ),
    createsOutput("rstCurrentBurnCumulative", "RasterLayer", "Cumulative number of times a pixel has burned"),
    createsOutput("sppEquiv", "data.table",
                  desc = paste("Same as input, but with new column, LandMine"))
  )
))

doEvent.LandMine <- function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    ### check for more detailed object dependencies:
    ### (use `checkObject` or similar)

    # do stuff for this event
    #  ff package, used in SpaDES.tools::spread2, doesn't always set this correctly.
    options(fftempdir = tempdir())
    sim <- EstimateTruncPareto(sim)
    sim <- Init(sim)

    # schedule future event(s)
    sim <- scheduleEvent(sim, P(sim)$burnInitialTime, "LandMine", "Burn", 2.5)
    sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "LandMine", "plot")
    sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "LandMine", "save")
  } else if (eventType == "plot") {
    if (!is.na(P(sim)$.plotInitialTime) && (P(sim)$.plotInitialTime == time(sim)))
      mod$LandMineDevice <- max(dev.list()) + 1

    devCur <- dev.cur()
    quickPlot::dev(mod$LandMineDevice, width = 18, height = 12)
    sim <- plotFn(sim)
    dev(devCur)
    sim <- scheduleEvent(sim, P(sim)$.plotInterval, "LandMine", "plot")
  } else if (eventType == "Burn") {
    sim <- Burn(sim)
    sim <- scheduleEvent(sim, time(sim) + P(sim)$fireTimestep, "LandMine", "Burn", 2.5)
  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}

### initialization
EstimateTruncPareto <- function(sim, verbose = getOption("LandR.verbose", TRUE)) {
  if (verbose > 0)
    message("Estimate Truncated Pareto parameters")

  findK_upper <- function(params = c(0.4), upper1) {
    fs <- round(VGAM::rtruncpareto(1e6, 1, upper = upper1, shape = params[1]))
    #meanFS <- meanTruncPareto(k = params[1], lower = 1, upper = upper1, alpha = 1)
    #diff1 <- abs(quantile(fs, 0.95) - meanFS)
    #abs(sum(fs[fs>quantile(fs, 0.95)])/sum(fs) - 0.9) # "90% of area is in 5% of fires" # from Dave rule of thumb

    # Eliot Adjustment because each year was too constant -- should create greater variation
    abs(sum(fs[fs > quantile(fs, 0.95)]) / sum(fs) - 0.95) # "95% of area (2nd term) is in 5% of fires (1st term)"
    # Eliot Adjustment Oct 23, 2018 because each year still too constant -- should create greater variation
    #abs(sum(fs[fs > quantile(fs, 0.90)]) / sum(fs) - 0.95) # "95% of area (2nd term) is in 10% of fires (1st term)"
  }

  sim$kBest <- Cache(optimize, interval = c(0.05, 0.99), f = findK_upper,
                     upper1 = P(sim)$biggestPossibleFireSizeHa,
                     cacheRepo = cachePath(sim))$minimum
  return(invisible(sim))
}

Init <- function(sim, verbose = getOption("LandR.verbose", TRUE)) {
  ## DEBUGGING: random seed issues
  #fseed <- file.path(outputPath(sim), "seed.txt")
  #writeEventInfo(sim, fseed, append = TRUE)
  #writeRNGInfo(fseed, append = TRUE)
  ## END DEBUGGING

  compareRaster(sim$rasterToMatch, sim$fireReturnInterval, sim$rstFlammable, sim$rstTimeSinceFire)

  sim$fireSizes <- list()

  if (!is.integer(sim$fireReturnInterval[]))
    sim$fireReturnInterval[] <- as.integer(sim$fireReturnInterval[])

  if (verbose > 0)
    message("Initializing fire maps")
  sim$fireTimestep <- P(sim)$fireTimestep
  sim$fireInitialTime <- P(sim)$burnInitialTime

  # check sim$fireReturnInterval should have no zeros
  zeros <- sim$fireReturnInterval[] == 0L
  if (any(zeros, na.rm = TRUE)) sim$fireReturnInterval[zeros] <- NA_integer_
  numPixelsPerPolygonNumeric <- Cache(freq, sim$fireReturnInterval, useNA = "no", cacheRepo = cachePath(sim)) %>%
    na.omit()
  colnames(numPixelsPerPolygonNumeric) <- c("fri", "count")
  numPixelsPerPolygonNumeric <- cbind(value = seq_len(NROW(numPixelsPerPolygonNumeric)), numPixelsPerPolygonNumeric)
  #numPixelsPerPolygonNumeric <- cbind(numPixelsPerPolygonNumeric, fri = raster::factorValues(sim$rasterToMatch,
  #                                                              numPixelsPerPolygonNumeric[, "value"],
  #                                                              att = "fireReturnInterval")[, 1])
  ordPolygons <- order(numPixelsPerPolygonNumeric[, "value"])
  numPixelsPerPolygonNumeric <- numPixelsPerPolygonNumeric[ordPolygons, , drop = FALSE]
  sim$fireReturnIntervalsByPolygonNumeric <- numPixelsPerPolygonNumeric[, "fri"]
  numPixelsPerPolygonNumeric <- numPixelsPerPolygonNumeric[, "count"]
  names(numPixelsPerPolygonNumeric) <- sim$fireReturnIntervalsByPolygonNumeric

  numHaPerPolygonNumeric <- numPixelsPerPolygonNumeric * (prod(res(sim$fireReturnInterval)) / 1e4)
  returnInterval <- sim$fireReturnIntervalsByPolygonNumeric

  if (verbose > 0)
    message("Determine mean fire size")
  meanFireSizeHa <- meanTruncPareto(k = sim$kBest, lower = 1,
                                    upper = P(sim)$biggestPossibleFireSizeHa,
                                    alpha = 1)
  numFiresByPolygonNumeric <- numHaPerPolygonNumeric / meanFireSizeHa
  sim$numFiresPerYear <- numFiresByPolygonNumeric / returnInterval

  sim$rstCurrentBurn <- raster(sim$fireReturnInterval) ## creates no-value raster
  sim$rstCurrentBurn[] <- 0L
  if (verbose > 0)
    message("6: ", Sys.time())

  mod$areaBurnedOverTime <- data.frame(time = numeric(0),
                                       nPixelsBurned = numeric(0),
                                       haBurned = numeric(0),
                                       FRI = numeric(0))

  mod$knownSpecies <- c(Pice_mar = "spruce", Pice_gla = "spruce",
                        Pinu_con = "pine", Pinu_ban = "pine",
                        Popu_tre = "decid", Betu_pap = "decid",
                        Abie_bal = "softwood", Abie_las = "softwood", Abie_sp = "softwood")
  sim$sppEquiv[, LandMine := mod$knownSpecies[LandR]]

  return(invisible(sim))
}

### plot events
plotFn <- compiler::cmpfun(function(sim) {
  if (time(sim) == P(sim)$.plotInitialTime) {
    friRast <- sim$fireReturnInterval
    friRast[] <- as.factor(sim$fireReturnInterval[])
    Plot(friRast, title = "Fire Return Interval", cols = c("pink", "darkred"), new = TRUE)
    sar <- sim$studyAreaReporting
    Plot(sar, addTo = "friRast", title = "", gp = gpar(col = "black", fill = 0))

    rstFlammable <- raster(sim$rstFlammable)
    rstFlammable[] <- getValues(sim$rstFlammable)
    Plot(rstFlammable, title = "Land Type (rstFlammable)", cols = c("mediumblue", "firebrick"), new = TRUE)
    Plot(sar, addTo = "rstFlammable", title = "", gp = gpar(col = "black", fill = 0))
  } else {
    firstPlot <- isTRUE(time(sim) == P(sim)$.plotInitialTime + P(sim)$.plotInterval)
    title1 <- if (firstPlot) "Current area burned (ha)" else ""
    abot <- mod$gg_areaBurnedOverTime
    Plot(abot, title = title1, new = TRUE, addTo = "areaBurnedOverTime")

    title2 <- if (firstPlot) "Cumulative Fire Map" else ""
    rcbc <- sim$rstCurrentBurnCumulative
    rcbc[!is.na(sim$rstCurrentBurn)] <- 0L
    Plot(rcbc, new = TRUE, title = title2, cols = c("pink", "red"), zero.color = "transparent")
    sar <- sim$studyAreaReporting
    Plot(sar, addTo = "rcbc", title = "", gp = gpar(col = "black", fill = 0))
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
})

### burn events
Burn <- compiler::cmpfun(function(sim, verbose = getOption("LandR.verbose", TRUE)) {
  ## DEBUGGING: random seed issues
  #fseed <- file.path(outputPath(sim), "seed.txt")
  #writeEventInfo(sim, fseed, append = TRUE)
  #writeRNGInfo(fseed, append = TRUE)
  ## END DEBUGGING

  sim$numFiresPerYear <- na.omit(sim$numFiresPerYear)
  NA_ids <- as.integer(attr(sim$numFiresPerYear, "na.action"))
  numFiresThisPeriod <- rnbinom(length(sim$numFiresPerYear),
                                mu = sim$numFiresPerYear * P(sim)$fireTimestep,
                                size = 1.3765) # Eliot lowered this from 1.8765 on Oct 23, 2018 because too constant

  thisYrStartCells <- data.table(pixel = seq(ncell(sim$fireReturnInterval)),
                                 polygonNumeric = sim$fireReturnInterval[],
                                 key = "polygonNumeric")

  thisYrStartCells <- thisYrStartCells[polygonNumeric %in% c(0, NA_ids), polygonNumeric := NA] %>%
    na.omit() %>%
    .[, SpaDES.tools:::resample(pixel, numFiresThisPeriod[.GRP]), by = polygonNumeric] %>%
    .$V1

  # If fire sizes are in hectares, must adjust based on resolution of maps
  #  NOTE: round causes fires < 0.5 pixels to NOT EXIST ... i.e., 3.25 ha fires are
  #  "not detectable" if resolution is 6.25 ha
  fireSizesThisPeriod <- VGAM::rtruncpareto(length(thisYrStartCells), lower = 1,
                                            upper = P(sim)$biggestPossibleFireSizeHa,
                                            shape = sim$kBest)

  # Because annual number of fires includes fires <6.25 ha, sometimes this will round down to 0 pixels.
  #   This calculation makes that probabilistic.
  fireSizesInPixels <- fireSizesThisPeriod / (prod(res(sim$rstFlammable)) / 1e4)
  ranDraws <- runif(length(fireSizesInPixels))
  truncVals <- trunc(fireSizesInPixels)
  decimalVals <- (unname(fireSizesInPixels - (truncVals))) > ranDraws

  fireSizesInPixels <- truncVals + decimalVals

  firesGT0 <- fireSizesInPixels > 0L
  thisYrStartCells <- thisYrStartCells[firesGT0]
  fireSizesInPixels <- fireSizesInPixels[firesGT0]

  ## Rate of Spread
  vegTypeMap <- vegTypeMapGenerator(sim$cohortData, sim$pixelGroupMap,
                                    P(sim)$vegLeadingProportion, mixedType = P(sim)$mixedType,
                                    sppEquiv = sim$sppEquiv, sppEquivCol = P(sim)$sppEquivCol,
                                    colors = sim$sppColorVect,
                                    doAssertion = P(sim)$.unitTest)

  ROSmap <- raster(sim$pixelGroupMap)
  ROSmap[] <- fireROS(sim, vegTypeMap = vegTypeMap)

  ## from DEoptim fitting - run in the LandMine.Rmd file
  spawnNewActive <- sns <- 10^c(-0.731520, -0.501823, -0.605968, -1.809726)
  spreadProb <- 0.9
  sizeCutoffs <- 10^c(2.202732, 4.696060)

  if (!all(is.na(thisYrStartCells)) & length(thisYrStartCells) > 0) {
    if (is.numeric(P(sim)$.useParallel)) {
      a <- data.table::setDTthreads(P(sim)$.useParallel)
      message("Burn should be using >100% CPU")
    } else {
      a <- data.table::setDTthreads(1)
    }
    on.exit(data.table::setDTthreads(a), add = TRUE)
    fires <- burn1(sim$fireReturnInterval,
                   startCells = thisYrStartCells,
                   fireSizes = fireSizesInPixels,
                   spreadProbRel = ROSmap,
                   sizeCutoffs = sizeCutoffs,
                   maxRetriesPerID = P(sim)$maxRetriesPerID,
                   spawnNewActive = spawnNewActive,
                   spreadProb = spreadProb)
    fa <- attr(fires, "spreadState")$clusterDT
    if (verbose > 0)
      print(fa[order(maxSize)][(.N - pmin(7, NROW(fa))):.N])

    fa1 <- fa[, list(numPixelsBurned = sum(size),
                     expectedNumBurned = sum(maxSize),
                     proportionBurned = sum(size) / sum(maxSize))]
    if (verbose > 0)
      print(fa1)

    fa[, maxSize := asInteger(maxSize)]
    sim$fireSizes[[round(time(sim) - P(sim)$burnInitialTime + 1, 0)]] <- fa[, c("size", "maxSize")]

    if (getOption("LandR.assertions", TRUE))
      if (any(tail(fa1$proportionBurned, 10)  < P(sim)$minPropBurn)) {
        mess <- "In 'LandMine' module 'Burn()': proportion area burned is less than 'minPropBurn'!"
        if (verbose > 0)
          message(crayon::red(mess))
        warning(mess, call. = FALSE)
      }

    sim$rstCurrentBurn[] <- 0L
    sim$rstCurrentBurn[fires$pixels] <- 1L #as.numeric(factor(fires$initialPixels))
  }

  if (is.null(sim$rstCurrentBurnCumulative)) {
    sim$rstCurrentBurnCumulative <- sim$rstCurrentBurn
    sim$rstCurrentBurnCumulative[!is.na(sim$rstCurrentBurnCumulative[])] <- 0
  } else {
    sim$rstCurrentBurnCumulative <- sim$rstCurrentBurn + sim$rstCurrentBurnCumulative
  }

  currBurn <- raster::mask(sim$rstCurrentBurn, sim$studyAreaReporting) %>% raster::stack()
  fris <- unique(na.omit(sim$fireReturnInterval[]))
  npix <- vapply(fris, function(x) {
    ids <- which(sim$fireReturnInterval[] == x)
    unname(table(currBurn[ids])[2])
  }, numeric(1)) %>% unname()
  npix[is.na(npix)] <- 0 # Show that zero pixels burned in a year with no pixels burned, rather than NA
  polys <- sim$fireReturnInterval
  burnedDF <- data.frame(time = as.numeric(times(sim)$current),
                         nPixelsBurned = npix,
                         haBurned = npix * prod(res(sim$rstCurrentBurn)) / 100^2, ## area in ha
                         FRI = as.factor(fris))
  mod$areaBurnedOverTime <- rbind(mod$areaBurnedOverTime, burnedDF)
  mod$gg_areaBurnedOverTime <- ggplot(mod$areaBurnedOverTime,
                                      aes(x = time, y = haBurned, fill = FRI, ymin = 0)) +
    #geom_line(size = 1.5) +
    geom_area() +
    theme(legend.text = element_text(size = 6))

  if (time(sim) == end(sim)) {
    figDir <- checkPath(file.path(outputPath(sim), "figures"), create = TRUE)
    ggsave(file.path(figDir, "areaBurnedOverTime.png"), mod$gg_areaBurnedOverTime)
  }

  return(invisible(sim))
})

.inputObjects <- function(sim) {
  ## DEBUGGING: random seed issues
  #fseed <- file.path(outputPath(sim), "seed.txt")
  #writeEventInfo(sim, fseed, append = TRUE)
  #writeRNGInfo(fseed, append = TRUE)
  ## END DEBUGGING

  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  if (getOption("LandR.verbose", TRUE) > 0)
    message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # Make random forest cover map
  nOT <- if (P(sim)$flushCachedRandomFRI) Sys.time() else NULL
  mod$numDefaultPixelGroups <- 20L
  mod$numDefaultPolygons <- 4L
  numDefaultSpeciesCodes <- 2L

  if (!suppliedElsewhere("studyArea", sim)) {
    if (getOption("LandR.verbose", TRUE) > 0)
      message("'studyArea' was not provided by user. Using a polygon in southwestern Alberta, Canada,")

    sim$studyArea <- randomStudyArea(seed = 1234, size = 1e9)
  }

  if (!suppliedElsewhere("studyAreaReporting", sim)) {
    if (getOption("LandR.verbose", TRUE) > 0)
      message("'studyAreaReporting' was not provided by user. Using the same as 'studyArea'.")
    sim$studyAreaReporting <- sim$studyArea
  }

  if (is.null(sim$rasterToMatch)) {
    if (!suppliedElsewhere("rasterToMatch", sim)) {
      sim$rasterToMatch <- raster(sim$studyArea, res = 100)
    }
  }

  if (!suppliedElsewhere("rasterToMatchReporting")) {
    sim$rasterToMatchReporting <- sim$rasterToMatch
  }

  if (!suppliedElsewhere("rstFlammable", sim)) {
    sim$rstFlammable <- sim$rasterToMatch
    sim$rstFlammable[] <- 1L  # 1 means flammable  ## TODO: use LandR::defineFlammable()
  }

  if (!suppliedElsewhere("fireReturnInterval", sim)) {
    sim$fireReturnInterval <- Cache(randomPolygons, sim$rasterToMatch,
                                    numTypes = mod$numDefaultPolygons)

    vals <- factor(sim$fireReturnInterval[],
                   levels = 1:mod$numDefaultPolygons,
                   labels = c(60, 100, 120, 250))
    sim$fireReturnInterval[] <- as.integer(as.character(vals)) ## TODO: need vals
  }

  if (!suppliedElsewhere(sim$ROSTable)) {
    sim$ROSTable <- rbindlist(list(
      list("immature_young", "decid", 6L),
      list("mature", "decid", 9L),
      list("immature_young", "mixed", 12L),
      list("mature", "mixed", 17L),
      list("immature", "pine", 14L),
      list("mature", "pine", 21L),
      list("young", "pine", 22L),
      list("immature_young", "softwood", 18L),
      list("mature", "softwood", 27L),
      list("immature_young", "spruce", 20L),
      list("mature", "spruce", 30L)
    ))
    setnames(sim$ROSTable, old = 1:3, new = c("age", "leading", "ros"))
  }

  # Upgrades to use suppliedElsewhere -- Eliot Oct 21 2018
  if (!suppliedElsewhere("pixelGroupMap", sim)) {
    sim$pixelGroupMap <- Cache(randomPolygons, sim$rasterToMatch,
                               numTypes = mod$numDefaultPixelGroups)
  }

  if (!suppliedElsewhere("rstTimeSinceFire", sim)) {
    sim$rstTimeSinceFire <- raster(sim$pixelGroupMap)
    sim$rstTimeSinceFire[] <- 200L
  }

  if (!suppliedElsewhere("species", sim)) {
    sim$species <- data.table(species = c("Pinu_sp", "Pice_gla"),
                              speciesCode = 1:numDefaultSpeciesCodes)
  }

  if (!suppliedElsewhere("cohortData", sim)) {
    sim$cohortData <- data.table(
      pixelGroup = seq(mod$numDefaultPixelGroups),
      speciesCode = factor(sample(sim$species$species, size = mod$numDefaultPixelGroups, replace = TRUE)),
      B = sample(10:20, size = mod$numDefaultPixelGroups, replace = TRUE) * 100,
      age = sample(5:20, size = mod$numDefaultPixelGroups, replace = TRUE) * 10
    )
  }

  if (!suppliedElsewhere("sppColorVect", sim)) {
    sim$sppColorVect <- c("Red", "Green")
    names(sim$sppColorVect) <- sim$species$species
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    data(sppEquivalencies_CA, envir = envir(sim))
    sim$sppEquiv <- sim$sppEquivalencies_CA
    rm("sppEquivalencies_CA", envir = envir(sim))
  }

  return(invisible(sim))
}

fireROS <- compiler::cmpfun(function(sim, vegTypeMap) {
  ROS <- rep(NA_integer_, ncell(vegTypeMap))

  vegType <- getValues(vegTypeMap)
  vegTypes <- data.table(raster::levels(vegTypeMap)[[1]]) # 2nd column in levels

  sppNames <- equivalentName(as.character(vegTypes[[2]]), sim$sppEquiv, P(sim)$sppEquivCol)
  suppressWarnings({
    onRaster <- rbindlist(list(
      list("mixed", which(is.na(sppNames))),
      list("spruce", grep(sppNames, pattern = "Pice")),
      list("pine", grep(sppNames, pattern = "Pinu")),
      list("decid", grep(sppNames, pattern = "Popu")),
      list("softwood", grep(sppNames, pattern = "Pice|Pinu|Popu", invert = TRUE))
    ))
  })
  # remove duplicates of softwood, which is NA
  onRaster <- na.omit(unique(onRaster, by = "V2"))
  setnames(onRaster, old = 1:2, new = c("leading", "pixelValue"))

  sppEquiv <- sim$sppEquiv[, c("LandMine", "LandR")][, leading := mod$knownSpecies[LandR]]
  sppEquiv <- na.omit(sppEquiv, on = "LandMine")
  sppEquiv <- unique(sppEquiv[onRaster, on = c("LandMine" = "leading")])

  sppEquivHere <- unique(na.omit(sppEquiv$LandR))
  haveAllKnown <- sppEquivHere %in% names(mod$knownSpecies)
  if (!all(haveAllKnown)) {
    stop("LandMine only has rate of spread burn rates for\n",
         paste(names(mod$knownSpecies), collapse = ", "),
         "\nMissing rate of spread for ", paste(sppEquivHere[!haveAllKnown], collapse = ", "))
  }
  sppEquiv <- unique(sppEquiv, by = c("LandMine", "leading", "pixelValue"))
  sppEquiv <- sppEquiv[sim$ROSTable, on = "leading", allow.cartesian = TRUE, nomatch = NULL]
  sppEquiv <- sppEquiv[, c("leading", "age", "ros", "pixelValue")]
  sppEquiv <- unique(sppEquiv, by = c("age", "leading", "pixelValue"))

  sppEquiv[, used := "no"]
  sppEquiv[(used == "no") & grepl("(^|_)mature", age), used := "mature"]
  sppEquiv[(used == "no") & grepl("(^|_)immature", age), used := "immature"]
  sppEquiv[(used == "no") & grepl("(^|_)young", age), used := "young"]
  setkeyv(sppEquiv, "used")

  # if there are no "mature_immature"
  cuts <- list()
  if (!any(grepl("_mature$|^mature_|_mature_", sppEquiv$age))) {
    cuts[[1]] <- sim$rstTimeSinceFire[] > 120
  } else {
    cuts[[1]] <- !is.na(sim$rstTimeSinceFire[])
  }

  if (!any(grepl("_immature$|^immature_|_immature_", sppEquiv$age))) {
    cuts[[2]] <- sim$rstTimeSinceFire[] > 40 & sim$rstTimeSinceFire[] <= 120
  } else {
    cuts[[2]] <- sim$rstTimeSinceFire[] <= 120
  }

  cuts[[3]] <- sim$rstTimeSinceFire[] <= 40

  # Now go through from mature through immature through young
  if (!all(sppEquiv["mature"]$pixelValue %in% vegTypes[[1]]))
    cuts[["mature"]] <- cuts[["mature"]] & vegType %in% sppEquiv["mature"]$pixelValue

  if (!all(sppEquiv["immature"]$pixelValue %in% vegTypes[[1]]))
    cuts[[2]] <- cuts[[2]] & vegType %in% sppEquiv["immature"]$pixelValue

  if (all(sppEquiv["young"]$pixelValue %in% vegTypes[[1]]))
    cuts[[3]] <- cuts[[3]] & vegType %in% sppEquiv["young"]$pixelValue

  mature <- which(cuts[[1]])
  immature <- which(cuts[[2]])
  young <- which(cuts[[3]])

  if (length(mature))
    ROS[mature] <- sppEquiv["mature"]$ros[match(vegType[mature], sppEquiv["mature"]$pixelValue)]
    #ROS[mature] <- plyr::mapvalues(vegType[mature], sppEquiv["mature"]$pixelValue, sppEquiv["mature"]$ros)
  if (length(immature))
    ROS[immature] <- sppEquiv["immature"]$ros[match(vegType[immature], sppEquiv["immature"]$pixelValue)]
    #ROS[immature] <- plyr::mapvalues(vegType[immature], sppEquiv["immature"]$pixelValue, sppEquiv["immature"]$ros)
  if (length(young))
    ROS[young] <- sppEquiv["young"]$ros[match(vegType[young], sppEquiv["young"]$pixelValue)]
    #ROS[young] <- plyr::mapvalues(vegType[young], sppEquiv["young"]$pixelValue, sppEquiv["young"]$ros)

  if (getOption("LandR.assertions", TRUE)) {
    names(cuts) <- c("mature", "immature", "young")
    dt <- data.table(
      ROS = ROS,
      pixelValue = vegType,
      age = cut(sim$rstTimeSinceFire[], breaks = c(0, 40, 120, 999),
                labels = c("young", "immature", "mature")),
      as.data.table(cuts)
    )
    dt <- na.omit(dt, cols = c("ROS", "age"))
    dtSumm <- dt[, list(derivedROS = unique(ROS)), by = c("pixelValue", "age")]
    dtSumm <- dtSumm[sppEquiv, on = c("pixelValue", "age" = "used"), nomatch = NULL]
    if (!(identical(dtSumm$derivedROS, dtSumm$ros))) {
      stop("fireROS failed its test")
    }
  }

  ## Other vegetation that can burn -- e.g., grasslands, lichen, shrub
  ## The original default value is the same as that of mature spruce stands (30L)
  matureSpruceROS <- sim$ROSTable[leading == "spruce" & age == "mature", ros]
  assertthat::assert_that(
    isTRUE(inRange(P(sim)$ROSother, min(sim$ROSTable$ros), max(sim$ROSTable$ros))),
    isTRUE(inRange(P(sim)$ROSother, 0.95*matureSpruceROS, 1.05*matureSpruceROS))
  )
  ROS[sim$rstFlammable[] == 1L & is.na(ROS)] <- as.integer(P(sim)$ROSother)
  ROS[sim$rstFlammable[] == 0L | is.na(sim$rstFlammable[])] <- NA ## non-flammable pixels

  return(ROS)
})
