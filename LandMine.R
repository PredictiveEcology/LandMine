defineModule(sim, list(
  name = "LandMine",
  description = "Reimplementation of Andison (1999) LandMine fire model",
  keywords = c("Fire", "Landscape", "Percolation", "Pixel-based"),
  authors = c(
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@canada.ca", role = c("aut", "cre")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@friresearch.ca", role = c("ctb"))
  ),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.3.9009", LandMine = numeric_version("0.0.1")),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "LandMine.Rmd"),
  reqdPkgs = list("data.table", "grDevices", "magrittr", "raster", "RColorBrewer", "VGAM",
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
    defineParameter("maxRetriesPerID", "integer", 10L, 0L, 20L,
                    "Minimum proportion burned pixels to use when triggering warnings about simulated fires."),
    defineParameter("ROStype", "character", "original", NA, NA,
                    paste("How to modify the 'rate of spread' parameters for different veg types.",
                          "One of 'equal', 'log', or 'original'.")),
    defineParameter("useSeed", "integer", NULL, NA, NA,
                    paste("Only used for creating a starting cohortData dataset.",
                          "If NULL, then it will be randomly generated;",
                          "If non-NULL, will pass this value to set.seed and be deterministic and identical each time.",
                          "WARNING: setting the seed to a specific value will cause all simulations to be identical!")),
    defineParameter("vegLeadingProportion", "numeric", 0.8, 0, 1,
                    "a number that define whether a species is leading for a given pixel"),
    defineParameter(".plotInitialTime", "numeric", start(sim, "year") + 1, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant"),
    defineParameter(".unitTest", "logical", TRUE, NA, NA,
                    "Some functions can have internal testing. This will turn those on or off, if any exist"),
    defineParameter(".useParallel", "numeric", parallel::detectCores(), NA, NA,
                    "Used in burning. Will be passed to data.table::setDTthreads")
  ),
  inputObjects = bind_rows(
    expectsInput("cohortData", "data.table",
                 desc = "Columns: B, pixelGroup, speciesCode, Indicating several features about ages and current vegetation of stand"),
    expectsInput("fireReturnInterval","Raster",
                 desc = "A raster layer that is a factor raster, with at least 1 column called fireReturnInterval, representing the fire return interval in years"),
    expectsInput("pixelGroupMap", "RasterLayer",
                 desc = "Pixels with identical values share identical stand features"),
    expectsInput("rasterToMatch", "RasterLayer",
                 #desc = "this raster contains two pieces of information: Full study area with fire return interval attribute",
                 desc = "DESCRIPTION NEEDED", # TODO: is this correct?
                 sourceURL = NA),
    expectsInput("rasterToMatchReporting", "RasterLayer",
                 desc = paste("Raster layer of study area used for plotting and reporting only.",
                              "Defaults to the kNN biomass map masked with `studyArea`"),
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"),
    expectsInput("rstFlammable", "Raster",
                 desc = paste("A raster layer, with 0, 1 and NA, where 1 indicates areas",
                              "that are flammable, 0 not flammable (e.g., lakes)",
                              "and NA not applicable (e.g., masked)")),
    expectsInput("rstTimeSinceFire", "Raster",
                 desc = "a time since fire raster layer", NA),
    expectsInput("species", "data.table",
                 desc = "Columns: species, speciesCode, Indicating several features about species"),
    expectsInput("sppColorVect", "character",
                 desc = "named character vector of hex colour codes corresponding to each species",
                 sourceURL = ""),
    expectsInput("studyAreaReporting", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (typically smaller/unbuffered than studyArea) to use for plotting/reporting.",
                              "Defaults to an area in Southwestern Alberta, Canada."),
                 sourceURL = NA)
  ),
  outputObjects = bind_rows(
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
      "A Raster map showing the fire return interval. THis is created from the rstCurrentBurn.")
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
    createsOutput("rstCurrentBurnCumulative", "RasterLayer", "Cumulative number of times a pixel has burned")
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

    mod$LandMineDevice <- max(dev.list()) + 1

    # schedule future event(s)
    sim <- scheduleEvent(sim, P(sim)$burnInitialTime, "LandMine", "Burn", 2.5)
    sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "LandMine", "plot")
    sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "LandMine", "save")
  } else if (eventType == "plot") {
    # ! ----- EDIT BELOW ----- ! #
    # do stuff for this event

    devCur <- dev.cur()
    quickPlot::dev(mod$LandMineDevice, width= 18, height = 12)
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
    fs <- round(rtruncpareto(1e6, 1, upper = upper1, shape = params[1]))
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
  sim$fireSizes <- list()

  if (!suppliedElsewhere("cohortData", sim)) {
    if (!is.null(P(sim)$useSeed)) {
      set.seed(P(sim)$useSeed)
    }

    sampleV <- Vectorize(sample, "size", SIMPLIFY = TRUE)
    repV <- Vectorize(rep.int, c("x","times"))
    numCohortsPerPG <- sample(1:2, replace = TRUE, mod$numDefaultPixelGroups)
    sim$cohortData <- data.table(speciesCode = unlist(sampleV(1:2, numCohortsPerPG)),
                                 B = runif(sum(numCohortsPerPG), 100, 1000),
                                 pixelGroup = unlist(repV(1:mod$numDefaultPixelGroups, times = numCohortsPerPG)))
  }

  if (verbose > 0)
    message("Initializing fire maps")
  sim$fireTimestep <- P(sim)$fireTimestep
  sim$fireInitialTime <- P(sim)$burnInitialTime

  # check sim$fireReturnInterval should have no zeros
  zeros <- sim$fireReturnInterval[] == 0
  if (any(zeros, na.rm = TRUE)) sim$fireReturnInterval[zeros] <- NA
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

  if (verbose > 0)
    message("Write fire return interval map to disk")

  #sim$fireReturnInterval <- raster(sim$rasterToMatch)
  #sim$fireReturnInterval[] <- raster::factorValues(sim$rasterToMatch, sim$rasterToMatch[], att = "fireReturnInterval")[, 1]
  #fireReturnIntFilename <- file.path(cachePath(sim), "rasters/fireReturnInterval.tif")
  #sim$fireReturnInterval <- writeRaster(sim$fireReturnInterval, filename = fireReturnIntFilename,
  #                                      datatype = "INT2U", overwrite = TRUE)

  sim$rstCurrentBurn <- raster(sim$fireReturnInterval) ## creates no-value raster
  sim$rstCurrentBurn[] <- 0L
  if (verbose > 0)
    message("6: ", Sys.time())

  mod$areaBurnedOverTime <- data.frame(time = numeric(0),
                                       nPixelsBurned = numeric(0),
                                       haBurned = numeric(0),
                                       FRI = numeric(0))
  return(invisible(sim))
}

### plot events
plotFn <- function(sim) {
  if (is.null(sim$rstCurrentBurnCumulative)) {
    sim$rstCurrentBurnCumulative <- raster(sim$rstCurrentBurn)
  }
  if (time(sim) == P(sim)$.plotInitialTime) {
    friRast <- sim$fireReturnInterval
    friRast[] <- as.factor(sim$fireReturnInterval[])
    Plot(friRast, title = "Fire Return Interval", cols = c("pink", "darkred"), new = TRUE)
    sar <- sim$studyAreaReporting
    Plot(sar, addTo = "friRast", title = "",
         gp = gpar(col = "black", fill = 0))

    sim$rstCurrentBurnCumulative[!is.na(sim$rstCurrentBurn)] <- 0L

    rstFlammable <- raster(sim$rstFlammable)
    rstFlammable[] <- getValues(sim$rstFlammable)
    Plot(rstFlammable, title = "Land Type (rstFlammable)", cols = c("mediumblue", "firebrick"), new = TRUE)
    Plot(sar, addTo = "rstFlammable", title = "",
         gp = gpar(col = "black", fill = 0))
  }

  currBurn <- raster::mask(sim$rstCurrentBurn, sim$studyAreaReporting) %>% stack()
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

  if (length(unique(mod$areaBurnedOverTime$time)) > 1) {

    gg_areaBurnedOverTime <- ggplot(mod$areaBurnedOverTime,
                                    aes(x = time, y = haBurned, fill = FRI, ymin = 0)) +
      #geom_line(size = 1.5) +
      geom_area() +
      theme(legend.text = element_text(size = 6))

    firstPlot <- isTRUE(time(sim) == P(sim)$.plotInitialTime + P(sim)$.plotInterval)
    title1 <- if (firstPlot) "Current area burned (ha)" else ""
    Plot(gg_areaBurnedOverTime, title = title1, new = TRUE, addTo = "areaBurnedOverTime")

    sim$rstCurrentBurnCumulative <- sim$rstCurrentBurn + sim$rstCurrentBurnCumulative
    title2 <- if (firstPlot) "Cumulative Fire Map" else ""
    rcbc <- sim$rstCurrentBurnCumulative
    Plot(rcbc, new = TRUE,
         title = title2,
         cols = c("pink", "red"), zero.color = "transparent")
    sar <- sim$studyAreaReporting
    Plot(sar, addTo = "rcbc", title = "",
         gp = gpar(col = "black", fill = 0))
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### burn events
Burn <- function(sim, verbose = getOption("LandR.verbose", TRUE)) {
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

  thisYrStartCells <- data.table(pixel = 1:ncell(sim$fireReturnInterval),
                                 polygonNumeric = sim$fireReturnInterval[],
                                 key = "polygonNumeric")

  thisYrStartCells <- thisYrStartCells[polygonNumeric %in% c(0, NA_ids), polygonNumeric := NA] %>%
    na.omit() %>%
    .[, SpaDES.tools:::resample(pixel, numFiresThisPeriod[.GRP]), by = polygonNumeric] %>%
    .$V1

  # If fire sizes are in hectares, must adjust based on resolution of maps
  #  NOTE: round causes fires < 0.5 pixels to NOT EXIST ... i.e., 3.25 ha fires are
  #  "not detectable" if resolution is 6.25 ha
  fireSizesThisPeriod <- rtruncpareto(length(thisYrStartCells), lower = 1,
                                      upper = P(sim)$biggestPossibleFireSizeHa,
                                      shape = sim$kBest)

  # Because annual number of fires includes fires <6.25 ha, sometimes this will round down to 0 pixels.
  #   This calculation makes that probabilistic.
  fireSizesInPixels <- fireSizesThisPeriod / (prod(res(sim$rstFlammable)) / 1e4)
  ranDraws <- runif(length(fireSizesInPixels))
  truncVals <- trunc(fireSizesInPixels)
  decimalVals <- (unname(fireSizesInPixels - (truncVals))) > ranDraws

  fireSizesInPixels <- truncVals + decimalVals

  firesGT0 <- fireSizesInPixels > 0
  thisYrStartCells <- thisYrStartCells[firesGT0]
  fireSizesInPixels <- fireSizesInPixels[firesGT0]

  ## Rate of Spread
  vegTypeMap <- vegTypeMapGenerator(sim$cohortData, sim$pixelGroupMap,
                                    P(sim)$vegLeadingProportion,
                                    colors = sim$sppColorVect,
                                    unitTest = P(sim)$.unitTest)

  ROSmap <- raster(sim$pixelGroupMap)
  ROSmap[] <- fireROS(sim, type = P(sim)$ROStype, vegTypeMap = vegTypeMap)

  # From DEoptim fitting - run in the LandMine.Rmd file
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
                   #spawnNewActive = c(0.65, 0.6, 0.2, 0.2),
                   sizeCutoffs = sizeCutoffs,
                   maxRetriesPerID = P(sim)$maxRetriesPerID,
                   spawnNewActive = spawnNewActive,
                   #spawnNewActive = c(0.76, 0.45, 1.0, 0.00),
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

  return(invisible(sim))
}

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
  emptyRas <- raster(extent(0, 2e4, 0, 2e4), res = 250)

  if (is.null(sim$rasterToMatch)) {
    if (!suppliedElsewhere("rasterToMatch", sim)) {
      stop("There is no 'rasterToMatch' supplied")
    }
  }

  if (!suppliedElsewhere("rasterToMatchReporting")) {
    sim$rasterToMatchReporting <- sim$rasterToMatch
  }

  if (!suppliedElsewhere("studyArea", sim)) {
    if (getOption("LandR.verbose", TRUE) > 0)
      message("'studyArea' was not provided by user. Using a polygon in southwestern Alberta, Canada,")

    sim$studyArea <- randomStudyArea(seed = 1234)
  }

  if (!suppliedElsewhere("studyAreaReporting", sim)) {
    if (getOption("LandR.verbose", TRUE) > 0)
      message("'studyAreaReporting' was not provided by user. Using the same as 'studyArea'.")
    sim$studyAreaReporting <- sim$studyArea
  }

  if (!suppliedElsewhere("rstFlammable", sim)) {
    sim$rstFlammable <- sim$rasterToMatch
    sim$rstFlammable[] <- 1L  # 1 means flammable
  }

  if (!suppliedElsewhere("fireReturnInterval", sim)) {
    sim$fireReturnInterval <- Cache(randomPolygons, emptyRas,
                                    numTypes = mod$numDefaultPolygons, notOlderThan = nOT,
                                    cacheRepo = cachePath(sim))

    #vals <- factor(sim$fireReturnInterval[],
    #               levels = 1:mod$numDefaultPolygons,
    #               labels = c(60, 100, 120, 250))
    sim$fireReturnInterval[] <- as.integer(as.character(vals)) ## TODO: need vals
  }

  # Upgrades to use suppliedElsewhere -- Eliot Oct 21 2018
  if (!suppliedElsewhere("pixelGroupMap", sim)) {
    sim$pixelGroupMap <- Cache(randomPolygons, emptyRas, numTypes = mod$numDefaultPixelGroups,
                               notOlderThan = nOT, cacheRepo = cachePath(sim))
  }

  if (!suppliedElsewhere("rstTimeSinceFire", sim)) {
    #if (is.null(sim$rstTimeSinceFire)) {
    sim$rstTimeSinceFire <- raster(sim$pixelGroupMap)
    sim$rstTimeSinceFire[] <- 200L
  }

  if (!suppliedElsewhere("species", sim)) {
    #if (is.null(sim$species)) {
    sim$species <- data.table(species = c("Pinu_sp", "Pice_gla"),
                              speciesCode = 1:numDefaultSpeciesCodes)
  }

  # if (!suppliedElsewhere("rstCurrentBurnCumulative))", sim)) {
  #   #if (is.null(sim$rstCurrentBurnCumulative)) {
  #   sim$rstCurrentBurnCumulative <- raster(sim$pixelGroupMap)
  #   sim$rstCurrentBurnCumulative[sim$rstTimeSinceFire[] == 0] <- 1
  # }

  # see https://github.com/PredictiveEcology/SpaDES.tools/issues#17 for discussion about this
  # meta <- depends(sim)@dependencies
  # mods <- unlist(modules(sim))
  # if (all(names(meta) %in% mods)) {
  #   # means there is more than just this module in the simList
  #
  #   outputs <- lapply(meta, function(x) {x@outputObjects$objectName})
  #   otherMods <- mods[!(mods %in% currentModule(sim))]
  #
  #   # is it or will it be supplied by another module, if yes, don't load a default here
  #   if (!("rstCurrentBurnCumulative" %in% unlist(outputs[otherMods]))) {
  #     if (is.null(sim$rstCurrentBurnCumulative)) {
  #       sim$rstCurrentBurnCumulative <- raster(sim$pixelGroupMap)
  #     }
  #   }
  # }

  return(invisible(sim))
}

fireROS <- function(sim, type = "original", vegTypeMap) {
  vegType <- getValues(vegTypeMap)
  vegTypes <- data.frame(raster::levels(vegTypeMap)[[1]][, 2, drop = FALSE]) # 2nd column in levels
  #vegTypes <- factorValues(vegTypeMap, seq_len(NROW(levels(vegTypeMap)[[1]]))) # [vegType, "Factor"]

  ## note: these are defined differently than in LandWeb, and that's ok?
  mature <- sim$rstTimeSinceFire[] > 120
  immature <- (sim$rstTimeSinceFire[] > 40) & !mature
  young <- !immature & !mature

  ROS <- rep(NA_integer_, NROW(vegType))
  mixed <- grep(tolower(vegTypes$Factor), pattern = "mix")
  spruce <- grep(tolower(vegTypes$Factor), pattern = "spruce")
  pine <- grep(tolower(vegTypes$Factor), pattern = "pine")
  decid <- grep(tolower(vegTypes$Factor), pattern = "deci")
  softwood <- grep(tolower(vegTypes$Factor), pattern = "soft")

  ROS[!mature & vegType %in% decid] <- 6L
  ROS[mature & vegType %in% decid] <- 9L

  ROS[!mature & vegType %in% mixed] <- 12L
  ROS[mature & vegType %in% mixed] <- 17L

  ROS[immature & vegType %in% pine] <- 14L
  ROS[mature & vegType %in% pine] <- 21L
  ROS[young & vegType %in% pine] <- 22L

  ROS[!mature & vegType %in% softwood] <- 18L
  ROS[mature & vegType %in% softwood] <- 27L

  ROS[!mature & vegType %in% spruce] <- 20L
  ROS[mature & vegType %in% spruce] <- 30L

  # Other vegetation that can burn -- e.g., grasslands, lichen, shrub
  ROS[sim$rstFlammable[] == 1L & is.na(ROS)] <- 30L

  if (type == "equal") {
    ## equal rates of spread
    ROS[young & vegType %in% c(mixed, spruce, pine, decid, softwood)] <- 1L
    ROS[immature & vegType %in% c(mixed, spruce, pine, decid, softwood)] <- 1L
    ROS[mature & vegType %in% c(mixed, spruce, pine, decid, softwood)] <- 1L
    ROS[sim$rstFlammable[] == 1L & is.na(ROS)] <- 1L
  }

  ## log(rates of spread), which maintains relationships but makes more equal
  if (type == "log") {
    ROS[] <- log(ROS[])
  }

  ROS
}
