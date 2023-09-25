defineModule(sim, list(
  name = "LandMine",
  description = "Reimplementation of Andison (1999) LandMine fire model",
  keywords = c("Fire", "Landscape", "Percolation", "Pixel-based"),
  authors = c(
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  childModules = character(0),
  version = list(LandMine = numeric_version("0.0.3")),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "LandMine.Rmd"),
  reqdPkgs = list("assertthat", "data.table", "fpCompare", "ggplot2", "grDevices", "gridExtra",
                  "magrittr", "raster", "RColorBrewer", "stats", "VGAM",
                  "quickPlot", "fasterize",
                  "PredictiveEcology/LandR@development (>= 1.1.0.9003)",
                  "PredictiveEcology/LandWebUtils@development (>= 0.1.7)",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/SpaDES.tools@development"),
  parameters = rbind(
    defineParameter("biggestPossibleFireSizeHa", "numeric", 1e6, 1e4, 2e6,
                    "An upper limit, in hectares, of the truncated Pareto distribution of fire sizes"),
    defineParameter("burnInitialTime", "numeric", start(sim, "year") + 1, NA, NA,
                    "This describes the simulation time at which the first burn event should occur"),
    defineParameter("fireTimestep", "numeric", 1, NA, NA,
                    "This describes the simulation time interval between burn events"),
    defineParameter("maxReburns", "integer", c(1L, 20L), 1L, 20L,
                    paste("Number of attempts to burn fires that don't reach their target fire size.",
                          "Reburning occurs in two phases, hence accepting a parameter value of length 2.",
                          "In the first phase, fires that did not reach their target size are reignited",
                          "from new pixels within the FRI zone, so they are less likely to continue",
                          "being stuck in a region with sinuous fires or discontinuous fuels.",
                          "If, after `maxReburns[1]` attempts, there are still fires that haven't reached",
                          "their target size, the second reburn phase is attempted.",
                          "After recording the the pixels that *did* burn in phase one,",
                          "*new* fires are ignited, whose target sizes are set equal the difference",
                          "between the previous target and the previously burned area. Repeats up to `maxReburns[2]` times.",
                          "This results in additional (smaller) fires, but since the purpose of LandMine",
                          "is to replicate area burned per year to achieve LTHFC, this is an acceptable compromise.")),
    defineParameter("maxRetriesPerID", "integer", 4L, 0L, 299L,
                    paste("Number of additional attempts ('jumps') that will be made per firelet ID, before abandoning.",
                          "See `?SpaDES.tools::spread2`.",
                          "NOTE: increasing this value results in longer simulation times when firelets get 'stuck',",
                          "but higher values are needed to achive larger fire sizes with discontinuous fuels.")),
    defineParameter("minPropBurn", "numeric", 0.90, 0.00, 1.00,
                    "Minimum proportion burned pixels to use when triggering warnings about simulated fires."),
    defineParameter("mixedType", "numeric", 2, 1, 2,
                    paste("How to define mixed stands: 1 for any species admixture;",
                          "2 for deciduous > conifer. See ?vegTypeMapGenerator.")),
    defineParameter("mode", "character", "single", NA, NA,
                    paste("use 'single' to run part of a landscape simulation;",
                          "use 'multi' to run as part of postprocessing multiple simulation runs.")),
    defineParameter("optimParsRowID", "integer", 1L, 1L, NA,
                    paste("which set of optimization parameter values to use for simulating fire spread,",
                          "specified by row number of the `LandMine_DEoptim_params.csv` file.",
                          "`1L` specifies the original 2018 values at 100m pixels;",
                          "all other rows were calculated using 250m pixels.")),
    defineParameter("reps", "integer", NA_integer_, 1L, NA_integer_,
                    paste("number of replicates/runs per study area when running in 'multi' mode.")),
    defineParameter("ROSother", "integer", 30L, NA, NA,
                    paste0("default ROS value for non-forest vegetation classes.",
                           "this is needed when passing a modified `ROSTable`, e.g. using log-transformed values.")),
    defineParameter("ROStype", "character", "default", NA, NA,
                    "One of 'burny', 'equal', 'log', or 'default'."),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    "The column in `sim$specieEquivalency` data.table to use as a naming convention."),
    defineParameter("useSeed", "integer", NULL, NA, NA,
                    paste("Only used for creating a starting `cohortData` dataset.",
                          "If `NULL`, then it will be randomly generated;",
                          "If non-`NULL`, will pass this value to `set.seed()` and be deterministic and identical each time.",
                          "WARNING: setting the seed to a specific value will cause all simulations to be identical!")),
    defineParameter("vegLeadingProportion", "numeric", 0.8, 0, 1,
                    "a number that define whether a species is leading for a given pixel"),
    defineParameter(".plotInitialTime", "numeric", start(sim, "year") + 1, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".plots", "character", c("png", "screen"), NA, NA,
                    paste("Passed to `types` in `Plots` (see `?Plots`). There are a few plots that are made within this module, if set.",
                          "Note that plots (or their data) saving will ONLY occur at `end(sim)`.",
                          "If `NA`, plotting is turned off completely (this includes plot saving).")),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events"),
    defineParameter(".studyAreaName", "character", "test", NA, NA,
                    "Human-readable name for the study area used - e.g., a hash of the study",
                    "area obtained using `reproducible::studyAreaName()`"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should this entire module be run with caching activated?",
                          "This is generally intended for data-type modules,",
                          "where stochasticity and time are not relevant")),
    defineParameter(".unitTest", "logical", getOption("LandR.assertions", TRUE), NA, NA,
                    "Some functions can have internal testing. This will turn those on or off, if any exist."),
    defineParameter(".useParallel", "numeric", 2, NA, NA,
                    paste("Used in burning. Will be passed to `data.table::setDTthreads()`.",
                          "NOTE: use `.useParallel <= 2` as the additonal RAM overhead too high given marginal speedup."))
  ),
  inputObjects = bindrows(
    expectsInput("cohortData", "data.table",
                 desc = paste("Columns: B, pixelGroup, speciesCode (as a factor of the names), age.",
                              "indicating several features about the current vegetation of stand."),
                 sourceURL = NA),
    expectsInput("fireReturnInterval", "Raster",
                 desc = paste("A raster layer that is a factor raster, with at least 1 column called",
                              "`fireReturnInterval`, representing the fire return interval in years."),
                 sourceURL = NA),
    expectsInput("pixelGroupMap", "RasterLayer",
                 desc = "Pixels with identical values share identical stand features",
                 sourceURL = NA),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = paste("a raster of the `studyArea` to use as a template raster",
                              "(resolution, projection, etc.) for all other rasters in the simulation."),
                 sourceURL = NA),
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
                              "Default taken from `LandR::sppEquivalencies_CA` which has names for",
                              "species of trees in Canada"),
                 sourceURL = NA),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon, typically buffered around an area of interest",
                              "(i.e., `studyAreaReporting`) to use for simulation.",
                              "Defaults to an area in Southwestern Alberta, Canada."),
                 sourceURL = NA),
    expectsInput("studyAreaReporting", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (typically smaller/unbuffered than `studyArea`)",
                              "to use for plotting/reporting.",
                              "Defaults to an area in Southwestern Alberta, Canada."),
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput("fireInitialTime", "numeric", paste(
      "The initial event time of the burn event.",
      "This is simply a reassignment from `P(sim)$burnInitialTime`.")
    ),
    createsOutput("fireSizes", "list", paste(
      "A list of data.tables, one per burn event, each with two columns, `size` and `maxSize`.",
      "These indicate the actual sizes and expected sizes burned, respectively.",
      "These can be put into a single data.table with `rbindlist(sim$fireSizes, idcol = 'year')`")
    ),
    createsOutput("fireReturnInterval", "RasterLayer", paste(
      "A `Raster` map showing the fire return interval. This is created from the `rstCurrentBurn`.")
    ),
    createsOutput("fireReturnIntervalsByPolygonNumeric", "numeric", paste(
      "A vector of the fire return intervals, ordered by the numeric representation of polygon ID")
    ),
    createsOutput("fireTimestep", "numeric", paste(
      "The number of time units between successive fire events in a fire module.")
    ),
    createsOutput("friSummary", "data.table", "summary fire return interval table"),
    createsOutput("kBest", "numeric", paste(
      "A numeric scalar that is the optimal value of `K` in the",
      "Truncated Pareto distribution (`rtruncpareto`)")
    ),
    createsOutput("numFiresPerYear", "numeric", paste(
      "The average number of fires per year, by fire return interval level on `rstCurrentBurn`.")
    ),
    createsOutput("rstCurrentBurn", "RasterLayer", paste(
      "A raster layer, produced at each timestep, where each",
      "pixel is either 1 or 0 indicating burned or not burned.")
    ),
    createsOutput("rstCurrentBurnCumulative", "RasterLayer", "Cumulative number of times a pixel has burned"),
    createsOutput("sppEquiv", "data.table", paste("Same as input, but with new column, `LandMine`."))
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
    if (P(sim)$mode == "single") {
      sim <- scheduleEvent(sim, P(sim)$burnInitialTime, "LandMine", "Burn", 2.5)
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "LandMine", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "LandMine", "save")
      sim <- scheduleEvent(sim, end(sim), "LandMine", "summarySingle")
    } else if (P(sim)$mode == "multi") {
      sim <- scheduleEvent(sim, start(sim), "LandMine", "summaryMulti")
    }
  } else if (eventType == "plot") {
    ## TODO: allow plot to file
    if (anyPlotting(P(sim)$.plots) && any(P(sim)$.plots == "screen")) {

      if (is.null(mod$LandMineDevice)) {
        dl <- dev.list()
        quickPlot::dev.useRSGD(FALSE)
        # if the device was already "this" size, meaning probably made here
        needDev <- FALSE
        desiredDims <- c(width = 14.3, height = 9.5)
        if (is.null(dl)) {
          needDev <- TRUE
        } else {
          if (!all(abs(dev.size() - desiredDims) < 0.5))
            needDev <- TRUE
        }
        if (needDev) {
          newDev <- max(dl) + 1
          do.call(quickPlot::dev, append(list(newDev), as.list(desiredDims)))
        }
        mod$LandMineDevice <- dev.cur()
      }
      quickPlot::dev(mod$LandMineDevice)
      sim <- plotFn(sim)

      sim <- scheduleEvent(sim, P(sim)$.plotInterval, "LandMine", "plot")
    }
  } else if (eventType == "Burn") {
    sim <- Burn(sim)
    sim <- scheduleEvent(sim, time(sim) + P(sim)$fireTimestep, "LandMine", "Burn", 2.5)
  } else if (eventType == "summarySingle") {
    sim <- SummarizeFRIsingle(sim)
  } else if (eventType == "summaryMulti") {
    sim <- SummarizeFRImulti(sim)
  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}

### initialization
EstimateTruncPareto <- function(sim, verbose = getOption("LandR.verbose", TRUE)) {
  if (verbose > 0) {
    message("Estimate Truncated Pareto parameters")
  }

  findK_upper <- function(params = c(0.4), upper1) {
    fs <- round(VGAM::rtruncpareto(1e6, 1, upper = upper1, shape = params[1]))
    # meanFS <- meanTruncPareto(k = params[1], lower = 1, upper = upper1, alpha = 1)
    # diff1 <- abs(quantile(fs, 0.95) - meanFS)

    ## "90% of area is in 5% of fires" - Dave rule of thumb
    # abs(sum(fs[fs>quantile(fs, 0.95)])/sum(fs) - 0.95)

    ## Eliot's adjustment because each year was too constant; should create greater variation.
    abs(sum(fs[fs > quantile(fs, 0.95)]) / sum(fs) - 0.95) ## "95% of area (2nd term) is in 5% of fires (1st term)"

    ## 2018-110-23: Eliot's adjustment because each year still too constant; need greater variation.
    # abs(sum(fs[fs > quantile(fs, 0.90)]) / sum(fs) - 0.95) # "95% of area (2nd term) is in 10% of fires (1st term)"
  }

  sim$kBest <- Cache(
    optimize,
    interval = c(0.05, 0.99),
    f = findK_upper,
    upper1 = P(sim)$biggestPossibleFireSizeHa,
    cacheRepo = cachePath(sim),
    useCache = FALSE
  )$minimum

  return(invisible(sim))
}

Init <- function(sim, verbose = getOption("LandR.verbose", TRUE)) {
  ## DEBUGGING: random seed issues
  #fseed <- file.path(outputPath(sim), "seed.txt")
  #writeEventInfo(sim, fseed, append = TRUE)
  #writeRNGInfo(fseed, append = TRUE)
  ## END DEBUGGING

  if (is.null(P(sim)$maxReburns) || is.null(P(sim)$maxRetriesPerID)) {
    stop("maxReburns and maxRetries must be integer values and cannot be NULL.")
  }

  if (length(P(sim)$maxReburns) == 1) {
    P(sim, "maxReburns", "LandMine") <- rep(P(sim)$maxReburns, 2)
  }

  P(sim, "maxReburns", "LandMine") <- as.integer(P(sim, "maxReburns", "LandMine"))
  P(sim, "maxRetriesPerID", "LandMine") <- as.integer(P(sim, "maxRetriesPerID", "LandMine"))

  compareRaster(sim$rasterToMatch, sim$fireReturnInterval, sim$rstFlammable, sim$rstTimeSinceFire)

  ## from DEoptim fitting, run in the LandMine.Rmd file
  optimPars <- read.csv(file.path(dataPath(sim), "LandMine_DEoptim_params.csv"))
  optimPars <- optimPars[P(sim)$optimParsRowID, grepl("^par", colnames(optimPars))]
  optimPars <- unlist(unname(optimPars))

  mod$spawnNewActive <- 10^c(optimPars[1], optimPars[2], optimPars[3], optimPars[4])
  mod$sizeCutoffs <- 10^c(optimPars[5], optimPars[6])

  mod$spreadProb <- !is.na(sim$fireReturnInterval)
  mod$spreadProb[mod$spreadProb[] == 0] <- NA_real_
  mod$spreadProb[mod$spreadProb[] == 1] <- optimPars[7]

  sim$fireSizes <- list()

  if (!is.integer(sim$fireReturnInterval[]))
    sim$fireReturnInterval[] <- as.integer(sim$fireReturnInterval[])

  if (verbose > 0)
    message("Initializing fire maps")
  sim$fireTimestep <- P(sim)$fireTimestep
  sim$fireInitialTime <- P(sim)$burnInitialTime

  ## fireReturnInterval should have no zeros
  zeros <- sim$fireReturnInterval[] == 0L
  if (any(zeros, na.rm = TRUE)) {
    sim$fireReturnInterval[zeros] <- NA_integer_
  }

  ## 2023-09: exclude non-flammable pixels for FRI calculations
  nonFlammable <- which(is.na(sim[["rstFlammable"]][]) | sim[["rstFlammable"]][] == 0)
  if (length(nonFlammable) > 0) {
    sim$fireReturnInterval[nonFlammable] <- NA
  }

  numPixelsPerPolygonNumeric <- Cache(freq, sim$fireReturnInterval, useNA = "no", cacheRepo = cachePath(sim)) |>
    na.omit()
  colnames(numPixelsPerPolygonNumeric) <- c("fri", "count")
  numPixelsPerPolygonNumeric <- cbind(value = seq_len(NROW(numPixelsPerPolygonNumeric)), numPixelsPerPolygonNumeric)
  ordPolygons <- order(numPixelsPerPolygonNumeric[, "value"])
  numPixelsPerPolygonNumeric <- numPixelsPerPolygonNumeric[ordPolygons, , drop = FALSE]
  sim$fireReturnIntervalsByPolygonNumeric <- numPixelsPerPolygonNumeric[, "fri"]
  numPixelsPerPolygonNumeric <- numPixelsPerPolygonNumeric[, "count"]
  names(numPixelsPerPolygonNumeric) <- sim$fireReturnIntervalsByPolygonNumeric

  numHaPerPolygonNumeric <- numPixelsPerPolygonNumeric * (prod(res(sim$fireReturnInterval)) / 1e4)
  returnInterval <- sim$fireReturnIntervalsByPolygonNumeric

  if (verbose > 0) {
    message("Determine mean fire size...")
  }
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
    Plot(sar, addTo = "friRast", title = "", cols = "transparent")

    rstFlammable <- raster(sim$rstFlammable)
    rstFlammable[] <- getValues(sim$rstFlammable)
    Plot(rstFlammable, title = "Land Type (rstFlammable)", cols = c("mediumblue", "firebrick"), new = TRUE)
    Plot(sar, addTo = "rstFlammable", title = "", cols = "transparent")
  } else {
    firstPlot <- isTRUE(time(sim) == P(sim)$.plotInitialTime + P(sim)$.plotInterval)
    title1 <- if (firstPlot) "Current area burned (ha)" else ""
    abot <- mod$gg_areaBurnedOverTime
    Plot(abot, title = title1, new = TRUE, addTo = "areaBurnedOverTime")

    title2 <- if (firstPlot) "Cumulative Fire Map" else ""
    rcbc <- sim$rstCurrentBurnCumulative
    rcbc[!is.na(sim$rstCurrentBurn) & sim$rstCurrentBurn == 0 &
           sim$rstCurrentBurnCumulative == 0] <- 0L
    Plot(rcbc, new = TRUE, title = title2, cols = c("pink", "red"), zero.color = "transparent")
    sar <- sim$studyAreaReporting
    Plot(sar, addTo = "rcbc", title = "", cols = "transparent")
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
  thisYrStartCellsDT <- data.table(pixel = seq(ncell(sim$fireReturnInterval)),
                                   polygonNumeric = sim$fireReturnInterval[],
                                   key = "polygonNumeric")

  ## August 2022: reburn fires that did not meet their target size

  ## Rate of Spread
  vegTypeMap <- vegTypeMapGenerator(sim$cohortData,
                                    pixelGroupMap = sim$pixelGroupMap,
                                    vegLeadingProportion = P(sim)$vegLeadingProportion,
                                    mixedType = P(sim)$mixedType,
                                    sppEquiv = sim$sppEquiv,
                                    sppEquivCol = P(sim)$sppEquivCol,
                                    colors = sim$sppColorVect,
                                    doAssertion = P(sim)$.unitTest)

  ROSmap <- raster(sim$pixelGroupMap)
  ROSmap[] <- fireROS(sim, vegTypeMap = vegTypeMap)

  spreadProbThisStep <- mod$spreadProb

  ## If fire sizes are in hectares, must adjust based on resolution of maps
  ##  NOTE: round causes fires < 0.5 pixels to NOT EXIST ... i.e., 3.25 ha fires are
  ##  "not detectable" if resolution is 6.25 ha
  fireSizesThisPeriod <- VGAM::rtruncpareto(sum(numFiresThisPeriod), lower = 1,
                                            upper = P(sim)$biggestPossibleFireSizeHa,
                                            shape = sim$kBest)

  ## Because annual number of fires includes fires <6.25 ha, sometimes this will round down to 0 pixels.
  ##   This calculation makes that probabilistic.
  fireSizesInPixels <- fireSizesThisPeriod / (prod(res(sim$rstFlammable)) / 1e4)
  ranDraws <- runif(length(fireSizesInPixels))
  truncVals <- trunc(fireSizesInPixels)
  decimalVals <- (unname(fireSizesInPixels - (truncVals))) > ranDraws

  fireSizesInPixels <- truncVals + decimalVals

  firesList <- fireSizes <- list()
  maxOrder <- 0L
  iter <- 1L

  ## 2023-09: after maxReburns, if not reaching fire size, take the last burn,
  ## and start new fire(s) to burn the remaining area until the target is achieved.
  ## Should be OK b/c LandMine replicates FRIs (i.e., area burned each year), not number of fires
  while (sum(numFiresThisPeriod) > 0 && (iter <= sum(P(sim)$maxReburns))) {
    thisYrStartCells <- thisYrStartCellsDT[polygonNumeric %in% c(0, NA_ids), polygonNumeric := NA] %>%
      na.omit() %>%
      .[, SpaDES.tools:::resample(pixel, numFiresThisPeriod[.GRP]), by = polygonNumeric] %>%
      .$V1

    firesGT0 <- fireSizesInPixels > 0L
    thisYrStartCells <- thisYrStartCells[firesGT0]
    fireSizesInPixels <- fireSizesInPixels[firesGT0]

    if (!all(is.na(thisYrStartCells)) && length(thisYrStartCells) > 0) {
      if (iter > 1 && iter <= P(sim)$maxReburns[1]) {
        message("Some fires did not reach their target size; reburning these fires (", iter, "/", P(sim)$maxReburns[1], ")")
      } else if (iter > P(sim)$maxReburns[1] && iter <= sum(P(sim)$maxReburns)) {
        message("Some fires did not reach their target size; starting additional fires (",
                iter - P(sim)$maxReburns[1], "/", P(sim)$maxReburns[2], ")")
      }

      if (is.numeric(P(sim)$.useParallel)) {
        a <- data.table::setDTthreads(P(sim)$.useParallel)
        message(sprintf("Burn should be using >100%% CPU (useParallel = %s)", as.character(P(sim)$.useParallel)))
      } else {
        a <- data.table::setDTthreads(1)
      }
      on.exit(data.table::setDTthreads(a), add = TRUE)

      fires <- landmine_burn1(sim$fireReturnInterval,
                              startCells = thisYrStartCells,
                              fireSizes = fireSizesInPixels,
                              spreadProbRel = ROSmap,
                              sizeCutoffs = mod$sizeCutoffs,
                              maxRetriesPerID = P(sim)$maxRetriesPerID,
                              spawnNewActive = mod$spawnNewActive,
                              spreadProb = spreadProbThisStep)

      ## occasionally, `order` col drops from fires, but it's not supposed to (SpaDES.tools#74)
      if (!"order" %in% colnames(fires)) {
        fires[, order := 1:nrow(fires)]
      }
      fires[, order := order + maxOrder]

      ## occasionally, `numNeighs` col appears in fires, but it's not supposed to (SpaDES.tools#74)
      if ("numNeighs" %in% colnames(fires)) {
        set(fires, NULL, "numNeighs", NULL)
      }

      fa <- attr(fires, "spreadState")$clusterDT
      fa1 <- fa[, list(numPixelsBurned = sum(size),
                       expectedNumBurned = sum(maxSize),
                       proportionBurned = sum(size) / sum(maxSize))]

      if (verbose > 0) {
        print(fa[order(maxSize)][(.N - pmin(7, NROW(fa))):.N])
        print(fa1)
      }

      fa[, maxSize := asInteger(maxSize)]

      tooSmall <- which(fa$size != fa$maxSize)
      if (length(tooSmall)) {
        tooSmallDT <- fa[tooSmall, c("initialPixels", "maxSize")]
        tooSmallByPoly <- thisYrStartCellsDT[tooSmallDT, on = c(pixel = "initialPixels")]
        friByPolyDT <- data.table(polygonNumeric = sim$fireReturnIntervalsByPolygonNumeric)

        if (iter <= P(sim)$maxReburns[1]) {
          firesOK <- fires[!initialPixels %in% tooSmallDT$initialPixels, ]
          firesList <- append(firesList, list(firesOK))
          fireSizes <- append(fireSizes, list(fa[!tooSmall, c("size", "maxSize")]))
          maxOrder <- max(fires$order)

          polysNeedMoreFires <- tooSmallByPoly[, N := .N, by = polygonNumeric]
          polysNeedMoreFires <- polysNeedMoreFires[friByPolyDT, on = "polygonNumeric"]
          polysNeedMoreFires[is.na(N), N := 0]
          set(polysNeedMoreFires, NULL, "pixel", NULL)

          numFiresThisPeriod <- polysNeedMoreFires[, N[1], by = "polygonNumeric"]$V1
          fireSizesInPixels <- na.omit(polysNeedMoreFires)$maxSize
          spreadProbThisStep[firesOK$pixels] <- NA_real_
        } else {
          firesTooSmall <- fires[initialPixels %in% tooSmallDT$initialPixels, ]

          fa2 <- fa[tooSmall, c("size", "maxSize")] ## track the fires that did burn
          fa3 <- copy(fa2)                          ## track what's left to burn

          fa2[, maxSize := size] ## consider the area that did burn as having reached target

          fa3[, maxSize2 := maxSize - size]
          fa3[, size := 0]
          fa3[, maxSize := maxSize2]
          set(fa3, NULL, "maxSize2", NULL)

          firesList <- append(firesList, list(firesTooSmall))
          fireSizes <- append(fireSizes, list(fa2))
          maxOrder <- max(fires$order)

          polysNeedMoreFires <- tooSmallByPoly[, N := .N, by = polygonNumeric]
          polysNeedMoreFires <- polysNeedMoreFires[, maxSize := fa3$maxSize] ## update what's left to burn
          polysNeedMoreFires <- polysNeedMoreFires[friByPolyDT, on = "polygonNumeric"]
          polysNeedMoreFires[is.na(N), N := 0]
          set(polysNeedMoreFires, NULL, "pixel", NULL)

          numFiresThisPeriod <- polysNeedMoreFires[, N[1], by = "polygonNumeric"]$V1
          fireSizesInPixels <- na.omit(polysNeedMoreFires)$maxSize
          spreadProbThisStep[firesTooSmall$pixels] <- NA_real_
        }
      } else {
        firesList <- append(firesList, list(fires))
        fireSizes <- append(fireSizes, list(fa[, c("size", "maxSize")]))

        if (length(tooSmall) == 0) {
          assertthat::assert_that(fa1$proportionBurned %==% 1)
          numFiresThisPeriod <- rep(0L, length(numFiresThisPeriod))
        } else if (isTRUE(P(sim)$.unitTest)) {
          if (any(tail(fa1$proportionBurned, 10)  < P(sim)$minPropBurn)) {
            mess <- "In 'LandMine' module 'Burn()': proportion area burned is less than 'minPropBurn'!"
            if (verbose > 0)
              message(crayon::red(mess))
            warning(mess, call. = FALSE)
          }
        }
      }
    }

    iter <- iter + 1L
  }

  fires <- rbindlist(firesList)
  sim$fireSizes[[round(time(sim) - P(sim)$burnInitialTime + 1, 0)]] <- rbindlist(fireSizes)
  sim$rstCurrentBurn[] <- 0L
  sim$rstCurrentBurn[fires$pixels] <- 1L #as.numeric(factor(fires$initialPixels))

  if (is.null(sim$rstCurrentBurnCumulative)) {
    sim$rstCurrentBurnCumulative <- sim$rstCurrentBurn # keeps 1s
    sim$rstCurrentBurnCumulative[!is.na(sim$rstCurrentBurnCumulative[])
                                 & sim$rstCurrentBurnCumulative[] == 0] <- 0
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

  burnedDF <- data.frame(time = as.numeric(times(sim)$current),
                         nPixelsBurned = npix,
                         haBurned = npix * prod(res(sim$rstCurrentBurn)) / 100^2, ## area in ha
                         FRI = as.factor(fris))
  mod$areaBurnedOverTime <- rbind(mod$areaBurnedOverTime, burnedDF)
  mod$gg_areaBurnedOverTime <- landmine_plot_areaBurnedOverTime(mod$areaBurnedOverTime)

  if (time(sim) == end(sim)) {
    fgg_areaBurnedOverTime <- file.path(figurePath(sim), "LandMine_areaBurnedOverTime.png")
    ggsave(fgg_areaBurnedOverTime, mod$gg_areaBurnedOverTime)
  }

  return(invisible(sim))
})

### summary events
SummarizeFRIsingle <- function(sim) {
  studyArea <- P(sim)$.studyAreaName

  flammableMap <- sim[["rstFlammable"]]   ## RasterLayer
  lthfc <- sim[["fireReturnInterval"]]    ## RasterLayer
  pixelRes <- res(sim[["rasterToMatch"]]) ## c(250, 250)

  meanAnnualCumulBurnMap <- sim[["rstCurrentBurnCumulative"]] / (end(sim) - start(sim))

  ## sanity check
  compareRaster(flammableMap, lthfc, meanAnnualCumulBurnMap, res = TRUE, orig = TRUE)

  nonFlammable <- which(is.na(flammableMap[]) | flammableMap[] == 0)
  if (length(nonFlammable) > 0) {
    flammableMap[nonFlammable] <- NA
    lthfc[nonFlammable] <- NA
    meanAnnualCumulBurnMap[nonFlammable] <- NA
  }

  expFRIs <- raster::getValues(lthfc) |>
    unique() |>
    na.omit() |>
    sort()

  simFRIs <- vapply(expFRIs, function(fri) {
    pixIds <- which(raster::getValues(lthfc) == fri)
    1 / (sum(meanAnnualCumulBurnMap[pixIds]) / (length(pixIds)))
  }, numeric(1))

  sim$friSummary <- data.table(
    studyArea = studyArea,
    LTHFC = expFRIs,
    FRI = simFRIs,
    stringsAsFactors = FALSE
  )

  f <- file.path(outputPath(sim), paste0("LandMine_FRI_summary.csv"))
  fwrite(sim$friSummary, f) ## TODO: add this file to list of outputs

  ## LTHFC/FRI polygons
  ggFriPolys <- landmine_plot_LTHFC(lthfc, studyArea)

  if ("png" %in% P(sim)$.plots) {
    fggFriPolys <- file.path(figurePath(sim), "LandMine_LTHFC_map.png")
    png(fggFriPolys, height = 1000, width = 1000)
    print(ggFriPolys)
    dev.off()
  }

  ## expected vs simulated fire return intervals
  ggFriExpVsSim <- landmine_plot_FRI(sim$friSummary)

  if ("png" %in% P(sim)$.plots) {
    fggFriExpVsSim <- file.path(figurePath(sim), "LandMine_FRI_exp_vs_sim.png")
    ggsave(fggFriExpVsSim, ggFriExpVsSim, height = 10, width = 10) ## NOTE: keep square aspect ratio
  }

  if ("screen" %in% P(sim)$.plots) {
    mod$summaryDevice <- max(dev.list()) + 1
    quickPlot::dev(mod$summaryDevice, width = 12)
    gridExtra::grid.arrange(ggFriPolys, ggFriExpVsSim, nrow = 1, ncol = 2)
  }

  return(invisible(sim))
}

SummarizeFRImulti <- function(sim) {
  studyArea <- P(sim)$.studyAreaName

  allReps <- P(sim)$reps
  flammableMap <- NULL
  lthfc <- NULL
  pixelRes <- NULL

  burnMaps <- lapply(allReps, function(rep) {
    fsim <- findSimFile(outputPath(sim), rep)

    tmpSim <- loadSimList(fsim)

    if (rep == 1L) {
      ## all reps have same flammable + LTHFC maps
      flammableMap <<- tmpSim[["rstFlammable"]]   ## RasterLayer
      lthfc <<- tmpSim[["fireReturnInterval"]]    ## RasterLayer
      pixelRes <<- res(tmpSim[["rasterToMatch"]]) ## c(250, 250)

      ## sanity check
      compareRaster(tmpSim[["fireReturnInterval"]],
                    tmpSim[["rstFlammable"]],
                    tmpSim[["rstCurrentBurnCumulative"]],
                    res = TRUE, orig = TRUE)
    }

    ## mean annual cumulative burn map
    tmpSim[["rstCurrentBurnCumulative"]] / (end(tmpSim) - start(tmpSim))
  }) |> raster::stack() |>
    raster::calc(sum, na.rm = TRUE)

  meanAnnualCumulBurnMap <- burnMaps / length(allReps)

  nonFlammable <- which(is.na(flammableMap[]) | flammableMap[] == 0)
  if (length(nonFlammable) > 0) {
    flammableMap[nonFlammable] <- NA
    lthfc[nonFlammable] <- NA
    meanAnnualCumulBurnMap[nonFlammable] <- NA
  }

  expFRIs <- raster::getValues(lthfc) |>
    unique() |>
    na.omit() |>
    sort()

  simFRIs <- vapply(expFRIs, function(fri) {
    pixIds <- which(raster::getValues(lthfc) == fri)
    1 / (sum(meanAnnualCumulBurnMap[pixIds]) / (length(pixIds)))
  }, numeric(1))

  sim$friSummary <- data.table(
    studyArea = studyArea,
    LTHFC = expFRIs,
    FRI = simFRIs,
    stringsAsFactors = FALSE
  )

  f <- file.path(outputPath(sim), paste0("LandMine_FRI_summary_multi.csv"))
  fwrite(sim$friSummary, f) ## TODO: add this file to list of outputs

  ## LTHFC/FRI polygons
  ggFriPolys <- landmine_plot_LTHFC_raster(lthfc)

  if ("png" %in% P(sim)$.plots) {
    fggFriPolys <- file.path(figurePath(sim), "LandMine_LTHFC_map.png")
    png(fggFriPolys, height = 1000, width = 1000)
    print(ggFriPolys)
    dev.off()
  }

  ## expected vs simulated fire return intervals
  ggFriExpVsSim <- landmine_plot_compare_FRI(sim$friSummary) +
    geom_smooth(method = "lm")

  if ("png" %in% P(sim)$.plots) {
    fggFriExpVsSim <- file.path(figurePath(sim), "LandMine_FRI_exp_vs_sim.png")
    ggsave(fggFriExpVsSim, ggFriExpVsSim, height = 10, width = 10) ## NOTE: keep square aspect ratio
  }

  if ("screen" %in% P(sim)$.plots) {
    clearPlot()
    gridExtra::grid.arrange(ggFriPolys, fggFriExpVsSim, nrow = 1, ncol = 2)
  }

  return(invisible(sim))
}

## .inputObjects
.inputObjects <- function(sim) {
  ## DEBUGGING: random seed issues
  #fseed <- file.path(outputPath(sim), "seed.txt")
  #writeEventInfo(sim, fseed, append = TRUE)
  #writeRNGInfo(fseed, append = TRUE)
  ## END DEBUGGING

  #cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  if (getOption("LandR.verbose", TRUE) > 0)
    message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # Make random forest cover map
  mod$numDefaultPixelGroups <- 20L
  mod$numDefaultPolygons <- 4L
  numDefaultSpeciesCodes <- 2L

  if (!suppliedElsewhere("studyArea", sim)) {
    if (getOption("LandR.verbose", TRUE) > 0)
      message("'studyArea' was not provided by user. Using a polygon in southwestern Alberta, Canada,")

    sim$studyArea <- randomStudyArea(seed = 1234, size = 1e9)
  }

  if (!is(sim$studyArea, "Spatial"))
    sim$studyArea <- as(sim$studyArea, "Spatial")

  if (!suppliedElsewhere("studyAreaReporting", sim)) {
    if (getOption("LandR.verbose", TRUE) > 0)
      message("'studyAreaReporting' was not provided by user. Using the same as 'studyArea'.")
    sim$studyAreaReporting <- sim$studyArea
  }

  if (is.null(sim$rasterToMatch)) {
    if (!suppliedElsewhere("rasterToMatch", sim)) {
      sim$rasterToMatch <- raster(sim$studyArea, res = 100)
      sim$rasterToMatch <- fasterize::fasterize(sf::st_as_sf(sim$studyArea), sim$rasterToMatch)
    }
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
    sim$fireReturnInterval[] <- as.integer(as.character(vals))
  }

  ## 2023-09: ensure fireReturnInterval map has non-flammable pixels removed
  nonFlammable <- which(is.na(sim[["rstFlammable"]][]) | sim[["rstFlammable"]][] == 0)
  if (length(nonFlammable) > 0) {
    sim$fireReturnInterval[nonFlammable] <- NA
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
    sim$sppEquiv <- LandR::sppEquivalencies_CA
    sppNames <- LandR::equivalentName(sim$species$species, sim$sppEquiv, column = Par$sppEquivCol)
    sim$sppEquiv <- sim$sppEquiv[get(Par$sppEquivCol) %in% sppNames]
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
  if (length(immature))
    ROS[immature] <- sppEquiv["immature"]$ros[match(vegType[immature], sppEquiv["immature"]$pixelValue)]
  if (length(young))
    ROS[young] <- sppEquiv["young"]$ros[match(vegType[young], sppEquiv["young"]$pixelValue)]

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
  ## 2023-02: discontinuous fuels (e.g., shield) requires increasing spread --
  ##          use same value as young deciduous (6L), per Dave's text messages
  ROSother <- switch(P(sim)$ROStype,
                     burny = sim$ROSTable[leading == "decid" & age == "immature_young", ros],
                     sim$ROSTable[leading == "spruce" & age == "mature", ros])

  assertthat::assert_that(
    isTRUE(inRange(P(sim)$ROSother, min(sim$ROSTable$ros), max(sim$ROSTable$ros))),
    isTRUE(inRange(P(sim)$ROSother, 0.95*ROSother, 1.05*ROSother)) ## TODO: tweak this to allow greater range
  )
  ROS[sim$rstFlammable[] == 1L & is.na(ROS)] <- as.integer(P(sim)$ROSother)
  ROS[sim$rstFlammable[] == 0L | is.na(sim$rstFlammable[])] <- NA ## non-flammable pixels

  return(ROS)
})

## older version of SpaDES.core used here doesn't have this function
if (packageVersion("SpaDES.core") < "2.0.2.9001") {
  figurePath <- function(sim) {
    file.path(outputPath(sim), "figures", current(sim)[["moduleName"]]) |>
      checkPath(create = TRUE)
  }
}
