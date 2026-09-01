defineModule(sim, list(
  name = "LandMine",
  description = "Reimplementation of Andison (1999) LandMine fire model",
  keywords = c("Fire", "Landscape", "Percolation", "Pixel-based"),
  authors = c(
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb", "cre"))
  ),
  childModules = character(0),
  version = list(LandMine = numeric_version("1.0.9")),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "LandMine.Rmd"),
  reqdPkgs = list(
    "assertthat", "cli", "data.table", "fpCompare", "ggplot2",
    "RColorBrewer", "stats", "terra", "tidyterra", "VGAM",
    "PredictiveEcology/LandR@development (>= 1.1.0.9003)",
    "PredictiveEcology/LandWebUtils@development (>= 1.0.3.9032)",
    "PredictiveEcology/pemisc@development",
    "PredictiveEcology/SpaDES.tools@development (>= 2.1.2.9000)"
  ),
  parameters = rbind(
    defineParameter("biggestPossibleFireSizeHa", "numeric", 1e6, 1e4, 2e6,
                    "An upper limit, in hectares, of the truncated Pareto distribution of fire sizes"),
    defineParameter("burnInitialTime", "numeric", start(sim, "year") + 1, NA, NA,
                    "This describes the simulation time at which the first burn event should occur"),
    defineParameter("fireTimestep", "integer", 1L, NA, NA,
                    "This describes the simulation time interval between burn events"),
    defineParameter("maxReburns", "integer", c(1L, 20L), 1L, 500L,
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
                          "Row 1 is the 120 m fit, and is the default at EVERY resolution:",
                          "it scores as well as the 240 m fit on the 240 m landscape (0.70 vs",
                          "0.81, 0.33 pooled sd) with a fifth of the variance, while the 240 m",
                          "fit is much worse at 120 m (1.78 vs 0.61, 4.91 sd). The 120 m fit is",
                          "also far better determined: repeat runs land within 4.6% of the",
                          "search interval, against 24.8% at 240 m.",
                          "`2L` is the 240 m fit, retained for reference. Rows fitted before",
                          "2026-08 were removed: they used a different objective (a size penalty",
                          "that never fired, a constant perimeter:area target, and fire sizes in",
                          "pixels rather than hectares) and are not comparable. They remain in",
                          "git history.")),
    defineParameter("reps", "integer", NA_integer_, 1L, NA_integer_,
                    paste("number of replicates/runs per study area when running in 'multi' mode.")),
    defineParameter("ROSother", "integer", 30L, NA, NA,
                    paste0("default ROS value for non-forest vegetation classes.",
                           "this is needed when passing a modified `ROSTable`, e.g. using log-transformed values.")),
    defineParameter("ROStype", "character", "default", NA, NA,
                    "One of 'default' or 'burny'."),
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
    expectsInput("fireReturnInterval", "SpatRaster",
                 desc = paste("A raster layer that is a factor raster, with at least 1 column called",
                              "`fireReturnInterval`, representing the fire return interval in years."),
                 sourceURL = NA),
    expectsInput("pixelGroupMap", "SpatRaster",
                 desc = "Pixels with identical values share identical stand features",
                 sourceURL = NA),
    expectsInput("rasterToMatch", "SpatRaster",
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
    expectsInput("flammableMap", "SpatRaster",
                 desc = paste("A raster layer, with 0, 1 and NA, where 1 indicates areas",
                              "that are flammable, 0 not flammable (e.g., lakes)",
                              "and NA not applicable (e.g., masked)")),
    expectsInput("rstTimeSinceFire", "SpatRaster",
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
    expectsInput("studyArea", "SpatVector",
                 desc = paste("multipolygon, typically buffered around an area of interest",
                              "(i.e., `studyAreaReporting`) to use for simulation.",
                              "Defaults to an area in Southwestern Alberta, Canada."),
                 sourceURL = NA),
    expectsInput("studyAreaReporting", "sf",
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
    createsOutput("fireReturnInterval", "SpatRaster", paste(
      "A `Raster` map showing the fire return interval. This is created from the `rstCurrentBurn`.")
    ),
    createsOutput("fireReturnIntervalsByPolygonNumeric", "numeric", paste(
      "A vector of the fire return intervals, ordered by the numeric representation of polygon ID")
    ),
    createsOutput("fireTimestep", "numeric", paste(
      "The number of time units between successive fire events in a fire module.")
    ),
    createsOutput("friSummary", "data.table", "summary fire return interval table"),
    createsOutput("friDiagnostics", "data.table", paste(
      "per-FRI-zone attainment diagnostics: achieved/target ratio over the BURNABLE area,",
      "the same ratio computed the way `friSummary` does it, the share of each zone that",
      "lies inside the study area, and the structural measurements that tell the failure",
      "modes apart (flammable fraction, patch structure, cohort coverage)."
    )),
    createsOutput("kBest", "numeric", paste(
      "A numeric scalar that is the optimal value of `K` in the",
      "Truncated Pareto distribution (`rtruncpareto`)")
    ),
    createsOutput("numFiresPerYear", "numeric", paste(
      "The average number of fires per year, by fire return interval level on `rstCurrentBurn`.")
    ),
    createsOutput("rstCurrentBurn", "SpatRaster", paste(
      "A raster layer, produced at each timestep, where each",
      "pixel is either 1 or 0 indicating burned or not burned.")
    ),
    createsOutput("burnMap", "SpatRaster", "Cumulative number of times a pixel has burned"),
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
      sim <- scheduleEvent(sim, end(sim), "LandMine", "plot")
      sim <- scheduleEvent(sim, end(sim), "LandMine", "summarySingle")
    } else if (P(sim)$mode == "multi") {
      sim <- scheduleEvent(sim, start(sim), "LandMine", "summaryMulti")
    }
  } else if (eventType == "plot") {
    sim <- plotFn(sim)
    sim <- scheduleEvent(sim, P(sim)$.plotInterval, "LandMine", "plot")
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
    # meanFS <- LandWebUtils::meanTruncPareto(k = params[1], lower = 1, upper = upper1, alpha = 1)
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
    cachePath = cachePath(sim),
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

  terra::compareGeom(sim$rasterToMatch, sim$fireReturnInterval, sim$flammableMap, sim$rstTimeSinceFire)

  ## from DEoptim fitting; see `LandWebUtils::landmine_optim_calibrate()` and LandMine.Rmd.
  ## The CSV schema and the `10^` parameter convention are defined once, in LandWebUtils.
  optimPars <- LandWebUtils::landmine_optim_unpack(
    LandWebUtils::landmine_optim_params_read(
      file.path(dataPath(sim), "LandMine_DEoptim_params.csv"),
      rowID = P(sim)$optimParsRowID
    )
  )

  mod$spawnNewActive <- optimPars$spawnNewActive
  mod$sizeCutoffs <- optimPars$sizeCutoffs

  ## 2024-10: use `sim$flammableMap` instead of `sim$fireReturnInterval` raster
  ##          to ensure coverage across entire studyArea (e.g., boundaries along grasslands)
  mod$spreadProb <- terra::rast(sim$flammableMap)
  mod$spreadProb[sim$flammableMap[] == 1] <- optimPars$spreadProb ## assign spreadProb to flammable pixels
  mod$spreadProb[is.na(sim$flammableMap[]) | sim$flammableMap[] == 0] <- switch(
    P(sim)$ROStype,
    burny = optimPars$spreadProb, ## non-flammable pixels *can* spread fire (but won't count as burned pixels for fire stats)
    NA_real_ ## non-flammable pixels don't have a spreadProb value (i.e., can't spread)
  )
  mod$spreadProb <- terra::mask(mod$spreadProb, sim$studyArea) ## mask only; don't crop

  sim$fireSizes <- list()

  if (!is.integer(sim$fireReturnInterval[])) {
    sim$fireReturnInterval[] <- as.integer(sim$fireReturnInterval[])
  }

  if (verbose > 0) {
    message("Initializing fire maps")
  }
  sim$fireTimestep <- P(sim)$fireTimestep
  sim$fireInitialTime <- P(sim)$burnInitialTime

  if (verbose > 0) {
    message("Determine mean fire size and per-zone ignition budget...")
  }
  ## Promoted to LandWebUtils: masks zero-FRI and non-flammable pixels out of the raster,
  ## tabulates pixels per FRI zone, and converts to expected fires/yr
  ## (`(area / meanFireSize) / FRI`). Unit-tested there for the ORDERING CONTRACT that the
  ## reburn loop below consumes POSITIONALLY via `numFiresThisPeriod[.GRP]`: `terra::freq()`
  ## sorts ascending, the NA row is appended last, and that entry is NA-VALUED so `na.omit()`
  ## drops it (an NA *name* alone would not). Break any one and zones silently receive each
  ## other's fire counts. Masking here is also why non-flammable pixels -- and now pixels outside
  ## the study area -- can never be ignition locations: the start-cell pool is built from this
  ## same raster.
  ## `studyArea` matters as much as `flammableMap` here, and for the same reason. `ROSmap` and
  ## `mod$spreadProb` are both masked to it, so a fire ignited outside burns its own start cell
  ## and spreads no further -- yet `fireReturnInterval` is not masked and overhangs the polygon,
  ## so those pixels were inflating each zone's area and therefore its expected fires per year,
  ## AND entering the start-cell pool. On WesternAlbertaUpland that allocated 2,918 fires/yr
  ## against a correct 2,015: 31% of every year's ignitions aimed at unburnable ground, with one
  ## zone (99.95% outside) over-allocated 2,008-fold.
  ignitionBudget <- LandWebUtils::landmine_ignition_budget(
    fireReturnInterval = sim$fireReturnInterval,
    flammableMap = sim[["flammableMap"]],
    kBest = sim$kBest,
    biggestPossibleFireSizeHa = P(sim)$biggestPossibleFireSizeHa,
    studyArea = sim$studyArea
  )
  sim$fireReturnInterval <- ignitionBudget$fireReturnInterval
  sim$fireReturnIntervalsByPolygonNumeric <- ignitionBudget$fireReturnIntervalsByPolygonNumeric
  sim$numFiresPerYear <- ignitionBudget$numFiresPerYear

  sim$rstCurrentBurn <- terra::rast(sim$fireReturnInterval) ## creates no-value raster
  sim$rstCurrentBurn[] <- 0L
  if (verbose > 0) {
    message("6: ", Sys.time())
  }

  mod$areaBurnedOverTime <- data.frame(
    time = numeric(0),
    nPixelsBurned = numeric(0),
    haBurned = numeric(0),
    FRI = numeric(0)
  )

  ## knownSpecies needs to use 'LandWeb' column, not 'LandR'!
  mod$knownSpecies <- LandWebUtils::landmine_known_species()
  sim$sppEquiv[, LandMine := mod$knownSpecies[LandWeb]]

  return(invisible(sim))
}

### plot events
plotFn <- compiler::cmpfun(function(sim) {
  if (time(sim) == P(sim)$.plotInitialTime) {
    gg_fri <- ggplot() +
      tidyterra::geom_spatraster(data = terra::as.factor(sim$fireReturnInterval)) +
      tidyterra::scale_fill_terrain_d() +
      tidyterra::geom_spatvector(data = sim$studyAreaReporting, fill = NA, linewidth = 1.5) +
      ggtitle("Fire Return Interval (i.e., LTHFC)") +
      theme_minimal()

    gg_flm <- ggplot() +
      tidyterra::geom_spatraster(data = terra::as.factor(sim$flammableMap)) +
      tidyterra::scale_fill_whitebox_d(palette = "bl_yl_rd") +
      tidyterra::geom_spatvector(data = sim$studyAreaReporting, fill = NA, linewidth = 1.5) +
      ggtitle("Landscape flammability") +
      theme_minimal()

    if ("png" %in% P(sim)$.plots) {
      f_gg_fri <- file.path(figurePath(sim), "LandMine_fireReturnInterval.png")
      ggsave(f_gg_fri, gg_fri)
      sim <- registerOutputs(f_gg_fri)

      f_gg_flm <- file.path(figurePath(sim), "LandMine_flammableMap.png")
      ggsave(f_gg_flm, gg_flm)
      sim <- registerOutputs(f_gg_flm)
    }

    if ("screen" %in% P(sim)$.plots) {
      print(cowplot::plot_grid(gg_fri, gg_flm))
    }
  } else {
    gg_abot <- mod$gg_areaBurnedOverTime +
      ggtitle("Current area burned (ha)")

    rcbc <- sim$burnMap
    rcbc[!is.na(sim$rstCurrentBurn) & sim$rstCurrentBurn == 0 &
           sim$burnMap == 0] <- 0L

    gg_cbc <- ggplot() +
      tidyterra::geom_spatraster(data = rcbc) +
      tidyterra::scale_fill_princess_c(palette = "maori") +
      tidyterra::geom_spatvector(data = sim$studyAreaReporting, fill = NA, linewidth = 1.5) +
      ggtitle(sprintf("Cumulative Fire Map (t = %d)", time(sim))) +
      theme_minimal()

    if ("png" %in% P(sim)$.plots) {
      f_gg_cbc <- file.path(figurePath(sim), sprintf("LandMine_cumulative_burn_map_%04d.png", time(sim)))
      ggsave(f_gg_cbc, gg_cbc)
      sim <- registerOutputs(f_gg_cbc)
    }

    if ("screen" %in% P(sim)$.plots) {
      print(cowplot::plot_grid(gg_abot, gg_cbc))
    }
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
  ## Promoted to LandWebUtils (unit-tested: names survive the draw -- `rnbinom()` drops them,
  ## and the counts are indexed by zone below -- and the overdispersion still gives a
  ## low-rate zone its long runs of zero-fire years).
  numFiresThisPeriod <- LandWebUtils::landmine_draw_num_fires(
    sim$numFiresPerYear,
    fireTimestep = P(sim)$fireTimestep
  )
  thisYrStartCellsDT <- data.table(
    pixel = seq(terra::ncell(sim$fireReturnInterval)),
    polygonNumeric = terra::values(sim$fireReturnInterval, mat = FALSE),
    key = "polygonNumeric"
  )

  ## August 2022: reburn fires that did not meet their target size

  ## Rate of Spread
  vegTypeMap <- LandR::vegTypeMapGenerator(
    sim$cohortData,
    pixelGroupMap = sim$pixelGroupMap,
    vegLeadingProportion = P(sim)$vegLeadingProportion,
    mixedType = P(sim)$mixedType,
    sppEquiv = sim$sppEquiv,
    sppEquivCol = P(sim)$sppEquivCol,
    colors = sim$sppColorVect,
    doAssertion = P(sim)$.unitTest
  )
  ROSmap <- terra::rast(sim$pixelGroupMap) ## empty raster as template
  ROSmap[] <- LandWebUtils::landmine_fire_ros(
    vegTypeMap = vegTypeMap,
    rstTimeSinceFire = sim$rstTimeSinceFire,
    flammableMap = sim$flammableMap,
    ROSTable = sim$ROSTable,
    sppEquiv = sim$sppEquiv,
    sppEquivCol = P(sim)$sppEquivCol,
    ROSother = P(sim)$ROSother,
    knownSpecies = mod$knownSpecies,
    ROStype = P(sim)$ROStype
  )
  ROSmap <- terra::mask(ROSmap, sim$studyArea)
  ## PERFORMANCE: for the same reason as `spreadProbThisStep` below, a raster
  ## `spreadProbRel` is re-materialised by `spread2()` on *every* spread step.
  ## Numeric `spreadProbRel` requires SpaDES.tools >= 2.1.2.9000 (SpaDES.tools#106).
  ROSvals <- terra::values(ROSmap, mat = FALSE)

  ## PERFORMANCE: `spread2()` re-materialises a `SpatRaster` `spreadProb` on *every*
  ## spread step (`as.vector(spreadProb)` in SpaDES.tools' spread2.R). Because
  ## `landmine_burn1()` calls `spread2(iterations = 1L)` in a while loop, that is one
  ## O(ncell) read per step -- ~36 MB per step on a 4.5 Mpix study area. Passing a
  ## numeric vector instead is ~1.5x faster overall and bit-identical (same values, so
  ## the same cells burn). It also makes the `[<-` NA-marking below an O(k) vector
  ## write rather than an O(ncell) terra raster write.
  spreadProbThisStep <- terra::values(mod$spreadProb, mat = FALSE)

  ## If fire sizes are in hectares, must adjust based on resolution of maps
  ##  NOTE: round causes fires < 0.5 pixels to NOT EXIST ...
  ##        e.g., 3.25 ha fires are "not detectable" if resolution is 6.25 ha
  fireSizesThisPeriod <- VGAM::rtruncpareto(
    n = sum(numFiresThisPeriod),
    lower = 1,
    upper = P(sim)$biggestPossibleFireSizeHa,
    shape = sim$kBest
  )

  ## Because annual number of fires includes fires <6.25 ha, sometimes this will round down to 0 pixels.
  ##   This calculation makes that probabilistic.
  ## Promoted to LandWebUtils (unit-tested: the rounding preserves E[pixels], which is the
  ## entire reason it is probabilistic, and exact pixel multiples never round up).
  fireSizesInPixels <- LandWebUtils::landmine_sizes_to_pixels(
    fireSizesThisPeriod,
    pixelAreaHa = prod(res(sim$flammableMap)) / 1e4
  )

  firesList <- fireSizes <- list()
  maxOrder <- 0L
  iter <- 1L

  ## 2024-10: normally, non-flammable pixels are NA in ROSvals and spreadProbThisStep;
  ##          except in 'burny' scenarios, where fires *can* spread through those pixels,
  ##          but aren't counted towards burn stats/summaries.
  ## PERFORMANCE: hoisted out of the reburn loop below -- `sim$flammableMap` is not
  ## modified within this event, so this is loop-invariant. It was previously recomputed
  ## every reburn iteration (two O(ncell) raster reads + a `which()` per round, and
  ## burny runs average ~6 rounds/yr).
  omitPixels <- switch(
    P(sim)$ROStype,
    burny = which(is.na(sim$flammableMap[]) | (sim$flammableMap[] == 0)),
    NULL
  )

  ## PERFORMANCE: the eligible start-cell pool is loop-invariant -- the `:=` filter is
  ## idempotent after the first pass and `thisYrStartCellsDT` is not otherwise modified --
  ## so build it ONCE rather than re-filtering and re-`na.omit()`ing every reburn round.
  ## Only the per-polygon resample has to be redone. On a 4.5 Mpix study area this is
  ## ~306 ms -> ~24 ms per round, and WesternAlbertaUpland averages ~8 rounds/yr.
  ## Verified to give the same start cells, in the same order, for the same seed: order
  ## matters because `fireSizesInPixels` is matched positionally to `thisYrStartCells`
  ## and `numFiresThisPeriod[.GRP]` depends on group order.
  ## `polygonNumeric` is read from `sim$fireReturnInterval`, which Init() has ALREADY set to NA
  ## on both zero-FRI and non-flammable pixels -- so `na.omit()` alone drops every ineligible
  ## start cell, and non-flammable pixels are never ignition locations. The `> 0` filter is a
  ## belt-and-braces guard for a caller-supplied raster that skipped Init's zero handling.
  ##
  ## This previously also tested `polygonNumeric %in% NA_ids`, where `NA_ids` came from
  ## `attr(na.omit(numFiresPerYear), "na.action")` -- i.e. POSITIONS compared against FRI
  ## VALUES. It was harmless only by accident (no study area has an FRI equal to one of those
  ## positions). Had it ever matched, that zone would have been dropped from `startCellPool`
  ## while remaining in `numFiresPerYear`, and the positional `.GRP` indexing below would then
  ## have silently handed every subsequent zone another zone's fire count.
  startCellPool <- na.omit(thisYrStartCellsDT)[polygonNumeric > 0]
  data.table::setkeyv(startCellPool, "polygonNumeric")

  ## The `numFiresThisPeriod[.GRP]` indexing below is POSITIONAL: it assumes the groups of
  ## `startCellPool` (keyed, so ascending `polygonNumeric`) line up 1:1 and in order with
  ## `sim$numFiresPerYear` (named by FRI, ascending from `terra::freq()`). Nothing enforces
  ## that, and a mismatch would not error -- it would quietly simulate the wrong fire regime.
  ## Assert it, so a future change to the masking, the FRI values, or the freq() ordering
  ## fails loudly here instead of producing plausible-looking wrong output.
  assertthat::assert_that(
    identical(
      as.character(startCellPool[, unique(polygonNumeric)]),
      names(sim$numFiresPerYear)
    ),
    msg = paste(
      "LandMine: start-cell pool zones do not match `numFiresPerYear` -- the positional",
      "`.GRP` indexing of `numFiresThisPeriod` would assign fire counts to the wrong FRI",
      "zones. Pool:", paste(startCellPool[, unique(polygonNumeric)], collapse = ","),
      "| numFiresPerYear:", paste(names(sim$numFiresPerYear), collapse = ",")
    )
  )

  ## 2023-09: after maxReburns, if not reaching fire size, take the last burn,
  ## and start new fire(s) to burn the remaining area until the target is achieved.
  ## Should be OK b/c LandMine replicates FRIs (i.e., area burned each year), not number of fires
  polysNeedMoreFires <- NULL ## set inside the loop; NULL when the loop body never reburns
  while (sum(numFiresThisPeriod) > 0 && (iter <= sum(P(sim)$maxReburns))) {
    thisYrStartCells <- startCellPool[
      , SpaDES.tools:::resample(pixel, numFiresThisPeriod[.GRP]), by = polygonNumeric
    ]$V1

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

      fires <- LandWebUtils::landmine_burn1(
        sim$fireReturnInterval,
        startCells = thisYrStartCells,
        fireSizes = fireSizesInPixels,
        spreadProbRel = ROSvals,
        sizeCutoffs = mod$sizeCutoffs,
        maxRetriesPerID = P(sim)$maxRetriesPerID,
        spawnNewActive = mod$spawnNewActive,
        spreadProb = spreadProbThisStep,
        omitPixels = omitPixels
      )

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

        if (iter <= P(sim)$maxReburns[1]) {
          firesOK <- fires[!initialPixels %in% tooSmallDT$initialPixels, ]
          firesList <- append(firesList, list(firesOK))
          fireSizes <- append(fireSizes, list(fa[!tooSmall, c("size", "maxSize")]))
          maxOrder <- max(fires$order)

          ## Promoted to LandWebUtils. Phase 1: each too-small fire keeps its FULL original
          ## target and is re-ignited from a fresh start cell.
          reburn <- LandWebUtils::landmine_reburn_budget(
            tooSmallByPoly, sim$fireReturnIntervalsByPolygonNumeric
          )
          polysNeedMoreFires <- reburn$polysNeedMoreFires
          numFiresThisPeriod <- reburn$numFiresThisPeriod
          fireSizesInPixels <- reburn$fireSizesInPixels
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

          ## Promoted to LandWebUtils. Phase 2: NEW fires sized to the REMAINING shortfall
          ## (`fa3$maxSize`), assigned positionally onto the too-small rows -- both descend
          ## from the same `fa[tooSmall]` subset. The promoted function is equivalence-tested
          ## against this exact logic over randomised multi-zone inputs.
          reburn <- LandWebUtils::landmine_reburn_budget(
            tooSmallByPoly, sim$fireReturnIntervalsByPolygonNumeric,
            remainingSize = fa3$maxSize
          )
          polysNeedMoreFires <- reburn$polysNeedMoreFires
          numFiresThisPeriod <- reburn$numFiresThisPeriod
          fireSizesInPixels <- reburn$fireSizesInPixels
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
              message(cli::col_red(mess))
            warning(mess, call. = FALSE)
          }
        }
      }
    }

    iter <- iter + 1L
  }

  ## The loop above exits EITHER because every fire reached its target (converged) OR because
  ## it ran out of reburn attempts. Those two are indistinguishable in the log otherwise, and
  ## the second one means this year's area-burned target was abandoned unmet. Record WHERE:
  ## `polygonNumeric` is the fire-return-interval zone, so this says which FRI zones the
  ## reburn ceiling is failing to satisfy -- the thing a `maxReburns[2]` change has to fix.
  ## One machine-parseable line per zone; only emitted when the ceiling actually binds.
  ## The tabulation is promoted to LandWebUtils (unit-tested: zero-outstanding zones and the
  ## NA-FRI row are excluded, and NULL/empty input gives zero rows rather than an error --
  ## a silently wrong diagnostic would mislead the next `maxReburns` decision).
  if (sum(numFiresThisPeriod) > 0) {
    unmet <- LandWebUtils::landmine_reburn_ceiling(polysNeedMoreFires, year = time(sim))
    for (i in seq_len(NROW(unmet))) {
      message(sprintf(
        "reburn-ceiling year=%g FRI=%g nFires=%d pixelsShort=%d",
        unmet$year[i], unmet$FRI[i], unmet$nFires[i], unmet$pixelsShort[i]
      ))
    }
  }

  firesDT <- rbindlist(firesList)
  fireSizesDT <- rbindlist(fireSizes)

  ## TODO: how to best deal with no fires and their impacts on FRI & area burned calculations?
  if (nrow(firesDT) == 0) {
    message(cli::col_yellow("no fires this period!"))
    firesDT <- data.table(
      initialPixels = integer(0),
      pixels = integer(0),
      state = character(0),
      order = integer(0)
    )
  }

  if (nrow(fireSizesDT) == 0) {
    fireSizesDT <- data.table(size = 0L, maxSize = 0L)
  }

  sim$fireSizes[[round(time(sim) - P(sim)$burnInitialTime + 1, 0)]] <- fireSizesDT

  sim$rstCurrentBurn[] <- 0L
  if (nrow(firesDT) > 0) {
    sim$rstCurrentBurn[firesDT$pixels] <- 1L # as.numeric(factor(firesDT$initialPixels))
  }

  if (is.null(sim$burnMap)) {
    sim$burnMap <- terra::deepcopy(sim$rstCurrentBurn) ## keeps 1s
    sim$burnMap[!is.na(sim$burnMap[])
                                 & sim$burnMap[] == 0] <- 0
  } else {
    sim$burnMap <- sim$rstCurrentBurn + sim$burnMap
  }

  currBurn <- terra::mask(sim$rstCurrentBurn, sim$studyAreaReporting)
  ## NOTE: this previously counted burned pixels with `table(currBurn[ids])[2]`, which returns
  ## NA when a zone burns COMPLETELY (only one level in the table) -- and the following
  ## `npix[is.na(npix)] <- 0` then recorded that total burn as ZERO ha. Promoted to
  ## LandWebUtils, where it counts with `sum(vals == 1)` and is unit-tested for the
  ## fully-burned, partly-burned and unburned cases.
  burnedDF <- LandWebUtils::landmine_area_burned_by_zone(
    currentBurn = currBurn,
    fireReturnInterval = sim$fireReturnInterval,
    time = as.numeric(times(sim)$current),
    pixelAreaHa = prod(res(sim$rstCurrentBurn)) / 100^2
  )
  mod$areaBurnedOverTime <- rbind(mod$areaBurnedOverTime, burnedDF)
  mod$gg_areaBurnedOverTime <- landmine_plot_areaBurnedOverTime(mod$areaBurnedOverTime)

  if (time(sim) == end(sim)) {
    fgg_areaBurnedOverTime <- file.path(figurePath(sim), "LandMine_areaBurnedOverTime.png")
    ggsave(fgg_areaBurnedOverTime, mod$gg_areaBurnedOverTime)
    sim <- registerOutputs(fgg_areaBurnedOverTime)
  }

  return(invisible(sim))
})

### summary events
SummarizeFRIsingle <- function(sim) {
  studyAreaName <- P(sim)$.studyAreaName

  flammableMap <- sim[["flammableMap"]]   ## SpatRaster
  lthfc <- sim[["fireReturnInterval"]]    ## SpatRaster
  pixelRes <- res(sim[["rasterToMatch"]]) ## c(240, 240)

  meanAnnualCumulBurnMap <- sim[["burnMap"]] / (end(sim) - start(sim))

  ## sanity check
  compareGeom(flammableMap, lthfc, meanAnnualCumulBurnMap, res = TRUE)

  ## Promoted to LandWebUtils: this block was duplicated VERBATIM in the single- and
  ## multi-mode summaries, so the two could silently diverge. It is now one unit-tested
  ## function (covering the non-flammable masking, the never-burned -> Inf case, and the
  ## implicit contract that the NA masks of `lthfc` and the burn map agree).
  ##
  ## `studyArea` is what keeps unburnable pixels OUT of the denominator: `ROSmap` and
  ## `mod$spreadProb` are both masked to it, so a pixel outside cannot be ignited or spread
  ## into, while `fireReturnInterval` is not masked and overhangs the polygon. Without this,
  ## every such pixel inflates its zone's achieved interval -- on WesternAlbertaUpland that was
  ## 29.4% of flammable zone pixels, and made two zones look as though they under-burned by
  ## 3.5x and 12.2x when in fact every zone was within 0.90-1.07 of its target.
  sim$friSummary <- LandWebUtils::landmine_fri_summary(
    lthfc = lthfc,
    flammableMap = flammableMap,
    meanAnnualCumulBurnMap = meanAnnualCumulBurnMap,
    studyAreaName = studyAreaName,
    studyArea = sim$studyArea
  )

  ## per-zone attainment diagnostics: a CSV a developer can diff between runs, two figures,
  ## and a one-line verdict, so nobody has to read this log to learn a zone missed its target.
  sim$friDiagnostics <- LandWebUtils::landmine_fri_metrics(
    lthfc = lthfc,
    flammableMap = flammableMap,
    meanAnnualCumulBurnMap = meanAnnualCumulBurnMap,
    studyAreaName = studyAreaName,
    pixelAreaHa = prod(pixelRes) / 1e4,
    studyArea = sim$studyArea,
    nYears = end(sim) - start(sim)
  )

  fDiag <- file.path(outputPath(sim), "LandMine_FRI_diagnostics.csv")
  fwrite(sim$friDiagnostics, fDiag)
  sim <- registerOutputs(fDiag)

  message(LandWebUtils::landmine_fri_verdict(sim$friDiagnostics))

  if ("png" %in% P(sim)$.plots) {
    fggFriZones <- file.path(figurePath(sim), "LandMine_FRI_zones_diagnostic.png")
    ggsave(fggFriZones, LandWebUtils::landmine_plot_fri_zones(
      lthfc, sim$friDiagnostics, studyAreaName
    ), height = 7, width = 12)
    sim <- registerOutputs(fggFriZones)

    fggFriDrivers <- file.path(figurePath(sim), "LandMine_FRI_drivers.png")
    ggsave(fggFriDrivers, LandWebUtils::landmine_plot_fri_drivers(
      sim$friDiagnostics, studyAreaName
    ), height = 5, width = 12)
    sim <- registerOutputs(fggFriDrivers)
  }

  f <- file.path(outputPath(sim), paste0("LandMine_FRI_summary.csv"))
  fwrite(sim$friSummary, f)
  sim <- registerOutputs(f)

  ## LTHFC/FRI polygons
  ggFriPolys <- landmine_plot_LTHFC(lthfc, studyAreaName) ## rasterVis::levelplot

  if ("png" %in% P(sim)$.plots) {
    fggFriPolys <- file.path(figurePath(sim), "LandMine_LTHFC_map.png")

    ## NOTE: this is a rasterVis::levelplot (not ggplot)
    png(fggFriPolys, height = 1000, width = 1000)
    print(ggFriPolys)
    dev.off()

    sim <- registerOutputs(fggFriPolys)
  }

  ## expected vs simulated fire return intervals
  ggFriExpVsSim <- landmine_plot_FRI(sim$friSummary)

  if ("png" %in% P(sim)$.plots) {
    fggFriExpVsSim <- file.path(figurePath(sim), "LandMine_FRI_exp_vs_sim.png")
    ggsave(fggFriExpVsSim, ggFriExpVsSim, height = 10, width = 10) ## NOTE: keep square aspect ratio
    sim <- registerOutputs(fggFriExpVsSim)
  }

  return(invisible(sim))
}

SummarizeFRImulti <- function(sim) {
  studyAreaName <- P(sim)$.studyAreaName

  allReps <- P(sim)$reps
  flammableMap <- NULL
  lthfc <- NULL
  pixelRes <- NULL

  burnMaps <- lapply(allReps, function(rep) {
    message(paste("loading simulation data for rep", rep, "..."))
    fsim <- findSimFile(outputPath(sim), rep)

    tmpSim <- suppressMessages(loadSimList(fsim))

    if (rep == 1L) {
      ## all reps have same flammable + LTHFC maps
      flammableMap <<- tmpSim[["flammableMap"]]   ## SpatRaster
      lthfc <<- tmpSim[["fireReturnInterval"]]    ## SpatRaster
      pixelRes <<- res(tmpSim[["rasterToMatch"]]) ## c(240, 240)

      ## sanity check
      compareGeom(tmpSim[["fireReturnInterval"]],
                  tmpSim[["flammableMap"]],
                  tmpSim[["burnMap"]],
                  res = TRUE)
    }

    ## mean annual cumulative burn map
    tmpSim[["burnMap"]] / (end(tmpSim) - start(tmpSim))
  }) |>
    terra::c() |>
    terra::app(sum, na.rm = TRUE) ## TODO: confirm c() and app() replaces stack() and calc()

  meanAnnualCumulBurnMap <- burnMaps / length(allReps)

  ## Promoted to LandWebUtils: this block was duplicated VERBATIM in the single- and
  ## multi-mode summaries, so the two could silently diverge. It is now one unit-tested
  ## function (covering the non-flammable masking, the never-burned -> Inf case, and the
  ## implicit contract that the NA masks of `lthfc` and the burn map agree).
  ##
  ## `studyArea` is what keeps unburnable pixels OUT of the denominator: `ROSmap` and
  ## `mod$spreadProb` are both masked to it, so a pixel outside cannot be ignited or spread
  ## into, while `fireReturnInterval` is not masked and overhangs the polygon. Without this,
  ## every such pixel inflates its zone's achieved interval -- on WesternAlbertaUpland that was
  ## 29.4% of flammable zone pixels, and made two zones look as though they under-burned by
  ## 3.5x and 12.2x when in fact every zone was within 0.90-1.07 of its target.
  sim$friSummary <- LandWebUtils::landmine_fri_summary(
    lthfc = lthfc,
    flammableMap = flammableMap,
    meanAnnualCumulBurnMap = meanAnnualCumulBurnMap,
    studyAreaName = studyAreaName,
    studyArea = sim$studyArea
  )

  ## per-zone attainment diagnostics: a CSV a developer can diff between runs, two figures,
  ## and a one-line verdict, so nobody has to read this log to learn a zone missed its target.
  sim$friDiagnostics <- LandWebUtils::landmine_fri_metrics(
    lthfc = lthfc,
    flammableMap = flammableMap,
    meanAnnualCumulBurnMap = meanAnnualCumulBurnMap,
    studyAreaName = studyAreaName,
    pixelAreaHa = prod(pixelRes) / 1e4,
    studyArea = sim$studyArea,
    nYears = end(sim) - start(sim)
  )

  fDiag <- file.path(outputPath(sim), "LandMine_FRI_diagnostics_multi.csv")
  fwrite(sim$friDiagnostics, fDiag)
  sim <- registerOutputs(fDiag)

  message(LandWebUtils::landmine_fri_verdict(sim$friDiagnostics))

  if ("png" %in% P(sim)$.plots) {
    fggFriZones <- file.path(figurePath(sim), "LandMine_FRI_zones_diagnostic_multi.png")
    ggsave(fggFriZones, LandWebUtils::landmine_plot_fri_zones(
      lthfc, sim$friDiagnostics, studyAreaName
    ), height = 7, width = 12)
    sim <- registerOutputs(fggFriZones)

    fggFriDrivers <- file.path(figurePath(sim), "LandMine_FRI_drivers_multi.png")
    ggsave(fggFriDrivers, LandWebUtils::landmine_plot_fri_drivers(
      sim$friDiagnostics, studyAreaName
    ), height = 5, width = 12)
    sim <- registerOutputs(fggFriDrivers)
  }

  f <- file.path(outputPath(sim), paste0("LandMine_FRI_summary_multi.csv"))
  fwrite(sim$friSummary, f)
  sim <- registerOutputs(f)

  ## LTHFC/FRI polygons
  ggFriPolys <- landmine_plot_LTHFC(lthfc, studyAreaName)

  if ("png" %in% P(sim)$.plots) {
    fggFriPolys <- file.path(figurePath(sim), "LandMine_LTHFC_map.png")
    png(fggFriPolys, height = 1000, width = 1000)
    print(ggFriPolys)
    dev.off()
    sim <- registerOutputs(fggFriPolys)
  }

  ## expected vs simulated fire return intervals
  ggFriExpVsSim <- landmine_plot_FRI(sim$friSummary) +
    geom_smooth(method = "lm")

  if ("png" %in% P(sim)$.plots) {
    fggFriExpVsSim <- file.path(figurePath(sim), "LandMine_FRI_exp_vs_sim.png")
    ggsave(fggFriExpVsSim, ggFriExpVsSim, height = 10, width = 10) ## NOTE: keep square aspect ratio
    sim <- registerOutputs(fggFriExpVsSim)
  }

  return(invisible(sim))
}

## .inputObjects
.inputObjects <- function(sim) {
  ## DEBUGGING: random seed issues
  # fseed <- file.path(outputPath(sim), "seed.txt")
  # writeEventInfo(sim, fseed, append = TRUE)
  # writeRNGInfo(fseed, append = TRUE)
  ## END DEBUGGING

  dPath <- asPath(inputPath(sim), 1)

  # Make random forest cover map
  mod$numDefaultPixelGroups <- 20L
  mod$numDefaultPolygons <- 4L
  numDefaultSpeciesCodes <- 2L

  if (!suppliedElsewhere("studyArea", sim)) {
    if (getOption("LandR.verbose", TRUE) > 0) {
      message("'studyArea' was not provided by user. Using a polygon in southwestern Alberta, Canada.")
    }

    sim$studyArea <- randomStudyArea(seed = 1234, size = 1e9)
  }

  if (!suppliedElsewhere("studyAreaReporting", sim)) {
    if (getOption("LandR.verbose", TRUE) > 0) {
      message("'studyAreaReporting' was not provided by user. Using the same as 'studyArea'.")
    }
    sim$studyAreaReporting <- sim$studyArea
  }

  if (is.null(sim$rasterToMatch)) {
    if (!suppliedElsewhere("rasterToMatch", sim)) {
      sim$rasterToMatch <- raster(sim$studyArea, res = 100)
      sim$rasterToMatch <- fasterize::fasterize(sf::st_as_sf(sim$studyArea), sim$rasterToMatch)
    }
  }
  ## TODO: use LandR::defineFlammable() below
  if (!suppliedElsewhere("flammableMap", sim)) {
    sim$flammableMap <- sim$rasterToMatch
    sim$flammableMap[] <- 1L  # 1 means flammable  ## TODO: use LandR::defineFlammable()
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
  nonFlammable <- which(is.na(sim[["flammableMap"]][]) | sim[["flammableMap"]][] == 0)
  if (length(nonFlammable) > 0) {
    sim$fireReturnInterval[nonFlammable] <- NA
  }

  if (!suppliedElsewhere(sim$ROSTable)) {
    ## ROS classes and values from Table 3.2 of Andison 1996
    ## - omitting 'water', 'non-productive brush', and 'non-productive black spruce' classes;
    ## - typo in Andison 1996: 'young mixed wood = 6' is really 'young hardwood = 6'.
    sim$ROSTable <- data.table::rbindlist(list(
      list("immature_young", "decid", 6L), ## aka hardwood
      list("mature", "decid", 9L), ## aka hardwood
      list("immature_young", "mixed", 12L),
      list("immature", "pine", 14L),
      list("mature", "mixed", 17L),
      list("immature_young", "softwood", 18L),
      list("immature_young", "spruce", 20L),
      list("mature", "pine", 21L),
      list("young", "pine", 22L),
      list("mature", "softwood", 27L),
      list("mature", "spruce", 30L)
    )) |>
      data.table::setnames(old = 1:3, new = c("age", "leading", "ros"))
  }

  if (!suppliedElsewhere("pixelGroupMap", sim)) {
    sim$pixelGroupMap <- Cache(randomPolygons, sim$rasterToMatch,
                               numTypes = mod$numDefaultPixelGroups)
  }

  if (!suppliedElsewhere("rstTimeSinceFire", sim)) {
    sim$rstTimeSinceFire <- terra::rast(sim$pixelGroupMap)
    sim$rstTimeSinceFire[] <- 200L
  }

  if (!suppliedElsewhere("species", sim)) {
    sim$species <- data.table(
      species = c("Pinu_sp", "Pice_gla"),
      speciesCode = 1:numDefaultSpeciesCodes
    )
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
    sppNames <- LandR::equivalentName(sim$species$species, sim$sppEquiv, column = P(sim)$sppEquivCol)
    sim$sppEquiv <- sim$sppEquiv[get(P(sim)$sppEquivCol) %in% sppNames]
  }

  return(invisible(sim))
}

