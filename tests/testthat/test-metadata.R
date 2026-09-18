## The module's metadata is its public contract: a project using LandMine binds to these
## object names and classes, and `reqdPkgs` states what it needs to run at all. These are
## CHARACTERIZATION tests -- they pin today's contract so a change to it has to be
## deliberate, rather than describing behaviour that did not exist before.
##
## When a change is intended, update this file in the same commit and bump the module
## version to match: removed, renamed or retyped is a MAJOR bump.

test_that("module metadata parses", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_type(md, "list")
  expect_identical(md$name, moduleName)
})

test_that("inputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  inputs <- stats::setNames(md$inputObjects$objectClass, md$inputObjects$objectName)
  expect_identical(
    inputs[order(names(inputs))],
    c(cohortData         = "data.table",
      fireReturnInterval = "SpatRaster",
      flammableMap       = "SpatRaster",
      pixelGroupMap      = "SpatRaster",
      rasterToMatch      = "SpatRaster",
      ROSTable           = "data.table",
      rstTimeSinceFire   = "SpatRaster",
      species            = "data.table",
      sppColorVect       = "character",
      sppEquiv           = "data.table",
      studyArea          = "SpatVector",
      studyAreaReporting = "sf")
  )
})

test_that("outputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  outputs <- stats::setNames(md$outputObjects$objectClass, md$outputObjects$objectName)
  expect_identical(
    outputs[order(names(outputs))],
    c(burnMap                             = "SpatRaster",
      fireInitialTime                     = "numeric",
      fireReturnInterval                  = "SpatRaster",
      fireReturnIntervalsByPolygonNumeric = "numeric",
      fireSizes                           = "list",
      fireTimestep                        = "numeric",
      friDiagnostics                      = "data.table",
      friSummary                          = "data.table",
      kBest                               = "numeric",
      numFiresPerYear                     = "numeric",
      rstCurrentBurn                      = "SpatRaster",
      sppEquiv                            = "data.table")
  )
})

test_that("parameters are the expected names", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_identical(
    sort(md$parameters$paramName),
    sort(c(".plotInitialTime", ".plotInterval", ".plots", ".saveInitialTime",
           ".saveInterval", ".studyAreaName", ".unitTest", ".useCache", ".useParallel",
           "biggestPossibleFireSizeHa", "burnInitialTime", "fireTimestep", "maxReburns",
           "maxRetriesPerID", "minPropBurn", "mixedType", "mode", "optimParsRowID",
           "reps", "ROSother", "ROStype", "sppEquivCol", "useSeed",
           "vegLeadingProportion"))
  )
})

test_that("reqdPkgs pins a minimum LandWebUtils", {
  ## LandMine's real logic lives in LandWebUtils, so the module is only correct when
  ## paired with a new enough copy of it -- `landmine_attach_identity()` and
  ## `landmine_fire_attainment()` arrived in 1.0.3.9038. A floor that silently loses its
  ## version bound is the failure this guards: the module would then resolve against
  ## whatever LandWebUtils happened to be installed.
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  pkgs <- unlist(md$reqdPkgs)
  lwu <- grep("LandWebUtils", pkgs, value = TRUE)

  expect_length(lwu, 1L)
  expect_match(lwu, "\\(>=\\s*[0-9.]+\\)", info = "LandWebUtils must carry a version floor")
  expect_gte(
    package_version(sub(".*\\(>=\\s*([0-9.]+)\\).*", "\\1", lwu)),
    package_version("1.0.3.9038")
  )
})
