#' Core Burn function for Andison Fire Module
#'
#' @description The main function for the Andison Fire Module. See details.
#'
#' @param landscape      A \code{RasterLayer}. This only provides the extent and
#'                       resolution for the fire spread algorithm.
#' @param startCells     A numeric vector indicating the indices on the \code{landscape}
#'                       where the fires will start with 100\% certainty.
#' @param fireSizes      A numeric vector indicating the final size of each of the fires.
#'                       Must be same length as \code{startCells}.
#' @param nActiveCells1  A numeric vector of length 2. These are cutoffs above and
#'                       below each of which different values of \code{spawnNewActive}
#'                       are used. See details.
#' @param spawnNewActive A numeric vector of length 4. These are the probabilities
#'                       of creating spreading to 2 neighbours instead of the 1
#'                       default neighbour, each time step.
#'                       The 4 values are for 4 different fire size conditions. See details.
#' @param sizeCutoffs    A numeric vector of length 2, in hectares.
#'                       These are 2 size thresholds that affect which \code{spawnNewActive}
#'                       probabilities are used. See details.
#'
#' @details
#' This algorithm is a modified contagious cellular automaton.
#'
#' @section Algorithm:
#'
#' \subsection{Core}{
#' Each fire starts at a single pixel, \code{startCells} and will spread,
#' i.e., visit and convert from a 0 to the fire id number.
#' It will iteratively spread until the number of cells visited is equal to \code{floor(fireSizes)}.
#' }
#' \subsection{Adjustments due to current fire size and number of active pixels}{
#'
#' That can vary too, but it gets a bit complicated, so use that for now.
#' Spawning probability was originally set at 13%, but created problems with very
#' large and very small fires, so over time has been adjusted to vary depending on
#' a) number of active "firelets" (NF); and b) fire size (FS), such that:
#'   -         If 10 <= NF < 36 and FS < 20,000 ha then P = 20%
#'   -         If NF > 36 and FS < 8,000 ha, P = 11%
#'   -         If NF < 36 and FS > 20,000 ha, P = 26%
#'   -         If NF < 10 then P = 46%
#' These rule create more heterogeneity in the pattern of burning.
#' }
#'
#' \subsection{Fire jumping}:
#' If the fire has not reached its target size, it will try to pick new neighbours among
#' the 8 immediate neighbours up to 4 times. If it still did not find enough neighbours, then
#' it will jump or "spot" up to 4 pixels away. It will then repeat the previous 2 stages again
#' once (i.e., 4 neighbours, 1 jump, 4 neighbours, 1 jump), then it will stop, unable to achieve
#' the desired fireSize.
#'
#' @return A \code{data.table} with 4 columns
#'
#' @export
#'
# burn <- compiler::cmpfun(function(landscape, startCells, fireSizes = 5, nActiveCells1 = c(10, 36),
#                  spawnNewActive = c(0.46, 0.2, 0.26, 0.11),
#                  sizeCutoffs = c(8e3, 2e4), spreadProbRel = 0.23) {
#   a <- spread(landscape, loci = startCells, spreadProbRel = spreadProbRel, persistence = 0,
#               neighProbs = c(1-spawnNewActive[1], spawnNewActive[1]), iterations = 1,
#               mask = NULL, maxSize = fireSizes, directions = 8, returnIndices = TRUE,
#               id = TRUE, plot.it = FALSE, exactSizes = TRUE)
#
#   while(sum(a$active) > 0) {
#     b <- a[,list(numActive = sum(active), size = .N), by = id]
#     ## This is undescribed in Andison -- NF >36 & FS >8,000 ha --
#     ## They look too circular, without this, so make this zero, no new spawn
#     set(b, NULL, "pSpawnNewActive", 0)
#     #set(b, NULL, "pSpawnNewActive", spawnNewActive[1])
#     b[numActive < nActiveCells1[1], pSpawnNewActive := spawnNewActive[1]]
#     b[numActive >= nActiveCells1[1] & numActive<nActiveCells1[2] & size < sizeCutoffs[2], pSpawnNewActive:=spawnNewActive[2]]
#     b[numActive > nActiveCells1[2] & size < sizeCutoffs[1], pSpawnNewActive := spawnNewActive[4]]
#     b[numActive < nActiveCells1[2] & size > sizeCutoffs[2], pSpawnNewActive := spawnNewActive[3]]
#     set(b, , "pNoNewSpawn", 1-b$pSpawnNewActive)
#
#     ## spawnNewActive must be joined sent in here as list...
#     b <- b[a]
#     a <- spread(landscape, spreadProbRel = spreadProbRel, spreadProb = spreadProb,
#                 spreadState = a, persistence = 0,
#                 neighProbs = transpose(as.list(b[active == TRUE, c("pNoNewSpawn", "pSpawnNewActive")])),
#                 iterations = 1, quick = TRUE,
#                 mask = NULL, maxSize = fireSizes, directions = 8, returnIndices = TRUE,
#                 id = TRUE, plot.it = FALSE, exactSizes = TRUE)
#   }
#   return(a)
# })


#' @param sizeCutoffs Numeric vector, length 2, indicating in the units of the map (e.g., m), not
#'                    pixels, the cutoff between the different \code{spawnNewActive} levels.
#'
#' @export
#' @importFrom data.table set
#' @importFrom purrr transpose
#' @importFrom raster res
#' @importFrom SpaDES.tools spread2
burn1 <- compiler::cmpfun(function(landscape, startCells, fireSizes = 5, nActiveCells1 = c(10, 36),
                                   spawnNewActive = c(0.46, 0.2, 0.26, 0.11), maxRetriesPerID = 10,
                                   sizeCutoffs = c(8e3, 2e4), spreadProbRel = spreadProbRel, spreadProb = 0.77) {
  # convert to pixels
  sizeCutoffs <- sizeCutoffs / (prod(res(landscape)) / 1e4)
  # add a little bit of stochasticity
  #spawnNewActive <- 10^rnorm(length(spawnNewActive), log10(spawnNewActive), 0.1)

  a <- spread2(landscape, start = startCells, spreadProb = 1,  # initial step can have spreadProb 1 so guarantees something
               spreadProbRel = spreadProbRel, #persistence = 0,
               neighProbs = c(1 - spawnNewActive[1], spawnNewActive[1]), iterations = 1,
               #mask=NULL,
               maxRetriesPerID = maxRetriesPerID,
               asRaster = FALSE, exactSize = fireSizes, directions = 8, #returnIndices = TRUE,
               #id = TRUE, plot.it = FALSE
               )
  whActive <- attr(a, "spreadState")$whActive#a$state=="activeSource"
  while (any(whActive)) {
    #a2 <- data.table::copy(a)
    # microbenchmark(times = 6, a = {
    #   a <- data.table::copy(a2)
    #   b <- a[,list(numActive = sum(state=="activeSource"), size = .N),by=initialPixels]
    #
    #  }, b = {
    #   a <- data.table::copy(a2)
       set(a, NULL, "numActive", 0L)
       a[whActive, numActive := .N, by = initialPixels]
       b <- attr(a, "spreadState")$cluster
       b <- a[b, mult = "last"]
       set(b, NULL, c("numRetries", "maxSize", "exactSize", "id", "state", "pixels"), NULL)
       set(a, NULL, c("numActive"), NULL)
    #
    # }, d = {
    #   a <- data.table::copy(a2)
    #   b <- a[,list(numActive = sum(state=="activeSource")),by=initialPixels]
    #   set(b, , "size", attr(a, "spreadState")$cluster$size)
    # }, f = {
    #   a <- data.table::copy(a2)
    #   b <- a[,list(numActive = sum(state=="activeSource")),by=initialPixels]
    #
    # })
    #
    set(b, NULL, "pSpawnNewActive", spawnNewActive[1])

    b[numActive >= nActiveCells1[1] & numActive < nActiveCells1[2] &
        size < sizeCutoffs[2], pSpawnNewActive := spawnNewActive[2]]
    b[numActive > nActiveCells1[2] & size < sizeCutoffs[1], pSpawnNewActive := spawnNewActive[4]]
    b[numActive < nActiveCells1[2] & size > sizeCutoffs[2], pSpawnNewActive := spawnNewActive[3]]
    set(b, NULL, "pNoNewSpawn", 1 - b$pSpawnNewActive)
    set(b, NULL, c("numActive"), NULL)

    # spawnNewActive must be joined sent in here as list...
    b <- b[a]#, on="initialPixels"]
    a <- spread2(landscape, spreadProbRel = spreadProbRel, spreadProb = spreadProb,
                 start = a, #persistence = 0,
                 neighProbs = transpose(as.list(b[state == "activeSource", c("pNoNewSpawn", "pSpawnNewActive")])),
                 iterations = 1, skipChecks = TRUE, asRaster = FALSE,
                 exactSize = attr(a, "spreadState")$cluster$maxSize,
                 #mask = NULL,
                 #maxSize = fireSizes,
                 maxRetriesPerID = maxRetriesPerID,
                 directions = 8, #returnIndices = TRUE,
                 #id = TRUE, plot.it = FALSE
                 )
    #message("max size:", max(attr(a, "spreadState")$clusterDT$maxSize), ", current size:",max(attr(a, "spreadState")$clusterDT$size))
    set(a, NULL, "order", seq_len(NROW(a)))
    whActive <- attr(a, "spreadState")$whActive#a$state=="activeSource"
  }
  return(a)
})
