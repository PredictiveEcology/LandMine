---
title: "LandMine Manual"
subtitle: "v.1.0.0"
date: "18 May 2018; updated 6 Sep 2022"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    theme: sandstone
    number_sections: false
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
  bibliography: citations/references_LandMine.bib
citation-style: citations/ecology-letters.csl
link-citations: true
always_allow_html: true
---

# *LandMine* Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:LandMine) *LandMine*



:::{.rmdwarning}
This documentation is work in progress.
Please report any discrepancies or omissions at <https://github.com/PredictiveEcology/LandMine/issues>.
:::

#### Authors:

Eliot J B McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut], Alex M. Chubaty <achubaty@for-cast.ca> [ctb, cre]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

Landmine is a model created for simulating the natural range of variation for landscapes in the boreal forest [@Andison:1996; @Andison:1998].
It has been widely used by the public and the private sector for various purposes.
This SpaDES module is a rewrite of the fire component in native R.

### Model Differences

The current version has not yet been fully tested and compared with the original version, but there are currently several known differences:

1. Fire sizes are taken from a Truncated Pareto distribution, resulting in numerous very small fires, and few large fires;
2. Parameters have been fitted to the landscapes that are under study in the LandWeb project.

### Known Species

Landmine requires the following codes as inputs (the genus and species codes below), which converts and groups species as follows.
Each of the species groups has its own Rate of Spread (ROS) for fire spreading:


Table: (\#tab:species-table)LandMine species codes.

|Species                  |Code     |Group                   |
|:------------------------|:--------|:-----------------------|
|Jack pine                |Pinu_ban |Pine (PINU)             |
|Lodgepole pine           |Pinu_con |Pine (PINU)             |
|Unspecified pine species |Pinu_spp |Pine (PINU)             |
|Paper birch              |Betu_pap |Deciduous (DECI)        |
|Balsam poplar            |Popu_bal |Deciduous (DECI)        |
|Trembling aspen          |Popu_tre |Deciduous (DECI)        |
|Larch/Tamarack           |Lari_lar |Deciduous (DECI)        |
|Black spruce             |Pice_mar |Black spruce (PICE_MAR) |
|White spruce             |Pice_gla |White spruce (PICE_GLA) |
|Fir and other softwoods  |Abie_spp |Fir (ABIE)              |

\newpage

### Module inputs and parameters

Table \@ref(tab:moduleInputs-LandMine) shows the full list of module inputs.

<table class="table" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-LandMine)(\#tab:moduleInputs-LandMine)List of (ref:LandMine) input objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
   <th style="text-align:left;"> sourceURL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cohortData </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Columns: B, pixelGroup, speciesCode (as a factor of the names), age. indicating several features about the current vegetation of stand. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireReturnInterval </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> A raster layer that is a factor raster, with at least 1 column called `fireReturnInterval`, representing the fire return interval in years. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pixelGroupMap </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Pixels with identical values share identical stand features </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatch </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> a raster of the `studyArea` to use as a template raster (resolution, projection, etc.) for all other rasters in the simulation. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ROSTable </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> A data.table with 3 columns, 'age', 'leading', and 'ros'. The values under the 'age' column can be 'mature', 'immature', 'young' and compound versions of these, e.g., 'immature_young' which can be used when 2 or more age classes share same 'ros'. 'leading' should be vegetation type. 'ros' gives the rate of spread values for each age and type. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstFlammable </td>
   <td style="text-align:left;"> Raster </td>
   <td style="text-align:left;"> A raster layer, with 0, 1 and NA, where 1 indicates areas that are flammable, 0 not flammable (e.g., lakes) and NA not applicable (e.g., masked) </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstTimeSinceFire </td>
   <td style="text-align:left;"> Raster </td>
   <td style="text-align:left;"> a time since fire raster layer </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> species </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Columns: species, speciesCode, Indicating several features about species </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppColorVect </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> named character vector of hex colour codes corresponding to each species </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquiv </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Multi-columned data.table indicating species name equivalencies. Default taken from `LandR::sppEquivalencies_CA` which has names for species of trees in Canada </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyArea </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> multipolygon, typically buffered around an area of interest (i.e., `studyAreaReporting`) to use for simulation. Defaults to an area in Southwestern Alberta, Canada. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaReporting </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> multipolygon (typically smaller/unbuffered than `studyArea`) to use for plotting/reporting. Defaults to an area in Southwestern Alberta, Canada. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Summary of user-visible parameters (Table \@ref(tab:moduleParams-LandMine)).

<table class="table" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-LandMine)(\#tab:moduleParams-LandMine)List of (ref:LandMine) parameters and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> paramName </th>
   <th style="text-align:left;"> paramClass </th>
   <th style="text-align:left;"> default </th>
   <th style="text-align:left;"> min </th>
   <th style="text-align:left;"> max </th>
   <th style="text-align:left;"> paramDesc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> biggestPossibleFireSizeHa </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1e+06 </td>
   <td style="text-align:left;"> 10000 </td>
   <td style="text-align:left;"> 2e+06 </td>
   <td style="text-align:left;"> An upper limit, in hectares, of the truncated Pareto distribution of fire sizes </td>
  </tr>
  <tr>
   <td style="text-align:left;"> burnInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time at which the first burn event should occur </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireTimestep </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time interval between burn events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> maxReburns </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 1, 20 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 20 </td>
   <td style="text-align:left;"> Number of attempts to burn fires that don't reach their target fire size. Reburning occurs in two phases, hence accepting a parameter value of length 2. In the first phase, fires that did not reach their target size are reignited from new pixels within the FRI zone, so they are less likely to continue being stuck in a region with sinuous fires or discontinuous fuels. If, after `maxReburns[1]` attempts, there are still fires that haven't reached their target size, the second reburn phase is attempted. After recording the the pixels that did burn in phase one, new fires are ignited, whose target sizes are set equal the difference between the previous target and the previously burned area. Repeats up to `maxReburns[2]` times. This results in additional (smaller) fires, but since the purpose of LandMine is to replicate area burned per year to achieve LTHFC, this is an acceptable compromise. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> maxRetriesPerID </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 299 </td>
   <td style="text-align:left;"> Number of additional attempts ('jumps') that will be made per firelet ID, before abandoning. See `?SpaDES.tools::spread2`. NOTE: increasing this value results in longer simulation times when firelets get 'stuck', but higher values are needed to achive larger fire sizes with discontinuous fuels. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> minPropBurn </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.9 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> Minimum proportion burned pixels to use when triggering warnings about simulated fires. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mixedType </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> How to define mixed stands: 1 for any species admixture; 2 for deciduous &gt; conifer. See ?vegTypeMapGenerator. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mode </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> single </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> use 'single' to run part of a landscape simulation; use 'multi' to run as part of postprocessing multiple simulation runs. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> optimParsRowID </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> which set of optimization parameter values to use for simulating fire spread, specified by row number of the `LandMine_DEoptim_params.csv` file. `1L` specifies the original 2018 values at 100m pixels; all other rows were calculated using 250m pixels. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reps </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> number of replicates/runs per study area when running in 'multi' mode. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ROSother </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 30 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> default ROS value for non-forest vegetation classes.this is needed when passing a modified `ROSTable`, e.g. using log-transformed values. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ROStype </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> default </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> One of 'burny', 'equal', 'log', or 'default'. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquivCol </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> LandR </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The column in `sim$specieEquivalency` data.table to use as a naming convention. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> useSeed </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Only used for creating a starting `cohortData` dataset. If `NULL`, then it will be randomly generated; If non-`NULL`, will pass this value to `set.seed()` and be deterministic and identical each time. WARNING: setting the seed to a specific value will cause all simulations to be identical! </td>
  </tr>
  <tr>
   <td style="text-align:left;"> vegLeadingProportion </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.8 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> a number that define whether a species is leading for a given pixel </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time at which the first plot event should occur </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time interval between plot events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plots </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> png, screen </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Passed to `types` in `Plots` (see `?Plots`). There are a few plots that are made within this module, if set. Note that plots (or their data) saving will ONLY occur at `end(sim)`. If `NA`, plotting is turned off completely (this includes plot saving). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time at which the first save event should occur </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time interval between save events </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .studyAreaName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> test </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Human-readable name for the study area used - e.g., a hash of the studyarea obtained using `reproducible::studyAreaName()` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .unitTest </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Some functions can have internal testing. This will turn those on or off, if any exist. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useParallel </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Used in burning. Will be passed to `data.table::setDTthreads()`. NOTE: use `.useParallel &lt;= 2` as the additonal RAM overhead too high given marginal speedup. </td>
  </tr>
</tbody>
</table>

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-LandMine)).

<table class="table" style="color: black; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-LandMine)(\#tab:moduleOutputs-LandMine)List of (ref:LandMine) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> fireInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> The initial event time of the burn event. This is simply a reassignment from `P(sim)$burnInitialTime`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSizes </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> A list of data.tables, one per burn event, each with two columns, `size` and `maxSize`. These indicate the actual sizes and expected sizes burned, respectively. These can be put into a single data.table with `rbindlist(sim$fireSizes, idcol = 'year')` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireReturnInterval </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> A `Raster` map showing the fire return interval. This is created from the `rstCurrentBurn`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireReturnIntervalsByPolygonNumeric </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> A vector of the fire return intervals, ordered by the numeric representation of polygon ID </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireTimestep </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> The number of time units between successive fire events in a fire module. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> friSummary </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> summary fire return interval table </td>
  </tr>
  <tr>
   <td style="text-align:left;"> kBest </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> A numeric scalar that is the optimal value of `K` in the Truncated Pareto distribution (`rtruncpareto`) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> numFiresPerYear </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> The average number of fires per year, by fire return interval level on `rstCurrentBurn`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstCurrentBurn </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> A raster layer, produced at each timestep, where each pixel is either 1 or 0 indicating burned or not burned. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstCurrentBurnCumulative </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Cumulative number of times a pixel has burned </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquiv </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Same as input, but with new column, `LandMine`. </td>
  </tr>
</tbody>
</table>

\newpage

## Usage

To run this Landmine module alone (*i.e.*, for fitting), the following should work (*iff* raster inputs for `studyArea` and `rasterToMatch` are available), assuming all R packages are available.

**NB:** Paths will have be changed for a different user.

### Package dependencies

To determine which packages are used by `LandMine`, use:


``` r
SpaDES.core::packages(modules = "LandMine", paths = "..")[[1]]
```

```
##  [1] "SpaDES.core"                                               
##  [2] "assertthat"                                                
##  [3] "data.table"                                                
##  [4] "fpCompare"                                                 
##  [5] "ggplot2"                                                   
##  [6] "ggspatial"                                                 
##  [7] "grDevices"                                                 
##  [8] "gridExtra"                                                 
##  [9] "magrittr"                                                  
## [10] "PredictiveEcology/LandR@development (>= 1.1.0.9003)"       
## [11] "PredictiveEcology/LandWebUtils@development (>= 1.0.3.9002)"
## [12] "PredictiveEcology/pemisc@development"                      
## [13] "PredictiveEcology/SpaDES.tools@development"                
## [14] "RColorBrewer"                                              
## [15] "stats"                                                     
## [16] "terra"                                                     
## [17] "VGAM"
```

### Module usage

First, define a study area and create a template raster.


``` r
studyArea <- SpaDES.tools::randomStudyArea(seed = 1234, size = 1e10)
rasterToMatch <- terra::rast(studyArea, resolution = 240)
```

Next, set up the simulation to run for 13 timesteps using default module parameters.


``` r
times <- list(start = 0, end = 13)

parameters <- list()

modules <- list("LandMine")

objects <- list(
  studyArea = studyArea,
  rasterToMatch = rasterToMatch
)

paths <- list(
  cachePath = cacheDir,
  modulePath = moduleDir,
  inputPath = inputDir,
  outputPath = outputDir
)
```


``` r
mySim <- simInit(
  times = times,
  params = parameters,
  modules = modules,
  objects = objects,
  paths = paths
)

dev()
mySimOut <- spades(mySim, .plotInitialTime = times$start, debug = TRUE)
```

\newpage

## Testing the burn algorithm


``` r
Require(c("data.table", "DEoptim", "parallel"))

s <- simInit(
  times = times, 
  params = parameters, 
  modules = modules,
  objects = objects, 
  paths = paths
)
```

\newpage

## Optimizing parameters

The following code chunk tries to find values of `spawnNewActive` that creates "reasonable" fire shapes at all sizes.


``` r
pixelSize <- 240
ros <- terra::rast(
  terra::ext(0, pixelSize * 1e3, 0, pixelSize * 1e3),
  resolution = pixelSize,
  vals = 0
)
ros <- ros == 0
fireSize <- 1e5

maxRetriesPerID <- 4 ## 4 retries (5 attempts total)
spreadProb <- 0.9
spawnNewActive <- c(0.46, 0.2, 0.26, 0.11)
sizeCutoffs <- c(8e3, 2e4)

NineCorners <- terra::cellFromRowCol(
  ros,
  row = nrow(ros) / 4 * rep(1:3, 3),
  col = ncol(ros) / 4 * rep(1:3, each = 3)
)

centreCell <- terra::cellFromRowCol(
  ros,
  row = nrow(ros) / 2,
  col = ncol(ros) / 2
)

## Set variables
objs <- c("ros", "centreCell", "fireSize", "spawnNewActive", "sizeCutoffs", "spreadProb")
pkgs <- c("data.table", "LandWebUtils", "raster", "SpaDES.tools")
```


``` r
### SET UP CLUSTER FOR PARALLEL
wantParallel <- TRUE

## numCores should be >= 70 and needs to be multiple of number of params to be fit (7)
numCores <- (parallelly::availableCores(constraints = c("connections")) %/% 7) * 7

machineName <- strsplit(Sys.info()["nodename"], "[.]")[[1]][1]

clNames <- switch(machineName,
                  pinus = c(
                    rep("localhost", 0),
                    rep("picea.for-cast.ca", 20),
                    rep("pseudotsuga.for-cast.ca", 82)
                  ),
                  picea = rep("localhost", numCores),
                  rep("localhost", numCores))

if (wantParallel) {
  cl <- landmine_optim_clusterSetup(nodes = clNames)
} else {
  cl <- NULL
}

if (!inherits(cl[[1]], "forknode")) {
  landmine_optim_clusterExport(cl, objs = objs, pkgs = pkgs)
}
```

### Visual examination


``` r
opt_sn <- DEoptim(
  landmine_optim_fitSN,
  lower = c(-2, -3, -3, -3, 1, 3.5, 0.75), 
  upper = c(-0.1, -0.5, -0.5, -1, 3.5, 5, 1),
  control = DEoptim.control(
    VTR = 0.001,
    NP = as.integer(length(clNames)),
    itermax = 200,
    cluster = cl,
    strategy = 6
  ),
  ros = ros,
  centreCell = centreCell,
  fireSizes = c(10, 100, 1000, 10000, 100000),
  desiredPerimeterArea = 0.003
)

opt_sn$optim$bestmem ## best param values

## assign with suffix to facilitate multiple DEoptim runs
optimParams <- data.frame(
  date = format(Sys.time(), "%Y-%m-%d"),
  pixelSize = pixelSize
) |>
  cbind(rbind(opt_sn$optim$bestmem))

fDEoptim <- file.path(moduleDir, "LandMine", "data", "LandMine_DEoptim_params.csv")
prevVals <- read.csv(fDEoptim)
write.csv(rbind(prevVals, optimParams), fDEoptim, row.names = FALSE)

parallel::stopCluster(cl)
cl <- NULL
```


``` r
fs_sn <- c(10, 100, 1000, 10000, 100000)
fit_sn <- landmine_optim_fitSN(
  sna = c(-1, -1, -1, -2, 2, 4, 0.9),
  ros = ros,
  centreCell = centreCell,
  fireSizes = fs_sn,
  desiredPerimeterArea = 0.003
)
bfs1_sn <- purrr::transpose(attr(fit_sn, "bfs1"))
LM_sn <- do.call(rbind, bfs1_sn$LM)
plot(log10(fs_sn), LM_sn[, "perim.area.ratio"]) ## NOTE: visual inspection - not too round; not too sinuous
```

## Alternate optimization

A second (alternative) version tries the optimization using fewer parameters, to test whether a simpler version gets better/different results.
Although this version was not used for the final module, we preserve it here for posterity.


``` r
fs_optim2 <- c(0.2, 1:8) * 10000
opt_sn2 <- DEoptim(
  landmine_optim_fitSN2,
  lower = c(1, -1, 1, 3, 4),
  upper = c(3, -0.3, 3, 4, 5), 
  control = DEoptim.control(
    VTR = 0.001, 
    itermax = 40, 
    cluster = cl,
    strategy = 6
  ), 
  ros = ros, 
  centreCell = centreCell, 
  fireSizes = fs_optim2, 
  desiredPerimeterArea = 0.003
)
```

### Visual examination


``` r
fs_sn2 <- round(runif(10, 10, 4000))
fit_sn2 <- landmine_optim_fitSN2(
  par = c(2, -0.63333, 1, 3.2, 4.4),
  ros = ros,
  centreCell = centreCell,
  fireSizes = fs_sn2,
  desiredPerimeterArea = 0.003,
  spreadProb = 0.9
)
bfs1_sn2 <- purrr::transpose(attr(fit_sn2, "bfs1"))
LM_sn2 <- do.call(rbind, bfs1_sn2$LM)

## NOTE: visual inspection - not too round; not too sinuous
plot(log10(fs_sn2), LM_sn2[, "perim.area.ratio"]) 
```

\newpage

## Manual inspection of optimization results

### Original (2018) version

The original version was run using 100m pixels, despite the simulations being run using 250m pixels.
This was corrected and rerun below.


``` r
## 10,000 hectares burns gave this
spawnNewActive[2:3] <- c(0.0235999945606232, 0.0263073265505955)

#100,000 hectare burns gave this
#spawnNewActive <- 10^c(-1.264780,   -1.970946,   -1.288213,   -2.202580)
spawnNewActive <- 10^c(-0.586503645288758, -1.08108837273903,
                       -2.14391896536108, -1.00221184641123)
sizeCutoffs <- 10^c(3.37711253212765, 4.52040993282571)

sns <- c(-1.733262, -0.933318, -2.562183, -2.493687, 3.064458, 4.812305)
spawnNewActive <- 10^sns[1:4]
sizeCutoffs <- 10^sns[5:6]
#spawnNewActive <- 10^c(-1.646419, -1.815395, -2.809013, -2.613337)
#sizeCutoffs <- 10^c(3.888317,    4.641961)

## 100,000 pixel fires -- the next worked, but I think we can get better
# sns <- structure(
#   c(-1.652459,-0.962121,-0.964879,-2.304902, 3.522345, 4.173242),
#   .Names = c("par1", 
sns <- structure(
  c(-1.641197,-1.152821,-0.697335,-1.751917, 3.720378, 4.034059),
  .Names = c("par1", "par2", "par3", "par4", "par5", "par6")
)
spawnNewActive <- 10^sns[1:4]
sizeCutoffs <- 10^(sns[5:6])
fireSize <- 30000

## 100
sns <- c(-0.77716149196811, -0.769325340166688, -1.2772046867758,
         -1.99332102853805, 3.14260408212431, 4.46155184064992) 

## 1000
sns <- c(-0.775107,-1.031760,-0.599669,-1.958105, 3.048958, 4.275831)

## seemed good for 100,000, pretty good for 1e3
sns <- c(-1.54885, -0.97052, -1.38305, -1.93759, 3.20379, 4.13237)

## good for 100 000, 10 000 ha -- too sinuous for 1000 and 100 ha
sns <- c(-1.537203,-1.462981,-0.524957,-1.002567, 3.642046, 4.501754)

## good for 100 000, 10 000 ha (except some fires @ 1e5 don't make it to full size)
## -- too sinuous for smaller
sns <- c(-1.484338,-1.220440,-2.948275,-2.155940, 3.945281, 4.904893)

sns <- c(-1.495247,-0.800494,-1.582350,-2.270646, 3.530671, 4.663245)

## final optimization after 75 iterations, Good: 1e5, 1e4
sns <- c(-1.47809, -0.86224, -1.34532, -1.93568, 3.27149, 4.20741)

## based on equal weights 10^(1:5)
sns <- c(-0.923528, -1.804549, -1.760455, -1.793594,  1.683355,  4.466668)

## With spreadProb = 0.9 # Pretty GOOD!
sns <- c(-0.731520, -0.501823, -0.605968, -1.809726,  2.202732,  4.696060, 0.9) ## used in module
optimParams <- data.frame(
  date = "2018-05-18",
  pixelSize = 100
) |>
  cbind(rbind(sns))
write.csv(optimParams, fDEoptim, row.names = FALSE)

## With spreadProb = 0.9 # Optimal
sns <- c(-0.978947, -0.540946, -0.790736, -1.583039,  2.532013,  4.267547,  0.946730)

spawnNewActive <- 10^sns[1:4]
sizeCutoffs <- 10^(sns[5:6])
if (length(sns) == 7) spreadProb <- sns[7]

# from linear model version
par <- c(1.548899,-0.396904, 2.191424, 3.903082, 4.854002)
sizeCutoffs <- 10^c(par[4], par[5])
sna <- min(-0.15, par[1] + par[2]*log10(fireSize))
sna <- 10^c(sna*par[3], sna*2*par[3], sna*3*par[3], sna*4*par[3])
spawnNewActive <- sna

###########################
dev.new()

for (i in 1:5) {
  fireSize <- 10^i
  dim <- round(sqrt(fireSize)*5 * 250)
  ros <- raster(extent(0, dim,0,dim), res = 250, vals = 1)
  centreCell <- cellFromRowCol(ros,
                               rownr = nrow(ros) / 2,
                               colnr = ncol(ros) / 2)
  reps <- paste0("rep", 1:4 + (log10(fireSize) - 1)*4)
  burnedMapList <- landmine_optim_clusterWrap(cl = cl, nodes = clNames, reps = reps, objs = objs, pkgs = pkgs)
  names(burnedMapList$out) <- reps
  burnedMapList <- purrr::transpose(burnedMapList$out)
  cl <- burnedMapList$cl
  do.call(rbind, burnedMapList$LM)
  
  terra::plot(
    burnedMapList$burnedMap,
    col = "red",
    colNA = "white",
    legend = FALSE,
    main = paste0("Fire size: ", fireSize)
  )
}
```


``` r
reps <- paste0("rep", 1:1)
perims <- list()
perm <- list()
mod <- list()

dev.new()

fireSizes <-  10^(4)

for (fs in fireSizes) {
  for (i in 1:1) {
    ros <- terra::rast(terra::ext(0, 2e5, 0, 2e5), resolution = 240, vals = 1)
    NineCorners <- terra::cellFromRowCol(
      ros,
      rownr = nrow(ros) / 4 * rep(1:3, 3),
      colnr = ncol(ros) / 4 * rep(1:3, each = 3)
    )
    centreCell <- NineCorners
    ran <- runif(4, -3, -1)
    spawnNewActive <- 10^ran
    # spawnNewActive <- 10^c(-0.1, -0.75, -1.2, ran*2.5)
    fireSize = rep(fs, length(centreCell))
    sizeCutoffs <- 10^c(1,3)
    burnedMapList <- landmine_optim_clusterWrap(cl = cl, nodes = clNames, reps = reps, objs = objs, pkgs = pkgs)
    names(burnedMapList$out) <- reps
    burnedMapList <- purrr::transpose(burnedMapList$out)
    cl <- burnedMapList$cl
    do.call(rbind, burnedMapList$LM)
    terra::plot(burnedMapList$burnedMap)
    perims[[i]] <- data.frame(
      perim = burnedMapList$LM$rep1$perim.area.ratio,
      spawnNewActive = mean(spawnNewActive),
      others = t(spawnNewActive)
    )
  }
  fsChar <- as.character(fs)
  perm[[fsChar]] <- rbindlist(perims)
  perm[[fsChar]]$perim <- log10(perm[[fsChar]]$perim)
  perm[[fsChar]]$spawnNewActive <- perm[[fsChar]]$spawnNewActive
  plot(perm[[fsChar]]$perim, perm[[fsChar]]$spawnNewActive, add = TRUE)
  mod[[fsChar]] <- lm(spawnNewActive ~ perim, data = perm[[fsChar]])
}

(predict(mod[["10"]], data.frame(perim = log10(0.003))))
(predict(mod[["100"]], data.frame(perim = log10(0.003))))
log10(predict(mod[["1000"]], data.frame(perim = log10(0.003))))
```

### Current (2025) version

The original version was run using 100m pixels, despite the simulations being run using 250m pixels.

This version uses 240m pixels to work with LandWeb v3.


``` r
pixelSize <- 240

## final optimization after 200 iterations (landmine_optim_fitSN)
optimParams <- read.csv(fDEoptim)
sns <- optimParams[, grepl("^par", colnames(optimParams))] |> tail(1) ## take last row

spawnNewActive <- 10^sns[1:4]
sizeCutoffs <- 10^(sns[5:6])
spreadProb <- sns[7]

## from linear model version
par <- c(1.548899,-0.396904, 2.191424, 3.903082, 4.854002)
sizeCutoffs <- 10^c(par[4], par[5])
sna <- min(-0.15, par[1] + par[2]*log10(fireSize))
sna <- 10^c(sna*par[3], sna*2*par[3], sna*3*par[3], sna*4*par[3])
spawnNewActive <- sna

for (i in 1:5) {
  fireSize <- 10^i
  dim <- round(sqrt(fireSize)*5 * pixelSize)
  ros <- terra::rast(terra::ext(0, dim, 0, dim), resolution = pixelSize, vals = 1)
  centreCell <- terra::cellFromRowCol(
    ros,
    row = nrow(ros) / 2,
    col = ncol(ros) / 2
  )
  reps <- paste0("rep", 1:4 + (log10(fireSize) - 1)*4)
  burnedMapList <- landmine_optim_clusterWrap(
    cl = cl, 
    nodes = clNames, 
    reps = reps, 
    objs = objs, 
    pkgs = pkgs
  )
  names(burnedMapList$out) <- reps
  burnedMapList <- purrr::transpose(burnedMapList$out)
  cl <- burnedMapList$cl
  do.call(rbind, burnedMapList$LM)
  
  terra::plot(
    burnedMapList$burnedMap,
    col = "red",
    colNA = "white",
    legend = FALSE,
    main = paste0("Fire size: ", fireSize)
  )
}
```

### Cleaning up


``` r
parallel::stopCluster(cl)
unlink(cacheDir, recursive = TRUE)
```

### Code and data availability

Code available from <https://github.com/PredictiveEcology/LandMine>.

### Links to other modules

Originally developed as part of the [LandWeb](https://github.com/PredictiveEcology/LandWeb) project.
