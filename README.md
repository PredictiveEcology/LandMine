
# Overview

Landmine is a model created for simulating the natural range of variation for landscapes in the boreal forest.
Written in the 1990s (Andison 1996), it has been widely used by the public and the private sector for various purposes.
This `SpaDES` module is a rewrite of the fire component in native R.

## Current Differences

The current version has not been fully tested and compared with the original version, and there are currently several known differences:

- Fire sizes are taken from a Truncated Pareto distribution, resulting in numerous very small fires, and few large fires.
- Parameters have been not been fitted to the landscapes that are under study in the LandWeb project

# Known species

Landmine requires the following codes as inputs (the genus_spec form below), and converts and groups species as follows. Each of the species groups has its own Rate of Spread (ROS) for fire spreading:

* PINU

    * "Pinu_ban"
    * "Pinu_con" 
    * "Pinu_sp" 
* DECI

    * Betu_pap
    * Popu_bal
    * Popu_tre
    * Lari_lar
    
* PICE_MAR

    * Pice_mar
    
* PICE_GLA

    * Pice_gla
    
* ABIE

    * Abie_sp
