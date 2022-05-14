# Imaging-EPP-with-VLF

Code used for the paper "A Method for Imaging Energetic Particle Precipitation with Subionospheric VLF Signals" (in review).

This repository contains source code and scripts in the [Julia](https://julialang.org/) programming language. Several registered, unregistered-but-public, and private repositories are dependencies for running the scripts and installing this package. Ideally the entire code base would be under MIT license, but I have not obtained permission to make GPILowerIonosphere.jl public. Contact me at Forrest.Gasdia@colorado.edu for additional information. Furthermore, this code calls a custom parallelized deployment of [LWPC](http://www.dtic.mil/docs/citations/ADA350375).

Although this code cannot be installed and run as-is by the general user, I hope the source code itself indicates how the LETKF algorithm can be applied to measurements made by an array of VLF receivers to image EPP. Output files with data used to make the plots in the paper are available at XXXXXXXXX. The base package is Imaging (src/Imaging.jl) and the scripts are scripts/epp_ampphase.jl and scripts/letkf_scenarios.jl. 
