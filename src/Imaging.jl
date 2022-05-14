module Imaging

using Dates, DelimitedFiles, LinearAlgebra
using StableRNGs, Parameters, Distributions, JLD2
using AxisKeys, ProgressMeter, Interpolations
using Proj4, Contour, GeoStats, ScatteredInterpolation
using LongwaveModePropagator
using LongwaveModePropagator: QE, ME, waitsparameter

# Unregistered packages. Install in this order.
# EPPIonization has a dependency on GPILowerIonosphere, which is not publically available.
using GeographicLib  # using Pkg; Pkg.add("https://github.com/anowacki/GeographicLib.jl")
using EPPIonization  # Pkg.add("https://github.com/fgasdia/EPPIonization.jl")
using PropagationModelPrep  # Pkg.add("https://github.com/fgasdia/PropagationModelPrep")
using LMPTools  # Pkg.add("https://github.com/fgasdia/LMPTools.jl")
using SubionosphericVLFInversionAlgorithms  # Pkg.add("https://github.com/fgasdia/SubionosphericVLFInversionAlgorithms.jl")
import SubionosphericVLFInversionAlgorithms as SIA


# Path at which to write output.
const RESDIR = Ref(normpath(joinpath(@__DIR__, "..", "results")))
resdir() = RESDIR[]
resdir(d) = joinpath(resdir(), d)
resdir!(d) = isdir(d) ? RESDIR[] = d : throw(ArgumentError("`d` must be a directory"))

"LWPC sometimes fails with low β. Clip the minimum β value to forward models at MIN_BETA."
const MIN_BETA = 0.22
"Sets size of simulated data - effectively sets the maximum number of LETKF iterations."
const DATALENGTH = 10

# consistent random numbers across Julia versions
reset_rng() = StableRNG(1234)
reset_rng(seed) = StableRNG(seed)

include("common.jl")

include("plots_base.jl")
include("letkf_plots.jl")

include("letkf.jl")

end  # module