using Imaging
using Imaging: buildpaths, segmentediono, GaussianRegion
using Dates, DelimitedFiles, Printf, Statistics
using Parameters, JLD2, Proj4, ProgressMeter
using LongwaveModePropagator
using LMPTools
import SubionosphericVLFInversionAlgorithms as SIA

include("truth_scenarios.jl")

resdir!("results")

function fluxpath(wpts, patch)
    fluxes = Vector{Float64}(undef, length(wpts))
    for i in eachindex(wpts)
        wpt = wpts[i]
        flux = patch(wpt.lon, wpt.lat)
        fluxes[i] = flux
    end
    return fluxes
end

function fluxampphase(tx, rx, background, epp, fluxes, βs; force=false)
    @unpack dt, pathstep = background

    if force || !isfile(resdir("epp_ampphase_normflux.gp.csv"))
        # Normalize flux profile for the plot
        _, wpts = SIA.pathpts(tx, rx; dist=pathstep)
        fluxmat = Matrix{Any}(undef, length(wpts)+1, 2)
        fluxmat[1,:] .= ("dist", "flux")
        fluxmat[2:end,1] .= getindex.(wpts, :dist)/1e3
        fluxvec = fluxpath(wpts, epp.patch)
        fluxmat[2:end,2] .= fluxvec ./ maximum(fluxvec)
        writedlm(resdir("epp_ampphase_normflux.gp.csv"), fluxmat, ',')
    end

    gs = GroundSampler(0:5e3:round(range(tx, rx)+10e3, digits=-4, RoundUp), Fields.Ez)

    mat = Matrix{Any}(undef, 1+length(gs.distance), 1+2*length(fluxes))
    mat[1,1] = "dist"
    mat[2:end,1] .= gs.distance/1e3
    @showprogress for j in eachindex(βs)
        ee = EnergeticElectrons(epp.ee.energy, exp.(-epp.ee.energy/βs[j]), epp.ee.pitchangle,
            epp.ee.pitchangledis)
        
        fullmat = copy(mat)
        if force || !isfile(resdir(@sprintf("epp_ampphase_%s_%d.gp.csv", dt, βs[j]/1e3)))
            for i in eachindex(fluxes)
                patch = GaussianRegion(fluxes[i], epp.patch.x0, epp.patch.y0, epp.patch.σx²,
                    epp.patch.σy², epp.patch.θ, epp.patch.p)

                _, a, p = segmentediono(tx, rx, gs, dt, patch, ee; background.elonly, pathstep)

                aidx = 1 + i
                pidx = aidx + length(fluxes)

                fullmat[1,aidx] = @sprintf("a_%d", log10(fluxes[i]))
                fullmat[1,pidx] = @sprintf("p_%d", log10(fluxes[i]))
                fullmat[2:end,aidx] = a
                fullmat[2:end,pidx] = rad2deg.(mod2pi.(p))
            end
            sname = resdir(@sprintf("epp_ampphase_%s_%d.gp.csv", dt, βs[j]/1e3))
            sname = replace(sname, ":"=>"-")
            writedlm(sname, fullmat, ',')
        end
    end
end

paths = buildpaths()
nml_wh = paths[12]
tx, rx = nml_wh[1], nml_wh[2]
gs = GroundSampler(0:5e3:round(range(tx, rx)+10e3, digits=-4, RoundUp), Fields.Ez)

fluxampphase(tx, rx, day1(), eppb(), (1e3, 1e4, 1e5, 1e6), 100e3:50e3:300e3; force=true)

_, a, p = segmentediono(tx, rx, gs, day1().dt; elonly=false, pathstep=day1().pathstep)
mat = Matrix{Any}(undef, 1+length(gs.distance), 3)
mat[1,:] .= ("dist", "a", "p")
mat[2:end,1] .= gs.distance/1e3
mat[2:end,2] .= a
mat[2:end,3] .= rad2deg.(mod2pi.(p))
sname = resdir("epp_ampphase_"*string(day1().dt)*"_0.gp.csv")
sname = replace(sname, ":"=>"-")
writedlm(sname, mat, ',')


fluxampphase(tx, rx, night1(), eppb(), (1e3, 1e4, 1e5, 1e6), 100e3:50e3:300e3; force=true)

_, a, p = segmentediono(tx, rx, gs, night1().dt; elonly=false, pathstep=night1().pathstep)
mat = Matrix{Any}(undef, 1+length(gs.distance), 3)
mat[1,:] .= ("dist", "a", "p")
mat[2:end,1] .= gs.distance/1e3
mat[2:end,2] .= a
mat[2:end,3] .= rad2deg.(mod2pi.(p))
sname = resdir("epp_ampphase_"*string(night1().dt)*"_0.gp.csv")
sname = replace(sname, ":"=>"-")
writedlm(sname, mat, ',')
