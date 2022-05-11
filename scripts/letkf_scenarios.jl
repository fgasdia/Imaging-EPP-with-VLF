import Imaging
using Imaging: common_simulation, buildpaths,
    buildmaps, buildeppmaps, map_labelled, ensembleestimate, plotensemble,
    buildionoprofile, plotprofiles, residuals, plotresiduals, resdir, resdir!
using Dates, Printf
using Parameters, JLD2, ScatteredInterpolation, Proj4
using SubionosphericVLFInversionAlgorithms, LMPTools

include("truth_scenarios.jl")
resdir!("results")


"""
    letkf1(background, prior=nothing)

Return a NamedTuple LETKF parameter set:
    - ens_size
    - ntimes (number of LETKF iterations)
    - numexe (number of parallel LWPC instances)
    - modelsteps (grid resolution and spatial covariance)
    - itp (interpolation object)
    - hB, bB (prior standard deviation of h′ and β)
    - h0, b0 (prior mean of h′ and β)
    - datatypes (amplitude and phase)
    - x_grid, y_grid (projected grid on which h′ and β are defined)
    - localization (localization object between grid and propagation paths)
"""
function letkf1(background, prior=nothing)
    @unpack west, east, south, north, modelproj = common_simulation()
    @unpack dt = background()

    ens_size = 100
    ntimes = 6

    numexe = 8
    datatypes = (:amp, :phase)

    # Setup grid
    dr = 300e3
    lengthscale = 600e3
    modelsteps = ((;dr, lengthscale),)

    x_grid, y_grid = build_xygrid(west, east, south, north, wgs84(), modelproj; dr)
    xy_grid = collect(densify(x_grid, y_grid))
    lola = permutedims(transform(modelproj, wgs84(), permutedims(xy_grid)))

    paths = buildpaths()
    localization = obs2grid_distance(lola, paths; r=lengthscale)
    filterbounds!(localization, lola, west, east, south, north)

    # Build interpolant (needed for plots)
    itppts = build_xygrid(anylocal(localization), x_grid, y_grid)
    itp = ScatteredInterpolant(GeneralizedPolyharmonic(1,1), modelproj, itppts)

    ncells = size(lola, 2)
    if isnothing(prior)
        hB = fill(2, ncells)  # σ_h′
        bB = fill(0.04, ncells)  # σ_β

        hb0 = [ferguson(lola[2,i], zenithangle(lola[2,i], lola[1,i], dt), dt) for i in axes(lola,2)]
        h0 = getindex.(hb0, 1)
        b0 = getindex.(hb0, 2)
    else
        @unpack hB, bB, h0, b0 = prior()
    end
    @assert length(h0) == length(hB) == ncells

    return (;ens_size, ntimes, numexe, modelsteps, itp, hB, bB, h0, b0, datatypes,
        x_grid, y_grid, localization)
end

"""
    letkf4(background, prior=nothing)

Like `letkf1`, except uses `b0 = 0.9` unless overridden by `prior`.
"""
function letkf4(background, prior=nothing)
    # Like letkf1 but with high beta prior
    @unpack hB, bB, h0, b0 = letkf1(background, prior)
    
    if isnothing(prior)
        fill!(b0, 0.9)
    else
        @unpack hB, bB, h0, b0 = prior()
    end

    return merge(letkf1(background, prior), (;hB, bB, h0, b0))
end

"""
    runandplotfullday(scenarios, letkffcn, force=false)

For each of `scenarios` (no EPP), apply the `letkffcn` to estimate the ionosphere and
generate data for plotting.

By default, this skips running `letkffcn` if it's already been run. Set `force=true` to
run `letkffcn` regardless.

See also: [runandplotfulldayepp](@ref)
"""
function runandplotfullday(scenarios, letkffcn, force=false)
    htype = (ARG3="h", ARG4="h′ (km)", ARG5=68.0, ARG6=90.0, ARG8="amp.pal")
    btype = (ARG3="b", ARG4="β (km^{-1})", ARG5=0.2, ARG6=0.8, ARG8="tempo.pal")
    hstdtype = (ARG3="hstd", ARG4="std h′ (km)", ARG5=0.0, ARG6=2.0, ARG8="dense.pal")        
    htruetype = (ARG3="htrue", ARG4="h′ (km)", ARG5=68.0, ARG6=90.0, ARG8="amp.pal")
    herrtype = (ARG3="herr", ARG4="err h′ (km)", ARG5=-2.0, ARG6=2.0, ARG8="coolwarm.pal")
    btruetype = (ARG3="btrue", ARG4="β (km^{-1})", ARG5=0.2, ARG6=0.8, ARG8="tempo.pal")
    berrtype = (ARG3="berr", ARG4="err β (km^{-1})", ARG5=-0.1, ARG6=0.1, ARG8="coolwarm.pal")
    fluxtype = (ARG3="flux", ARG4="flux", ARG5=0, ARG6=6, ARG8="devon.pal")

    for f in scenarios
        s = string(f)
        scenario = s*"_"*string(letkffcn)
        isdir(resdir(scenario)) || mkdir(resdir(scenario))
        if :hbfcn in keys(f())
            data = Imaging.waittruth(s, f())
        else
            data = Imaging.truth(s, f())
        end
        params() = merge(f(), letkffcn(f), (;scenario, data))

        # Skip running LETKF if it's already been run
        if force || !isfile(joinpath(resdir(scenario), scenario*".jld2"))
            try
                Imaging.runletkf(params)
            catch e
                error(e)
                nothing
            end
        end
        
        state, data = load(joinpath(resdir(scenario), scenario*".jld2"), "state", "data")
        buildmaps(state, params)
        if :hbfcn in keys(params())
            [map_labelled(params, t) for t in (htype, btype, hstdtype, htruetype, btruetype, herrtype, berrtype)]
        else
            [map_labelled(params, t) for t in (htype, btype, hstdtype)]
            map_labelled(params, fluxtype, 98)
        end

        xycoords = [(state.x[12], state.y[4]), (state.x[9], state.y[6]), (state.x[7], state.y[9])]
        ensembleestimate(state, xycoords, params)
        plotensemble(params, htype, btype)

        buildionoprofile(state, xycoords, params)
        plotprofiles(params)

        try
            residuals(state, data, params)
            plotresiduals(params, 7)
        catch e
            @error e
            nothing
        end
    end
end

"""
    runandplotfulldayepp(scenarios, letkffcn, epp, force=false; hBbB=nothing, b0=nothing)

For each of `scenarios` with `epp`, apply the `letkffcn` to estimate the ionosphere and
generate data for plotting.

Optionally provide `hBbB`, a tuple with the prior standard deviation of h′ and β, and `b0`,
the prior mean of β. Otherwise, the prior used is the estimate generated for the `background`
scenario without EPP.

See also: [runandplotfullday](@ref)
"""
function runandplotfulldayepp(scenarios, letkffcn, epp, force=false; hBbB=nothing, b0=nothing)
    htype = (ARG3="h", ARG4="h′ (km)", ARG5=60.0, ARG6=90.0, ARG8="amp.pal")
    btype = (ARG3="b", ARG4="β (km^{-1})", ARG5=0.2, ARG6=0.8, ARG8="tempo.pal")
    hstdtype = (ARG3="hstd", ARG4="std h′ (km)", ARG5=0.0, ARG6=2.0, ARG8="dense.pal")        
    htruetype = (ARG3="htrue", ARG4="h′ (km)", ARG5=60.0, ARG6=90.0, ARG8="amp.pal")
    herrtype = (ARG3="herr", ARG4="err h′ (km)", ARG5=-2.0, ARG6=2.0, ARG8="coolwarm.pal")
    btruetype = (ARG3="btrue", ARG4="β (km^{-1})", ARG5=0.2, ARG6=0.8, ARG8="tempo.pal")
    berrtype = (ARG3="berr", ARG4="err β (km^{-1})", ARG5=-0.1, ARG6=0.1, ARG8="coolwarm.pal")
    fluxtype = (ARG3="flux", ARG4="flux", ARG5=0, ARG6=6, ARG8="devon_r.pal")

    vartype = (ARG3="var", ARG4="var", ARG5=0, ARG6=1, ARG8="dense.pal")

    for f in scenarios
        s = string(f)
        backgroundscenario = s*"_"*string(letkffcn)
        if isnothing(hBbB)
            scenario = s*"_"*string(epp)*"_"*string(letkffcn)
        else
            if isnothing(b0)
                scenario = s*"_"*string(epp)*"_"*string(letkffcn)*"_"*@sprintf("%.1f-%.2f", hBbB...)
            else
                scenario = s*"_"*string(epp)*"_"*string(letkffcn)*"_"*@sprintf("%.2f", b0)
            end
        end
        isdir(resdir(scenario)) || mkdir(resdir(scenario))

        if :hbfcn in keys(f())
            data = Imaging.waittruth(s*"_"*string(epp), f(), epp())
        else
            data = Imaging.truth(s*"_"*string(epp), f(), epp())
        end        

        # Load prior from background estimate
        state = load(joinpath(resdir(backgroundscenario), backgroundscenario*".jld2"), "state")
        tmax = maximum(state.t)
        while all(isnan, state(:h)(t=tmax))
            tmax -= 1
        end

        hb0 = mean(state(t=tmax), dims=:ens)
        if isnothing(b0)
            h0 = vec(hb0(:h))
            b0 = vec(hb0(:b))
        else
            h0 = vec(hb0(:h))
            b0 = fill(b0, length(hb0(:b)))
        end

        if isnothing(hBbB)
            # Let's try using the same hB, bB as the background...
            @unpack hB, bB = letkffcn(f)
        else
            hB = fill(hBbB[1], length(h0))
            bB = fill(hBbB[2], length(b0))
        end

        prior() = (;h0, b0, hB, bB)
        params() = merge(f(), epp(), letkffcn(f, prior), (;scenario, data))

        # Skip running LETKF if it's already been run
        if force || !isfile(joinpath(resdir(scenario), scenario*".jld2"))
            try
                if isnothing(b0)
                    Imaging.runletkf(params)
                else
                    Imaging.runletkf_honly(params, b0)
                end
            catch e
                error(e)
                nothing
            end
        end

        state, data = load(joinpath(resdir(scenario), scenario*".jld2"),
            "state", "data")
        buildmaps(state, params)

        if :hbfcn in keys(params())
            [map_labelled(params, t) for t in (htype, btype, hstdtype, htruetype, btruetype, herrtype, berrtype)]
        else
            buildeppmaps(params)
            [map_labelled(params, t) for t in (htype, btype, hstdtype)]
            map_labelled(params, fluxtype, 98)

            map_labelled(params, vartype, 0)
        end

        xycoords = [(state.x[12], state.y[4]), (state.x[9], state.y[6]), (state.x[7], state.y[9])]
        ensembleestimate(state, xycoords, params)
        plotensemble(params, htype, btype)

        buildionoprofile(state, xycoords, params)
        plotprofiles(params)

        try
            residuals(state, data, params)
            plotresiduals(params, 7)
        catch e
            @error e
            nothing
        end
    end
end

########
# Generate simulated truth observations and estimate the ionosphere using the following
# scenarios.
#
# NOTE: These functions call LWPC using a custom parallelized deployment and use non-public
# code repositories. These functions cannot be run as-is on general installs.

# No EPP, Wait and Spies and realistic daytime ionospheres
# runandplotfullday((day1, waitday1), letkf1)

# No EPP, realistic nighttime ionosphere
runandplotfullday((night1,), letkf1)

# EPP scenarios with realistic daytime ionosphere
runandplotfulldayepp((day1,), letkf1, eppa)
# runandplotfulldayepp((day1,), letkf1, eppb)
# runandplotfulldayepp((day1,), letkf1, eppc)
# runandplotfulldayepp((day1,), letkf1, eppd)

# Additional realistic nighttime ionospheres, no EPP
# runandplotfullday((night2, night3), letkf1)

# Alternative nighttime prior with β = 0.9 km⁻¹
# runandplotfullday((night1,), letkf4)
