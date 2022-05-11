function buildmaps(state, parameters)
    params = parameters()
    if :hbfcn in keys(params)
        truthfcn = true
        @unpack scenario, ntimes, pathstep, dt, itp, modelsteps, hbfcn = params
    else
        truthfcn = false
        @unpack scenario, ntimes, pathstep, dt, itp, modelsteps = params
    end
    @unpack modelproj = common_simulation()

    dr, lengthscale = only(modelsteps)

    paths = buildpaths()
    mapproj = epsg_102002()

    # Get map boundaries based off of `state`
    xy_grid = densify(state.x, state.y)
    (xmin, xmax), (ymin, ymax) = extrema(xy_grid; dims=2)
    mapx, mapy = build_xygrid(xmin, xmax, ymin, ymax, modelproj, mapproj; dr=20e3)
    mapxy_grid = permutedims(densify(mapx, mapy))

    buildbasemapfiles(scenario, basemapdir(), paths, mapproj; overwrite=false)

    if truthfcn
        if :patch in keys(params)
            @unpack patch = params
            function totalhb(x, y, dt)
                h, b = hbfcn(x, y, dt)
                Δh, Δb = patch(x, y)
                return h+Δh, b+Δb
            end
            htrue, btrue = buildmapfiles(totalhb, dt, mapproj, mapx, mapy)
        else
            htrue, btrue = buildmapfiles(hbfcn, dt, mapproj, mapx, mapy)
        end
    end

    ctrlpts = permutedims(itp.coords)

    # Convert WGS84 range of 1200km to equivalent grid distance 
    lonlat = permutedims(transform(modelproj, wgs84(), ctrlpts))
    truedr = mediandr(lonlat)
    modelprojscalar = dr/truedr

    # Variogram mask
    # According to https://www.goldensoftware.com/variogramTutorial.pdf
    # the variogram range is often close to the average size of physical anomalies
    varmap = krigingmask(itp, paths, mapproj, mapx, mapy; pathstep, range=600e3*modelprojscalar)
    varmask = replace(x->x > 0.2^2 ? NaN : 1.0, varmap)

    for t = 0:ntimes
        hprimes = dropdims(mean(state(:h)(t=t); dims=:ens); dims=:ens)
        betas = dropdims(mean(state(:b)(t=t); dims=:ens); dims=:ens)

        if all(isnan, hprimes)  # ens didn't run for all of ntimes
            break
        end

        hprimestd = dropdims(std(state(:h)(t=t); dims=:ens); dims=:ens)
        betastd = dropdims(std(state(:b)(t=t); dims=:ens); dims=:ens)

        hmap = buildmapfile(itp, hprimes, mapproj, mapx, mapy)
        bmap = buildmapfile(itp, betas, mapproj, mapx, mapy)

        hstdmap = buildmapfile(itp, hprimestd, mapproj, mapx, mapy)
        bstdmap = buildmapfile(itp, betastd, mapproj, mapx, mapy)

        hmap .*= varmask
        bmap .*= varmask
        hstdmap .*= varmask
        bstdmap .*= varmask

        if truthfcn
            htrue .*= varmask
            btrue .*= varmask

            herr = hmap - htrue
            berr = bmap - btrue

            hcntr = errcontour(herr, mapx, mapy, [-0.5, 0.5])
            bcntr = errcontour(berr, mapx, mapy, [-0.02, 0.02])

            writedlm(joinpath(resdir(scenario), "herrcntr_$(t).gp.csv"), hcntr, ',')
            writedlm(joinpath(resdir(scenario), "berrcntr_$(t).gp.csv"), bcntr, ',')

            mat = Matrix{Any}(undef, length(hmap)+1, 11)
            mat[1,:] .= ("x", "y", "h", "b", "hstd", "bstd", "htrue", "btrue", "herr", "berr", "var")
            mat[2:end,:] .= [mapxy_grid[:,1] mapxy_grid[:,2] vec(hmap) vec(bmap) vec(hstdmap) vec(bstdmap) vec(htrue) vec(btrue) vec(herr) vec(berr) vec(varmap)]
            writedlm(joinpath(resdir(scenario), "maps_$(t).gp.csv"), mat, ',')
        else
            mat = Matrix{Any}(undef, length(hmap)+1, 7)
            mat[1,:] .= ("x", "y", "h", "b", "hstd", "bstd", "var")
            mat[2:end,:] .= [mapxy_grid[:,1] mapxy_grid[:,2] vec(hmap) vec(bmap) vec(hstdmap) vec(bstdmap) vec(varmap)]
            writedlm(joinpath(resdir(scenario), "maps_$(t).gp.csv"), mat, ',')
        end

        t_ctrlpts = transform(modelproj, mapproj, ctrlpts)
        writedlm(joinpath(resdir(scenario), "ctrlpts_$(t).gp.csv"), t_ctrlpts, ',')
    end
    return varmap
end

function residuals(state, data, parameters)
    @unpack itp, dt, pathstep, scenario = parameters()
    paths = buildpaths()

    amps = Matrix{Float64}(undef, length(paths), length(state.t))
    phases = similar(amps)
    i = 1
    for t in state.t
        if all(isnan, state(t=t)(:h))  # ens didn't run for all of ntimes
            amps[:,i] .= NaN
            phases[:,i] .= NaN
        else
            a, p = model(itp, dropdims(mean(state(t=t); dims=:ens); dims=:ens), paths, dt; pathstep, lwpc=true)
            amps[:,i] = a
            phases[:,i] = p
        end

        i += 1
    end
    # Post-fit residuals except for the prior
    ampresid = [amps[:,1].-data(:amp)(t=1) amps[:,2:end].-data(:amp)(t=1:maximum(state.t))]
    phaseresid = [rad2deg.(phasediff.(phases[:,1], data(:phase)(t=1))) rad2deg.(phasediff.(phases[:,2:end], data(:phase)(t=1:maximum(state.t))))]

    writedlm(joinpath(resdir(scenario), "amp_resid.csv"), ampresid, ',')
    writedlm(joinpath(resdir(scenario), "phase_resid.csv"), phaseresid, ',')

    return ampresid, phaseresid
end


"""
    ensembleerror(state, xycoords, parameters)

Compute `state` error over all times at vector of tuples `xycoords` in `modelproj` from `parameters`.
"""
function ensembleerror(state, xycoords, parameters)
    @unpack scenario, modelproj, truthfcn, dt = parameters()

    xh = -20:0.01:20
    xb = -0.5:0.001:0.5

    open(joinpath(scenario, "hens_gaus.csv"), "w") do g
    open(joinpath(scenario, "hens.csv"), "w") do f
        for (xc, yc) in xycoords
            xy_grid = collect(densify(xc, yc))
            lola = permutedims(transform(modelproj, wgs84(), permutedims(xy_grid)))
            @info lola
            
            trueh = truthfcn(lola[1], lola[2], dt)[1]
            xerr = state(y=yc,x=xc) .- trueh
            m = mean(state(y=yc,x=xc)(:h); dims=:ens) .- trueh
            s = std(state(y=yc,x=xc)(:h); dims=:ens)

            tmp = Matrix{Float64}(undef, length(xh), length(state.t)+1)
            tmp[:,1] .= collect(xh)
            for i in eachindex(state.t)
                t = state.t[i]
                println(f, t, ",", join(m(t=t), ","), join(s(t=t), ","), join(xerr(:h)(t=t), ","))
                @. tmp[:,i+1] = (1/(only(s(t=t))*sqrt(2π)))*exp(-(xh-only(m(t=t)))^2/(2*only(s(t=t))^2))
            end
            writedlm(g, tmp, ',')
            println(f, "\n")  # double blank line
            println(g, "\n")
        end
    end
    end

    open(joinpath(scenario, "bens_gaus.csv"), "w") do g
    open(joinpath(scenario, "bens.csv"), "w") do f
        for (xc, yc) in xycoords
            xy_grid = collect(densify(xc, yc))
            lola = permutedims(transform(modelproj, wgs84(), permutedims(xy_grid)))
            
            trueb = truthfcn(lola[1], lola[2], dt)[2]
            xerr = state(y=yc,x=xc) .- trueb
            m = mean(state(y=yc,x=xc)(:b); dims=:ens) .- trueb
            s = std(state(y=yc,x=xc)(:b); dims=:ens)

            tmp = Matrix{Float64}(undef, length(xb), length(state.t)+1)
            tmp[:,1] .= collect(xb)
            for i in eachindex(state.t)
                t = state.t[i]
                println(f, t, ",", join(m(t=t), ","), join(s(t=t), ","), join(xerr(:b)(t=t), ","))
                @. tmp[:,i+1] = (1/(only(s(t=t))*sqrt(2π)))*exp(-(xb-only(m(t=t)))^2/(2*only(s(t=t))^2))
            end
            writedlm(g, tmp, ',')
            println(f, "\n")  # double blank line
            println(g, "\n")
        end
    end
    end

    mapproj = epsg_102002()
    t_xycoords = transform(modelproj, mapproj, [getindex.(xycoords,1) getindex.(xycoords,2)])
    mat = Matrix{Any}(undef, length(xycoords)+1, 3)
    mat[1,:] .= ("i", "x", "y")
    mat[2:end,1] .= collect(1:length(xycoords))
    mat[2:end,2:3] .= t_xycoords
    writedlm(joinpath(scenario, "ens_xycoords.csv"), mat, ',')
end


"""
    ensembleestimate(state, xycoords, parameters)

Fit gaussian to `state` over all times at vector of tuples `xycoords` in `modelproj` from `parameters`.
"""
function ensembleestimate(state, xycoords, parameters)
    @unpack scenario, dt = parameters()
    @unpack modelproj = common_simulation()

    xh = 50:0.01:100
    xb = 0.1:0.001:2.0

    open(joinpath(resdir(scenario), "hens_gaus.csv"), "w") do g
    open(joinpath(resdir(scenario), "hens.csv"), "w") do f
        for (xc, yc) in xycoords
            xy_grid = collect(densify(xc, yc))
            lola = permutedims(transform(modelproj, wgs84(), permutedims(xy_grid)))
            @info lola
            
            x = state(y=yc,x=xc)
            m = mean(x(:h); dims=:ens)
            s = std(x(:h); dims=:ens)

            tmp = Matrix{Float64}(undef, length(xh), length(state.t)+1)
            tmp[:,1] .= collect(xh)
            for i in eachindex(state.t)
                t = state.t[i]
                println(f, t, ",", join(m(t=t), ","), join(s(t=t), ","), join(x(:h)(t=t), ","))
                @. tmp[:,i+1] = (1/(only(s(t=t))*sqrt(2π)))*exp(-(xh-only(m(t=t)))^2/(2*only(s(t=t))^2))
            end
            writedlm(g, tmp, ',')
            println(f, "\n")  # double blank line
            println(g, "\n")
        end
    end
    end

    open(joinpath(resdir(scenario), "bens_gaus.csv"), "w") do g
    open(joinpath(resdir(scenario), "bens.csv"), "w") do f
        for (xc, yc) in xycoords           
            x = state(y=yc,x=xc)
            m = mean(x(:b); dims=:ens)
            s = std(x(:b); dims=:ens)

            tmp = Matrix{Float64}(undef, length(xb), length(state.t)+1)
            tmp[:,1] .= collect(xb)
            for i in eachindex(state.t)
                t = state.t[i]
                println(f, t, ",", join(m(t=t), ","), join(s(t=t), ","), join(x(:b)(t=t), ","))
                @. tmp[:,i+1] = (1/(only(s(t=t))*sqrt(2π)))*exp(-(xb-only(m(t=t)))^2/(2*only(s(t=t))^2))
            end
            writedlm(g, tmp, ',')
            println(f, "\n")  # double blank line
            println(g, "\n")
        end
    end
    end

    mapproj = epsg_102002()
    t_xycoords = transform(modelproj, mapproj, [getindex.(xycoords,1) getindex.(xycoords,2)])
    mat = Matrix{Any}(undef, length(xycoords)+1, 3)
    mat[1,:] .= ("i", "x", "y")
    mat[2:end,1] .= collect(1:length(xycoords))
    mat[2:end,2:3] .= t_xycoords
    writedlm(joinpath(resdir(scenario), "ens_xycoords.csv"), mat, ',')
end

function plotensemble(params, htype, btype)
    @unpack scenario = params()
    
    hlow, hhi = htype.ARG5, htype.ARG6
    blow, bhi = btype.ARG5, btype.ARG6

    basedir = joinpath(resdir(scenario),"")
    outfile = joinpath(resdir(scenario), "ens.png")

    Base.run(`gnuplot -c $(joinpath("src", "letkf_ens.gp")) $basedir $outfile $hlow $hhi $blow $bhi --verbose`)
end

function get_Wr_height(hp, b; zs=0:110e3)
    bfield = BField(50e-6, π/2, 0.0)  # not used
    freq = Frequency(25e3)  # approx
    species = Species(QE, ME, z->waitprofile(z, hp, b; cutoff_low=40e3), electroncollisionfrequency)

    wrs = waitsparameter.(zs, (freq,), (bfield,), (species,))
    idx = argmin(abs.(wrs .- freq.ω))
    return zs[idx]
end


"""
    ionoprofile(lat, lon, dt, patch=nothing, ee=nothing; z=0:110) → (n1, n2, n3, n4, n5)
    
Return electron density and neutral species density profiles as a function of `z` in km at
`lat`, `lon`, `dt`, including EPP if `ee` and `patch` is not nothing.

Usually used for plotting ionosphere profile.
"""
function ionoprofile(lat, lon, dt, patch=nothing, ee=nothing; z=0:110)
    wdcpath = joinpath(@__DIR__, "..", "wdc")

    if !isnothing(patch) && !isnothing(ee)
        flux = patch(lon, lat)
        background, perturbed = chargeprofiles(flux, lat, lon, ee, z, dt; datafilepath=wdcpath)
        n1 = Interpolations.interpolate(z, perturbed[:,1], FritschButlandMonotonicInterpolation())
        n2 = Interpolations.interpolate(z, perturbed[:,2], FritschButlandMonotonicInterpolation())
        n3 = Interpolations.interpolate(z, perturbed[:,3], FritschButlandMonotonicInterpolation())
        n4 = Interpolations.interpolate(z, perturbed[:,4], FritschButlandMonotonicInterpolation())
        n5 = Interpolations.interpolate(z, perturbed[:,5], FritschButlandMonotonicInterpolation())
    else
        background = chargeprofiles(lat, lon, z, dt; datafilepath=wdcpath)
        n1 = Interpolations.interpolate(z, background[:,1], FritschButlandMonotonicInterpolation())
        n2 = Interpolations.interpolate(z, background[:,2], FritschButlandMonotonicInterpolation())
        n3 = Interpolations.interpolate(z, background[:,3], FritschButlandMonotonicInterpolation())
        n4 = Interpolations.interpolate(z, background[:,4], FritschButlandMonotonicInterpolation())
        n5 = Interpolations.interpolate(z, background[:,5], FritschButlandMonotonicInterpolation())
    end

    return n1, n2, n3, n4, n5
end

function buildionoprofile(state, xycoords, parameters)
    alt = 0:0.1:110
    zs = alt*1e3

    waittruth = false
    params = parameters()
    if :hbfcn in keys(params)
        waittruth = true
        @unpack hbfcn, dt, scenario = params
        if :patch in keys(params)
            function totalhb(x, y, dt)
                h, b = hbfcn(x, y, dt)
                Δh, Δb = params.patch(x, y)
                return h+Δh, b+Δb
            end
        else
            totalhb = hbfcn
        end
    elseif :ee in keys(params) && :patch in keys(params)
        @unpack ee, patch, dt, scenario = params
    else
        @unpack dt, scenario = params
        ee, patch = nothing, nothing
    end

    # Find last valid state t estimate
    tmax = maximum(state.t)
    while all(isnan, state(:h)(t=tmax))
        tmax -= 1
    end

    # ionoprofile needs lat/lon (not y,x)
    lola = permutedims(transform(common_simulation().modelproj, wgs84(),
        permutedims(reinterpret(reshape, Float64, xycoords))))

    hbmat = Matrix{Any}(undef, 1+length(xycoords), 2)
    hbmat[1,:] .= ("h", "b")

    mat = Matrix{Any}(undef, length(alt)+1, 1+3*length(xycoords))
    mat[1,1] = "alt"
    mat[2:end,1] .= alt
    for (i, (xc, yc)) in enumerate(xycoords)
        x = state(y=yc, x=xc, t=tmax)
        h = mean(x(:h); dims=:ens)
        b = mean(x(:b); dims=:ens)

        wrcol = 2 + length(xycoords)*(i - 1)
        truecol = wrcol + 1
        waitcol = truecol + 1
        mat[1,wrcol] = "Wr$i"
        mat[1,truecol] = "true$i"
        mat[1,waitcol] = "wait$i"
        mat[2:end,wrcol] .= get_Wr_height(only(h), only(b); zs)/1000 # convert to km

        if waittruth
            truthprofile = waitprofile.(zs, totalhb(lola[1,i], lola[2,i], dt)...)
        else
            n1, n2, n3, n4, n5 = ionoprofile(lola[2,i], lola[1,i], dt, patch, ee; z=alt)
            truthprofile = n1.(alt)
        end

        hbmat[i+1,:] .= (only(h), only(b))
        
        mat[2:end,truecol] .= truthprofile
        mat[2:end,waitcol] .= replace(waitprofile.(zs, h, b; cutoff_low=40e3), 0=>NaN)
    end

    writedlm(joinpath(resdir(scenario), "ionoprofiles_hb.csv"), hbmat, ',')
    writedlm(joinpath(resdir(scenario), "ionoprofiles.csv"), mat, ',')
end

function plotprofiles(params)
    @unpack scenario = params()

    basedir = joinpath(resdir(scenario),"")
    outfile = joinpath(resdir(scenario), "profiles.png")

    Base.run(`gnuplot -c $(joinpath("src", "ionoprofiles.gp")) $basedir $outfile --verbose`)
end
