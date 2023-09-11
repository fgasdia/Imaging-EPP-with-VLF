function runletkf(parameters)
    @unpack scenario, ens_size, ntimes, dt, pathstep, x_grid, y_grid, modelsteps,
        datatypes, h0, b0, hB, bB, rng, σamp, σphase, data, itp, localization, numexe = parameters()
    @unpack modelproj = common_simulation()

    lengthscale = only(modelsteps).lengthscale

    paths = buildpaths()
    npaths = length(paths)
    
    gridshape = (length(y_grid), length(x_grid))  # useful later
    
    CI = CartesianIndices(gridshape)

    xy_grid = densify(x_grid, y_grid)
    trans = Proj.Transformation(modelproj, wgs84())
    lola = trans.(parent(parent(xy_grid)))

    distarr = lonlatgrid_dists(lola)
    gc = gaspari1999_410(distarr, compactlengthscale(lengthscale))

    # We wrap in Symmetric because we know it is. Without it, eigvals of β can be small
    # (but positive) and still fail the `isposdef` check
    hdistribution = MvNormal(h0, Symmetric(Diagonal(hB)*gc*Diagonal(hB)))  # yes, w/ matrix argument we need variance
    bdistribution = MvNormal(b0, Symmetric(Diagonal(bB)*gc*Diagonal(bB)))

    # Initial ensemble
    h_init = reshape(rand(rng, hdistribution, ens_size), gridshape..., ens_size)
    b_init = reshape(rand(rng, bdistribution, ens_size), gridshape..., ens_size)
    replace!(x->x < MIN_BETA ? MIN_BETA : x, b_init)

    locmask = anylocal(localization)
    h_init[CI[.!locmask],:] .= NaN
    b_init[CI[.!locmask],:] .= NaN

    state = KeyedArray(fill(NaN, 2, length(y_grid), length(x_grid), ens_size, ntimes+1),
            field=[:h, :b], y=y_grid, x=x_grid, ens=1:ens_size, t=0:ntimes)
    state(:h)(t=0) .= h_init
    state(:b)(t=0) .= b_init

    # Generate measurements
    R = [fill(σamp^2, npaths); fill(σphase^2, npaths)]

    # Run model
    ym = KeyedArray(Array{Float64,4}(undef, 2, npaths, ens_size, ntimes+1);
        field=[:amp, :phase], path=pathname.(paths), ens=state.ens, t=0:ntimes)
    H!(x,t) = ensemble_model!(ym(t=t), z->model(itp, z, paths, dt; pathstep, numexe, lwpc=true), x)

    for i in 1:ntimes
        start_time = Dates.now()
        @info "Iteration" i=i start_time
        xnew = LETKF_measupdate(x->H!(x,i-1), state(t=i-1), data(t=i), R;
            localization=localization, datatypes=datatypes)

        # Floor to β = 0.16
        xnew(:b)[xnew(:b) .< MIN_BETA] .= MIN_BETA
        state(t=i) .= xnew

        @info "Elapsed" Δt=canonicalize(Dates.now() - start_time)
        jldsave(joinpath(resdir(scenario), "$scenario.jld2"); state, data, ym)
    end

    # Compute ym with final estimate
    H!(state(t=ntimes), ntimes)

    jldsave(joinpath(resdir(scenario), "$scenario.jld2"); state, data, ym)

    return state, data, ym
end
