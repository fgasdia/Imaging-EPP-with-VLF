"""
    common_simulation()

Return parameters common to all scenarios.
"""
function common_simulation()
    modelproj = esri_102010()

    # Simulation boundaries
    west, east = -136.5, -91
    south, north = 46, 64

    return (;modelproj, west, east, south, north)
end

"""
    buildpaths()

Return a vector of `(Transmitter, Receiver)` propagation paths used in the scenarios.
"""
function buildpaths()
    transmitters = [TRANSMITTER[:NLK], TRANSMITTER[:NML]]
    receivers = [
        Receiver("Whitehorse", 60.724, -135.043, 0.0, VerticalDipole()),
        Receiver("Churchill", 58.74, -94.085, 0.0, VerticalDipole()),
        Receiver("Stony Rapids", 59.253, -105.834, 0.0, VerticalDipole()),
        Receiver("Fort Smith", 60.006, -111.92, 0.0, VerticalDipole()),
        Receiver("Bella Bella", 52.1675508, -128.1545219, 0.0, VerticalDipole()),
        Receiver("Nahanni Butte", 61.0304412, -123.3926734, 0.0, VerticalDipole()),
        Receiver("Juneau", 58.32, -134.41, 0.0, VerticalDipole()),
        Receiver("Ketchikan", 55.35, -131.673, 0.0, VerticalDipole()),
        Receiver("Winnipeg", 49.8822, -97.1308, 0.0, VerticalDipole()),
        Receiver("IslandLake", 53.8626, -94.6658, 0.0, VerticalDipole()),
        Receiver("Gillam", 56.3477, -94.7093, 0.0, VerticalDipole())
    ]
    paths = [(tx, rx) for tx in transmitters for rx in receivers]

    return paths
end

"Localize grid points to propagation paths using `obs2grid_distance`."
struct Localizer
    lengthscale::Float64
    west::Float64
    east::Float64
    south::Float64
    north::Float64
end
(l::Localizer)(lonlat) = anylocal(filterbounds!(
    obs2grid_distance(lonlat, buildpaths(); r=l.lengthscale), lonlat, l.west, l.east, l.south, l.north
))

function localize!(v, x_grid, y_grid, localizationfcn, modelproj)
    xy_grid = densify(x_grid, y_grid)
    trans = Proj.Transformation(modelproj, wgs84())
    lola = trans.(parent(parent(xy_grid)))
    locmask = localizationfcn(lola)

    v[.!locmask] .= NaN

    return v
end

"GaussianRegion used to represent EPP patches."
struct GaussianRegion
    peak::Float64
    x0::Float64
    y0::Float64
    σx²::Float64
    σy²::Float64
    θ::Float64  # rad, rotation CW
    p::Float64
end
GaussianRegion(peak, x0, y0, σx²,σy², θ) = GaussianRegion(peak, x0, y0, σx², σy², θ, 1)

function (r::GaussianRegion)(x, y)
    x0, y0 = r.x0, r.y0
    θ, σx², σy² = r.θ, r.σx², r.σy²
    p = r.p

    s, c = sincos(θ)
    s², c² = s^2, c^2
    s2θ = sin(2θ)

    a = c²/(2σx²) + s²/(2σy²)
    b = -s2θ/(4σx²) + s2θ/(4σy²)
    c = s²/(2σx²) + c²/(2σy²)

    f = r.peak*exp(-(a*(x - x0)^2 + 2b*(x - x0)*(y - y0) + c*(y - y0)^2)^p)

    return f
end

"""
    waittruth(name, background, epp=nothing; force=false)

Generate simulated truth observations when the truth ionosphere has a Wait and Spies
exponential profile.
"""
function waittruth(name, background, epp=nothing; force=false)
    paths = buildpaths()

    @unpack dt, pathstep, hbfcn, rng, σamp, σphase = background

    obs_fname = resdir(name*".jld2")
    if force || !isfile(obs_fname)
        if isnothing(epp)
            obsamp, obsphase = waitbatchpropagate(paths, dt, hbfcn; pathstep)
        else
            @unpack patch = epp
            obsamp, obsphase = waitbatchpropagate(paths, dt, hbfcn, patch; pathstep)
        end

        jldsave(obs_fname; obsamp, obsphase)
    else
        f = jldopen(obs_fname, "r")
        obsamp, obsphase = f["obsamp"], f["obsphase"]
    end

    npaths = length(paths)
    data = KeyedArray(Array{Float64,3}(undef, 4, npaths, DATALENGTH);
        field=[:amp, :phase, :amp_noiseless, :phase_noiseless], path=pathname.(paths), t=1:DATALENGTH)
    data(:amp_noiseless) .= obsamp
    data(:phase_noiseless) .= obsphase
    data(:amp) .= obsamp .+ σamp.*randn(rng, npaths, DATALENGTH)
    data(:phase) .= obsphase .+ σphase.*randn(rng, npaths, DATALENGTH)

    return data
end

function waitsegmentediono(tx, rx, dt, hbfcn, patch=nothing; pathstep=100e3)
    _, wpts = SIA.pathpts(tx, rx; dist=pathstep)
    geoaz = inverse(tx.longitude, tx.latitude, rx.longitude, rx.latitude).azi

    wvgs = Vector{HomogeneousWaveguide{Species}}(undef, length(wpts))
    for i in eachindex(wpts)
        wpt = wpts[i]
        
        ground = GROUND[LMPTools.get_groundcode(wpt.lat, wpt.lon)]
        bfield = igrf(geoaz, wpt.lat, wpt.lon, year(dt))

        h, b = hbfcn(wpt.lon, wpt.lat, dt)
        if !isnothing(patch)
            Δh, Δb = patch(wpt.lon, wpt.lat)
            h += Δh  # Δh is usually a negative value for EPP
            b += Δb
        end

        species = Species(QE, ME, z->waitprofile(z, h, b), electroncollisionfrequency)
        
        wvgs[i] = HomogeneousWaveguide(bfield, species, ground, wpt.dist)
    end

    wvg = SegmentedWaveguide(wvgs)

    gs = GroundSampler(range(tx, rx), Fields.Ez)
    propagate(wvg, tx, gs)
end

function waitbatchpropagate(paths, dt, hbfcn, patch=nothing; pathstep)
    amps = Vector{Float64}(undef, length(paths))
    phases = similar(amps)
    @showprogress for i in eachindex(paths)
        tx, rx = paths[i]
        
        _, a, p = waitsegmentediono(tx, rx, dt, hbfcn, patch; pathstep)
        
        amps[i] = a
        phases[i] = p
    end

    return amps, phases
end

"""
    truth(name, background, epp=nothing; force=false)

Generate simulated truth observations with realistic ionosphere profiles.

# Example

```julia
data = truth("obs_day1", day1(), epp1())
```
"""
function truth(name, background, epp=nothing; force=false)
    paths = buildpaths()

    @unpack dt, pathstep, elonly, rng, σamp, σphase = background

    obs_fname = resdir(name*".jld2")
    if force || !isfile(obs_fname)
        if isnothing(epp)
            obsamp, obsphase = batchpropagate(paths, dt; elonly, pathstep)
        else
            @unpack patch, ee = epp
            obsamp, obsphase = batchpropagate(paths, dt, patch, ee; elonly, pathstep)
        end

        jldsave(obs_fname; obsamp, obsphase)
    else
        f = jldopen(obs_fname, "r")
        obsamp, obsphase = f["obsamp"], f["obsphase"]
    end

    npaths = length(paths)
    data = KeyedArray(Array{Float64,3}(undef, 4, npaths, DATALENGTH);
        field=[:amp, :phase, :amp_noiseless, :phase_noiseless], path=pathname.(paths), t=1:DATALENGTH)
    data(:amp_noiseless) .= obsamp
    data(:phase_noiseless) .= obsphase
    data(:amp) .= obsamp .+ σamp.*randn(rng, npaths, DATALENGTH)
    data(:phase) .= obsphase .+ σphase.*randn(rng, npaths, DATALENGTH)

    return data
end

function segmentediono(tx, rx, gs, dt::DateTime, patch=nothing, ee=nothing; elonly::Bool, pathstep=100e3)
    z = 0:110
    h = z*1000

    _, wpts = SIA.pathpts(tx, rx; dist=pathstep)

    geoaz = inverse(tx.longitude, tx.latitude, rx.longitude, rx.latitude).azi

    if elonly
        wvgs = Vector{HomogeneousWaveguide{Species}}(undef, length(wpts))
    else
        wvgs = Vector{HomogeneousWaveguide{NTuple{5, Species}}}(undef, length(wpts))
    end

    for i in eachindex(wpts)
        wpt = wpts[i]
        
        ground = GROUND[LMPTools.get_groundcode(wpt.lat, wpt.lon)]
        bfield = igrf(geoaz, wpt.lat, wpt.lon, year(dt))

        if !isnothing(patch) && !isnothing(ee)
            flux = patch(wpt.lon, wpt.lat)
            background, perturbed = chargeprofiles(flux, wpt.lat, wpt.lon, ee, z, dt)
            n1 = Interpolations.interpolate(h, perturbed[:,1], FritschButlandMonotonicInterpolation())
            n2 = Interpolations.interpolate(h, perturbed[:,2], FritschButlandMonotonicInterpolation())
            n3 = Interpolations.interpolate(h, perturbed[:,3], FritschButlandMonotonicInterpolation())
            n4 = Interpolations.interpolate(h, perturbed[:,4], FritschButlandMonotonicInterpolation())
            n5 = Interpolations.interpolate(h, perturbed[:,5], FritschButlandMonotonicInterpolation())
        else
            background = chargeprofiles(wpt.lat, wpt.lon, z, dt)
            n1 = Interpolations.interpolate(h, background[:,1], FritschButlandMonotonicInterpolation())
            n2 = Interpolations.interpolate(h, background[:,2], FritschButlandMonotonicInterpolation())
            n3 = Interpolations.interpolate(h, background[:,3], FritschButlandMonotonicInterpolation())
            n4 = Interpolations.interpolate(h, background[:,4], FritschButlandMonotonicInterpolation())
            n5 = Interpolations.interpolate(h, background[:,5], FritschButlandMonotonicInterpolation())
        end

        if elonly
            species = Species(QE, ME, n1, electroncollisionfrequency)
        else
            species = (
                Species(QE, ME, n1, electroncollisionfrequency),  # Ne
                Species(QE, 58000ME, n2, ioncollisionfrequency),  # N-
                Species(-QE, 120000ME, n3, ioncollisionfrequency),  # Nx+
                Species(QE, 120000ME, n4, ioncollisionfrequency),  # Nx-
                Species(-QE, 58000ME, n5, ioncollisionfrequency),  # N+
            )
        end

        wvgs[i] = HomogeneousWaveguide(bfield, species, ground, wpt.dist)
    end

    wvg = SegmentedWaveguide(wvgs)

    propagate(wvg, tx, gs)
end

function segmentediono(tx, rx, dt::DateTime, patch=nothing, ee=nothing; elonly::Bool, pathstep=100e3)
    gs = GroundSampler(range(tx, rx), Fields.Ez)
    segmentediono(tx, rx, gs, dt, patch, ee; elonly, pathstep)
end

function batchpropagate(paths, dt, patch=nothing, ee=nothing; elonly, pathstep)
    amps = Vector{Float64}(undef, length(paths))
    phases = similar(amps)
    @showprogress for i in eachindex(paths)
        tx, rx = paths[i]
        
        _, a, p = segmentediono(tx, rx, dt, patch, ee; elonly, pathstep)
        
        amps[i] = a
        phases[i] = p
    end

    return amps, phases
end
