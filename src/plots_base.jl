# Canada Lambert Conformal Conic, for plotting
# epsg_102002() = "+proj=lcc +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
epsg_102002() = "ESRI:102002"

function basemapdir()
    if Sys.iswindows()
        basemapdir = "C:\\gis\\basemaps"
    elseif Sys.islinux()
        basemapdir = "/home/forrest/gis/basemaps"
    end
    return basemapdir
end

"""
    shp2csv(file; crs=nothing, overwrite=false) → newfile

Convert shape files to a CSV called `newfile` using `ogr2ogr`.

If a `crs` is provided, also transform to `crs`. If `overwrite`, then always generate the
new file. Otherwise, if it already exists, just return `newfile`. 

See http://www.gnuplotting.org/code/shape2txt.
"""
function shp2csv(file; crs=nothing, overwrite=false)
    fdir, fname = splitdir(file)
    fname, _ = splitext(fname)
    newfile = joinpath(fdir, fname*".csv")

    overwrite || isfile(newfile) && return newfile

    crsstring = string(crs)

    if Sys.iswindows()
        isfile("C:\\OSGeo4W64\\OSGeo4W.bat") || error("You must install OSGeo4W")

        if isnothing(crs)
            Base.run(`cmd /c C:\\OSGeo4W64\\OSGeo4W.bat ogr2ogr -f CSV -sql "select OGR_GEOMETRY from $fname" $newfile $file -lco GEOMETRY=AS_WKT`)
        else
            Base.run(`cmd /c C:\\OSGeo4W64\\OSGeo4W.bat ogr2ogr -f CSV -sql "select OGR_GEOMETRY from $fname" -t_srs $crsstring $newfile $file -lco GEOMETRY=AS_WKT`)
        end
    else
        if isnothing(crs)
            Base.run(`ogr2ogr -f CSV -sql "select OGR_GEOMETRY from $fname" $newfile $file -lco GEOMETRY=AS_WKT`)
        else
            Base.run(`ogr2ogr -f CSV -sql "select OGR_GEOMETRY from $fname" -t_srs $crsstring $newfile $file -lco GEOMETRY=AS_WKT`)
        end
    end

    return newfile
end

"""
    gpformat(infile, outfile)

Reformat a GIS `infile` to be read by Gnuplot as `outfile`.

See http://www.gnuplotting.org/code/shape2txt.
"""
function gpformat(infile, outfile)
    open(infile, "r") do f
        # Ignore first line
        readline(f)

        open(outfile, "w") do fo
            while !eof(f)
                line = readline(f, keep=true)  # keep newline in string

                # Remove '"LINESTRING (' or '"MULTILINESTRING ((' or '"MULTIPOLYGON (((' at beginning of every line
                line = replace(line, r"^\"LINESTRING \(" => "")
                line = replace(line, r"^\"MULTILINESTRING \(\(" => "")
                line = replace(line, r"^\"POLYGON \(\(" => "")
                line = replace(line, r"^\"MULTIPOLYGON \(\(\(" => "")

                # Remove '"),0,Coastline,d.d' at the end of every line (or any other word)
                # (\r\n|\r|\n) matches `\r` or "\r\n" which we're keeping at the end of the line
                # but then we have to make sure to put a newline back
                line = replace(line, r"\)\"\,[0-9]+\,[A-Z][a-z]*\,[0-9]\.[0-9](\r\n|\r|\n)$" => "\n")

                # Replace '),(' with '\n'
                line = replace(line, r"\)\,\(" => "\n")

                # Double space the file
                line = replace(line, r"\n" => "\n\n")

                # Replace ',' with '\n' for gnuplot
                line = replace(line, r"\," => "\n")

                # Remove remaining POLYGON or MULTIPOLYGON
                line = replace(line, "MULTIPOLYGON" => "")
                line = replace(line, "POLYGON" => "")

                # Remove remaining (, ), or "
                line = replace(line, r"(\)*\")?(\)*)?(\")?(\(*)?" => "")

                # Not in original script, but instead of space separated columns, I prefer
                # commas
                line = replace(line, r"(-?[0-9]+\.[0-9]+) (-?[0-9]+\.[0-9]+)" => s"\g<1>, \g<2>")

                write(fo, line)
            end
        end
    end
end

"""
    gnuplotformat(filename; crs=nothing, overwrite=false) → gpfile

Read the shapefile `filename`, transforming to `crs` if not `nothing`, and then process it
to be read by Gnuplot.

If `overwrite` is false, then only generate `gpfile` if it doesn't already exist.

See also: [`shp2csv`](@ref), [`gpformat`](@ref)
"""
function gnuplotformat(filename; crs=nothing, overwrite=false)
    csvfile = shp2csv(filename; crs, overwrite)

    fdir, _ = splitdir(filename)
    fname, _ = splitext(filename)
    gpfile = fname*".gp.csv"

    overwrite || isfile(gpfile) && return gpfile

    @assert fname*".csv" == joinpath(fdir, csvfile)
    gpformat(fname*".csv", gpfile)

    return gpfile
end

"""
    buildbasemapfiles(scenario, basemapdir, mapproj, paths, mapproj)

Construct:

    - ne_50m_coastline
    - states ne_50m_admin_1_states_provinces_lakes provinces
    - ne_50m_graticules_10

located in `basemapdir` for Gnuplot.

Construct:

    - transmitters
    - receivers
    - great circle paths

saved to the `scenario` directory for Gnuplot using `mapproj` projection.
"""
function buildbasemapfiles(scenario, basemapdir, paths, mapproj; overwrite=false)
    files = (
        "ne_50m_coastline",
        "ne_50m_admin_1_states_provinces_lakes",
        "ne_50m_graticules_10"
    )

    for file in files
        filename = joinpath(basemapdir, file, file*".shp")
        gnuplotformat(filename; crs=mapproj, overwrite)
    end

    transmitters = unique(getindex.(paths, 1))
    receivers = unique(getindex.(paths, 2))

    maptrans = Proj.Transformation(wgs84(), mapproj)
    t_transmitters = maptrans.([(tx.longitude, tx.latitude) for tx in transmitters])
    t_receivers = maptrans.([(rx.longitude, rx.latitude) for rx in receivers])

    t_gcp = Vector{Vector{Tuple{Float64, Float64}}}(undef, length(paths))
    for i in eachindex(paths)
        _, gcp_wpts = SIA.pathpts(paths[i][1], paths[i][2]; dist=25e3)
        lola = [(w.lon, w.lat) for w in gcp_wpts]
        t_gcp_wpts = maptrans.(lola)
        t_gcp[i] = t_gcp_wpts
    end

    writedlm(joinpath(resdir(scenario), "transmitters.gp.csv"), t_transmitters, ',')
    writedlm(joinpath(resdir(scenario), "receivers.gp.csv"), t_receivers, ',')

    rm(joinpath(resdir(scenario), "gcps.gp.csv"); force=true)
    open(joinpath(resdir(scenario), "gcps.gp.csv"), "a") do io
        for i in eachindex(t_gcp)
            writedlm(io, t_gcp[i], ',')
            print(io, "\n")
        end
    end
end

"""
    buildmapfile(itp, values, mapproj, mapx, mapy)

Compute dense map of `values` on a dense grid of `mapx`, `mapy` defined on the projection
`mapproj` using the interpolator `itp`.
"""
function buildmapfile(itp::ScatteredInterpolant, values, mapproj, mapx, mapy)
    vitp = ScatteredInterpolation.interpolate(itp.method, itp.coords, filter(!isnan, values))

    mapxy = densify(mapx, mapy)  # n × 2
    maptrans = Proj.Transformation(mapproj, itp.projection)
    pts = PointSet(maptrans.(parent(parent(mapxy))))

    vmap = Matrix{Float64}(undef, length(mapy), length(mapx))
    for i in eachindex(pts)
        vmap[i] = only(ScatteredInterpolation.evaluate(vitp, pts[i].coords))
    end

    return vmap
end

function buildmapfile(itp::GeoStatsInterpolant, values, mapproj, mapx, mapy)
    geox = georef((v=vec(values),), PointSet(itp.coords))

    mapxy = densify(mapx, mapy)  # n × 2
    maptrans = Proj.Transformation(mapproj, itp.projection)
    pts = maptrans.(parent(parent(mapxy)))

    problem = EstimationProblem(geox, PointSet(pts), :v)
    solution = solve(problem, itp.method)

    vmap = Matrix{Float64}(undef, length(mapy), length(mapx))
    for i in eachindex(pts)
        vmap[i] = solution.v[i]
    end

    return vmap
end

function buildmapfiles(hbfcn, datetime, mapproj, mapx, mapy)
    mapxy = densify(mapx, mapy)  # n × 2
    maptrans = Proj.Transformation(mapproj, wgs84())
    lola = maptrans.(parent(parent(mapxy)))

    hmap = Matrix{Float64}(undef, length(mapy), length(mapx))
    bmap = similar(hmap)

    for i in eachindex(lola)
        h, b = hbfcn(lola[i][1], lola[i][2], datetime)
        hmap[i] = h
        bmap[i] = b
    end

    return hmap, bmap
end

function buildtruthmaps(parameters)
    @unpack modelproj = common_simulation()
    params = parameters()
    @unpack scenario, dt, hbfcn, itp, x_grid, y_grid, modelsteps, pathstep = params

    if :patch in keys(params)
        @unpack patch = params
        epp = true
    else
        epp = false
    end

    paths = buildpaths()
    mapproj = epsg_102002()

    dr, _ = modelsteps[1]

    # Get map boundaries based off of `state`
    xy_grid = densify(x_grid, y_grid)
    (xmin, xmax), (ymin, ymax) = extrema(xy_grid; dims=2)
    mapx, mapy = build_xygrid(xmin, xmax, ymin, ymax, modelproj, mapproj; dr=20e3)
    mapxy_grid = permutedims(densify(mapx, mapy))

    buildbasemapfiles(resdir(scenario), basemapdir(), paths, mapproj; overwrite=false)

    if epp
        function totalhb(x, y, dt)
            h, b = hbfcn(x, y, dt)
            Δh, Δb = patch(x, y)
            return h+Δh, b+Δb
        end
    else
        totalhb = hbfcn
    end
    htrue, btrue = buildmapfiles(totalhb, dt, mapproj, mapx, mapy)

    ctrlpts = [(v[1], v[2]) for v in eachcol(itp.coords)]

    # Convert WGS84 range of 1200km to equivalent grid distance 
    trans = Proj.Transformation(modelproj, wgs84())
    lonlat = trans.(ctrlpts)
    truedr = mediandr(lonlat)
    modelprojscalar = dr/truedr

    # Variogram mask
    # According to https://www.goldensoftware.com/variogramTutorial.pdf
    # the variogram range is often close to the average size of physical anomalies
    varmap = krigingmask(itp, paths, mapproj, mapx, mapy; pathstep, range=600e3*modelprojscalar)
    varmask = replace(x->x > 0.2^2 ? NaN : 1.0, varmap)

    htrue .*= varmask
    btrue .*= varmask

    t = 99
    mat = Matrix{Any}(undef, length(htrue)+1, 4)
    mat[1,:] .= ("x", "y", "htrue", "btrue")
    mat[2:end,:] .= [mapxy_grid[:,1] mapxy_grid[:,2] vec(htrue) vec(btrue)]
    writedlm(joinpath(resdir(scenario), "maps_$t.gp.csv"), mat, ',')

    maptrans = Proj.Transformation(modelproj, mapproj)
    t_ctrlpts = maptrans.(ctrlpts)
    writedlm(joinpath(resdir(scenario), "ctrlpts_$t.gp.csv"), t_ctrlpts, ',')
end

function buildeppmaps(parameters; modelstepindex=1)
    @unpack scenario, patch, modelsteps = parameters()
    @unpack modelproj, west, east, south, north = common_simulation()

    paths = buildpaths()
    mapproj = epsg_102002()

    # Get map boundaries
    dr, lengthscale = modelsteps[modelstepindex]
    x_grid, y_grid = build_xygrid(west, east, south, north, wgs84(), modelproj; dr)
    xy_grid = densify(x_grid, y_grid)
    (xmin, xmax), (ymin, ymax) = extrema(xy_grid; dims=2)
    mapx, mapy = build_xygrid(xmin, xmax, ymin, ymax, modelproj, mapproj; dr=20e3)
    mapxy_grid = densify(mapx, mapy)

    # densexy_grid = permutedims(transform(mapproj, modelproj, mapxy_grid))
    maptrans = Proj.Transformation(mapproj, wgs84())
    densexy_grid = maptrans.(mapxy_grid)

    buildbasemapfiles(scenario, basemapdir(), paths, mapproj; overwrite=false)

    # NOTE TEMP: We're taking log, but we can potentially to a special log plot instead
    fluxmap = Matrix{Float64}(undef, length(mapy), length(mapx))
    for i in eachindex(densexy_grid)
        flux = patch(densexy_grid[i]...)
        fluxmap[i] = log10(flux)
    end

    t = 98
    mat = Matrix{Any}(undef, length(fluxmap)+1, 3)
    mat[1,:] .= ("x", "y", "flux")
    mat[2:end,:] .= [mapxy_grid[:,1] mapxy_grid[:,2] vec(fluxmap)]
    writedlm(joinpath(resdir(scenario), "maps_$t.gp.csv"), mat, ',')

    fcntr = errcontour(fluxmap, mapx, mapy, [3])
    writedlm(joinpath(resdir(scenario), "cntr_$t.gp.csv"), fcntr, ',')
end

function buildfluxprofile(parameters)
    @unpack scenario, dt, patch = parameters()

    # tx = TRANSMITTER[:NML]
    # rx = Receiver("Whitehorse", 60.724, -135.043, 0.0, VerticalDipole())
    # _, wpts = SIA.pathpts(tx, rx; dist=pathstep)
    # wpt = wpts[length(wpts)÷2]
    # lat, lon = wpt.lat, wpt.lon

    lat, lon = patch.r.y0, patch.r.x0
    flux = patch(lon, lat)
    
    z, e1, e2, e3, e4, e5 = epp_profile(patch, lat, lon, dt; epp=true)
    z, n1, n2, n3, n4, n5 = epp_profile(patch, lat, lon, dt; epp=false)
    h = z*1000

    energy = 90e3:1e4:2.2e6  # eV; 90 keV to 2.2 MeV every 10 keV
    energydis = exp.(-energy/2e5)  # f(E) ∝ exp(-E/β) where β ranges from 100 to 300 keV
    pitchangle = 0:90
    pitchdis = ones(length(pitchangle))

    # `md` must be defined at kilometer intervals, but the ionization profile can be at
    # finer `z` steps
    zstepped = first(z):last(z)
    neutraltable = neutralprofiles(lat, lon, z, dt)
    p = GPI.Profiles(neutraltable)
    md = EPPIonization.massdensity.((p,), zstepped)  # g/cm³
    S = ionizationprofile(z, energy, energydis, pitchangle, pitchdis, md/1000)*1e6
    S *= flux  # ionization rate in pairs/m³/s

    mat = Matrix{Any}(undef, length(e1)+1, 12)
    mat[1,:] .= ("alt", "flux", "e1", "e2", "e3", "e4", "e5", "n1", "n2", "n3", "n4", "n5")
    mat[2:end,:] .= [z S e1.(h) e2.(h) e3.(h) e4.(h) e5.(h) n1.(h) n2.(h) n3.(h) n4.(h) n5.(h)]
    writedlm(joinpath(scenario, "flux.csv"), mat, ',')
end

function krigingmask(itp, paths, mapproj, mapx, mapy; pathstep=100e3, range=600e3)
    allwpts = Vector{Tuple{Float64,Float64}}()
    for i in eachindex(paths)
        _, wpts = SIA.pathpts(paths[i][1], paths[i][2]; dist=pathstep)
        append!(allwpts, [(w.lon, w.lat) for w in wpts])
    end

    # There are duplicate entries of some points - let's remove them
    # (otherwise GeoStats doesn't work)
    uidx = unique(x->allwpts[x], 1:length(allwpts))

    trans = Proj.Transformation(wgs84(), itp.projection)
    wptpts = PointSet(trans.(allwpts[uidx]))  # in model proj

    # weights (more samples at same point has more weight)
    # f = Vector{Float64}(undef, length(uidx))
    # for i in eachindex(uidx)
    #     f[i] = count(==(allwpts[uidx[i]]), allwpts)
    # end

    # GeoStats problem setup in model proj
    # geox = georef((f=f,), wptpts)
    geox = georef((f=zeros(length(wptpts)),), wptpts)
    solver = Kriging(:f => (variogram=GaussianVariogram(range=range, sill=1.0), degree=0))

    mapxy = densify(mapx, mapy)  # n × 2
    maptrans = Proj.Transformation(mapproj, itp.projection)
    pts = maptrans.(parent(parent(mapxy)))
    problem = EstimationProblem(geox, PointSet(pts), :f)
    solution = solve(problem, solver)
    # XXX This solution is broken - returning all NaN

    varmap = Matrix{Float64}(undef, length(mapy), length(mapx))
    for i in eachindex(varmap, solution.f_variance)
        varmap[i] = solution.f_variance[i]
    end

    return varmap
end

function errcontour(v, mapx, mapy, bnds)
    c = contours(mapy, mapx, v, bnds)
    allxs = Float64[]
    allys = similar(allxs)
    for cl in Contour.levels(c)
        for line in lines(cl)
            ys, xs = Contour.coordinates(line)
            append!(allxs, xs)
            append!(allys, ys)
            push!(allxs, NaN)
            push!(allys, NaN)
        end
    end
    return [allxs allys]
end

function map_labelled(parameters, maptype, times=0:parameters().ntimes)
    @unpack scenario = parameters()
    @unpack ARG3, ARG4, ARG5, ARG6, ARG8 = maptype

    ARG1 = resdir(scenario)*"/"
    for t in times
        ARG2 = string(t)
        ARG7 = joinpath(ARG1, ARG3*"_"*ARG2*".png")
        c = `gnuplot -c $(joinpath("src", "map_labelled.gp")) $ARG1 $ARG2 $ARG3 $ARG4 $ARG5 $ARG6 $ARG7 $ARG8 --verbose`
        Base.run(c)
    end
end

function plotresiduals(params, finaliter)
    @unpack scenario = params()

    basedir = joinpath(resdir(scenario), "")
    outfile = joinpath(resdir(scenario), "resids.png")

    Base.run(`gnuplot -c $(joinpath("src", "resids.gp")) $basedir $outfile $finaliter`)
end
