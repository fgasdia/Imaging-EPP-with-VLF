import Imaging
using Dates
using LMPTools

# NOT a public repository. Contact forrest.gasdia@colorado.edu for further information.
# GPI, an EPPIonization dependency, is also NOT a public repository.
using EPPIonization

########
# "Background" scenarios, without EPP

function day1()
    dt = DateTime(2020, 3, 1, 20, 00)  # day
    elonly = false

    pathstep = 50e3  # defined in WGS84!

    σamp, σphase = 0.1, deg2rad(1.0)

    rng = Imaging.reset_rng()

    return (;dt, elonly, pathstep, σamp, σphase, rng)
end

function waitday1()
    dt = DateTime(2020, 3, 1, 20, 00)  # day

    σamp, σphase = 0.1, deg2rad(1.0)
    
    pathstep = 50e3

    rng = Imaging.reset_rng()

    function hbfcn(lon, lat, dt)
        sza = zenithangle(lat, lon, dt)
        coeffh = (2.35, 0.98, -0.17, -0.28, 0.1)
        coeffb = (0.03, 0.01, 0.008, -0.002, -0.008)
        h0, b0 = ferguson(lat, sza, dt)
        hpert = fourierperturbation(sza, coeffh)
        bpert = fourierperturbation(sza, coeffb)
        return h0+hpert, b0+bpert
    end

    return (;dt, hbfcn, pathstep, σamp, σphase, rng)
end

function night1()
    dt = DateTime(2020, 3, 2, 5)  # early night
    elonly = false

    σamp, σphase = 0.1, deg2rad(1.0)
    
    pathstep = 50e3

    rng = Imaging.reset_rng()

    return (;dt, elonly, pathstep, σamp, σphase, rng)
end
night2() = merge(night1(), (;dt=DateTime(2020, 3, 2, 8, 00)))  # West coast midnight
night3() = merge(night1(), (;dt=DateTime(2020, 3, 2, 11, 00)))  # pre dawn

function waitnight1()
    dt = DateTime(2020, 3, 2, 5)  # early night

    σamp, σphase = 0.1, deg2rad(1.0)
    
    pathstep = 50e3

    rng = Imaging.reset_rng()

    function hbfcn(lon, lat, dt)
        sza = zenithangle(lat, lon, dt)
        coeffh = (2.35, 0.98, -0.17, -0.28, 0.1)
        coeffb = (0.03, 0.01, 0.008, -0.002, -0.008)
        h0, b0 = ferguson(lat, sza, dt)
        hpert = fourierperturbation(sza, coeffh)
        bpert = fourierperturbation(sza, coeffb)
        return h0+hpert, b0+bpert
    end

    return (;dt, hbfcn, pathstep, σamp, σphase, rng)
end

########
# EPP scenarios

function eppa()
    # Define patch in wgs84
    patch_lon = -110
    patch_lat = 55
    patch_width, patch_height = 3, 1.5
    peakflux = 1e5
    patch_xy = [patch_lon patch_lat]
    patch_θ = deg2rad(0)
    patch = Imaging.GaussianRegion(peakflux, patch_xy[1], patch_xy[2], patch_width^2, patch_height^2, patch_θ, 4)

    energy = 90e3:1e4:2e6  # 90 keV to 2.2 MeV every 10 keV
    energydis = exp.(-energy/2e5)  # f(E) = exp(-E/β) where β = 200 keV
    pitchangle = 0:90
    pitchangledis = ones(length(pitchangle))
    ee = EnergeticElectrons(energy, energydis, pitchangle, pitchangledis)

    return (;patch, ee)
end

function eppb()
    # Define patch in wgs84
    patch_lon = -120
    patch_lat = 55
    patch_width, patch_height = 11, 1.3
    peakflux = 1e5
    patch_xy = [patch_lon patch_lat]
    patch_θ = deg2rad(1.5)
    patch = Imaging.GaussianRegion(peakflux, patch_xy[1], patch_xy[2], patch_width^2, patch_height^2, patch_θ, 4)

    energy = 90e3:1e4:2e6  # 90 keV to 2.2 MeV every 10 keV
    energydis = exp.(-energy/2e5)  # f(E) = exp(-E/β) where β = 200 keV
    pitchangle = 0:90
    pitchangledis = ones(length(pitchangle))
    ee = EnergeticElectrons(energy, energydis, pitchangle, pitchangledis)

    return (;patch, ee)
end

function eppc()
    # Define patch in wgs84
    patch_lon = -125
    patch_lat = 58
    patch_width, patch_height = 12, 1.5
    peakflux = 1e4
    patch_xy = [patch_lon patch_lat]
    patch_θ = deg2rad(-3)
    patch = Imaging.GaussianRegion(peakflux, patch_xy[1], patch_xy[2], patch_width^2, patch_height^2, patch_θ, 4)

    energy = 90e3:1e4:2e6  # 90 keV to 2.2 MeV every 10 keV
    energydis = exp.(-energy/2e5)  # f(E) = exp(-E/β) where β = 200 keV
    pitchangle = 0:90
    pitchangledis = ones(length(pitchangle))
    ee = EnergeticElectrons(energy, energydis, pitchangle, pitchangledis)

    return (;patch, ee)
end

function eppd()
    # Define patch in wgs84
    patch_lon = -100
    patch_lat = 52
    patch_width, patch_height = 6.0, 0.6
    peakflux = 1e4
    patch_xy = [patch_lon patch_lat]
    patch_θ = deg2rad(0)
    patch = Imaging.GaussianRegion(peakflux, patch_xy[1], patch_xy[2], patch_width^2, patch_height^2, patch_θ, 4)

    energy = 90e3:1e4:2e6  # 90 keV to 2.2 MeV every 10 keV
    energydis = exp.(-energy/2e5)  # f(E) = exp(-E/β) where β = 200 keV
    pitchangle = 0:90
    pitchangledis = ones(length(pitchangle))
    ee = EnergeticElectrons(energy, energydis, pitchangle, pitchangledis)

    return (;patch, ee)
end
