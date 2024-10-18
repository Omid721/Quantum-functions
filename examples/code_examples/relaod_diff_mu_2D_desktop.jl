using ITensors
using HDF5

# Open the HDF5 file and load the data
h5open("results.h5", "r") do file
    # Read kx_range and ky_range
    kx_range = read(file, "kx_range")
    ky_range = read(file, "ky_range")

    # Read Gk2_shift
    Gk2_shift = read(file, "Gk2_shift")

    # Read mps_shifted
    mps_shifted = read(file, "mps_shifted", MPS)

    # Read details group
    details_grp = file["details"]
    details = Dict()
    for key in keys(details_grp)
        details[key] = read(details_grp, key)
    end

    
    println("Details:")
    for (key, value) in details
        println(key, ": ", value)
    end
end
