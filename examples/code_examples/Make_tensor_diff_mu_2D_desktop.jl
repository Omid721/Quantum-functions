using ITensors
import TensorCrossInterpolation as TCI
using TCIITensorConversion
using QuanticsGrids
using Plots
using HDF5
include("functions.jl")



Ln = 20 # for n
Lpx = 10 # for px
Lpy = 10 # for py
Lkx = 8 # for kx
Lky = 8 # for ky
Lμ = 3 # for μ


###########################
#  kx  px  n   py  ky  Lμ
#  O---O---O---O---O---O
#  |   |   |   |   |   |
###########################

L = Int(Ln+ Lpx + Lpy + Lkx + Lky + Lμ)

β = 100
p_max = 2π
k_max = 2π
ω_max = 0
μ_max = 3
#kx = k_max * cos(π/3)
#ky = k_max * sin(π/3)
T = 1/β
m = 0 #The


epsilon(kx:: Real, ky::Real, μ:: Real) = sqrt(abs(3 + 2 * cos(sqrt(3) * ky) + 4 * cos(sqrt(3) / 2 * ky) * cos(3 / 2 * kx))) - μ #sqrt(px^2+ py^2) - μ

G(iω, kx, ky, μ) = 1/(iω - epsilon(kx, ky, μ))



function f(v)
    # Extract sections of the vector `v` in the order: kx, px, n, py, ky

    x3 = b2c(v[1:Lkx])  # kx is the first segment
    x1 = b2c(v[Lkx + 1:Lkx + Lpx])  # px follows kx
    x = b2c(v[Lkx + Lpx + 1:Lkx + Lpx + Ln])  # n follows px
    x2 = b2c(v[Lkx + Lpx + Ln + 1:Lkx + Lpx + Ln + Lpy])  # py follows n
    x4 = b2c(v[Lkx + Lpx + Ln + Lpy + 1:Lkx + Lpx + Ln + Lpy + Lky])  # ky follows py
    x5 = b2c(v[Lkx + Lpx + Ln + Lpy + 1:Lkx + Lpx + Ln + Lpy + Lky + Lμ]) #μ follows ky

    # Compute intermediate values
    n = 2^(Ln + 1) * (x - 1/2)
    iω = (2*n + 1) * π * T * im
    px = p_max * (x1 - 1/2)
    py = p_max * (x2 - 1/2)
    kx = k_max * (x3 - 1/2)
    ky = k_max * (x4 - 1/2)
    μ = μ_max * x5

    # Return the Green's function G
    return G(iω, kx + px, ky + py, μ)
end

localdims = fill(2, L)  # Fill the tensor with dimensions of size 2

tolerance = 1e-8

tci, ranks, errors = TCI.crossinterpolate2(ComplexF64, f, localdims; tolerance=tolerance, maxbonddim = 1500)

mps_shifted = MPS(tci)

println("mps_shifted calculation is done")

h5open("mps_shifted_kx8_px10_n20_py10_ky8_mu.h5", "w") do file
    write(file, "MPS", mps_shifted)
end
println("file saved successfully")


# Save kx_range, ky_range, and Gk2_shift to an HDF5 file


h5open("MPS_shift_results.h5", "w") do file
    # Convert LinRange to array using collect() before saving
    write(file, "mps_shifted", mps_shifted)

    # Save the details of the network as attributes or in a dataset named "details"
    details = Dict(
        "Ln" => Ln,
        "Lpx" => Lpx,
        "Lpy" => Lpy,
        "Lkx" => Lkx,
        "Lky" => Lky,
        "Lμ" => Lμ,
        "β" => β,
        "p_max" => p_max,
        "k_max" => k_max,
        "μ_max" => μ_max,
        "network_structure" => """
            ###########################
            #  Lμ  kx  px  n   py  ky
            #  O---O---O---O---O---O
            #  |   |   |   |   |   |
            ###########################
        """
    )

    # Write the details dictionary to a group named "details"
    grp = create_group(file, "details")
    for (key, value) in details
        write(grp, key, value)
    end
end

println("Done!")
