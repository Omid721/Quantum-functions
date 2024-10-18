using ITensors
import TensorCrossInterpolation as TCI
using TCIITensorConversion
using QuanticsGrids
using Plots
using HDF5
include("functions.jl")


mps_shifted = h5open("mps_shifted_kx8_px10_n20_py10_ky8_mu.h5", "r") do file
    read(file, "MPS", MPS)
end



Ln = 20 # for n
Lpx = 10 # for px
Lpy = 10 # for py
Lkx = 8 # for kx
Lky = 8 # for ky
Lμ = 5 # for n


###########################
#  Lμ  kx  px  n   py  ky
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


# Define k_max and the ranges for kx and ky
kx_range = LinRange(-k_max/2, k_max/2, 100)
ky_range = LinRange(-k_max/2, k_max/2, 100)
n =  0
qx = 2
qy = 1
μ = 1
qx_b = c2b(qx/k_max + 1/2, Lkx)
qy_b = c2b(qy/k_max + 1/2, Lky)

μ_b = c2b(μ/μ_max, Lμ)



x = c2b(1/2, Ln)
x1 = c2b(2/p_max + 1/2, Lpx)
x2 = c2b(1/p_max + 1/2, Lpy)
s = [qx_b; x1; x ;x2; qy_b; μ_b ]



# Initialize an empty matrix to store Πk values
Gk2_shift = zeros(ComplexF64, length(kx_range), length(ky_range))



# Loop over kx and ky to evaluate Πk
for i in 1:length(kx_range)
    for j in 1:length(ky_range)
        kx = kx_range[i]
        ky = ky_range[j]

        # Transformation functions
        x = c2b(1/2, Ln)
        x1 = c2b(kx/p_max + 1/2, Lpx)
        x2 = c2b(ky/p_max + 1/2, Lpy)
        s = [μ_b ; qx_b; x1; x ;x2; qy_b;]        # Evaluate Πk for the current kx and ky
        Gk2_shift[i, j] = evaluate(mps_shifted, s)

    end


    println(i/length(kx_range))

end


# Save kx_range, ky_range, and Gk2_shift to an HDF5 file


h5open("results.h5", "w") do file
    # Convert LinRange to array using collect() before saving
    write(file, "kx_range", collect(kx_range))
    write(file, "ky_range", collect(ky_range))
    write(file, "Gk2_shift", Gk2_shift)
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
