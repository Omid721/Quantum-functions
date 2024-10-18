using ITensors
import TensorCrossInterpolation as TCI
using TCIITensorConversion
using QuanticsGrids
using Plots
using HDF5
include("functions.jl")

using LaTeXStrings

# Set default font family (Times New Roman) and apply individual font sizes
default(fontfamily="Times New Roman")

# Check if there are enough arguments
if length(ARGS) < 1
    println("Please provide at least one argument.")
else
    # Convert the first argument to a number (if needed)
    first_arg = parse(Float64, ARGS[1])
    println("The number of Grids for k is : ", Int(first_arg))
end

Nk = Int(first_arg)

β = 100
# Load the MPS from the HDF5 file
I = h5open("mps_$(β)_data.h5", "r") do file
    read(file, "MPS", MPS)
end

truncate!(I ;maxdim=100,cutoff=1E-11)


χ = [i for i in linkdims(I)]
χmax = maximum(χ)
println("-"^30)
println("MPS loaded with max bond dimension: $χmax")

L = 75

L1 = 25 #for n
L2 = 25 #for p
L3 = 25 #for k
β = 1000
k_max = 3 # momentum of the polarization
p_max = 2000
T = 1/β
m = 0 #The




s = [rand([1,2]) for i in range(1,L3)]
#println("Done4")
#println(" Evaluating I for random number $(evaluate(I,s))")



k_range = LinRange(0.05,k_max, Nk)
#=
for k in k_range
    x3 = c2b(k/k_max,L3)
    push!(Πk , evaluate(I,x3))
end
=#

k_normalized = k_range ./ k_max
x3_values = c2b.(k_normalized, L3)

println("-"^30)
println("start evaluating the polarization with $Nk Grids")

Πk = [(1/4)*evaluate(I, x3) for x3 in x3_values]

Πk_exact = [(1/4)*(1/q)*log((2*q+q^2)^2/(2*q-q^2)^2) for q in k_range]


default(dpi=300)
# Assuming k_range and Πk are defined
p = scatter(k_range, real.(Πk), xlabel="k",ylabel="Re(Πk)",ylim = (0,2), label="Real(Πk)")
plot!(k_range, real.(Πk_exact), label="exact_Πk",  lw=1.5)
# Save the first plot
savefig(p, "Pi_VS_E.png")


# Plot link dimensions for MPS
χ = [i for i in linkdims(I)]
p2 = plot(χ, legend=false)

# Save the second plot
savefig(p2, "link_dimensions2.png")
println("-"^60)
println("Done")
println("-"^60)
