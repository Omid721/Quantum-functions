
using ITensors
import TensorCrossInterpolation as TCI
using TCIITensorConversion
using QuanticsGrids
using Plots
using HDF5


println("packages uploaded")
println("-"^40)



function integrate(M::MPS)
    I = ITensor(1.)
    for (m,s) in zip(M,siteinds(M))
      I *= (m*ITensor([1/2,1/2],s))
    end
    return scalar(I)
end


function partial_integrate(psi::MPS, trace_sites::Vector{Int})
  new_psi = copy(psi)
  site_num = [i for i in 1:length(psi) if !(i in trace_sites)]
  site_indx = [siteinds(psi)[i] for i in 1:length(psi) if !(i in trace_sites)]
  if isempty(site_indx)
    return integrate(new_psi)
  else
    I = MPS(site_indx)
    first_I_site = site_num[1]
    I[1] = psi[first_I_site]
    l = first_I_site - 1
    while l >= 1
      I[1] *= new_psi[l]* ITensor([1/2, 1/2], siteinds(psi)[l])
      l -= 1
    end
    site = first_I_site + 1
    i = 1

    for (m, s) in zip(psi[first_I_site + 1:end], siteinds(psi)[first_I_site + 1:end])
        if (site in trace_sites)
            I[i] *= m * ITensor([1/2, 1/2], s)
        else
          i += 1
          I[i] =  m
        end
        site += 1
    end
    return I
  end
end



function b2c(xs...)
  n = length(xs)
  x = 0.0
  for j=1:n
    x += (xs[j]-1)/2^j
  end
  return x
end

b2c(xs::Vector) = b2c(xs...)
b2c(xs::Tuple) = b2c(xs...)


function c2b(x::Real, n::Int)
b = fill(1,n)
val = 0.0
δ = 1.0
for j=1:n
  δ /= 2
  nextval = val+δ
  if nextval <= x
    b[j] = 2
    val = nextval
  end
end
return b
end


#############################################################################################################
# writing the evaluate function. in the evalute2 function there might be 'nothing' index due to the partial_integrate

"""
Evaluate the function represented by the MPS
given a collection of settings of the site indices.
"""
function evaluate(M::MPS, v)
  V = ITensor(1.)
  s = siteinds(M)
  for (j,sj) in enumerate(s)
    if sj !== nothing
      V *= M[j]*onehot(s[j]=>v[j])
    end
  end
  return scalar(V)
end



###############################
#  upgraded evaluate function #
###############################

function evaluate2(M::MPS, v)
  i = 1
  for (j, sj) in enumerate(s)
      if sj !== nothing
          V *= M[j] * onehot(s[j] => v[i])
          i += 1
      else
        V *= M[j]
      end
  end
  return scalar(V)
end

"""
Evaluate the function represented by the MPS
at the grid point corresponding to the value x ∈ [0,1).
"""
evaluate(M::MPS, x::Float64) = evaluate(M,c2b(x,length(M)))
#############################################################################################################




"""
Collect all values taken by the QTT and
return a vector of these.
"""
function evaluate_all(M::MPS; precision=length(M))

  s = siteinds(M)
  P1 = M[1]*onehot(s[1]=>1)
  P2 = M[1]*onehot(s[1]=>2)
  (length(M) == 1) && return [scalar(P1),scalar(P2)]
  L = MPS(M[2:end])
  R = MPS(M[2:end])
  L[1] *= P1
  R[1] *= P2
  return vcat(evaluate_all(L),evaluate_all(R))
end


#####################################################


L1 = 30 #for n
L2 = 20 #for p
L3 = 10 #for k

L = Int(L1 + L2 + L3)
β = 500
k_max = 2*π # momentum of the polarization
p_max = 2*π
T = 1/β
m = 0 #The


epsilon(p:: Real) = cos(p) - μ
μ = 1
G(iω, k) = 1/(iω - epsilon(k))







println("functions defined")


function f(v)
    iν = 2 * m * π * T  * im
    x = b2c(v[1:L1])
    x1 = b2c(v[L1+1:L1+L2])
    x2 = b2c(v[L1+L2+1:L])
    n = 2^(L1 + 1) * (x - 1/2)
    iω = (2*n + 1) * π * T *im
    p = p_max * (x1 -1/2)
    k = k_max * (x2)
    return - G(iω,p) * G(iν+iω,k + p)
end

localdims = fill(2, L)  # Fill the tensor with dimensions of size 2

tolerance = 1e-11
println("START THE TCI CALCULATION")
tci, ranks, errors = TCI.crossinterpolate2(ComplexF64, f, localdims; tolerance=tolerance, maxbonddim=600)

mps = MPS(tci)
println("Done")
println("-"^20)




integrate_list_indx = [i for i in range(1,L1+L2)]

println("Done2")
I = 2* p_max * 2^(L1+1) * (1/β) * partial_integrate(mps,integrate_list_indx)
println("Done3")
s = [rand([1,2]) for i in range(1,L3)]
println("Done4")
println(" Evaluating I for random number $(evaluate(I,s))")



k_range = LinRange(0.5,k_max,100)
Πk = []
for k in k_range
    x3 = c2b(k/k_max,L3)
    push!(Πk , evaluate(I,x3))
end
# Assuming k_range and Πk are defined
p = scatter(k_range, real.(Πk), xlabel="k",ylabel="Re(Πk)", label="Real(Πk)")
# Save the first plot
savefig(p, "Pi_VS_E.png")


# Plot link dimensions for MPS
χ = [i for i in linkdims(mps)]
p2 = plot(χ, legend=false)

# Save the second plot
savefig(p2, "link_dimensions2.png")


h5open("mps_$(β)_data.h5", "w") do file
    write(file, "MPS", I)
end
