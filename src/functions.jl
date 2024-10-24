using ITensors
import TensorCrossInterpolation as TCI
using TCIITensorConversion
using QuanticsGrids
using Plots
using HDF5


println("packages uploaded")
println("-"^30)

##############################################################################
#  You can visit the series of talks                                         #
# "Tensor network for machine learning applications" by miles stoudenmire    #
# to undrestand how the following functions are define.                      # 
##############################################################################


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
  #TODO: implement precision keyword argument
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


########################################
#    New fucitons
########################################
function MPS2MPO(mps::MPS)
    n = length(mps)
    C = MPO(n)
    i = 1
    for (m,s) in zip(mps, siteinds(mps))
        C[i] = m * delta(s,s',s'')
        i += 1
    end
    return C
end




function site_evaluate(M::MPS, v, s :: Vector{<:Index})
    M1 = copy(M)
    @assert length(v) == length(s)
    for (j,sj) in enumerate(s)
      M1[j] = M[j]*onehot(s[j]=>v[j])
    end
    return M1
end
function site_evaluate_MPO(M::MPO, v, s :: Vector{<:Index})
    M1 = copy(M)
    @assert length(v) == length(s)
    for (j,sj) in enumerate(s)
      M1[j] = M[j]*onehot(s[j]=>v[j])
    end
    return M1
end


##############################
#  This fucntion does the MPS contraction from left O-O-O-O-O
#                                                   | | | | |
#                                                   O-O-O-O-O-O-O
############################
function zipup_MPS_MPS_from_left(Mps1::MPS, Mps2::MPS)
    M1 = copy(Mps1)
    M2 = copy(Mps2)
    n1 = length(M1)
    n2 = length(M2)
    @assert n1 <= n2
    #changing the M1 index
    for j in range(1,n1)
        M1[j] = M1[j]* delta(siteinds(M1)[j],siteinds(M2)[j])
    end
    M = MPS(n2-n1)
    Clinked = M1[1]*M2[1]
    for i in range(2,n1)
        Clinked *= M1[i]
        Clinked *= M2[i]
    end
    M[1] = Clinked * M2[n1+1]
    for i in range(n1+2,n2)
        M[i-n1] = M2[i]
    end
    return M
end

println("functions uploaded")
println("-"^30)
