# define key default parameters

# default constants
DOFS_PER_WAVELENGTH = 20
DOFS_FOR_NONOSCILLATORS = 5

function vindex_to_scalar(M::Integer, ℓ::Integer, m::AbstractVector{<:Integer})
    n = 1
    for j=1:ℓ
        n += M^(ℓ-j) * (m[j]-1)
    end
    return n
end

include("basis.jl")
include("innerproducts.jl")
include("matrixreps.jl")
include("discreteoperators.jl")
include("projections.jl")