include("jacobimatrices.jl")
include("productquadrature.jl")
include("barycentrerule.jl")
include("gauss.jl")
include("senergy.jl")

# generic quadrature function:

function mapquadrule(μ::AbstractInvariantMeasure, m::AbstractVector{<:Integer}, X, W)
    for mᵢ in reverse(m)
        X = μ.supp.ifs[mᵢ].(X)
    end
    return X, prod(μ.weights[m]).*W
end