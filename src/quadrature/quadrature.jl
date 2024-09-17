# quadrature struct, to be compactly passed around inside other structs

struct QuadStruct{T<:AbstractArray, R<:AbstractArray}
    nodes::T
    weights::R
end

include("jacobimatrices.jl")
include("productquadrature.jl")
include("barycentrerule.jl")
include("chaos.jl")
include("gauss.jl")
include("senergy.jl")

# default parameters
QUAD_DEFAULT_GAUSS = 5
QUAD_EXTRA_LEVELS = 2

getdefault_quadwidth(Γ::AbstractAttractor) =
    Γ.diam * maximum(sₘ.ρ for sₘ in Γ.ifs)^QUAD_EXTRA_LEVELS

default_barywidth(μ::AbstractInvariantMeasure) = getdefault_quadwidth(μ.supp)

default_senergy_barywidth(μ₁, μ₂) = max(default_barywidth(μ₁), default_barywidth(μ₂))


getdefault_quad_premap(μ, h_mesh, h_quad, n_quad) =
    getdefault_quad(μ, h_quad*diam(μ)/h_mesh, n_quad)


function getdefault_quad(μ::AbstractInvariantMeasure{N},
                        h_quad,
                        n_quad) where N
        if h_quad > 0
            x, w = barycentre_quadrule(μ, h_quad)
        elseif N == 1
            if n_quad >= 1
                x, w = gauss_quadrule(μ, n_quad)
            else
                x, w = gauss_quadrule(μ, QUAD_DEFAULT_GAUSS)
            end
        else
            x, w = barycentre_quadrule(μ, default_barywidth(μ))
        end
    return QuadStruct(x, w)
end

getdefault_quad_premap(μ₁, μ₂, h_mesh = max(diam(μ₁), diam(μ₂)); h_quad = 0.0, N_quad = 0) =
    combine_quadrules(  getdefault_quad_premap(μ₁, h_mesh, h_quad, N_quad),
                        getdefault_quad_premap(μ₂, h_mesh, h_quad, N_quad))


# generic quadrature function:

function mapquadrule(μ::AbstractInvariantMeasure, m::AbstractVector{<:Integer}, X, W)
    for mᵢ in reverse(m)
        X = μ.supp.ifs[mᵢ].(X)
    end
    return X, prod(μ.weights[m]).*W
end

# two-dimensional analogue
function mapquadrule(μ₁::AbstractInvariantMeasure,
                    μ₂::AbstractInvariantMeasure,
                    m::AbstractVector{<:Integer},
                    m_::AbstractVector{<:Integer},
                    X::AbstractVector,
                    Y::AbstractVector,
                    W::AbstractVector)
    for mᵢ in reverse(m)
        X = μ₁.supp.ifs[mᵢ].(X)
    end

    for mᵢ in reverse(m_)
        Y = μ₂.supp.ifs[mᵢ].(Y)
    end

    return X, Y, prod(μ₁.weights[m]).*prod(μ₂.weights[m_]).*W
end