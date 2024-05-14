include("jacobimatrices.jl")
include("productquadrature.jl")
include("barycentrerule.jl")
include("gauss.jl")
include("senergy.jl")

# default parameters
QUAD_DEFAULT_GAUSS = 5
QUAD_EXTRA_LEVELS = 2

default_barywidth(μ::AbstractInvariantMeasure) =
    maximum(sₘ.ρ for sₘ in μ.supp.ifs)^QUAD_EXTRA_LEVELS

default_senergy_barywidth(μ₁, μ₂) = max(default_barywidth(μ₁), default_barywidth(μ₂))

function getdefault_quad(μ::AbstractInvariantMeasure;
                        h_quad::Real = 0.0,
                        N_quad::Integer = 0)
    if h_quad > 0
        x, w = barycentre_quadrule(μ, h_quad)
    elseif μ.supp.n == 1
        if N_quad >= 1
            x, w = gauss_quadrule(μ, N_quad)
        else
            x, w = gauss_quadrule(μ, QUAD_DEFAULT_GAUSS)
        end
    else
        x, w = barycentre_quadrule(μ, default_barywidth(μ))
    end
    return x, w
end

getdefault_quad(μ₁, μ₂; h_quad = 0.0, N_quad = 0) =
    combine_quadrules(  getdefault_quad(μ₁; h_quad = h_quad, N_quad = N_quad)...,
                        getdefault_quad(μ₂; h_quad = h_quad, N_quad = N_quad)...)

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