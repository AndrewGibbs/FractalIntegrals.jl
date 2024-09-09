# quadrature struct, to be compactly passed around inside other structs

struct QuadStruct{T<:AbstractArray, R<:AbstractArray}
    nodes::T
    weights::R
end

# putting this function here for now...
function subdivide_indices(Γ::AbstractAttractor, h::Real, max_num_indices = Inf)
    I = Vector{Int64}[]
    M = length(Γ.ifs)
    r = zeros(M)

    @assert (h>0 || max_num_indices<Inf
            ) "either meshwidth must be positive, or max_num_indices must be finite"

    if Γ.diam >= h
        subdiv = true
        for m in 1:M
            push!(I,[m])
            r[m] = Γ.ifs[m].ρ
        end
    else
        subdiv = false
    end

    while subdiv && (length(I)<max_num_indices)
        subdiv = false
        split_vecs = Int64[]
        for j in eachindex(I)
           if Γ.diam*prod(r[I[j]]) >= h
                subdiv = true
                push!(split_vecs,j)
            end
        end
        if subdiv
            new_vecs = Vector{Int64}[]
            for j in split_vecs
                for m in 1:M
                    push!(new_vecs, vcat(I[j],[m]))
                end
            end
            deleteat!(I,split_vecs)
            I = vcat(I,new_vecs)
        end
    end
    #quick bodge - this convention means we can keep the same type
        # and it's (more) consistent with the paper
    if isempty(I)
        I = [[0]]
    end
    return I
end

# include("jacobimatrices.jl")
include("productquadrature.jl")
include("barycentrerule.jl")
# include("gauss.jl")
# include("senergy.jl")

# default parameters
QUAD_DEFAULT_GAUSS = 5
QUAD_EXTRA_LEVELS = 2

getdefault_quadwidth(Γ::AbstractAttractor) =
    Γ.diam * maximum(sₘ.ρ for sₘ in Γ.ifs)^QUAD_EXTRA_LEVELS

default_barywidth(μ::AbstractInvariantMeasure) = getdefault_quadwidth(μ.supp)

default_senergy_barywidth(μ₁, μ₂) = max(default_barywidth(μ₁), default_barywidth(μ₂))

getdefault_quad_premap(μ₁, μ₂, h_mesh = max(diam(μ₁),diam(μ₂)); h_quad = 0.0, N_quad = 0) =
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