abstract type FractalBasisElement end

struct P0BasisElement{M<:AbstractInvariantMeasure,
                    T<:Number,
                    I<:Integer,
                    V<:AbstractVector{I}}
    measure :: M
    normalisation :: T
    index :: I
    vindex :: V
end

(ϕₙ::P0BasisElement)(::Any) = ϕₙ.normalisation
# (ϕₙ::P0BasisElement)(x::AbstractArray) = fill(ϕₙ.normalisation, size(x))

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

abstract type FractalBasis{M<:AbstractInvariantMeasure, E<:AbstractVector} end

struct P0Basis{ M <: AbstractInvariantMeasure,
                E <: AbstractVector{<:P0BasisElement}
                } <: FractalBasis{M, E}
    measure :: M
    elements :: E
end

Base.getindex(Vₙ::FractalBasis, j::Integer) = Vₙ.elements[j]
Base.length(Vₙ::FractalBasis) = length(Vₙ.elements)
Base.iterate(Vₙ::FractalBasis, state=1) = state > length(Vₙ) ? nothing : (Vₙ.elements[state], state+1)

function construct_p0basis(μ::AbstractInvariantMeasure, h::Real)
    Lₕ = subdivide_indices(μ.supp, h::Real)
    Vₕ = P0Basis(μ, [P0BasisElement(μ[m], 1.0, n, m) for (n,m) in enumerate(Lₕ)])
    return Vₕ
end

# quadrature type function - but needs to be defined after FractalBasis

# function mapquadrule_to_elements(Vₕ::FractalBasis, X, W)
function mapquadrule_to_elements(Vₕ::FractalBasis, q::QuadStruct)

    # allquads = Vector{typeof(X)}(undef, length(Vₕ))
    # allweights = Vector{typeof(W)}(undef, length(Vₕ))
    # quads = Vector{QuadStruct{typeof(X),typeof(W)}}(undef, length(Vₕ))
    # for (n, ϕₙ) in enumerate(Vₕ)
    #     allquads[n], allweights[n] = mapquadrule(Vₕ.measure, ϕₙ.vindex, X, W)
    # end

    # return allquads, allweights
    return [QuadStruct(mapquadrule(Vₕ.measure, ϕₙ.vindex, q.nodes, q.weights)...) for ϕₙ in Vₕ]
end

get_h_mesh(Vₕ::FractalBasis) = maximum(ϕₙ.measure.supp.diam for ϕₙ in Vₕ)