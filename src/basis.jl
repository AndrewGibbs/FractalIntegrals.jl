abstract type FractalBasisElement end

struct P0BasisElement{M<:AbstractInvariantMeasure,
                    T<:Number,
                    I<:AbstractVector{<:Integer}}
    measure :: M
    normalisation :: T
    vindex :: I
end

function subdivide_indices(Γ::AbstractAttractor, h::Real)
    I = Vector{Int64}[]
    M = length(Γ.ifs)
    r = zeros(M)

    if Γ.diam >= h
        subdiv = true
        for m in 1:M
            push!(I,[m])
            r[m] = Γ.ifs[m].ρ
        end
    else
        subdiv = false
    end

    while subdiv
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

abstract type FractalBasis{M<: AbstractInvariantMeasure, E<: AbstractVector} end

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
    Vₕ = P0Basis(μ, [P0BasisElement(μ[m], 1.0, m) for m in Lₕ])
    return Vₕ
end