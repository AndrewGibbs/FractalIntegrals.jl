# Efficient representation of vector indices of the form 𝐦 = [𝐦₁, 𝐦₂, 𝐦₃, 𝐦₄, ...]
struct VectorIndex{M, T<:Integer} <: AbstractVector{T}
    intdex::T
    ℓ::T
end

# some AbstractVector properties
Base.length(𝐦::VectorIndex{M, T}) where {M, T} = 𝐦.ℓ
Base.size(𝐦::VectorIndex) = (length(𝐦),)

# Convert to vector of type T
function Base.Vector(𝐦::VectorIndex{M, T}) where {M, T}
    𝐯_main = 𝐦.intdex !=0 ? T.(digits(𝐦.intdex - 1, base=M) .+ 1) : [zero(T)]
    𝐯_pad = ones(T, 𝐦.ℓ-length(𝐯_main))
    return vcat(𝐯_main, 𝐯_pad)
end

# Convert a vector to VectorIndex type
function VectorIndex{M}(𝐯::AbstractVector{T}) where {M, T}
    @assert (minimum(𝐯) > 0) || (𝐯 == [0]) "If vector is not trivial [0], all entries must be positive"
    @assert (maximum(𝐯) <= M) "Max entry must be ≤$M "
    VectorIndex{M, T}(T(sum( (𝐯 .- 1) .* (M .^ (eachindex(𝐯) .- 1)) ) + 1), length(𝐯))
end

# accessing individual entries of vector index without creating full vector
function Base.getindex(𝐦::VectorIndex{M, T}, i::Int) where {M, T}
    @boundscheck 0 < i ≤ 𝐦.ℓ
    div((𝐦.intdex-1) % M^i, M^(i-1)) + 1 #int_to_vindex(𝐦)[i]
end

# subdividing fractals corresponds to the following operation
split(𝐦::VectorIndex{M, T}) where {M, T} =
    𝐦.intdex == 0 ?
    [VectorIndex{M, T}(𝐦_, 1) for 𝐦_ in 1:M] : 
    [VectorIndex{M, T}(𝐦_, 𝐦.ℓ+1) for 𝐦_ in 𝐦.intdex .+ ((M^𝐦.ℓ) .* (0:(M-1)))]

# extend current vector indexing to access subcomponents of attractors and measures
Base.getindex(Γ::Union{AbstractAttractor, AbstractInvariantMeasure}, 𝐦::VectorIndex) = 
    getindex(Γ, Vector(𝐦))

# output of vector indices must be as one would write it, e.g. [1,4,2,1,1,3] etc
Base.show(io::IO, 𝐦::VectorIndex{M, T}) where {M, T} = Base.show(io, Vector(𝐦))

# subdivide attractor and get array of VectorIndex 
function subdivide_indices( Γ::AbstractAttractor{N, M, T},
                            h::Real;
                            max_num_indices = Inf
                            ) where {N, M, T}


    @assert (h>0 || max_num_indices<Inf
    ) "either meshwidth must be positive, or max_num_indices must be finite"

    Lₕ = [VectorIndex{M}([zero(T)])]
    ρs = [s.ρ for s in Γ.ifs]
    diams = [diam(Γ)]

    keep_subdividing = true
    while keep_subdividing && length(Lₕ)<max_num_indices
        split_vecs = Int64[]
        keep_subdividing = false
        for j in eachindex(Lₕ)
            if diams[j] ≥ h
                keep_subdividing = true
                # append new vector index onto end of index set
                append!(Lₕ, split(Lₕ[j]))
                # similar for diameter logging
                append!(diams, diams[j].*ρs)

                # log current indices to delete
                push!(split_vecs,j)
            end
        end

        # delete vectors and radii from split components
        deleteat!(Lₕ, split_vecs)
        deleteat!(diams, split_vecs)
    end

    return Lₕ
end