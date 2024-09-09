# abstract type AbstractVectorIndex{M} end

# Efficient representation of vector indices of the form 𝐦 = [𝐦₁, 𝐦₂, 𝐦₃, 𝐦₄, ...]
struct VectorIndex{M, T<:UInt} <: AbstractVectorIndex{M}
    intdex :: T
end

# Convert to vector of type T
Base.Vector(𝐦::VectorIndex{M, T}) where M =
    𝐦.index !=0 ? T.(digits(𝐦.intdex, base=M) .+ 1) : [zero(T)]

# Convert a vector to VectorIndex type
VectorIndex{M,T}(𝐯::AbstractVector{T<:UInt}) = VectorIndex{M, T}(T(sum( 𝐯 .* (M .^ eachindex(𝐯)) )))

# accessing individual entries of vector index without creating full vector
Base.getindex(𝐦::VectorIndex{M, T}, i) = div(𝐦.intdex % M^(i+1), M^i) #int_to_vindex(𝐦)[i]

# extend current vector indexing to access subcomponents of attractors and measures
Base.getindex(Γ::Union{AbstractAttractor, AbstractInvariantMeasure}, 𝐦::VectorIndex) = 
    getindex(Γ, Vector(𝐦))

# output of vector indices must be as one would write it, e.g. [1,4,2,1,1,3] etc
Base.show(io::IO, 𝐦::VectorIndex{M, T}) where {M, T} = Base.show(io, int_to_vindex(𝐦))