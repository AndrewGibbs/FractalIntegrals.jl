InvariantMap(δ, A) = Similarity(1.0, δ, A)

# get_identity_map(Γ::AbstractAttractor{R, <:AbstractArray}) where R =
#     TranslatingSimilarity(one(R), zeros(eltype{T}, Γ.n), one(R), one(R))
# get_identity_map(::AbstractAttractor{R, <:Number}) where R =
#     OneDimensionalSimilarity(one(R), zero(R), one(R), one(R))

Base.one(s::AbstractSimilarity{R, T}) where {R, T} =
    TranslatingSimilarity{R, T}(one(R), zeros(eltype(T), length(s.δ)), one(R), one(R))

Base.one(::AbstractSimilarity{R, T}) where {R, T<:Number} =
    OneDimensionalSimilarity{R, T}(one(R), zero(T), one(R), one(R))

# get_identity_map(R::Type, T::Type, n::Integer) = 
#     TranslatingSimilarity{R, T}(one(R), zeros(eltype(T), n), one(R), one(R))

# get_identity_map(R::Type, n::Integer) = 
#     TranslatingSimilarity{R, T}(one(R), zeros(eltype(T), n), one(R), one(R))

# get_identity_map(R::Type, T::Type, ::Integer) = 
#     OneDimensionalSimilarity{R, T}(one(R), zero(T), one(R), one(R))
# get_identity_map(n::Integer)  = get_identity_map(Float64, SVector{n, Float64}, n)

# get_identity_map(Γ::AbstractAttractor{R, T}) where {R,T} = get_identity_map(R, T, Γ.n)

# abstract type AbstractInvariantMap{T} <:AffineMap end

# struct InvariantMap{T<:SVector, M<:SMatrix} <: AbstractInvariantMap{T}
#     δ :: T
#     A :: M
# end

# # check this outer constructor fixes the problem!
# function InvariantMap(δ::AbstractVector{T1}, A::AbstractArray{T2}
#                         ) where {
#                     T1<:Number,
#                     T2<:Number}
#     T = promote_type(T1, T2)
#     n = length(δ)
#     δ_ = SVector{n,T}(δ)
#     A_ = SMatrix{n,n,T}(A)
#     return InvariantMap{typeof(δ_), typeof(A_)}(δ_, A_)
#     #InvariantMap{typeof(δ_), typeof(A_)}(δ_, A_)
# end

# struct OneDimensionalInvariantMap{T} <: AbstractInvariantMap{T}
#     δ :: T
#     A :: T
# end
# # OneDimensionalInvariantMap(δ::Number, A::Number) = OneDimensionalInvariantMap(promote(δ, A)...)

# Base.convert(::Type{OneDimensionalInvariantMap{Tout}},
#             m::OneDimensionalInvariantMap{Tin}
#             ) where {Tin <:Number, Tout <: Number} = 
#             OneDimensionalInvariantMap(convert(Tout, m.δ), convert(Tout, m.A))

# Base.promote_rule(::Type{OneDimensionalInvariantMap{T}}, ::Type{OneDimensionalInvariantMap{I}}
#                     ) where {T<:Number,I<:Number} = OneDimensionalInvariantMap{promote_type(T,I)}

# Base.one(::OneDimensionalSimilarity{T}) where {T<:Number} = 
#     OneDimensionalSimilarity(one(T), zero(T), one(T), zero(T))

# Base.one(::TranslatingSimilarity{T, V}) where {N, T<:Number, V<:SVector{N,T}} = 
#     TranslatingSimilarity(one(T), zero(V), one(T), zero(T))

# Base.one(::Similarity{T, V, M}) where {N, T<:Number, V<:SVector{N,T}, M<:SMatrix{N,N,T}} = 
#     Similarity(one(T), zero(V), one(V), zero(M))

# define the identity automorphism
# function IdentityInvariantMap(T::Type, n::Integer)
#     if n==1
#         I = OneDimensionalInvariantMap(zero(T), one(T))
#     else
#         I = InvariantMap(SVector{n,T}(zeros(n)), SMatrix{n,n,T}(Matrix(IdMat(n))))
#     end
#     return I
# end
# IdentityInvariantMap(n::Integer) = IdentityInvariantMap(Float64, n)

# function simcompsymmetries(s::AbstractSimilarity, symmetries::AbstractVector{T}) where T<:AbstractSimilarity
#     Achain = [s.A*symmetries[m].A/s.A for m in eachindex(symmetries)]
#     return [T((IdMat-Achain[m])*s.δ + s.ρA*symmetries[m].δ, Achain[m]) for m in eachindex(symmetries)]
# end

## Invariant/automoprhic maps/symmetries - focusing on 2d stuff for now

rotationmatrix2d(θ::T) where T = SMatrix{2,2,T}([cos(θ) -sin(θ); sin(θ) cos(θ)])
reflectionmatrix2d(θ::T) where T = SMatrix{2,2,T}([cos(2θ) sin(2θ); sin(2θ) -cos(2θ)])

function get_group_operations2d(num_rotations::Integer,
                                reflections::AbstractVector{T},
                                centre::AbstractVector{T}) where T<:Number
    δθ = T(2π) / num_rotations
    num_reflections = length(reflections)
    
    ivmaps = Vector{Similarity{T, SVector{2,T}, SMatrix{2,2,T,4}}
                    }(undef, num_rotations + num_reflections)#[IdentityMap(2) for _=1:(num_rotations+num_reflections)]

    counter = 0
    for θ in 0:δθ:(2π-δθ)
        counter += 1
        rotmat = rotationmatrix2d(θ)
        ivmaps[counter] = InvariantMap(centre - rotmat * centre, rotmat)
    end
    for ϑ in reflections
        counter += 1
        refmat = reflectionmatrix2d(ϑ)
        ivmaps[counter] = InvariantMap(centre - refmat * centre, refmat)
    end
    # the union of rotations and reflections gives all compositions
    return ivmaps
end

function DihedralGroup( T::Type,
                        n::Integer;
                        centre::AbstractVector{<:Number} = zeros(T, 2),
                        angle_offest::Number = zero(T))
    if mod(n,2) == 0 # for even n, want to double frequency of reflections, but restrict to [0,π]
        δθ = T(π)/n
        reflections = (0:δθ:(T(π)-δθ)) .+ angle_offest
    else
        δθ = 2π/n # for odd n, split [0,2π] into n segments
        reflections = (0:δθ:(2T(π)-δθ)) .+ angle_offest
    end
    return get_group_operations2d(n, reflections, centre)
end

DihedralGroup(n; varags...) = DihedralGroup(Float64, n; varags...)

function trivialgroup(  T::Type,
                        n::Integer)
    if n==1
        tg = OneDimensionalSimilarity(one(T), zero(T), one(T), zero(T))
    else
        tg = TranslatingSimilarity(one(T), SVector{n}(zeros(T,n)), one(T), zero(T))
    end
    return [tg]
end
trivialgroup(n::Integer) = trivialgroup(Float64, n)

# trivialgroup(args...) = [get_identity_map(args...)]
# D₂_in_1D(;centre::Float64=0.0) = [IdentityMap(1), OneDimensionalInvariantMap(-1.0, 2*centre)]