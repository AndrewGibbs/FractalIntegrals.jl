# guy here
# https://discourse.julialang.org/t/staticarrays-parametric-types-type-stability/59095
# has same issue as me

# abstract types

abstract type AffineMap end
abstract type AbstractSimilarity <: AffineMap end

# concrete types

struct Similarity{T<:Number, V<:AbstractVector{T}, M<:AbstractMatrix{T}} <: AbstractSimilarity
    ρ :: T
    δ :: V
    A :: M
    ρA :: M
end

struct TranslatingSimilarity{T<:Number, V<:AbstractVector{T}} <: AbstractSimilarity
    ρ :: T
    δ :: V
    A :: T
    ρA :: T
end

struct OneDimensionalSimilarity{T<:Number} <: AbstractSimilarity
    ρ :: T
    δ :: T
    A :: T
    ρA :: T
end

# outer constructors
function Similarity(ρ::T1, δ::AbstractVector{T2}) where {T1<:Number, T2<:Number}
    @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    T = promote_type(T1,T2)
    return TranslatingSimilarity(T(ρ), SVector{length(δ),T}(δ), one(T), T(ρ))
end

function Similarity(ρ::T1, δ::AbstractVector{T2}, A::AbstractMatrix{T3}) where {T1<:Number, T2<:Number, T3<:Number}
    @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    @assert abs(det(A)) ≈ 1 "Matrix (third argument) must have determinant one"
    N = length(δ)
    @assert (N,N) == size(A) "Matrix dimension must match translation vector"
    T = promote_type(T1,T2,T3)
    #{T,SVector{N,T},SMatrix{N,N,T,N^2}(A)}
    return Similarity(T(ρ), SVector{N,T}(δ), SMatrix{N,N,T}(A), SMatrix{N,N,T}(ρ*A))
end

function Similarity(ρ::T1, δ::T2, r::T3=one(T1)) where {T1<:Number, T2<:Number, T3<:Number}
    @assert abs(ρ) < 1 "Contraction (first argument) must be less than one"
    @assert abs(r) ≈ 1 "Reflection (third argument) must have length one"
    T = promote_type(T1,T2,T3)
    return OneDimensionalSimilarity(T(ρ), T(δ), T(r), convert(T,ρ*r))
end

# promotion of non-rotating map to rotating map (with identity rotation matrix)

function Base.convert(::Type{Similarity{T,V,M}}, s::TranslatingSimilarity{T,V}) where {T<:Number,V<:AbstractVector,M<:AbstractMatrix}
    A = convert(M,IdMat(length(s.δ)))
    return Similarity(s.ρ, s.δ, A, s.ρ*A)
end

# Base.promote_rule(Type{Similarity{T, V, M}}, Type{TranslatingSimilarity{T, V}}) where {T<:Number,V<:AbstractVector,M<:AbstractMatrix} = Similarity{T, V, M}
Base.promote_rule(::Type{Similarity{T, V, M}}, ::Type{TranslatingSimilarity{T, V}}) where {T<:Number,V<:AbstractVector,M<:AbstractMatrix} = Similarity{T, V, M}

# composition of two similarities
simcomp(s₁::S,s₂::S) where S<:AbstractSimilarity = S(s₁.ρ*s₂.ρ, s₁.δ+s₁.ρA*s₂.δ, s₁.A*s₂.A, s₁.ρA*s₂.ρA)
simcomp(s₁::TranslatingSimilarity,s₂::Similarity) = Similarity(s₁.ρ*s₂.ρ, s₁.δ+s₁.ρ*s₂.δ, s₂.A, s₂.ρA)
# (does not commute)
simcomp(s₁::Similarity,s₂::TranslatingSimilarity) = Similarity(s₁.ρ*s₂.ρ, s₁.δ+s₁.ρA*s₂.δ, s₁.A, s₁.ρA)
simcomp(s₁::TranslatingSimilarity,s₂::TranslatingSimilarity) = TranslatingSimilarity(s₁.ρ*s₂.ρ, s₁.δ+s₁.ρA*s₂.δ, one(typeof(s₁.ρ*s₂.ρ)), s₁.ρ*s₂.ρ)
# one dimensional case
# simcomp(s₁::OneDimensionalSimilarity, s₂::OneDimensionalSimilarity) = OneDimensionalSimilarity(s₁.ρ*s₂.ρ, s₁.δ+s₁.ρA*s₂.δ, s₁.A*s₂.A, s₁.ρA*s₂.ρA)

# Similarities as maps, optimised as far as possible

# similarity acting on an IFS
function simcompifs(s::T, ifs::AbstractVector{T}) where T<:AbstractSimilarity
    Achain = [s.A*ifs[m].A/s.A for m in eachindex(ifs)]
    rAchain = [ifs[m].ρ*Achain[m] for m in eachindex(ifs)]
    return [T(ifs[m].ρ, (IdMat-rAchain[m])*s.δ + s.ρA*ifs[m].δ, Achain[m], rAchain[m]) for m in eachindex(ifs)]
end

# simplified case without rotation
function simcompifs(s::TranslatingSimilarity, IFS::Vector{<:TranslatingSimilarity})
    return [TranslatingSimilarity(IFS[m].ρ, (IdMat-IFS[m].ρA)*s.δ + s.ρA*IFS[m].δ, IFS[m].A, IFS[m].ρA) for m in eachindex(IFS)]
end

(s::AbstractSimilarity)(x)  = s.ρA*x + s.δ
(s::TranslatingSimilarity)(x) = s.ρ*convert(typeof(s.δ),x) + s.δ

# nested composition:
function simmulticomp(IFS::Vector{<:AbstractSimilarity}, m::AbstractVector{<:Integer})
    s = IFS[m[1]]
        for j in 2:length(m)
            s = simcomp(s,IFS[m[j]])
        end
    return s
end

#inverse map
sim_map_inv(s::Similarity, x) = s.rA/(x-s.δ)

# determine how Similarity appears in the REPL

info_string(s::TranslatingSimilarity) = "x ↦ "*string(round.(s.ρ,digits=2))*"x + "*string(round.(s.δ,digits=2))
info_string(s::Similarity) = "x ↦ "*string(round.(s.ρ,digits=2))*string(round.(s.A,digits=2))*"x + "*string(round.(s.δ,digits=2))
info_string(s::OneDimensionalSimilarity) = "x ↦ "*string(round(s.ρA,digits=2))*"x + "*string(round(s.δ,digits=2))

Base.show(io::IO, s::AbstractSimilarity)  = print(io,'\n',typeof(s),':','\n', info_string(s))

# affine maps

abstract type AbstractInvariantMap <:AffineMap end

struct InvariantMap{T<:Number, V<:AbstractVector{T}, M<:AbstractMatrix{T}} <: AbstractInvariantMap
    δ :: V
    A :: M
end


# check this outer constructor fixes the problem!
function InvariantMap(δ::AbstractVector{T1}, A::AbstractArray{T2}
                        ) where {
                    T1<:Number,
                    T2<:Number}
    T = promote_type(T1, T2)
    n = length(δ)
    δ_ = SVector{n,T}(δ)
    A_ = SMatrix{n,n,T}(A)
    return InvariantMap{T, typeof(δ_), typeof(A_)}(δ_, A_)
end

struct OneDimensionalInvariantMap{T<:Number} <: AbstractInvariantMap
    δ :: T
    A :: T
end
OneDimensionalInvariantMap(δ::Number, A::Number) = OneDimensionalInvariantMap(promote(δ, A)...)

Base.convert(::Type{OneDimensionalInvariantMap{Tout}},
            m::OneDimensionalInvariantMap{Tin}
            ) where {Tin <:Number, Tout <: Number} = 
            OneDimensionalInvariantMap(convert(Tout, m.δ), convert(Tout, m.A))

Base.promote_rule(::Type{OneDimensionalInvariantMap{T}}, ::Type{OneDimensionalInvariantMap{I}}) where {T<:Number,I<:Number} = OneDimensionalInvariantMap{promote_type(T,I)}

# need to define identity similarity
# need to define 'one' equivalent for these three, which unfornately needs extra input
# Have removed function below because it is not type stable
# function IdentitySimilarity(T::Type, n::Integer)
#     if n==1
#         I = OneDimensionalSimilarity(one(T), zero(T), one(T), zero(T))
#     else
#         I = TranslatingSimilarity(one(T), # contraction
#                                 SVector{n,T}(zeros(n)), # translation
#                                 one(T), # rotation
#                                 one(T)) # rotation*contraction
#     end
#     return I
# end

Base.one(::OneDimensionalSimilarity{T}) where {T<:Number} = 
    OneDimensionalSimilarity(one(T), zero(T), one(T), zero(T))

Base.one(::TranslatingSimilarity{T, V}) where {N, T<:Number, V<:SVector{N,T}} = 
    TranslatingSimilarity(one(T), zero(V), one(T), zero(T))

Base.one(::Similarity{T, V, M}) where {N, T<:Number, V<:SVector{N,T}, M<:SMatrix{N,N,T}} = 
    Similarity(one(T), zero(V), one(V), zero(M))

# define the identity automorphism
function IdentityInvariantMap(T::Type, n::Integer)
    if n==1
        I = OneDimensionalInvariantMap(zero(T), one(T))
    else
        I = InvariantMap(SVector{n,T}(zeros(n)), SMatrix{n,n,T}(Matrix(IdMat(n))))
    end
    return I
end
IdentityInvariantMap(n::Integer) = IdentityInvariantMap(Float64, n)

function simcompsymmetries(s::AbstractSimilarity, symmetries::AbstractVector{T}) where T<:AbstractInvariantMap
    Achain = [s.A*symmetries[m].A/s.A for m in eachindex(symmetries)]
    # rAchain = [symmetries[m].ρ*Achain[m] for m in eachindex(symmetries)]
    return [T((IdMat-Achain[m])*s.δ + s.ρA*symmetries[m].δ, Achain[m]) for m in eachindex(symmetries)]
end

## Invariant/automoprhic maps/symmetries - focusing on 2d stuff for now

rotationmatrix2d(θ::T) where T = SMatrix{2,2,T}([cos(θ) -sin(θ); sin(θ) cos(θ)])
reflectionmatrix2d(θ::T) where T = SMatrix{2,2,T}([cos(2θ) sin(2θ); sin(2θ) -cos(2θ)])

function get_group_operations2d(num_rotations::Integer,
                                reflections::AbstractVector{<:Real},
                                centre::AbstractVector{T}) where T<:Number
    δθ = 2π / num_rotations
    num_reflections = length(reflections)
    
    ivmaps = Vector{InvariantMap{T, SVector{2,T}, SMatrix{2,2,T,4}}
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

function DihedralGroup( n::Integer;
                        centre::AbstractVector{<:Number} = [0.0,0.0],
                        angle_offest::Number = 0.0)
    if mod(n,2) == 0 # for even n, want to double frequency of reflections, but restrict to [0,π]
        δθ = π/n
        reflections = (0:δθ:(π-δθ)) .+ angle_offest
    else
        δθ = 2π/n # for odd n, split [0,2π] into n segments
        reflections = (0:δθ:(2π-δθ)) .+ angle_offest
    end
    return get_group_operations2d(n, reflections,centre)
end

trivialgroup(args...) = [IdentityInvariantMap(args...)]
# D₂_in_1D(;centre::Float64=0.0) = [IdentityMap(1), OneDimensionalInvariantMap(-1.0, 2*centre)]