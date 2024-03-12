abstract type Fractal end

abstract type AbstractAttractor end

struct HomogenousAttractor{ T<:Real,
                            S<:AbstractArray{<:AbstractSimilarity},
                            G<:AbstractArray{<:AbstractInvariantMap}
                            } <: AbstractAttractor
    ifs::S
    diam::T
    d::T
    n::Int64
    symmetries::G
    ρ::T
end

struct Attractor{   T<:Real,
                    S<:AbstractArray{<:AbstractSimilarity},
                    G<:AbstractArray{<:AbstractInvariantMap},
                    } <: AbstractAttractor
    ifs::S
    diam::T
    dimH::T
    n::Int64
    symmetries::G
end

function ifs_map!(  Sx::AbstractVector{<:T},
                    S::AbstractVector{<:AbstractSimilarity},
                    x::AbstractVector{<:T}) where T
    # M = length(S)
    N = length(x)
    # y = zeros(T, M*N)
    for m in eachindex(S)
        Sx[((m-1)*N+1) : (m*N)] .= S[m].(x)
    end
end

function ifs_map(   S::AbstractVector{<:AbstractSimilarity},
                    x::AbstractVector{T}) where T
    N = length(x)
    Sx = zeros(T, length(S)*N)
    for m in eachindex(S)
        @inbounds Sx[((m-1)*N+1) : (m*N)] .= S[m].(x)
    end
    return Sx
end

(Γ::AbstractAttractor)(x::AbstractVector) = ifs_map(Γ.ifs, x)

abstract type Measure{T<:Real, B, V<:AbstractVector{T}, A<:AbstractAttractor} end

# I'm wondering if its actuall worth defining different types of measures.
# But I think that HausdorffMeasure{HomogenousAttractor} will be quite useful
# also for symmetries, e.g. Koch is inhomogenous but has symmetries

struct HausdorffMeasure{T <: Real,
                        B, # probably want this to match the translation in A
                        V <: AbstractVector{T},
                        A <: AbstractAttractor
                        } <: Measure{T, B, V, A}
    supp::A
    barycentre::B
    suppmeasure::T
    weights::V
end

struct InvariantMeasure{T <: Real,
                        B, # probably want this to match the translation in A
                        V <: AbstractVector{T},
                        A <: AbstractAttractor
                        } <: Measure{T, B, V, A}
    supp::A
    barycentre::B
    suppmeasure::T
    weights::V
end

# is there much argument for defining the submeasure?
# The idea is that it points to the 'parent' anyway,
# and this is memory-cheap because the parent is an immutable struct.
# But V is an immutable struct anyway and μ is necessarilty different in each submeasure.
# So the only issue lies with A.
# dimH and n are the only major duplicates, and these are both tiny in size.
# I also assume that if I write tte constructor carefully, there will be no extra allocations.