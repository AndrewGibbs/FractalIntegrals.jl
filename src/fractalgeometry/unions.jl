# unions of attractors

struct AttractorUnion{A<:Tuple{Vararg{<:AbstractAttractor}}} <: Fractal
    attractors :: A
end

# Base.:∪(Γ₁::AbstractAttractor, Γ₂::AbstractAttractor) = AttractorUnion((Γ₁, Γ₂))
Base.:∪(Γs::AbstractAttractor...) = AttractorUnion(embed_into_same_dimension(tuple(Γs...)))
Base.:∪(Γunion::AttractorUnion, Γs::AbstractAttractor...) = ∪(Γunion.attractors..., Γs...)
Base.getindex(Γunion::AttractorUnion, m::Integer) = Γunion.attractors[m]
Base.length(Γ::AttractorUnion) = length(Γ.attractors)


# unions of measures

struct MeasureUnion{A<:Tuple{Vararg{<:AbstractInvariantMeasure}}} <: Measure
    measures :: A
end

# Base.:∪(μ₁::AbstractInvariantMeasure, μ₂::AbstractInvariantMeasure) = MeasureUnion((μ₁, μ₂))
Base.:∪(μs::AbstractInvariantMeasure...) = MeasureUnion(embed_into_same_dimension(tuple(μs...)))
Base.:∪(μnion::MeasureUnion, μs::AbstractInvariantMeasure...) = ∪(μnion.attractors..., μs...)
Base.getindex(μnion::MeasureUnion, m::Integer) = μnion.measures[m]
Base.length(μ::MeasureUnion) = length(μ.measures)