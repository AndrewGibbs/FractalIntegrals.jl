abstract type FractalBasisElement end

struct P0BasisElement{A<:AbstractInvariantMeasure,
                    T<:Number,
                    I<:Integer,
                    V<:VectorIndex}
    measure :: A
    normalisation :: T
    index :: I
    vindex :: V

    # inner constructor check here that M of vindex matches measure
    function P0BasisElement(measure::AbstractInvariantMeasure{<:Any, M1, <:Any, <:Any},
                    normalisation::Any,
                    index::Any,
                    vindex::VectorIndex{M2, <:Any}) where {M1, M2}
        if M1 != M2
            throw(ArgumentError("Parameter 'M' of measure and vindex must match"))
        end

        # Return the constructed object
        new{typeof(measure), typeof(normalisation), typeof(index), typeof(vindex)}(measure, normalisation, index, vindex)
    end
end

(ϕₙ::P0BasisElement)(::Any) = ϕₙ.normalisation

abstract type FractalBasis{M<:AbstractInvariantMeasure} <: AbstractVector{M} end

struct P0Basis{ M <: AbstractInvariantMeasure,
                E <: AbstractVector{<:P0BasisElement}
                } <: FractalBasis{M}
    measure :: M
    elements :: E
    uniform :: Bool
end

Base.getindex(Vₙ::FractalBasis, j::Integer) = Vₙ.elements[j]
Base.length(Vₙ::FractalBasis) = length(Vₙ.elements)
Base.size(Vₙ::FractalBasis) = size(Vₙ.elements)

function construct_p0basis(μ::AbstractInvariantMeasure, h::Real)
    Lₕ = subdivide_indices(μ.supp, h)
    return P0Basis(μ,
                [P0BasisElement(μ[𝐦], 1.0, n, 𝐦) for (n, 𝐦) in enumerate(Lₕ)],
                isa(μ.supp, AbstractHomogenousAttractor))
end

# default to Hausdorff measure if an attractor is passed as first arg
@hausdorffdefault construct_p0basis

# quadrature type function - but needs to be defined after FractalBasis

# function mapquadrule_to_elements
mapquadrule_to_elements(Vₕ::FractalBasis, q::QuadStruct) =
    [QuadStruct(mapquadrule(Vₕ.measure, ϕₙ.vindex, q.nodes, q.weights)...) for ϕₙ in Vₕ]

get_h_mesh(Vₕ::FractalBasis) = maximum(ϕₙ.measure.supp.diam for ϕₙ in Vₕ)