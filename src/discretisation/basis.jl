abstract type FractalBasisElement end
abstract type AbstractP0BasisElement <: FractalBasisElement end

struct P0BasisElement{A<:AbstractInvariantMeasure,
                    T<:Number,
                    I<:Integer,
                    V<:VectorIndex} <: AbstractP0BasisElement
    measure :: A
    normalisation :: T
    index :: I
    vindex :: V
end

struct PreQuadP0BasisElement{A<:AbstractInvariantMeasure,
                            T<:Number,
                            I<:Integer,
                            V<:VectorIndex,
                            Q<:QuadStruct} <: AbstractP0BasisElement
    measure :: A
    normalisation :: T
    index :: I
    vindex :: V
    quad :: Q

    # # inner constructor check here that M of vindex matches measure
    # function P0BasisElement(measure::AbstractInvariantMeasure{<:Any, M1, <:Any, <:Any},
    #                 normalisation::Any,
    #                 index::Any,
    #                 vindex::VectorIndex{M2, <:Any}) where {M1, M2}
    #     if M1 != M2
    #         throw(ArgumentError("Parameter 'M' of measure and vindex must match"))
    #     end

    #     # Return the constructed object
    #     new{typeof(measure), typeof(normalisation), typeof(index), typeof(vindex)}(measure, normalisation, index, vindex)
    # end
end

(ϕₙ::AbstractP0BasisElement)(::Any) = ϕₙ.normalisation

abstract type FractalBasis{M<:AbstractInvariantMeasure} <: AbstractVector{M} end

struct QuasiUniformBasis{ M <: AbstractInvariantMeasure,
                E <: FractalBasisElement,
                Q <: QuadStruct
                } <: FractalBasis{M}
    measure :: M
    elements :: Vector{E}
    parent_quadrule :: Q
    # possible h parameter here... still considering this?
end

Base.getindex(Vₙ::FractalBasis, j::Integer) = Vₙ.elements[j]
Base.length(Vₙ::FractalBasis) = length(Vₙ.elements)
Base.size(Vₙ::FractalBasis) = size(Vₙ.elements)

function construct_quasiuniform_prebary_p0basis(μ::AbstractInvariantMeasure,
                                                h_mesh::Real,
                                                h_quad::Real)
    Lₕ = subdivide_indices(μ.supp, h_mesh)
    h_quad_mod = diam(μ) * h_quad / h_mesh
    quad_nodes, quad_weights = barycentre_quadrule(μ, h_quad_mod)
    return QuasiUniformBasis(μ,
                [PreQuadP0BasisElement(μ[𝐦], # sub-measure
                                1.0, # normalisation
                                n, # scalar index
                                𝐦, # vector index
                                QuadStruct(mapquadrule(μ, 𝐦, quad_nodes, quad_weights)...)#QuadStruct(barycentre_quadrule(μ[𝐦], h_quad)...)
                                ) for (n, 𝐦) in enumerate(Lₕ)],
                QuadStruct(quad_nodes, quad_weights)
                )
end

# default to Hausdorff measure if an attractor is passed as first arg
@hausdorffdefault construct_quasiuniform_prebary_p0basis

# quadrature type function - but needs to be defined after FractalBasis

# function mapquadrule_to_elements
mapquadrule_to_elements(Vₕ::FractalBasis, q::QuadStruct) =
    [QuadStruct(mapquadrule(Vₕ.measure, ϕₙ.vindex, q.nodes, q.weights)...) for ϕₙ in Vₕ]

get_h_mesh(Vₕ::FractalBasis) = maximum(ϕₙ.measure.supp.diam for ϕₙ in Vₕ)