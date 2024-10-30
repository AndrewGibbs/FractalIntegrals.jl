# ------------------ basis element ----------------------------------#

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
    quadrule :: Q

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

# ------------------ basis -------------------------------------- #
abstract type FractalBasis{M <: Measure} <: AbstractVector{M} end
abstract type InvariantMeasureBasis{M<:AbstractInvariantMeasure} <: FractalBasis{M} end

struct QuasiUniformBasis{ M <: AbstractInvariantMeasure,
                E <: FractalBasisElement,
                Q <: QuadStruct
                } <: InvariantMeasureBasis{M}
    measure :: M
    elements :: Vector{E}
    parent_quadrule :: Q
    # possible h parameter here... still considering this?
end

Base.getindex(Vₙ::FractalBasis, j::Integer) = Vₙ.elements[j]
Base.length(Vₙ::FractalBasis) = length(Vₙ.elements)
Base.size(Vₙ::FractalBasis) = size(Vₙ.elements)


function construct_quasiuniform_p0basis(μ::AbstractInvariantMeasure,
                                        h_mesh::Real,
                                        quadrule::QuadStruct)
    Lₕ = subdivide_indices(μ.supp, h_mesh)
    return QuasiUniformBasis(μ,
                            [PreQuadP0BasisElement(
                                μ[𝐦], # sub-measure
                                1.0/sqrt(μ[𝐦].suppmeasure), # normalisation
                                n, # scalar index
                                𝐦, # vector index
                                QuadStruct(mapquadrule(μ, 𝐦, quadrule.nodes, quadrule.weights)...)
                                ) for (n, 𝐦) in enumerate(Lₕ)
                            ],
                            quadrule
                            )
end

construct_quasiuniform_prebary_p0basis(μ::AbstractInvariantMeasure,
                                                h_mesh::Real,
                                                h_quad::Real) =
    construct_quasiuniform_p0basis( μ,
                                    h_mesh,
                                    QuadStruct(barycentre_quadrule(μ, diam(μ) * h_quad / h_mesh)...))

# default to Hausdorff measure if an attractor is passed as first arg
@hausdorffdefault construct_quasiuniform_prebary_p0basis

# quadrature type function - but needs to be defined after FractalBasis

# function mapquadrule_to_elements
mapquadrule_to_elements(Vₕ::FractalBasis, q::QuadStruct) =
    [QuadStruct(mapquadrule(Vₕ.measure, ϕₙ.vindex, q.nodes, q.weights)...) for ϕₙ in Vₕ]

get_h_mesh(Vₕ::FractalBasis) = maximum(ϕₙ.measure.supp.diam for ϕₙ in Vₕ)

# --------------------- define union of bases -------------------------------- #
struct UnionBasis{  M<:MeasureUnion,
                    B<:Tuple{Vararg{FractalBasis}}} <: FractalBasis{M}
    measures :: M
    bases :: B
    basis_lengths :: Vector{Int64}
end

Base.:∪(a::InvariantMeasureBasis, b::InvariantMeasureBasis) =
    UnionBasis(∪(a.measure, b.measure), (a, b), length.([a, b]))

function Base.getindex(b::UnionBasis, n::Integer)
    @assert n>0 "basis index must be greater than zero"
    cum_lens = vcat([0], cumsum([length.(b.bases)...]))
    # initialise as final basis entry
    basis_index = cum_lens[end] + length(b.bases[end])
    el_index = 0
    for m in 1:(length(cum_lens)-1)
        if cum_lens[m] < n ≤ cum_lens[m+1]
            basis_index = m
            el_index = n - cum_lens[m]
            break
        end
    end
    return b.bases[basis_index][el_index]
end

Base.length(b::UnionBasis) = sum(length.(b.bases))
Base.size(b::UnionBasis) = (length(b), )