
# split this into a barycentric-centric version and a more general version
struct DiscreteFractalOperator{FO<:FractalOperator,
                                IP<:AbstractInnerProduct,
                                B<:FractalBasis,
                                M<:AbstractMatrix
                                } <: FractalOperator
    op :: FO
    ip :: IP
    basis :: B
    galerkinmatrix :: M
end

function discretise(sio::AbstractSingularIntegralOperator;
                    h_mesh::Number = sio.measure.supp.diam/5,
                    h_quad::Number = h_mesh/5)
    Vₕ = construct_p0basis(sio.measure, h_mesh)
    bip = BarycentreHomogInnerProduct(sio, Vₕ, h_mesh, h_quad)
    return discretise(sio, bip, Vₕ)
end

function discretise(sio::AbstractSingularIntegralOperator,
                    ip::AbstractInnerProduct,
                    Vₕ::FractalBasis
                    )
    N = length(Vₕ)
    galerkinmatrix = zeros(typeof(sio.singularconst), N, N)
    @sync for n in 1:N
        @spawn begin
            for m in n:N
                @inbounds galerkinmatrix[m, n] = sesquilinearform(ip, Vₕ[m], Vₕ[n])
            end
        end
    end

    # store matrix as symmetric:
    # most efficient way to invert matrix, also allows to do everything column major

    return DiscreteFractalOperator(sio, ip, Vₕ, Symmetric(galerkinmatrix, :L))
end

struct projection{B<:FractalBasis, V<:AbstractVector}
    basis::B
    coeffs::V
end

function project(bip::BarycentreHomogInnerProduct,
                Vₕ::FractalBasis, 
                f::Function)
    # non-iterate version:
    #Vₕ[n])/Vₕ[n].normalisation for n in eachindex[Vₕ]
    coeffs = [innerproduct(bip, f, ϕₕ)/ϕₕ.normalisation for ϕₕ ∈ Vₕ]
    return projection(Vₕ, coeffs)
end

function Base.:\(op::DiscreteFractalOperator, f::Function)
    fₕ = project(   op.ip, # inner product
                    op.basis, # basis
                    f)
    coeffs = op.galerkinmatrix \ fₕ.coeffs
    return projection(op.basis, coeffs)
end