
# split this into a barycentric-centric version and a more general version

function discretise(sio::AbstractSingularIntegralOperator;
                    h_mesh::Number = sio.measure.supp.diam/5,
                    h_quad::Number = h_mesh/5)
    Vₕ = construct_p0basis(sio.measure, h_mesh)
    bip = BarycentreHomogInnerProduct(sio, Vₕ, h_mesh, h_quad)
    return discretise(sio, bip, Vₕ)
end

function discretise(sio::AbstractSingularIntegralOperator,
                    bip::BarycentreHomogInnerProduct,
                    Vₕ::FractalBasis
                    )
    N = length(Vₕ)
    galerkinmatrix = zeros(typeof(sio.singularconst), N, N)
    @inbounds for n in 1:N
        @inbounds for m in n:N
            galerkinmatrix[m, n] = sesquilinearform(bip, Vₕ[m], Vₕ[n])
        end
    end

    # most efficient way to invert matrix, also allows to do everything column major
    return Symmetric(galerkinmatrix, :L)
end