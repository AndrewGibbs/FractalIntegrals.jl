abstract type DiscreteFractalOperator end

struct DiscreteGalerkinOperator{
    K <: FractalOperator,
    B <: FractalBasis,
    M <: AbstractMatrix
    } <: DiscreteFractalOperator
    op :: K
    basis :: B
    stiffness_matrix :: M
end

getdefault_meshwidth(sio::OscillatorySeparableIntegralOperator) =
    2π / abs(DOFS_PER_WAVELENGTH*sio.wavenumber)

getdefault_meshwidth(sio::IntegralOperator) =
    sio.measure.supp.diam / DOFS_FOR_NONOSCILLATORS

# new style which allows us to try different quadrature rules
function discretise(sio::AbstractSingularIntegralOperator;
            h_mesh::Real = getdefault_meshwidth(sio),#sio.measure.supp.diam/5,
            h_quad::Real = 0.0, # quick option for Barycentre rule
            N_quad::Integer = 0,
            quadrule::QuadStruct =
            getdefault_quad_premap(sio.measure, h_mesh, h_quad, N_quad),
            kwargs...)
    Vₕ = construct_p0basis(sio.measure, h_mesh)
    return discretise(sio, Vₕ; kwargs...)
end

function convert_quad_to_tuple(Q::QuadStruct)
    (; nodes, weights) = Q
    return (nodes, weights)
end
