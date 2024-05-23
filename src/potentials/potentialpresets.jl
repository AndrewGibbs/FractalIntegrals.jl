function singlelayer_potential_helmholtz(ϕ::Projection,
                                        k::Number;
                                        ambient_dimension::Integer = ϕ.basis.measure.supp.n,
                                        h_quad::Real = 0.0,
                                        N_quad::Integer = 0,
                                        quadrule::Tuple{AbstractVector, AbstractVector} =
                                            getdefault_quad(ϕ.basis.measure,
                                                            h_quad = h_quad,
                                                            N_quad = N_quad
                                                            ),
                                        )
    # select appropriate kernel for ambient dimension
    if ambient_dimension == 2
        Φ = (x,y) -> helmholtzkernel2d(k, x, y)
    elseif ambient_dimension == 3
        Φ = (x,y) -> helmholtzkernel3d(k, x, y)
    else
        error("Haven't coded single layer potential for this many dimensions")
    end

    # create and return instance of potential type
    return Potential(ϕ, Φ, quadrule)
end

function singlelayer_potential_laplace(ϕ::Projection,
                                        ambient_dimension::Integer = ϕ.basis.measure.supp.n,
                                        h_quad::Real = 0.0,
                                        N_quad::Integer = 0,
                                        quadrule::Tuple{AbstractVector, AbstractVector} =
                                            getdefault_quad(ϕ.basis.measure,
                                                            h_quad = h_quad,
                                                            N_quad = N_quad
                                                            ),
                                        )
    # select appropriate kernel for ambient dimension
    @assert ambient_dimension <= 1 "Ambeint dimension must be greater than one"

    n = ambient_dimension + 1
    n == 2 ? Aₙ = -1/(2π) :  Aₙ = 2π^((n-1)/2)/gammafn((n-1)/2)

    Φ = (x,y) -> 1/((n-2)*Aₙ)*energykernel(n-2, x, y)

    # create and return instance of potential type
    return Potential(ϕ, Φ, quadrule)
end