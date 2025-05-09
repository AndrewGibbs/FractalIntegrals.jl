function singlelayer_potential_helmholtz(
    ϕ::Projection,
    k::Number;
    ambient_dimension::Integer = get_ambient_dimension(ϕ.basis.measure),#,
    # h_quad::Real = 0.0,
    # N_quad::Integer = 0,
    # quadrule::QuadStruct =
    #     getdefault_quad_premap(ϕ.basis.measure,
    #                     get_h_mesh(ϕ.basis),
    #                     h_quad,
    #                     N_quad
    #                     )
)
    ambient_dimension = check_ambient_dimension(ambient_dimension)

    # select appropriate kernel for ambient dimension
    if ambient_dimension == 2
        if get_ambient_dimension(ϕ.basis.measure) == 1
            Φ = (x, y) -> helmholtzkernel2d_screenpot(k, x, y)
        else
            Φ = (x, y) -> helmholtzkernel2d(k, x, y)
        end
    elseif ambient_dimension == 3
        if get_ambient_dimension(ϕ.basis.measure) == 2
            Φ = (x, y) -> helmholtzkernel3d_screenpot(k, x, y)
        else
            Φ = (x, y) -> helmholtzkernel3d(k, x, y)
        end
    else
        error("Haven't coded single layer potential for this many dimensions")
    end

    # # map (X,W) quadrature rule to each basis element
    # quadrules_mapped = mapquadrule_to_elements(ϕ.basis, quadrule)

    # create and return instance of potential type
    return Potential(ϕ, Φ)
end

function farfield_pattern(
    ϕ::Projection,
    k::Number;
    ambient_dimension::Integer = get_ambient_dimension(ϕ.basis.measure),#,
    # h_quad::Real = 0.0,
    # N_quad::Integer = 0,
    # quadrule::QuadStruct =
    #     getdefault_quad_premap(ϕ.basis.measure,
    #                     get_h_mesh(ϕ.basis),
    #                     h_quad,
    #                     N_quad
    #                     ),
)

    # const kwave = k

    ambient_dimension = check_ambient_dimension(ambient_dimension)

    if ambient_dimension == 2
        if get_ambient_dimension(ϕ.basis.measure) == 1
            ffkernel = (θ, y) -> -sqrt(1im / (8π * k)) * cis.(-k * (cos(θ) * y))
        else
            ffkernel =
                (θ, Y) ->
                    -sqrt(1im / (8π * k)) *
                    cis.(-k * (cos(θ) * [y[1] for y in Y] + sin(θ) * [y[2] for y in Y]))
        end
    elseif ambient_dimension == 3
        if get_ambient_dimension(ϕ.basis.measure) == 2
            ffkernel =
                ((θ, ψ), Y) ->
                    -(1 / (4π)) *
                    cis.(
                        -k * (
                            sin(θ) * cos(ψ) * [y[1] for y in Y] + sin(θ) * sin(ψ) * [y[2] for y in Y]
                        )
                    )
        else
            ffkernel =
                ((θ, ψ), Y) ->
                    -(1 / (4π)) *
                    cis.(
                        -k * (
                            sin(θ) * cos(ψ) * [y[1] for y in Y] +
                            sin(θ) * sin(ψ) * [y[2] for y in Y] +
                            cos(θ) * [y[3] for y in Y]
                        )
                    )
        end
    else
        @error "Haven't coded ffp for this many dimensions... does it even make sense?"
    end

    # # map (X,W) quadrature rule to each basis element
    # quadrules_mapped = mapquadrule_to_elements(ϕ.basis, quadrule)

    # create and return instance of potential type
    return Potential(ϕ, ffkernel)
end

function singlelayer_potential_laplace(
    ϕ::Projection;
    ambient_dimension::Integer = get_ambient_dimension(ϕ.basis.measure),#,
    # h_quad::Real = 0.0,
    # N_quad::Integer = 0,
    # quadrule::QuadStruct =
    #     getdefault_quad_premap(ϕ.basis.measure,
    #                     get_h_mesh(ϕ.basis),
    #                     h_quad,
    #                     N_quad
    #                     ),
)
    ambient_dimension = check_ambient_dimension(ambient_dimension)

    n = ambient_dimension + 1
    n == 2 ? Aₙ = -1 / (2π) : Aₙ = 2π^((n - 1) / 2) / gammafn((n - 1) / 2)

    Φ = (x, y) -> 1 / ((n - 2) * Aₙ) * energykernel(n - 2, x, y)

    # # map (X,W) quadrature rule to each basis element
    # quadrules_mapped = mapquadrule_to_elements(ϕ.basis, quadrule)

    # create and return instance of potential type
    return Potential(ϕ, Φ)
end
