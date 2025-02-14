# ---------------------------------------------------------------------------------------- #
# -------------------- popular 2D fractals of measure zero ------------------------------- #
# ---------------------------------------------------------------------------------------- #

function cantorset(T = Float64; ρ = 1 / T(3))
    @assert 0 < ρ <= 0.5 "Contraction must satisfy ρ∈(0,0.5]"
    @assert isa(ρ, T) "Contraction ρ must be same type as first input"
    # note that automatic conversion T(ρ) will truncate ρ in some cases, so is best avoided

    # get iterated function system
    ifs = [Similarity(ρ, zero(T)), Similarity(ρ, 1 - ρ)]

    # get connectedness matrix, depends on if subcomponents are touching
    ρ == 1 / 2 ? connectedness = ones(Bool, 2, 2) : connectedness = Matrix(IdMat(2))

    # get symmetry group
    reflectiongroup1d = [one(ifs[1]), InvariantMap(T(1), T(-1))]

    # return fractal attractor
    return Attractor(
        ifs;
        d = log(T(2)) / log(1 / ρ),
        diam = one(T),
        symmetries = reflectiongroup1d,
        connectedness = connectedness,
    )
end

function cantordust(T = Float64; ρ = 1 / T(3))
    @assert 0 < ρ <= 0.5 "Contraction must satisfy ρ∈(0,0.5]"
    @assert isa(ρ, T) "Contraction ρ must be same type as first input"
    # note that automatic conversion T(ρ) will truncate ρ in some cases, so is best avoided

    # get iterated function system
    ifs = [
        Similarity(ρ, [zero(T), zero(T)]),
        Similarity(ρ, [1 - ρ, zero(T)]),
        Similarity(ρ, [zero(T), 1 - ρ]),
        Similarity(ρ, [1 - ρ, 1 - ρ]),
    ]

    # get connectedness matrix, depends on if subcomponents are touching
    ρ == 1 / 2 ? connectedness = ones(Bool, 4, 4) : connectedness = Matrix(IdMat(4))

    # return fractal attractor
    return Attractor(
        ifs;
        d = log(T(4)) / log(1 / ρ),
        diam = sqrt(T(2)),
        symmetries = DihedralGroup(T, 4),
        connectedness = connectedness,
    )
end

function vicsek(T = Float64; ρ = 1 / T(3))
    @assert 0 < ρ <= 1 / T(3) "Contraction must satisfy ρ∈(0,1/3]"
    @assert isa(ρ, T) "Contraction ρ must be same type as first input"
    # note that automatic conversion T(ρ) will truncate ρ in some cases, so is best avoided

    # get iterated function system
    ifs = [
        Similarity(ρ, [zero(T), zero(T)]),
        Similarity(ρ, [1 - ρ, zero(T)]),
        Similarity(ρ, [zero(T), 1 - ρ]),
        Similarity(ρ, [1 - ρ, 1 - ρ]),
        Similarity(ρ, [1 / T(2) - ρ / 2, 1 / T(2) - ρ / 2]),
    ] # keep centre around [0.5,0.5]

    # get connectedness matrix, depends on if subcomponents are touching
    if ρ == 1 / 3
        connectedness = Bool[
            1 0 0 0 0
            0 1 0 0 0
            0 0 1 0 0
            0 0 0 1 0
            1 1 1 1 1
        ]
    else
        connectedness = Matrix(IdMat(5))
    end

    # return fractal attractor
    return Attractor(
        ifs;
        d = log(T(5)) / log(1 / ρ),
        diam = sqrt(T(2)),
        symmetries = DihedralGroup(T, 4),
        connectedness = connectedness,
    )
end

function sierpinskitriangle(T = Float64; ρ = T(1 / 2))
    @assert 0 < ρ <= 0.5 "Contraction must satisfy ρ∈(0,0.5]"
    @assert isa(ρ, T) "Contraction ρ must be same type as first input"
    # note that automatic conversion T(ρ) will truncate ρ in some cases, so is best avoided

    # get iterated function system
    courage = Similarity(ρ, [zero(T), zero(T)])
    wisdom = Similarity(ρ, [1 / T(2), zero(T)])
    power = Similarity(ρ, [1 / T(4), sqrt(T(3)) / 4])
    ifs = [courage, wisdom, power]

    # get connectedness matrix, depends on if subcomponents are touching
    ρ == 1 / 2 ? connectedness = ones(Bool, 3, 3) : connectedness = Matrix(IdMat(3))

    # return fractal attractor
    return Attractor(
        ifs;
        d = log(T(5)) / log(1 / T(ρ)),
        diam = one(T),
        symmetries = DihedralGroup(T, 3),
        connectedness = connectedness,
    )
end

function kochcurve(T = Float64)

    # get iterated function system:
    ifs = [
        Similarity(1 / T(3), [0, 0]),
        Similarity(1 / T(3), [1 / T(3), 0], rotationmatrix2d(T(π) / T(3))),
        Similarity(1 / T(3), [1 / 2, sqrt(T(3)) / 6], rotationmatrix2d(-T(π) / T(3))),
        Similarity(1 / T(3), [2 / T(3), 0]),
    ]

    # get connectedness matrix:
    connectedness = [
        true true false false
        true true true false
        false true true true
        false false true true
    ]

    # get reflective symmetry group:

    horizontal_reflection_group =
        [one(ifs[1]), InvariantMap(Vector{T}([1, 0]), Matrix{T}([-1 0; 0 1]))]

    # return attractor
    return Attractor(
        ifs;
        diam = one(T),
        connectedness = connectedness,
        symmetries = horizontal_reflection_group,
    )
end

function sierpinskicarpet(T = Float64; ρ = T(1 / 3))
    @assert 0 < ρ <= 1 / 3 "Contraction must satisfy ρ∈(0,1/3]"
    @assert isa(ρ, T) "Contraction ρ must be same type as first input"

    ifs = [
        Similarity(ρ, [0, 0]),
        Similarity(ρ, [ρ, 0]),
        Similarity(ρ, [2ρ, 0]),
        Similarity(ρ, [2ρ, ρ]),
        Similarity(ρ, [2ρ, 2ρ]),
        Similarity(ρ, [ρ, 2ρ]),
        Similarity(ρ, [0, 2ρ]),
        Similarity(ρ, [0, ρ]),
    ]

    # get connectedness matrix
    connectedness = zeros(Bool, 8, 8)
    for m = 1:8
        for m_ = 1:8
            if m == m_
                connectedness[m, m] = true
            end

            # edges of subcomponents touching
            if abs(m - m_) == 1 || abs(m - m_) == 7
                connectedness[m, m_] = true
            end

            # corners of subcomponents touching
            if mod(m, 2) == 0 && mod(m_, 2) == 0 && (abs(m_ - m) == 2 || abs(m - m_) == 6)
                connectedness[m, m_] = true
            end
        end
    end

    # return attractor
    return Attractor(
        ifs;
        d = log(T(8)) / log(T(3)),
        diam = sqrt(T(2)),
        symmetries = DihedralGroup(T, 4),
        connectedness = connectedness,
    )
end
