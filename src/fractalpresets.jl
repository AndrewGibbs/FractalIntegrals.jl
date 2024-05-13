# ---------------------------------------------------------------------------------------- #
# -------------------- popular 2D fractals of measure zero ------------------------------- #
# ---------------------------------------------------------------------------------------- #

function cantorset(T = Float64; ρ = 1/T(3))
    @assert 0 < ρ <= 0.5 "Contraction must satisfy ρ∈(0,0.5]"
    @assert isa(ρ, T) "Contraction ρ must be same type as first input"
    # note that automatic conversion T(ρ) will truncate ρ in some cases, so is best avoided

    # get iterated function system
    ifs = [Similarity(ρ, zero(T)), Similarity(ρ, 1-ρ)]

    # get connectedness matrix, depends on if subcomponents are touching
    ρ == 1/2 ? connectedness = ones(Bool,2,2) : connectedness = Matrix(IdMat(2))

    # get symmetry group
    reflectiongroup1d = [IdentityInvariantMap(T, 1), OneDimensionalInvariantMap(T(1), T(-1))]

    # return fractal attractor
    return Attractor(   ifs,
                        d = log(T(2))/log(1/ρ),
                        diam = one(T),
                        symmetries = reflectiongroup1d,
                        connectedness = connectedness)
end

function cantordust(T = Float64; ρ = 1/T(3))
    @assert 0 < ρ <= 0.5 "Contraction must satisfy ρ∈(0,0.5]"
    @assert isa(ρ, T) "Contraction ρ must be same type as first input"
    # note that automatic conversion T(ρ) will truncate ρ in some cases, so is best avoided

    # get iterated function system
    ifs = [ Similarity(ρ, [zero(T), zero(T)]),
            Similarity(ρ, [1-ρ, zero(T)]),
            Similarity(ρ, [zero(T), 1-ρ]),
            Similarity(ρ, [1-ρ, 1-ρ])]

    # get connectedness matrix, depends on if subcomponents are touching
    ρ == 1/2 ? connectedness = ones(Bool,4,4) : connectedness = Matrix(IdMat(4))

    # return fractal attractor
    return Attractor(   ifs,
                        d = log(T(4))/log(1/ρ),
                        diam = sqrt(T(2)),
                        symmetries = DihedralGroup(4),
                        connectedness = connectedness)
end

function vicsek(T = Float64; ρ = 1/T(3))
        @assert 0 < ρ <= 1/3 "Contraction must satisfy ρ∈(0,1/3]"
        @assert isa(ρ, T) "Contraction ρ must be same type as first input"
        # note that automatic conversion T(ρ) will truncate ρ in some cases, so is best avoided

        # get iterated function system
        ifs = [ Similarity(ρ, [zero(T), zero(T)]),
                Similarity(ρ, [1-ρ, zero(T)]),
                Similarity(ρ, [zero(T), 1-ρ]),
                Similarity(ρ, [1-ρ, 1-ρ]),
                Similarity(ρ, [1/T(2) - ρ/2, 1/T(2) - ρ/2])] # keep centre around [0.5,0.5]

        
        # get connectedness matrix, depends on if subcomponents are touching
        if ρ == 1/3
            connectedness = Matrix{Bool}([  1 0 0 0 0;
                                            0 1 0 0 0;
                                            0 0 1 0 0;
                                            0 0 0 1 0
                                            1 1 1 1 1])
        else
            connectedness = Matrix(IdMat(5))
        end

        # return fractal attractor
        return Attractor(ifs,
                        d = log(T(5))/log(1/ρ),
                        diam = sqrt(T(2)),
                        symmetries = DihedralGroup(4),
                        connectedness = connectedness)
end

function sierpinskitriangle(T = Float64; ρ = T(1/2))
    @assert 0<ρ<=0.5 "Contraction must satisfy ρ∈(0,0.5]"
    @assert isa(ρ, T) "Contraction ρ must be same type as first input"
    # note that automatic conversion T(ρ) will truncate ρ in some cases, so is best avoided
    
    # get iterated function system
    courage = Similarity(T, [zero(T), zero(T)])
    wisdom  = Similarity(T, [1/T(2), zero(T)])
    power   = Similarity(T, [1/T(4), sqrt(T(3))/4])
    ifs = [courage, wisdom, power]

    # get connectedness matrix, depends on if subcomponents are touching
    ρ == 1/2 ? connectedness = ones(Bool,3,3) : connectedness = Matrix(IdMat(3))

    # return fractal attractor
    return Attractor(ifs,
                    d = log(T(5))/log(1/T(ρ)),
                    diam = one(T),
                    symmetries = DihedralGroup(3),
                    connectedness = connectedness)
end

# ---------------------------------------------------------------------------------------- #
# --------------------- main functions to export fractals  ------------------------------- #
# ---------------------------------------------------------------------------------------- #
# the following dict will be used to get fractals without exporting all of these functions
fractaldict = Dict(
                    :cantorset => cantorset,
                    :cantordust => cantordust,
                    :sierpinskitriangle => sierpinskitriangle,
                    :sierpinski => sierpinskitriangle,
                    :sierpinskigasket => sierpinskitriangle,
                    :vicsek => vicsek
                    )

function getfractal(T::Type, fractalname::Symbol; vargs...)
    if haskey(fractaldict, fractalname)
        return fractaldict[fractalname](T; vargs...)
    else
        str_fractalname = String(fractalname)
        error("Invalid fractal name: $str_fractalname")
    end
end

# make option to pass String instead of Symbol
getfractal(fractalname::String; vargs...) = getfractal(Symbol(fractalname); vargs...)

# define default type
getfractal(fractalname; vargs...) = getfractal(Float64, fractalname; vargs...)