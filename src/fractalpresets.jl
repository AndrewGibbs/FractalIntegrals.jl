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
    ρ == 1/2 ? connectedness = ones(Bool, 2, 2) : connectedness = Matrix(IdMat(2))

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
            connectedness = Bool[  1 0 0 0 0;
                                    0 1 0 0 0;
                                    0 0 1 0 0;
                                    0 0 0 1 0
                                    1 1 1 1 1]
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

function kochcurve(T = Float64)

    # get iterated function system:
    ifs = [ Similarity(1/T(3), [0, 0]),
            Similarity(1/T(3), [1/T(3),0], T(π)/T(3)),
            Similarity(1/T(3), [1/2, sqrt(T(3))/6], rotationmatrix2d(-T(π)/T(3))),
            Similarity(1/T(3), [2/T(3),0])]

    # get connectedness matrix:
    connectedness = [  true true false false;
                        true true true false;
                        false true true true;
                        false false true true]

    # get reflective symmetry group:
    G = [IdentityInvariantMap(T, 1), OneDimensionalInvariantMap(T(1), T(-1))]

    # return attractor
    return Attractor(ifs,
                    diam = one(T),
                    connectedness = connectedness,
                    symmetries = G);
end

function sierpinskicarpet(T = Float64; ρ = T(1/3))
    @assert 0<ρ<=1/3 "Contraction must satisfy ρ∈(0,1/3]"
    @assert isa(ρ, T) "Contraction ρ must be same type as first input"

    ifs = [ Similarity(ρ, [0,0]),
            Similarity(ρ, [ρ,0]),
            Similarity(ρ, [2ρ,0]),
            Similarity(ρ, [2ρ,ρ]),
            Similarity(ρ, [2ρ,2ρ]),
            Similarity(ρ, [ρ,2ρ]),
            Similarity(ρ, [0,2ρ]),
            Similarity(ρ, [0,ρ])
            ]

    # get connectedness matrix
    connectedness = zeros(Bool, 8, 8)
    for m=1:8
        for m_=1:8
            if m == m_
                connectedness[m,m] = true
            end
            
            # edges of subcomponents touching
            if abs(m-m_)==1 || abs(m-m_)==7
                connectedness[m,m_] = true
            end

            # corners of subcomponents touching
            if mod(m,2) == 0 && mod(m_,2) == 0 && (abs(m_-m)==2 || abs(m-m_)==6)
                connectedness[m,m_] = true
            end
        end
    end

    # return attractor
    return Attractor(ifs,
                    d = log(T(8))/log(T(3)),
                    diam = sqrt(T(2)),
                    symmetries = DihedralGroup(4),
                    connectedness = connectedness)

end

# ---------------------------------------------------------------------------------------- #
# -------------- popular 2D fractals/dragons of dimension two ---------------------------- #
# ---------------------------------------------------------------------------------------- #

function heighwaydragon(T = Float64)
    # get the iterated function system
    ifs = [ Similarity(1/sqrt(T(2)), Vector{T}([0.0, 0.0]), rotation2d(T(π)/4)),
            Similarity(1/sqrt(T(2)), Vector{T}([1.0, 0.0]), rotation2d(3T(π)/4))]

    # construct the connectedness matrix
    connectedness = Bool[  1 0 0 0;
                            0 1 1 0;
                            0 1 1 0;
                            0 0 0 1]

    # no symmetries or analytic formula for the diameter
    InvariantMeasure(ifs, d = 2, connectedness = connectedness)
end

function kochsnowflake(T = Float64; rescale = sqrt(T(3))/3)
    
    # construct the iterated function system for the Koch
    ifs = [ Similarity(sqrt(1/T(3)), rescale*[0, 0], T(π)/6),
            Similarity(1/3, rescale*[0, 2/3]),
            Similarity(1/3, rescale*[-1/sqrt(T(3)), 1/3]),
            Similarity(1/3, rescale*[-1/sqrt(T(3)), -1/3]),
            Similarity(1/3, rescale*[0,-2/3]),
            Similarity(1/3, rescale*[1/sqrt(T(3)), -1/3]),
            Similarity(1/3, rescale*[1/sqrt(T(3)), 1/3])
            ]

    # create connectedness matrix
    connectedness = zeros(Bool, 7^2, 7^2)
    connectedness_level_one = Bool[1  1  1  1  1  1  1;
                                    1  1  1  0  0  0  1;
                                    1  1  1  1  0  0  0;
                                    1  0  1  1  1  0  0;
                                    1  0  0  1  1  1  0;
                                    1  0  0  0  1  1  1;
                                    1  1  0  0  0  1  1]
    # make diagonal blocks
    for m in 1:7
        diag_block_range = ((m-1)*7+1):(7*m)
        connectedness[diag_block_range, diag_block_range] .= connectedness_level_one
    end

    # now list indices which are also true
    Λ = [([2,5],[1,2]),
        ([2,4],[1,2]),
        ([3,7],[1,2]),
        ([3,6],[1,2]),
        ([3,7],[2,4]), # all smalls around [1,2] done
        ([3,6],[1,3]),
        ([3,5],[1,3]),
        ([4,2],[1,3]),
        ([4,7],[1,3]),
        ([3,5],[4,2]), # all around [1,3] done
        ([4,7],[1,4]),
        ([4,6],[1,4]),
        ([5,3],[1,4]),
        ([5,2],[1,4]),
        ([4,6],[5,3]), # all around [1,4] done
        ([1,5],[5,2]),
        ([1,5],[5,7]),
        ([1,5],[6,4]),
        ([1,5],[6,3]),
        ([5,7],[6,4]),# all around [1,5] done
        ([1,6],[6,3]),
        ([1,6],[6,2]),
        ([1,6],[7,5]),
        ([1,6],[7,4]),
        ([7,5],[6,2]),# all around [1,6] done
        ([1,7],[7,4]),
        ([1,7],[7,3]),
        ([1,7],[2,6]),
        ([1,7],[2,5]),
        ([7,3],[2,6]),# all around [1,7] done
        ([1,1],[2,5]),
        ([1,1],[3,6]),
        ([1,1],[4,7]),
        ([1,1],[5,2]),
        ([1,1],[6,3]),
        ([1,1],[7,4]),# all the new level 2 singularities done
        ([2,1],[1,2]),
        ([2,1],[1,7]),
        ([3,1],[1,2]),
        ([3,1],[1,3]),
        ([4,1],[1,3]),
        ([4,1],[1,4]),
        ([5,1],[1,4]),
        ([5,1],[1,5]),
        ([6,1],[1,5]),
        ([6,1],[1,6]),
        ([7,1],[1,6]),
        ([7,1],[1,7])]#
    
    for (m,m_) ∈ Λ
        connectedness[(m[1]-1)*7+m[2], (m_[1]-1)*7+m_[2]] = true
        connectedness[(m_[1]-1)*7+m_[2], (m[1]-1)*7+m[2]] = true
    end
    area =  2*3*sqrt(T(3))/5 * rescale^2
    
    InvariantMeasure(ifs, d = 2, connectedness = connectedness, selfmeasure = area)
end

# ---------------------------------------------------------------------------------------- #
# --------------------- main functions to export fractals  ------------------------------- #
# ---------------------------------------------------------------------------------------- #
# the following dict will be used to get fractals without exporting all of these functions
fractaldict = Dict( 
                    :cantorset => cantorset,
                    # measure zero Γ ⊂ ℝ² attractors
                    :cantordust => cantordust,
                    :sierpinskitriangle => sierpinskitriangle,
                    :sierpinski => sierpinskitriangle, # alt name for sierpinskitriangle
                    :sierpinskigasket => sierpinskitriangle, # alt name for sierpinskitriangle
                    :vicsek => vicsek,
                    # Lebesgue measurable d=2 attractors / dragons
                    :heighwaydragon => heighwaydragon,
                    :kochsnowflake => kochsnowflake,
                    :kochflake => :kochsnowflake, # alt name for kochsnowflake
                    :kochcurve => kochcurve,
                    :sierpinskicarpet => sierpinskicarpet
                    :carpet => sierpinskicarpet # alt name for sierpinskicarpet
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