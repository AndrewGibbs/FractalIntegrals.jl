# ---------------------------------------------------------------------------------------- #
# -------------- popular 2D fractals/dragons of dimension two ---------------------------- #
# ---------------------------------------------------------------------------------------- #

function heighwaydragon(T = Float64)
    # get the iterated function system
    ifs = [ Similarity(1/sqrt(T(2)), Vector{T}([0.0, 0.0]), rotationmatrix2d(T(π)/4)),
            Similarity(1/sqrt(T(2)), Vector{T}([1.0, 0.0]), rotationmatrix2d(3T(π)/4))]

    # construct the connectedness matrix
    connectedness = Bool[  1 0 0 0;
                            0 1 1 0;
                            0 1 1 0;
                            0 0 0 1]

    # no symmetries or analytic formula for the diameter
    Attractor(ifs, d = 2, connectedness = connectedness)
end

function kochsnowflake(T = Float64; rescale = sqrt(T(3))/3)
    
    # construct the iterated function system for the Koch
    ifs = [ Similarity(sqrt(1/T(3)), rescale*[0, 0], rotationmatrix2d(T(π)/6)),
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
    
    Attractor(  ifs,
                d = 2,
                connectedness = connectedness,
                symmetries = DihedralGroup(T, 6))
end