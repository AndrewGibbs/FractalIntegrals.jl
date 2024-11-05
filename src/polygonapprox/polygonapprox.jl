struct Polygon{T}
    nodes::Vector{SVector{2, T}}
end

function ∪(poly1::Polygon, poly2::Polygon, pres = 8)
    mag = 1 # need to adjust this later
    path1 = [IntPoint(node[1], node[2], mag, pres) for node in poly1.nodes]
    path2 = [IntPoint(node[1], node[2], mag, pres) for node in poly2.nodes]
    c = Clip()
    add_path!(c, path1, PolyTypeSubject, true)
    add_path!(c, path2, PolyTypeClip, true)
    result, polys = execute(c, ClipTypeUnion, PolyFillTypeEvenOdd, PolyFillTypeEvenOdd)
    @assert result "Clipper unification failed"
    @assert length(polys) == 1 "Clipper union has returned two disjoint polygons"
    # poly_float = tofloat.(polys[1], mag, pres)
    return Polygon(SVector{2}([tofloat(node, mag, pres)...]) for node in polys[1])
end

# could generalise below fn to hull of fixed points for homogenous attractors
# also, this should be for a subtype of attractors only, open sets, or something
function get_bounding_square(Γ::AbstractAttractor,
                            r = 1 + 1e-4*rand() # random stretch
                            )
    centre = get_boundingball_centre(Γ)
    hsl = diam(Γ) # half side length
    return Polygon([centre .+ r*[-hsl[1], -hsl[2]],
                    centre .+ r*[-hsl[1], hsl[2]],
                    centre .+ r*[hsl[1], hsl[2]],
                    centre .+ r*[hsl[1], -hsl[2]]])
end

function polygon_approx(Γ::AbstractAttractor,
                        levels = 2 #prefractal levels required
                        )
    union_poly = get_bounding_square(Γ)

    for j in 1:levels

        # apply each ifs map to 
        mapped_polys = [s(union_poly) for s in Γ.ifs]

        # initialise union
        union_poly = mapped_polys[1]

        # unify with all other mapped polygons
        for poly in union_poly[2:end]
            union_poly = union_poly ∪ poly # take union of union with next polygon
        end
    end

    return union_poly
end

Shape(p::Polygon) = Shape([x[1] for x in p.nodes], [x[2] for x in p.nodes])