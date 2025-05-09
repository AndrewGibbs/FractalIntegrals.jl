struct Polygon{T}
    nodes::Vector{SVector{2,T}}
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
    return Polygon([SVector{2}([tofloat(node, mag, pres)...]) for node in polys[1]])
end

# could generalise below fn to hull of fixed points for homogenous attractors
# also, this should be for a subtype of attractors only, open sets, or something
function get_bounding_square(Γ::AbstractAttractor, r = 0.6 + 1e-1 * rand())
    centre = get_boundingball_centre(Γ)
    hsl = diam(Γ) # half side length
    return Polygon(
        SVector{
            2,
        }.([
            centre .+ r * [-hsl, -hsl],
            centre .+ r * [-hsl, hsl],
            centre .+ r * [hsl, hsl],
            centre .+ r * [hsl, -hsl],
        ]),
    )
end

# define Similarity map applied to polygon
(s::Similarity)(p::Polygon) = Polygon(s.(p.nodes))

function polygon_approx(Γ::AbstractAttractor; levels = 2)
    union_poly = get_bounding_square(Γ)

    for j = 1:levels

        # apply each ifs map to 
        mapped_polys = [s(union_poly) for s in Γ.ifs]

        # initialise union
        union_poly = mapped_polys[1]

        # unify with all other mapped polygons
        for poly in mapped_polys[2:end]
            union_poly = union_poly ∪ poly # take union of union with next polygon
        end
    end

    return union_poly
end
