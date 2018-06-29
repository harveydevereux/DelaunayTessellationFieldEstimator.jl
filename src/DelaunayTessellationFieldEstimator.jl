module DelaunayTessellationFieldEstimator

# Implements the DTFE algoritm as found in
# Schaap, W.E., 2007. The delaunay tessellation field estimator. Ph. D. Thesis.
# https://www.rug.nl/research/portal/files/2816076/c2.pdf
# Author Harvey Devereux (h.devereux@warwick.ac.uk)
# License: MIT
# Any collaboration welcomed

# university wesite https://warwick.ac.uk/fac/sci/mathsys/people/students/2017intake/devereux/
# personal website https://harveydevereux.github.io/

# still very much a work in progress!

# TODO triangle counter and density
# TODO Interpolation for inbetween points
# TODO Wrap it all in the DTFEMap function
# TODO Parallel support?
# TODO Gif visualisation

using VoronoiDelaunay
import VoronoiDelaunay: delaunayedges 
using GeometricalPredicates

export
Coordinates, GeomToDataAreaScalingFactor, BuildDelaunay!, FindTriangles,
ContiguousVoronoiCellArea, DelaunayTessellation2D, DelaunayTessellation,
Point,Point2D,getplotxyd,delaunayedges

const min_lim = VoronoiDelaunay.min_coord
const max_lim = VoronoiDelaunay.max_coord

"""Holds the original data and the rescaled geometric data to
    comform the VoronoiDelaunay\'s range requirements

    NB: GeometricalPredicates requires points ∈ [1,2] but
        VoronoiDelaunay requires points in [1+eps(float),2-2eps(float)]
"""
mutable struct Coordinates{T <: Array{<:Real,2}}
    DataCoordinates::T
    GeometryCoordinates::T
    DataLims::T
    GeomLims::T
end

"""Constructor to deal which automatically calculates the geometry scale
   data and holds the limits for future scaling.
"""
function Coordinates{T <: Array{<:Real,2}}(Points::T)
    Glim = [min_lim max_lim; min_lim max_lim]
    Dlim = [minimum(Points[:,1]) maximum(Points[:,1]); minimum(Points[:,2]) maximum(Points[:,2])]
    GPoints = Scaling(Points, Glim)
    return Coordinates(Points,GPoints,Dlim,Glim)
end

"""Min max scaling to the VoronoiDelaunay range requirement
"""
function Scaling{T <: Array{<:Real,2}}(coords::T,lim::T)
    Rescale = zeros(size(coords))
    for d in 1:2
        min = minimum(coords[:,d])
        max = maximum(coords[:,d])
        a = minimum(lim[d,:])
        b = maximum(lim[d,:])
        Rescale[:,d] = [ ((coords[i,d] - min)/(max-min))*(b-a)+a for i in 1:size(coords,1)]
    end
    return Rescale
end

"""Returns the scaling factor required to calculate the area in the original
   point space.

   If the scaling is x* = [(x-min)/(max-min)](M-m)+m by writing
   out (a₁-b₁)(a₂-b₂) in the form of x* you\'ll find
   (a₁-b₁)(a₂-b₂)*SF⁻¹
   where SF  = ∏ᵢ(M-m)/(maxᵢ-minᵢ)
"""
function GeomToDataAreaScalingFactor(C::Coordinates)
    prod = 0
    D = C.DataLims
    G = C.GeomLims
    for d in 1:size(D,1)
        prod+=((maximum(D[d,:])-minimum(D[d,:]))/(maximum(G[d,:])-minimum(G[d,:])))^2
    end
    return prod
end

"""Wrapper to build a Delaunay tessellation from an
   array of coordinates

   returns the points, the tess object is modified so
   no return necessary
"""
# TODO same function but with the coordinates object?
function BuildDelaunay!{T <: VoronoiDelaunay.DelaunayTessellation2D,
                       A <: Array{<: Real,2}}(Tess::T,Coords::A)
    # doesnt't seem to exist
    # sizehint(tess,size(Coords,1))
    Points = Point2D[Point(Coords[i,1],Coords[i,2]) for i in 1:size(Coords,1)]
    push!(Tess, Points)
    return Points
end

"""Finds all the delaunay triangles with the common vertex of \'point\'
   returns an array of these triangles.
"""
function FindTriangles{T <: VoronoiDelaunay.DelaunayTessellation2D,
                       P <: GeometricalPredicates.Point2D}(tess::T,point::P)
    Triangles = Array{VoronoiDelaunay.DelaunayTriangle,1}()
    for triangle in tess
        if point == triangle._a || point == triangle._b || point == triangle._c
            push!(Triangles,triangle)
        end
    end
    return Triangles
end

"""The contiguous voronoi cell at a point i (Schaap, 2007 PhD) is the union of
   all delaunay triangles with the coomon vertex i

   returns the Geometry scale area of the union of the triangles
   in \'Triangles\', use FindTriangles to get these for a point.
"""
function ContiguousVoronoiCellArea{T <: VoronoiDelaunay.DelaunayTriangle}(Triangles::Array{T,1})
    Area=0
    for triangle in Triangles
        Area+=GeometricalPredicates.area(triangle)
    end
    return Area
end

end # module DelaunayTessellationFieldEstimator
