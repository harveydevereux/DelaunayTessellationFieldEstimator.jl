module DelaunayTessellationFieldEstimator

# Implements the DTFE algoritm as found in
# Schaap, W.E., 2007. The delaunay tessellation field estimator. Ph. D. Thesis.
# https://www.rug.nl/research/portal/files/2816076/c2.pdf
# Author Harvey Devereux (h.devereux@warwick.ac.uk)
# License: MIT
# Any collaboration welcomed

# university wesite https://warwick.ac.uk/fac/sci/mathsys/people/students/2017intake/devereux/
# personal website https://harveydevereux.github.io/

# Works now but is slower than python (possibly due to pyimport overhead?)

# TODO See if GridInterpolations.jl SimplexGrid is faster than Py-importing
# TODO 3D -> N-D
# TODO Gif visualisation

using VoronoiDelaunay
import VoronoiDelaunay: delaunayedges
using GeometricalPredicates
using NearestNeighbors
import NearestNeighbors: knn, KDTree

# for the simplex based interpolation
using PyCall
@pyimport scipy.interpolate as si

export
DTFE,
Coordinates, GeomToDataAreaScalingFactor, BuildDelaunay!, FindTriangles,
ContiguousVoronoiCellArea, DelaunayTessellation2D, DelaunayTessellation,
Point, Point2D, getplotxyd, delaunayedges, UniquePoints, GeomToDataAreaScalingFactor,
DelaunayVertexDensity, DelaunayVertexGeometryDensity, Point2DToArray, DTFEMap,
Grid, meshgrid

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

   Here M,m and the new min and max (assumed the same for all dimensions)
   and min,max are the dimensions maximum and minimum.

   If the scaling is x* = [(x-min)/(max-min)](M-m)+m by writing
   out (a₁-b₁)(a₂-b₂) in the form of x* you\'ll find
   (a₁-b₁)(a₂-b₂)*SF⁻¹
   where SF  = ∏ᵢ(M-m)/(maxᵢ-minᵢ)
"""
function GeomToDataAreaScalingFactor(C::Coordinates)
    prod = 1
    D = C.DataLims
    G = C.GeomLims
    for d in 1:size(D,1)
        prod *= (maximum(D[d,:])-minimum(D[d,:]))/(maximum(G[d,:])-minimum(G[d,:]))
    end
    return prod
end

"""Wrapper to build a Delaunay tessellation from an
   array of coordinates

   returns the points, the tess object is modified so
   no return necessary
"""
# TODO same function but with the coordinates object?
function BuildDelaunay!{A <: Array{<: Real,2}}(Coords::A)
    # doesnt't seem to exist
    # sizehint(tess,size(Coords,1))
    # This because if you try to push
    # to an old tess the function hangs
    # in an infinite loop
    # there is an issue for this in
    # VoronoiDelaunay.jl
    Tess = DelaunayTessellation()
    Points = Point2D[Point(Coords[i,1],Coords[i,2]) for i in 1:size(Coords,1)]
    push!(Tess, Points)
    return Tess,Points
end

"""Usefull for quikcly getting th triangles vertices
all at once, returns a T <: Point2D
"""
function UnpackVertices{T <: VoronoiDelaunay.DelaunayTriangle}(Triangle::T)
    return geta(Triangle), getb(Triangle), getc(Triangle)
end

"""Usefull for quikcly getting th triangles vertices
all at once, returns a T <: Array
"""
function UnpackVerticesAsArray{T <: VoronoiDelaunay.DelaunayTriangle}(Triangle::T)
    a = [getx(geta(Triangle)),gety(geta(Triangle))]
    b = [getx(getb(Triangle)),gety(getb(Triangle))]
    c = [getx(getc(Triangle)),gety(getc(Triangle))]
    return a,b,c
end

"""Finds all the delaunay triangles with the common vertex of \'point\'
   returns an array of these triangles.
"""
function FindTriangles{T <: VoronoiDelaunay.DelaunayTessellation2D,
                       P <: GeometricalPredicates.Point2D}(tess::T,point::P)
    Triangles = Array{VoronoiDelaunay.DelaunayTriangle,1}()
    for triangle in tess
        a, b, c = UnpackVertices(triangle)
        if point == a || point == b || point == c
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

"""Finds the set of unique points within a contigous
   Voronoi cell

   returns the array of unique points
"""
function UniquePoints{T <: VoronoiDelaunay.DelaunayTriangle}(Triangles::Array{T,1})
    Uniques = []
    for triangle in Triangles
        for v in UnpackVertices(triangle)
            if !(v in Uniques)
                push!(Uniques,v)
            end
        end
    end
    return Uniques
end

"""Finds the density in the geometry space (no scaling)
   for a single Delaunay vertex point
"""
function DelaunayVertexGeometryDensity{T <: VoronoiDelaunay.DelaunayTessellation2D,
                                       P <: GeometricalPredicates.Point2D}(tess::T,point::P,mass::Float64=1.0)
    Triangles = FindTriangles(tess,point)
    A = ContiguousVoronoiCellArea(Triangles)
    ρ = mass*length(UniquePoints(Triangles))/A
    return ρ
end

"""Wraps around DelaunayVertexGeometryDensity to calculate
   the density in data space for every delaunay vertex point

   Masses can be set to some field value other than mass. The
   default is a uniform unit mass.

   Returns an array of densities scaled to the data space
"""
function DelaunayVertexDensity{T <: VoronoiDelaunay.DelaunayTessellation2D}(tess::T,C::Coordinates,Masses=ones(length(C.GeometryCoordinates)))
    G = C.GeometryCoordinates
    ρ = zeros(size(G,1))
    for i in 1:size(G,1)
        P = Point2D(G[i,1],G[i,2])
        ρ[i] = DelaunayVertexGeometryDensity(tess,P,Masses[i])
    end
    return ρ./GeomToDataAreaScalingFactor(C)
end

meshgrid(x,y) = (repmat(x',length(y),1),repmat(y,1,length(x)))

"""Computes the interpolated DTFE grid from the given points
   with resolution given by Size
"""
function DTFE(Points,Size=500)
    C = Coordinates(Points)
    tess,P = BuildDelaunay!(C.GeometryCoordinates)
    x,y = meshgrid(linspace(C.DataLims[1,:]...,Size),linspace(C.DataLims[2,:]...,Size))

    D = DelaunayTessellationFieldEstimator.DelaunayVertexDensity(tess,C)
    grid = si.griddata(C.DataCoordinates,D,(x,y),method="linear",fill_value=0)
    return grid
end

end # module DelaunayTessellationFieldEstimator
