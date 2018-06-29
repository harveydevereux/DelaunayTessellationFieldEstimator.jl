include("../src/DelaunayTessellationFieldEstimator.jl")
using DelaunayTessellationFieldEstimator
using Plots

Coords = rand(0.0:500.0,10,2)

tess = DelaunayTessellation()

C = Coordinates(Coords)

P = BuildDelaunay!(tess,C.GeometryCoordinates)
x,y = VoronoiDelaunay.getplotxy(delaunayedges(tess))
plot(x,y)
scatter!(C.GeometryCoordinates[:,1],C.GeometryCoordinates[:,2])

t = FindTriangles(tess,P[1])
A = ContiguousVoronoiCellArea(t)

# yet to get a triangle counter
(A*GeomToDataAreaScalingFactor(C))

u=UniquePoints(t)
DelaunayVertexDensity(tess,P[1])

length(u)/(A*GeomToDataAreaScalingFactor(C))

DelaunayVertexDensity(tess,C)
