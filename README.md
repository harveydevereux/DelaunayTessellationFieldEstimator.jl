# DelaunayTessellationFieldEstimator.jl

### What is it?

This package is an implementation the DTFE algorithm [(Schaap, 2007 PhD)](https://www.rug.nl/research/portal/files/2816076/c2.pdf)
for calculating an estimate of a particle density field (can be used for other fields too if the 'mass' is some another physical 
attribute).

### Dependent on

- [GeometriclPredicates.jl](https://github.com/JuliaGeometry/GeometricalPredicates.jl)
- [VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl)


### Current limitations

- Currently 2D since [VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl) is and everything
  is based on this.

- 'Due to numerical restrictions the point coordinates must be within min_coord <= x <= max_coord where min_coord=1.0+eps(Float64) and max_coord=2.0-2eps(Float64). Note this is a bit different than what is required by the GeometricalPredicates package.' - Quote
[VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl)

### TODO 

- Nice and easy plotting + movie making :)
- 3D, would need [VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl) updated.

### Example of what this does

[Notebook](https://github.com/harveydevereux/DelaunayTessellationFieldEstimator.jl/blob/master/src/Example.ipynb)
