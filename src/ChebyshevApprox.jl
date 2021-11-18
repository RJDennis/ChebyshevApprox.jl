module ChebyshevApprox

using LinearAlgebra
using ThreadPools

include("structures.jl")
include("chebyshev_nodes.jl")
include("chebyshev_extrema.jl")
include("chebyshev_extended.jl")
include("vertesi_nodes.jl")
include("chebyshev_quadrature.jl")
include("normalize_node.jl")
include("chebyshev_polynomial.jl")
include("chebyshev_weights.jl")
include("chebyshev_evaluate.jl")
include("chebyshev_derivative.jl")
include("clenshaw_evaluate_generated.jl")

export ChebPoly,
       ChebInterpRoots,
       ChebInterpExtrema,
       ChebInterpExtended,
       ChebInterpVertesi

export chebyshev_nodes,
       chebyshev_extrema,
       chebyshev_extended,
       vertesi_nodes,
       chebyshev_quadrature,
       chebyshev_polynomial,
       chebyshev_weights,
       chebyshev_weights_extrema,
       chebyshev_weights_extended,
       chebyshev_weights_vertesi,
       chebyshev_weights_threaded,
       chebyshev_weights_extrema_threaded,
       chebyshev_weights_extended_threaded,
       chebyshev_weights_vertesi_threaded,
       chebyshev_evaluate,
       cheb_interp,
       chebyshev_derivative,
       chebyshev_gradient,
       clenshaw_evaluate

end