module ChebyshevApprox

using LinearAlgebra
using ThreadPools

include("chebyshev_approx_functions")
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