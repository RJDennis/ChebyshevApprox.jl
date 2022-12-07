module ChebyshevApprox

using LinearAlgebra
using ThreadPools

include("chebyshev_approx_functions.jl")

export ChebRoots,
       ChebExtrema,
       ChebExtended,
       VertesiNodes,
       LegendreNodes,
       Grid,
       ChebPoly,
       CApproxPlan,
       CApproxPlanPoly


export chebyshev_nodes,
       chebyshev_extrema,
       chebyshev_extended,
       vertesi_nodes,
       legendre_nodes,
       nodes,
       chebyshev_polynomial,
       chebyshev_polynomial_deriv,
       chebyshev_polynomial_sec_deriv,
       chebyshev_weights,
       chebyshev_weights_extrema,
       chebyshev_weights_extended,
       chebyshev_weights_vertesi,
       chebyshev_weights_legendre,
       chebyshev_weights_threaded,
       chebyshev_weights_extrema_threaded,
       chebyshev_weights_extended_threaded,
       chebyshev_weights_vertesi_threaded,
       chebyshev_weights_legendre_threaded,
       chebyshev_evaluate,
       chebyshev_interp,
       chebyshev_interp_threaded,
       chebyshev_derivative,
       chebyshev_gradient,
       chebyshev_gradient_threaded,
       chebyshev_hessian,
       chebyshev_hessian_threaded

end