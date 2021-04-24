module ChebyshevApprox

using ThreadPools

include("structures.jl")
include("chebyshev_nodes.jl")
include("chebyshev_extrema.jl")
include("normalize_node.jl")
include("chebyshev_polynomial.jl")
include("chebyshev_weights.jl")
include("chebyshev_evaluate.jl")
include("clenshaw_evaluate_generated.jl")
include("chebyshev_derivative.jl")

export ChebPolyTensor,
       ChebPolyComplete
       ChebInterpTensor,
       ChebInterpComplete

export chebyshev_nodes,
       chebyshev_extrema,
       chebyshev_polynomial,
       chebyshev_weights,
       chebyshev_evaluate,
       cheb_interp,
       chebyshev_derivative,
       chebyshev_gradient,
       clenshaw_evaluate

end
