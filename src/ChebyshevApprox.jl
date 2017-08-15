module ChebyshevApprox

include("chebyshev_nodes.jl")
include("normalize_node.jl")
include("chebyshev_polynomial.jl")
#include("chebyshev_weights.jl")
include("chebyshev_weights_generated.jl")
#include("chebyshev_evaluate.jl")
include("chebyshev_evaluate_generated.jl")
include("clenshaw_evaluate.jl")
include("chebyshev_polynomial_derivative.jl")
include("chebyshev_derivative.jl")

export chebyshev_nodes,
       chebyshev_polynomial,
       chebyshev_weights,
       chebyshev_evaluate,
       clenshaw_evaluate,
       chebyshev_derivative

end
