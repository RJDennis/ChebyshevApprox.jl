module ChebyshevApprox

include("chebyshev_nodes.jl")
include("normalize_node.jl")
include("chebyshev_polynomial.jl")
include("chebyshev_weights_generated.jl")
include("chebyshev_evaluate_generated.jl")
include("clenshaw_evaluate_generated.jl")
include("chebyshev_derivative_generated.jl")


export chebyshev_nodes,
       chebyshev_polynomial,
       chebyshev_weights,
       chebyshev_evaluate,
       clenshaw_evaluate,
       chebyshev_derivative

end
