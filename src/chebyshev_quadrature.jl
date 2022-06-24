function chebyshev_quadrature(node_type::Function,n::S,domain::Array{T,1} = [1.0,-1.0]) where {S<:Integer,T<:AbstractFloat}

  if n <= 0
      error("The number of nodes must be positive.")
  end
      
  nodes   = Array{T,1}(undef,n)
  weights = Array{T,1}(undef,n)

  if node_type == chebyshev_nodes # weighting function is (1.0 .- nodes.^2).^(-1/2) -- (first kind)

      nodes .= chebyshev_nodes(n,domain)
      for i = 1:n
          weights[i] = pi/n
      end

      return nodes, weights

  elseif node_type == chebyshev_extrema # weighting function is (1.0 .- nodes.^2).^(1/2) -- (second kind)

      nodes .= chebyshev_extrema(n,domain)
      for i = 1:n
          weights[i] = (pi/(n-1))*sin(((i-1)/(n-1))*pi)^2
      end

      return nodes, weights
    
  else
    
      error("This node type is not supported")

  end

end