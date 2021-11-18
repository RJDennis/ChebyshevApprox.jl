function chebyshev_quadrature(node_type::Function,n::S,domain = [1.0,-1.0]) where {S <: Integer}

  if n <= 0

    error("The number of nodes must be positive.")
  
  end
      
  nodes   = zeros(n)
  weights = zeros(n)

  if node_type == chebyshev_nodes # weighting function is (1.0 .- nodes.^2).^(-1/2) -- (first kind)

    for i = 1:n
      nodes[i] = (domain[1]+domain[2])/2.0 - cos((i-0.5)*pi/n)*(domain[1]-domain[2])/2.0
      weights[i] = pi/n
    end

    return nodes, weights

  elseif node_type == chebyshev_extrema # weighting function is (1.0 .- nodes.^2).^(1/2) -- (second kind)

    for i = 1:n
      nodes[i] = (domain[1]+domain[2])/2.0 - cos((i-1)*pi/(n-1))*(domain[1]-domain[2])/2.0
      weights[i] = (pi/(n-1))*sin(((i-1)/(n-1))*pi)^2
    end

    return nodes, weights
    
  else
    
    error("This node type is not supported")

  end

end
    