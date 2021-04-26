function chebyshev_nodes(n::S,domain = [1.0,-1.0]) where {S <: Integer}

  if n <= 0
    error("The number of nodes must be positive.")
  end

  nodes = zeros(n)

  if isodd(n)
    nodes[Int((n-1)/2)+1] = (domain[1]+domain[2])/2.0
  end

  for i = 1:div(n,2)
    x = -cos((i-0.5)*pi/n)*(domain[1]-domain[2])/2.0
    nodes[i]       = (domain[1]+domain[2])/2.0 + x
    nodes[end-i+1] = (domain[1]+domain[2])/2.0 - x
  end

  return nodes

end
