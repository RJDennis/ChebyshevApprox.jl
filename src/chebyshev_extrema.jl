function chebyshev_extrema(n::S,domain = [1.0,-1.0]) where S <: Integer

  # Construct the nodes on the [-1.0,1.0] interval

  nodes = Array{typeof(1.0)}(undef,n)

  if n >= 2
    for i = 0:(n-1)
      nodes[i+1] = domain[2] + (1.0 - cos(i*pi/(n-1)))*(domain[1]-domain[2])/2.0
    end
    return nodes
  else
    error("chebyshev_extrema: the number of nodes must be at least 2")
  end

end
