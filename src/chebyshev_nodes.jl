function chebyshev_nodes(n::S,domain = [1.0,-1.0]) where {S <: Integer}

  nodes = Array{typeof(1.0)}(undef,n)

  for i = 1:n
    nodes[i] = domain[2] + ( 1.0 - cos((2.0*i-1.0)*pi/(2.0*n)) )*(domain[1]-domain[2])/2.0
  end

  return nodes

end
