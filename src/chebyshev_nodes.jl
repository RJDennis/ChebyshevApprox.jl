function chebyshev_nodes(n::S,domain::Array{T,1} = [1.0,-1.0]) where {S<:Integer,T<:AbstractFloat}

  if n <= 0
      error("The number of nodes must be positive.")
  end

  nodes = Array{T,1}(undef,n)

  if isodd(n)
      nodes[div(n-1,2)+1] = (domain[1]+domain[2])*0.5
  end

  for i = 1:div(n,2)
      x = -cos((i-0.5)*pi/n)*(domain[1]-domain[2])*0.5
      nodes[i]       = (domain[1]+domain[2])*0.5 + x
      nodes[end-i+1] = (domain[1]+domain[2])*0.5 - x
  end

  return nodes

end
