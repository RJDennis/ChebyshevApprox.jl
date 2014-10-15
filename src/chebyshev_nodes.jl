function chebyshev_nodes{S<:Integer}(n::S)

  nodes = Array(typeof(1.0),n)

  for i = 1:n

    nodes[i] = -cos((2*i-1.0)*pi/(2*n))

  end

  return nodes

end

function chebyshev_nodes{T<:FloatingPoint,S<:Integer}(n::S,range::Array{T,1})

  nodes = Array(T,n)

  for i = 1:n

    nodes[i] = -cos((2*i-1.0)*pi/(2*n))

  end

  nodes = range[2] .+ (1.0 .+ nodes)*(range[1]-range[2])/2

  return nodes

end
