function chebyshev_nodes{S<:Integer}(n::S)

  nodes = Array{Float64}(n)

  for i = 1:n
    nodes[i] = -cos((2.0*i-1.0)*pi/(2.0*n))
  end

  return nodes

end

function chebyshev_nodes{T<:AbstractFloat,S<:Integer}(n::S,domain::Array{T,1})

  nodes = Array{T}(n)

  for i = 1:n
    nodes[i] = domain[2] + ( 1.0 - cos((2.0*i-1.0)*pi/(2.0*n)) )*(domain[1]-domain[2])/2.0
  end

  return nodes

end
