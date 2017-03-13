function normalize_node{T<:AbstractFloat}(node::T,domain::Array{T,1})

  if domain[1] == domain[2]
    node = (domain[1]+domain[2])/2.0
  else
    node = 2.0*(node-domain[2])/(domain[1]-domain[2])-1.0
  end

  return node

end

function normalize_node{T<:AbstractFloat}(node::Array{T,1},domain::Array{T,1})

  for i = 1:length(node)
    node[i] = normalize_node(node[i],domain)
  end

  return node

end
