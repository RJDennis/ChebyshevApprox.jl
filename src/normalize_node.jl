function normalize_node{T<:AbstractFloat}(node::T,domain::Array{T,1})

  normalized_nodes = similar(node)
  if domain[1] == domain[2]
    normalized_nodes = (domain[1]+domain[2])/2.0
  else
    normalized_nodes = 2.0*(node-domain[2])/(domain[1]-domain[2])-1.0
  end

  return normalized_nodes

end

function normalize_node{T<:AbstractFloat}(node::Array{T,1},domain::Array{T,1})

  normalized_nodes = similar(node)
  for i = 1:length(node)
    normalized_nodes[i] = normalize_node(node[i],domain)
  end

  return normalized_nodes

end
