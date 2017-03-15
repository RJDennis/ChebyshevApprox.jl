function normalize_node{T<:AbstractFloat}(node::T,domain::Array{T,1})

  if domain[1] == domain[2]
    norm_node = (domain[1]+domain[2])/2.0
    return norm_node
  else
    norm_node = 2.0*(node-domain[2])/(domain[1]-domain[2])-1.0
    return norm_node
  end

end

function normalize_node{T<:AbstractFloat}(node::Array{T,1},domain::Array{T,1})

  norm_nodes = similar(node)
  for i = 1:length(node)
    norm_nodes[i] = normalize_node(node[i],domain)
  end

  return norm_nodes

end
