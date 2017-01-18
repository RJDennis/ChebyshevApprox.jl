function normalize_node{T<:AbstractFloat}(node::T,range::Array{T,1})

  if range[1] == range[2]
    node = (range[1]+range[2])/2.0
  else
    node = 2.0*(node-range[2])/(range[1]-range[2])-1.0
  end

  return node

end

function normalize_node{T<:AbstractFloat}(node::Array{T,1},range::Array{T,1})

  for i = 1:length(node)
    node[i] = normalize_node(node[i],range)
  end

  return node

end
