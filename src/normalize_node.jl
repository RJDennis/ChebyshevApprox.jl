function normalize_node{T<:FloatingPoint}(node::T,range::Array{T,1})

  if range[1] == range[2]
    node = (range[1]+range[2])/2
  else
    node = 2*(node-range[2])/(range[1]-range[2])-one(T)
  end

  return node

end

function normalize_node{T<:FloatingPoint}(node::Array{T,1},range::Array{T,1})

  if range[1] == range[2]
    node = ones(T,length(node))*(range[1]+range[2])/2
  else
    node = 2*(node-range[2])/(range[1]-range[2])-one(T)
  end

  return node

end

