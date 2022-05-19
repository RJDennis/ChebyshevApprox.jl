function normalize_node(node::R,domain::Array{T,1}) where {R<:Number,T<:AbstractFloat}

  if domain[1] == domain[2]
    norm_node = zero(T)
    return norm_node
  else
    norm_node = 2.0*(node-domain[2])/(domain[1]-domain[2])-1.0
    return norm_node
  end
  
end
  
function normalize_node(node::AbstractArray{R,1},domain::Array{T,1}) where {R<:Number,T<:AbstractFloat}
  
  norm_nodes = similar(node)
  for i in eachindex(node)
    norm_nodes[i] = normalize_node(node[i],domain)
  end
  
  return norm_nodes
  
end
  
function normalize_node!(node::AbstractArray{R,1},domain::Array{T,1}) where {R<:Number,T<:AbstractFloat}
  
  for i in eachindex(node)
    node[i] = normalize_node(node[i],domain)
  end
  
end