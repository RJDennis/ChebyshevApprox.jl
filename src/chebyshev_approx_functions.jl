abstract type ChebInterp end

struct ChebPoly{T<:AbstractFloat, S<:Integer, N}

  weights::Array{T,N}
  order::Union{S,Array{S,1}}
  domain::Union{Array{T,1},Array{T,2}}

end

struct ChebInterpRoots{T<:AbstractFloat, S<:Integer, N} <: ChebInterp

  data::Array{T,N}
  nodes::NTuple{N,Array{T,1}}
  order::Union{S,Array{S,1}}
  domain::Union{Array{T,1},Array{T,2}}

end

struct ChebInterpExtrema{T<:AbstractFloat, S<:Integer, N} <: ChebInterp

  data::Array{T,N}
  nodes::NTuple{N,Array{T,1}}
  order::Union{S,Array{S,1}}
  domain::Union{Array{T,1},Array{T,2}}

end

struct ChebInterpExtended{T<:AbstractFloat, S<:Integer, N} <: ChebInterp

  data::Array{T,N}
  nodes::NTuple{N,Array{T,1}}
  order::Union{S,Array{S,1}}
  domain::Union{Array{T,1},Array{T,2}}

end

struct ChebInterpVertesi{T<:AbstractFloat, S<:Integer, N} <: ChebInterp

  data::Array{T,N}
  nodes::NTuple{N,Array{T,1}}
  order::Union{S,Array{S,1}}
  domain::Union{Array{T,1},Array{T,2}}

end

function chebyshev_nodes(n::S,domain = [1.0,-1.0]) where {S <: Integer}

  if n <= 0
    error("The number of nodes must be positive.")
  end
  
  nodes = zeros(n)
  
  if isodd(n)
    nodes[div(n-1,2)+1] = (domain[1]+domain[2])/2.0
  end
  
  for i = 1:div(n,2)
    x = -cos((i-0.5)*pi/n)*(domain[1]-domain[2])/2.0
    nodes[i]       = (domain[1]+domain[2])/2.0 + x
    nodes[end-i+1] = (domain[1]+domain[2])/2.0 - x
  end
  
  return nodes
  
end

function chebyshev_extrema(n::S,domain = [1.0,-1.0]) where S <: Integer

  if n <= 0
    error("The number of nodes must be positive.")
  end
  
  nodes = zeros(n)
  
  if isodd(n)
    nodes[div(n-1,2)+1] = (domain[1]+domain[2])/2.0
  end
  
  for i = 1:div(n,2)
    x = -cos((i-1)*pi/(n-1))*(domain[1]-domain[2])/2.0
    nodes[i]       = (domain[1]+domain[2])/2.0  + x
    nodes[end-i+1] = (domain[1]+domain[2])/2.0  - x
  end
  
  return nodes
  
end

function chebyshev_extended(n::S,domain = [1.0,-1.0]) where {S <: Integer}

  if n <= 0
    error("The number of nodes must be positive.")
  end
  
  nodes = zeros(n)
  
  if isodd(n)
    nodes[div(n-1,2)+1] = 0.0
  end
  
  for i = 1:div(n,2)
    x = -cos((i-0.5)*pi/n)
    nodes[i]       = x
    nodes[end-i+1] = -x
  end
  
  return (domain[1]+domain[2])*0.5 .+ nodes*((domain[1]-domain[2])*0.5)/cos(π/(2N))
  
end

function vertesi_nodes(n::S,domain = [1.0,-1.0]) where {S <: Integer}

  if n <= 0
    error("The number of nodes must be positive.")
  end
        
  nodes = zeros(n)
        
  if isodd(n)               
    nodes[div(n-1,2)+1] = (domain[1]+domain[2])/2.0
  end
  
  if n == 1

    return nodes
  
  elseif n == 2
  
    nodes[1] = domain[2]
    nodes[n] = domain[1]
  
    return nodes
  
  else
        
    nodes[1] = domain[2]
    nodes[n] = domain[1]
  
    for i = 2:div(n,2)
      x = -cos((pi/2)*(2(i-1)+1.0)/(n))/cos((pi/(2*(n)))*(1.0+1/(4*log(n))))*(domain[1]-domain[2])/2.0
      nodes[i]     = (domain[1]+domain[2])/2.0 + x
      nodes[n-i+1] = (domain[1]+domain[2])/2.0 - x
    end
        
    return nodes
      
  end
      
end

function legendre_nodes(N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = Array{Float64,1}(undef, N)
  update = Array{Float64,1}(undef, N)

  p = zeros(N, N)

  if isodd(N)
    points[div(N - 1, 2)+1] = 0.0
  end

  if N == 1
    point[1] = (domain[1] + domain[2]) * 0.5
    return points
  else
    points[begin] = -1.0
    points[end] = 1.0

    @inbounds for i = 2:div((N - 1), 2)
      points[i] = -cos(π * (i - 1) / (N - 1))
      points[end-i+1] = cos(π * (i - 1) / (N - 1))
    end

    p[:, 1] .= 1.0

    len = Inf
    @views while len > eps(1.0)
      p[:, 2] .= points
      @inbounds for i = 2:(N-1)
        p[:, i+1] .= ((2i - 1) * points .* p[:, i] .- (i - 1) * p[:, i-1]) / i
      end
      update .= (points .* p[:, end] .- p[:, end-1]) ./ (N * p[:, end])
      points .-= update
      len = maximum(abs, update)
    end

    points .= (domain[1] + domain[2]) * 0.5 .+ points * ((domain[1] - domain[2]) * 0.5)

    return points
  end

end

function chebyshev_quadrature(node_type::Function,n::S,domain = [1.0,-1.0]) where {S <: Integer}

  if n <= 0
    error("The number of nodes must be positive.")
  end
        
  nodes   = zeros(n)
  weights = zeros(n)
  
  if node_type == chebyshev_nodes # weighting function is (1.0 .- nodes.^2).^(-1/2) -- (first kind)
  
    for i = 1:n
      nodes[i] = (domain[1]+domain[2])/2.0 - cos((i-0.5)*pi/n)*(domain[1]-domain[2])/2.0
      weights[i] = pi/n
    end
  
    return nodes, weights
  
  elseif node_type == chebyshev_extrema # weighting function is (1.0 .- nodes.^2).^(1/2) -- (second kind)
  
    for i = 1:n
      nodes[i] = (domain[1]+domain[2])/2.0 - cos((i-1)*pi/(n-1))*(domain[1]-domain[2])/2.0
      weights[i] = (pi/(n-1))*sin(((i-1)/(n-1))*pi)^2
    end
  
    return nodes, weights
      
  else
      
    error("This node type is not supported")
  
  end
  
end

function normalize_node(node::R,domain::Array{T,1}) where {R<:Number,T<:AbstractFloat}

  if domain[1] == domain[2]
    norm_node = zero(T)
    return norm_node
  else
    norm_node = 2.0*(node-domain[2])/(domain[1]-domain[2])-1.0
    return norm_node
  end
  
end
  
function normalize_node(node::Array{R,1},domain::Array{T,1}) where {R<:Number,T<:AbstractFloat}
  
  norm_nodes = similar(node)
  for i in eachindex(node)
    norm_nodes[i] = normalize_node(node[i],domain)
  end
  
  return norm_nodes
  
end
  
function normalize_node!(node::Array{R,1},domain::Array{T,1}) where {R<:Number,T<:AbstractFloat}
  
  for i in eachindex(node)
    node[i] = normalize_node(node[i],domain)
  end
  
end

function chebyshev_polynomial(order::S,x::T) where {T<:Number,S<:Integer}

  polynomial    = Array{T}(undef,1,order+1)
  polynomial[1] = one(T)
  
  for i = 2:order+1
    if i == 2
      polynomial[i] = x
    else
      polynomial[i] = 2*x*polynomial[i-1]-polynomial[i-2]
    end
  end
  
  return polynomial
  
end
  
function chebyshev_polynomial(order::S,x::AbstractArray{T,1}) where {T<:Number,S<:Integer}
  
  polynomial      = Array{T}(undef,length(x),order+1)
  polynomial[:,1] = ones(T,length(x))
  
  for i = 2:order+1
    for j in eachindex(x)
      if i == 2
        polynomial[j,i] = x[j]
      else
        polynomial[j,i] = 2*x[j]*polynomial[j,i-1]-polynomial[j,i-2]
      end
    end
  end
  
  return polynomial
  
end
  
function derivative_of_chebyshev_polynomial(order::S,x::T) where {T<:Number,S<:Integer}
  
  polynomial    = Array{T}(undef,1,order+1)
  poly_deriv    = Array{T}(undef,1,order+1)
  polynomial[1] = one(T)
  poly_deriv[1] = zero(T)
  
  for i = 2:order+1
    if i == 2
      polynomial[i] = x
      poly_deriv[i] = one(T)
    else
      polynomial[i] = 2*x*polynomial[i-1]-polynomial[i-2]
      poly_deriv[i] = 2*polynomial[i-1]+2*x*poly_deriv[i-1]-poly_deriv[i-2] 
    end
  end
  
  return poly_deriv
  
end

function derivative_of_chebyshev_polynomial(order::S,x::AbstractArray{T,1}) where {S<:Integer,T<:Number}

  polynomial    = Array{T}(undef,length(x),order+1)
  poly_deriv    = Array{T}(undef,length(x),order+1)
  polynomial[:,1] = ones(T,length(x))
  poly_deriv[:,1] = zeros(T,length(x))

  for i = 2:order+1
      for j in eachindex(x)
          if i == 2
              polynomial[j,i] = x[j]
              poly_deriv[j,i] = one(T)
          else
              polynomial[j,i] = 2*x[j]*polynomial[j,i-1]-polynomial[j,i-2]
              poly_deriv[j,i] = 2*polynomial[j,i-1]+2*x[j]*poly_deriv[j,i-1]-poly_deriv[j,i-2] 
          end
      end
  end

  return poly_deriv

end

function chebyshev_weights(f::Array{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i],normalize_node(nodes[i],domain[:,i]))
  end
  
  weights = Array{T,N}(undef,Tuple(order.+1))
    
  @inbounds for i in CartesianIndices(weights)
  
    numerator   = zero(T)
    denominator = zero(T)
  
    @inbounds for s in CartesianIndices(f)
  
      product = one(T)
      @inbounds for j = 1:N
        product *= poly[j][s[j],i[j]]
      end
  
      numerator   += f[s]*product
      denominator += product^2
  
    end
  
    weights[i] = numerator/denominator
  
  end
  
  return weights
  
end
  
function chebyshev_weights_extrema(f::Array{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  n = size(f)
  
  poly = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i],normalize_node(nodes[i],domain[:,i]))
  end
  
  weights = Array{T,N}(undef,Tuple(order.+1))
    
  @inbounds for i in CartesianIndices(weights)
  
    numerator   = zero(T)
    denominator = zero(T)
  
    @inbounds for s in CartesianIndices(f)
  
      num = f[s]
      den = one(T)
      @inbounds for j = 1:N
        if s[j] === 1 || s[j] === n[j]
          scale = 0.5  
        else
          scale = 1.0
        end
        temp = poly[j][s[j],i[j]]
        num *= temp*scale
        den *= (temp^2)*scale
      end
  
      numerator   += num
      denominator += den
  
    end
  
    weights[i] = numerator/denominator
  
  end
  
  return weights
  
end
  
function chebyshev_weights_extended(f::Array{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  poly       = Array{Array{T,2},1}(undef,N)
  complement = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    poly[i]       = chebyshev_polynomial(order[i],normalize_node(nodes[i],domain[:,i]))
    complement[i] = pinv(poly[i])'
  end
  
  weights = Array{T,N}(undef,Tuple(order.+1))
    
  @inbounds for i in CartesianIndices(weights)
  
    numerator   = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)
  
      num = f[s]
      den = one(T)
      @inbounds for j = 1:N
        temp = complement[j][s[j],i[j]]
        num *= temp
        den *= temp*poly[j][s[j],i[j]]
      end
  
      numerator   += num
      denominator += den
  
    end
  
    weights[i] = numerator/denominator
  
  end
  
  return weights
  
end
  
function chebyshev_weights(f::Array{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = Array{T,N}(undef,Tuple(order.+1))
    
  @inbounds for i in CartesianIndices(weights)
  
    numerator   = zero(T)
    denominator = zero(T)
  
    @inbounds for s in CartesianIndices(f)
  
      product = one(T)
      @inbounds for j = 1:N
        product *= poly[j][s[j],i[j]]
      end
  
      numerator   += f[s]*product
      denominator += product^2
  
    end
  
    weights[i] = numerator/denominator
  
  end
  
  return weights
  
end
  
function chebyshev_weights_extrema(f::Array{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  n = size(f)
  
  weights = Array{T,N}(undef,Tuple(order.+1))
    
  @inbounds for i in CartesianIndices(weights)
  
    numerator   = zero(T)
    denominator = zero(T)
  
    @inbounds for s in CartesianIndices(f)
  
      num = f[s]
      den = one(T)
      @inbounds for j = 1:N
        if s[j] == 1 || s[j] == n[j]
          scale = 0.5  
        else
          scale = 1.0
        end
        temp = poly[j][s[j],i[j]]
        num *= temp*scale
        den *= (temp^2)*scale
      end
  
      numerator   += num
      denominator += den
  
    end
  
    weights[i] = numerator/denominator
  
  end
  
  return weights
  
end
  
function chebyshev_weights_extended(f::Array{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  complement = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end
  
  weights = Array{T,N}(undef,Tuple(order.+1))
    
  @inbounds for i in CartesianIndices(weights)
  
    numerator   = zero(T)
    denominator = zero(T)
  
    @inbounds for s in CartesianIndices(f)
  
      num = f[s]
      den = one(T)
      @inbounds for j = 1:N
        temp = complement[j][s[j],i[j]]
        num *= temp
        den *= temp*poly[j][s[j],i[j]]
      end
  
      numerator   += num
      denominator += den
  
    end
  
    weights[i] = numerator/denominator
  
  end
  
  return weights
  
end
  
function chebyshev_weights(f::Array{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  poly = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order,normalize_node(nodes[i],domain[:,i]))
  end
  
  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
  
  weights = Array{T,N}(undef,ord.+1)
  
  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
  
      numerator   = zero(T)
      denominator = zero(T)
  
      @inbounds for s in CartesianIndices(f)
  
        product = one(T)
        @inbounds for j = 1:N
          product *= poly[j][s[j],i[j]]
        end
    
        numerator   += f[s]*product
        denominator += product^2
    
      end
  
      weights[i] = numerator/denominator
      
    else
      weights[i] = zero(T)
    end
  
  end
  
  return weights
  
end
  
function chebyshev_weights_extrema(f::Array{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  n = size(f)
    
  poly = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order,normalize_node(nodes[i],domain[:,i]))
  end
  
  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
  
  weights = Array{T,N}(undef,ord.+1)
    
  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
 
      numerator   = zero(T)
      denominator = zero(T)
    
      @inbounds for s in CartesianIndices(f)
    
        num = f[s]
        den = one(T)
        @inbounds for j = 1:N
          if s[j] == 1 || s[j] == n[j]
            scale = 0.5  
          else
            scale = 1.0
          end
          temp = poly[j][s[j],i[j]]
          num *= temp*scale
          den *= (temp^2)*scale
        end
    
        numerator   += num
        denominator += den
    
      end
    
      weights[i] = numerator/denominator
        
    else
      weights[i] = zero(T)
    end
    
  end
    
  return weights
    
end
  
function chebyshev_weights_extended(f::Array{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  poly       = Array{Array{T,2},1}(undef,N)
  complement = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    poly[i]       = chebyshev_polynomial(order,normalize_node(nodes[i],domain[:,i]))
    complement[i] = pinv(poly[i])'
  end
  
  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
  
  weights = Array{T,N}(undef,ord.+1)
    
  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
  
      numerator   = zero(T)
      denominator = zero(T)
    
      @inbounds for s in CartesianIndices(f)
    
        num = f[s]
        den = one(T)
        @inbounds for j = 1:N
          temp = complement[j][s[j],i[j]]
          num *= temp
          den *= temp*poly[j][s[j],i[j]]
        end
    
        numerator   += num
        denominator += den
    
      end
    
      weights[i] = numerator/denominator
        
    else
      weights[i] = zero(T)
    end
    
  end
    
  return weights
    
end
  
function chebyshev_weights(f::Array{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
    
  weights = Array{T,N}(undef,ord.+1)
    
  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
  
      numerator   = zero(T)
      denominator = zero(T)
    
      @inbounds for s in CartesianIndices(f)
    
        product = one(T)
        @inbounds for j = 1:N
          product *= poly[j][s[j],i[j]]
        end
    
        numerator   += f[s]*product
        denominator += product^2
      
      end
    
      weights[i] = numerator/denominator
        
    else
      weights[i] = zero(T)
    end
    
  end
    
  return weights
    
end
  
function chebyshev_weights_extrema(f::Array{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  n = size(f)
  
  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
    
  weights = Array{T,N}(undef,ord.+1)
    
  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
  
      numerator   = zero(T)
      denominator = zero(T)
    
      @inbounds for s in CartesianIndices(f)
    
        num = f[s]
        den = one(T)
        @inbounds for j = 1:N
          if s[j] == 1 || s[j] == n[j]
            scale = 0.5  
          else
            scale = 1.0
          end
          temp = poly[j][s[j],i[j]]
          num *= temp*scale
          den *= (temp^2)*scale
        end
    
        numerator   += num
        denominator += den
    
      end
    
      weights[i] = numerator/denominator
        
    else
      weights[i] = zero(T)
    end
    
  end
    
  return weights
    
end
  
function chebyshev_weights_extended(f::Array{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  complement = Array{Array{T,2},1}(undef,N)
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end
  
  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
  
  weights = Array{T,N}(undef,ord.+1)
    
  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N

      numerator   = zero(T)
      denominator = zero(T)
    
      @inbounds for s in CartesianIndices(f)
    
        num = f[s]
        den = one(T)
        @inbounds for j = 1:N
          temp = complement[j][s[j],i[j]]
          num *= temp
          den *= temp*poly[j][s[j],i[j]]
        end
    
        numerator   += num
        denominator += den
    
      end
    
      weights[i] = numerator/denominator
        
    else
      weights[i] = zero(T)
    end
    
  end
    
  return weights
    
end
  
function chebyshev_weights_threaded(f::Array{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  poly = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i],normalize_node(nodes[i],domain[:,i]))
  end
  
  weights = zeros(Tuple(order.+1))
    
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
  
    numerator   = zero(T)
    denominator = zero(T)
  
    @inbounds for s in CartesianIndices(f)
  
      product = one(T)
      @inbounds for j = 1:N
        product *= poly[j][s[j],i[j]]
      end
  
      numerator   += f[s]*product
      denominator += product^2
  
    end
  
    weights[i] = numerator/denominator
  
  end
  
  return weights
  
end
  
function chebyshev_weights_extrema_threaded(f::Array{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  n = size(f)
  
  poly = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i],normalize_node(nodes[i],domain[:,i]))
  end
  
  weights = zeros(Tuple(order.+1))
    
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
  
    numerator   = zero(T)
    denominator = zero(T)
  
    @inbounds for s in CartesianIndices(f)
  
      num = f[s]
      den = one(T)
      @inbounds for j = 1:N
        if s[j] == 1 || s[j] == n[j]
          scale = 0.5  
        else
          scale = 1.0
        end
        temp = poly[j][s[j],i[j]]
        num *= temp*scale
        den *= (temp^2)*scale
      end
  
      numerator   += num
      denominator += den
  
    end
  
    weights[i] = numerator/denominator
  
  end
  
  return weights
  
end
  
function chebyshev_weights_extended_threaded(f::Array{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  poly       = Array{Array{T,2},1}(undef,N)
  complement = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    poly[i]       = chebyshev_polynomial(order[i],normalize_node(nodes[i],domain[:,i]))
    complement[i] = pinv(poly[i])'
  end
  
  weights = zeros(Tuple(order.+1))
    
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
  
    numerator   = zero(T)
    denominator = zero(T)
  
    @inbounds for s in CartesianIndices(f)
  
      num = f[s]
      den = one(T)
      @inbounds for j = 1:N
        temp = complement[j][s[j],i[j]]
        num *= temp
        den *= temp*poly[j][s[j],i[j]]
      end
  
      numerator   += num
      denominator += den
  
    end
  
    weights[i] = numerator/denominator
  
  end
  
  return weights
  
end
  
function chebyshev_weights_threaded(f::Array{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = Array{T,N}(undef,Tuple(order.+1))
    
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
  
    numerator   = zero(T)
    denominator = zero(T)
  
    @inbounds for s in CartesianIndices(f)
  
      product = one(T)
      @inbounds for j = 1:N
        product *= poly[j][s[j],i[j]]
      end
  
      numerator   += f[s]*product
      denominator += product^2
  
    end
  
    weights[i] = numerator/denominator
  
  end
  
  return weights
  
end
  
function chebyshev_weights_extrema_threaded(f::Array{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  n = size(f)
  
  weights = Array{T,N}(undef,Tuple(order.+1))
    
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
  
    numerator   = zero(T)
    denominator = zero(T)
  
    @inbounds for s in CartesianIndices(f)
  
      num = f[s]
      den = one(T)
      @inbounds for j = 1:N
        if s[j] == 1 || s[j] == n[j]
          scale = 0.5  
        else
          scale = 1.0
        end
        temp = poly[j][s[j],i[j]]
        num *= temp*scale
        den *= (temp^2)*scale
      end
  
      numerator   += num
      denominator += den
  
    end
  
    weights[i] = numerator/denominator
  
  end
  
  return weights
  
end
  
function chebyshev_weights_extended_threaded(f::Array{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  complement = Array{Array{T,2},1}(undef,N)
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end
  
  weights = Array{T,N}(undef,Tuple(order.+1))
    
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
  
    numerator   = zero(T)
    denominator = zero(T)
  
    @inbounds for s in CartesianIndices(f)
  
      num = f[s]
      den = one(T)
      @inbounds for j = 1:N
        temp = complement[j][s[j],i[j]]
        num *= temp
        den *= temp*poly[j][s[j],i[j]]
      end
  
      numerator   += num
      denominator += den
  
    end
  
    weights[i] = numerator/denominator
  
  end
  
  return weights
  
end
  
function chebyshev_weights_threaded(f::Array{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  poly = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order,normalize_node(nodes[i],domain[:,i]))
  end
  
  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
    
  weights = Array{T,N}(undef,ord.+1)
    
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
  
      numerator   = zero(T)
      denominator = zero(T)
    
      @inbounds for s in CartesianIndices(f)
    
        product = one(T)
        @inbounds for j = 1:N
          product *= poly[j][s[j],i[j]]
        end
    
        numerator   += f[s]*product
        denominator += product^2
      
      end
    
      weights[i] = numerator/denominator
        
    else
      weights[i] = zero(T)
    end
    
  end
    
  return weights
    
end
  
function chebyshev_weights_extrema_threaded(f::Array{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  n = size(f)
  
  poly = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order,normalize_node(nodes[i],domain[:,i]))
  end
  
  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
    
  weights = Array{T,N}(undef,ord.+1)
    
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
  
      numerator   = zero(T)
      denominator = zero(T)
    
      @inbounds for s in CartesianIndices(f)
    
        num = f[s]
        den = one(T)
        @inbounds for j = 1:N
          if s[j] == 1 || s[j] == n[j]
            scale = 0.5  
          else
            scale = 1.0
          end
          temp = poly[j][s[j],i[j]]
          num *= temp*scale
          den *= (temp^2)*scale
        end
    
        numerator   += num
        denominator += den
    
      end
    
      weights[i] = numerator/denominator
        
    else
      weights[i] = zero(T)
    end
    
  end
    
  return weights
    
end
  
function chebyshev_weights_extended_threaded(f::Array{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  poly       = Array{Array{T,2},1}(undef,N)
  complement = Array{Array{T,2},1}(undef,N)
  
  @inbounds for i = 1:N
    poly[i]       = chebyshev_polynomial(order,normalize_node(nodes[i],domain[:,i]))
    complement[i] = pinv(poly[i])'
  end
  
  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
    
  weights = Array{T,N}(undef,ord.+1)
    
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
  
      numerator   = zero(T)
      denominator = zero(T)
    
      @inbounds for s in CartesianIndices(f)
    
        num = f[s]
        den = one(T)
        @inbounds for j = 1:N
          temp = complement[j][s[j],i[j]]
          num *= temp
          den *= temp*poly[j][s[j],i[j]]
        end
    
        numerator   += num
        denominator += den
    
      end
    
      weights[i] = numerator/denominator
        
    else
      weights[i] = zero(T)
    end
    
  end
   
  return weights
    
end
  
function chebyshev_weights_threaded(f::Array{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
      
  weights = Array{T,N}(undef,ord.+1)
    
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
  
      numerator   = zero(T)
      denominator = zero(T)
    
      @inbounds for s in CartesianIndices(f)
    
        product = one(T)
        @inbounds for j = 1:N
          product *= poly[j][s[j],i[j]]
        end
    
        numerator   += f[s]*product
        denominator += product^2
      
      end
    
      weights[i] = numerator/denominator
        
    else
      weights[i] = zero(T)
    end
    
  end
    
  return weights
    
end
  
function chebyshev_weights_extrema_threaded(f::Array{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  n = size(f)
  
  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
      
  weights = Array{T,N}(undef,ord.+1)
    
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
  
      numerator   = zero(T)
      denominator = zero(T)
    
      @inbounds for s in CartesianIndices(f)
    
        num = f[s]
        den = one(T)
        @inbounds for j = 1:N
          if s[j] == 1 || s[j] == n[j]
            scale = 0.5  
          else
            scale = 1.0
          end
          temp = poly[j][s[j],i[j]]
          num *= temp*scale
          den *= (temp^2)*scale
        end
    
        numerator   += num
        denominator += den
    
      end
    
      weights[i] = numerator/denominator
        
    else
      weights[i] = zero(T)
    end
    
  end
    
  return weights
    
end
  
function chebyshev_weights_extended_threaded(f::Array{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  complement = Array{Array{T,2},1}(undef,N)
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end
  
  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
      
  weights = Array{T,N}(undef,ord.+1)
    
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
  
      numerator   = zero(T)
      denominator = zero(T)
    
      @inbounds for s in CartesianIndices(f)
    
        num = f[s]
        den = one(T)
        @inbounds for j = 1:N
          temp = complement[j][s[j],i[j]]
          num *= temp
          den *= temp*poly[j][s[j],i[j]]
        end
    
        numerator   += num
        denominator += den
    
      end
    
      weights[i] = numerator/denominator
        
    else
      weights[i] = zero(T)
    end
    
  end
    
  return weights
    
end

function chebyshev_weights(cheb::ChebInterpRoots)
  
  weights = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  
  return weights
  
end
  
function chebyshev_weights(cheb::ChebInterpExtrema)
  
  weights = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  
  return weights
  
end
  
function chebyshev_weights(cheb::ChebInterpExtended)
  
  weights = chebyshev_weights_extended_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  
  return weights
  
end

function chebyshev_weights(cheb::ChebInterpVertesi)
  
  weights = chebyshev_weights_vertesi_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  
  return weights
  
end

# Functions for the one-variable case where the nodes are a vector
  
function chebyshev_weights(f::Array{T,1},nodes::Array{T,1},order::Union{S,Array{S,1}},domain=[one(T);-one(T)]) where {T<:AbstractFloat,S<:Integer}
  
  weights = chebyshev_weights(f,(nodes,),order,domain)
    
  return weights
    
end
  
function chebyshev_weights_extrema(f::Array{T,1},nodes::Array{T,1},order::Union{S,Array{S,1}},domain=[one(T);-one(T)]) where {T<:AbstractFloat,S<:Integer}
  
  weights = chebyshev_weights_extrema(f,(nodes,),order,domain)
    
  return weights
    
end
  
function chebyshev_weights_extended(f::Array{T,1},nodes::Array{T,1},order::Union{S,Array{S,1}},domain=[one(T);-one(T)]) where {T<:AbstractFloat,S<:Integer}
  
  weights = chebyshev_weights_extended(f,(nodes,),order,domain)
    
  return weights
    
end
  
function chebyshev_weights(f::Array{T,1},poly::Array{T,2},order::Union{S,Array{S,1}}) where {T<:AbstractFloat,S<:Integer}
    
  weights = chebyshev_weights(f,(poly,),order)
    
  return weights
    
end
    
function chebyshev_weights_extrema(f::Array{T,1},poly::Array{T,2},order::Union{S,Array{S,1}}) where {T<:AbstractFloat,S<:Integer}
    
  weights = chebyshev_weights_extrema(f,(poly,),order)
    
  return weights
    
end
  
function chebyshev_weights_extended(f::Array{T,1},poly::Array{T,2},order::Union{S,Array{S,1}}) where {T<:AbstractFloat,S<:Integer}
    
  weights = chebyshev_weights_extended(f,(poly,),order)
    
  return weights
    
end
  
# Functions that allow the nodes to be in an array of arrays
  
function chebyshev_weights(f::Array{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights(f,tuple(nodes...),order,domain)
  
  return weights
  
end
  
function chebyshev_weights_extrema(f::Array{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
 
  weights = chebyshev_weights_extrema(f,tuple(nodes...),order,domain)
  
  return weights
  
end
  
function chebyshev_weights_extended(f::Array{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extended(f,tuple(nodes...),order,domain)
  
  return weights
  
end
  
function chebyshev_weights(f::Array{T,N},poly::Array{Array{T,2},1},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights(f,tuple(poly...),order)
  
  return weights
  
end
  
function chebyshev_weights_extrema(f::Array{T,N},poly::Array{Array{T,2},1},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extrema(f,tuple(poly...),order)
  
  return weights
  
end
  
function chebyshev_weights_extended(f::Array{T,N},poly::Array{Array{T,2},1},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extended(f,tuple(poly...),order)
  
  return weights
  
end
  
function chebyshev_weights(f::Array{T,N},nodes::Array{Array{T,1},1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights(f,tuple(nodes...),order,domain)
  
  return weights
  
end
  
function chebyshev_weights_extrema(f::Array{T,N},nodes::Array{Array{T,1},1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extrema(f,tuple(nodes...),order,domain)
  
  return weights
  
end
  
function chebyshev_weights_extended(f::Array{T,N},nodes::Array{Array{T,1},1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extended(f,tuple(nodes...),order,domain)
  
  return weights
  
end
  
function chebyshev_weights(f::Array{T,N},poly::Array{Array{T,2},1},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights(f,tuple(poly...),order)
  
  return weights
  
end
  
function chebyshev_weights_extrema(f::Array{T,N},poly::Array{Array{T,2},1},order::S) where {T<:AbstractFloat,N,S<:Integer}
 
  weights = chebyshev_weights_extrema(f,tuple(poly...),order)
  
  return weights
  
end
  
function chebyshev_weights_extended(f::Array{T,N},poly::Array{Array{T,2},1},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extended(f,tuple(poly...),order)
  
  return weights
  
end
  
# Threaded functions
  
function chebyshev_weights_threaded(f::Array{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_threaded(f,tuple(nodes...),order,domain)
  
  return weights
  
end
  
function chebyshev_weights_extrema_threaded(f::Array{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extrema_threaded(f,tuple(nodes...),order,domain)
  
  return weights

end
  
function chebyshev_weights_extended_threaded(f::Array{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extended_threaded(f,tuple(nodes...),order,domain)
  
  return weights
  
end
  
function chebyshev_weights_threaded(f::Array{T,N},poly::Array{Array{T,2},1},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_threaded(f,tuple(poly...),order)
  
  return weights
  
end
  
function chebyshev_weights_extrema_threaded(f::Array{T,N},poly::Array{Array{T,2},1},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extrema_threaded(f,tuple(poly...),order)
  
  return weights
 
end
  
function chebyshev_weights_extended_threaded(f::Array{T,N},poly::Array{Array{T,2},1},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extended_threaded(f,tuple(poly...),order)
  
  return weights
  
end
  
function chebyshev_weights_threaded(f::Array{T,N},nodes::Array{Array{T,1},1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_threaded(f,tuple(nodes...),order,domain)
  
  return weights
  
end
  
function chebyshev_weights_extrema_threaded(f::Array{T,N},nodes::Array{Array{T,1},1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extrema_threaded(f,tuple(nodes...),order,domain)
  
  return weights
  
end
  
function chebyshev_weights_extended_threaded(f::Array{T,N},nodes::Array{Array{T,1},1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extended_threaded(f,tuple(nodes...),order,domain)
  
  return weights
  
end
  
function chebyshev_weights_threaded(f::Array{T,N},poly::Array{Array{T,2},1},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_threaded(f,tuple(poly...),order)
  
  return weights
  
end
  
function chebyshev_weights_extrema_threaded(f::Array{T,N},poly::Array{Array{T,2},1},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extrema_threaded(f,tuple(poly...),order)
  
  return weights
  
end
  
function chebyshev_weights_extended_threaded(f::Array{T,N},poly::Array{Array{T,2},1},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extended_threaded(f,tuple(poly...),order)
  
  return weights
  
end

const chebyshev_weights_vertesi          = chebyshev_weights_extended
const chebyshev_weights_vertesi_threaded = chebyshev_weights_extended_threaded

# Regular functions for evaluating Chebyshev polynominals

function chebyshev_evaluate(weights::Array{T,N},x::AbstractArray{R,1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  poly = Array{Array{R,2},1}(undef,N)
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i],normalize_node(x[i],domain[:,i]))
  end
  
  yhat = zero(T)
  @inbounds for i in CartesianIndices(weights)
    poly_product = poly[1][i[1]]
    @inbounds for j = 2:N
      poly_product *= poly[j][i[j]]
    end
    yhat += weights[i]*poly_product
  end
  
  return yhat
  
end
  
function chebyshev_evaluate(weights::Array{T,N},x::AbstractArray{R,1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}
  
  poly = Array{Array{R,2},1}(undef,N)
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order,normalize_node(x[i],domain[:,i]))
  end
  
  yhat = zero(T)
  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
      poly_product = poly[1][i[1]]
      @inbounds for j = 2:N
        poly_product *= poly[j][i[j]]
      end
      yhat += weights[i]*poly_product
    end
  end
  
  return yhat
  
end
  
function chebyshev_evaluate(cheb_poly::ChebPoly,x::AbstractArray{R,1}) where {R<:Number}
  
  yhat = chebyshev_evaluate(cheb_poly.weights,x,cheb_poly.order,cheb_poly.domain) 
      
  return yhat
    
end
  
function cheb_interp(cheb::ChebInterpRoots,x::AbstractArray{R,1}) where {R<:Number}
  
  weights = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  yhat    = chebyshev_evaluate(weights,x,cheb.order,cheb.domain) 
     
  return yhat
    
end
  
function cheb_interp(cheb::ChebInterpExtrema,x::AbstractArray{R,1}) where {R<:Number}
  
  weights = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  yhat    = chebyshev_evaluate(weights,x,cheb.order,cheb.domain) 
      
  return yhat
    
end
  
function cheb_interp(cheb::ChebInterpExtended,x::AbstractArray{R,1}) where {R<:Number}
  
  weights = chebyshev_weights_extended_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  yhat    = chebyshev_evaluate(weights,x,cheb.order,cheb.domain) 
      
  return yhat
    
end

function cheb_interp(cheb::ChebInterpVertesi,x::AbstractArray{R,1}) where {R<:Number}
  
  weights = chebyshev_weights_vertesi_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  yhat    = chebyshev_evaluate(weights,x,cheb.order,cheb.domain) 
        
  return yhat
      
end
  
function chebyshev_evaluate(weights::Array{T,N},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}
  
  function chebeval(x::AbstractArray{R,1}) where {R <: Number}
  
    return chebyshev_evaluate(weights,x,order,domain)
  
  end
  
  return chebeval
  
end
  
function chebyshev_evaluate(weights::Array{T,N},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}
  
  function chebeval(x::AbstractArray{R,1}) where {R <: Number}
  
    return chebyshev_evaluate(weights,x,order,domain)
  
  end
  
  return chebeval
  
end
  
function chebyshev_evaluate(cheb_poly::ChebPoly)
  
  function chebeval(x::AbstractArray{R,1}) where {R <: Number}
  
    return chebyshev_evaluate(cheb_poly,x)
  
  end
  
  return chebeval

end
  
function cheb_interp(cheb::ChebInterpRoots)
 
  weights = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
    
  function chebeval(x::AbstractArray{R,1}) where {R <: Number}
  
    return chebyshev_evaluate(weights,x,cheb.order,cheb.domain)
  
  end
  
  return chebeval
  
end
  
function cheb_interp(cheb::ChebInterpExtrema)
  
  weights = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
    
  function chebeval(x::AbstractArray{R,1}) where {R <: Number}
  
    return chebyshev_evaluate(weights,x,cheb.order,cheb.domain)
  
  end
  
  return chebeval
  
end
  
function cheb_interp(cheb::ChebInterpExtended)
  
  weights = chebyshev_weights_extended_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
    
  function chebeval(x::AbstractArray{R,1}) where {R <: Number}
  
    return chebyshev_evaluate(weights,x,cheb.order,cheb.domain)
  
  end
  
  return chebeval
  
end

function cheb_interp(cheb::ChebInterpVertesi)
  
  weights = chebyshev_weights_vertesi_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
      
  function chebeval(x::AbstractArray{R,1}) where {R <: Number}
    
    return chebyshev_evaluate(weights,x,cheb.order,cheb.domain)
    
  end
    
  return chebeval
    
end
  
# Regular functions for evaluating derivatives and gradients of Chebyshev polynomials

# Regular functions for derivatives

function chebyshev_derivative(weights::Array{T,N},x::AbstractArray{R,1},pos::S,order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  poly = Array{Array{R,2},1}(undef,N)
  @inbounds for i = 1:N
    if i === pos
      poly[i] = derivative_of_chebyshev_polynomial(order[i],normalize_node(x[i],domain[:,i]))
    else
      poly[i] = chebyshev_polynomial(order[i],normalize_node(x[i],domain[:,i]))
    end
  end
  
  derivative = zero(T)
  @inbounds for i in CartesianIndices(weights)
    poly_product = poly[1][i[1]]
    @inbounds for j = 2:N
      poly_product *= poly[j][i[j]]
    end
    derivative += weights[i]*poly_product
  end
  
  return derivative*(2.0/(domain[1,pos]-domain[2,pos]))
  
end
  
function chebyshev_derivative(weights::Array{T,N},x::AbstractArray{R,1},pos::S,order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}
  
  poly = Array{Array{R,2},1}(undef,N)
  @inbounds for i = 1:N
    if i === pos
      poly[i] = derivative_of_chebyshev_polynomial(order,normalize_node(x[i],domain[:,i]))
    else
      poly[i] = chebyshev_polynomial(order,normalize_node(x[i],domain[:,i]))
    end
  end
  
  derivative = zero(T)
  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
      poly_product = poly[1][i[1]]
      @inbounds for j = 2:N
        poly_product *= poly[j][i[j]]
      end
      derivative += weights[i]*poly_product
    end
  end
  
  return derivative*(2.0/(domain[1,pos]-domain[2,pos]))
  
end
  
function chebyshev_derivative(cheb_poly::ChebPoly,x::AbstractArray{R,1},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  derivative = chebyshev_derivative(cheb_poly.weights,x,pos,cheb_poly.order,cheb_poly.domain) 
      
  return derivative
  
end
  
function chebyshev_derivative(cheb::ChebInterpRoots,x::AbstractArray{R,1},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  weights    = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  derivative = chebyshev_derivative(weights,x,pos,cheb.order,cheb.domain) 
      
  return derivative
  
end
  
function chebyshev_derivative(cheb::ChebInterpExtrema,x::AbstractArray{R,1},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  weights    = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  derivative = chebyshev_derivative(weights,x,pos,cheb.order,cheb.domain) 
    
  return derivative

end
  
function chebyshev_derivative(cheb::ChebInterpExtended,x::AbstractArray{R,1},pos::S) where {T<:AbstractFloat,R<:Number,S<:Integer}
  
  weights    = chebyshev_weights_extended_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  derivative = chebyshev_derivative(weights,x,pos,cheb.order,cheb.domain) 
      
  return derivative
  
end
  
# Regular functions for gradients
  
function chebyshev_gradient(weights::Array{T,N},x::AbstractArray{R,1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}
  
  gradient = Array{R,2}(undef,1,N)
  
  @inbounds for i = 1:N
    gradient[i] = chebyshev_derivative(weights,x,i,order,domain)
  end
  
  return gradient
  
end
  
function chebyshev_gradient(weights::Array{T,N},x::AbstractArray{R,1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}
  
  gradient = Array{R,2}(undef,1,N)
  
  @inbounds for i = 1:N
    gradient[i] = chebyshev_derivative(weights,x,i,order,domain)
  end
  
  return gradient
 
end
  
function chebyshev_gradient(cheb_poly::ChebPoly,x::AbstractArray{R,1}) where {R<:Number}
  
  gradient = chebyshev_gradient(cheb_poly.weights,x,cheb_poly.order,cheb_poly.domain) 
      
  return gradient
  
end
  
function chebyshev_gradient(cheb::ChebInterpRoots,x::AbstractArray{R,1}) where {R<:Number}
  
  weights  = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  gradient = chebyshev_gradient(weights,x,cheb.order,cheb.domain) 
      
  return gradient
  
end
  
function chebyshev_gradient(cheb::ChebInterpExtrema,x::AbstractArray{R,1}) where {R<:Number}
  
  weights  = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  gradient = chebyshev_gradient(weights,x,cheb.order,cheb.domain) 
      
  return gradient
  
end
  
function chebyshev_gradient(cheb::ChebInterpExtended,x::AbstractArray{R,1}) where {R<:Number}
  
  weights  = chebyshev_weights_extended_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  gradient = chebyshev_gradient(weights,x,cheb.order,cheb.domain) 
      
  return gradient
  
end

function chebyshev_gradient(cheb::ChebInterpVertesi,x::AbstractArray{R,1}) where {R<:Number}
  
  weights  = chebyshev_weights_vertesi_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  gradient = chebyshev_gradient(weights,x,cheb.order,cheb.domain) 
        
  return gradient
    
end
  
function chebyshev_derivative(weights::Array{T,N},pos::S,order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}
  
  function chebderiv(x::AbstractArray{R,1}) where {R<:Number}
  
    return chebyshev_derivative(weights,x,pos,order,domain)
  
  end
  
  return chebderiv
  
end
  
function chebyshev_derivative(weights::Array{T,N},pos::S,order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}
  
  function chebderiv(x::AbstractArray{R,1}) where {R<:Number}
  
    return chebyshev_derivative(weights,x,pos,order,domain)
  
  end
  
  return chebderiv
  
end
  
function chebyshev_derivative(cheb_poly::ChebPoly,pos::S) where {S <: Integer}
  
  function chebderiv(x::AbstractArray{R,1}) where {R<:Number}
  
    return chebyshev_derivative(cheb_poly,x,pos)
  
  end
  
  return chebderiv
  
end
  
function chebyshev_derivative(cheb::ChebInterpRoots,pos::S) where {S <: Integer}
  
  weights = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  
  function chebderiv(x::AbstractArray{R,1}) where {R<:Number}
  
    return chebyshev_derivative(weights,x,pos,cheb.order,cheb.domain)
  
  end
  
  return chebderiv
  
end
  
function chebyshev_derivative(cheb::ChebInterpExtrema,pos::S) where {S <: Integer}
  
  weights = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  
  function chebderiv(x::AbstractArray{R,1}) where {R<:Number}
  
    return chebyshev_derivative(weights,x,pos,cheb.order,cheb.domain)
  
  end
  
  return chebderiv

end
  
function chebyshev_derivative(cheb::ChebInterpExtended,pos::S) where {S <: Integer}
  
  weights = chebyshev_weights_extended_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  
  function chebderiv(x::AbstractArray{R,1}) where {R<:Number}
  
    return chebyshev_derivative(weights,x,pos,cheb.order,cheb.domain)
  
  end
  
  return chebderiv
  
end

function chebyshev_derivative(cheb::ChebInterpVertesi,pos::S) where {S <: Integer}
  
  weights = chebyshev_weights_vertesi_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
    
  function chebderiv(x::AbstractArray{R,1}) where {R<:Number}
    
    return chebyshev_derivative(weights,x,pos,cheb.order,cheb.domain)
    
  end
    
  return chebderiv
    
end
  
function chebyshev_gradient(weights::Array{T,N},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}
  
  function chebgrad(x::AbstractArray{R,1}) where {R<:Number}
  
    return chebyshev_gradient(weights,x,pos,order,domain)
  
  end
  
  return chebgrad
  
end
  
function chebyshev_gradient(weights::Array{T,N},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}
  
  function chebgrad(x::AbstractArray{R,1}) where {R<:Number}
  
    return chebyshev_gradient(weights,x,pos,order,domain)
  
  end
  
  return chebgrad
  
end
  
function chebyshev_gradient(cheb_poly::ChebPoly) where {S <: Integer}
  
  function chebgrad(x::AbstractArray{R,1}) where {R<:Number}
  
    return chebyshev_gradient(cheb_poly,x)
  
  end
  
  return chebgrad
  
end
  
function chebyshev_gradient(cheb::ChebInterpRoots) where {S <: Integer}
  
  weights = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  
  function chebgrad(x::AbstractArray{R,1}) where {R<:Number}
  
    return chebyshev_gradient(weights,x,cheb.order,cheb.domain)
  
  end
  
  return chebgrad
 
end
  
function chebyshev_gradient(cheb::ChebInterpExtrema) where {S <: Integer}
  
  weights = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)

  function chebgrad(x::AbstractArray{R,1}) where {R<:Number}
  
    return chebyshev_gradient(weights,x,cheb.order,cheb.domain)
  
  end
  
  return chebgrad
  
end
  
function chebyshev_gradient(cheb::ChebInterpExtended) where {S <: Integer}
  
  weights = chebyshev_weights_extended_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  
  function chebgrad(x::AbstractArray{R,1}) where {R<:Number}
  
    return chebyshev_gradient(weights,x,cheb.order,cheb.domain)
  
  end
  
  return chebgrad
  
end

function chebyshev_gradient(cheb::ChebInterpVertesi) where {S <: Integer}
  
  weights = chebyshev_weights_vertesi_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
    
  function chebgrad(x::AbstractArray{R,1}) where {R<:Number}
    
    return chebyshev_gradient(weights,x,cheb.order,cheb.domain)
    
  end
    
  return chebgrad
    
end