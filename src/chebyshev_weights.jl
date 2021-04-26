function chebyshev_weights(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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
      dom = one(T)
      @inbounds for j = 1:N
        num *= poly[j][s[j],i[j]]
        dom *= poly[j][s[j],i[j]]
      end

      numerator   += num
      denominator += dom^2

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights_extrema(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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
      dom = one(T)
      @inbounds for j = 1:N
        pow = (sum(s[j]==1)+sum(s[j]==length(nodes[j])))
        num *= poly[j][s[j],i[j]]/(2^pow)
        dom *= poly[j][s[j],i[j]]*poly[j][s[j],i[j]]/(2^pow)
      end

      numerator   += num
      denominator += dom

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef,Tuple(order.+1))
  
  @inbounds for i in CartesianIndices(weights)

    numerator   = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      num = f[s]
      dom = one(T)
      @inbounds for j = 1:N
        num *= poly[j][s[j],i[j]]
        dom *= poly[j][s[j],i[j]]
      end

      numerator   += num
      denominator += dom^2

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights_extrema(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef,Tuple(order.+1))
  
  @inbounds for i in CartesianIndices(weights)

    numerator   = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      num = f[s]
      dom = one(T)
      @inbounds for j = 1:N
        pow = (sum(s[j]==1)+sum(s[j]==length(nodes[j])))
        num *= poly[j][s[j],i[j]]/(2^pow)
        dom *= poly[j][s[j],i[j]]*poly[j][s[j],i[j]]/(2^pow)
      end

      numerator   += num
      denominator += dom

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order,normalize_node(nodes[i],domain[:,i]))
  end

  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end

  weights = zeros(ord.+1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N

      numerator   = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        num = f[s]
        dom = one(T)
        @inbounds for j = 1:N
          num *= poly[j][s[j],i[j]]
          dom *= poly[j][s[j],i[j]]
        end

        numerator   += num
        denominator += dom^2

      end

      weights[i] = numerator/denominator
    
    end

  end

  return weights

end

function chebyshev_weights_extrema(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order,normalize_node(nodes[i],domain[:,i]))
  end

  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end

  weights = zeros(ord.+1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N

      numerator   = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        num = f[s]
        dom = one(T)
        @inbounds for j = 1:N
          pow = (sum(s[j]==1)+sum(s[j]==length(nodes[j])))
          num *= poly[j][s[j],i[j]]/(2^pow)
          dom *= poly[j][s[j],i[j]]*poly[j][s[j],i[j]]/(2^pow)
        end

        numerator   += num
        denominator += dom

      end

      weights[i] = numerator/denominator
    
    end

  end

  return weights

end

function chebyshev_weights(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
  
  weights = zeros(Tuple(ord.+1))
  
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N

      numerator   = zero(T)
      denominator = zero(T)
  
      @inbounds for s in CartesianIndices(f)
  
        num = f[s]
        dom = one(T)
        @inbounds for j = 1:N
          num *= poly[j][s[j],i[j]]
          dom *= poly[j][s[j],i[j]]
        end
  
        numerator   += num
        denominator += dom^2
  
      end
  
      weights[i] = numerator/denominator
      
    end
  
  end
  
  return weights
  
end

function chebyshev_weights_extrema(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
  
  weights = zeros(Tuple(ord.+1))
  
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N

      numerator   = zero(T)
      denominator = zero(T)
  
      @inbounds for s in CartesianIndices(f)
  
        num = f[s]
        dom = one(T)
        @inbounds for j = 1:N
          pow = (sum(s[j]==1)+sum(s[j]==length(nodes[j])))
          num *= poly[j][s[j],i[j]]/(2^pow)
          dom *= poly[j][s[j],i[j]]*poly[j][s[j],i[j]]/(2^pow)
        end
  
        numerator   += num
        denominator += dom
  
      end
  
      weights[i] = numerator/denominator
      
    end
  
  end
  
  return weights
  
end

function chebyshev_weights_threaded(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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
      dom = one(T)
      @inbounds for j = 1:N
        num *= poly[j][s[j],i[j]]
        dom *= poly[j][s[j],i[j]]
      end

      numerator   += num
      denominator += dom^2

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights_extrema_threaded(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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
      dom = one(T)
      @inbounds for j = 1:N
        pow = (sum(s[j]==1)+sum(s[j]==length(nodes[j])))
        num *= poly[j][s[j],i[j]]/(2^pow)
        dom *= poly[j][s[j],i[j]]*poly[j][s[j],i[j]]/(2^pow)
      end

      numerator   += num
      denominator += dom

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights_threaded(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef,Tuple(order.+1))
  
  @inbounds @sync @qthreads for i in CartesianIndices(weights)

    numerator   = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      num = f[s]
      dom = one(T)
      @inbounds for j = 1:N
        num *= poly[j][s[j],i[j]]
        dom *= poly[j][s[j],i[j]]
      end

      numerator   += num
      denominator += dom^2

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights_extrema_threaded(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef,Tuple(order.+1))
  
  @inbounds @sync @qthreads for i in CartesianIndices(weights)

    numerator   = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      num = f[s]
      dom = one(T)
      @inbounds for j = 1:N
        pow = (sum(s[j]==1)+sum(s[j]==length(nodes[j])))
        num *= poly[j][s[j],i[j]]/(2^pow)
        dom *= poly[j][s[j],i[j]]*poly[j][s[j],i[j]]/(2^pow)
      end

      numerator   += num
      denominator += dom

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights_threaded(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i],normalize_node(nodes[i],domain[:,i]))
  end

  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
  
  weights = zeros(Tuple(ord.+1))
  
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N

      numerator   = zero(T)
      denominator = zero(T)
  
      @inbounds for s in CartesianIndices(f)
  
        num = f[s]
        dom = one(T)
        @inbounds for j = 1:N
          num *= poly[j][s[j],i[j]]
          dom *= poly[j][s[j],i[j]]
        end
  
        numerator   += num
        denominator += dom^2
  
      end
  
      weights[i] = numerator/denominator
      
    end
  
  end
  
  return weights
  
end

function chebyshev_weights_extrema_threaded(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i],normalize_node(nodes[i],domain[:,i]))
  end

  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
  
  weights = zeros(Tuple(ord.+1))
  
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N

      numerator   = zero(T)
      denominator = zero(T)
  
      @inbounds for s in CartesianIndices(f)
  
        num = f[s]
        dom = one(T)
        @inbounds for j = 1:N
          pow = (sum(s[j]==1)+sum(s[j]==length(nodes[j])))
          num *= poly[j][s[j],i[j]]/(2^pow)
          dom *= poly[j][s[j],i[j]]*poly[j][s[j],i[j]]/(2^pow)
        end
  
        numerator   += num
        denominator += dom
  
      end
  
      weights[i] = numerator/denominator
      
    end
  
  end
  
  return weights
  
end

function chebyshev_weights_threaded(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
    
  weights = zeros(Tuple(ord.+1))
  
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N

      numerator   = zero(T)
      denominator = zero(T)
  
      @inbounds for s in CartesianIndices(f)
  
        num = f[s]
        dom = one(T)
        @inbounds for j = 1:N
          num *= poly[j][s[j],i[j]]
          dom *= poly[j][s[j],i[j]]
        end
  
        numerator   += num
        denominator += dom^2
  
      end
  
      weights[i] = numerator/denominator
      
    end
  
  end
  
  return weights
  
end

function chebyshev_weights_extrema_threaded(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = (order,)
  for i = 2:N
    ord = (ord...,order)
  end
    
  weights = zeros(Tuple(ord.+1))
  
  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N

      numerator   = zero(T)
      denominator = zero(T)
  
      @inbounds for s in CartesianIndices(f)
  
        num = f[s]
        dom = one(T)
        @inbounds for j = 1:N
          pow = (sum(s[j]==1)+sum(s[j]==length(nodes[j])))
          num *= poly[j][s[j],i[j]]/(2^pow)
          dom *= poly[j][s[j],i[j]]*poly[j][s[j],i[j]]/(2^pow)
        end
  
        numerator   += num
        denominator += dom
  
      end
  
      weights[i] = numerator/denominator
      
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

# Functions for the one-variable case where the nodes are a vector

function chebyshev_weights(f::AbstractArray{T,N},nodes::Array{T,1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights(f,(nodes,),order,domain)
  
  return weights
  
end

function chebyshev_weights_extrema(f::AbstractArray{T,N},nodes::Array{T,1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extrema(f,(nodes,),order,domain)
  
  return weights
  
end

function chebyshev_weights(f::AbstractArray{T,N},poly::Array{T,2},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights(f,(poly,),order)
  
  return weights
  
end
  
function chebyshev_weights_extrema(f::AbstractArray{T,N},poly::Array{T,2},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extrema(f,(poly,),order)
  
  return weights
  
end

function chebyshev_weights(f::AbstractArray{T,N},nodes::Array{T,1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights(f,(nodes,),order,domain)
  
  return weights
  
end

function chebyshev_weights_extrema(f::AbstractArray{T,N},nodes::Array{T,1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extrema(f,(nodes,),order,domain)
  
  return weights
  
end

function chebyshev_weights(f::AbstractArray{T,N},poly::Array{T,2},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights(f,(poly,),order)
  
  return weights
  
end

function chebyshev_weights_extrema(f::AbstractArray{T,N},poly::Array{T,2},order::S) where {T<:AbstractFloat,N,S<:Integer}
  
  weights = chebyshev_weights_extrema(f,(poly,),order)
  
  return weights
  
end

# Functions that allow the nodes to be in an array of arrays

function chebyshev_weights(f::AbstractArray{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights(f,tuple(nodes...),order,domain)

  return weights

end

function chebyshev_weights_extrema(f::AbstractArray{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extrema(f,tuple(nodes...),order,domain)

  return weights

end

function chebyshev_weights(f::AbstractArray{T,N},poly::Array{Array{T,2},1},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights(f,tuple(poly...),order)

  return weights

end

function chebyshev_weights_extrema(f::AbstractArray{T,N},poly::Array{Array{T,2},1},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extrema(f,tuple(poly...),order)

  return weights

end

function chebyshev_weights(f::AbstractArray{T,N},nodes::Array{Array{T,1},1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights(f,tuple(nodes...),order,domain)

  return weights

end

function chebyshev_weights_extrema(f::AbstractArray{T,N},nodes::Array{Array{T,1},1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extrema(f,tuple(nodes...),order,domain)

  return weights

end

function chebyshev_weights(f::AbstractArray{T,N},poly::Array{Array{T,2},1},order::S) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights(f,tuple(poly...),order)

  return weights

end

function chebyshev_weights_extrema(f::AbstractArray{T,N},poly::Array{Array{T,2},1},order::S) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extrema(f,tuple(poly...),order)

  return weights

end

# Threaded functions

function chebyshev_weights_threaded(f::AbstractArray{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_threaded(f,tuple(nodes...),order,domain)

  return weights

end

function chebyshev_weights_extrema_threaded(f::AbstractArray{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extrema_threaded(f,tuple(nodes...),order,domain)

  return weights

end

function chebyshev_weights_threaded(f::AbstractArray{T,N},poly::Array{Array{T,2},1},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_threaded(f,tuple(poly...),order)

  return weights

end

function chebyshev_weights_extrema_threaded(f::AbstractArray{T,N},poly::Array{Array{T,2},1},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extrema_threaded(f,tuple(poly...),order)

  return weights

end

function chebyshev_weights_threaded(f::AbstractArray{T,N},nodes::Array{Array{T,1},1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_threaded(f,tuple(nodes...),order,domain)

  return weights

end

function chebyshev_weights_extrema_threaded(f::AbstractArray{T,N},nodes::Array{Array{T,1},1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extrema_threaded(f,tuple(nodes...),order,domain)

  return weights

end

function chebyshev_weights_threaded(f::AbstractArray{T,N},poly::Array{Array{T,2},1},order::S) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_threaded(f,tuple(poly...),order)

  return weights

end

function chebyshev_weights_extrema_threaded(f::AbstractArray{T,N},poly::Array{Array{T,2},1},order::S) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extrema_threaded(f,tuple(poly...),order)

  return weights

end