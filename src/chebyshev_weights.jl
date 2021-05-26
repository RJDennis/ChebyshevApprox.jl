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

function chebyshev_weights_extrema(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extended(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extrema(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extended(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extrema(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extended(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extrema(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extended(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extrema_threaded(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extended_threaded(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_threaded(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extrema_threaded(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extended_threaded(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_threaded(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extrema_threaded(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extended_threaded(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_threaded(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extrema_threaded(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights_extended_threaded(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T<:AbstractFloat,N,S<:Integer}

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

function chebyshev_weights(f::AbstractArray{T,1},nodes::Array{T,1},order::Union{S,Array{S,1}},domain=[one(T);-one(T)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes,),order,domain)
  
  return weights
  
end

function chebyshev_weights_extrema(f::AbstractArray{T,1},nodes::Array{T,1},order::Union{S,Array{S,1}},domain=[one(T);-one(T)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights_extrema(f,(nodes,),order,domain)
  
  return weights
  
end

function chebyshev_weights_extended(f::AbstractArray{T,1},nodes::Array{T,1},order::Union{S,Array{S,1}},domain=[one(T);-one(T)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights_extended(f,(nodes,),order,domain)
  
  return weights
  
end

function chebyshev_weights(f::AbstractArray{T,1},poly::Array{T,2},order::Union{S,Array{S,1}}) where {T<:AbstractFloat,S<:Integer}
  
  weights = chebyshev_weights(f,(poly,),order)
  
  return weights
  
end
  
function chebyshev_weights_extrema(f::AbstractArray{T,1},poly::Array{T,2},order::Union{S,Array{S,1}}) where {T<:AbstractFloat,S<:Integer}
  
  weights = chebyshev_weights_extrema(f,(poly,),order)
  
  return weights
  
end

function chebyshev_weights_extended(f::AbstractArray{T,1},poly::Array{T,2},order::Union{S,Array{S,1}}) where {T<:AbstractFloat,S<:Integer}
  
  weights = chebyshev_weights_extended(f,(poly,),order)
  
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

function chebyshev_weights_extended(f::AbstractArray{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extended(f,tuple(nodes...),order,domain)

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

function chebyshev_weights_extended(f::AbstractArray{T,N},poly::Array{Array{T,2},1},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extended(f,tuple(poly...),order)

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

function chebyshev_weights_extended(f::AbstractArray{T,N},nodes::Array{Array{T,1},1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extended(f,tuple(nodes...),order,domain)

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

function chebyshev_weights_extended(f::AbstractArray{T,N},poly::Array{Array{T,2},1},order::S) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extended(f,tuple(poly...),order)

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

function chebyshev_weights_extended_threaded(f::AbstractArray{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extended_threaded(f,tuple(nodes...),order,domain)

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

function chebyshev_weights_extended_threaded(f::AbstractArray{T,N},poly::Array{Array{T,2},1},order::Array{S,1}) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extended_threaded(f,tuple(poly...),order)

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

function chebyshev_weights_extended_threaded(f::AbstractArray{T,N},nodes::Array{Array{T,1},1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extended_threaded(f,tuple(nodes...),order,domain)

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

function chebyshev_weights_extended_threaded(f::AbstractArray{T,N},poly::Array{Array{T,2},1},order::S) where {T<:AbstractFloat,N,S<:Integer}

  weights = chebyshev_weights_extended_threaded(f,tuple(poly...),order)

  return weights

end

const chebyshev_weights_vertesi          = chebyshev_weights_extended
const chebyshev_weights_vertesi_threaded = chebyshev_weights_extended_threaded
