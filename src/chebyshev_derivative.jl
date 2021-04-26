# Regular functions for evaluating derivatives and gradients of Chebyshev polynomials

# Regular functions for derivatives

function chebyshev_derivative(weights::Array{T,N},x::Array{T,1},pos::S,order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)
  @inbounds for i = 1:N
    if i === pos
      poly[i] = derivative_of_chebyshev_polynomial(order[i],normalize_node(x[i],domain[:,i]))
    else
      poly[i] = chebyshev_polynomial(order[i],normalize_node(x[i],domain[:,i]))
    end
  end

  derivative = 0.0
  @inbounds for i in CartesianIndices(weights)
    poly_product = poly[1][i[1]]
    @inbounds for j = 2:N
      poly_product *= poly[j][i[j]]
    end
    derivative += weights[i]*poly_product
  end

  return derivative*(2.0/(domain[1,pos]-domain[2,pos]))

end

function chebyshev_derivative(weights::Array{T,N},x::Array{T,1},pos::S,order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)
  @inbounds for i = 1:N
    if i === pos
      poly[i] = derivative_of_chebyshev_polynomial(order,normalize_node(x[i],domain[:,i]))
    else
      poly[i] = chebyshev_polynomial(order,normalize_node(x[i],domain[:,i]))
    end
  end

  derivative = 0.0
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

function chebyshev_derivative(cheb_poly::ChebPoly,x::Array{T,1},pos::S) where {T<:AbstractFloat,S<:Integer}

  derivative = chebyshev_derivative(cheb_poly.weights,x,pos,cheb_poly.order,cheb_poly.domain) 
    
  return derivative

end

function chebyshev_derivative(cheb::ChebInterpRoots,x::Array{T,1},pos::S) where {T<:AbstractFloat,S<:Integer}

  weights    = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  derivative = chebyshev_derivative(weights,x,pos,cheb.order,cheb.domain) 
    
  return derivative

end

function chebyshev_derivative(cheb::ChebInterpExtrema,x::Array{T,1},pos::S) where {T<:AbstractFloat,S<:Integer}

  weights    = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  derivative = chebyshev_derivative(weights,x,pos,cheb.order,cheb.domain) 
    
  return derivative

end

# Regular functions for gradients

function chebyshev_gradient(weights::Array{T,N},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  gradient = Array{T,2}(undef,1,N)

  @inbounds for i = 1:N
    gradient[i] = chebyshev_derivative(weights,x,i,order,domain)
  end

  return gradient

end

function chebyshev_gradient(weights::Array{T,N},x::Array{T,1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  gradient = Array{T,2}(undef,1,N)

  @inbounds for i = 1:N
    gradient[i] = chebyshev_derivative(weights,x,i,order,domain)
  end

  return gradient

end

function chebyshev_gradient(cheb_poly::ChebPoly,x::Array{T,1}) where {T<:AbstractFloat}

  gradient = chebyshev_gradient(cheb_poly.weights,x,cheb_poly.order,cheb_poly.domain) 
    
  return gradient

end

function chebyshev_gradient(cheb::ChebInterpRoots,x::Array{T,1}) where {T<:AbstractFloat}

  weights  = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  gradient = chebyshev_gradient(weights,x,cheb.order,cheb.domain) 
    
  return gradient

end

function chebyshev_gradient(cheb::ChebInterpExtrema,x::Array{T,1}) where {T<:AbstractFloat}

  weights  = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  gradient = chebyshev_gradient(weights,x,cheb.order,cheb.domain) 
    
  return gradient

end

function chebyshev_derivative(weights::Array{T,N},pos::S,order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function chebderiv(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_derivative(weights,x,pos,order,domain)

  end

  return chebderiv

end

function chebyshev_derivative(weights::Array{T,N},pos::S,order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function chebderiv(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_derivative(weights,x,pos,order,domain)

  end

  return chebderiv

end

function chebyshev_derivative(cheb_poly::ChebPoly,pos::S) where {S <: Integer}

  function chebderiv(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_derivative(cheb_poly,x,pos)

  end

  return chebderiv

end

function chebyshev_derivative(cheb::ChebInterpRoots,pos::S) where {S <: Integer}

  weights = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)

  function chebderiv(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_derivative(weights,x,pos,cheb.order,cheb.domain)

  end

  return chebderiv

end

function chebyshev_derivative(cheb::ChebInterpExtrema,pos::S) where {S <: Integer}

  weights = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)

  function chebderiv(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_derivative(weights,x,pos,cheb.order,cheb.domain)

  end

  return chebderiv

end

function chebyshev_gradient(weights::Array{T,N},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function chebgrad(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_gradient(weights,x,pos,order,domain)

  end

  return chebgrad

end

function chebyshev_gradient(weights::Array{T,N},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function chebgrad(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_gradient(weights,x,pos,order,domain)

  end

  return chebgrad

end

function chebyshev_gradient(cheb_poly::ChebPoly) where {S <: Integer}

  function chebgrad(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_gradient(cheb_poly,x)

  end

  return chebgrad

end

function chebyshev_gradient(cheb::ChebInterpRoots) where {S <: Integer}

  weights = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)

  function chebgrad(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_gradient(weights,x,cheb.order,cheb.domain)

  end

  return chebgrad

end

function chebyshev_gradient(cheb::ChebInterpExtrema) where {S <: Integer}

  weights = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)

  function chebgrad(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_gradient(weights,x,cheb.order,cheb.domain)

  end

  return chebgrad

end