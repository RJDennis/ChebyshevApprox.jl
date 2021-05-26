# Regular functions for evaluating Chebyshev polynominals

function chebyshev_evaluate(weights::Array{T,N},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)
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

function chebyshev_evaluate(weights::Array{T,N},x::Array{T,1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)
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

function chebyshev_evaluate(cheb_poly::ChebPoly,x::Array{T,1}) where {T<:AbstractFloat}

  yhat = chebyshev_evaluate(cheb_poly.weights,x,cheb_poly.order,cheb_poly.domain) 
    
  return yhat
  
end

function cheb_interp(cheb::ChebInterpRoots,x::Array{T,1}) where {T<:AbstractFloat}
  
  weights = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  yhat    = chebyshev_evaluate(weights,x,cheb.order,cheb.domain) 
     
  return yhat
    
end
  
function cheb_interp(cheb::ChebInterpExtrema,x::Array{T,1}) where {T<:AbstractFloat}
  
  weights = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  yhat    = chebyshev_evaluate(weights,x,cheb.order,cheb.domain) 
      
  return yhat
    
end
  
function cheb_interp(cheb::ChebInterpExtended,x::Array{T,1}) where {T<:AbstractFloat}
  
  weights = chebyshev_weights_extended_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  yhat    = chebyshev_evaluate(weights,x,cheb.order,cheb.domain) 
      
  return yhat
    
end

function cheb_interp(cheb::ChebInterpVertesi,x::Array{T,1}) where {T<:AbstractFloat}
  
  weights = chebyshev_weights_vertesi_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
  yhat    = chebyshev_evaluate(weights,x,cheb.order,cheb.domain) 
        
  return yhat
      
end
  
function chebyshev_evaluate(weights::Array{T,N},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}
  
  function chebeval(x::Array{T,1}) where {T <: AbstractFloat}
  
    return chebyshev_evaluate(weights,x,order,domain)
  
  end
  
  return chebeval
  
end
  
function chebyshev_evaluate(weights::Array{T,N},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}
  
  function chebeval(x::Array{T,1}) where {T <: AbstractFloat}
  
    return chebyshev_evaluate(weights,x,order,domain)
  
  end
  
  return chebeval
  
end
  
function chebyshev_evaluate(cheb_poly::ChebPoly)
  
  function chebeval(x::Array{T,1}) where {T <: AbstractFloat}
  
    return chebyshev_evaluate(cheb_poly,x)
  
  end
  
  return chebeval

end
  
function cheb_interp(cheb::ChebInterpRoots)
 
  weights = chebyshev_weights_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
    
  function chebeval(x::Array{T,1}) where {T <: AbstractFloat}
  
    return chebyshev_evaluate(weights,x,cheb.order,cheb.domain)
  
  end
  
  return chebeval
  
end
  
function cheb_interp(cheb::ChebInterpExtrema)
  
  weights = chebyshev_weights_extrema_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
    
  function chebeval(x::Array{T,1}) where {T <: AbstractFloat}
  
    return chebyshev_evaluate(weights,x,cheb.order,cheb.domain)
  
  end
  
  return chebeval
  
end
  
function cheb_interp(cheb::ChebInterpExtended)
  
  weights = chebyshev_weights_extended_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
    
  function chebeval(x::Array{T,1}) where {T <: AbstractFloat}
  
    return chebyshev_evaluate(weights,x,cheb.order,cheb.domain)
  
  end
  
  return chebeval
  
end

function cheb_interp(cheb::ChebInterpVertesi)
  
  weights = chebyshev_weights_vertesi_threaded(cheb.data,cheb.nodes,cheb.order,cheb.domain)
      
  function chebeval(x::Array{T,1}) where {T <: AbstractFloat}
    
    return chebyshev_evaluate(weights,x,cheb.order,cheb.domain)
    
  end
    
  return chebeval
    
end
