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
    for j = 1:length(x)
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
