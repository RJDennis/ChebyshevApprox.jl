function chebyshev_polynomial_derivative(order::S,x::T) where {T<:AbstractFloat,S<:Integer}

  polynomial = ones(T,1,order+1)
  poly_deriv = zeros(T,1,order+1)

  for i = 2:order+1
    if i == 2
      polynomial[i] = x
      poly_deriv[i] = 1.0
    else
      polynomial[i] = 2*x*polynomial[i-1]-polynomial[i-2]
      poly_deriv[i] = ((i-1)*polynomial[i-1]-(i-1)*x*polynomial[i])/(1-x^2)
    end
  end

  return poly_deriv

end

function chebyshev_polynomial_derivative(order::S,x::Array{T,1}) where {T<:AbstractFloat,S<:Integer}

  polynomial = ones(T,length(x),order+1)
  poly_deriv = zeros(T,length(x),order+1)

  for i = 1:length(x)
    poly_deriv[i,:] = chebyshev_polynomial_derivative(order,x[i])
  end

#=

  for i = 2:order+1
	for j = 1:length(x)
      if i == 2
        polynomial[j,i] = x[j]
        poly_deriv[j,i] = 1.0
      else
        polynomial[j,i] = 2*x[j]*polynomial[j,i-1]-polynomial[j,i-2]
        poly_deriv[j,i] = ((i-1)*polynomial[j,i-1]-(i-1)*x[j]*polynomial[j,i])/(1-x[j]^2)
      end
    end
  end

=#

  return poly_deriv

end
