function chebyshev_polynomial{T<:AbstractFloat,S<:Integer}(order::S,x::T)

  polynomial    = Array(T,1,order+1)
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

function chebyshev_polynomial{T<:AbstractFloat,S<:Integer}(order::S,x::Array{T,1})

  polynomial      = Array(T,length(x),order+1)
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
