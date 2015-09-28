function chebyshev_polynomial{T<:AbstractFloat,S<:Integer}(order::S,x::T)

  polynomial = Array(T,1,order+1)

  for i = 1:order+1

    if i == 1

      polynomial[i] = one(T)

    elseif i == 2

      polynomial[i] = x

    elseif i == 3

      polynomial[i] = 2*x*x-one(T)

    else

      polynomial[i] = 2*x*polynomial[i-1]-polynomial[i-2]

    end

  end

  return polynomial

end

function chebyshev_polynomial{T<:AbstractFloat,S<:Integer}(order::S,x::Array{T,1})

  polynomial = Array(T,length(x),order+1)

  for i = 1:order+1
	  for j = 1:length(x)

      if i == 1

        polynomial[j,i] = one(T)

      elseif i == 2

        polynomial[j,i] = x[j]

      elseif i == 3

        polynomial[j,i] = 2*x[j]*x[j]-one(T)

      else

        polynomial[j,i] = 2*x[j]*polynomial[j,i-1]-polynomial[j,i-2]

      end

    end
  end

  return polynomial

end
