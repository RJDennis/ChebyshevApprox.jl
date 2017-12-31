# Differentiating tensor-product polynomials

function chebyshev_derivative(weights::Array{T,1},x::Array{T,1},order::Array{S,1},domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain)

  poly_deriv_1 = chebyshev_polynomial_derivative(order[1],[x1])

  evaluated_derivative = zero(T)

  for i = 1:order[1]+1

    evaluated_derivative += (2.0/(domain[1]-domain[2]))*poly_deriv_1[i]*weights[i]

  end

  return evaluated_derivative

end

function chebyshev_derivative(weights::Array{T,2},x::Array{T,1},order::Array{S,1},domain=[ones(1,2);-ones(1,2)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])

  polynomial_1_data = zeros(2,order[1]+1)
  polynomial_1_data[1,:] = chebyshev_polynomial(order[1],x1)
  polynomial_1_data[2,:] = chebyshev_polynomial_derivative(order[1],x1)

  polynomial_2_data = zeros(2,order[2]+1)
  polynomial_2_data[1,:] = chebyshev_polynomial(order[2],x2)
  polynomial_2_data[2,:] = chebyshev_polynomial_derivative(order[2],x2)

  evaluated_derivative = zeros(1,2)

  for j = 1:order[2]+1
    for i = 1:order[1]+1

      evaluated_derivative[1] += (2.0/(domain[1,1]-domain[2,1]))*polynomial_1_data[2,i]*polynomial_2_data[1,j]*weights[i,j]
      evaluated_derivative[2] += (2.0/(domain[1,2]-domain[2,2]))*polynomial_1_data[1,i]*polynomial_2_data[2,j]*weights[i,j]

    end
  end

  return evaluated_derivative

end

function chebyshev_derivative(weights::Array{T,3},x::Array{T,1},order::Array{S,1},domain=[ones(1,3);-ones(1,3)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])
  x3 = normalize_node(x[3],domain[:,3])

  polynomial_1_data = zeros(2,order[1]+1)
  polynomial_1_data[1,:] = chebyshev_polynomial(order[1],x1)
  polynomial_1_data[2,:] = chebyshev_polynomial_derivative(order[1],x1)

  polynomial_2_data = zeros(2,order[2]+1)
  polynomial_2_data[1,:] = chebyshev_polynomial(order[2],x2)
  polynomial_2_data[2,:] = chebyshev_polynomial_derivative(order[2],x2)

  polynomial_3_data = zeros(2,order[3]+1)
  polynomial_3_data[1,:] = chebyshev_polynomial(order[3],x3)
  polynomial_3_data[2,:] = chebyshev_polynomial_derivative(order[3],x3)

  evaluated_derivative = zeros(1,3)

  for k = 1:order[3]+1
    for j = 1:order[2]+1
      for i = 1:order[1]+1

        evaluated_derivative[1] += (2.0/(domain[1,1]-domain[2,1]))*polynomial_1_data[2,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*weights[i,j,k]
        evaluated_derivative[2] += (2.0/(domain[1,2]-domain[2,2]))*polynomial_1_data[1,i]*polynomial_2_data[2,j]*polynomial_3_data[1,k]*weights[i,j,k]
        evaluated_derivative[3] += (2.0/(domain[1,3]-domain[2,3]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[2,k]*weights[i,j,k]

      end
    end
  end

  return evaluated_derivative

end

function chebyshev_derivative(weights::Array{T,4},x::Array{T,1},order::Array{S,1},domain=[ones(1,4);-ones(1,4)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])
  x3 = normalize_node(x[3],domain[:,3])
  x4 = normalize_node(x[4],domain[:,4])

  polynomial_1_data = zeros(2,order[1]+1)
  polynomial_1_data[1,:] = chebyshev_polynomial(order[1],x1)
  polynomial_1_data[2,:] = chebyshev_polynomial_derivative(order[1],x1)

  polynomial_2_data = zeros(2,order[2]+1)
  polynomial_2_data[1,:] = chebyshev_polynomial(order[2],x2)
  polynomial_2_data[2,:] = chebyshev_polynomial_derivative(order[2],x2)

  polynomial_3_data = zeros(2,order[3]+1)
  polynomial_3_data[1,:] = chebyshev_polynomial(order[3],x3)
  polynomial_3_data[2,:] = chebyshev_polynomial_derivative(order[3],x3)

  polynomial_4_data = zeros(2,order[4]+1)
  polynomial_4_data[1,:] = chebyshev_polynomial(order[4],x4)
  polynomial_4_data[2,:] = chebyshev_polynomial_derivative(order[4],x4)

  evaluated_derivative = zeros(1,4)

  for l = 1:order[4]+1
    for k = 1:order[3]+1
      for j = 1:order[2]+1
        for i = 1:order[1]+1

          evaluated_derivative[1] += (2.0/(domain[1,1]-domain[2,1]))*polynomial_1_data[2,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*weights[i,j,k,l]
          evaluated_derivative[2] += (2.0/(domain[1,2]-domain[2,2]))*polynomial_1_data[1,i]*polynomial_2_data[2,j]*polynomial_3_data[2,k]*polynomial_4_data[2,l]*weights[i,j,k,l]
          evaluated_derivative[3] += (2.0/(domain[1,3]-domain[2,3]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[2,k]*polynomial_4_data[1,l]*weights[i,j,k,l]
          evaluated_derivative[4] += (2.0/(domain[1,4]-domain[2,4]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[2,l]*weights[i,j,k,l]

  			 end
			end
    end
  end

  return evaluated_derivative

end

function chebyshev_derivative(weights::Array{T,5},x::Array{T,1},order::Array{S,1},domain=[ones(1,5);-ones(1,5)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])
  x3 = normalize_node(x[3],domain[:,3])
  x4 = normalize_node(x[4],domain[:,4])
  x5 = normalize_node(x[5],domain[:,5])

  polynomial_1_data = zeros(2,order[1]+1)
  polynomial_1_data[1,:] = chebyshev_polynomial(order[1],x1)
  polynomial_1_data[2,:] = chebyshev_polynomial_derivative(order[1],x1)

  polynomial_2_data = zeros(2,order[2]+1)
  polynomial_2_data[1,:] = chebyshev_polynomial(order[2],x2)
  polynomial_2_data[2,:] = chebyshev_polynomial_derivative(order[2],x2)

  polynomial_3_data = zeros(2,order[3]+1)
  polynomial_3_data[1,:] = chebyshev_polynomial(order[3],x3)
  polynomial_3_data[2,:] = chebyshev_polynomial_derivative(order[3],x3)

  polynomial_4_data = zeros(2,order[4]+1)
  polynomial_4_data[1,:] = chebyshev_polynomial(order[4],x4)
  polynomial_4_data[2,:] = chebyshev_polynomial_derivative(order[4],x4)

  polynomial_5_data = zeros(2,order[5]+1)
  polynomial_5_data[1,:] = chebyshev_polynomial(order[5],x5)
  polynomial_5_data[2,:] = chebyshev_polynomial_derivative(order[5],x5)

  evaluated_derivative = zeros(1,5)

  for m = 1:order[5]+1
	  for l = 1:order[4]+1
      for k = 1:order[3]+1
        for j = 1:order[2]+1
          for i = 1:order[1]+1

            evaluated_derivative[1] += (2.0/(domain[1,1]-domain[2,1]))*polynomial_1_data[2,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*weights[i,j,k,l,m]
            evaluated_derivative[2] += (2.0/(domain[1,2]-domain[2,2]))*polynomial_1_data[1,i]*polynomial_2_data[2,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*weights[i,j,k,l,m]
            evaluated_derivative[3] += (2.0/(domain[1,3]-domain[2,3]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[2,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*weights[i,j,k,l,m]
            evaluated_derivative[4] += (2.0/(domain[1,4]-domain[2,4]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[2,l]*polynomial_5_data[1,m]*weights[i,j,k,l,m]
            evaluated_derivative[5] += (2.0/(domain[1,5]-domain[2,5]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[2,m]*weights[i,j,k,l,m]

					end
				end
      end
    end
  end

  return evaluated_derivative

end

function chebyshev_derivative(weights::Array{T,6},x::Array{T,1},order::Array{S,1},domain=[ones(1,6);-ones(1,6)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])
  x3 = normalize_node(x[3],domain[:,3])
  x4 = normalize_node(x[4],domain[:,4])
  x5 = normalize_node(x[5],domain[:,5])
  x6 = normalize_node(x[6],domain[:,6])

  polynomial_1_data = zeros(2,order[1]+1)
  polynomial_1_data[1,:] = chebyshev_polynomial(order[1],x1)
  polynomial_1_data[2,:] = chebyshev_polynomial_derivative(order[1],x1)

  polynomial_2_data = zeros(2,order[2]+1)
  polynomial_2_data[1,:] = chebyshev_polynomial(order[2],x2)
  polynomial_2_data[2,:] = chebyshev_polynomial_derivative(order[2],x2)

  polynomial_3_data = zeros(2,order[3]+1)
  polynomial_3_data[1,:] = chebyshev_polynomial(order[3],x3)
  polynomial_3_data[2,:] = chebyshev_polynomial_derivative(order[3],x3)

  polynomial_4_data = zeros(2,order[4]+1)
  polynomial_4_data[1,:] = chebyshev_polynomial(order[4],x4)
  polynomial_4_data[2,:] = chebyshev_polynomial_derivative(order[4],x4)

  polynomial_5_data = zeros(2,order[5]+1)
  polynomial_5_data[1,:] = chebyshev_polynomial(order[5],x5)
  polynomial_5_data[2,:] = chebyshev_polynomial_derivative(order[5],x5)

  polynomial_6_data = zeros(2,order[6]+1)
  polynomial_6_data[1,:] = chebyshev_polynomial(order[6],x6)
  polynomial_6_data[2,:] = chebyshev_polynomial_derivative(order[6],x6)

  evaluated_derivative = zeros(1,6)

  for n = 1:order[6]+1
	  for m = 1:order[5]+1
	    for l = 1:order[4]+1
        for k = 1:order[3]+1
          for j = 1:order[2]+1
            for i = 1:order[1]+1

              evaluated_derivative[1] += (2.0/(domain[1,1]-domain[2,1]))*polynomial_1_data[2,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*polynomial_6_data[1,n]*weights[i,j,k,l,m,n]
              evaluated_derivative[2] += (2.0/(domain[1,2]-domain[2,2]))*polynomial_1_data[1,i]*polynomial_2_data[2,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*polynomial_6_data[1,n]*weights[i,j,k,l,m,n]
              evaluated_derivative[3] += (2.0/(domain[1,3]-domain[2,3]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[2,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*polynomial_6_data[1,n]*weights[i,j,k,l,m,n]
              evaluated_derivative[4] += (2.0/(domain[1,4]-domain[2,4]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[2,l]*polynomial_5_data[1,m]*polynomial_6_data[1,n]*weights[i,j,k,l,m,n]
              evaluated_derivative[5] += (2.0/(domain[1,5]-domain[2,5]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[2,m]*polynomial_6_data[1,n]*weights[i,j,k,l,m,n]
              evaluated_derivative[6] += (2.0/(domain[1,6]-domain[2,6]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*polynomial_6_data[2,n]*weights[i,j,k,l,m,n]

            end
          end
				end
      end
    end
  end

  return evaluated_derivative

end

# Differentiating complete polynomials

function chebyshev_derivative(weights::Array{T,1},x::Array{T,1},order::S,domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain)

  poly_deriv_1 = chebyshev_polynomial_derivative(order,[x1])

  evaluated_derivative = zero(T)

  for i = 1:order+1

    evaluated_derivative += (2.0/(domain[1,1]-domain[2,1]))*poly_deriv_1[i]*weights[i]

  end

  return evaluated_derivative

end

function chebyshev_derivative(weights::Array{T,2},x::Array{T,1},order::S,domain=[ones(1,2);-ones(1,2)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])

  polynomial_1_data = zeros(2,order+1)
  polynomial_1_data[1,:] = chebyshev_polynomial(order,x1)
  polynomial_1_data[2,:] = chebyshev_polynomial_derivative(order,x1)

  polynomial_2_data = zeros(2,order+1)
  polynomial_2_data[1,:] = chebyshev_polynomial(order,x2)
  polynomial_2_data[2,:] = chebyshev_polynomial_derivative(order,x2)

  evaluated_derivative = zeros(1,2)

  for j = 1:order+1
    for i = 1:order+1

      if (i+j <= order+2)

        evaluated_derivative[1] += (2.0/(domain[1,1]-domain[2,1]))*polynomial_1_data[2,i]*polynomial_2_data[1,j]*weights[i,j]
        evaluated_derivative[2] += (2.0/(domain[1,2]-domain[2,2]))*polynomial_1_data[1,i]*polynomial_2_data[2,j]*weights[i,j]

      end

    end
  end

  return evaluated_derivative

end

function chebyshev_derivative(weights::Array{T,3},x::Array{T,1},order::S,domain=[ones(1,3);-ones(1,3)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])
  x3 = normalize_node(x[3],domain[:,3])

  polynomial_1_data = zeros(2,order+1)
  polynomial_1_data[1,:] = chebyshev_polynomial(order,x1)
  polynomial_1_data[2,:] = chebyshev_polynomial_derivative(order,x1)

  polynomial_2_data = zeros(2,order+1)
  polynomial_2_data[1,:] = chebyshev_polynomial(order,x2)
  polynomial_2_data[2,:] = chebyshev_polynomial_derivative(order,x2)

  polynomial_3_data = zeros(2,order+1)
  polynomial_3_data[1,:] = chebyshev_polynomial(order,x3)
  polynomial_3_data[2,:] = chebyshev_polynomial_derivative(order,x3)

  evaluated_derivative = zeros(1,3)

  for k = 1:order+1
    for j = 1:order+1
      for i = 1:order+1

        if (i+j+k <= order+3)

          evaluated_derivative[1] += (2.0/(domain[1,1]-domain[2,1]))*polynomial_1_data[2,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*weights[i,j,k]
          evaluated_derivative[2] += (2.0/(domain[1,2]-domain[2,2]))*polynomial_1_data[1,i]*polynomial_2_data[2,j]*polynomial_3_data[1,k]*weights[i,j,k]
          evaluated_derivative[3] += (2.0/(domain[1,3]-domain[2,3]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[2,k]*weights[i,j,k]

        end

      end
    end
  end

  return evaluated_derivative

end

function chebyshev_derivative(weights::Array{T,4},x::Array{T,1},order::S,domain=[ones(1,4);-ones(1,4)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])
  x3 = normalize_node(x[3],domain[:,3])
  x4 = normalize_node(x[4],domain[:,4])

  polynomial_1_data = zeros(2,order+1)
  polynomial_1_data[1,:] = chebyshev_polynomial(order,x1)
  polynomial_1_data[2,:] = chebyshev_polynomial_derivative(order,x1)

  polynomial_2_data = zeros(2,order+1)
  polynomial_2_data[1,:] = chebyshev_polynomial(order,x2)
  polynomial_2_data[2,:] = chebyshev_polynomial_derivative(order,x2)

  polynomial_3_data = zeros(2,order+1)
  polynomial_3_data[1,:] = chebyshev_polynomial(order,x3)
  polynomial_3_data[2,:] = chebyshev_polynomial_derivative(order,x3)

  polynomial_4_data = zeros(2,order+1)
  polynomial_4_data[1,:] = chebyshev_polynomial(order,x4)
  polynomial_4_data[2,:] = chebyshev_polynomial_derivative(order,x4)

  evaluated_derivative = zeros(1,4)

  for l = 1:order+1
    for k = 1:order+1
      for j = 1:order+1
        for i = 1:order+1

          if (i+j+k+l <= order+4)

            evaluated_derivative[1] += (2.0/(domain[1,1]-domain[2,1]))*polynomial_1_data[2,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*weights[i,j,k,l]
            evaluated_derivative[2] += (2.0/(domain[1,2]-domain[2,2]))*polynomial_1_data[1,i]*polynomial_2_data[2,j]*polynomial_3_data[2,k]*polynomial_4_data[2,l]*weights[i,j,k,l]
            evaluated_derivative[3] += (2.0/(domain[1,3]-domain[2,3]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[2,k]*polynomial_4_data[1,l]*weights[i,j,k,l]
            evaluated_derivative[4] += (2.0/(domain[1,4]-domain[2,4]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[2,l]*weights[i,j,k,l]

          end

        end
			end
    end
  end

  return evaluated_derivative

end

function chebyshev_derivative(weights::Array{T,5},x::Array{T,1},order::S,domain=[ones(1,5);-ones(1,5)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])
  x3 = normalize_node(x[3],domain[:,3])
  x4 = normalize_node(x[4],domain[:,4])
  x5 = normalize_node(x[5],domain[:,5])

  polynomial_1_data = zeros(2,order+1)
  polynomial_1_data[1,:] = chebyshev_polynomial(order,x1)
  polynomial_1_data[2,:] = chebyshev_polynomial_derivative(order,x1)

  polynomial_2_data = zeros(2,order+1)
  polynomial_2_data[1,:] = chebyshev_polynomial(order,x2)
  polynomial_2_data[2,:] = chebyshev_polynomial_derivative(order,x2)

  polynomial_3_data = zeros(2,order+1)
  polynomial_3_data[1,:] = chebyshev_polynomial(order,x3)
  polynomial_3_data[2,:] = chebyshev_polynomial_derivative(order,x3)

  polynomial_4_data = zeros(2,order+1)
  polynomial_4_data[1,:] = chebyshev_polynomial(order,x4)
  polynomial_4_data[2,:] = chebyshev_polynomial_derivative(order,x4)

  polynomial_5_data = zeros(2,order+1)
  polynomial_5_data[1,:] = chebyshev_polynomial(order,x5)
  polynomial_5_data[2,:] = chebyshev_polynomial_derivative(order,x5)

  evaluated_derivative = zeros(1,5)

  for m = 1:order+1
	  for l = 1:order+1
      for k = 1:order+1
        for j = 1:order+1
          for i = 1:order+1

            if (i+j+k+l+m <= order+5)

              evaluated_derivative[1] += (2.0/(domain[1,1]-domain[2,1]))*polynomial_1_data[2,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*weights[i,j,k,l,m]
              evaluated_derivative[2] += (2.0/(domain[1,2]-domain[2,2]))*polynomial_1_data[1,i]*polynomial_2_data[2,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*weights[i,j,k,l,m]
              evaluated_derivative[3] += (2.0/(domain[1,3]-domain[2,3]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[2,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*weights[i,j,k,l,m]
              evaluated_derivative[4] += (2.0/(domain[1,4]-domain[2,4]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[2,l]*polynomial_5_data[1,m]*weights[i,j,k,l,m]
              evaluated_derivative[5] += (2.0/(domain[1,5]-domain[2,5]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[2,m]*weights[i,j,k,l,m]

            end

					end
				end
      end
    end
  end

  return evaluated_derivative

end

function chebyshev_derivative(weights::Array{T,6},x::Array{T,1},order::S,domain=[ones(1,6);-ones(1,6)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])
  x3 = normalize_node(x[3],domain[:,3])
  x4 = normalize_node(x[4],domain[:,4])
  x5 = normalize_node(x[5],domain[:,5])
  x6 = normalize_node(x[6],domain[:,6])

  polynomial_1_data = zeros(2,order+1)
  polynomial_1_data[1,:] = chebyshev_polynomial(order,x1)
  polynomial_1_data[2,:] = chebyshev_polynomial_derivative(order,x1)

  polynomial_2_data = zeros(2,order+1)
  polynomial_2_data[1,:] = chebyshev_polynomial(order,x2)
  polynomial_2_data[2,:] = chebyshev_polynomial_derivative(order,x2)

  polynomial_3_data = zeros(2,order+1)
  polynomial_3_data[1,:] = chebyshev_polynomial(order,x3)
  polynomial_3_data[2,:] = chebyshev_polynomial_derivative(order,x3)

  polynomial_4_data = zeros(2,order+1)
  polynomial_4_data[1,:] = chebyshev_polynomial(order,x4)
  polynomial_4_data[2,:] = chebyshev_polynomial_derivative(order,x4)

  polynomial_5_data = zeros(2,order+1)
  polynomial_5_data[1,:] = chebyshev_polynomial(order,x5)
  polynomial_5_data[2,:] = chebyshev_polynomial_derivative(order,x5)

  polynomial_6_data = zeros(2,order+1)
  polynomial_6_data[1,:] = chebyshev_polynomial(order,x6)
  polynomial_6_data[2,:] = chebyshev_polynomial_derivative(order,x6)

  evaluated_derivative = zeros(1,6)

  for n = 1:order+1
	  for m = 1:order+1
	    for l = 1:order+1
        for k = 1:order+1
          for j = 1:order+1
            for i = 1:order+1

              if (i+j+k+l+m+n <= order+6)

                evaluated_derivative[1] += (2.0/(domain[1,1]-domain[2,1]))*polynomial_1_data[2,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*polynomial_6_data[1,n]*weights[i,j,k,l,m,n]
                evaluated_derivative[2] += (2.0/(domain[1,2]-domain[2,2]))*polynomial_1_data[1,i]*polynomial_2_data[2,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*polynomial_6_data[1,n]*weights[i,j,k,l,m,n]
                evaluated_derivative[3] += (2.0/(domain[1,3]-domain[2,3]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[2,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*polynomial_6_data[1,n]*weights[i,j,k,l,m,n]
                evaluated_derivative[4] += (2.0/(domain[1,4]-domain[2,4]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[2,l]*polynomial_5_data[1,m]*polynomial_6_data[1,n]*weights[i,j,k,l,m,n]
                evaluated_derivative[5] += (2.0/(domain[1,5]-domain[2,5]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[2,m]*polynomial_6_data[1,n]*weights[i,j,k,l,m,n]
                evaluated_derivative[6] += (2.0/(domain[1,6]-domain[2,6]))*polynomial_1_data[1,i]*polynomial_2_data[1,j]*polynomial_3_data[1,k]*polynomial_4_data[1,l]*polynomial_5_data[1,m]*polynomial_6_data[2,n]*weights[i,j,k,l,m,n]

              end
              
            end
          end
				end
      end
    end
  end

  return evaluated_derivative

end
