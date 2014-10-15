function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,1},x::Array{T,1},order::Array{S,1},range::Array{T,1})

  x1 = normalize_node(x[1],range)

  polynomial_1 = chebyshev_polynomial(order[1],[x1])

  evaluated_polynomial = zero(T)

  for i = 1:order[1]+1

    evaluated_polynomial += polynomial_1[i]*weights[i]

  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,2},x::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])

  polynomial_1 = chebyshev_polynomial(order[1],[x1])
  polynomial_2 = chebyshev_polynomial(order[2],[x2])

  evaluated_polynomial = zero(T)

  for j = 1:order[2]+1
    for i = 1:order[1]+1

      evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*weights[i,j]

    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,3},x::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])

  polynomial_1 = chebyshev_polynomial(order[1],[x1])
  polynomial_2 = chebyshev_polynomial(order[2],[x2])
  polynomial_3 = chebyshev_polynomial(order[3],[x3])

  evaluated_polynomial = zero(T)

  for k = 1:order[3]+1
    for j = 1:order[2]+1
      for i = 1:order[1]+1

        evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*polynomial_3[k]*weights[i,j,k]

      end
    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,4},x::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])
  x4 = normalize_node(x[4],range[:,4])

  polynomial_1 = chebyshev_polynomial(order[1],[x1])
  polynomial_2 = chebyshev_polynomial(order[2],[x2])
  polynomial_3 = chebyshev_polynomial(order[3],[x3])
  polynomial_4 = chebyshev_polynomial(order[4],[x4])

  evaluated_polynomial = zero(T)

  for l = 1:order[4]+1
    for k = 1:order[3]+1
      for j = 1:order[2]+1
        for i = 1:order[1]+1

          evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*polynomial_3[k]*polynomial_4[l]*weights[i,j,k,l]

  			 end
			end
    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,5},x::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])
  x4 = normalize_node(x[4],range[:,4])
  x5 = normalize_node(x[5],range[:,5])

  polynomial_1 = chebyshev_polynomial(order[1],[x1])
  polynomial_2 = chebyshev_polynomial(order[2],[x2])
  polynomial_3 = chebyshev_polynomial(order[3],[x3])
  polynomial_4 = chebyshev_polynomial(order[4],[x4])
  polynomial_5 = chebyshev_polynomial(order[5],[x5])

  evaluated_polynomial = zero(T)

	for m = 1:order[5]+1
	  for l = 1:order[4]+1
      for k = 1:order[3]+1
        for j = 1:order[2]+1
          for i = 1:order[1]+1

            evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*polynomial_3[k]*polynomial_4[l]*polynomial_5[m]*weights[i,j,k,l,m]

					end
				end
      end
    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,1},x::Array{T,1},order::S,range::Array{T,1})

  x1 = normalize_node(x[1],range)

  polynomial_1 = chebyshev_polynomial(order,[x1])

  evaluated_polynomial = zero(T)

  for i = 1:order+1

    evaluated_polynomial += polynomial_1[i]*weights[i]

  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,2},x::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])

  polynomial_1 = chebyshev_polynomial(order,[x1])
  polynomial_2 = chebyshev_polynomial(order,[x2])

  evaluated_polynomial = zero(T)

  for j = 1:order+1
    for i = 1:order+1

      if (i+j <= order+2)

        evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*weights[i,j]

      end

    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,3},x::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])

  polynomial_1 = chebyshev_polynomial(order,[x1])
  polynomial_2 = chebyshev_polynomial(order,[x2])
  polynomial_3 = chebyshev_polynomial(order,[x3])

  evaluated_polynomial = zero(T)

  for k = 1:order+1
    for j = 1:order+1
      for i = 1:order+1

        if (i+j+k <= order+3)

          evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*polynomial_3[k]*weights[i,j,k]

        end

      end
    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,4},x::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])
  x4 = normalize_node(x[4],range[:,4])

  polynomial_1 = chebyshev_polynomial(order,[x1])
  polynomial_2 = chebyshev_polynomial(order,[x2])
  polynomial_3 = chebyshev_polynomial(order,[x3])
  polynomial_4 = chebyshev_polynomial(order,[x4])

  evaluated_polynomial = zero(T)

  for l = 1:order+1
    for k = 1:order+1
      for j = 1:order+1
        for i = 1:order+1

          if (i+j+k+l <= order+4)

            evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*polynomial_3[k]*polynomial_4[l]*weights[i,j,k,l]

          end

  			end
			end
    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,5},x::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])
  x4 = normalize_node(x[4],range[:,4])
  x5 = normalize_node(x[5],range[:,5])

  polynomial_1 = chebyshev_polynomial(order,[x1])
  polynomial_2 = chebyshev_polynomial(order,[x2])
  polynomial_3 = chebyshev_polynomial(order,[x3])
  polynomial_4 = chebyshev_polynomial(order,[x4])
  polynomial_5 = chebyshev_polynomial(order,[x5])

  evaluated_polynomial = zero(T)

  for m = 1:order+1
    for l = 1:order+1
      for k = 1:order+1
	      for j = 1:order+1
          for i = 1:order+1

            if (i+j+k+l+m <= order+5)

              evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*polynomial_3[k]*polynomial_4[l]*polynomial_5[m]*weights[i,j,k,l,m]

            end

					end
				end
      end
    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,1},x::Array{T,1},order::Array{S,1})

  polynomial_1 = chebyshev_polynomial(order[1],[x[1]])

  evaluated_polynomial = zero(T)

  for i = 1:order[1]+1

    evaluated_polynomial += polynomial_1[i]*weights[i]

  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,2},x::Array{T,1},order::Array{S,1})

  polynomial_1 = chebyshev_polynomial(order[1],[x[1]])
  polynomial_2 = chebyshev_polynomial(order[2],[x[2]])

  evaluated_polynomial = zero(T)

  for j = 1:order[2]+1
    for i = 1:order[1]+1

      evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*weights[i,j]

    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,3},x::Array{T,1},order::Array{S,1})

  polynomial_1 = chebyshev_polynomial(order[1],[x[1]])
  polynomial_2 = chebyshev_polynomial(order[2],[x[2]])
  polynomial_3 = chebyshev_polynomial(order[3],[x[3]])

  evaluated_polynomial = zero(T)

  for k = 1:order[3]+1
    for j = 1:order[2]+1
      for i = 1:order[1]+1

        evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*polynomial_3[k]*weights[i,j,k]

      end
    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,4},x::Array{T,1},order::Array{S,1})

  polynomial_1 = chebyshev_polynomial(order[1],[x[1]])
  polynomial_2 = chebyshev_polynomial(order[2],[x[2]])
  polynomial_3 = chebyshev_polynomial(order[3],[x[3]])
  polynomial_4 = chebyshev_polynomial(order[4],[x[4]])

  evaluated_polynomial = zero(T)

  for l = 1:order[4]+1
    for k = 1:order[3]+1
      for j = 1:order[2]+1
        for i = 1:order[1]+1

          evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*polynomial_3[k]*polynomial_4[l]*weights[i,j,k,l]

  			 end
			end
    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,5},x::Array{T,1},order::Array{S,1})

  polynomial_1 = chebyshev_polynomial(order[1],[x[1]])
  polynomial_2 = chebyshev_polynomial(order[2],[x[2]])
  polynomial_3 = chebyshev_polynomial(order[3],[x[3]])
  polynomial_4 = chebyshev_polynomial(order[4],[x[4]])
  polynomial_5 = chebyshev_polynomial(order[5],[x[5]])

  evaluated_polynomial = zero(T)

	for m = 1:order[5]+1
	  for l = 1:order[4]+1
      for k = 1:order[3]+1
        for j = 1:order[2]+1
          for i = 1:order[1]+1

            evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*polynomial_3[k]*polynomial_4[l]*polynomial_5[m]*weights[i,j,k,l,m]

					end
				end
      end
    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,1},x::Array{T,1},order::S)

  polynomial_1 = chebyshev_polynomial(order,[x[1]])

  evaluated_polynomial = zero(T)

  for i = 1:order+1

    evaluated_polynomial += polynomial_1[i]*weights[i]

  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,2},x::Array{T,1},order::S)

  polynomial_1 = chebyshev_polynomial(order,[x[1]])
  polynomial_2 = chebyshev_polynomial(order,[x[2]])

  evaluated_polynomial = zero(T)

  for j = 1:order+1
    for i = 1:order+1

      if (i+j <= order+2)

        evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*weights[i,j]

      end

    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,3},x::Array{T,1},order::S)

  polynomial_1 = chebyshev_polynomial(order,[x[1]])
  polynomial_2 = chebyshev_polynomial(order,[x[2]])
  polynomial_3 = chebyshev_polynomial(order,[x[3]])

  evaluated_polynomial = zero(T)

  for k = 1:order+1
    for j = 1:order+1
      for i = 1:order+1

        if (i+j+k <= order+3)

          evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*polynomial_3[k]*weights[i,j,k]

        end

      end
    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,4},x::Array{T,1},order::S)

  polynomial_1 = chebyshev_polynomial(order,[x[1]])
  polynomial_2 = chebyshev_polynomial(order,[x[2]])
  polynomial_3 = chebyshev_polynomial(order,[x[3]])
  polynomial_4 = chebyshev_polynomial(order,[x[4]])

  evaluated_polynomial = zero(T)

  for l = 1:order+1
    for k = 1:order+1
      for j = 1:order+1
        for i = 1:order+1

          if (i+j+k+l <= order+4)

            evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*polynomial_3[k]*polynomial_4[l]*weights[i,j,k,l]

          end

        end
			end
    end
  end

  return evaluated_polynomial

end

function chebyshev_evaluate{T<:FloatingPoint,S<:Integer}(weights::Array{T,5},x::Array{T,1},order::S)

  polynomial_1 = chebyshev_polynomial(order,[x[1]])
  polynomial_2 = chebyshev_polynomial(order,[x[2]])
  polynomial_3 = chebyshev_polynomial(order,[x[3]])
  polynomial_4 = chebyshev_polynomial(order,[x[4]])
  polynomial_5 = chebyshev_polynomial(order,[x[5]])

  evaluated_polynomial = zero(T)

  for m = 1:order+1
    for l = 1:order+1
      for k = 1:order+1
	      for j = 1:order+1
          for i = 1:order+1

            if (i+j+k+l+m <= order+5)

              evaluated_polynomial += polynomial_1[i]*polynomial_2[j]*polynomial_3[k]*polynomial_4[l]*polynomial_5[m]*weights[i,j,k,l,m]

            end

					end
				end
      end
    end
  end

  return evaluated_polynomial

end
