function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,1},nodes_1::Array{T,1},order::Array{S,1},range::Array{T,1})

  x1 = normalize_node(nodes_1,range)

  polynomial_1 = chebyshev_polynomial(order[1],x1)

  weights = Array(T,order[1]+1)

  for i = 1:order[1]+1

    numerator   = zero(T)
    denominator = zero(T)

    for s1 = 1:length(nodes_1)

      numerator   += f[s1]*polynomial_1[s1,i]
      denominator += (polynomial_1[s1,i])^2

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,2},nodes_1::Array{T,1},nodes_2::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])

  polynomial_1 = chebyshev_polynomial(order[1],x1)
  polynomial_2 = chebyshev_polynomial(order[2],x2)

  weights = Array(T,order[1]+1,order[2]+1)

  for j = 1:order[2]+1
    for i = 1:order[1]+1

      numerator   = zero(T)
      denominator = zero(T)

      for s2 = 1:length(nodes_2)
        for s1 = 1:length(nodes_1)

          numerator   += f[s1,s2]*polynomial_1[s1,i]*polynomial_2[s2,j]
          denominator += (polynomial_1[s1,i]*polynomial_2[s2,j])^2

        end
      end

      weights[i,j] = numerator/denominator

    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,3},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])
  x3 = normalize_node(nodes_3,range[:,3])

  polynomial_1 = chebyshev_polynomial(order[1],x1)
  polynomial_2 = chebyshev_polynomial(order[2],x2)
  polynomial_3 = chebyshev_polynomial(order[3],x3)

  weights = Array(T,order[1]+1,order[2]+1,order[3]+1)

  for k = 1:order[3]+1
    for j = 1:order[2]+1
      for i = 1:order[1]+1

        numerator   = zero(T)
        denominator = zero(T)

        for s3 = 1:length(nodes_3)
          for s2 = 1:length(nodes_2)
            for s1 = 1:length(nodes_1)

              numerator   += f[s1,s2,s3]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]
              denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k])^(2)

            end
          end
        end

        weights[i,j,k] = numerator/denominator

      end
    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,4},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])
  x3 = normalize_node(nodes_3,range[:,3])
  x4 = normalize_node(nodes_4,range[:,4])

  polynomial_1 = chebyshev_polynomial(order[1],x1)
  polynomial_2 = chebyshev_polynomial(order[2],x2)
  polynomial_3 = chebyshev_polynomial(order[3],x3)
  polynomial_4 = chebyshev_polynomial(order[4],x4)

  weights = Array(T,order[1]+1,order[2]+1,order[3]+1,order[4]+1)

  for l = 1:order[4]+1
    for k = 1:order[3]+1
      for j = 1:order[2]+1
        for i = 1:order[1]+1

          numerator   = zero(T)
          denominator = zero(T)

          for s4 = 1:length(nodes_4)
            for s3 = 1:length(nodes_3)
              for s2 = 1:length(nodes_2)
                for s1 = 1:length(nodes_1)

                  numerator   += f[s1,s2,s3,s4]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]
                  denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l])^(2.0)

                end
              end
            end
				  end

          weights[i,j,k,l] = numerator/denominator

        end
      end
    end
	 end

   return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,5},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])
  x3 = normalize_node(nodes_3,range[:,3])
  x4 = normalize_node(nodes_4,range[:,4])
  x5 = normalize_node(nodes_5,range[:,5])

  polynomial_1 = chebyshev_polynomial(order[1],x1)
  polynomial_2 = chebyshev_polynomial(order[2],x2)
  polynomial_3 = chebyshev_polynomial(order[3],x3)
  polynomial_4 = chebyshev_polynomial(order[4],x4)
  polynomial_5 = chebyshev_polynomial(order[5],x5)

  weights = Array(T,order[1]+1,order[2]+1,order[3]+1,order[4]+1,order[5]+1)

	for m = 1:order[5]+1
	  for l = 1:order[4]+1
      for k = 1:order[3]+1
        for j = 1:order[2]+1
          for i = 1:order[1]+1

            numerator   = zero(T)
            denominator = zero(T)

            for s5 = 1:length(nodes_5)
              for s4 = 1:length(nodes_4)
                for s3 = 1:length(nodes_3)
                  for s2 = 1:length(nodes_2)
                    for s1 = 1:length(nodes_1)

                      numerator   += f[s1,s2,s3,s4,s5]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m]
                      denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m])^(2.0)

                    end
                  end
                end
              end
						end

            weights[i,j,k,l,m] = numerator/denominator

          end
        end
      end
    end
	end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,6},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},nodes_6::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])
  x3 = normalize_node(nodes_3,range[:,3])
  x4 = normalize_node(nodes_4,range[:,4])
  x5 = normalize_node(nodes_5,range[:,5])
  x6 = normalize_node(nodes_6,range[:,6])

  polynomial_1 = chebyshev_polynomial(order[1],x1)
  polynomial_2 = chebyshev_polynomial(order[2],x2)
  polynomial_3 = chebyshev_polynomial(order[3],x3)
  polynomial_4 = chebyshev_polynomial(order[4],x4)
  polynomial_5 = chebyshev_polynomial(order[5],x5)
  polynomial_6 = chebyshev_polynomial(order[6],x6)

  weights = Array(T,order[1]+1,order[2]+1,order[3]+1,order[4]+1,order[5]+1,order[6]+1)

	for n = 1:order[6]+1
	  for m = 1:order[5]+1
	    for l = 1:order[4]+1
        for k = 1:order[3]+1
          for j = 1:order[2]+1
            for i = 1:order[1]+1

              numerator   = zero(T)
              denominator = zero(T)

              for s6 = 1:length(nodes_6)
                for s5 = 1:length(nodes_5)
                  for s4 = 1:length(nodes_4)
                    for s3 = 1:length(nodes_3)
                      for s2 = 1:length(nodes_2)
                        for s1 = 1:length(nodes_1)

                          numerator   += f[s1,s2,s3,s4,s5,s6]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m]*polynomial_6[s6,n]
                          denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m]*polynomial_6[s6,n])^(2.0)

                        end
                      end
                    end
                  end
                end
						  end

              weights[i,j,k,l,m,n] = numerator/denominator

            end
          end
        end
      end
    end
	end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,1},nodes_1::Array{T,1},order::S,range::Array{T,1})

  x1 = normalize_node(nodes_1,range)

  polynomial_1 = chebyshev_polynomial(order,x1)

  weights = Array(T,order+1)

  for i = 1:order+1

    numerator   = zero(T)
    denominator = zero(T)

    for s1 = 1:length(nodes_1)

      numerator   += f[s1]*polynomial_1[s1,i]
      denominator += (polynomial_1[s1,i])^2

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,2},nodes_1::Array{T,1},nodes_2::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])

  polynomial_1 = chebyshev_polynomial(order,x1)
  polynomial_2 = chebyshev_polynomial(order,x2)

  weights = zeros(order+1,order+1)

  for j = 1:order+1
    for i = 1:order+1

      numerator   = zero(T)
      denominator = zero(T)

      if (i+j <= order+2)

        for s2 = 1:length(nodes_2)
          for s1 = 1:length(nodes_1)

            numerator   += f[s1,s2]*polynomial_1[s1,i]*polynomial_2[s2,j]
            denominator += (polynomial_1[s1,i]*polynomial_2[s2,j])^(2)

          end
        end

        weights[i,j] = numerator/denominator

      end

    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,3},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])
  x3 = normalize_node(nodes_3,range[:,3])

  polynomial_1 = chebyshev_polynomial(order,x1)
  polynomial_2 = chebyshev_polynomial(order,x2)
  polynomial_3 = chebyshev_polynomial(order,x3)

  weights = zeros(order+1,order+1,order+1)

  for k = 1:order+1
    for j = 1:order+1
      for i = 1:order+1

        numerator   = zero(T)
        denominator = zero(T)

        if (i+j+k <= order+3)

          for s3 = 1:length(nodes_3)
            for s2 = 1:length(nodes_2)
              for s1 = 1:length(nodes_1)

                numerator   += f[s1,s2,s3]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]
                denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k])^2

              end
            end
          end

          weights[i,j,k] = numerator/denominator

        end

      end
    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,4},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])
  x3 = normalize_node(nodes_3,range[:,3])
  x4 = normalize_node(nodes_4,range[:,4])

  polynomial_1 = chebyshev_polynomial(order,x1)
  polynomial_2 = chebyshev_polynomial(order,x2)
  polynomial_3 = chebyshev_polynomial(order,x3)
  polynomial_4 = chebyshev_polynomial(order,x4)

  weights = zeros(order+1,order+1,order+1,order+1)

  for l = 1:order+1
    for k = 1:order+1
      for j = 1:order+1
        for i = 1:order+1

          numerator   = zero(T)
          denominator = zero(T)

          if (i+j+k+l <= order+4)

            for s4 = 1:length(nodes_4)
              for s3 = 1:length(nodes_3)
                for s2 = 1:length(nodes_2)
                  for s1 = 1:length(nodes_1)

                    numerator   += f[s1,s2,s3,s4]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]
                    denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l])^2

                  end
                end
              end
            end

            weights[i,j,k,l] = numerator/denominator

          end

        end
      end
    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,5},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])
  x3 = normalize_node(nodes_3,range[:,3])
  x4 = normalize_node(nodes_4,range[:,4])
  x5 = normalize_node(nodes_5,range[:,5])

  polynomial_1 = chebyshev_polynomial(order,x1)
  polynomial_2 = chebyshev_polynomial(order,x2)
  polynomial_3 = chebyshev_polynomial(order,x3)
  polynomial_4 = chebyshev_polynomial(order,x4)
  polynomial_5 = chebyshev_polynomial(order,x5)

  weights = zeros(order+1,order+1,order+1,order+1,order+1)

  for m = 1:order+1
    for l = 1:order+1
      for k = 1:order+1
        for j = 1:order+1
          for i = 1:order+1

            numerator   = zero(T)
            denominator = zero(T)

            if (i+j+k+l+m <= order+5)

              for s5 = 1:length(nodes_5)
                for s4 = 1:length(nodes_4)
                  for s3 = 1:length(nodes_3)
                    for s2 = 1:length(nodes_2)
                      for s1 = 1:length(nodes_1)

                        numerator   += f[s1,s2,s3,s4,s5]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m]
                        denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m])^2

                      end
                    end
                  end
                end
              end

              weights[i,j,k,l,m] = numerator/denominator

            end

          end
        end
      end
    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,6},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},nodes_6::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])
  x3 = normalize_node(nodes_3,range[:,3])
  x4 = normalize_node(nodes_4,range[:,4])
  x5 = normalize_node(nodes_5,range[:,5])
  x6 = normalize_node(nodes_6,range[:,6])

  polynomial_1 = chebyshev_polynomial(order,x1)
  polynomial_2 = chebyshev_polynomial(order,x2)
  polynomial_3 = chebyshev_polynomial(order,x3)
  polynomial_4 = chebyshev_polynomial(order,x4)
  polynomial_5 = chebyshev_polynomial(order,x5)
  polynomial_6 = chebyshev_polynomial(order,x6)

  weights = zeros(order+1,order+1,order+1,order+1,order+1,order+1)

  for n = 1:order+1
    for m = 1:order+1
      for l = 1:order+1
        for k = 1:order+1
          for j = 1:order+1
            for i = 1:order+1

              numerator   = zero(T)
              denominator = zero(T)

              if (i+j+k+l+m+n <= order+6)

                for s6 = 1:length(nodes_6)
                  for s5 = 1:length(nodes_5)
                    for s4 = 1:length(nodes_4)
                      for s3 = 1:length(nodes_3)
                        for s2 = 1:length(nodes_2)
                          for s1 = 1:length(nodes_1)

                            numerator   += f[s1,s2,s3,s4,s5,s6]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m]*polynomial_6[s6,n]
                            denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m]*polynomial_6[s6,n])^2.0

                          end
                        end
                      end
                    end
                  end
                end

                weights[i,j,k,l,m,n] = numerator/denominator

              end

            end
          end
        end
      end
    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,1},nodes_1::Array{T,1},order::Array{S,1})

  polynomial_1 = chebyshev_polynomial(order[1],nodes_1)

  weights = Array(T,order[1]+1)

  for i = 1:order[1]+1

    numerator   = zero(T)
    denominator = zero(T)

    for s1 = 1:length(nodes_1)

      numerator   += f[s1]*polynomial_1[s1,i]
      denominator += (polynomial_1[s1,i])^2

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,2},nodes_1::Array{T,1},nodes_2::Array{T,1},order::Array{S,1})

  polynomial_1 = chebyshev_polynomial(order[1],nodes_1)
  polynomial_2 = chebyshev_polynomial(order[2],nodes_2)

  weights = Array(T,order[1]+1,order[2]+1)

  for j = 1:order[2]+1
    for i = 1:order[1]+1

      numerator   = zero(T)
      denominator = zero(T)

      for s2 = 1:length(nodes_2)
        for s1 = 1:length(nodes_1)

          numerator   += f[s1,s2]*polynomial_1[s1,i]*polynomial_2[s2,j]
          denominator += (polynomial_1[s1,i]*polynomial_2[s2,j])^2

        end
      end

      weights[i,j] = numerator/denominator

    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,3},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},order::Array{S,1})

  polynomial_1 = chebyshev_polynomial(order[1],nodes_1)
  polynomial_2 = chebyshev_polynomial(order[2],nodes_2)
  polynomial_3 = chebyshev_polynomial(order[3],nodes_3)

  weights = Array(T,order[1]+1,order[2]+1,order[3]+1)

  for k = 1:order[3]+1
    for j = 1:order[2]+1
      for i = 1:order[1]+1

        numerator   = zero(T)
        denominator = zero(T)

        for s3 = 1:length(nodes_3)
          for s2 = 1:length(nodes_2)
            for s1 = 1:length(nodes_1)

              numerator   += f[s1,s2,s3]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]
              denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k])^(2)

            end
          end
        end

        weights[i,j,k] = numerator/denominator

      end
    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,4},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},order::Array{S,1})

  polynomial_1 = chebyshev_polynomial(order[1],nodes_1)
  polynomial_2 = chebyshev_polynomial(order[2],nodes_2)
  polynomial_3 = chebyshev_polynomial(order[3],nodes_3)
  polynomial_4 = chebyshev_polynomial(order[4],nodes_4)

  weights = Array(T,order[1]+1,order[2]+1,order[3]+1,order[4]+1)

  for l = 1:order[4]+1
    for k = 1:order[3]+1
      for j = 1:order[2]+1
        for i = 1:order[1]+1

          numerator   = zero(T)
          denominator = zero(T)

          for s4 = 1:length(nodes_4)
            for s3 = 1:length(nodes_3)
              for s2 = 1:length(nodes_2)
                for s1 = 1:length(nodes_1)

                  numerator   += f[s1,s2,s3,s4]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]
                  denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l])^(2.0)

                end
              end
            end
				  end

          weights[i,j,k,l] = numerator/denominator

        end
      end
    end
	 end

   return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,5},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},order::Array{S,1})

  polynomial_1 = chebyshev_polynomial(order[1],nodes_1)
  polynomial_2 = chebyshev_polynomial(order[2],nodes_2)
  polynomial_3 = chebyshev_polynomial(order[3],nodes_3)
  polynomial_4 = chebyshev_polynomial(order[4],nodes_4)
  polynomial_5 = chebyshev_polynomial(order[5],nodes_5)

  weights = Array(T,order[1]+1,order[2]+1,order[3]+1,order[4]+1,order[5]+1)

	for m = 1:order[5]+1
	  for l = 1:order[4]+1
      for k = 1:order[3]+1
        for j = 1:order[2]+1
          for i = 1:order[1]+1

            numerator   = zero(T)
            denominator = zero(T)

            for s5 = 1:length(nodes_5)
              for s4 = 1:length(nodes_4)
                for s3 = 1:length(nodes_3)
                  for s2 = 1:length(nodes_2)
                    for s1 = 1:length(nodes_1)

                      numerator   += f[s1,s2,s3,s4,s5]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m]
                      denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m])^2

                    end
                  end
                end
              end
						end

            weights[i,j,k,l,m] = numerator/denominator

          end
        end
      end
    end
	end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,6},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},nodes_6::Array{T,1},order::Array{S,1})

  polynomial_1 = chebyshev_polynomial(order[1],nodes_1)
  polynomial_2 = chebyshev_polynomial(order[2],nodes_2)
  polynomial_3 = chebyshev_polynomial(order[3],nodes_3)
  polynomial_4 = chebyshev_polynomial(order[4],nodes_4)
  polynomial_5 = chebyshev_polynomial(order[5],nodes_5)
  polynomial_6 = chebyshev_polynomial(order[6],nodes_6)

  weights = Array(T,order[1]+1,order[2]+1,order[3]+1,order[4]+1,order[5]+1,order[6]+1)

	for n = 1:order[6]+1
	  for m = 1:order[5]+1
	    for l = 1:order[4]+1
        for k = 1:order[3]+1
          for j = 1:order[2]+1
            for i = 1:order[1]+1

              numerator   = zero(T)
              denominator = zero(T)

              for s6 = 1:length(nodes_6)
                for s5 = 1:length(nodes_5)
                  for s4 = 1:length(nodes_4)
                    for s3 = 1:length(nodes_3)
                      for s2 = 1:length(nodes_2)
                        for s1 = 1:length(nodes_1)

                          numerator   += f[s1,s2,s3,s4,s5,s6]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m]*polynomial_6[s6,n]
                          denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m]*polynomial_6[s6,n])^2.0

                        end
                      end
                    end
                  end
                end
						  end

              weights[i,j,k,l,m,n] = numerator/denominator

            end
          end
        end
      end
    end
	end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,1},nodes_1::Array{T,1},order::S)

  polynomial_1 = chebyshev_polynomial(order,nodes_1)

  weights = Array(T,order+1)

  for i = 1:order+1

    numerator   = zero(T)
    denominator = zero(T)

    for s1 = 1:length(nodes_1)

      numerator   += f[s1]*polynomial_1[s1,i]
      denominator += (polynomial_1[s1,i])^2

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,2},nodes_1::Array{T,1},nodes_2::Array{T,1},order::S)

  polynomial_1 = chebyshev_polynomial(order,nodes_1)
  polynomial_2 = chebyshev_polynomial(order,nodes_2)

  weights = zeros(order+1,order+1)

  for j = 1:order+1
    for i = 1:order+1

      numerator   = zero(T)
      denominator = zero(T)

      if (i+j <= order+2)

        for s2 = 1:length(nodes_2)
          for s1 = 1:length(nodes_1)

            numerator   += f[s1,s2]*polynomial_1[s1,i]*polynomial_2[s2,j]
            denominator += (polynomial_1[s1,i]*polynomial_2[s2,j])^(2)

          end
        end

        weights[i,j] = numerator/denominator

      end

    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,3},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},order::S)

  polynomial_1 = chebyshev_polynomial(order,nodes_1)
  polynomial_2 = chebyshev_polynomial(order,nodes_2)
  polynomial_3 = chebyshev_polynomial(order,nodes_3)

  weights = zeros(order+1,order+1,order+1)

  for k = 1:order+1
    for j = 1:order+1
      for i = 1:order+1

        numerator   = zero(T)
        denominator = zero(T)

        if (i+j+k <= order+3)

          for s3 = 1:length(nodes_3)
            for s2 = 1:length(nodes_2)
              for s1 = 1:length(nodes_1)

              numerator   += f[s1,s2,s3]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]
              denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k])^2

              end
            end
          end

          weights[i,j,k] = numerator/denominator

        end

      end
    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,4},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},order::S)

  polynomial_1 = chebyshev_polynomial(order,nodes_1)
  polynomial_2 = chebyshev_polynomial(order,nodes_2)
  polynomial_3 = chebyshev_polynomial(order,nodes_3)
  polynomial_4 = chebyshev_polynomial(order,nodes_4)

  weights = zeros(order+1,order+1,order+1,order+1)

  for l = 1:order+1
    for k = 1:order+1
      for j = 1:order+1
        for i = 1:order+1

          numerator   = zero(T)
          denominator = zero(T)

          if (i+j+k+l <= order+4)

            for s4 = 1:length(nodes_4)
              for s3 = 1:length(nodes_3)
                for s2 = 1:length(nodes_2)
                  for s1 = 1:length(nodes_1)

                    numerator   += f[s1,s2,s3,s4]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]
                    denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l])^2

                  end
                end
              end
            end

            weights[i,j,k,l] = numerator/denominator

          end

        end
      end
    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,5},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},order::S)

  polynomial_1 = chebyshev_polynomial(order,nodes_1)
  polynomial_2 = chebyshev_polynomial(order,nodes_2)
  polynomial_3 = chebyshev_polynomial(order,nodes_3)
  polynomial_4 = chebyshev_polynomial(order,nodes_4)
  polynomial_5 = chebyshev_polynomial(order,nodes_5)

  weights = zeros(order+1,order+1,order+1,order+1,order+1)

  for m = 1:order+1
    for l = 1:order+1
      for k = 1:order+1
        for j = 1:order+1
          for i = 1:order+1

            numerator   = zero(T)
            denominator = zero(T)

            if (i+j+k+l+m <= order+5)

              for s5 = 1:length(nodes_5)
                for s4 = 1:length(nodes_4)
                  for s3 = 1:length(nodes_3)
                    for s2 = 1:length(nodes_2)
                      for s1 = 1:length(nodes_1)

                        numerator   += f[s1,s2,s3,s4,s5]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m]
                        denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m])^2

                      end
                    end
                  end
                end
              end

              weights[i,j,k,l,m] = numerator/denominator

            end

          end
        end
      end
    end
  end

  return weights

end

function chebyshev_weights{T<:FloatingPoint,S<:Integer}(f::Array{T,6},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},nodes_6::Array{T,1},order::S)

  polynomial_1 = chebyshev_polynomial(order,nodes_1)
  polynomial_2 = chebyshev_polynomial(order,nodes_2)
  polynomial_3 = chebyshev_polynomial(order,nodes_3)
  polynomial_4 = chebyshev_polynomial(order,nodes_4)
  polynomial_5 = chebyshev_polynomial(order,nodes_5)
  polynomial_6 = chebyshev_polynomial(order,nodes_6)

  weights = zeros(order+1,order+1,order+1,order+1,order+1,order+1)

  for n = 1:order+1
    for m = 1:order+1
      for l = 1:order+1
        for k = 1:order+1
          for j = 1:order+1
            for i = 1:order+1

              numerator   = zero(T)
              denominator = zero(T)

              if (i+j+k+l+m+n <= order+6)

                for s6 = 1:length(nodes_6)
                  for s5 = 1:length(nodes_5)
                    for s4 = 1:length(nodes_4)
                      for s3 = 1:length(nodes_3)
                        for s2 = 1:length(nodes_2)
                          for s1 = 1:length(nodes_1)

                            numerator   += f[s1,s2,s3,s4,s5,s6]*polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m]*polynomial_6[s6,n]
                            denominator += (polynomial_1[s1,i]*polynomial_2[s2,j]*polynomial_3[s3,k]*polynomial_4[s4,l]*polynomial_5[s5,m]*polynomial_6[s6,n])^2.0

                          end
                        end
                      end
                    end
                  end
                end

                weights[i,j,k,l,m,n] = numerator/denominator

              end

            end
          end
        end
      end
    end
  end

  return weights

end
