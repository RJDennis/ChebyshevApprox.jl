# Functions that compute weights for complete polynomials with up to six dimensions

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::Array{T,1},nodes_1::Array{T,1},order::S,range::Array{T,1})

  x1 = normalize_node(nodes_1,range)

  polynomial_1 = chebyshev_polynomial(order,[x1;])

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::Array{T,2},nodes_1::Array{T,1},nodes_2::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])

  polynomial_1 = chebyshev_polynomial(order,[x1;])
  polynomial_2 = chebyshev_polynomial(order,[x2;])

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::Array{T,3},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])
  x3 = normalize_node(nodes_3,range[:,3])

  polynomial_1 = chebyshev_polynomial(order,[x1;])
  polynomial_2 = chebyshev_polynomial(order,[x2;])
  polynomial_3 = chebyshev_polynomial(order,[x3;])

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::Array{T,4},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])
  x3 = normalize_node(nodes_3,range[:,3])
  x4 = normalize_node(nodes_4,range[:,4])

  polynomial_1 = chebyshev_polynomial(order,[x1;])
  polynomial_2 = chebyshev_polynomial(order,[x2;])
  polynomial_3 = chebyshev_polynomial(order,[x3;])
  polynomial_4 = chebyshev_polynomial(order,[x4;])

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::Array{T,5},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])
  x3 = normalize_node(nodes_3,range[:,3])
  x4 = normalize_node(nodes_4,range[:,4])
  x5 = normalize_node(nodes_5,range[:,5])

  polynomial_1 = chebyshev_polynomial(order,[x1;])
  polynomial_2 = chebyshev_polynomial(order,[x2;])
  polynomial_3 = chebyshev_polynomial(order,[x3;])
  polynomial_4 = chebyshev_polynomial(order,[x4;])
  polynomial_5 = chebyshev_polynomial(order,[x5;])

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::Array{T,6},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},nodes_6::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(nodes_1,range[:,1])
  x2 = normalize_node(nodes_2,range[:,2])
  x3 = normalize_node(nodes_3,range[:,3])
  x4 = normalize_node(nodes_4,range[:,4])
  x5 = normalize_node(nodes_5,range[:,5])
  x6 = normalize_node(nodes_6,range[:,6])

  polynomial_1 = chebyshev_polynomial(order,[x1;])
  polynomial_2 = chebyshev_polynomial(order,[x2;])
  polynomial_3 = chebyshev_polynomial(order,[x3;])
  polynomial_4 = chebyshev_polynomial(order,[x4;])
  polynomial_5 = chebyshev_polynomial(order,[x5;])
  polynomial_6 = chebyshev_polynomial(order,[x6;])

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::Array{T,1},nodes_1::Array{T,1},order::S)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::Array{T,2},nodes_1::Array{T,1},nodes_2::Array{T,1},order::S)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::Array{T,3},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},order::S)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::Array{T,4},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},order::S)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::Array{T,5},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},order::S)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::Array{T,6},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},nodes_6::Array{T,1},order::S)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::AbstractArray{T,1},polynomial_1::Array{T,2},order::S)

  weights = Array(T,order+1)

  for i = 1:order+1

    numerator   = zero(T)
    denominator = zero(T)

    for s1 = 1:size(polynomial_1,1)

      numerator   += f[s1]*polynomial_1[s1,i]
      denominator += (polynomial_1[s1,i])^2

    end

    weights[i] = numerator/denominator

  end

  return weights

end

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::AbstractArray{T,2},polynomial_1::Array{T,2},polynomial_2::Array{T,2},order::S)

  weights = zeros(order+1,order+1)

  for j = 1:order+1
    for i = 1:order+1

      numerator   = zero(T)
      denominator = zero(T)

      if (i+j <= order+2)

        for s2 = 1:size(polynomial_2,1)
          for s1 = 1:size(polynomial_1,1)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::AbstractArray{T,3},polynomial_1::Array{T,2},polynomial_2::Array{T,2},polynomial_3::Array{T,2},order::S)

  weights = zeros(order+1,order+1,order+1)

  for k = 1:order+1
    for j = 1:order+1
      for i = 1:order+1

        numerator   = zero(T)
        denominator = zero(T)

        if (i+j+k <= order+3)

          for s3 = 1:size(polynomial_3,1)
            for s2 = 1:size(polynomial_2,1)
              for s1 = 1:size(polynomial_1,1)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::AbstractArray{T,4},polynomial_1::Array{T,2},polynomial_2::Array{T,2},polynomial_3::Array{T,2},polynomial_4::Array{T,2},order::S)

  weights = zeros(order+1,order+1,order+1,order+1)

  for l = 1:order+1
    for k = 1:order+1
      for j = 1:order+1
        for i = 1:order+1

          numerator   = zero(T)
          denominator = zero(T)

          if (i+j+k+l <= order+4)

            for s4 = 1:size(polynomial_4,1)
              for s3 = 1:size(polynomial_3,1)
                for s2 = 1:size(polynomial_2,1)
                  for s1 = 1:size(polynomial_1,1)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::AbstractArray{T,5},polynomial_1::Array{T,2},polynomial_2::Array{T,2},polynomial_3::Array{T,2},polynomial_4::Array{T,2},polynomial_5::Array{T,2},order::S)

  weights = zeros(order+1,order+1,order+1,order+1,order+1)

  for m = 1:order+1
    for l = 1:order+1
      for k = 1:order+1
        for j = 1:order+1
          for i = 1:order+1

            numerator   = zero(T)
            denominator = zero(T)

            if (i+j+k+l+m <= order+5)

              for s5 = 1:size(polynomial_5,1)
                for s4 = 1:size(polynomial_4,1)
                  for s3 = 1:size(polynomial_3,1)
                    for s2 = 1:size(polynomial_2,1)
                      for s1 = 1:size(polynomial_1,1)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::AbstractArray{T,6},polynomial_1::Array{T,2},polynomial_2::Array{T,2},polynomial_3::Array{T,2},polynomial_4::Array{T,2},polynomial_5::Array{T,2},polynomial_6::Array{T,2},order::S)

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

                for s6 = 1:size(polynomial_6,1)
                  for s5 = 1:size(polynomial_5,1)
                    for s4 = 1:size(polynomial_4,1)
                      for s3 = 1:size(polynomial_3,1)
                        for s2 = 1:size(polynomial_2,1)
                          for s1 = 1:size(polynomial_1,1)

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

# Generated function to compute Chebyshev weights using Chebyshev regression

@generated function chebyshev_weights{T,N,S}(f::Array{T,N},nodes::Array{Array{T,1},1},order::Array{S,1},range::Array{T,2})

  chebyshev_polynomials = :( poly = Array{T,2}[];
                             for k = 1:N;
                               orderk = order[k];
                               xk = nodes[k];

                               # normalize nodes

                               if range[1,k] == range[2,k];
                                 fill!(xk,(range[1,k]+range[2,k])/2);
                               else;
                                 xk = 2*(xk-range[2,k])/(range[1,k]-range[2,k])-one(T);
                               end;

                               polynomial = Array(T,length(xk),orderk+1);
                               for i = 1:orderk+1;
                                 for j = 1:length(xk);
                                   if i == 1;
                                     polynomial[j,i] = one(T);
                                   elseif i == 2;
                                     polynomial[j,i] = xk[j];
                                   elseif i == 3;
                                     polynomial[j,i] = 2*xk[j]*xk[j]-one(T);
                                   else;
                                     polynomial[j,i] = 2*xk[j]*polynomial[j,i-1]-polynomial[j,i-2];
                                   end;
                                 end;
                               end;
                               push!(poly,polynomial);
                             end )

  i_vars = Array(Symbol,N)
  s_vars = Array(Symbol,N)
  for i = 1:N
    i_vars[i] = symbol("i$i")
    s_vars[i] = symbol("s$i")
  end

  # Construct denominator term

  inner_prod_denominator = :( poly[1][s1,i1] )
  for i = 2:N
    inner_prod_denominator = :( poly[$i][$(s_vars[i]),$(i_vars[i])]*$inner_prod_denominator )
  end
  denominator_term = :( $inner_prod_denominator*$inner_prod_denominator )

  # Construct numerator term

  inner_prod_numerator = :( poly[1][s1,i1]*f[$(s_vars...)] )
  for i = 2:N
      inner_prod_numerator = :( poly[$i][$(s_vars[i]),$(i_vars[i])]*$inner_prod_numerator )
  end
  numerator_term = inner_prod_numerator

  # Construct weights term

  weights_term = :( weights[$(i_vars...)] = numerator/denominator )

  # Construct the inner loops that compute numerator and denominator

  inner = :( numerator += $numerator_term;
             denominator += $denominator_term )
  outer  = inner

  for i = 1:N
    outer = :(
      for $(s_vars[i]) = 1:size(f,$i)
        $outer
      end
    )
  end

  # Construct the outer loops that compute the weights

  new_inner = :( numerator = zero(T);
                 denominator = zero(T);
                 $outer;
                 $weights_term )
  new_outer = new_inner

  for i = 1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order[$i]+1)
        $new_outer
      end
    )
  end

  # Initialize the weight Array

  initial_weight = string("weights = Array(T,order[1]+1,",)
  for i = 2:N
    initial_weight = string(initial_weight,"order[",i,"]+1,")
  end
  initial_weight = parse(string(initial_weight,")"))

    # Put it all together to compute the weights array

  final = :( $chebyshev_polynomials;
             $initial_weight;
             $new_outer;
             return weights )

  return final

end

@generated function chebyshev_weights{T,N,S}(f::Array{T,N},poly::Array{Array{T,2},1},order::Array{S,1})

  i_vars = Array(Symbol,N)
  s_vars = Array(Symbol,N)
  for i = 1:N
    i_vars[i] = symbol("i$i")
    s_vars[i] = symbol("s$i")
  end

  # Construct denominator term

  inner_prod_denominator = :( poly[1][s1,i1] )
  for i = 2:N
    inner_prod_denominator = :( poly[$i][$(s_vars[i]),$(i_vars[i])]*$inner_prod_denominator )
  end
  denominator_term = :( $inner_prod_denominator*$inner_prod_denominator )

  # Construct numerator term

  inner_prod_numerator = :( poly[1][s1,i1]*f[$(s_vars...)] )
  for i = 2:N
      inner_prod_numerator = :( poly[$i][$(s_vars[i]),$(i_vars[i])]*$inner_prod_numerator )
  end
  numerator_term = inner_prod_numerator

  # Construct weights term

  weights_term = :( weights[$(i_vars...)] = numerator/denominator )

  # Construct the inner loops that compute numerator and denominator

  inner = :( numerator += $numerator_term;
             denominator += $denominator_term )
  outer  = inner

  for i = 1:N
    outer = :(
      for $(s_vars[i]) = 1:size(f,$i)
        $outer
      end
    )
  end

  # Construct the outer loops that compute the weights

  new_inner = :( numerator = zero(T);
                 denominator = zero(T);
                 $outer;
                 $weights_term )
  new_outer = new_inner

  for i = 1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order[$i]+1)
        $new_outer
      end
    )
  end

  # Initialize the weight Array

  initial_weight = string("weights = Array(T,order[1]+1,",)
  for i = 2:N
    initial_weight = string(initial_weight,"order[",i,"]+1,")
  end
  initial_weight = parse(string(initial_weight,")"))

    # Put it all together to compute the weights array

  final = :( $initial_weight;
             $new_outer;
             return weights )

  return final

end
