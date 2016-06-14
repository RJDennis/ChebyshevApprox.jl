function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::AbstractArray{T,1},polynomial_1::Array{T,2},order::Array{S,1})

  weights = Array(T,order[1]+1)

  for i = 1:order[1]+1

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::AbstractArray{T,2},polynomial_1::Array{T,2},polynomial_2::Array{T,2},order::Array{S,1})

  weights = Array(T,order[1]+1,order[2]+1)

  for j = 1:order[2]+1
    for i = 1:order[1]+1

      numerator   = zero(T)
      denominator = zero(T)

      for s2 = 1:size(polynomial_2,1)
        for s1 = 1:size(polynomial_1,1)

          numerator   += f[s1,s2]*polynomial_1[s1,i]*polynomial_2[s2,j]
          denominator += (polynomial_1[s1,i]*polynomial_2[s2,j])^2

        end
      end

      weights[i,j] = numerator/denominator

    end
  end

  return weights

end

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::AbstractArray{T,3},polynomial_1::Array{T,2},polynomial_2::Array{T,2},polynomial_3::Array{T,2},order::Array{S,1})

  weights = Array(T,order[1]+1,order[2]+1,order[3]+1)

  for k = 1:order[3]+1
    for j = 1:order[2]+1
      for i = 1:order[1]+1

        numerator   = zero(T)
        denominator = zero(T)

        for s3 = 1:size(polynomial_3,1)
          for s2 = 1:size(polynomial_2,1)
            for s1 = 1:size(polynomial_1,1)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::AbstractArray{T,4},polynomial_1::Array{T,2},polynomial_2::Array{T,2},polynomial_3::Array{T,2},polynomial_4::Array{T,1},order::Array{S,1})

  weights = Array(T,order[1]+1,order[2]+1,order[3]+1,order[4]+1)

  for l = 1:order[4]+1
    for k = 1:order[3]+1
      for j = 1:order[2]+1
        for i = 1:order[1]+1

          numerator   = zero(T)
          denominator = zero(T)

          for s4 = 1:size(polynomial_4,1)
            for s3 = 1:size(polynomial_3,1)
              for s2 = 1:size(polynomial_2,1)
                for s1 = 1:size(polynomial_1,1)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::AbstractArray{T,5},polynomial_1::Array{T,2},polynomial_2::Array{T,2},polynomial_3::Array{T,2},polynomial_4::Array{T,2},polynomial_5::Array{T,2},order::Array{S,1})

  weights = Array(T,order[1]+1,order[2]+1,order[3]+1,order[4]+1,order[5]+1)

	for m = 1:order[5]+1
	  for l = 1:order[4]+1
      for k = 1:order[3]+1
        for j = 1:order[2]+1
          for i = 1:order[1]+1

            numerator   = zero(T)
            denominator = zero(T)

            for s5 = 1:size(polynomial_5,1)
              for s4 = 1:size(polynomial_4,1)
                for s3 = 1:size(polynomial_3,1)
                  for s2 = 1:size(polynomial_2,1)
                    for s1 = 1:size(polynomial_1,1)

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

function chebyshev_weights{T<:AbstractFloat,S<:Integer}(f::AbstractArray{T,6},polynomial_1::Array{T,2},polynomial_2::Array{T,2},polynomial_3::Array{T,2},polynomial_4::Array{T,2},polynomial_5::Array{T,2},polynomial_6::Array{T,2},order::Array{S,1})

  weights = Array(T,order[1]+1,order[2]+1,order[3]+1,order[4]+1,order[5]+1,order[6]+1)

	for n = 1:order[6]+1
	  for m = 1:order[5]+1
	    for l = 1:order[4]+1
        for k = 1:order[3]+1
          for j = 1:order[2]+1
            for i = 1:order[1]+1

              numerator   = zero(T)
              denominator = zero(T)

              for s6 = 1:size(polynomial_6,1)
                for s5 = 1:size(polynomial_5,1)
                  for s4 = 1:size(polynomial_4,1)
                    for s3 = 1:size(polynomial_3,1)
                      for s2 = 1:size(polynomial_2,1)
                        for s1 = 1:size(polynomial_1,1)

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

@generated function chebyshev_weights{T,N,S}(f::Array{T,N},poly::Array{Array{T,2},1},order::Array{S,1})

  # Construct numerator term

  inner_prod_numerator = "f["
  for i = 1:N
    if i == N
      inner_prod_numerator = string("poly[",i,"][s",i,",i",i,"]*",inner_prod_numerator,"s",i)
    else
      inner_prod_numerator = string("poly[",i,"][s",i,",i",i,"]*",inner_prod_numerator,"s",i,",")
    end
  end
  inner_prod_numerator = string(inner_prod_numerator,"]")
  numerator_term = parse(inner_prod_numerator)

  # Construct denominator term

  inner_prod_denominator = string("poly[",1,"][s",1,",i",1,"]")
  for i = 2:N
    inner_prod_denominator = string("poly[",i,"][s",i,",i",i,"]*",inner_prod_denominator)
  end
  inner_prod_denominator = string("(",inner_prod_denominator,")^2.0")
  denominator_term = parse(inner_prod_denominator)

  # Construct weights term

  weights_term = string("weights[i",1)
  for i = 2:N
    weights_term = string(weights_term,",i",i)
  end
  weights_term = parse(string(weights_term,"] = numerator/denominator"))

  # Construct the inner loops that compute numerator and denominator

  inner = :( numerator += $numerator_term;
             denominator += $denominator_term )
  outer  = inner

  for i = 1:N
    var = symbol("s$i")
    outer = :(
      for $var = 1:size(f,$i)
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
    var = symbol("i$i")
    new_outer = :(
      for $var = 1:(order[$i]+1)
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
