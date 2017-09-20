# Using Clenshaw's method to evaulate tensor-product Chebyshev polynomials

@generated function clenshaw_recursion(weights::Array{T,N},x::T) where {T,N}

  first_term = :( p = Array{T}(Base.tail(size(weights))) )

  i_vars = Array{Symbol}(N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
  end

  inner_term = :( z = zeros(T,size(weights,1)+2);
                  for $(i_vars[1]) in size(weights,1):-1:1;
                    z[$(i_vars[1])] = weights[$(i_vars...)]+2*x*z[$(i_vars[1])+1]-z[$(i_vars[1])+2];
                  end;
                  p[$(i_vars[2:N]...)] = z[1]-x*z[2]
                  )

  outer_term = inner_term
  for i = N:-1:2
    outer_term = :( for $(i_vars[i]) = 1:size(weights,$i);
                      $outer_term;
                    end
                    )
  end

  final = :( $first_term;
             $outer_term;
             return p
             )

  return final

end

function clenshaw_evaluate(weights::Array{T,N},x::Array{T,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N}

  for i = 1:N
    normalized_x = normalize_node(x[i],domain[:,i])
    weights = clenshaw_recursion(weights,normalized_x)
  end

  return weights[1]

end
