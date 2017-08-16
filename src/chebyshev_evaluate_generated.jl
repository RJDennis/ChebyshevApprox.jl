# Generated functions for evaluating Chebyshev polynomials

@generated function chebyshev_evaluate(weights::Array{T,N},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T,N,S}

  chebyshev_polynomials = :( poly = Array{T,2}[];
                             for i = 1:size(x,1);
                               xi = x[i];

                               # Normalize nodes to [-1,1]

                               if domain[1,i] == domain[2,i];
                                 xi = (domain[1,i]+domain[2,i])/2;
                               else;
                                 xi = 2*(xi-domain[2,i])/(domain[1,i]-domain[2,i])-one(T);
                               end;

                               polynomial = Array{T}(1,order[i]+1);
                               for j = 1:order[i]+1;
                                 if j == 1;
                                   polynomial[j] = one(T);
                                 elseif j == 2;
                                   polynomial[j] = xi;
                                 else;
                                   polynomial[j] = 2*xi*polynomial[j-1]-polynomial[j-2];
                                 end;
                               end;
                               push!(poly,polynomial);
                             end
                             )

  i_vars = Array{Symbol}(N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
  end

  inner_prod = :( poly[1][i1]*weights[$(i_vars...)] )
  for i = 2:N
    inner_prod = :( poly[$i][$(i_vars[i])]*$inner_prod )
  end

  inner = :( evaluated_polynomial += $inner_prod )
  outer = inner
  for i = 1:N
    outer = :(
      for $(i_vars[i]) = 1:size(weights,$i)
        $outer
      end
    )
  end

  final = :( $chebyshev_polynomials;
             evaluated_polynomial = zero(T);
             $outer;
             return evaluated_polynomial
             )

  return final

end

@generated function chebyshev_evaluate(weights::Array{T,N},x::Array{T,1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T,N,S}

  chebyshev_polynomials = :( poly = Array{T,2}[];
                             for i = 1:size(x,1);
                               xi = x[i];

                               # Normalize nodes to [-1,1]

                               if domain[1,i] == domain[2,i];
                                 xi = (domain[1,i]+domain[2,i])/2;
                               else;
                                 xi = 2*(xi-domain[2,i])/(domain[1,i]-domain[2,i])-one(T);
                               end;

                               polynomial = Array{T}(1,order+1);
                               for j = 1:order+1;
                                 if j == 1;
                                   polynomial[j] = one(T);
                                 elseif j == 2;
                                   polynomial[j] = xi;
                                 else;
                                   polynomial[j] = 2*xi*polynomial[j-1]-polynomial[j-2];
                                 end;
                               end;
                               push!(poly,polynomial);
                             end
                             )

  i_vars = Array{Symbol}(N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
  end

  inner_prod = :( poly[1][i1]*weights[$(i_vars...)] )
  for i = 2:N
    inner_prod = :( poly[$i][$(i_vars[i])]*$inner_prod )
  end

  inner = :( if sum(tuple($(i_vars...))) <= order+N;
               evaluated_polynomial += $inner_prod;
             end )
  outer = inner
  for i = 1:N
    outer = :(
      for $(i_vars[i]) = 1:size(weights,$i)
        $outer
      end
    )
  end

  final = :( $chebyshev_polynomials;
             evaluated_polynomial = zero(T);
             $outer;
             return evaluated_polynomial
          )

  return final

end
