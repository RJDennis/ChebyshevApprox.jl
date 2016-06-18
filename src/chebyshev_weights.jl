# Generated functions for tensor-product polynomials

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

@generated function chebyshev_weights{T,N,S}(f::Array{T,N},nodes::Array{Array{T,1},1},order::Array{S,1})

  chebyshev_polynomials = :( poly = Array{T,2}[];
                             for k = 1:N;
                               orderk = order[k];
                               xk = nodes[k];
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

# Generated functions for complete polynomials

@generated function chebyshev_weights{T,N,S}(f::Array{T,N},nodes::Array{Array{T,1},1},order::S,range::Array{T,2})

  chebyshev_polynomials = :( poly = Array{T,2}[];
                             for k = 1:N;
                               xk = nodes[k];

                               # normalize nodes

                               if range[1,k] == range[2,k];
                                 fill!(xk,(range[1,k]+range[2,k])/2);
                               else;
                                 xk = 2*(xk-range[2,k])/(range[1,k]-range[2,k])-one(T);
                               end;

                               polynomial = Array(T,length(xk),order+1);
                               for i = 1:order+1;
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
                 if sum([$(i_vars...)]) <= order+N;
                   $outer;
                   $weights_term;
                 end )
  new_outer = new_inner

  for i = 1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order+1)
        $new_outer
      end
    )
  end

  # Initialize the weight Array

  initial_weight = string("weights = zeros(T,order+1,",)
  for i = 2:N
    initial_weight = string(initial_weight,"order+1,")
  end
  initial_weight = parse(string(initial_weight,")"))

    # Put it all together to compute the weights array

  final = :( $chebyshev_polynomials;
             $initial_weight;
             $new_outer;
             return weights )

  return final

end

@generated function chebyshev_weights{T,N,S}(f::Array{T,N},nodes::Array{Array{T,1},1},order::S)

  chebyshev_polynomials = :( poly = Array{T,2}[];
                             for k = 1:N;
                               xk = nodes[k];
                               polynomial = Array(T,length(xk),order+1);
                               for i = 1:order+1;
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
                 if sum([$(i_vars...)]) <= order+N;
                   $outer;
                   $weights_term;
                 end )
  new_outer = new_inner

  for i = 1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order+1)
        $new_outer
      end
    )
  end

  # Initialize the weight Array

  initial_weight = string("weights = zeros(T,order+1,",)
  for i = 2:N
    initial_weight = string(initial_weight,"order+1,")
  end
  initial_weight = parse(string(initial_weight,")"))

    # Put it all together to compute the weights array

  final = :( $chebyshev_polynomials;
             $initial_weight;
             $new_outer;
             return weights )

  return final

end

@generated function chebyshev_weights{T,N,S}(f::Array{T,N},poly::Array{Array{T,2},1},order::S)

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
                 if sum([$(i_vars...)]) <= order+N;
                   $outer;
                   $weights_term;
                 end )
  new_outer = new_inner

  for i = 1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order+1)
        $new_outer
      end
    )
  end

  # Initialize the weight Array

  initial_weight = string("weights = zeros(T,order+1,",)
  for i = 2:N
    initial_weight = string(initial_weight,"order+1,")
  end
  initial_weight = parse(string(initial_weight,")"))

    # Put it all together to compute the weights array

  final = :( $initial_weight;
             $new_outer;
             return weights )

  return final

end
