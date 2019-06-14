# Generated functions for tensor-product polynomials where nodes are in an array of arrays

@generated function chebyshev_weights(f::AbstractArray{T,N},nodes::Union{Array{Array{T,1},1},NTuple{N,Array{T,1}}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T,N,S}

  chebyshev_polynomials = :( poly = Array{Array{T,2},1}(undef,N);
                             for k = 1:N;
                               orderk = order[k];
                               xk = nodes[k];

                               # normalize nodes

                               if domain[1,k] == domain[2,k];
                                 fill!(xk,(domain[1,k].+domain[2,k])/2);
                               else;
                                 xk = 2*(xk.-domain[2,k])/(domain[1,k].-domain[2,k]).-one(T);
                               end;

                               polynomial = Array{T}(undef,length(xk),orderk+1);
                               for i = 1:orderk+1;
                                 for j = 1:length(xk);
                                   if i == 1;
                                     polynomial[j,i] = one(T);
                                   elseif i == 2;
                                     polynomial[j,i] = xk[j];
                                   else;
                                     polynomial[j,i] = 2*xk[j]*polynomial[j,i-1]-polynomial[j,i-2];
                                   end;
                                 end;
                               end;
                               poly[i] = polynomial;#push!(poly,polynomial);
                             end )

  i_vars = Array{Symbol}(undef,N)
  s_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
    s_vars[i] = Symbol("s$i")
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

  for i = N:-1:1#1:N
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

  for i = N:-1:1#1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order[$i]+1)
        $new_outer
      end
    )
  end

  # Put it all together to compute the weights array

  final = :( $chebyshev_polynomials;
             weights = Array{T}(undef,(order+1)...);
             $new_outer;
             return weights )

  return final

end

@generated function chebyshev_weights(f::AbstractArray{T,N},poly::Union{Array{Array{T,2},1},NTuple{N,Array{T,2}}},order::Array{S,1}) where {T,N,S}

  i_vars = Array{Symbol}(undef,N)
  s_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
    s_vars[i] = Symbol("s$i")
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

  for i = N:-1:1#1:N
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

  for i = N:-1:1#1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order[$i]+1)
        $new_outer
      end
    )
  end

    # Put it all together to compute the weights array

  final = :( weights = Array{T}(undef,(order.+1)...);
             $new_outer;
             return weights )

  return final

end

# Generated functions for tensor-product polynomials where nodes are in a tuple

#=

@generated function chebyshev_weights(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T,N,S}

  chebyshev_polynomials = :( poly = Array{Array{T,2},1}(undef,N);
                             for k = 1:N;
                               orderk = order[k];
                               xk = nodes[k];

                               # normalize nodes

                               if domain[1,k] == domain[2,k];
                                 fill!(xk,(domain[1,k].+domain[2,k])/2);
                               else;
                                 xk = 2*(xk.-domain[2,k])/(domain[1,k].-domain[2,k]).-one(T);
                               end;

                               polynomial = Array{T}(undef,length(xk),orderk+1);
                               for i = 1:orderk+1;
                                 for j = 1:length(xk);
                                   if i == 1;
                                     polynomial[j,i] = one(T);
                                   elseif i == 2;
                                     polynomial[j,i] = xk[j];
                                   else;
                                     polynomial[j,i] = 2*xk[j]*polynomial[j,i-1]-polynomial[j,i-2];
                                   end;
                                 end;
                               end;
                               poly[k] = polynomial;#push!(poly,polynomial);
                             end )

  i_vars = Array{Symbol}(undef,N)
  s_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
    s_vars[i] = Symbol("s$i")
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

  for i = N:-1:1#1:N
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

  for i = N:-1:1#1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order[$i]+1)
        $new_outer
      end
    )
  end

  # Put it all together to compute the weights array

  final = :( $chebyshev_polynomials;
             weights = Array{T}(undef,(order.+1)...);
             $new_outer;
             return weights )

  return final

end

@generated function chebyshev_weights(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::Array{S,1}) where {T,N,S}

  i_vars = Array{Symbol}(undef,N)
  s_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
    s_vars[i] = Symbol("s$i")
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

  for i = N:-1:1#1:N
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

  for i = N:-1:1#1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order[$i]+1)
        $new_outer
      end
    )
  end

  # Put it all together to compute the weights array

  final = :( weights = Array{T}(undef,(order.+1)...);
             $new_outer;
             return weights )

  return final

end

=#

# Generated functions for complete polynomials where nodes are in an array of arrays

@generated function chebyshev_weights(f::AbstractArray{T,N},nodes::Union{Array{Array{T,1},1},NTuple{N,Array{T,1}}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T,N,S}

  chebyshev_polynomials = :( poly = Array{Array{T,2},1}(undef,N);
                             for k = 1:N;
                               xk = nodes[k];

                               # normalize nodes

                               if domain[1,k] == domain[2,k];
                                 fill!(xk,(domain[1,k].+domain[2,k])/2);
                               else;
                                 xk = 2*(xk.-domain[2,k])/(domain[1,k].-domain[2,k]).-one(T);
                               end;

                               polynomial = Array{T}(undef,length(xk),order+1);
                               for i = 1:order+1;
                                 for j = 1:length(xk);
                                   if i == 1;
                                     polynomial[j,i] = one(T);
                                   elseif i == 2;
                                     polynomial[j,i] = xk[j];
                                   else;
                                     polynomial[j,i] = 2*xk[j]*polynomial[j,i-1]-polynomial[j,i-2];
                                   end;
                                 end;
                               end;
                               poly[k] = polynomial;#push!(poly,polynomial);
                             end )

  i_vars = Array{Symbol}(undef,N)
  s_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
    s_vars[i] = Symbol("s$i")
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

  for i = N:-1:1#1:N
    outer = :(
      for $(s_vars[i]) = 1:size(f,$i)
        $outer
      end
    )
  end

  # Construct the outer loops that compute the weights

  new_inner = :( numerator = zero(T);
                 denominator = zero(T);
                 if sum(tuple($(i_vars...))) <= order+N;
                   $outer;
                   $weights_term;
                 end )
  new_outer = new_inner

  for i = N:-1:1#1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order+1)
        $new_outer
      end
    )
  end

  # Initialize the weight Array

#=

  initial_weight = string("weights = zeros(order+1,",)
  for i = 2:N
    initial_weight = string(initial_weight,"order+1,")
  end
  initial_weight = parse(string(initial_weight,")"))

=#

  # Put it all together to compute the weights array

  final = :( $chebyshev_polynomials;
             weights = Array{T}(undef,(order.+1)...);
#             $initial_weight;
             $new_outer;
             return weights )

  return final

end

@generated function chebyshev_weights(f::AbstractArray{T,N},poly::Union{Array{Array{T,2},1},NTuple{N,Array{T,2}}},order::S) where {T,N,S}

  i_vars = Array{Symbol}(undef,N)
  s_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
    s_vars[i] = Symbol("s$i")
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

  for i = N:-1:1#1:N
    outer = :(
      for $(s_vars[i]) = 1:size(f,$i)
        $outer
      end
    )
  end

  # Construct the outer loops that compute the weights

  new_inner = :( numerator = zero(T);
                 denominator = zero(T);
                 if sum(tuple($(i_vars...))) <= order+N;
                   $outer;
                   $weights_term;
                 end )
  new_outer = new_inner

  for i = N:-1:1#1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order+1)
        $new_outer
      end
    )
  end

  # Initialize the weight Array

#=

  initial_weight = string("weights = zeros(order+1,",)
  for i = 2:N
    initial_weight = string(initial_weight,"order+1,")
  end
  initial_weight = parse(string(initial_weight,")"))

=#

  # Put it all together to compute the weights array

  final = :( weights = Array{T}(undef,(order.+1)...);
#             $initial_weight;
             $new_outer;
             return weights )

  return final

end

#=

# Generated functions for complete polynomials where nodes are in a tuple

@generated function chebyshev_weights(f::AbstractArray{T,N},nodes::NTuple{N,Array{T,1}},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T,N,S}

  chebyshev_polynomials = :( poly = Array{Array{T,2},1}(undef,N);
                             for k = 1:N;
                               xk = nodes[k];

                               # normalize nodes

                               if domain[1,k] == domain[2,k];
                                 fill!(xk,(domain[1,k].+domain[2,k])/2);
                               else;
                                 xk = 2*(xk.-domain[2,k])/(domain[1,k].-domain[2,k]).-one(T);
                               end;

                               polynomial = Array{T}(undef,length(xk),order+1);
                               for i = 1:order+1;
                                 for j = 1:length(xk);
                                   if i == 1;
                                     polynomial[j,i] = one(T);
                                   elseif i == 2;
                                     polynomial[j,i] = xk[j];
                                   else;
                                     polynomial[j,i] = 2*xk[j]*polynomial[j,i-1]-polynomial[j,i-2];
                                   end;
                                 end;
                               end;
                               poly[k] = polynomial;#push!(poly,polynomial);
                             end )

  i_vars = Array{Symbol}(undef,N)
  s_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
    s_vars[i] = Symbol("s$i")
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

  for i = N:-1:1#1:N
    outer = :(
      for $(s_vars[i]) = 1:size(f,$i)
        $outer
      end
    )
  end

  # Construct the outer loops that compute the weights

  new_inner = :( numerator = zero(T);
                 denominator = zero(T);
                 if sum(tuple($(i_vars...))) <= order+N;
                   $outer;
                   $weights_term;
                 end )
  new_outer = new_inner

  for i = N:-1:1#1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order+1)
        $new_outer
      end
    )
  end

  # Initialize the weight Array

  initial_weight = string("weights = zeros(order+1,",)
  for i = 2:N
    initial_weight = string(initial_weight,"order+1,")
  end
  initial_weight = Meta.parse(string(initial_weight,")"))

  # Put it all together to compute the weights array

  final = :( $chebyshev_polynomials;
#             weights = Array{T}(undef,(order.+1)...);
             $initial_weight;
             $new_outer;
             return weights )

  return final

end

@generated function chebyshev_weights(f::AbstractArray{T,N},poly::NTuple{N,Array{T,2}},order::S) where {T,N,S}

  i_vars = Array{Symbol}(undef,N)
  s_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
    s_vars[i] = Symbol("s$i")
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

  for i = N:-1:1#1:N
    outer = :(
      for $(s_vars[i]) = 1:size(f,$i)
        $outer
      end
    )
  end

  # Construct the outer loops that compute the weights

  new_inner = :( numerator = zero(T);
                 denominator = zero(T);
                 if sum(tuple($(i_vars...))) <= order+N;
                   $outer;
                   $weights_term;
                 end )
  new_outer = new_inner

  for i = N:-1:1#1:N
    new_outer = :(
      for $(i_vars[i]) = 1:(order+1)
        $new_outer
      end
    )
  end

  # Initialize the weight Array

  initial_weight = string("weights = zeros(order+1,",)
  for i = 2:N
    initial_weight = string(initial_weight,"order+1,")
  end
  initial_weight = Meta.parse(string(initial_weight,")"))

  # Put it all together to compute the weights array

  final = :( #weights = Array{T}(undef,(order.+1)...);
             $initial_weight;
             $new_outer;
             return weights )

  return final

end

=#

# Functions to ensure backward compatibility with an older API

# Tensor product polynomials

function chebyshev_weights(f::AbstractArray{T,1},nodes_1::Array{T,1},order::Array{S,1},domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes_1,),order,reshape(domain,2,1))

  return weights

end

function chebyshev_weights(f::AbstractArray{T,2},nodes_1::Array{T,1},nodes_2::Array{T,1},order::Array{S,1},domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes_1,nodes_2),order,domain)

  return weights

end

function chebyshev_weights(f::AbstractArray{T,3},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},order::Array{S,1},domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes_1,nodes_2,nodes_3),order,domain)

  return weights

end

function chebyshev_weights(f::AbstractArray{T,4},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},order::Array{S,1},domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes_1,nodes_2,nodes_3,nodes_4),order,domain)

  return weights

end

function chebyshev_weights(f::AbstractArray{T,5},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},order::Array{S,1},domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes_1,nodes_2,nodes_3,nodes_4,nodes_5),order,domain)

  return weights

end

function chebyshev_weights(f::AbstractArray{T,6},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},nodes_6::Array{T,1},order::Array{S,1},domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes_1,nodes_2,nodes_3,nodes_4,nodes_5,nodes_6),order,domain)

  return weights

end

# Complete polynomials

function chebyshev_weights(f::AbstractArray{T,1},nodes_1::Array{T,1},order::S,domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes_1,),order,reshape(domain,2,1))

  return weights

end

function chebyshev_weights(f::AbstractArray{T,2},nodes_1::Array{T,1},nodes_2::Array{T,1},order::S,domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes_1,nodes_2),order,domain)

  return weights

end

function chebyshev_weights(f::AbstractArray{T,3},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},order::S,domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes_1,nodes_2,nodes_3),order,domain)

  return weights

end

function chebyshev_weights(f::AbstractArray{T,4},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},order::S,domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes_1,nodes_2,nodes_3,nodes_4),order,domain)

  return weights

end

function chebyshev_weights(f::AbstractArray{T,5},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},order::S,domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes_1,nodes_2,nodes_3,nodes_4,nodes_5),order,domain)

  return weights

end

function chebyshev_weights(f::AbstractArray{T,6},nodes_1::Array{T,1},nodes_2::Array{T,1},nodes_3::Array{T,1},nodes_4::Array{T,1},nodes_5::Array{T,1},nodes_6::Array{T,1},order::S,domain=[ones(1);-ones(1)]) where {T<:AbstractFloat,S<:Integer}

  weights = chebyshev_weights(f,(nodes_1,nodes_2,nodes_3,nodes_4,nodes_5,nodes_6),order,domain)

  return weights

end
