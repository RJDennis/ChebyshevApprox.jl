# Generated functions for evaluating Chebyshev polynomials

# Tensor product polynomials

@generated function chebyshev_derivative(weights::Array{T,N},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T,N,S}

  chebyshev_polynomials = :( poly        = Array{Array{T,2},1}(undef,N);
                             poly_derivs = Array{Array{T,2},1}(undef,N);
                             for i = 1:size(x,1);
                               xi = x[i];

                               # Normalize nodes to [-1,1]

                               if domain[1,i] == domain[2,i];
                                 xi = (domain[1,i].+domain[2,i])/2;
                               else;
                                 xi = 2*(xi.-domain[2,i])/(domain[1,i].-domain[2,i]).-one(T);
                               end;

                               polynomial = ones(T,1,order[i]+1);
                               poly_deriv = zeros(T,1,order[i]+1);
                               for j = 2:order[i]+1;
                                 if j == 2;
                                   polynomial[j] = xi;
                                   poly_deriv[j] = one(T);
                                 else;
                                   polynomial[j] = 2*xi*polynomial[j-1]-polynomial[j-2];
                                   poly_deriv[j] = ((j-1)*polynomial[j-1]-(j-1)*xi*polynomial[j])/(1-xi^2)
                                 end;
                               end;
                               poly[i]        = polynomial;#push!(poly,polynomial);
                               poly_derivs[i] = poly_deriv;#push!(poly_derivs,poly_deriv);
                             end
                             )

  i_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
  end

  j = Symbol("j")

  inner_prod = :( (($j!==1)*poly[1][i1]+($j==1)*poly_derivs[1][i1])*weights[$(i_vars...)] )
  for i = 2:N
    inner_prod = :( (($j!==$i)*poly[$i][$(i_vars[i])]+($j==$i)*poly_derivs[$i][$(i_vars[i])])*$inner_prod )
  end

  inner = :( evaluated_derivatives[$j] += (2.0/(domain[1,$j].-domain[2,$j]))*$inner_prod )
  outer = inner
  for i = N:-1:1#1:N
    outer = :(
      for $(i_vars[i]) = 1:size(weights,$i)
        $outer
      end
    )
  end

  new_outer = :(
               for j = 1:N
                 $outer
               end
               )

  final = :( $chebyshev_polynomials;
             evaluated_derivatives = zeros(1,N);
             $new_outer;
             return evaluated_derivatives
             )

  return final

end

# Complete polynomials

@generated function chebyshev_derivative(weights::Array{T,N},x::Array{T,1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T,N,S}

  chebyshev_polynomials = :( poly        = Array{Array{T,2},1}(undef,N);
                             poly_derivs = Array{Array{T,2},1}(undef,N);
                             for i = 1:size(x,1);
                               xi = x[i];

                               # Normalize nodes to [-1,1]

                               if domain[1,i] == domain[2,i];
                                 xi = (domain[1,i].+domain[2,i])/2;
                               else;
                                 xi = 2*(xi.-domain[2,i])/(domain[1,i].-domain[2,i]).-one(T);
                               end;

                               polynomial = ones(T,1,order+1);
                               poly_deriv = zeros(T,1,order+1);
                               for j = 2:order+1;
                                 if j == 2;
                                   polynomial[j] = xi;
                                   poly_deriv[j] = one(T);
                                 else;
                                   polynomial[j] = 2*xi*polynomial[j-1]-polynomial[j-2];
                                   poly_deriv[j] = ((j-1)*polynomial[j-1]-(j-1)*xi*polynomial[j])/(1-xi^2)
                                 end;
                               end;
                               poly[i]        = polynomial;#push!(poly,polynomial);
                               poly_derivs[i] = poly_deriv;#push!(poly_derivs,poly_deriv);
                             end
                             )

  i_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
  end

  j = Symbol("j")

  inner_prod = :( (($j!==1)*poly[1][i1]+($j==1)*poly_derivs[1][i1])*weights[$(i_vars...)] )
  for i = 2:N
    inner_prod = :( (($j!==$i)*poly[$i][$(i_vars[i])]+($j==$i)*poly_derivs[$i][$(i_vars[i])])*$inner_prod )
  end

  inner = :( if sum(tuple($(i_vars...))) <= order+N;
               evaluated_derivatives[$j] += (2.0/(domain[1,$j].-domain[2,$j]))*$inner_prod;
             end )

  outer = inner
  for i = N:-1:1#1:N
    outer = :(
      for $(i_vars[i]) = 1:size(weights,$i)
        $outer
      end
    )
  end

  new_outer = :(
               for j = 1:N
                 $outer
               end
               )

 final = :( $chebyshev_polynomials;
            evaluated_derivatives = zeros(1,N);
            $new_outer;
            return evaluated_derivatives
          )

  return final

end

# Generated functions for differentiating with respect to specific variables

# Tensor product polynomials

@generated function chebyshev_derivative(weights::Array{T,N},x::Array{T,1},order::Array{S,1},pos::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T,N,S}

  chebyshev_polynomials = :( poly        = Array{Array{T,2},1}(undef,N);
                             poly_derivs = Array{Array{T,2},1}(undef,N);
                             for i = 1:size(x,1);
                               xi = x[i];

                               # Normalize nodes to [-1,1]

                               if domain[1,i] == domain[2,i];
                                 xi = (domain[1,i].+domain[2,i])/2;
                               else;
                                 xi = 2*(xi.-domain[2,i])/(domain[1,i].-domain[2,i]).-one(T);
                               end;

                               polynomial = ones(T,1,order[i]+1);
                               poly_deriv = zeros(T,1,order[i]+1);
                               for j = 2:order[i]+1;
                                 if j == 2;
                                   polynomial[j] = xi;
                                   poly_deriv[j] = one(T);
                                 else;
                                   polynomial[j] = 2*xi*polynomial[j-1]-polynomial[j-2];
                                   poly_deriv[j] = ((j-1)*polynomial[j-1]-(j-1)*xi*polynomial[j])/(1-xi^2)
                                 end;
                               end;
                               poly[i]        = polynomial;#push!(poly,polynomial);
                               poly_derivs[i] = poly_deriv;#push!(poly_derivs,poly_deriv);
                             end
                             )

  i_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
  end

  j = Symbol("j")

  inner_prod = :( ((pos[$j]!==1)*poly[1][i1]+(pos[$j]==1)*poly_derivs[1][i1])*weights[$(i_vars...)] )
  for i = 2:N
    inner_prod = :( ((pos[$j]!==$i)*poly[$i][$(i_vars[i])]+(pos[$j]==$i)*poly_derivs[$i][$(i_vars[i])])*$inner_prod )
  end

  inner = :( evaluated_derivatives[$j] += (2.0/(domain[1,pos[$j]].-domain[2,pos[$j]]))*$inner_prod )
  outer = inner
  for i = N:-1:1#1:N
    outer = :(
      for $(i_vars[i]) = 1:size(weights,$i)
        $outer
      end
    )
  end

  new_outer = :(
               for j = 1:length(pos)
                 $outer
               end
               )

  final = :( $chebyshev_polynomials;
             evaluated_derivatives = zeros(1,length(pos));
             $new_outer;
             return evaluated_derivatives
             )

  return final

end

# Complete polynomial

@generated function chebyshev_derivative(weights::Array{T,N},x::Array{T,1},order::S,pos::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T,N,S}

  chebyshev_polynomials = :( poly        = Array{Array{T,2},1}(undef,N);
                             poly_derivs = Array{Array{T,2},1}(undef,N);
                             for i = 1:size(x,1);
                               xi = x[i];

                               # Normalize nodes to [-1,1]

                               if domain[1,i] == domain[2,i];
                                 xi = (domain[1,i].+domain[2,i])/2;
                               else;
                                 xi = 2*(xi.-domain[2,i])/(domain[1,i].-domain[2,i]).-one(T);
                               end;

                               polynomial = ones(T,1,order+1);
                               poly_deriv = zeros(T,1,order+1);
                               for j = 2:order+1;
                                 if j == 2;
                                   polynomial[j] = xi;
                                   poly_deriv[j] = one(T);
                                 else;
                                   polynomial[j] = 2*xi*polynomial[j-1]-polynomial[j-2];
                                   poly_deriv[j] = ((j-1)*polynomial[j-1]-(j-1)*xi*polynomial[j])/(1-xi^2)
                                 end;
                               end;
                               poly[i]        = polynomial;#push!(poly,polynomial);
                               poly_derivs[i] = poly_deriv;#push!(poly_derivs,poly_deriv);
                             end
                             )

  i_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
  end

  j = Symbol("j")

  inner_prod = :( ((pos[$j]!==1)*poly[1][i1]+(pos[$j]==1)*poly_derivs[1][i1])*weights[$(i_vars...)] )
  for i = 2:N
    inner_prod = :( ((pos[$j]!==$i)*poly[$i][$(i_vars[i])]+(pos[$j]==$i)*poly_derivs[$i][$(i_vars[i])])*$inner_prod )
  end

  inner = :( if sum(tuple($(i_vars...))) <= order+N;
               evaluated_derivatives[$j] += (2.0/(domain[1,pos[$j]].-domain[2,pos[$j]]))*$inner_prod;
             end )

  outer = inner
  for i = N:-1:1#1:N
    outer = :(
      for $(i_vars[i]) = 1:size(weights,$i)
        $outer
      end
    )
  end

  new_outer = :(
               for j = 1:length(pos)
                 $outer
               end
               )

 final = :( $chebyshev_polynomials;
            evaluated_derivatives = zeros(1,length(pos));
            $new_outer;
            return evaluated_derivatives
          )

  return final

end

function chebyshev_derivative(weights::Array{T,N},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_derivative(weights,x,order,domain)

  end

  return goo

end

function chebyshev_derivative(weights::Array{T,N},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_derivative(weights,x,order,domain)

  end

  return goo

end

function chebyshev_derivative(weights::Array{T,N},order::Array{S,1},pos::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_derivative(weights,x,order,pos,domain)

  end

  return goo

end

function chebyshev_derivative(weights::Array{T,N},order::S,pos::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function goo(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_derivative(weights,x,order,pos,domain)

  end

  return goo

end
