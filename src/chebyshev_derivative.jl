# Regular functions for evaluating derivatives and gradients of Chebyshev polynomials

# Regular functions for derivatives

function chebyshev_derivative(weights::Array{T,N},x::Array{T,1},pos::S,order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)
  @inbounds for i = 1:N
    if i === pos
      poly[i] = derivative_of_chebyshev_polynomial(order[i],normalize_node(x[i],domain[:,i]))
    else
      poly[i] = chebyshev_polynomial(order[i],normalize_node(x[i],domain[:,i]))
    end
  end

  derivative = 0.0
  @inbounds for i in CartesianIndices(weights)
    poly_product = poly[1][i[1]]
    @inbounds for j = 2:N
      poly_product *= poly[j][i[j]]
    end
    derivative += weights[i]*poly_product
  end

  return derivative*(2.0/(domain[1,pos]-domain[2,pos]))

end

function chebyshev_derivative(weights::Array{T,N},x::Array{T,1},pos::S,order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)
  @inbounds for i = 1:N
    if i === pos
      poly[i] = derivative_of_chebyshev_polynomial(order,normalize_node(x[i],domain[:,i]))
    else
      poly[i] = chebyshev_polynomial(order,normalize_node(x[i],domain[:,i]))
    end
  end

  derivative = 0.0
  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
      poly_product = poly[1][i[1]]
      @inbounds for j = 2:N
        poly_product *= poly[j][i[j]]
      end
      derivative += weights[i]*poly_product
    end
  end

  return derivative*(2.0/(domain[1,pos]-domain[2,pos]))

end

function chebyshev_derivative(cheb_poly::ChebyshevPolynomial,x::Array{T,1},pos::S) where {T<:AbstractFloat,S<:Integer}

  derivative = chebyshev_derivative(cheb_poly.weights,x,pos,cheb_poly.order,cheb_poly.domain) 
    
  return derivative

end

# Regular functions for gradients

function chebyshev_gradient(weights::Array{T,N},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  gradient = Array{T,2}(undef,1,N)

  @inbounds for i = 1:N
    gradient[i] = chebyshev_derivative(weights,x,i,order,domain)
  end

  return gradient

end

function chebyshev_gradient(weights::Array{T,N},x::Array{T,1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  gradient = Array{T,2}(undef,1,N)

  @inbounds for i = 1:N
    gradient[i] = chebyshev_derivative(weights,x,i,order,domain)
  end

  return gradient

end

function chebyshev_gradient(cheb_poly::ChebyshevPolynomial,x::Array{T,1}) where {T<:AbstractFloat}

  gradient = chebyshev_gradient(cheb_poly.weights,x,cheb_poly.order,cheb_poly.domain) 
    
  return gradient

end

function chebyshev_derivative(weights::Array{T,N},pos::S,order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function chebderiv(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_derivative(weights,x,pos,order,domain)

  end

  return chebderiv

end

function chebyshev_derivative(weights::Array{T,N},pos::S,order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function chebderiv(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_derivative(weights,x,pos,order,domain)

  end

  return chebderiv

end

function chebyshev_derivative(cheb_poly::ChebyshevPolynomial,pos::S) where {S <: Integer}

  function chebderiv(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_derivative(cheb_poly,x,pos)

  end

  return chebderiv

end

function chebyshev_gradient(weights::Array{T,N},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function chebgrad(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_gradient(weights,x,pos,order,domain)

  end

  return chebgrad

end

function chebyshev_gradient(weights::Array{T,N},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function chebgrad(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_gradient(weights,x,pos,order,domain)

  end

  return chebgrad

end

function chebyshev_gradient(cheb_poly::ChebyshevPolynomial) where {S <: Integer}

  function chebgrad(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_gradient(cheb_poly,x)

  end

  return chebgrad

end

##########################################################

# Generated functions for evaluating Chebyshev polynomials

# Tensor product polynomials

#=

@generated function chebyshev_derivative(weights::Array{T,N},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T,N,S}

  chebyshev_polynomials = :( poly        = Array{Array{T,2},1}(undef,N);
                             poly_derivs = Array{Array{T,2},1}(undef,N);
                             for i = 1:size(x,1);
                               xi = x[i];

                               # Normalize nodes to [-1,1]

                               if domain[1,i] == domain[2,i];
                                 xi = zero(T);
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
                               poly[i]        = polynomial;
                               poly_derivs[i] = poly_deriv;
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
                                 xi = zero(T);
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
                               poly[i]        = polynomial;
                               poly_derivs[i] = poly_deriv;
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
                                 xi = zero(T);
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
                               poly[i]        = polynomial;
                               poly_derivs[i] = poly_deriv;
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
  for i = N:-1:1
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
                                 xi = zero(T);
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
                               poly[i]        = polynomial;
                               poly_derivs[i] = poly_deriv;
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
  for i = N:-1:1
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

=#