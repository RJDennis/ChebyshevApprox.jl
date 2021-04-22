# Regular functions for evaluating Chebyshev polynominals

function chebyshev_evaluate(weights::Array{T,N},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i],normalize_node(x[i],domain[:,i]))
  end

  yhat = 0.0
  @inbounds for i in CartesianIndices(weights)
    poly_product = poly[1][i[1]]
    @inbounds for j = 2:N
      poly_product *= poly[j][i[j]]
    end
    yhat += weights[i]*poly_product
  end

  return yhat

end

function chebyshev_evaluate(weights::Array{T,N},x::Array{T,1},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef,N)
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order,normalize_node(x[i],domain[:,i]))
  end

  yhat = 0.0
  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order+N
      poly_product = poly[1][i[1]]
      @inbounds for j = 2:N
        poly_product *= poly[j][i[j]]
      end
      yhat += weights[i]*poly_product
    end
  end

  return yhat

end

function chebyshev_evaluate(cheb_poly::ChebyshevPolynomial,x::Array{T,1}) where {T<:AbstractFloat}

  yhat = chebyshev_evaluate(cheb_poly.weights,x,cheb_poly.order,cheb_poly.domain) 
    
  return yhat
  
end

function chebyshev_evaluate(weights::Array{T,N},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function chebeval(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_evaluate(weights,x,order,domain)

  end

  return chebeval

end

function chebyshev_evaluate(weights::Array{T,N},order::S,domain=[ones(T,1,N);-ones(T,1,N)]) where {T <: AbstractFloat,N,S <: Integer}

  function chebeval(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_evaluate(weights,x,order,domain)

  end

  return chebeval

end

function chebyshev_evaluate(cheb_poly::ChebyshevPolynomial)

  function chebeval(x::Array{T,1}) where {T <: AbstractFloat}

    return chebyshev_evaluate(cheb_poly,x)

  end

  return chebeval

end

##########################################################

# Generated functions for evaluating Chebyshev polynomials

# In the complete polynomial case the @generated function is faster than the regular function

#=

@generated function chebyshev_evaluate(weights::Array{T,N},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,N);-ones(T,1,N)]) where {T,N,S}

  chebyshev_polynomials = :( poly = Array{Array{T,2},1}(undef,N);
                             for i = 1:size(x,1);
                               xi = x[i];

                               # Normalize nodes to [-1,1]

                               if domain[1,i] == domain[2,i];
                                 xi = zero(T)
                               else;
                                 xi = 2*(xi.-domain[2,i])/(domain[1,i].-domain[2,i]).-one(T);
                               end;

                               polynomial    = Array{T}(undef,1,order[i]+1);
                               polynomial[1] = one(T);
                               for j = 2:order[i]+1;
                                 if j == 2;
                                   polynomial[j] = xi;
                                 else;
                                   polynomial[j] = 2*xi*polynomial[j-1]-polynomial[j-2];
                                 end;
                               end;
                               poly[i] = polynomial;
                             end
                             )

  i_vars = Array{Symbol}(undef,N)
  for i = 1:N
    i_vars[i] = Symbol("i$i")
  end

  inner_prod = :( poly[1][i1]*weights[$(i_vars...)] )
  for i = 2:N
    inner_prod = :( poly[$i][$(i_vars[i])]*$inner_prod )
  end

  inner = :( evaluated_polynomial += $inner_prod )
  outer = inner
  for i = N:-1:1
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

  chebyshev_polynomials = :( poly = Array{Array{T,2},1}(undef,N);
                             for i = 1:size(x,1);
                               xi = x[i];

                               # Normalize nodes to [-1,1]

                               if domain[1,i] == domain[2,i];
                                 xi = zero(T)
                               else;
                                 xi = 2*(xi.-domain[2,i])/(domain[1,i].-domain[2,i]).-one(T);
                               end;

                               polynomial    = Array{T}(undef,1,order+1);
                               polynomial[1] = one(T);
                               for j = 2:order+1;
                                 if j == 2;
                                   polynomial[j] = xi;
                                 else;
                                   polynomial[j] = 2*xi*polynomial[j-1]-polynomial[j-2];
                                 end;
                               end;
                               poly[i] = polynomial;
                             end
                             )

  i_vars = Array{Symbol}(undef,N)
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
  for i = N:-1:1
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

=#