abstract type Nodes end
abstract type CApproximationPlan end

struct ChebRoots{T<:AbstractFloat,R<:AbstractFloat} <: Nodes

  points::Array{R,1}
  domain::Array{T,1}

end

struct ChebExtrema{T<:AbstractFloat,R<:AbstractFloat} <: Nodes

  points::Array{R,1}
  domain::Array{T,1}

end

struct ChebExtended{T<:AbstractFloat,R<:AbstractFloat} <: Nodes

  points::Array{R,1}
  domain::Array{T,1}

end

struct Grid{G<:Nodes}

  grid::Tuple{Vararg{G}}

end

struct ChebPoly{T<:AbstractFloat} # Holds polynomials as well as their first- and second-derivatives

  poly::Union{Array{T,2},Transpose{T,Array{T,2}}}
  nodetype::DataType

end

struct CApproxPlan{G<:Grid,S<:Integer,T<:AbstractFloat} <: CApproximationPlan

  grid::G
  order::Union{S,Tuple{Vararg{S}}}
  domain::Union{Array{T,1},Array{T,2}}

end

struct CApproxPlanPoly{N,P<:ChebPoly,S<:Integer,T<:AbstractFloat} <: CApproximationPlan

  polys::NTuple{N,P}
  order::Union{S,Tuple{Vararg{S}}}
  domain::Union{Array{T,1},Array{T,2}}

end

"""
Compute 'N' roots of the Chebyshev polynomial and scale the roots to the interval given in 'domain'.

N --- an integer specifying the number of Chebyshev roots.

domain --- (with with default [1.0,-1.0]) is a 2-element vector specifying the upper and lower bounds on the domain.
"""
function chebyshev_nodes(N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = fill((domain[1]+domain[2])*0.5,N)

  @inbounds for i = 1:div(N,2)
    x = -cos((i - 0.5)*π/N)*(domain[1] - domain[2])*0.5
    points[i]     += x
    points[N-i+1] -= x
  end

  return points

end

"""
Compute 'N' extrema of the Chebyshev polynomial and scale the roots to the interval given in 'domain'.

N --- an integer specifying the number of Chebyshev extrema.

domain --- (with with default [1.0,-1.0]) is a 2-element vector specifying the upper and lower bounds on the domain.
"""
function chebyshev_extrema(N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = fill((domain[1]+domain[2])*0.5,N)

  @inbounds for i = 1:div(N,2)
    x = -cos((i - 1)*π/(N - 1))*(domain[1] - domain[2])*0.5
    points[i]     += x
    points[N-i+1] -= x
  end

  return points

end

"""
Compute 'N' extended Chebyshev points and scale the points to the interval given in 'domain'.

N --- an integer specifying the number of extended Chebyshev points.

domain --- (with with default [1.0,-1.0]) is a 2-element vector specifying the upper and lower bounds on the domain.
"""
function chebyshev_extended(N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = fill((domain[1]+domain[2])*0.5,N)

  @inbounds for i = 1:div(N,2)
    x = -cos((i - 0.5)*π/N)*((domain[1] - domain[2])*0.5)/cos(π/(2N))
    points[i]     += x
    points[N-i+1] -= x
  end

  return points

end

"""
Computes 'N' points according to 'node_generator' and scales those points to the interval given in 'domain'.

N --- an integer specifying the number of points.

node_generator --- a function that computes the desired points.

domain --- (with with default [1.0,-1.0]) is a 2-element vector specifying the upper and lower bounds on the domain.
"""
function nodes(N::S, node_generator::Symbol, domain=[1.0, -1.0],) where {S<:Integer}

  if node_generator == :chebyshev_nodes
    return ChebRoots(chebyshev_nodes(N, domain), domain)
  elseif node_generator == :chebyshev_extrema
    return ChebExtrema(chebyshev_extrema(N, domain), domain)
  elseif node_generator == :chebyshev_extended
    return ChebExtended(chebyshev_extended(N, domain), domain)
  end

end

"""
Normalizes a point, 'node', in 'domain' to the [1.0,-1.0] interval.

node --- a number.

domain --- (with with default [1.0,-1.0]) is a 2-element vector specifying the upper and lower bounds on the domain.
""" 
function normalize_node(node::R, domain::Array{T,1}) where {R<:Number,T<:AbstractFloat}

  if domain[1] == domain[2]
    norm_node = zero(T)
    return norm_node
  else
    norm_node = 2*(node - domain[2])/(domain[1] - domain[2]) - 1
    return norm_node
  end

end

"""
Normalizes a vector of points, 'node', with in element in 'domain' to the [1.0,-1.0] interval.

node --- a vector of numbers.

domain --- (with with default [1.0,-1.0]) is a 2-element vector specifying the upper and lower bounds on the domain.
""" 
function normalize_node(nodes::Array{R,1}, domain::Array{T,1}) where {R<:Number,T<:AbstractFloat}

  norm_nodes = map(x -> normalize_node(x, domain), nodes)

  return norm_nodes

end

"""
Normalizes the points in 'node' to reside on the [1.0,-1.0] interval.

node --- A Nodes type.
"""
function normalize_node(nodes::G) where {G<:Nodes}

  norm_nodes = map(x -> normalize_node(x, nodes.domain), nodes.points)

  return norm_nodes

end

"""
Compute a Chebyshev polynomial of order 'order' at point, 'x'.

order --- an integer specifying the order of the polynomial.

x --- a number.
"""
function chebyshev_polynomial(order::S, x::R) where {S<:Integer,R<:Number}

  # x must reside in [-1,1]

  poly = Array{R}(undef, 1, order + 1)
  poly[1] = one(R)

  @inbounds for i = 2:order+1
    if i == 2
      poly[i] = x
    else
      poly[i] = 2*x*poly[i-1] - poly[i-2]
    end
  end

  return poly

end

"""
Compute a Chebyshev polynomial of order 'order' at points, 'x'.

order --- an integer specifying the order of the polynomial.

x --- a vector of numbers.
"""
function chebyshev_polynomial(order::S, x::AbstractArray{R,1}) where {S<:Integer,R<:Number}

  # Elements of x must reside in [-1,1]

  poly = Array{R}(undef, length(x), order + 1)
  poly[:, 1] .= ones(R, length(x))

  @inbounds for i = 2:order+1
    for j in eachindex(x)
      if i == 2
        poly[j, i] = x[j]
      else
        poly[j, i] = 2*x[j]*poly[j,i-1] - poly[j,i-2]
      end
    end
  end

  return poly

end

"""
Compute a Chebyshev polynomial of order 'order' from the structure, 'g'.

order --- an integer specifying the order of the polynomial.

g --- a Nodes structure.
"""
function chebyshev_polynomial(order::S, nodes::G) where {S<:Integer,G<:Nodes}

  T = eltype(nodes.points)

  poly = Array{T}(undef, length(nodes.points), order + 1)
  poly[:, 1] .= ones(T, length(nodes.points))

  @inbounds for i = 2:order+1
    for j in eachindex(nodes.points)
      if i == 2
        poly[j, i] = normalize_node(nodes.points[j],nodes.domain)
      else
        poly[j, i] = 2*normalize_node(nodes.points[j],nodes.domain)*poly[j,i-1] - poly[j,i-2]
      end
    end
  end

  return ChebPoly(poly, G)

end

"""
Compute the derivatives of a Chebyshev polynomial of order 'order' at point, 'x'.

order --- an integer specifying the order of the polynomial.

x --- a number.
"""
function chebyshev_polynomial_deriv(order::S, x::R) where {S<:Integer,R<:Number}

  poly_deriv = Array{R}(undef, 1, order + 1)
  poly_deriv[1] = zero(R)

  p   = one(R)
  pl  = NaN
  pll = NaN

  @inbounds for i = 2:order+1
    if i == 2
      pl, p = p, x
      poly_deriv[i] = one(R)
    else
      pll, pl = pl, p
      p = 2*x*pl - pll
      poly_deriv[i] = 2*pl + 2*x*poly_deriv[i-1] - poly_deriv[i-2]
    end
  end

  return poly_deriv

end

"""
Compute the derivatives of a Chebyshev polynomial of order 'order' at points, 'x'.

order --- an integer specifying the order of the polynomial.

x --- a vector of numbers.
"""
function chebyshev_polynomial_deriv(order::S, x::AbstractArray{R,1}) where {S<:Integer,R<:Number}

  poly_deriv = Array{R}(undef, order + 1, length(x))
  poly_deriv[1, :] .= zeros(R, length(x))

  @inbounds for j in eachindex(x)
    p = one(R)
    pl = NaN
    pll = NaN
    for i = 2:order+1
      if i == 2
        pl, p = p, x[j]
        poly_deriv[i,j] = one(R)
      else
        pll, pl = pl, p
        p = 2*x[j]*pl - pll
        poly_deriv[i,j] = 2*pl + 2*x[j]*poly_deriv[i-1,j] - poly_deriv[i-2,j]
      end
    end
  end

  return transpose(poly_deriv)# <: AbstractArray

end

"""
Compute the derivatives of a Chebyshev polynomial of order 'order' from the structure, 'g'.

order --- an integer specifying the order of the polynomial.

g --- a Nodes structure.
"""
function chebyshev_polynomial_deriv(order::S, nodes::G) where {S<:Integer,G<:Nodes}

  T = eltype(nodes.points)

  poly_deriv = Array{T}(undef, order + 1,length(nodes.points))
  poly_deriv[1,:] .= zeros(T,length(nodes.points))

  @inbounds for j in eachindex(nodes.points)
    p = one(T)
    pl = NaN
    pll = NaN
    for i = 2:order+1
      if i == 2
        pl, p = p, normalize_node(nodes.points[j],nodes.domain)
        poly_deriv[i,j] = one(T)
      else
        pll, pl = pl, p
        p = 2*normalize_node(nodes.points[j],nodes.domain)*pl - pll
        poly_deriv[i,j] = 2*pl + 2*normalize_node(nodes.points[j],nodes.domain)*poly_deriv[i-1,j] - poly_deriv[i-2,j]
      end
    end
  end

  return ChebPoly(transpose(poly_deriv), G)

end

"""
Compute the second derivatives of a Chebyshev polynomial of order 'order' at point, 'x'.

order --- an integer specifying the order of the polynomial.

x --- a number.
"""
function chebyshev_polynomial_sec_deriv(order::S, x::T) where {T<:Number,S<:Integer}

  poly_sec_deriv = Array{T}(undef, 1, order + 1)
  poly_sec_deriv[1] = zero(T)

  p = one(T)
  pl = NaN
  pll = NaN
  pd = zero(T)
  pdl = NaN
  pdll = NaN

  @inbounds for i = 2:order+1
    if i == 2
      pl, p = p, x
      pdl, pd = pd, one(T)
      poly_sec_deriv[i] = zero(T)
    else
      pll, pl = pl, p
      p = 2*x*pl - pll
      pdll, pdl = pdl, pd
      pd = 2*pl + 2*x*pdl - pdll
      poly_sec_deriv[i] = 2*x*poly_sec_deriv[i-1] + 4*pdl - poly_sec_deriv[i-2]
    end
  end

  return poly_sec_deriv

end

"""
Compute the second derivatives of a Chebyshev polynomial of order 'order' at points, 'x'.

order --- an integer specifying the order of the polynomial.

x --- a vector of numbers.
"""
function chebyshev_polynomial_sec_deriv(order::S, x::AbstractArray{T,1}) where {S<:Integer,T<:Number}

  poly_sec_deriv = Array{T}(undef, order + 1, length(x))
  poly_sec_deriv[1,:] .= zeros(T, length(x))

  @inbounds for j in eachindex(x)
    p = one(T)
    pl = NaN
    pll = NaN
    pd = zero(T)
    pdl = NaN
    pdll = NaN
    for i = 2:order+1
      if i == 2
        pl, p = p, x[j]
        pdl, pd = pd, one(T)
        poly_sec_deriv[i,j] = zero(T)
      else
        pll, pl = pl, p
        p = 2*x[j]*pl - pll
        pdll, pdl = pdl, pd
        pd = 2*pl + 2*x[j]*pdl - pdll
        poly_sec_deriv[i,j] = 2*x[j]*poly_sec_deriv[i-1,j] + 4*pdl - poly_sec_deriv[i-2,j]
      end
    end
  end

  return transpose(poly_sec_deriv)

end

"""
Compute the second derivatives of a Chebyshev polynomial of order 'order' from the structure, 'g'.

order --- an integer specifying the order of the polynomial.

g --- a Nodes structure.
"""
function chebyshev_polynomial_sec_deriv(order::S,nodes::G) where {G<:Nodes,S<:Integer}

  T = eltype(nodes.points)

  poly_sec_deriv = Array{T}(undef, order + 1, length(nodes.points))
  poly_sec_deriv[1,:] .= zeros(T, length(nodes.points))

  @inbounds for j in eachindex(nodes.points)
    p = one(T)
    pl = NaN
    pll = NaN
    pd = zero(T)
    pdl = NaN
    pdll = NaN
    for i = 2:order+1
      if i == 2
        pl, p = p, normalize_node(nodes.points[j], nodes.domain)
        pdl, pd = pd, one(T)
        poly_sec_deriv[i,j] = zero(T)
      else
        pll, pl = pl, p
        p = 2*normalize_node(nodes.points[j],nodes.domain)*pl - pll
        pdll, pdl = pdl, pd
        pd = 2 * pl + 2 * normalize_node(nodes.points[j], nodes.domain) * pdl - pdll
        poly_sec_deriv[i,j] = 2*normalize_node(nodes.points[j],nodes.domain)*poly_sec_deriv[i-1,j] + 4*pdl - poly_sec_deriv[i-2,j]
      end
    end
  end

  return ChebPoly(transpose(poly_sec_deriv), G)

end

"""
Computes the Chebyshev weights given the approximation sample, 'y', and the approximation plan, 'plan'.

y --- a N-dimensional array of data evaluated at the approximation points.

plan --- a structure specifying how the approximation should be undertaken.
""" 
function chebyshev_weights(y::AbstractArray{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  if typeof(plan) <: CApproxPlan

    nodes = Array{Array{T,1},1}(undef, length(plan.grid.grid))
    for i in eachindex(plan.grid.grid)
      nodes[i] = plan.grid.grid[i].points
    end

    if eltype(plan.grid.grid) <: ChebRoots
      return chebyshev_weights(y, Tuple(nodes), plan.order, plan.domain)
    elseif eltype(plan.grid.grid) <: ChebExtrema
      return chebyshev_weights_extrema(y, Tuple(nodes), plan.order, plan.domain)
    elseif eltype(plan.grid.grid) <: ChebExtended
      return chebyshev_weights_extended(y, Tuple(nodes), plan.order, plan.domain)
    end

  elseif typeof(plan) <: CApproxPlanPoly

    polynomials = Array{Array{T,2},1}(undef, length(plan.polys))
    for i in eachindex(plan.polys)
      polynomials[i] = plan.polys[i].poly
    end

    if plan.polys[1].nodetype == ChebRoots{T,T}
      return chebyshev_weights(y, Tuple(polynomials), plan.order)
    elseif plan.polys[1].nodetype == ChebExtrema{T,T}
      return chebyshev_weights_extrema(y, Tuple(polynomials), plan.order)
    elseif plan.polys[1].nodetype == ChebExtended{T,T}
      return chebyshev_weights_extended(y, Tuple(polynomials), plan.order)
    end

  end

end

"""
Computes the Chebyshev weights using multi-threading given the approximation sample, 'y', and the approximation plan, 'plan'.

y --- a N-dimensional array of data evaluated at the approximation points.

plan --- a structure specifying how the approximation should be undertaken.
""" 
function chebyshev_weights_threaded(y::AbstractArray{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  if typeof(plan) <: CApproxPlan

    nodes = Array{Array{T,1},1}(undef, length(plan.grid.grid))
    for i in eachindex(plan.grid.grid)
      nodes[i] = plan.grid.grid[i].points
    end

    if eltype(plan.grid.grid) <: ChebRoots
      return chebyshev_weights_threaded(y, Tuple(nodes), plan.order, plan.domain)
    elseif eltype(plan.grid.grid) <: ChebExtrema
      return chebyshev_weights_extrema_threaded(y, Tuple(nodes), plan.order, plan.domain)
    elseif eltype(plan.grid.grid) <: ChebExtended
      return chebyshev_weights_extended_threaded(y, Tuple(nodes), plan.order, plan.domain)
    end

  elseif typeof(plan) <: CApproxPlanPoly

    polynomials = Array{Array{T,2},1}(undef, length(plan.polys))
    for i in eachindex(plan.polys)
      polynomials[i] = plan.polys[i].poly
    end

    if plan.polys[1].nodetype == ChebRoots{T,T}
      return chebyshev_weights_threaded(y, Tuple(polynomials), plan.order)
    elseif plan.polys[1].nodetype == ChebExtrema{T,T}
      return chebyshev_weights_extrema_threaded(y, Tuple(polynomials), plan.order)
    elseif plan.polys[1].nodetype == ChebExtended{T,T}
      return chebyshev_weights_extended_threaded(y, Tuple(polynomials), plan.order)
    end

  end

end

"""
Computes the Chebyshev weights in a tensor-product polynomial given the data sample, 'y', the Chebyshev roots, 'nodes', the order of the polynomial, 'order', and the domain for the 
sampling points, 'domain'.

f --- an array containing the function evaluated on the approximation grid.

nodes --- a tuple of vectors containing the Chebyshev nodes used for approximation along each spacial dimension.

order --- a tuple of a vector specifying the polynomial's order for each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_weights(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
  end

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(y)

      product = one(T)
      @inbounds for j = 1:N
        product *= poly[j][s[j], i[j]]
      end

      numerator += y[s] * product
      denominator += product^2

    end

    weights[i] = numerator / denominator

  end

  return weights

end

"""
Computes the Chebyshev weights in a tensor-product polynomial given the data sample, 'y', the Chebyshev extrema, 'nodes', the order of the polynomial, 'order', and the domain for the 
sampling points, 'domain'.

f --- an array containing the function evaluated on the approximation grid.

nodes --- a tuple of vectors containing the Chebyshev extrema used for approximation along each spacial dimension.

order --- a tuple of a vector specifying the polynomial's order for each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_weights_extrema(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
  end

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(y)

      num = y[s]
      den = one(T)
      @inbounds for j = 1:N
        if s[j] === 1 || s[j] === n[j]
          scale = 0.5
        else
          scale = 1.0
        end
        temp = poly[j][s[j], i[j]]
        num *= temp * scale
        den *= (temp^2) * scale
      end

      numerator += num
      denominator += den

    end

    weights[i] = numerator / denominator

  end

  return weights

end

"""
Computes the Chebyshev weights in a tensor-product polynomial given the data sample, 'y', the extended Chebyshev points, 'nodes', the order of the polynomial, 'order', and the domain for the 
sampling points, 'domain'.

f --- an array containing the function evaluated on the approximation grid.

nodes --- a tuple of vectors containing the extended Chebyshev points used for approximation along each spacial dimension.

order --- a tuple of a vector specifying the polynomial's order for each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_weights_extended(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)
  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
    complement[i] = pinv(poly[i])'
  end

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(y)

      num = y[s]
      den = one(T)
      @inbounds for j = 1:N
        temp = complement[j][s[j], i[j]]
        num *= temp
        den *= temp * poly[j][s[j], i[j]]
      end

      numerator += num
      denominator += den

    end

    weights[i] = numerator / denominator

  end

  return weights

end

"""
Computes the Chebyshev weights in a tensor-product polynomial given the data sample, 'y', the Chebyshev polynomials evaluated at the Chebyshev roots, 'poly', and the order of the 
polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the Chebyshev roots.

order --- a tuple of a vector specifying the polynomial's order for each spacial dimension.
"""
function chebyshev_weights(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(y)

      product = one(T)
      @inbounds for j = 1:N
        product *= poly[j][s[j], i[j]]
      end

      numerator += y[s] * product
      denominator += product^2

    end

    weights[i] = numerator / denominator

  end

  return weights

end

"""
Computes the Chebyshev weights in a tensor-product polynomial given the data sample, 'y', the Chebyshev polynomials evaluated at the Chebyshev extrema, 'poly', and the order of the 
polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the Chebyshev extrema.

order --- a tuple of a vector specifying the polynomial's order for each spacial dimension.
"""
function chebyshev_weights_extrema(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(y)

      num = y[s]
      den = one(T)
      @inbounds for j = 1:N
        if s[j] == 1 || s[j] == n[j]
          scale = 0.5
        else
          scale = 1.0
        end
        temp = poly[j][s[j], i[j]]
        num *= temp * scale
        den *= (temp^2) * scale
      end

      numerator += num
      denominator += den

    end

    weights[i] = numerator / denominator

  end

  return weights

end

"""
Computes the Chebyshev weights in a tensor-product polynomial given the data sample, 'y', the Chebyshev polynomials evaluated at the extended Chebyshev points, 'poly', and the order of the 
polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the extended Chebyshev points.

order --- a tuple of a vector specifying the polynomial's order for each spacial dimension.
"""
function chebyshev_weights_extended(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(y)

      num = y[s]
      den = one(T)
      @inbounds for j = 1:N
        temp = complement[j][s[j], i[j]]
        num *= temp
        den *= temp * poly[j][s[j], i[j]]
      end

      numerator += num
      denominator += den

    end

    weights[i] = numerator / denominator

  end

  return weights

end

"""
Computes the Chebyshev weights in a complete polynomial given the data sample, 'y', the Chebyshev roots, 'nodes', the order of the polynomial, 'order', and the domain for the 
sampling points, 'domain'.

f --- an array containing the function evaluated on the approximation grid.

nodes --- a tuple of vectors containing the Chebyshev nodes used for approximation along each spacial dimension.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_weights(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
  end

  ord = Tuple([order for _ in 1:N])

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(y)

        product = one(T)
        @inbounds for j = 1:N
          product *= poly[j][s[j], i[j]]
        end

        numerator += y[s] * product
        denominator += product^2

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

"""
Computes the Chebyshev weights in a complete polynomial given the data sample, 'y', the Chebyshev extrema, 'nodes', the order of the polynomial, 'order', and the domain for the 
sampling points, 'domain'.

f --- an array containing the function evaluated on the approximation grid.

nodes --- a tuple of vectors containing the Chebyshev extrema used for approximation along each spacial dimension.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_weights_extrema(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
  end

  ord = Tuple([order for _ in 1:N])

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(y)

        num = y[s]
        den = one(T)
        @inbounds for j = 1:N
          if s[j] == 1 || s[j] == n[j]
            scale = 0.5
          else
            scale = 1.0
          end
          temp = poly[j][s[j], i[j]]
          num *= temp * scale
          den *= (temp^2) * scale
        end

        numerator += num
        denominator += den

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

"""
Computes the Chebyshev weights in a complete polynomial given the data sample, 'y', the extended Chebyshev points, 'nodes', the order of the polynomial, 'order', and the domain for the 
sampling points, 'domain'.

f --- an array containing the function evaluated on the approximation grid.

nodes --- a tuple of vectors containing the extended Chebyshev points used for approximation along each spacial dimension.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_weights_extended(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)
  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
    complement[i] = pinv(poly[i])'
  end

  ord = Tuple([order for _ in 1:N])

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(y)

        num = y[s]
        den = one(T)
        @inbounds for j = 1:N
          temp = complement[j][s[j], i[j]]
          num *= temp
          den *= temp * poly[j][s[j], i[j]]
        end

        numerator += num
        denominator += den

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

"""
Computes the Chebyshev weights in a complete polynomial given the data sample, 'y', the Chebyshev polynomials evaluated at the Chebyshev roots, 'poly', and the order of the 
polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the Chebyshev roots.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.
"""
function chebyshev_weights(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple([order for _ in 1:N])

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(y)

        product = one(T)
        @inbounds for j = 1:N
          product *= poly[j][s[j], i[j]]
        end

        numerator += y[s] * product
        denominator += product^2

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

"""
Computes the Chebyshev weights in a complete polynomial given the data sample, 'y', the Chebyshev polynomials evaluated at the Chebyshev extrema, 'poly', and the order of the 
polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the Chebyshev extrema.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.
"""
function chebyshev_weights_extrema(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  ord = Tuple([order for _ in 1:N])

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(y)

        num = y[s]
        den = one(T)
        @inbounds for j = 1:N
          if s[j] == 1 || s[j] == n[j]
            scale = 0.5
          else
            scale = 1.0
          end
          temp = poly[j][s[j], i[j]]
          num *= temp * scale
          den *= (temp^2) * scale
        end

        numerator += num
        denominator += den

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

"""
Computes the Chebyshev weights in a complete polynomial given the data sample, 'y', the Chebyshev polynomials evaluated at the extended Chebyshev points, 'poly', and the order of the 
polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the extended Chebyshev points.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.
"""
function chebyshev_weights_extended(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  complement = Array{Array{T,2},1}(undef, N)
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end

  ord = Tuple([order for _ in 1:N])

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(y)

        num = y[s]
        den = one(T)
        @inbounds for j = 1:N
          temp = complement[j][s[j], i[j]]
          num *= temp
          den *= temp * poly[j][s[j], i[j]]
        end

        numerator += num
        denominator += den

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

"""
Computes the Chebyshev weights using multi-threading in a complete polynomial given the data sample, 'y', the Chebyshev roots, 'nodes', the order of the polynomial, 'order', and the domain for the 
sampling points, 'domain'.

f --- an array containing the function evaluated on the approximation grid.

nodes --- a tuple of vectors containing the Chebyshev roots used for approximation along each spacial dimension.

order --- a vector of integers specifying the polynomial's order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_weights_threaded(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
  end

  weights = zeros(Tuple(order .+ 1))

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(y)

      product = one(T)
      @inbounds for j = 1:N
        product *= poly[j][s[j], i[j]]
      end

      numerator += y[s] * product
      denominator += product^2

    end

    weights[i] = numerator / denominator

  end

  return weights

end

"""
Computes the Chebyshev weights using multi-threading in a complete polynomial given the data sample, 'y', the Chebyshev extrema, 'nodes', the order of the polynomial, 'order', and the domain for the 
sampling points, 'domain'.

f --- an array containing the function evaluated on the approximation grid.

nodes --- a tuple of vectors containing the Chebyshev extrema used for approximation along each spacial dimension.

order --- a vector of integers specifying the polynomial's order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_weights_extrema_threaded(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
  end

  weights = zeros(Tuple(order .+ 1))

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(y)

      num = y[s]
      den = one(T)
      @inbounds for j = 1:N
        if s[j] == 1 || s[j] == n[j]
          scale = 0.5
        else
          scale = 1.0
        end
        temp = poly[j][s[j], i[j]]
        num *= temp * scale
        den *= (temp^2) * scale
      end

      numerator += num
      denominator += den

    end

    weights[i] = numerator / denominator

  end

  return weights

end

"""
Computes the Chebyshev weights using multi-threading in a complete polynomial given the data sample, 'y', the extended Chebyshev points, 'nodes', the order of the polynomial, 'order', and the domain for the 
sampling points, 'domain'.

f --- an array containing the function evaluated on the approximation grid.

nodes --- a tuple of vectors containing the extended Chebyshev points used for approximation along each spacial dimension.

order --- a vector of integers specifying the polynomial's order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_weights_extended_threaded(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)
  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
    complement[i] = pinv(poly[i])'
  end

  weights = zeros(Tuple(order .+ 1))

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(y)

      num = y[s]
      den = one(T)
      @inbounds for j = 1:N
        temp = complement[j][s[j], i[j]]
        num *= temp
        den *= temp * poly[j][s[j], i[j]]
      end

      numerator += num
      denominator += den

    end

    weights[i] = numerator / denominator

  end

  return weights

end

"""
Computes the Chebyshev weights using multi-threading in a tensor-product polynomial given the data sample, 'y', the Chebyshev polynomials 
evaluated at the Chebyshev roots, 'poly', and the order of the polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the Chebyshev roots.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.
"""
function chebyshev_weights_threaded(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(y)

      product = one(T)
      @inbounds for j = 1:N
        product *= poly[j][s[j], i[j]]
      end

      numerator += y[s] * product
      denominator += product^2

    end

    weights[i] = numerator / denominator

  end

  return weights

end

"""
Computes the Chebyshev weights using multi-threading in a tensor-product polynomial given the data sample, 'y', the Chebyshev polynomials 
evaluated at the Chebyshev extrema, 'poly', and the order of the polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the Chebyshev extrema.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.
"""
function chebyshev_weights_extrema_threaded(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(y)

      num = y[s]
      den = one(T)
      @inbounds for j = 1:N
        if s[j] == 1 || s[j] == n[j]
          scale = 0.5
        else
          scale = 1.0
        end
        temp = poly[j][s[j], i[j]]
        num *= temp * scale
        den *= (temp^2) * scale
      end

      numerator += num
      denominator += den

    end

    weights[i] = numerator / denominator

  end

  return weights

end

"""
Computes the Chebyshev weights using multi-threading in a tensor-product polynomial given the data sample, 'y', the Chebyshev polynomials 
evaluated at the extended Chebyshev points, 'poly', and the order of the polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the Chebyshev roots.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.
"""
function chebyshev_weights_extended_threaded(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  complement = Array{Array{T,2},1}(undef, N)
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(y)

      num = y[s]
      den = one(T)
      @inbounds for j = 1:N
        temp = complement[j][s[j], i[j]]
        num *= temp
        den *= temp * poly[j][s[j], i[j]]
      end

      numerator += num
      denominator += den

    end

    weights[i] = numerator / denominator

  end

  return weights

end

"""
Computes the Chebyshev weights using multi-threading in a complete polynomial given the data sample, 'y', the Chebyshev polynomials 
evaluated at the Chebyshev roots, 'poly', and the order of the polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the Chebyshev roots.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.
"""
function chebyshev_weights_threaded(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
  end

  ord = Tuple([order for _ in 1:N])

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(y)

        product = one(T)
        @inbounds for j = 1:N
          product *= poly[j][s[j], i[j]]
        end

        numerator += y[s] * product
        denominator += product^2

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

"""
Computes the Chebyshev weights using multi-threading in a complete polynomial given the data sample, 'y', the Chebyshev polynomials 
evaluated at the Chebyshev extrema, 'poly', and the order of the polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the Chebyshev extrema.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.
"""
function chebyshev_weights_extrema_threaded(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
  end

  ord = Tuple([order for _ in 1:N])

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(y)

        num = y[s]
        den = one(T)
        @inbounds for j = 1:N
          if s[j] == 1 || s[j] == n[j]
            scale = 0.5
          else
            scale = 1.0
          end
          temp = poly[j][s[j], i[j]]
          num *= temp * scale
          den *= (temp^2) * scale
        end

        numerator += num
        denominator += den

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

"""
Computes the Chebyshev weights using multi-threading in a complete polynomial given the data sample, 'y', the Chebyshev polynomials 
evaluated at the extended Chebyshev points, 'poly', and the order of the polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the extended Chebyshev points.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.
"""
function chebyshev_weights_extended_threaded(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)
  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
    complement[i] = pinv(poly[i])'
  end

  ord = Tuple([order for _ in 1:N])

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(y)

        num = y[s]
        den = one(T)
        @inbounds for j = 1:N
          temp = complement[j][s[j], i[j]]
          num *= temp
          den *= temp * poly[j][s[j], i[j]]
        end

        numerator += num
        denominator += den

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

"""
Computes the Chebyshev weights using multi-threading in a complete polynomial given the data sample, 'y', the Chebyshev polynomials 
evaluated at the Chebyshev roots, 'poly', and the order of the polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the Chebyshev roots.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.
"""
function chebyshev_weights_threaded(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple([order for _ in 1:N])

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(y)

        product = one(T)
        @inbounds for j = 1:N
          product *= poly[j][s[j], i[j]]
        end

        numerator += y[s] * product
        denominator += product^2

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

"""
Computes the Chebyshev weights using multi-threading in a complete polynomial given the data sample, 'y', the Chebyshev polynomials 
evaluated at the Chebyshev extrema, 'poly', and the order of the polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the Chebyshev extrema.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.
"""
function chebyshev_weights_extrema_threaded(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  ord = Tuple([order for _ in 1:N])

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(y)

        num = y[s]
        den = one(T)
        @inbounds for j = 1:N
          if s[j] == 1 || s[j] == n[j]
            scale = 0.5
          else
            scale = 1.0
          end
          temp = poly[j][s[j], i[j]]
          num *= temp * scale
          den *= (temp^2) * scale
        end

        numerator += num
        denominator += den

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

"""
Computes the Chebyshev weights using multi-threading in a complete polynomial given the data sample, 'y', the Chebyshev polynomials 
evaluated at the extended Chebyshev points, 'poly', and the order of the polynomial, 'order'.

f --- an array containing the function evaluated on the approximation grid.

poly --- a tuple of matrices containing the Chebyshev polynomials evaluated at the extended Chebyshev points.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.
"""
function chebyshev_weights_extended_threaded(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  complement = Array{Array{T,2},1}(undef, N)
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end

  ord = Tuple([order for _ in 1:N])

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(y)

        num = y[s]
        den = one(T)
        @inbounds for j = 1:N
          temp = complement[j][s[j], i[j]]
          num *= temp
          den *= temp * poly[j][s[j], i[j]]
        end

        numerator += num
        denominator += den

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

# Functions for the one-variable case where the nodes are a vector

chebyshev_weights(y::Array{T,1}, nodes::Array{T,1}, order::Union{S,Array{S,1}}, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights(y, (nodes,), order, domain)
chebyshev_weights_extrema(y::Array{T,1}, nodes::Array{T,1}, order::Union{S,Array{S,1}}, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extrema(y, (nodes,), order, domain)
chebyshev_weights_extended(y::Array{T,1}, nodes::Array{T,1}, order::Union{S,Array{S,1}}, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extended(y, (nodes,), order, domain)
chebyshev_weights(y::Array{T,1}, poly::Array{T,2}, order::Union{S,Array{S,1}}) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights(y, (poly,), order)
chebyshev_weights_extrema(y::Array{T,1}, poly::Array{T,2}, order::Union{S,Array{S,1}}) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extrema(y, (poly,), order)
chebyshev_weights_extended(y::Array{T,1}, poly::Array{T,2}, order::Union{S,Array{S,1}}) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extended(y, (poly,), order)

# Functions that allow the nodes to be in an array of arrays

# Serial functions

chebyshev_weights(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights(y, Tuple(nodes), order, domain)
chebyshev_weights_extrema(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema(y, Tuple(nodes), order, domain)
chebyshev_weights_extended(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended(y, Tuple(nodes), order, domain)
chebyshev_weights(y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights(y, Tuple(poly), order)
chebyshev_weights_extrema(y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema(y, Tuple(poly), order)
chebyshev_weights_extended(y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended(y, Tuple(poly), order)
chebyshev_weights(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights(y, Tuple(nodes), order, domain)
chebyshev_weights_extrema(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema(y, Tuple(nodes), order, domain)
chebyshev_weights_extended(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended(y, Tuple(nodes), order, domain)
chebyshev_weights(y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights(y, Tuple(poly), order)
chebyshev_weights_extrema(y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema(y, Tuple(poly), order)
chebyshev_weights_extended(y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended(y, Tuple(poly), order)

# Threaded functions

chebyshev_weights_threaded(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded(y, Tuple(nodes), order, domain)
chebyshev_weights_extrema_threaded(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded(y, Tuple(nodes), order, domain)
chebyshev_weights_extended_threaded(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded(y, Tuple(nodes), order, domain)
chebyshev_weights_threaded(y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded(y, Tuple(poly), order)
chebyshev_weights_extrema_threaded(y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded(y, Tuple(poly), order)
chebyshev_weights_extended_threaded(y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded(y, Tuple(poly), order)
chebyshev_weights_threaded(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded(y, Tuple(nodes), order, domain)
chebyshev_weights_extrema_threaded(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded(y, Tuple(nodes), order, domain)
chebyshev_weights_extended_threaded(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded(y, Tuple(nodes), order, domain)
chebyshev_weights_threaded(y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded(y, Tuple(poly), order)
chebyshev_weights_extrema_threaded(y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded(y, Tuple(poly), order)
chebyshev_weights_extended_threaded(y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded(y, Tuple(poly), order)

# Functions to evaluate Chebyshev polynominals

"""
Evaluate a tensor-product Chebyshev polynomial at point, 'x', given the polynomial weights, 'weights', the order of the polynomial, 'order', and the domain, 'domain'.

weights --- an N-dimensional array containing the polynomial weights.

x --- a vector of Numbers representing the evaluation point.

order --- a tuple or vector specifying the polynomial's order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_evaluate(weights::AbstractArray{T,N}, x::AbstractArray{R,1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  if length(x) != N
    error("A value for 'x' is needed for each spacial dimension.")
  end

  poly = Array{Array{R,2},1}(undef, N)
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(x[i], domain[:, i]))
  end

  yhat = zero(T)
  @inbounds for i in CartesianIndices(weights)
    poly_product = poly[1][i[1]]
    @inbounds for j = 2:N
      poly_product *= poly[j][i[j]]
    end
    yhat += weights[i] * poly_product
  end

  return yhat

end

"""
Evaluate a complete Chebyshev polynomial at point, 'x', given the polynomial weights, 'weights', the order of the polynomial, 'order', and the domain, 'domain'.

weights --- an N-dimensional array containing the polynomial weights.

x --- a vector of Numbers representing the evaluation point.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_evaluate(weights::AbstractArray{T,N}, x::AbstractArray{R,1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  if length(x) != N
    error("A value for 'x' is needed for each spacial dimension.")
  end

  poly = Array{Array{R,2},1}(undef, N)
  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(x[i], domain[:, i]))
  end

  yhat = zero(T)
  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N
      poly_product = poly[1][i[1]]
      @inbounds for j = 2:N
        poly_product *= poly[j][i[j]]
      end
      yhat += weights[i] * poly_product
    end
  end

  return yhat

end

"""
Creates an interpolating function given an approximation plan, 'plan', that can evaluate the Chebyshev polynomial at any point in the state-space.
  
y --- an N-dimensional array containing the function evaluated at each point of the approximation grid.

plan --- a structure specifying how the approximation should be undertaken.
"""
function chebyshev_interp(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights(y, plan)

  function interp(x::AbstractArray{R,1}) where {R<:Number}

    yhat = chebyshev_evaluate(w, x, plan.order, plan.domain)

    return yhat

  end

  return interp

end

"""
Creates an interpolating function given an approximation plan, 'plan', that can evaluate using multi-threading the Chebyshev polynomial at any point in the state-space.
  
y --- an N-dimensional array containing the function evaluated at each point of the approximation grid.

plan --- a structure specifying how the approximation should be undertaken.
"""
function chebyshev_interp_threaded(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights_threaded(y, plan)

  function interp(x::AbstractArray{R,1}) where {R<:Number}

    yhat = chebyshev_evaluate(w, x, plan.order, plan.domain)

    return yhat

  end

  return interp

end

# Functions for derivatives

"""
Computes the first derivative of a tensor-product Chebyshev polynomial with respect to the i'th variable, 'pos', evaluated at 'x', given the 
polynomial weights, 'weights', the polynomial's order, 'order', and the domain, 'domain'.

weights --- an N-dimensional array containinbg the polynomial weights.

x --- a vector containing the point at which to evaluate the derivative.

pos --- an integer specifying the variable that the function is being differentiated with respect to.

order --- a tuple or vector specifying the polynomial's order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_derivative(weights::Array{T,N}, x::AbstractArray{R,1}, pos::S, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  if length(x) != N
    error("A value for 'x' is needed for each spacial dimension.")
  end

  poly = Array{Array{R,2},1}(undef, N)
  @inbounds for i = 1:N
    if i === pos
      poly[i] = chebyshev_polynomial_deriv(order[i], normalize_node(x[i], domain[:, i]))
    else
      poly[i] = chebyshev_polynomial(order[i], normalize_node(x[i], domain[:, i]))
    end
  end

  derivative = zero(T)
  @inbounds for i in CartesianIndices(weights)
    poly_product = poly[1][i[1]]
    @inbounds for j = 2:N
      poly_product *= poly[j][i[j]]
    end
    derivative += weights[i] * poly_product
  end

  return derivative * (2.0 / (domain[1, pos] - domain[2, pos]))

end

"""
Computes the first derivative of a complete Chebyshev polynomial with respect to the i'th variable, 'pos', evaluated at 'x', given the 
polynomial weights, 'weights', the polynomial's order, 'order', and the domain, 'domain'.

weights --- an N-dimensional array containinbg the polynomial weights.

x --- a vector containing the point at which to evaluate the derivative.

pos --- an integer specifying the variable that the function is being differentiated with respect to.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_derivative(weights::Array{T,N}, x::AbstractArray{R,1}, pos::S, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  if length(x) != N
    error("A value for 'x' is needed for each spacial dimension.")
  end

  poly = Array{Array{R,2},1}(undef, N)
  @inbounds for i = 1:N
    if i === pos
      poly[i] = chebyshev_polynomial_deriv(order, normalize_node(x[i], domain[:, i]))
    else
      poly[i] = chebyshev_polynomial(order, normalize_node(x[i], domain[:, i]))
    end
  end

  derivative = zero(T)
  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N
      poly_product = poly[1][i[1]]
      @inbounds for j = 2:N
        poly_product *= poly[j][i[j]]
      end
      derivative += weights[i] * poly_product
    end
  end

  return derivative * (2.0 / (domain[1, pos] - domain[2, pos]))

end

# Functions for gradients

"""
Computes the gradient of a tensor-product Chebyshev polynomial evaluated at 'x', given the 
polynomial weights, 'weights', the polynomial's order, 'order', and the domain, 'domain'.

weights --- an N-dimensional array containinbg the polynomial weights.

x --- a vector containing the point at which to evaluate the derivative.

order --- a tuple or vector specifying the polynomial's order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_gradient(weights::Array{T,N}, x::AbstractArray{R,1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  if length(x) != N
    error("A value for 'x' is needed for each spacial dimension.")
  end

  gradient = Array{R,2}(undef, 1, N)

  @inbounds for i = 1:N
    gradient[i] = chebyshev_derivative(weights, x, i, order, domain)
  end

  return gradient

end

"""
Computes the gradient of a complete Chebyshev polynomial evaluated at 'x', given the 
polynomial weights, 'weights', the polynomial's order, 'order', and the domain, 'domain'.

weights --- an N-dimensional array containinbg the polynomial weights.

x --- a vector containing the point at which to evaluate the derivative.

order --- an integer specifying the polynomial's maximal order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_gradient(weights::Array{T,N}, x::AbstractArray{R,1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  if length(x) != N
    error("A value for 'x' is needed for each spacial dimension.")
  end

  gradient = Array{R,2}(undef, 1, N)

  @inbounds for i = 1:N
    gradient[i] = chebyshev_derivative(weights, x, i, order, domain)
  end

  return gradient

end

"""
Creates a function that evaluates the gradient of a Chebyshev polynomial evaluated at 'x', given an approximation plan, 'plan'.

y --- an N-dimensional array containing function evaluated on the approximation grid.

plan --- a structure specifying how the approximation should be undertaken.
"""
function chebyshev_gradient(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights(y, plan)

  function cheb_grad(x::AbstractArray{R,1}) where {R<:Number}

    grad = chebyshev_gradient(w, x, plan.order, plan.domain)

    return grad

  end

  return cheb_grad

end

"""
Creates a function that evaluates using multi-threading the gradient of a Chebyshev polynomial evaluated at 'x', given an approximation plan, 'plan'.

y --- an N-dimensional array containing function evaluated on the approximation grid.

plan --- a structure specifying how the approximation should be undertaken.
"""
function chebyshev_gradient_threaded(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights_threaded(y, plan)

  function cheb_grad(x::AbstractArray{R,1}) where {R<:Number}

    grad = chebyshev_gradient(w, x, plan.order, plan.domain)

    return grad

  end

  return cheb_grad

end

# Functions for hessians

"""
Computes the hessian of a tensor-product Chebyshev polynomial evaluated at 'x', given the 
polynomial weights, 'weights', the polynomial's order, 'order', and the domain, 'domain'.

weights --- an N-dimensional array containinbg the polynomial weights.

x --- a vector containing the point at which to evaluate the derivative.

order --- a tuple or vector specifying the polynomial's order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_hessian(weights::Array{T,N}, x::AbstractArray{R,1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  if length(x) != N
    error("A value for 'x' is needed for each spacial dimension.")
  end

  hess = Array{T,2}(undef, N, N)
  poly = Array{Array{R,2},1}(undef, N)

  @inbounds for j in CartesianIndices(hess)
    @inbounds for i = 1:N
      if i == j[1] == j[2]
        poly[i] = chebyshev_polynomial_sec_deriv(order[i], normalize_node(x[i], domain[:, i]))
      elseif i == j[1] || i == j[2]
        poly[i] = chebyshev_polynomial_deriv(order[i], normalize_node(x[i], domain[:, i]))
      else
        poly[i] = chebyshev_polynomial(order[i], normalize_node(x[i], domain[:, i]))
      end
    end

    deriv = zero(T)
    @inbounds for i in CartesianIndices(weights)
      poly_product = poly[1][i[1]]
      @inbounds for j = 2:N
        poly_product *= poly[j][i[j]]
      end
      deriv += weights[i] * poly_product
    end

    hess[j] = deriv * (2.0 / (domain[1, j[1]] - domain[2, j[1]])) * (2.0 / (domain[1, j[2]] - domain[2, j[2]]))

  end

  return hess

end

"""
Computes the hessian of a complete Chebyshev polynomial evaluated at 'x', given the 
polynomial weights, 'weights', the polynomial's order, 'order', and the domain, 'domain'.

weights --- an N-dimensional array containinbg the polynomial weights.

x --- a vector containing the point at which to evaluate the derivative.

order --- a tuple or vector specifying the polynomial's order along each spacial dimension.

domain --- a matrix containing the upper and lower bounds on the domain for each spacial dimension.
"""
function chebyshev_hessian(weights::Array{T,N}, x::AbstractArray{R,1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  if length(x) != N
    error("A value for 'x' is needed for each spacial dimension.")
  end

  hess = Array{T,2}(undef, N, N)
  poly = Array{Array{R,2},1}(undef, N)

  @inbounds for j in CartesianIndices(hess)
    @inbounds for i = 1:N
      if i == j[1] == j[2]
        poly[i] = chebyshev_polynomial_sec_deriv(order, normalize_node(x[i], domain[:, i]))
      elseif i == j[1] || i == j[2]
        poly[i] = chebyshev_polynomial_deriv(order, normalize_node(x[i], domain[:, i]))
      else
        poly[i] = chebyshev_polynomial(order, normalize_node(x[i], domain[:, i]))
      end
    end

    deriv = zero(T)
    @inbounds for i in CartesianIndices(weights)
      poly_product = poly[1][i[1]]
      @inbounds for j = 2:N
        poly_product *= poly[j][i[j]]
      end
      deriv += weights[i] * poly_product
    end

    hess[j] = deriv * (2.0 / (domain[1, j[1]] - domain[2, j[1]])) * (2.0 / (domain[1, j[2]] - domain[2, j[2]]))

  end

  return hess

end

"""
Creates a function that evaluates the hessian of a Chebyshev polynomial evaluated at 'x', given an approximation plan, 'plan'.

y --- an N-dimensional array containing function evaluated on the approximation grid.

plan --- a structure specifying how the approximation should be undertaken.
"""
function chebyshev_hessian(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights(y, plan)

  function cheb_hess(x::AbstractArray{R,1}) where {R<:Number}

    hess = chebyshev_hessian(w, x, plan.order, plan.domain)

    return hess

  end

  return cheb_hess

end

"""
Creates a function that uses multi-threading to evaluate the hessian of a Chebyshev polynomial evaluated at 'x', given an approximation plan, 'plan'.

y --- an N-dimensional array containing function evaluated on the approximation grid.

plan --- a structure specifying how the approximation should be undertaken.
"""
function chebyshev_hessian_threaded(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights_threaded(y, plan)

  function cheb_hess(x::AbstractArray{R,1}) where {R<:Number}

    hess = chebyshev_hessian(w, x, plan.order, plan.domain)

    return hess

  end

  return cheb_hess

end