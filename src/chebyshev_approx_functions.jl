abstract type Nodes end
abstract type CApproximationPlan end

struct ChebRoots{T<:AbstractFloat} <: Nodes

  points::Array{T,1}
  domain::Array{T,1}

end

struct ChebExtrema{T<:AbstractFloat} <: Nodes

  points::Array{T,1}
  domain::Array{T,1}

end

struct ChebExtended{T<:AbstractFloat} <: Nodes

  points::Array{T,1}
  domain::Array{T,1}

end

struct VertesiNodes{T<:AbstractFloat} <: Nodes

  points::Array{T,1}
  domain::Array{T,1}

end

struct LegendreNodes{T<:AbstractFloat} <: Nodes

  points::Array{T,1}
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

struct CApproxPlanPoly{N,T<:AbstractFloat,S<:Integer,G<:ChebPoly} <: CApproximationPlan

  polys::NTuple{N,G}
  order::Union{S,Tuple{Vararg{S}}}
  domain::Union{Array{T,1},Array{T,2}}

end

function chebyshev_nodes(N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = Array{Float64,1}(undef, N)

  if isodd(N)
    points[div(N - 1, 2)+1] = (domain[1] + domain[2]) * 0.5
  end

  @inbounds for i = 1:div(N, 2)
    x = -cos((i - 0.5) * π / N) * (domain[1] - domain[2]) * 0.5
    points[i] = (domain[1] + domain[2]) * 0.5 + x
    points[end-i+1] = (domain[1] + domain[2]) * 0.5 - x
  end

  return points

end

function chebyshev_extrema(N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = Array{Float64,1}(undef, N)

  if isodd(N)
    points[div(N - 1, 2)+1] = (domain[1] + domain[2]) * 0.5
  end

  @inbounds for i = 1:div(N, 2)
    x = -cos((i - 1) * π / (N - 1)) * (domain[1] - domain[2]) * 0.5
    points[i] = (domain[1] + domain[2]) * 0.5 + x
    points[end-i+1] = (domain[1] + domain[2]) * 0.5 - x
  end

  return points

end

function chebyshev_extended(N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = Array{Float64,1}(undef, N)

  if isodd(N)
    points[div(N - 1, 2)+1] = 0.0
  end

  @inbounds for i = 1:div(N, 2)
    x = -cos((i - 0.5) * π / N)
    points[i] = x
    points[end-i+1] = -x
  end

  points .= (domain[1] + domain[2]) * 0.5 .+ points * ((domain[1] - domain[2]) * 0.5) / cos(π / (2N))

  return points

end

function vertesi_nodes(N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = Array{Float64,1}(undef, N)

  if isodd(N)
    points[div(N - 1, 2)+1] = (domain[1] + domain[2]) * 0.5
  end

  if N == 1

    return points

  elseif N == 2

    points[begin] = domain[2]
    points[end] = domain[1]

    return points

  else

    points[begin] = domain[2]
    points[end] = domain[1]

    @inbounds for i = 2:div(N, 2)
      x = -cos((π * 0.5) * (2i - 1) / N) / cos((π / (2N)) * (1 + 1 / (4 * log(N)))) * (domain[1] - domain[2]) * 0.5
      points[i] = (domain[1] + domain[2]) * 0.5 + x
      points[end-i+1] = (domain[1] + domain[2]) * 0.5 - x
    end

    return points

  end

end

function legendre_nodes(N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = Array{Float64,1}(undef, N)
  update = Array{Float64,1}(undef, N)

  p = zeros(N, N)

  if isodd(N)
    points[div(N - 1, 2)+1] = 0.0
  end

  if N == 1
    point[1] = (domain[1] + domain[2]) * 0.5
    return points
  else
    points[begin] = -1.0
    points[end] = 1.0

    @inbounds for i = 2:div(N, 2)
      points[i] = -cos(π * (i - 1) / (N - 1))
      points[end-i+1] = cos(π * (i - 1) / (N - 1))
    end

    p[:, 1] .= 1.0

    len = Inf
    @views while len > eps(1.0)
      p[:, 2] .= points
      @inbounds for i = 2:(N-1)
        p[:, i+1] .= ((2i - 1) * points .* p[:, i] .- (i - 1) * p[:, i-1]) / i
      end
      update .= (points .* p[:, end] .- p[:, end-1]) ./ (N * p[:, end])
      points .-= update
      len = maximum(abs, update)
    end

    points .= (domain[1] + domain[2]) * 0.5 .+ points * ((domain[1] - domain[2]) * 0.5)

    return points
  end

end

function nodes(N::S, node_generator::Symbol, domain=[1.0, -1.0],) where {S<:Integer}

  if node_generator == :chebyshev_nodes
    return ChebRoots(chebyshev_nodes(N, domain), domain)
  elseif node_generator == :chebyshev_extrema
    return ChebExtrema(chebyshev_extrema(N, domain), domain)
  elseif node_generator == :chebyshev_extended
    return ChebExtended(chebyshev_extended(N, domain), domain)
  elseif node_generator == :vertesi_nodes
    return VertesiNodes(vertesi_nodes(N, domain), domain)
  elseif node_generator == :legendre_nodes
    return LegendreNodes(legendre_nodes(N, domain), domain)
  end
end

function normalize_node(node::R, domain::Array{T,1}) where {R<:Number,T<:AbstractFloat}

  if domain[1] == domain[2]
    norm_node = zero(T)
    return norm_node
  else
    norm_node = 2 * (node - domain[2]) / (domain[1] - domain[2]) - 1
    return norm_node
  end

end

function normalize_node(node::Array{R,1}, domain::Array{T,1}) where {R<:Number,T<:AbstractFloat}

  norm_nodes = map(x -> normalize_node(x, domain), node)
  return norm_nodes

end

function normalize_node(node::G) where {G<:Nodes}

  norm_nodes = map(x -> normalize_node(x, node.domain), node.points)
  return norm_nodes

end

function chebyshev_polynomial(order::S, x::R) where {R<:Number,S<:Integer}

  poly = Array{R}(undef, 1, order + 1)
  poly[1] = one(R)

  @inbounds for i = 2:order+1
    if i == 2
      poly[i] = x
    else
      poly[i] = 2 * x * poly[i-1] - poly[i-2]
    end
  end

  return poly

end

function chebyshev_polynomial(order::S, x::AbstractArray{R,1}) where {R<:Number,S<:Integer}

  poly = Array{R}(undef, length(x), order + 1)
  poly[:, 1] .= ones(R, length(x))

  @inbounds for i = 2:order+1
    for j in eachindex(x)
      if i == 2
        poly[j, i] = x[j]
      else
        poly[j, i] = 2 * x[j] * poly[j, i-1] - poly[j, i-2]
      end
    end
  end

  return poly

end

function chebyshev_polynomial(order::S, g::G) where {G<:Nodes,S<:Integer}

  T = eltype(g.points)

  poly = Array{T}(undef, length(g.points), order + 1)
  poly[:, 1] .= ones(T, length(g.points))

  @inbounds for i = 2:order+1
    for j in eachindex(g.points)
      if i == 2
        poly[j, i] = normalize_node(g.points[j], g.domain)
      else
        poly[j, i] = 2 * normalize_node(g.points[j], g.domain) * poly[j, i-1] - poly[j, i-2]
      end
    end
  end

  return ChebPoly(poly, typeof(g))

end

function chebyshev_polynomial_deriv(order::S, x::T) where {T<:Number,S<:Integer}

  poly_deriv = Array{T}(undef, 1, order + 1)
  poly_deriv[1] = zero(T)
  p = one(T)
  pl = NaN
  pll = NaN

  @inbounds for i = 2:order+1
    if i == 2
      pl, p = p, x
      poly_deriv[i] = one(T)
    else
      pll, pl = pl, p
      p = 2 * x * pl - pll
      poly_deriv[i] = 2 * pl + 2 * x * poly_deriv[i-1] - poly_deriv[i-2]
    end
  end

  return poly_deriv

end

function chebyshev_polynomial_deriv(order::S, x::AbstractArray{T,1}) where {S<:Integer,T<:Number}

  poly_deriv = Array{T}(undef, order + 1, length(x))
  poly_deriv[1, :] .= zeros(T, length(x))

  @inbounds for j in eachindex(x)
    p = one(T)
    pl = NaN
    pll = NaN
    for i = 2:order+1
      if i == 2
        pl, p = p, x[j]
        poly_deriv[i, j] = one(T)
      else
        pll, pl = pl, p
        p = 2 * x[j] * pl - pll
        poly_deriv[i, j] = 2 * pl + 2 * x[j] * poly_deriv[i-1, j] - poly_deriv[i-2, j]
      end
    end
  end

  return transpose(poly_deriv)# <: AbstractArray

end

function chebyshev_polynomial_deriv(order::S, g::G) where {G<:Nodes,S<:Integer}

  T = eltype(g.points)

  poly_deriv = Array{T}(undef, order + 1, length(g.points))
  poly_deriv[1, :] .= zeros(T, length(g.points))

  @inbounds for j in eachindex(g.points)
    p = one(T)
    pl = NaN
    pll = NaN
    for i = 2:order+1
      if i == 2
        pl, p = p, normalize_node(g.points[j], g.domain)
        poly_deriv[i, j] = one(T)
      else
        pll, pl = pl, p
        p = 2 * normalize_node(g.points[j], g.domain) * pl - pll
        poly_deriv[i, j] = 2 * pl + 2 * normalize_node(g.points[j], g.domain) * poly_deriv[i-1, j] - poly_deriv[i-2, j]
      end
    end
  end

  return ChebPoly(transpose(poly_deriv), typeof(g))

end

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
      p = 2 * x * pl - pll
      pdll, pdl = pdl, pd
      pd = 2 * pl + 2 * x * pdl - pdll
      poly_sec_deriv[i] = 2 * x * poly_sec_deriv[i-1] + 4 * pdl - poly_sec_deriv[i-2]
    end
  end

  return poly_sec_deriv

end

function chebyshev_polynomial_sec_deriv(order::S, x::AbstractArray{T,1}) where {S<:Integer,T<:Number}

  poly_sec_deriv = Array{T}(undef, order + 1, length(x))
  poly_sec_deriv[1, :] .= zeros(T, length(x))

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
        poly_sec_deriv[i, j] = zero(T)
      else
        pll, pl = pl, p
        p = 2 * x[j] * pl - pll
        pdll, pdl = pdl, pd
        pd = 2 * pl + 2 * x[j] * pdl - pdll
        poly_sec_deriv[i, j] = 2 * x[j] * poly_sec_deriv[i-1, j] + 4 * pdl - poly_sec_deriv[i-2, j]
      end
    end
  end

  return transpose(poly_sec_deriv)

end

function chebyshev_polynomial_sec_deriv(order::S, g::G) where {G<:Nodes,S<:Integer}

  T = eltype(g.points)

  poly_sec_deriv = Array{T}(undef, order + 1, length(g.points))
  poly_sec_deriv[1, :] .= zeros(T, length(g.points))

  @inbounds for j in eachindex(g.points)
    p = one(T)
    pl = NaN
    pll = NaN
    pd = zero(T)
    pdl = NaN
    pdll = NaN
    for i = 2:order+1
      if i == 2
        pl, p = p, normalize_node(g.points[j], g.domain)
        pdl, pd = pd, one(T)
        poly_sec_deriv[i, j] = zero(T)
      else
        pll, pl = pl, p
        p = 2 * normalize_node(g.points[j], g.domain) * pl - pll
        pdll, pdl = pdl, pd
        pd = 2 * pl + 2 * normalize_node(g.points[j], g.domain) * pdl - pdll
        poly_sec_deriv[i, j] = 2 * normalize_node(g.points[j], g.domain) * poly_sec_deriv[i-1, j] + 4 * pdl - poly_sec_deriv[i-2, j]
      end
    end
  end

  return ChebPoly(transpose(poly_sec_deriv), typeof(g))

end

function chebyshev_weights(y::AbstractArray{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  if typeof(plan) <: CApproxPlan

    nodes = Array{Array{T,1},1}(undef, length(plan.grid.grid))
    for i in eachindex(plan.grid.grid)
      nodes[i] = plan.grid.grid[i].points
    end

    if eltype(plan.grid.grid) <: ChebRoots
      return chebyshev_weights(y, tuple(nodes...), plan.order, plan.domain)
    elseif eltype(plan.grid.grid) <: ChebExtrema
      return chebyshev_weights_extrema(y, tuple(nodes...), plan.order, plan.domain)
    else # Extended_nodes, Vertesi_nodes, Legendre_nodes
      return chebyshev_weights_extended(y, tuple(nodes...), plan.order, plan.domain)
    end

  elseif typeof(plan) <: CApproxPlanPoly

    polynomials = Array{Array{T,2},1}(undef, length(plan.polys))
    for i in eachindex(plan.polys)
      polynomials[i] = plan.polys[i].poly
    end

    if plan.polys[1].nodetype == ChebRoots{T}
      return chebyshev_weights(y, tuple(polynomials...), plan.order)
    elseif plan.polys[1].nodetype == ChebExtrema{T}
      return chebyshev_weights_extrema(y, tuple(polynomials...), plan.order)
    else # Extended_nodes, Vertesi_nodes, Legendre_nodes
      return chebyshev_weights_extended(y, tuple(polynomials...), plan.order)
    end

  end

end

function chebyshev_weights_threaded(y::AbstractArray{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  if typeof(plan) <: CApproxPlan

    nodes = Array{Array{T,1},1}(undef, length(plan.grid.grid))
    for i in eachindex(plan.grid.grid)
      nodes[i] = plan.grid.grid[i].points
    end

    if eltype(plan.grid.grid) <: ChebRoots
      return chebyshev_weights_threaded(y, tuple(nodes...), plan.order, plan.domain)
    elseif eltype(plan.grid.grid) <: ChebExtrema
      return chebyshev_weights_extrema_threaded(y, tuple(nodes...), plan.order, plan.domain)
    else # Extended_nodes, Vertesi_nodes, Legendre_nodes
      return chebyshev_weights_extended_threaded(y, tuple(nodes...), plan.order, plan.domain)
    end

  elseif typeof(plan) <: CApproxPlanPoly

    polynomials = Array{Array{T,2},1}(undef, length(plan.polys))
    for i in eachindex(plan.polys)
      polynomials[i] = plan.polys[i].poly
    end

    if plan.polys[1].nodetype == ChebRoots{T}
      return chebyshev_weights_threaded(y, tuple(polynomials...), plan.order)
    elseif plan.polys[1].nodetype == ChebExtrema{T}
      return chebyshev_weights_extrema_threaded(y, tuple(polynomials...), plan.order)
    else # Extended_nodes, Vertesi_nodes, Legendre_nodes
      return chebyshev_weights_extended_threaded(y, tuple(polynomials...), plan.order)
    end

  end

end

function chebyshev_weights(f::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
  end

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      product = one(T)
      @inbounds for j = 1:N
        product *= poly[j][s[j], i[j]]
      end

      numerator += f[s] * product
      denominator += product^2

    end

    weights[i] = numerator / denominator

  end

  return weights

end

function chebyshev_weights_extrema(f::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  n = size(f)

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
  end

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      num = f[s]
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

function chebyshev_weights_extended(f::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

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

    @inbounds for s in CartesianIndices(f)

      num = f[s]
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

function chebyshev_weights(f::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      product = one(T)
      @inbounds for j = 1:N
        product *= poly[j][s[j], i[j]]
      end

      numerator += f[s] * product
      denominator += product^2

    end

    weights[i] = numerator / denominator

  end

  return weights

end

function chebyshev_weights_extrema(f::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  n = size(f)

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      num = f[s]
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

function chebyshev_weights_extended(f::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      num = f[s]
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

function chebyshev_weights(f::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
  end

  ord = (order,)
  for i = 2:N
    ord = (ord..., order)
  end

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        product = one(T)
        @inbounds for j = 1:N
          product *= poly[j][s[j], i[j]]
        end

        numerator += f[s] * product
        denominator += product^2

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

function chebyshev_weights_extrema(f::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  n = size(f)

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
  end

  ord = (order,)
  for i = 2:N
    ord = (ord..., order)
  end

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        num = f[s]
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

function chebyshev_weights_extended(f::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)
  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
    complement[i] = pinv(poly[i])'
  end

  ord = (order,)
  for i = 2:N
    ord = (ord..., order)
  end

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        num = f[s]
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

function chebyshev_weights(f::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = (order,)
  for i = 2:N
    ord = (ord..., order)
  end

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        product = one(T)
        @inbounds for j = 1:N
          product *= poly[j][s[j], i[j]]
        end

        numerator += f[s] * product
        denominator += product^2

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

function chebyshev_weights_extrema(f::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  n = size(f)

  ord = (order,)
  for i = 2:N
    ord = (ord..., order)
  end

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        num = f[s]
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

function chebyshev_weights_extended(f::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  complement = Array{Array{T,2},1}(undef, N)
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end

  ord = (order,)
  for i = 2:N
    ord = (ord..., order)
  end

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        num = f[s]
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

function chebyshev_weights_threaded(f::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
  end

  weights = zeros(Tuple(order .+ 1))

  @inbounds @sync @qthreads for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      product = one(T)
      @inbounds for j = 1:N
        product *= poly[j][s[j], i[j]]
      end

      numerator += f[s] * product
      denominator += product^2

    end

    weights[i] = numerator / denominator

  end

  return weights

end

function chebyshev_weights_extrema_threaded(f::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  n = size(f)

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
  end

  weights = zeros(Tuple(order .+ 1))

  @inbounds @sync @qthreads for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      num = f[s]
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

function chebyshev_weights_extended_threaded(f::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)
  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
    complement[i] = pinv(poly[i])'
  end

  weights = zeros(Tuple(order .+ 1))

  @inbounds @sync @qthreads for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      num = f[s]
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

function chebyshev_weights_threaded(f::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds @sync @qthreads for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      product = one(T)
      @inbounds for j = 1:N
        product *= poly[j][s[j], i[j]]
      end

      numerator += f[s] * product
      denominator += product^2

    end

    weights[i] = numerator / denominator

  end

  return weights

end

function chebyshev_weights_extrema_threaded(f::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  n = size(f)

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds @sync @qthreads for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      num = f[s]
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

function chebyshev_weights_extended_threaded(f::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  complement = Array{Array{T,2},1}(undef, N)
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end

  weights = Array{T,N}(undef, Tuple(order .+ 1))

  @inbounds @sync @qthreads for i in CartesianIndices(weights)

    numerator = zero(T)
    denominator = zero(T)

    @inbounds for s in CartesianIndices(f)

      num = f[s]
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

function chebyshev_weights_threaded(f::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
  end

  ord = (order,)
  for i = 2:N
    ord = (ord..., order)
  end

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        product = one(T)
        @inbounds for j = 1:N
          product *= poly[j][s[j], i[j]]
        end

        numerator += f[s] * product
        denominator += product^2

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

function chebyshev_weights_extrema_threaded(f::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  n = size(f)

  poly = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
  end

  ord = (order,)
  for i = 2:N
    ord = (ord..., order)
  end

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        num = f[s]
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

function chebyshev_weights_extended_threaded(f::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)
  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
    complement[i] = pinv(poly[i])'
  end

  ord = (order,)
  for i = 2:N
    ord = (ord..., order)
  end

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        num = f[s]
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

function chebyshev_weights_threaded(f::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = (order,)
  for i = 2:N
    ord = (ord..., order)
  end

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        product = one(T)
        @inbounds for j = 1:N
          product *= poly[j][s[j], i[j]]
        end

        numerator += f[s] * product
        denominator += product^2

      end

      weights[i] = numerator / denominator

    else
      weights[i] = zero(T)
    end

  end

  return weights

end

function chebyshev_weights_extrema_threaded(f::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  n = size(f)

  ord = (order,)
  for i = 2:N
    ord = (ord..., order)
  end

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        num = f[s]
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

function chebyshev_weights_extended_threaded(f::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  complement = Array{Array{T,2},1}(undef, N)
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end

  ord = (order,)
  for i = 2:N
    ord = (ord..., order)
  end

  weights = Array{T,N}(undef, ord .+ 1)

  @inbounds @sync @qthreads for i in CartesianIndices(weights)
    if sum(Tuple(i)) <= order + N

      numerator = zero(T)
      denominator = zero(T)

      @inbounds for s in CartesianIndices(f)

        num = f[s]
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

chebyshev_weights(f::Array{T,1}, nodes::Array{T,1}, order::Union{S,Array{S,1}}, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights(f, (nodes,), order, domain)
chebyshev_weights_extrema(f::Array{T,1}, nodes::Array{T,1}, order::Union{S,Array{S,1}}, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extrema(f, (nodes,), order, domain)
chebyshev_weights_extended(f::Array{T,1}, nodes::Array{T,1}, order::Union{S,Array{S,1}}, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extended(f, (nodes,), order, domain)
chebyshev_weights(f::Array{T,1}, poly::Array{T,2}, order::Union{S,Array{S,1}}) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights(f, (poly,), order)
chebyshev_weights_extrema(f::Array{T,1}, poly::Array{T,2}, order::Union{S,Array{S,1}}) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extrema(f, (poly,), order)
chebyshev_weights_extended(f::Array{T,1}, poly::Array{T,2}, order::Union{S,Array{S,1}}) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extended(f, (poly,), order)

# Functions that allow the nodes to be in an array of arrays

# Serial functions

chebyshev_weights(f::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights(f, tuple(nodes...), order, domain)
chebyshev_weights_extrema(f::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema(f, tuple(nodes...), order, domain)
chebyshev_weights_extended(f::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended(f, tuple(nodes...), order, domain)
chebyshev_weights(f::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights(f, tuple(poly...), order)
chebyshev_weights_extrema(f::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema(f, tuple(poly...), order)
chebyshev_weights_extended(f::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended(f, tuple(poly...), order)
chebyshev_weights(f::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights(f, tuple(nodes...), order, domain)
chebyshev_weights_extrema(f::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema(f, tuple(nodes...), order, domain)
chebyshev_weights_extended(f::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended(f, tuple(nodes...), order, domain)
chebyshev_weights(f::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights(f, tuple(poly...), order)
chebyshev_weights_extrema(f::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema(f, tuple(poly...), order)
chebyshev_weights_extended(f::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended(f, tuple(poly...), order)

# Threaded functions

chebyshev_weights_threaded(f::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded(f, tuple(nodes...), order, domain)
chebyshev_weights_extrema_threaded(f::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded(f, tuple(nodes...), order, domain)
chebyshev_weights_extended_threaded(f::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded(f, tuple(nodes...), order, domain)
chebyshev_weights_threaded(f::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded(f, tuple(poly...), order)
chebyshev_weights_extrema_threaded(f::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded(f, tuple(poly...), order)
chebyshev_weights_extended_threaded(f::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded(f, tuple(poly...), order)
chebyshev_weights_threaded(f::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded(f, tuple(nodes...), order, domain)
chebyshev_weights_extrema_threaded(f::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded(f, tuple(nodes...), order, domain)
chebyshev_weights_extended_threaded(f::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded(f, tuple(nodes...), order, domain)
chebyshev_weights_threaded(f::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded(f, tuple(poly...), order)
chebyshev_weights_extrema_threaded(f::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded(f, tuple(poly...), order)
chebyshev_weights_extended_threaded(f::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded(f, tuple(poly...), order)

const chebyshev_weights_vertesi = chebyshev_weights_extended
const chebyshev_weights_legendre = chebyshev_weights_extended

# Functions to evaluate Chebyshev polynominals

function chebyshev_evaluate(weights::AbstractArray{T,N}, x::AbstractArray{R,1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

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

function chebyshev_evaluate(weights::AbstractArray{T,N}, x::AbstractArray{R,1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

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

function chebyshev_interp(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights(y, plan)

  function interp(x::AbstractArray{R,1}) where {R<:Number}

    yhat = chebyshev_evaluate(w, x, plan.order, plan.domain)

    return yhat

  end

  return interp

end

function chebyshev_interp_threaded(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights_threaded(y, plan)

  function interp(x::AbstractArray{R,1}) where {R<:Number}

    yhat = chebyshev_evaluate(w, x, plan.order, plan.domain)

    return yhat

  end

  return interp

end

# Functions for derivatives

function chebyshev_derivative(weights::Array{T,N}, x::AbstractArray{R,1}, pos::S, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

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

function chebyshev_derivative(weights::Array{T,N}, x::AbstractArray{R,1}, pos::S, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

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

function chebyshev_gradient(weights::Array{T,N}, x::AbstractArray{R,1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  gradient = Array{R,2}(undef, 1, N)

  @inbounds for i = 1:N
    gradient[i] = chebyshev_derivative(weights, x, i, order, domain)
  end

  return gradient

end

function chebyshev_gradient(weights::Array{T,N}, x::AbstractArray{R,1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  gradient = Array{R,2}(undef, 1, N)

  @inbounds for i = 1:N
    gradient[i] = chebyshev_derivative(weights, x, i, order, domain)
  end

  return gradient

end

function chebyshev_gradient(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights(y, plan)

  function cheb_grad(x::AbstractArray{R,1}) where {R<:Number}

    grad = chebyshev_gradient(w, x, plan.order, plan.domain)

    return grad

  end

  return cheb_grad

end

function chebyshev_gradient_threaded(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights_threaded(y, plan)

  function cheb_grad(x::AbstractArray{R,1}) where {R<:Number}

    grad = chebyshev_gradient(w, x, plan.order, plan.domain)

    return grad

  end

  return cheb_grad

end

# Functions for hessians

function chebyshev_hessian(weights::Array{T,N}, x::AbstractArray{R,1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

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

function chebyshev_hessian(weights::Array{T,N}, x::AbstractArray{R,1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

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

function chebyshev_hessian(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights(y, plan)

  function cheb_hess(x::AbstractArray{R,1}) where {R<:Number}

    hess = chebyshev_hessian(w, x, plan.order, plan.domain)

    return hess

  end

  return cheb_hess

end

function chebyshev_hessian_threaded(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights_threaded(y, plan)

  function cheb_hess(x::AbstractArray{R,1}) where {R<:Number}

    hess = chebyshev_hessian(w, x, plan.order, plan.domain)

    return hess

  end

  return cheb_hess

end