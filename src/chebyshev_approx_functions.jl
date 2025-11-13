abstract type Nodes end
abstract type CApproximationPlan end

"""
ChebRoots is an immutable struct that contains two fields: the roots of the Chebyshev polynomial and the domain over
which they are located.
"""
struct ChebRoots{R<:AbstractFloat,T<:AbstractFloat} <: Nodes

  points::Array{R,1}
  domain::Array{T,1}

end

"""
ChebExtrema is an immutable struct that contains two fields: the extrema of the Chebyshev polynomial and the domain
over which they are located.
"""
struct ChebExtrema{R<:AbstractFloat,T<:AbstractFloat} <: Nodes

  points::Array{R,1}
  domain::Array{T,1}

end

"""
ChebExtended is an immutable struct that contains two fields: the extended oints of the Chebyshev polynomial and the
domain over which they are located.  The extended points are the Chebyshev roots scaled so that the first and last
points equal -1.0 and 1.0, respectively.
"""
struct ChebExtended{R<:AbstractFloat,T<:AbstractFloat} <: Nodes

  points::Array{R,1}
  domain::Array{T,1}

end

"""
Grid is an immutable struct that contains a single field: a tuple of Nodes subtypes (ChebRoots, ChebExtrema,
ChebExtended).
"""
struct Grid{G<:Nodes}

  grid::Tuple{Vararg{G}}

end

"""
ChebPoly is an immutable struct that contains two fields: an array that contains a Chebyshev polynomial of degree
```n```` evaluated over ```N``` points and a DataType that describes how the points were generated.  First- and
second-derivatives of Chebyshev polynomials are also represented using the ChebPoly type.
"""
struct ChebPoly{T<:AbstractFloat} # Holds polynomials as well as their first- and second-derivative

  poly::Array{T,2}
  nodetype::DataType

end

"""
CApproxPlan is an immutable struct that contains the information used to approximate a multi-dimensional function.
CApproxPlan has three fileds: A Grid type that summarizes the approximating points, the degree of the
approximating polynomials, and the approximation domain.
"""
struct CApproxPlan{G<:Grid,S<:Integer,T<:AbstractFloat} <: CApproximationPlan

  grid::G
  order::Union{S,Tuple{Vararg{S}}}
  domain::Union{Array{T,1},Array{T,2}}

end

"""
CApproxPoly is an immutable struct that contains the information used to approximate a multi-dimensional function.
CApproxPlan has three fileds: A Tuple that contains the univariate Chebyshev polynomials, the degree of the
approximating polynomials, and the approximation domain.
"""
struct CApproxPlanPoly{N,P<:ChebPoly,S<:Integer,T<:AbstractFloat} <: CApproximationPlan

  polys::NTuple{N,P}
  order::Union{S,Tuple{Vararg{S}}}
  domain::Union{Array{T,1},Array{T,2}}

end

"""
Computes ```N``` roots of the Chebyshev polynomial and scales the roots to the interval given in ```domain```
(defaults to [1.0,-1.0]).  Returns a vector containing the ```N``` Chebyshev roots located over the specified domain.

Signature
=========

points = chebyshev_nodes(N,domain)

Examples
========
```
julia> points = chebyshev_nodes(2)
[-0.7071067811865476
  0.7071067811865476]

julia> points = chebyshev_nodes(2,[3.0,-1.0])
[-0.41421356237309515
  2.414213562373095]
```
"""
function chebyshev_nodes(N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = fill((domain[1]+domain[2])*0.5,N)

  for i = 1:div(N,2)
    x = -cos((i - 0.5)*π/N)*(domain[1] - domain[2])*0.5
    points[i]     += x
    points[N-i+1] -= x
  end

  return points

end

"""
Computes ```N``` roots with floating point type ```T``` of the Chebyshev polynomial and scales the roots to the
interval given in ```domain``` (defaults to [1.0,-1.0]).  Returns a vector containing the ```N``` Chebyshev roots
with type ```T``` located over the specified domain.

Signature
=========

points = chebyshev_nodes(T,N,domain)

Examples
========
```
julia> using DoubleFloats
julia> points = chebyshev_nodes(Float64,2)
[-4.1421356237309504880168872420969394e-01
  2.41421356237309504880168872420969394]

julia> points = chebyshev_nodes(Double64,2,[3.0,-1.0])
[-4.1421356237309504880168872420969394e-01
  2.41421356237309504880168872420969394]

julia> points = chebyshev_nodes(BigFloat,2,[3.0,-1.0])
[-0.4142135623730950488016887242096980785696718753769480731766797379907324784621193
  2.414213562373095048801688724209698078569671875376948073176679737990732478462119]
```
"""
function chebyshev_nodes(T::DataType,N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = fill(T(domain[1]+domain[2])*0.5,N)

  for i = 1:div(N,2)
    x = -cos(T(i - 0.5)*π/N)*(domain[1] - domain[2])*0.5
    points[i]     += x
    points[N-i+1] -= x
  end

  return points

end

"""
Computes ```N``` extrema of the Chebyshev polynomial and scales the extrema to the interval given in ```domain```
(defaults to [1.0,-1.0]).  Returns a vector containing the ```N``` Chebyshev extrema located over the specified
domain.

Signature
=========

points = chebyshev_extrema(N,domain)

Examples
========
```
julia> points = chebyshev_extrema(4)
[-1.0
 -0.5000000000000001
  0.5000000000000001
  1.0]

julia> points = chebyshev_extrema(4,[3.0,-1.0])
[-1.0
 -2.220446049250313e-16
  2.0
  3.0]
```
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
Computes ```N``` extrema with floating point type ```T``` of the Chebyshev polynomial and scales the extrema to the
interval given in ```domain``` (defaults to [1.0,-1.0]).  Returns a vector containing the ```N``` Chebyshev extrema
with type ```T``` located over the specified domain.

Signature
=========

points = chebyshev_extrema(T,N,domain)

Examples
========
```
julia> using DoubleFloats
julia> points = chebyshev_extrema(Float64,4)
[-1.0
 -5.00000000000000000000000000000002311e-01
  5.00000000000000000000000000000002311e-01
  1.0]

julia> points = chebyshev_extrema(Double64,4,[3.0,-1.0])
[-1.0
 -4.622231866529366e-33
  2.00000000000000000000000000000000462
  3.0]

  julia> points = chebyshev_extrema(BigFloat,4,[3.0,-1.0])
[-1.0
 -1.727233711018888925077270372560079914223200072887256277004740694033718360632485e-77
  2.0
  3.0]
```
"""
function chebyshev_extrema(T::DataType, N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = fill(T(domain[1]+domain[2])*0.5,N)

  @inbounds for i = 1:div(N,2)
    x = -cos(T(i - 1)*π/(N - 1))*(domain[1] - domain[2])*0.5
    points[i]     += x
    points[N-i+1] -= x
  end

  return points

end

"""
Computes ```N``` Chebyshev extended points and scales the points to the interval given in ```domain``` (defaults to
[1.0,-1.0]).  Returns a vector containing the ```N``` Chebyshev extended points located over the specified domain.

Signature
=========

points = chebyshev_extended(N,domain)

Examples
========
```
julia> points = chebyshev_extended(4)
[-1.0
 -0.41421356237309515
  0.41421356237309515
  1.0]

julia> points = chebyshev_extended(4,[3.0,-1.0])
[-1.0
  0.1715728752538097
  1.8284271247461903
  3.0]
```
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
Computes ```N``` Chebyshev extended points with floating point type ```T``` and scales the points to the interval
given in ```domain``` (defaults to [1.0,-1.0]).  Returns a vector containing the ```N``` Chebyshev extended points
of type ```T``` located over the specified domain.

Signature
=========

points = chebyshev_extended(T,N,domain)

Examples
========
```
julia> using DoubleFloats
julia> points = chebyshev_extended(Float64,4)
[-1.00000000000000001909886133787806654
 -4.1421356237309505671269611624194948e-01
  4.1421356237309505671269611624194948e-01
  1.00000000000000001909886133787806654]

julia> points = chebyshev_extended(Double64,4,[3.0,-1.0])
[-1.00000000000000003819772267575613308
  1.71572875253809886574607767516101039e-01
  1.82842712474619011342539223248389896
  3.00000000000000003819772267575613308]

  julia> points = chebyshev_extended(BigFloat,4,[3.0,-1.0])
[-1.000000000000000038197722675756133178078150949228553435907665016103921602647826
  0.1715728752538098865746077675161036699979185261802141131475309447067879422084131
  1.828427124746190113425392232483896330002081473819785886852469055293212057791578
  3.000000000000000038197722675756133178078150949228553435907665016103921602647826]
```
"""
function chebyshev_extended(T::DataType, N::S, domain=[1.0, -1.0]) where {S<:Integer}

  points = fill(T(domain[1]+domain[2])*0.5,N)

  @inbounds for i = 1:div(N,2)
    x = -cos(T(i - 0.5)*π/N)*((domain[1] - domain[2])*0.5)/cos(π/(2N))
    points[i]     += x
    points[N-i+1] -= x
  end

  return points

end

"""
Computes ```N``` points according to ```node_generator``` and scales the points to the interval given in
```domain``` (defaults to [1.0,-1.0]).  Return the points in a struct (ChebRoots, ChebExtrema, ChebExtended)
according to the ```node_generator```.

Signature
=========

points = nodes(N,node_generator,domain)

Examples
========
```
julia> points = nodes(4,:chebyshev_nodes)
julia> points = nodes(4,:chebyshev_extrema)
julia> points = nodes(4,:chebyshev_extended)

julia> points = nodes(4,:chebyshev_nodes,[3.0,-1.0])
julia> points = nodes(4,:chebyshev_extrema,[3.0,-1.0])
julia> points = nodes(4,:chebyshev_extended,[3.0,-1.0])
```
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
Computes ```N``` points of type ```T``` according to ```node_generator``` and scales the points to the interval
given in ```domain``` (defaults to [1.0,-1.0]).  Return the points in a struct (ChebRoots, ChebExtrema,
ChebExtended) according to the ```node_generator```.

Signature
=========

points = nodes(T,N,node_generator,domain)

Examples
========
```
julia> using DoubleFloats
julia> points = nodes(Double64,4,:chebyshev_nodes)
julia> points = nodes(Double64,4,:chebyshev_extrema)
julia> points = nodes(Double64,4,:chebyshev_extended)

julia> points = nodes(Double64,4,:chebyshev_nodes,[3.0,-1.0])
julia> points = nodes(Double64,4,:chebyshev_extrema,[3.0,-1.0])
julia> points = nodes(Double64,4,:chebyshev_extended,[3.0,-1.0])
```
"""
function nodes(T::DataType, N::S, node_generator::Symbol, domain=[1.0, -1.0],) where {S<:Integer}

  if node_generator == :chebyshev_nodes
    return ChebRoots(chebyshev_nodes(T, N, domain), domain)
  elseif node_generator == :chebyshev_extrema
    return ChebExtrema(chebyshev_extrema(T, N, domain), domain)
  elseif node_generator == :chebyshev_extended
    return ChebExtended(chebyshev_extended(T, N, domain), domain)
  end

end

function normalize_node(node::R, domain::Array{T,1}) where {R<:Number,T<:AbstractFloat} # Internal function, not exported

  if domain[1] == domain[2]
    norm_node = zero(T)
    return norm_node
  else
    norm_node = 2*(node - domain[2])/(domain[1] - domain[2]) - 1.0
    return norm_node
  end

end

function normalize_node(nodes::Array{R,1}, domain::Array{T,1}) where {R<:Number,T<:AbstractFloat} # Internal function, not exported

  norm_nodes = map(x -> normalize_node(x, domain), nodes)

  return norm_nodes

end

function normalize_node(nodes::G) where {G<:Nodes} # Internal function, not exported

  norm_nodes = map(x -> normalize_node(x, nodes.domain), nodes.points)

  return norm_nodes

end

"""
Computes the Chebyshev polynomial of degree ```order``` at the scalar point ```x```, which must lie in [1.0,-1.0].
Returns a 2d array.  The element-type of the polynomial is given by the element-type of ```x```.

Signature
=========

P = chebyshev_polynomial(order,x)

Example
=======
```
julia> P = chebyshev_polynomial(3,0.6)
[1.0  0.6  -0.28  -0.936]
```
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
Computes the Chebyshev polynomial of degree ```order``` at end point in the vector ```x```. The elements of ```x```
must lie in [1.0,-1.0].  Returns a 2d array.  The element-type of the polynomial is given by the element-type of
```x```.

Signature
=========

P = chebyshev_polynomial(order,x)

Example
=======
```
julia> P = chebyshev_polynomial(3,[0.6,0.4])
[1.0  0.6  -0.28  -0.936
 1.0  0.4  -0.68  -0.944]
```
"""
function chebyshev_polynomial(order::S, x::AbstractArray{R,1}) where {S<:Integer,R<:Number}

  # Elements of x must reside in [-1,1]
  poly = Array{R}(undef, length(x), order + 1)
  poly[:, 1] .= one(R)

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
Computes the Chebyshev polynomial of degree ```order``` at the scalar point ```x```, which must lie in ```dom```.
Returns a 2d array.  The element-type of the polynomial is given by the element-type of ```x```.

Signature
=========

P = chebyshev_polynomial(order,x,dom)

Example
=======
```
julia> P = chebyshev_polynomial(3,0.6,[1.0,-1.0])
[1.0  0.6  -0.28  -0.936]
```
"""
function chebyshev_polynomial(order::S, x::R, dom::Array{T,1}) where {S<:Integer,R<:Number,T<:AbstractFloat}

  # x must reside in [-1,1]

  point = normalize_node(x,dom)

  poly = Array{R}(undef, 1, order + 1)
  poly[1] = one(R)

  @inbounds for i = 2:order+1
    if i == 2
      poly[i] = point
    else
      poly[i] = 2*point*poly[i-1] - poly[i-2]
    end
  end

  return poly

end

"""
Computes the Chebyshev polynomial of degree ```order``` at end point in the vector ```x```. The elements of ```x```
must lie in ```dom```.  Returns a 2d array.  The element-type of the polynomial is given by the element-type of
```x```.

Signature
=========

P = chebyshev_polynomial(order,x,dom)

Example
=======
```
julia> P = chebyshev_polynomial(3,[0.6,0.4],[1.0,-1.0])
[1.0  0.6  -0.28  -0.936
 1.0  0.4  -0.68  -0.944]
```
"""
function chebyshev_polynomial(order::S, x::AbstractArray{R,1}, dom::Array{T,1}) where {S<:Integer,R<:Number,T<:AbstractFloat}

  points = normalize_node(x,dom
  )
  # Elements of x must reside in [-1,1]
  poly = Array{R}(undef, length(x), order + 1)
  poly[:, 1] .= one(R)

  @inbounds for i = 2:order+1
    for j in eachindex(x)
      if i == 2
          poly[j, i] = points[j]
      else
        poly[j, i] = 2*points[j]*poly[j,i-1] - poly[j,i-2]
      end
    end
  end

  return poly

end

"""
Computes the Chebyshev polynomial of degree ```order``` at each point specified in the ```nodes``` struct.  The
elements of ```nodes.points``` must lie in [1.0,-1.0].  Returns a ChebPoly type.  The element-type of the polynomial
is given by the element-type of ```x```.

Signature
=========

P = chebyshev_polynomial(order,x)

Example
=======
```
julia> n = nodes(2,:chebyshev_nodes)
ChebRoots{Float64, Float64}([-0.7071067811865476, 0.7071067811865476], [1.0, -1.0])

julia> P = chebyshev_polynomial(3,n)
ChebPoly{Float64}([1.0 -0.7071067811865476 2.220446049250313e-16 0.7071067811865472; 1.0 0.7071067811865475 -2.220446049250313e-16 -0.7071067811865478], ChebRoots{Float64, Float64})
```
"""
function chebyshev_polynomial(order::S, nodes::G) where {S<:Integer,G<:Nodes}

  T = eltype(nodes.points)

  poly = Array{T}(undef, length(nodes.points), order + 1)
  poly[:, 1] .= one(T)

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
Computes the first derivative of the Chebyshev polynomial of ```order``` at the scalar ```x```, which must lie in
[1.0,-1.0].  Returns a 2d array.  The element-type of the polynomial is given by the element-type of ```x```.

Signature
=========

P = chebyshev_polynomial_deriv(order,x)

Example
=======
```
julia> P = chebyshev_polynomial_deriv(3,0.6)
[0.0  1.0  2.4  1.32]
```
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
Computes the first derivative of the Chebyshev polynomial of ```order``` at each point in the vector ```x```.  The
elements of ```x``` must lie in [1.0,-1.0].  Returns a 2d array.  The element-type of the polynomial is given by
the element-type of ```x```.

Signature
=========

P = chebyshev_polynomial_deriv(order,x)

Example
=======
```
julia> P = chebyshev_polynomial_deriv(3,[0.6,0.4])
[0.0  1.0  2.4   1.32
 0.0  1.0  1.6  -1.08]
```
"""
function chebyshev_polynomial_deriv(order::S, x::AbstractArray{R,1}) where {S<:Integer,R<:Number}

  poly_deriv = Array{R}(undef, length(x), order + 1)
  poly_deriv[:, 1] .= zero(R)

  @inbounds for j in eachindex(x)
    p = one(R)
    pl = NaN
    pll = NaN
    for i = 2:order+1
      if i == 2
        pl, p = p, x[j]
        poly_deriv[j,i] = one(R)
      else
        pll, pl = pl, p
        p = 2*x[j]*pl - pll
        poly_deriv[j,i] = 2*pl + 2*x[j]*poly_deriv[j,i-1] - poly_deriv[j,i-2]
      end
    end
  end

  return poly_deriv

end

"""
Computes the first derivative of the Chebyshev polynomial of ```order``` at each point specified in the ```nodes```
struct.  The elements of ```nodes.points``` must lie in [1.0,-1.0].  Returns a ChebPoly type.  The element-type of
the polynomial is given by the element-type of ```x```.

Signature
=========

P = chebyshev_polynomial_deriv(order,x)

Example
=======
```
julia> n = nodes(2,:chebyshev_nodes)
ChebRoots{Float64, Float64}([-0.7071067811865476, 0.7071067811865476], [1.0, -1.0])

julia> P = chebyshev_polynomial_deriv(3,n)
ChebPoly{Float64}([0.0 1.0 -2.8284271247461903 3.0000000000000018; 0.0 1.0 2.82842712474619 2.9999999999999987], ChebRoots{Float64, Float64})
```
"""
function chebyshev_polynomial_deriv(order::S, nodes::G) where {S<:Integer,G<:Nodes}

  T = eltype(nodes.points)

  poly_deriv = Array{T}(undef,length(nodes.points), order + 1)
  poly_deriv[:,1] .= zero(T)

  @inbounds for j in eachindex(nodes.points)
    p = one(T)
    pl = NaN
    pll = NaN
    for i = 2:order+1
      if i == 2
        pl, p = p, normalize_node(nodes.points[j],nodes.domain)
        poly_deriv[j,i] = one(T)
      else
        pll, pl = pl, p
        p = 2*normalize_node(nodes.points[j],nodes.domain)*pl - pll
        poly_deriv[j,i] = 2*pl + 2*normalize_node(nodes.points[j],nodes.domain)*poly_deriv[j,i-1] - poly_deriv[j,i-2]
      end
    end
  end

  return ChebPoly(poly_deriv, G)

end

"""
Computes the second derivative of the Chebyshev polynomial of ```order``` at the scalar ```x```, which must lie in
[1.0,-1.0].  Returns a 2d array.  The element-type of the polynomial is given by the element-type of ```x```.

Signature
=========

P = chebyshev_polynomial_sec_deriv(order,x)

Example
=======
```
julia> P = chebyshev_polynomial_sec_deriv(3,0.6)
[0.0  0.0  4.0  14.4]
```
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
Computes the second derivative of the Chebyshev polynomial of ```order``` at each poiunt in the vector ```x```. The
elements of ```x``` must lie in [1.0,-1.0].  Returns a 2d array.  The element-type of the polynomial is given by the
element-type of ```x```.

Signature
=========

P = chebyshev_polynomial_sec_deriv(order,x)

Example
=======
```
julia> P = chebyshev_polynomial_sec_deriv(3,[0.6,0.4])
[0.0  0.0  4.0  14.4
 0.0  0.0  4.0   9.6]
```
"""
function chebyshev_polynomial_sec_deriv(order::S, x::AbstractArray{T,1}) where {S<:Integer,T<:Number}

  poly_sec_deriv = Array{T}(undef, length(x), order + 1)
  poly_sec_deriv[:,1] .= zero(T)

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
        poly_sec_deriv[j,i] = zero(T)
      else
        pll, pl = pl, p
        p = 2*x[j]*pl - pll
        pdll, pdl = pdl, pd
        pd = 2*pl + 2*x[j]*pdl - pdll
        poly_sec_deriv[j,i] = 2*x[j]*poly_sec_deriv[j,i-1] + 4*pdl - poly_sec_deriv[j,i-2]
      end
    end
  end

  return poly_sec_deriv

end

"""
Computes the second derivative of the Chebyshev polynomial of ```order``` at each point specified in the ```nodes```
struct.  The elements of ```nodes.points``` must lie in [1.0,-1.0].  Returns a ChebPoly type.  The element-type of
the polynomial is given by the element-type of ```x```.

Signature
=========

P = chebyshev_polynomial_sec_deriv(order,x)

Example
=======
```
julia> n = nodes(2,:chebyshev_nodes)
ChebRoots{Float64, Float64}([-0.7071067811865476, 0.7071067811865476], [1.0, -1.0])

julia> P = chebyshev_polynomial_sec_deriv(3,n)
ChebPoly{Float64}([0.0 0.0 4.0 -16.970562748477143; 0.0 0.0 4.0 16.97056274847714], ChebRoots{Float64, Float64})
```
"""
function chebyshev_polynomial_sec_deriv(order::S,nodes::G) where {G<:Nodes,S<:Integer}

  T = eltype(nodes.points)

  poly_sec_deriv = Array{T}(undef, length(nodes.points), order + 1)
  poly_sec_deriv[:,1] .= zero(T)

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
        poly_sec_deriv[j,i] = zero(T)
      else
        pll, pl = pl, p
        p = 2*normalize_node(nodes.points[j],nodes.domain)*pl - pll
        pdll, pdl = pdl, pd
        pd = 2 * pl + 2 * normalize_node(nodes.points[j], nodes.domain) * pdl - pdll
        poly_sec_deriv[j,i] = 2*normalize_node(nodes.points[j],nodes.domain)*poly_sec_deriv[j,i-1] + 4*pdl - poly_sec_deriv[j,i-2]
      end
    end
  end

  return ChebPoly(poly_sec_deriv, G)

end

"""
Computes the weights in a Chebyshev polynomial given the approximation sample, ```y```, and the approximation plan,
```plan```.  Returns a multi-dimensional array containing the weights in the specified ```plan```.  The
element-type of the Chebyshev weights is given by the element-type of the approximating points specified in the
approximation ```plan```.

Signature
=========

w = chebyshev_weights(y,plan)

Examples
========
```
julia> nodes1 = nodes(5,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_nodes,[1.5,0.5])
julia> g = Grid((nodes1,nodes2))
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> plan = CApproxPlan(g,ord,dom)
julia> weights = chebyshev_weights(y,plan)
[1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]

julia> nodes1 = nodes(9,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(9,:chebyshev_nodes,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:9, j in 1:9]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> plan = CApproxPlanPoly((poly1,poly2),ord,dom)
julia> weights = chebyshev_weights(y,plan)
[1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]
```
"""
function chebyshev_weights(y::AbstractArray{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  if typeof(plan) <: CApproxPlan

    nodes = Tuple(plan.grid.grid[i].points for i in eachindex(plan.grid.grid))

    if eltype(plan.grid.grid) <: ChebRoots
      return chebyshev_weights(y, nodes, plan.order, plan.domain)
    elseif eltype(plan.grid.grid) <: ChebExtrema
      return chebyshev_weights_extrema(y, nodes, plan.order, plan.domain)
    elseif eltype(plan.grid.grid) <: ChebExtended
      return chebyshev_weights_extended(y, nodes, plan.order, plan.domain)
    end

  elseif typeof(plan) <: CApproxPlanPoly

    polynomials = Tuple(plan.polys[i].poly for i in eachindex(plan.polys))

    if plan.polys[1].nodetype <: ChebRoots
      return chebyshev_weights(y, polynomials, plan.order)
    elseif plan.polys[1].nodetype <: ChebExtrema
      return chebyshev_weights_extrema(y, polynomials, plan.order)
    elseif plan.polys[1].nodetype <: ChebExtended
      return chebyshev_weights_extended(y, polynomials, plan.order)
    end

  end

end

"""
Computes the weights in a tensor-product Chebyshev polynomial given the approximation sample, ```y```, the Chebyshev
roots, ```nodes```, the ```order``` of the polynomial, and the ```domain```.  Returns a multi-dimensional array
containing the weights.  The element-type of the Chebyshev weights is given by the element-type of the ```nodes```.

Signature
=========

w = chebyshev_weights(y,nodes,order,domain)

Example
=======
```
julia> nodes1 = chebyshev_nodes(5,[9.0,3.0])
julia> nodes2 = chebyshev_nodes(5,[1.5,0.5])
julia> y = [nodes1[i]^0.3*nodes2[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> weights = chebyshev_weights(y,(nodes1,nodes2),ord,dom)
[1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]
```
"""
function chebyshev_weights(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, (order .+ 1)...) # allocates
  return chebyshev_weights!(weights, y, nodes, order, domain)

end

"""
In-place variant of ```chebyshev_weights(y, nodes, order, domain)```.
See the documentation for the out-of-place version. **Warning**: for performance, no bounds checking is done, so be sure that weights and order are consistent with each other.
"""
function chebyshev_weights!(weights::Array{T,N}, y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = [chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i])) for i in 1:N] # allocates

  @inbounds for i in CartesianIndices(weights) # loop does not allocate

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
Computes the weights in a tensor-product Chebyshev polynomial given the approximation sample, ```y```, the Chebyshev
extrema, ```nodes```, the ```order``` of the polynomial, and the ```domain```.  Returns a multi-dimensional array
containing the weights.  The element-type of the Chebyshev weights is given by the element-type of the ```nodes```.

Signature
=========

w = chebyshev_weights_extrema(y,nodes,order,domain)

Example
=======
```
julia> nodes1 = chebyshev_extrema(5,[9.0,3.0])
julia> nodes2 = chebyshev_extrema(5,[1.5,0.5])
julia> y = [nodes1[i]^0.3*nodes2[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> weights = chebyshev_weights_extrema(y,(nodes1,nodes2),ord,dom)
[1.66416      0.598461   -0.02372       0.00280618
  0.263786     0.0948622  -0.00375987    0.000444808
 -0.024647    -0.0088635   0.000351305  -4.15608e-5
  0.00386243   0.001389   -5.5053e-5     6.513e-6]
```
"""
function chebyshev_weights_extrema(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, (order .+ 1)...) # allocates
  return chebyshev_weights_extrema!(weights, y, nodes, order, domain)

end

"""
In-place variant of ```chebyshev_weights_extrema(y, nodes, order, domain)```.
See the documentation for the out-of-place version. **Warning**: for performance, no bounds checking is done, so be sure that weights and order are consistent with each other.
"""
function chebyshev_weights_extrema!(weights::Array{T,N}, y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y) # allocates

  poly = [chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i])) for i in 1:N] # allocates

  @inbounds for i in CartesianIndices(weights)

  numerator = zero(T)
  denominator = zero(T)

  @inbounds for s in CartesianIndices(y)

    num = y[s]
    den = one(T)
    @inbounds for j = 1:N
      scale = 1.0
      if s[j] == 1 || s[j] == n[j] # allocates
        scale = 0.5
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
Computes the weights in a tensor-product Chebyshev polynomial given the approximation sample, ```y```, the Chebyshev
extended points, ```nodes```, the ```order``` of the polynomial, and the ```domain```.  Returns a multi-dimensional
array containing the weights.  The element-type of the Chebyshev weights is given by the element-type of the
```nodes```.

Signature
=========

w = chebyshev_weights_extended(y,nodes,order,domain)

Example
=======
```
julia> nodes1 = chebyshev_extended(5,[9.0,3.0])
julia> nodes2 = chebyshev_extended(5,[1.5,0.5])
julia> y = [nodes1[i]^0.3*nodes2[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> weights = chebyshev_weights_extended(y,(nodes1,nodes2),ord,dom)
[1.66393      0.598444    -0.0238922     0.00277041
  0.263837     0.0948909   -0.00378841    0.000439284
 -0.0249247   -0.00896435   0.000357891  -4.14991e-5
  0.00379791   0.00136595  -5.45338e-5    6.32345e-6]
```
"""
function chebyshev_weights_extended(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, (order .+ 1)...)
  return chebyshev_weights_extended!(weights, y, nodes, order, domain)

end

"""
In-place variant of ```chebyshev_weights_extended(y, nodes, order, domain)```.
 See the documentation for the out-of-place version. **Warning**: for performance, no bounds checking is done, so be sure that weights and order are consistent with each other.
"""
function chebyshev_weights_extended!(weights::Array{T,N}, y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)
  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
    complement[i] = pinv(poly[i])'
  end

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
Computes the weights in a tensor-product Chebyshev polynomial given the approximation sample, ```y```, the Chebyshev
polynomials, ```poly``` evaluated at the Chebyshev roots, and the ```order``` of the polynomial.  Returns a
multi-dimensional array containing the weights.  The element-type of the Chebyshev weights is given by the
element-type of ```poly```.

Signature
=========

w = chebyshev_weights(y,poly,order)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_nodes,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> weights = chebyshev_weights(y,(poly1.poly,poly2.poly),ord)
[1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]
```
"""
function chebyshev_weights(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, (order .+ 1)...)
  return chebyshev_weights!(weights, y, poly, order)

end

"""
In-place variant of ```chebyshev_weights(y, poly, order)```.
See the documentation for the out-of-place version. **Warning**: for performance, no bounds checking is done, so be sure that weights, poly, and order are consistent with each other.
"""
function chebyshev_weights!(weights::Array{T,N}, y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

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
Computes the weights in a tensor-product Chebyshev polynomial given the approximation sample, ```y```, the Chebyshev
polynomials, ```poly``` evaluated at the Chebyshev extrema, and the ```order``` of the polynomial.  Returns a
multi-dimensional array containing the weights.  The element-type of the Chebyshev weights is given by the
element-type of ```poly```.

Signature
=========

w = chebyshev_weights_extrema(y,poly,order)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_extrema,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_extrema,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> weights = chebyshev_weights_extrema(y,(poly1.poly,poly2.poly),ord)
[1.66416      0.598461   -0.02372       0.00280618
  0.263786     0.0948622  -0.00375987    0.000444808
 -0.024647    -0.0088635   0.000351305  -4.15608e-5
  0.00386243   0.001389   -5.5053e-5     6.513e-6]
```
"""
function chebyshev_weights_extrema(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, (order .+ 1)...)
  return chebyshev_weights_extrema!(weights, y, poly, order)

end

"""
In-place variant of ```chebyshev_weights_extrema(y, poly, order)```.
See the documentation for the out-of-place version. **Warning**: for performance, no bounds checking is done, so be sure that weights, poly, and order are consistent with each other.
"""
function chebyshev_weights_extrema!(weights::Array{T,N}, y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

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
Computes the weights in a tensor-product Chebyshev polynomial given the approximation sample, ```y```, the Chebyshev
polynomials, ```poly``` evaluated at the extended Chebyshev points, and the ```order``` of the polynomial.  Returns a
multi-dimensional array containing the weights.  The element-type of the Chebyshev weights is given by the
element-type of ```poly```.

Signature
=========

w = chebyshev_weights_extended(y,poly,order)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_extended,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_extended,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> weights = chebyshev_weights_extended(y,(poly1.poly,poly2.poly),ord)
[1.66393      0.598444    -0.0238922     0.00277041
  0.263837     0.0948909   -0.00378841    0.000439284
 -0.0249247   -0.00896435   0.000357891  -4.14991e-5
  0.00379791   0.00136595  -5.45338e-5    6.32345e-6]
```
"""
function chebyshev_weights_extended(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, (order .+ 1)...)
  return chebyshev_weights_extended!(weights, y, poly, order)

end

"""
In-place variant of ```chebyshev_weights_extended(y, poly, order)````.
See the documentation for the out-of-place version. **Warning**: for performance, no bounds checking is done, so be sure that weights, poly, and order are consistent with each other.
"""
function chebyshev_weights_extended!(weights::Array{T,N}, y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end

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
Computes the weights in a complete Chebyshev polynomial given the approximation sample, ```y```, the Chebyshev
roots, ```nodes```, the ```order``` of the polynomial, and the ```domain```.  Returns a multi-dimensional array
containing the weights.  The element-type of the Chebyshev weights is given by the element-type of the ```nodes```.

Signature
=========

w = chebyshev_weights(y,nodes,order,domain)

Example
=======
```
julia> nodes1 = chebyshev_nodes(5,[9.0,3.0])
julia> nodes2 = chebyshev_nodes(5,[1.5,0.5])
julia> y = [nodes1[i]^0.3*nodes2[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = 3
julia> dom = [9.0 1.5; 3.0 0.5]
julia> weights = chebyshev_weights(y,(nodes1,nodes2),ord,dom)
[1.66416      0.598458    -0.0237052   0.00272941
  0.26378      0.0948592   -0.00375743  0.0
 -0.0246176   -0.00885287   0.0         0.0
  0.00372291   0.0          0.0         0.0]
```
"""
function chebyshev_weights(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)
  weights = Array{T,N}(undef, ord .+ 1)
  return chebyshev_weights!(weights, y, nodes, order, domain)

end

"""
In-place variant for scalar-order complete Chebyshev weights (roots):
```chebyshev_weights!(weights, y, nodes, order, domain)```.
"""
function chebyshev_weights!(weights::Array{T,N}, y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = [chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i])) for i in 1:N] # allocates

  ord = Tuple(order for _ in 1:N)

  @inbounds for i in CartesianIndices(weights)
    if sum(i.I) <= order + N

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
Computes the weights in a complete Chebyshev polynomial given the approximation sample, ```y```, the Chebyshev
extrema, ```nodes```, the ```order``` of the polynomial, and the ```domain```.  Returns a multi-dimensional array
containing the weights.  The element-type of the Chebyshev weights is given by the element-type of the ```nodes```.

Signature
=========

w = chebyshev_weights(y,nodes,order,domain)

Example
=======
```
julia> nodes1 = chebyshev_extrema(5,[9.0,3.0])
julia> nodes2 = chebyshev_extrema(5,[1.5,0.5])
julia> y = [nodes1[i]^0.3*nodes2[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = 3
julia> dom = [9.0 1.5; 3.0 0.5]
julia> weights = chebyshev_weights_extrema(y,(nodes1,nodes2),ord,dom)
[1.66416      0.598461   -0.02372     0.00280618
  0.263786     0.0948622  -0.00375987  0.0
 -0.024647    -0.0088635   0.0         0.0
  0.00386243   0.0         0.0         0.0]
```
"""
function chebyshev_weights_extrema(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)
  weights = Array{T,N}(undef, ord .+ 1)
  return chebyshev_weights_extrema!(weights, y, nodes, order, domain)

end

"""
In-place variant for scalar-order complete Chebyshev weights (extrema):
```chebyshev_weights_extrema!(weights, y, nodes, order, domain)```.
"""
function chebyshev_weights_extrema!(weights::Array{T,N}, y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  poly = [chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i])) for i in 1:N] # allocates

  ord = Tuple(order for _ in 1:N)

  @inbounds for i in CartesianIndices(weights)
    if sum(i.I) <= order + N

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
Computes the weights in a complete Chebyshev polynomial given the approximation sample, ```y```, the extended
Chebyshev points, ```nodes```, the ```order``` of the polynomial, and the ```domain```.  Returns a multi-dimensional
array containing the weights.  The element-type of the Chebyshev weights is given by the element-type of the ```nodes```.

Signature
=========

w = chebyshev_weights_extended(y,nodes,order,domain)

Example
=======
```
julia> nodes1 = chebyshev_extended(5,[9.0,3.0])
julia> nodes2 = chebyshev_extended(5,[1.5,0.5])
julia> y = [nodes1[i]^0.3*nodes2[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = 3
julia> dom = [9.0 1.5; 3.0 0.5]
julia> weights = chebyshev_weights_extended(y,(nodes1,nodes2),ord,dom)
[1.66393      0.598444    -0.0238922   0.00277041
  0.263837     0.0948909   -0.00378841  0.0
 -0.0249247   -0.00896435   0.0         0.0
  0.00379791   0.0          0.0         0.0]
```
"""
function chebyshev_weights_extended(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)
  weights = Array{T,N}(undef, ord .+ 1)
  return chebyshev_weights_extended!(weights, y, nodes, order, domain)

end

"""
In-place variant for scalar-order complete Chebyshev weights (extended):
```chebyshev_weights_extended!(weights, y, nodes, order, domain)```.
"""
function chebyshev_weights_extended!(weights::Array{T,N}, y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)
  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
    complement[i] = pinv(poly[i])'
  end

  ord = Tuple(order for _ in 1:N)

  @inbounds for i in CartesianIndices(weights)
    if sum(i.I) <= order + N

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
Computes the weights in a complete Chebyshev polynomial given the approximation sample, ```y```, the Chebyshev
polynomials, ```poly``` evaluated at the Chebyshev roots, and the ```order``` of the polynomial.  Returns a
multi-dimensional array containing the weights.  The element-type of the Chebyshev weights is given by the
element-type of ```poly```.

Signature
=========

w = chebyshev_weights(y,poly,order)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_nodes,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = 3
julia> weights = chebyshev_weights(y,(poly1.poly,poly2.poly),ord)
[1.66416      0.598458    -0.0237052   0.00272941
  0.26378      0.0948592   -0.00375743  0.0
 -0.0246176   -0.00885287   0.0         0.0
  0.00372291   0.0          0.0         0.0]
```
"""
function chebyshev_weights(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)
  weights = Array{T,N}(undef, ord .+ 1)
  return chebyshev_weights!(weights, y, poly, order)

end

"""
In-place variant for scalar-order complete Chebyshev weights (roots, pre-computed poly):
```chebyshev_weights!(weights, y, poly, order)```.
"""
function chebyshev_weights!(weights::Array{T,N}, y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)

  @inbounds for i in CartesianIndices(weights)
    if sum(i.I) <= order + N

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
Computes the weights in a complete Chebyshev polynomial given the approximation sample, ```y```, the Chebyshev
polynomials, ```poly``` evaluated at the Chebyshev extrema, and the ```order``` of the polynomial.  Returns a
multi-dimensional array containing the weights.  The element-type of the Chebyshev weights is given by the
element-type of ```poly```.

Signature
=========

w = chebyshev_weights_extrema(y,poly,order)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_extrema,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_extrema,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = 3
julia> weights = chebyshev_weights_extrema(y,(poly1.poly,poly2.poly),ord)
[1.66416      0.598461   -0.02372     0.00280618
  0.263786     0.0948622  -0.00375987  0.0
 -0.024647    -0.0088635   0.0         0.0
  0.00386243   0.0         0.0         0.0]
```
"""
function chebyshev_weights_extrema(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)
  weights = Array{T,N}(undef, ord .+ 1)
  return chebyshev_weights_extrema!(weights, y, poly, order)

end

"""
In-place variant for scalar-order complete Chebyshev weights (extrema, pre-computed poly):
```chebyshev_weights_extrema!(weights, y, poly, order)```.
"""
function chebyshev_weights_extrema!(weights::Array{T,N}, y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  ord = Tuple(order for _ in 1:N)

  @inbounds for i in CartesianIndices(weights)
    if sum(i.I) <= order + N

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
Computes the weights in a complete Chebyshev polynomial given the approximation sample, ```y```, the Chebyshev
polynomials, ```poly``` evaluated at the extended Chebyshev points, and the ```order``` of the polynomial.
Returns a multi-dimensional array containing the weights.  The element-type of the Chebyshev weights is given
by the element-type of ```poly```.

Signature
=========

w = chebyshev_weights_extended(y,poly,order)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_extended,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_extended,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = 3
julia> weights = chebyshev_weights_extended(y,(poly1.poly,poly2.poly),ord)
[1.66393      0.598444    -0.0238922   0.00277041
  0.263837     0.0948909   -0.00378841  0.0
 -0.0249247   -0.00896435   0.0         0.0
  0.00379791   0.0          0.0         0.0]
```
"""
function chebyshev_weights_extended(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)
  weights = Array{T,N}(undef, ord .+ 1)
  return chebyshev_weights_extended!(weights, y, poly, order)

end

"""
In-place variant for scalar-order complete Chebyshev weights (extended, pre-computed poly):
```chebyshev_weights_extended!(weights, y, poly, order)```.
"""
function chebyshev_weights_extended!(weights::Array{T,N}, y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  complement = Array{Array{T,2},1}(undef, N)
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end

  ord = Tuple(order for _ in 1:N)

  @inbounds for i in CartesianIndices(weights)
    if sum(i.I) <= order + N

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
Computes the weights in a Chebyshev polynomial using multi-threading given the approximation sample, ```y```, and
the approximation plan, ```plan```.  Returns a multi-dimensional array containing the weights in the specified
```plan```.  The element-type of the Chebyshev weights is given by the element-type of the approximating points
specified in the  approximation ```plan```.

Signature
=========

w = chebyshev_weights_threaded(y,plan)

Examples
========
```
julia> nodes1 = nodes(5,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_nodes,[1.5,0.5])
julia> g = Grid((nodes1,nodes2))
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> plan = CApproxPlan(g,ord,dom)
julia> weights = chebyshev_weights_threaded(y,plan)
[1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]

julia> nodes1 = nodes(9,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(9,:chebyshev_nodes,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:9, j in 1:9]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> plan = CApproxPlanPoly((poly1,poly2),ord,dom)
julia> weights = chebyshev_weights_threaded(y,plan)
[1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]
```
"""
function chebyshev_weights_threaded(y::AbstractArray{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  if typeof(plan) <: CApproxPlan

    nodes = Tuple(plan.grid.grid[i].points for i in eachindex(plan.grid.grid))

    if eltype(plan.grid.grid) <: ChebRoots
      return chebyshev_weights_threaded(y, nodes, plan.order, plan.domain)
    elseif eltype(plan.grid.grid) <: ChebExtrema
      return chebyshev_weights_extrema_threaded(y, nodes, plan.order, plan.domain)
    elseif eltype(plan.grid.grid) <: ChebExtended
      return chebyshev_weights_extended_threaded(y, nodes, plan.order, plan.domain)
    end

  elseif typeof(plan) <: CApproxPlanPoly

    polynomials = Tuple(plan.polys[i].poly for i in eachindex(plan.polys))

    if plan.polys[1].nodetype <: ChebRoots
      return chebyshev_weights_threaded(y, polynomials, plan.order)
    elseif plan.polys[1].nodetype <: ChebExtrema
      return chebyshev_weights_extrema_threaded(y, polynomials, plan.order)
    elseif plan.polys[1].nodetype <: ChebExtended
      return chebyshev_weights_extended_threaded(y, polynomials, plan.order)
    end

  end

end

"""
Computes the weights in a tensor-product Chebyshev polynomial using multi-threading given the approximation sample,
```y```, the Chebyshev roots, ```nodes```, the ```order``` of the polynomial, and the ```domain```.  Returns a
multi-dimensional array containing the weights.  The element-type of the Chebyshev weights is given by the
element-type of the ```nodes```.

Signature
=========

w = chebyshev_weights_threaded(y,nodes,order,domain)

Example
=======
```
julia> nodes1 = chebyshev_nodes(5,[9.0,3.0])
julia> nodes2 = chebyshev_nodes(5,[1.5,0.5])
julia> y = [nodes1[i]^0.3*nodes2[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> weights = chebyshev_weights_threaded(y,(nodes1,nodes2),ord,dom)
[1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]
```
"""
function chebyshev_weights_threaded(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = zeros(T,(order .+ 1)...)
  return chebyshev_weights_threaded!(weights, y, nodes, order, domain)

end

"""
In-place variant for tensor-product Chebyshev weights (roots, threaded):
```chebyshev_weights_threaded!(weights, y, nodes, order, domain)```.
"""
function chebyshev_weights_threaded!(weights::Array{T,N}, y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = [chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i])) for i in 1:N] # allocates

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
Computes the weights in a tensor-product Chebyshev polynomial using multi-threading given the approximation sample,
```y```, the Chebyshev extrema, ```nodes```, the ```order``` of the polynomial, and the ```domain```.  Returns a
multi-dimensional array containing the weights.  The element-type of the Chebyshev weights is given by the
element-type of the ```nodes```.

Signature
=========

w = chebyshev_weights_extrema_threaded(y,nodes,order,domain)

Example
=======
```
julia> nodes1 = chebyshev_extrema(5,[9.0,3.0])
julia> nodes2 = chebyshev_extrema(5,[1.5,0.5])
julia> y = [nodes1[i]^0.3*nodes2[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> weights = chebyshev_weights_extrema_threaded(y,(nodes1,nodes2),ord,dom)
[1.66416      0.598461   -0.02372       0.00280618
  0.263786     0.0948622  -0.00375987    0.000444808
 -0.024647    -0.0088635   0.000351305  -4.15608e-5
  0.00386243   0.001389   -5.5053e-5     6.513e-6]
```
"""
function chebyshev_weights_extrema_threaded(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = zeros(T,(order .+ 1)...)
  return chebyshev_weights_extrema_threaded!(weights, y, nodes, order, domain)

end

"""
In-place variant for tensor-product Chebyshev weights (extrema, threaded):
```chebyshev_weights_extrema_threaded!(weights, y, nodes, order, domain)```.
"""
function chebyshev_weights_extrema_threaded!(weights::Array{T,N}, y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  poly = [chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i])) for i in 1:N] # allocates

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
Computes the weights in a tensor-product Chebyshev polynomial using multi-threading given the approximation sample,
```y```, the Chebyshev extended points, ```nodes```, the ```order``` of the polynomial, and the ```domain```.
Returns a multi-dimensional array containing the weights.  The element-type of the Chebyshev weights is given by the
element-type of the ```nodes```.

Signature
=========

w = chebyshev_weights_extended_threaded(y,nodes,order,domain)

Example
=======
```
julia> nodes1 = chebyshev_extended(5,[9.0,3.0])
julia> nodes2 = chebyshev_extended(5,[1.5,0.5])
julia> y = [nodes1[i]^0.3*nodes2[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> weights = chebyshev_weights_extended_threaded(y,(nodes1,nodes2),ord,dom)
[1.66393      0.598444    -0.0238922     0.00277041
  0.263837     0.0948909   -0.00378841    0.000439284
 -0.0249247   -0.00896435   0.000357891  -4.14991e-5
  0.00379791   0.00136595  -5.45338e-5    6.32345e-6]
```
"""
function chebyshev_weights_extended_threaded(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  weights = zeros(T,(order .+ 1)...)
  return chebyshev_weights_extended_threaded!(weights, y, nodes, order, domain)

end

"""
In-place variant for tensor-product Chebyshev weights (extended, threaded):
```chebyshev_weights_extended_threaded!(weights, y, nodes, order, domain)```.
"""
function chebyshev_weights_extended_threaded!(weights::Array{T,N}, y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)
  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order[i], normalize_node(nodes[i], domain[:, i]))
    complement[i] = pinv(poly[i])'
  end

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
Computes the weights in a tensor-product Chebyshev polynomial using multi-threading given the approximation sample,
```y```, the Chebyshev polynomials, ```poly``` evaluated at the Chebyshev roots, and the ```order``` of the
polynomial.  Returns a multi-dimensional array containing the weights.  The element-type of the Chebyshev weights
is given by the element-type of ```poly```.

Signature
=========

w = chebyshev_weights_threaded(y,poly,order)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_nodes,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> weights = chebyshev_weights_threaded(y,(poly1.poly,poly2.poly),ord)
[1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]
```
"""
function chebyshev_weights_threaded(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, (order .+ 1)...)
  return chebyshev_weights_threaded!(weights, y, poly, order)

end

"""
In-place variant for tensor-product Chebyshev weights (pre-computed poly, threaded):
```chebyshev_weights_threaded!(weights, y, poly, order)```.
"""
function chebyshev_weights_threaded!(weights::Array{T,N}, y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

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
Computes the weights in a tensor-product Chebyshev polynomial using multi-threading given the approximation sample,
```y```, the Chebyshev polynomials, ```poly``` evaluated at the Chebyshev extrema, and the ```order``` of the
polynomial.  Returns a multi-dimensional array containing the weights.  The element-type of the Chebyshev weights
is given by the element-type of ```poly```.

Signature
=========

w = chebyshev_weights_extrema_threaded(y,poly,order)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_extrema,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_extrema,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> weights = chebyshev_weights_extrema_threaded(y,(poly1.poly,poly2.poly),ord)
[1.66416      0.598461   -0.02372       0.00280618
  0.263786     0.0948622  -0.00375987    0.000444808
 -0.024647    -0.0088635   0.000351305  -4.15608e-5
  0.00386243   0.001389   -5.5053e-5     6.513e-6]
```
"""
function chebyshev_weights_extrema_threaded(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, (order .+ 1)...)
  return chebyshev_weights_extrema_threaded!(weights, y, poly, order)

end

"""
In-place variant for tensor-product Chebyshev weights (extrema, pre-computed poly, threaded):
```chebyshev_weights_extrema_threaded!(weights, y, poly, order)```.
"""
function chebyshev_weights_extrema_threaded!(weights::Array{T,N}, y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

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
Computes the weights in a tensor-product Chebyshev polynomial using multi-threadinggiven the approximation sample,
```y```, the Chebyshev polynomials, ```poly``` evaluated at the extended Chebyshev points, and the ```order``` of
the polynomial.  Returns a multi-dimensional array containing the weights.  The element-type of the Chebyshev
weights is given by the element-type of ```poly```.

Signature
=========

w = chebyshev_weights_extended_threaded(y,poly,order)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_extended,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_extended,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> weights = chebyshev_weights_extended_threaded(y,(poly1.poly,poly2.poly),ord)
[1.66393      0.598444    -0.0238922     0.00277041
  0.263837     0.0948909   -0.00378841    0.000439284
 -0.0249247   -0.00896435   0.000357891  -4.14991e-5
  0.00379791   0.00136595  -5.45338e-5    6.32345e-6]
```
"""
function chebyshev_weights_extended_threaded(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  weights = Array{T,N}(undef, (order .+ 1)...)
  return chebyshev_weights_extended_threaded!(weights, y, poly, order)

end

"""
In-place variant for tensor-product Chebyshev weights (extended, pre-computed poly, threaded):
```chebyshev_weights_extended_threaded!(weights, y, poly, order)```.
"""
function chebyshev_weights_extended_threaded!(weights::Array{T,N}, y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer}

  complement = Array{Array{T,2},1}(undef, N)
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end

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
Computes the weights in a complete Chebyshev polynomial using multi-threading given the approximation sample,
```y```, the Chebyshev roots, ```nodes```, the ```order``` of the polynomial, and the ```domain```.  Returns
a multi-dimensional array containing the weights.  The element-type of the Chebyshev weights is given by the
element-type of the ```nodes```.

Signature
=========

w = chebyshev_weights_threaded(y,nodes,order,domain)

Example
=======
```
julia> nodes1 = chebyshev_nodes(5,[9.0,3.0])
julia> nodes2 = chebyshev_nodes(5,[1.5,0.5])
julia> y = [nodes1[i]^0.3*nodes2[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = 3
julia> dom = [9.0 1.5; 3.0 0.5]
julia> weights = chebyshev_weights_threaded(y,(nodes1,nodes2),ord,dom)
[1.66416      0.598458    -0.0237052   0.00272941
  0.26378      0.0948592   -0.00375743  0.0
 -0.0246176   -0.00885287   0.0         0.0
  0.00372291   0.0          0.0         0.0]
```
"""
function chebyshev_weights_threaded(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)
  weights = Array{T,N}(undef, ord .+ 1)
  return chebyshev_weights_threaded!(weights, y, nodes, order, domain)

end

"""
In-place variant for scalar-order complete Chebyshev weights (roots, threaded):
```chebyshev_weights_threaded!(weights, y, nodes, order, domain)```.
"""
function chebyshev_weights_threaded!(weights::Array{T,N}, y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = [chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i])) for i in 1:N] # allocates

  ord = Tuple(order for _ in 1:N)

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)
    if sum(i.I) <= order + N

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
Computes the weights in a complete Chebyshev polynomial using multi-threading given the approximation sample,
```y```, the Chebyshev extrema, ```nodes```, the ```order``` of the polynomial, and the ```domain```.  Returns
a multi-dimensional array containing the weights.  The element-type of the Chebyshev weights is given by the
element-type of the ```nodes```.

Signature
=========

w = chebyshev_weights_threaded(y,nodes,order,domain)

Example
=======
```
julia> nodes1 = chebyshev_extrema(5,[9.0,3.0])
julia> nodes2 = chebyshev_extrema(5,[1.5,0.5])
julia> y = [nodes1[i]^0.3*nodes2[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = 3
julia> dom = [9.0 1.5; 3.0 0.5]
julia> weights = chebyshev_weights_extrema_threaded(y,(nodes1,nodes2),ord,dom)
[1.66416      0.598461   -0.02372     0.00280618
  0.263786     0.0948622  -0.00375987  0.0
 -0.024647    -0.0088635   0.0         0.0
  0.00386243   0.0         0.0         0.0]
```
"""
function chebyshev_weights_extrema_threaded(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)
  weights = Array{T,N}(undef, ord .+ 1)
  return chebyshev_weights_extrema_threaded!(weights, y, nodes, order, domain)

end

"""
In-place variant for scalar-order complete Chebyshev weights (extrema, threaded):
```chebyshev_weights_extrema_threaded!(weights, y, nodes, order, domain)```.
"""
function chebyshev_weights_extrema_threaded!(weights::Array{T,N}, y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  poly = [chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i])) for i in 1:N] # allocates

  ord = Tuple(order for _ in 1:N)

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)
    if sum(i.I) <= order + N

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
Computes the weights in a complete Chebyshev polynomial using multi-threading given the approximation sample,
```y```, the extended Chebyshev points, ```nodes```, the ```order``` of the polynomial, and the ```domain```.
Returns a multi-dimensional array containing the weights.  The element-type of the Chebyshev weights is given
by the element-type of the ```nodes```.

Signature
=========

w = chebyshev_weights_extended_threaded(y,nodes,order,domain)

Example
=======
```
julia> nodes1 = chebyshev_extended(5,[9.0,3.0])
julia> nodes2 = chebyshev_extended(5,[1.5,0.5])
julia> y = [nodes1[i]^0.3*nodes2[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = 3
julia> dom = [9.0 1.5; 3.0 0.5]
julia> weights = chebyshev_weights_extended_threaded(y,(nodes1,nodes2),ord,dom)
[1.66393      0.598444    -0.0238922   0.00277041
  0.263837     0.0948909   -0.00378841  0.0
 -0.0249247   -0.00896435   0.0         0.0
  0.00379791   0.0          0.0         0.0]
```
"""
function chebyshev_weights_extended_threaded(y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)
  weights = Array{T,N}(undef, ord .+ 1)
  return chebyshev_weights_extended_threaded!(weights, y, nodes, order, domain)

end

"""
In-place variant for scalar-order complete Chebyshev weights (extended, threaded):
```chebyshev_weights_extended_threaded!(weights, y, nodes, order, domain)```.
"""
function chebyshev_weights_extended_threaded!(weights::Array{T,N}, y::Array{T,N}, nodes::NTuple{N,Array{T,1}}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer}

  poly = Array{Array{T,2},1}(undef, N)
  complement = Array{Array{T,2},1}(undef, N)

  @inbounds for i = 1:N
    poly[i] = chebyshev_polynomial(order, normalize_node(nodes[i], domain[:, i]))
    complement[i] = pinv(poly[i])'
  end

  ord = Tuple(order for _ in 1:N)

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)
    if sum(i.I) <= order + N

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
Computes the weights in a complete Chebyshev polynomial using multi-threading given the approximation sample,
```y```, the Chebyshev polynomials, ```poly``` evaluated at the Chebyshev roots, and the ```order``` of the
polynomial.  Returns a multi-dimensional array containing the weights.  The element-type of the Chebyshev
weights is given by the element-type of ```poly```.

Signature
=========

w = chebyshev_weights_threaded(y,poly,order)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_nodes,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = 3
julia> weights = chebyshev_weights_threaded(y,(poly1.poly,poly2.poly),ord)
[1.66416      0.598458    -0.0237052   0.00272941
  0.26378      0.0948592   -0.00375743  0.0
 -0.0246176   -0.00885287   0.0         0.0
  0.00372291   0.0          0.0         0.0]
```
"""
function chebyshev_weights_threaded(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)
  weights = Array{T,N}(undef, ord .+ 1)
  return chebyshev_weights_threaded!(weights, y, poly, order)

end

"""
In-place variant for scalar-order complete Chebyshev weights (roots, pre-computed poly, threaded):
```chebyshev_weights_threaded!(weights, y, poly, order)```.
"""
function chebyshev_weights_threaded!(weights::Array{T,N}, y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)
    if sum(i.I) <= order + N

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
Computes the weights in a complete Chebyshev polynomial using multi-threading given the approximation sample,
```y```, the Chebyshev polynomials, ```poly``` evaluated at the Chebyshev extrema, and the ```order``` of the
polynomial.  Returns a multi-dimensional array containing the weights.  The element-type of the Chebyshev
weights is given by the element-type of ```poly```.

Signature
=========

w = chebyshev_weights_extrema_threaded(y,poly,order)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_extrema,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_extrema,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = 3
julia> weights = chebyshev_weights_extrema_threaded(y,(poly1.poly,poly2.poly),ord)
[1.66416      0.598461   -0.02372     0.00280618
  0.263786     0.0948622  -0.00375987  0.0
 -0.024647    -0.0088635   0.0         0.0
  0.00386243   0.0         0.0         0.0]
```
"""
function chebyshev_weights_extrema_threaded(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)
  weights = Array{T,N}(undef, ord .+ 1)
  return chebyshev_weights_extrema_threaded!(weights, y, poly, order)

end

"""
In-place variant for scalar-order complete Chebyshev weights (extrema, pre-computed poly, threaded):
```chebyshev_weights_extrema_threaded!(weights, y, poly, order)```.
"""
function chebyshev_weights_extrema_threaded!(weights::Array{T,N}, y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  n = size(y)

  ord = Tuple(order for _ in 1:N)

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)
    if sum(i.I) <= order + N

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
Computes the weights in a complete Chebyshev polynomial using multi-threading given the approximation sample,
```y```, the Chebyshev polynomials, ```poly``` evaluated at the extended Chebyshev points, and the ```order```
of the polynomial.  Returns a multi-dimensional array containing the weights.  The element-type of the
Chebyshev weights is given by the element-type of ```poly```.

Signature
=========

w = chebyshev_weights_extended_threaded(y,poly,order)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_extended,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_extended,[1.5,0.5])
julia> poly1 = chebyshev_polynomial(3,nodes1)
julia> poly2 = chebyshev_polynomial(3,nodes2)
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = 3
julia> weights = chebyshev_weights_extended_threaded(y,(poly1.poly,poly2.poly),ord)
[1.66393      0.598444    -0.0238922   0.00277041
  0.263837     0.0948909   -0.00378841  0.0
 -0.0249247   -0.00896435   0.0         0.0
  0.00379791   0.0          0.0         0.0]
```
"""
function chebyshev_weights_extended_threaded(y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  ord = Tuple(order for _ in 1:N)
  weights = Array{T,N}(undef, ord .+ 1)
  return chebyshev_weights_extended_threaded!(weights, y, poly, order)

end

"""
In-place variant for scalar-order complete Chebyshev weights (extended, pre-computed poly, threaded):
```chebyshev_weights_extended_threaded!(weights, y, poly, order)```.
"""
function chebyshev_weights_extended_threaded!(weights::Array{T,N}, y::Array{T,N}, poly::NTuple{N,Array{T,2}}, order::S) where {T<:AbstractFloat,N,S<:Integer}

  complement = Array{Array{T,2},1}(undef, N)
  @inbounds for i = 1:N
    complement[i] = pinv(poly[i])'
  end

  ord = Tuple(order for _ in 1:N)

  @inbounds @sync Threads.@threads for i in CartesianIndices(weights)
    if sum(i.I) <= order + N

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

# Serial functions

# Tensor-product case
chebyshev_weights(y::Array{T,1}, nodes::Array{T,1}, order::Array{S,1}, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights(y, (nodes,), order, domain)
chebyshev_weights_extrema(y::Array{T,1}, nodes::Array{T,1}, order::Array{S,1}, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extrema(y, (nodes,), order, domain)
chebyshev_weights_extended(y::Array{T,1}, nodes::Array{T,1}, order::Array{S,1}, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extended(y, (nodes,), order, domain)
chebyshev_weights(y::Array{T,1}, poly::Array{T,2}, order::Array{S,1}) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights(y, (poly,), order)
chebyshev_weights_extrema(y::Array{T,1}, poly::Array{T,2}, order::Array{S,1}) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extrema(y, (poly,), order)
chebyshev_weights_extended(y::Array{T,1}, poly::Array{T,2}, order::Array{S,1}) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extended(y, (poly,), order)

# Complete polynomial case defaults to tensor-product case
chebyshev_weights(y::Array{T,1}, nodes::Array{T,1}, order::S, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights(y, (nodes,), [order], domain)
chebyshev_weights_extrema(y::Array{T,1}, nodes::Array{T,1}, order::S, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extrema(y, (nodes,), [order], domain)
chebyshev_weights_extended(y::Array{T,1}, nodes::Array{T,1}, order::S, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extended(y, (nodes,), [order], domain)
chebyshev_weights(y::Array{T,1}, poly::Array{T,2}, order::S) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights(y, (poly,), order)
chebyshev_weights_extrema(y::Array{T,1}, poly::Array{T,2}, order::S) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extrema(y, (poly,), [order])
chebyshev_weights_extended(y::Array{T,1}, poly::Array{T,2}, order::S) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extended(y, (poly,), [order])

# Threaded functions

# Tensor-product case
chebyshev_weights_threaded(y::Array{T,1}, nodes::Array{T,1}, order::Array{S,1}, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_threaded(y, (nodes,), order, domain)
chebyshev_weights_extrema_threaded(y::Array{T,1}, nodes::Array{T,1}, order::Array{S,1}, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extrema_threaded(y, (nodes,), order, domain)
chebyshev_weights_extended_threaded(y::Array{T,1}, nodes::Array{T,1}, order::Array{S,1}, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extended_threaded(y, (nodes,), order, domain)
chebyshev_weights_threaded(y::Array{T,1}, poly::Array{T,2}, order::Array{S,1}) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_threaded(y, (poly,), order)
chebyshev_weights_extrema_threaded(y::Array{T,1}, poly::Array{T,2}, order::Array{S,1}) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extrema_threaded(y, (poly,), order)
chebyshev_weights_extended_threaded(y::Array{T,1}, poly::Array{T,2}, order::Array{S,1}) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extended_threaded(y, (poly,), order)

# Complete polynomial case defaults to tensor-product case
chebyshev_weights_threaded(y::Array{T,1}, nodes::Array{T,1}, order::S, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_threaded(y, (nodes,), [order], domain)
chebyshev_weights_extrema_threaded(y::Array{T,1}, nodes::Array{T,1}, order::S, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extrema_threaded(y, (nodes,), [order], domain)
chebyshev_weights_extended_threaded(y::Array{T,1}, nodes::Array{T,1}, order::S, domain=[one(T); -one(T)]) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extended_threaded(y, (nodes,), [order], domain)
chebyshev_weights_threaded(y::Array{T,1}, poly::Array{T,2}, order::S) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_threaded(y, (poly,), order)
chebyshev_weights_extrema_threaded(y::Array{T,1}, poly::Array{T,2}, order::S) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extrema_threaded(y, (poly,), [order])
chebyshev_weights_extended_threaded(y::Array{T,1}, poly::Array{T,2}, order::S) where {T<:AbstractFloat,S<:Integer} = chebyshev_weights_extended_threaded(y, (poly,), [order])

# Functions that allow the nodes to be in an array of arrays

# Serial functions

# Tensor-product case
chebyshev_weights(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights(y, Tuple(nodes), order, domain)
chebyshev_weights_extrema(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema(y, Tuple(nodes), order, domain)
chebyshev_weights_extended(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended(y, Tuple(nodes), order, domain)
chebyshev_weights(y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights(y, Tuple(poly), order)
chebyshev_weights_extrema(y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema(y, Tuple(poly), order)
chebyshev_weights_extended(y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended(y, Tuple(poly), order)

chebyshev_weights!(w::Array{T,N}, y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights!(w, y, Tuple(nodes), order, domain)
chebyshev_weights_extrema!(w::Array{T,N}, y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema!(w, y, Tuple(nodes), order, domain)
chebyshev_weights_extended!(w::Array{T,N}, y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended!(w, y, Tuple(nodes), order, domain)
chebyshev_weights!(w::Array{T,N}, y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights!(w, y, Tuple(poly), order)
chebyshev_weights_extrema!(w::Array{T,N}, y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema!(w, y, Tuple(poly), order)
chebyshev_weights_extended!(w::Array{T,N}, y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended!(w, y, Tuple(poly), order)

# Complete polynomial case
chebyshev_weights(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights(y, Tuple(nodes), order, domain)
chebyshev_weights_extrema(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema(y, Tuple(nodes), order, domain)
chebyshev_weights_extended(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended(y, Tuple(nodes), order, domain)
chebyshev_weights(y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights(y, Tuple(poly), order)
chebyshev_weights_extrema(y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema(y, Tuple(poly), order)
chebyshev_weights_extended(y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended(y, Tuple(poly), order)

chebyshev_weights!(w::Array{T,N}, y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights!(w, y, Tuple(nodes), order, domain)
chebyshev_weights_extrema!(w::Array{T,N}, y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema!(w, y, Tuple(nodes), order, domain)
chebyshev_weights_extended!(w::Array{T,N}, y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended!(w, y, Tuple(nodes), order, domain)
chebyshev_weights!(w::Array{T,N}, y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights!(w, y, Tuple(poly), order)
chebyshev_weights_extrema!(w::Array{T,N}, y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema!(w, y, Tuple(poly), order)
chebyshev_weights_extended!(w::Array{T,N}, y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended!(w, y, Tuple(poly), order)

# Threaded functions

# Tensor-product case
chebyshev_weights_threaded(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded(y, Tuple(nodes), order, domain)
chebyshev_weights_extrema_threaded(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded(y, Tuple(nodes), order, domain)
chebyshev_weights_extended_threaded(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded(y, Tuple(nodes), order, domain)
chebyshev_weights_threaded(y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded(y, Tuple(poly), order)
chebyshev_weights_extrema_threaded(y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded(y, Tuple(poly), order)
chebyshev_weights_extended_threaded(y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded(y, Tuple(poly), order)

chebyshev_weights_threaded!(w::Array{T,N}, y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded!(w, y, Tuple(nodes), order, domain)
chebyshev_weights_extrema_threaded!(w::Array{T,N}, y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded!(w, y, Tuple(nodes), order, domain)
chebyshev_weights_extended_threaded!(w::Array{T,N}, y::Array{T,N}, nodes::Array{Array{T,1},1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded!(w, y, Tuple(nodes), order, domain)
chebyshev_weights_threaded!(w::Array{T,N}, y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded!(w, y, Tuple(poly), order)
chebyshev_weights_extrema_threaded!(w::Array{T,N}, y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded!(w, y, Tuple(poly), order)
chebyshev_weights_extended_threaded!(w::Array{T,N}, y::Array{T,N}, poly::Array{Array{T,2},1}, order::Union{NTuple{N,S},Array{S,1}}) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded!(w, y, Tuple(poly), order)

# Complete polynomial case
chebyshev_weights_threaded(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded(y, Tuple(nodes), order, domain)
chebyshev_weights_extrema_threaded(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded(y, Tuple(nodes), order, domain)
chebyshev_weights_extended_threaded(y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded(y, Tuple(nodes), order, domain)
chebyshev_weights_threaded(y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded(y, Tuple(poly), order)
chebyshev_weights_extrema_threaded(y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded(y, Tuple(poly), order)
chebyshev_weights_extended_threaded(y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded(y, Tuple(poly), order)

chebyshev_weights_threaded!(w::Array{T,N}, y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded!(w, y, Tuple(nodes), order, domain)
chebyshev_weights_extrema_threaded!(w::Array{T,N}, y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded!(w, y, Tuple(nodes), order, domain)
chebyshev_weights_extended_threaded!(w::Array{T,N}, y::Array{T,N}, nodes::Array{Array{T,1},1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded!(w, y, Tuple(nodes), order, domain)
chebyshev_weights_threaded!(w::Array{T,N}, y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_threaded!(w, y, Tuple(poly), order)
chebyshev_weights_extrema_threaded!(w::Array{T,N}, y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extrema_threaded!(w, y, Tuple(poly), order)
chebyshev_weights_extended_threaded!(w::Array{T,N}, y::Array{T,N}, poly::Array{Array{T,2},1}, order::S) where {T<:AbstractFloat,N,S<:Integer} = chebyshev_weights_extended_threaded!(w, y, Tuple(poly), order)

# Functions to evaluate Chebyshev polynominals

"""
Evaluate a tensor-product Chebyshev polynomial at point, ```x````, given the ```weights``` the ```order``` of the polynomial, and the ```domain```.

Signature
=========

yhat = chebyshev_evaluate(weights,x,order,domain)

Example
=======
```
julia> weights = [1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]
julia> x = [5.5,0.9]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> yhat = chebyshev_evaluate(weights,x,ord,dom)
1.5500018804115083
```
"""
function chebyshev_evaluate(weights::AbstractArray{T,N}, x::AbstractArray{R,1}, order::Union{NTuple{N,S},Array{S,1}}, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  if length(x) != N
    error("A value for 'x' is needed for each spacial dimension.")
  end

  poly = [chebyshev_polynomial(order[i], normalize_node(x[i], domain[:, i])) for i in 1:N] # allocates

  yhat = chebyshev_evaluate(weights, poly)

  return yhat

end

"""
Evaluate a complete Chebyshev polynomial at point, ```x````, given the ```weights``` the ```order``` of the polynomial, and the ```domain```.

Signature
=========

yhat = chebyshev_evaluate(weights,x,order,domain)
Example
=======
```
julia> weights = [1.66416      0.598458    -0.0237052   0.00272941
  0.26378      0.0948592   -0.00375743  0.0
 -0.0246176   -0.00885287   0.0         0.0
  0.00372291   0.0          0.0         0.0]
julia> x = [5.5,0.9]
julia> ord = 3
julia> dom = [9.0 1.5; 3.0 0.5]
julia> yhat = chebyshev_evaluate(weights,x,ord,dom)
1.5498202486133335
```
"""
function chebyshev_evaluate(weights::AbstractArray{T,N}, x::AbstractArray{R,1}, order::S, domain=[ones(T, 1, N); -ones(T, 1, N)]) where {T<:AbstractFloat,R<:Number,N,S<:Integer}

  if length(x) != N
    error("A value for 'x' is needed for each spacial dimension.")
  end

  poly = [chebyshev_polynomial(order, normalize_node(x[i], domain[:, i])) for i in 1:N] # allocates

  yhat = chebyshev_evaluate(weights, poly, order)

  return yhat

end

"""
Evaluate a tensor-product Chebyshev polynomial, given the ```weights``` and pre-computed ```poly``` matrices.

Signature
=========

yhat = chebyshev_evaluate(weights,poly)

Example
=======
```
julia> weights = [1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]
julia> x = [5.5,0.9]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> poly = [chebyshev_polynomial(ord[i], ChebyshevApprox.normalize_node(x[i], dom[:, i])) for i in 1:2]
julia> yhat = chebyshev_evaluate(weights,poly)
1.5500018804115083
```
"""
function chebyshev_evaluate(weights::AbstractArray{T1,N}, poly::Union{NTuple{N,<:AbstractArray{T2,2}},AbstractArray{<:AbstractArray{T2,2},1}}) where {T1<:AbstractFloat,T2<:Number,N}

  yhat = zero(T1)
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
Evaluate a complete Chebyshev polynomial, given the ```weights```, pre-computed ```poly``` matrices, and ```order```.

Signature
=========

yhat = chebyshev_evaluate(weights,poly,order)

Example
=======
```
julia> weights = [1.66416      0.598458    -0.0237052   0.00272941
  0.26378      0.0948592   -0.00375743  0.0
 -0.0246176   -0.00885287   0.0         0.0
  0.00372291   0.0          0.0         0.0]
julia> x = [5.5,0.9]
julia> ord = 3
julia> dom = [9.0 1.5; 3.0 0.5]
julia> poly = [chebyshev_polynomial(ord, ChebyshevApprox.normalize_node(x[i], dom[:, i])) for i in 1:2]
julia> yhat = chebyshev_evaluate(weights,poly,ord)
1.5498202486133335
```
"""
function chebyshev_evaluate(weights::AbstractArray{T1,N}, poly::Union{NTuple{N,<:AbstractArray{T2,2}},AbstractArray{<:AbstractArray{T2,2},1}}, order::S) where {T1<:AbstractFloat,T2<:Number,N,S<:Integer}

  yhat = zero(T1)
  @inbounds for i in CartesianIndices(weights)
    if sum(i.I) <= order + N
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
Create an interpolating function, given the approximation sample, ```y``` and approximation ```plan```.

Signature
=========

f = chebyshev_interp(y,plan)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_nodes,[1.5,0.5])
julia> g = Grid((nodes1,nodes2))
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> plan = CApproxPlan(g,ord,dom)
julia> f = chebyshev_interp(y,plan)
julia> f([5.5,0.9])
1.550003750285832
```
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
Create an interpolating function using multi-threading, given the approximation sample, ```y``` and approximation ```plan```.

Signature
=========

f = chebyshev_interp_threaded(y,plan)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_nodes,[1.5,0.5])
julia> g = Grid((nodes1,nodes2))
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> plan = CApproxPlan(g,ord,dom)
julia> f = chebyshev_interp_threaded(y,plan)
julia> f([5.5,0.9])
1.550003750285832
```
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
Computes the first derivative of a tensor-product Chebyshev polynomial with respect to the variable in position ```pos```, at point ```x```, given the ```weights```, the polynomial ```order```, and the ```domain```.

Signature
=========

d = chebyshev_derivative(weights,x,pos,order,domain)

Example
=======
```
julia> weights = [1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]
julia> x = [5.5,0.9]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> pos = 1
julia> d = chebyshev_derivative(weights,x,pos,ord,dom)
0.08487312307427555
```
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
Computes the first derivative of a complete Chebyshev polynomial with respect to the variable in position ```pos```, at point ```x```, given the ```weights```, the polynomial ```order```, and the ```domain```.

Signature
=========

d = chebyshev_derivative(weights,x,pos,order,domain)

Example
=======
```
julia> weights = [1.66416      0.598458    -0.0237052   0.00272941
  0.26378      0.0948592   -0.00375743  0.0
 -0.0246176   -0.00885287   0.0         0.0
  0.00372291   0.0          0.0         0.0]
julia> x = [5.5,0.9]
julia> ord = 3
julia> dom = [9.0 1.5; 3.0 0.5]
julia> pos = 1
julia> d = chebyshev_derivative(weights,x,pos,ord,dom)
0.08452286208888889
```
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

  derivative = zero(R)
  @inbounds for i in CartesianIndices(weights)
    if sum(i.I) <= order + N
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
Computes the gradient of a tensor-product Chebyshev polynomial evaluated at ```x````, given the ```weights```, the polynomial ```order```, and the ```domain```.

Signature
=========

g = chebyshev_gradient(weights,x,order,domain)

Example
=======
```
julia> weights = [1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]
julia> x = [5.5,0.9]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> g = chebyshev_gradient(weights,x,ord,dom)
[0.0848731  1.2072]
```
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
Computes the gradient of a complete Chebyshev polynomial evaluated at ```x````, given the ```weights```, the polynomial ```order```, and the ```domain```.

Signature
=========

g = chebyshev_gradient(weights,x,order,domain)

Example
=======
```
julia> weights = [1.66416      0.598458    -0.0237052   0.00272941
  0.26378      0.0948592   -0.00375743  0.0
 -0.0246176   -0.00885287   0.0         0.0
  0.00372291   0.0          0.0         0.0]
julia> x = [5.5,0.9]
julia> ord = 3
julia> dom = [9.0 1.5; 3.0 0.5]
julia> g = chebyshev_gradient(weights,x,ord,dom)
[0.0845229  1.20519]
```
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
Creates a function that evaluates the gradient of a Chebyshev polynomial evaluated at ```x```, given the approximation sample ```y``` and an approximation ```plan```.

Signature
=========

g = chebyshev_gradient(y,plan)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_nodes,[1.5,0.5])
julia> g = Grid((nodes1,nodes2))
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> plan = CApproxPlan(g,ord,dom)
julia> grad = chebyshev_gradient(y,plan)
julia> grad([5.5,0.9])
[0.0848731  1.2072]
```
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
Creates a function that evaluates using multi-threading the gradient of a Chebyshev polynomial evaluated at ```x```, given the approximation sample ```y``` and an approximation ```plan```.

Signature
=========

g = chebyshev_gradient_threaded(y,plan)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_nodes,[1.5,0.5])
julia> g = Grid((nodes1,nodes2))
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> plan = CApproxPlan(g,ord,dom)
julia> grad = chebyshev_gradient_threaded(y,plan)
julia> grad([5.5,0.9])
[0.0848731  1.2072]
```
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
Computes the hessian of a tensor-product Chebyshev polynomial evaluated at ```x````, given the ```weights```, the polynomial ```order```, and the ```domain```.

Signature
=========

h = chebyshev_hessian(weights,x,order,domain)

Example
=======
```
julia> weights = [1.66416      0.598458    -0.0237052     0.00272941
  0.26378      0.0948592   -0.00375743    0.000432628
 -0.0246176   -0.00885287   0.000350667  -4.03756e-5
  0.00372291   0.00133882  -5.30312e-5    6.10598e-6]
julia> x = [5.5,0.9]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> h = chebyshev_hessian(weights,x,ord,dom)
[-0.0118667   0.0661025
  0.0661025  -0.42678]
```
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

""""
Computes the hessian of a complete Chebyshev polynomial evaluated at ```x````, given the ```weights```, the polynomial ```order```, and the ```domain```.

Signature
=========

h = chebyshev_hessian(weights,x,order,domain)

Example
=======
```
julia> weights = [1.66416      0.598458    -0.0237052   0.00272941
  0.26378      0.0948592   -0.00375743  0.0
 -0.0246176   -0.00885287   0.0         0.0
  0.00372291   0.0          0.0         0.0]
julia> x = [5.5,0.9]
julia> ord = 3
julia> dom = [9.0 1.5; 3.0 0.5]
julia> h = chebyshev_hessian(weights,x,ord,dom)
[-0.0118089   0.069178
  0.069178   -0.421668]
```
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
Creates a function that evaluates the hessian of a Chebyshev polynomial evaluated at ```x```, given the approximation sample ```y``` and an approximation ```plan```.

Signature
=========

h = chebyshev_hessian(y,plan)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_nodes,[1.5,0.5])
julia> g = Grid((nodes1,nodes2))
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> plan = CApproxPlan(g,ord,dom)
julia> h = chebyshev_hessian(y,plan)
julia> h([5.5,0.9])
[-0.0118667   0.0661025
  0.0661025  -0.42678]
```
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
Creates a function that evaluates using multi-threading the hessian of a Chebyshev polynomial evaluated at ```x```, given the approximation sample ```y``` and an approximation ```plan```.

Signature
=========

h = chebyshev_hessian_threaded(y,plan)

Example
=======
```
julia> nodes1 = nodes(5,:chebyshev_nodes,[9.0,3.0])
julia> nodes2 = nodes(5,:chebyshev_nodes,[1.5,0.5])
julia> g = Grid((nodes1,nodes2))
julia> y = [nodes1.points[i]^0.3*nodes2.points[j]^0.7 for i in 1:5, j in 1:5]
julia> ord = (3,3)
julia> dom = [9.0 1.5; 3.0 0.5]
julia> plan = CApproxPlan(g,ord,dom)
julia> h = chebyshev_hessian_threaded(y,plan)
julia> h([5.5,0.9])
[-0.0118667   0.0661025
  0.0661025  -0.42678]
```
"""
function chebyshev_hessian_threaded(y::Array{T,N}, plan::P) where {T<:AbstractFloat,P<:CApproximationPlan,N}

  w = chebyshev_weights_threaded(y, plan)

  function cheb_hess(x::AbstractArray{R,1}) where {R<:Number}

    hess = chebyshev_hessian(w, x, plan.order, plan.domain)

    return hess

  end

  return cheb_hess

end