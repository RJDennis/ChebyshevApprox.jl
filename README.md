**ChebyshevApprox**
===============

ChebyshevApprox is a Julia package for approximating continuous functions using Chebyshev polynomials.  The package's focus is on multivariate functions that depend on an arbitrary number of variables.  Both tensor-product polynomials and complete polynomials are implemented.  Working with complete polynomials often leads to a considerable decrease in computation time with little loss of accuracy.  The package allows the nodes to be either the roots of the Chebyshev polynomial (points of the first kind), the extrema of the Chebyshev polynomial (points of the second kind), Chebyshev extended points (Chebyshev roots normalized so that the boundry nodes equal -1.0 and 1.0), the Vertesi nodes, or the Legendre nodes.  In addition to approximating functions the package also uses the approximating polynomial to compute derivatives, gradients, and hessians.

Installation
------------

ChebyshevApprox.jl is a registered package.  To install it simply type in the REPL:

```julia
using Pkg
Pkg.add("ChebyshevApprox")
```

Nodes
-----

The package contains functions for computing servel different types of approximating points.  Which points you use may depend on your application.

To compute the Chebyshev roots within the [1.0, -1.0] interval use:

```julia
n = 11
points = nodes(n,:chebyshev_nodes)
```

where `n`, an integer, is the number of nodes and `:chebyshev_nodes` is a symbol indicated the type of nodes to be produced.  Alternatives to `:chebyshev_nodes` are: `chebyshev_extrema`, `:chebyshev_extended`, `:vertesi_nodes`, and `:legendre_nodes`. 

To compute nodes over bounded domains other than the [1.0,-1.0] interval, the `nodes` function accepts an optional third argument containing the domain in the form of a 1D array (a vector) containing two elements, where the first element is the upper bound on the interval and the second is the lower bound.  For example,

```julia
domain = [3.5,0.5]
points = nodes(n,:chebyshev_extrema,domain)
```

would compute `n` extrema of the Chebyshev polynomial and scale those points to the [3.5,0.5] interval.

Grids
-----

Once the approximating points have been computed, these points can be collected together in a `Grid` structure:

```julia
p1 = nodes(11,:chebyshev_nodes)
p2 = nodes(15,:chebyshev_nodes)

g = Grid((p1,p2))
```

Approximation Plans
-------------------

The approximate a function, the package makes use of an approximation plan, which is a stucture containing the information needed to produce an approximating function.

```julia
A_plan = ApproxPlan(g,order,domain)
```

where `g` is a Grid, `order` is an integer (complete polynomial approximation) or a tuple (tensor-product approximation), and `domain` is an array containing the upper and lower limits for each variable.  For example,

```julia
dom_1 = [5.0,3.0]
dom_2 = [1.5,-0.5]
dom = [dom_1 dom_2]

p1 = nodes(11,:chebyshev_nodes,dom_1)
p2 = nodes(15,:chebyshev_nodes,dom_2)

g = Grid((p1,p2))
order = (6,6)
A_plan = ApproxPlan(g,order,dom)
```

Function approximation
----------------------

To approximate a function, you use the `chebyshev_interp` function, which itself returns a function.  If the data on the function you wish to approximate, sampled on the Grid, `g`, are contained in the Array, `y`, then the approximating function is generated using:

```julia
f_approx = chebyshev_interp(y,A_plan)
```

which can then be evaluated at a point, `x`, in the domain according to:

```julia
x = [2.0,0.1]
y_hat = f_approx(x)
```

Functions to approximate gradients and hessians are produced similarly:

```julia
f_grad = chebyshev_gradient(y,A_plan)
f_hess = chebyshev_hessian(y,A_plan)
grad_hat = f_grad(x)
hess_hat = f_hess(x)
```

There are multi-threaded versions of these functions: `chebyshev_interp_threaded()`, `chebyshev_gradient_threaded()`, and `chebyshev_hessian_threaded()`.

Under the hood
--------------

The functions documented above should cover a lot of use cases.  However, the machinery used to produce the chebyshev polynomials, find the polynomial weights, and evaluate polynomials at given points---which took center stage in releases before version 0.3.0---is still there and can be used as before.  As such this 0.3.x release is (mostly) not expected to be code-breaking (some stuctures have gone, but I don't think they were used much anyway).

Polynomials
-----------

Chebyshev polynomials are constructed using the `chebyshev_polynomial()` function, which takes two arguments.  The first argument is an integer representing the order of the polynomial.  The second argument is the point in the [1.0,-1.0] interval at which the polynominal is evaluated.  This second argument can be a scalar or a 1D array.  For example,

```julia
order = 5
x = 0.5
p = chebyshev_polynomial(order,x)
```

will return a 2D array containing the Chebyshev polynomials of orders 0---5 evaluated at the point `x`.  If `x` is a 1D array of points, as in:

```julia
order = 5
x = chebyshev_nodes(11)
p = chebyshev_polynomial(order,x)
```

then `p` will be a 2D array (11 $\times$ 6) containing the Chebyshev polynomials of orders 0---5 evaluated at each element in `x`.

Weights
-------

ChebyshevApprox.jl uses Chebyshev regression to compute the weights in the Chebyshev polynomial.  The central function for computing Chebyshev weights is the following:

```julia
w = chebyshev_weights(y,nodes,order,domain)
```

where `y` is a n-D array containing the function evaluations at `nodes`, `nodes` is a tuple of 1D arrays containing Chebyshev-roots, `order` is a tuple (tensor-product polynomial) or an integer (complete polynomial) specifying the order of the polynomial in each dimension, and `domain` is a 2D array containing the upper and lower bounds on the approximating interval in each dimension.  So,

```julia
order_x1  = 5
nodes_x1  = chebyshev_nodes(11)
domain_x1 = [3.5,0.5]

order_x2  = 7
nodes_x2  = chebyshev_nodes(15)
domain_x2 = [1.7,-0.3]

order  = (order_x1,order_x2)
nodes  = (nodes_x1,nodes_x2)
domain = [domain_x1 domain_x2]

w = chebyshev_weights(y,nodes,order,domain)
```

would compute the weights, `w`, (a 2D array in this example) in a tensor-product polynomial.  The domain-argument is optional, needed only if one or more variable does not have domain [1.0,-1.0].  The nodes-argument can be an array-of-arrays (instead of a tuple).  Alternatively, the polynominals can be computed and entered directly into the `chebyshev_weights()` function:

```julia
poly_1 = chebyshev_polynomial(order_x1,nodes_x1)
poly_2 = chebyshev_polynomial(order_x2,nodes_x2)
poly   = (poly_1,poly_2)

w = chebyshev_weights(y,poly,order)
```

The `poly`-argument can be an array-of-arrays (instead of a tuple).  The weights, `w` are returned in a (multi-dimensional) array.

If the solution nodes are instead the Chebyshev-extrema, then the analogue to the above is the use the chebyshev_weights_extrema() function.  For example:

```julia
order_x1  = 5
nodes_x1  = chebyshev_extrema(11)
domain_x1 = [3.5,0.5]

order_x2  = 7
nodes_x2  = chebyshev_extrema(15)
domain_x2 = [1.7,-0.3]

order  = [order_x1,order_x2]
nodes  = (nodes_x1,nodes_x2)
domain = [domain_x1 domain_x2]

w = chebyshev_weights_extrema(y,nodes,order,domain)
```

Other possibilities are to use `chebyshev_weights_extended()`, `chebyshev_weights_vertesi()`, or `chebyshev_weights_legendre()`.

Function evaluation
-------------------

ChebyshevApprox uses the `chebyshev_evaluate()` function, which accommodates several methods, to evaluate Chebyshev polynomials.  If `x` is a 1D array representing the point at which the polynomial is to be evaluated, then:

```julia
yhat = chebyshev_evaluate(w,x,order,domain)
```

For the case where a complete polynomial rather than a tensor-product polynomial is to be evaluated, the `order` variable is now simply an integer rather than a tuple of integers.

Derivatives
-----------

The `chebyshev_derivative()` function can be used to approximate the partial derivative of a function with respect to a designated variable.  For example, the partial derivative with respect to the 3'rd variable evaluated at point `x` can be computed by:

```julia
deriv = chebyshev_derivative(w,x,3,order,domain)
```

where `deriv` is a floating point number.

Gradients
---------

Gradients are computed using the `chebyshev_gradient()` function.

```julia
grad = chebyshev_gradient(w,x,order,domain)
```

where `grad` is a 2D array with one row.

Multi-threading
---------------

Computing the weights in a multivariate Chebyshev polynomial can be time-consuming for functions whose dimensions are large, or where the number of nodes and/or the order of the polynomals is large.  For this reason, multi-threaded functions for computing the weights are provided.  For example, if the nodes are Chebyshev-roots:

```julia
w = chebyshev_weights_threaded(y,nodes,order,domain)
```

etc.

For completeness, there are also:

```julia
w = chebyshev_weights_extrema_threaded(y,nodes,order,domain)
w = chebyshev_weights_extended_threaded(y,nodes,order,domain)
w = chebyshev_weights_vertesi_threaded(y,nodes,order,domain)
w = chebyshev_weights_legendre_threaded(y,nodes,order,domain)
```

Related packages
----------------

If you are looking to approximate a function of one variable, then there is:

- ApproxFun.jl

For multivariate functions, there is:

- SmolyakApprox.jl
- HyperbolicCrossApprox.jl
- PiecewiseLinearApprox.jl

Finally, if you are approximating complex-valued multivariate functions, then a package insipired by ChebyshevApprox.jl and offering similar functionality is:

- FastChebInterp.jl
