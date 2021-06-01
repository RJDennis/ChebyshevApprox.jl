ChebyshevApprox
===============

ChebyshevApprox is a Julia package for approximating continuous functions using Chebyshev polynomials.  The package's focus is on multivariate functions that depend on an arbitrary number of variables.  Both tensor-product polynomials and complete polynomials are implemented.  Working with complete polynomials often leads to a considerable decrease in computation time with little loss of accuracy.  The package allows the nodes to be either the roots of the Chebyshev polynomial (points of the first kind), the extrema of the Chebyshev polynomial (points of the second kind), Chebyshev extended points (Chebyshev roots normalized so that the boundry nodes equal -1.0 and 1.0), or the Vertesi nodes.  In addition to approximating functions the package also uses the approximating polynomial to compute derivatives and gradients.

Installation
------------

ChebyshevApprox is a registered package.  To install it simply type in the REPL:

```
using Pkg
Pkg.add("ChebyshevApprox")
````

Nodes
-----

The package contains functions for computing both the roots of the Chebyshev polynomial and the extrema of the Chebyshev polynominal.  Depending of the application, you may wish to use one or the other.

To compute the Chebyshev roots within the [1.0, -1.0] interval use:

```
nodes = chebyshev_nodes(n)
```

where `n`, an integer, is the number of nodes.  Similarly, to compute the Chebyshev extrema within the [1.0,-1.0] interval use:

```
nodes = chebyshev_extrema(n)
```

The extended nodes and the Vertesi nodes over the [1.0,-1.0] interval are obtained from:

```
nodes = chebyshev_extended(n)
nodes = vertesi_nodes(n)
```

respectively.

To compute nodes over bounded domains other than the [1.0,-1.0] interval, all three functions accept a second argument containing the domain in the form of a 1D array (a vector) containing two elements, where the first element is the upper bound on the interval and the second is the lower bound.  For example,

```
domain = [3.5,0.5]
nodes = chebyshev_nodes(n,domain)
```

would compute `n` roots of the Chebyshev polynomial and scale those roots to the [3.5,0.5] interval.

Polynomials
-----------

Chebyshev polynomials are constructed using the chebyshev_polynomial() function, which takes two arguments.  The first argument is an integer representing the order of the polynomial.  The second argument is the point in the [1.0,-1.0] interval at which the polynominal is evaluated.  This second argument can be a scalar or a 1D array.  For example,

```
order = 5
x = 0.5
p = chebyshev_polynomial(order,x)
```

will return a 2D array containing the Chebyshev polynomials of orders 0---5 evaluated at the point `x`.  If `x` is a 1D array of points, as in:

```
order = 5
x = chebyshev_nodes(11)
p = chebyshev_polynomial(order,x)
```

then `p` will be a 2D array (11*6) containing the Chebyshev polynomials of orders 0---5 evaluated at each element in `x`.

Structures
----------

ChebyshevApprox contains five structures that can make your life easier.  The first contains the information needed to evaluate a polynomial at a point.  I.e.:

```
chebpoly = ChebPoly(w,order,domain)
```

where `order` would be an integer (complete polynomial) or a 1D array of integers (tensor-product polynomial).

The remaining four structures are interpolation objects, which are created as follows:

```
cheby = ChebInterpRoots(y,nodes,order,domain)
cheby = ChebInterpExtrema(y,nodes,order,domain)
cheby = ChebInterpExtended(y,nodes,order,domain)
cheby = ChebInterpVertesi(y,nodes,order,domain)
```

where `y` is an n-D array containing function-values evaluated at the nodes, `nodes` is a tuple, `order` would be an integer or a 1D array of integers, and `nodes` would be Chebyshev-roots in the first, Chebyshev-extrema in the second, Chebyshev-extended-roots in the third, and Vertesi nodes in the fourth.

Weights
-------

We focus here on the case where the solution nodes are Chebyshev-roots and cover the cases where they are Chebyshev-extrema, Chebyshev-extended-roots, and Vertesi-nodes subsequently.

ChebyshevApprox uses Chebyshev regression to compute the weights in the Chebyshev polynomial.  The central function for computing Chebyshev weights is the following:

```
w = chebyshev_weights(y,nodes,order,domain)
```

where `y` is a n-D array containing the function evaluations at `nodes`, `nodes` is a tuple of 1D arrays containing Chebyshev-roots, 'order' is a 1D array (tensor-product polynomial) or an integer (complete polynomial) specifying the order of the polynomial in each dimension, and `domain` is a 2D array containing the upper and lower bounds on the approximating interval in each dimension.  So,

```
order_x1  = 5
nodes_x1  = chebyshev_nodes(11)
domain_x1 = [3.5,0.5]

order_x2  = 7
nodes_x2  = chebyshev_nodes(15)
domain_x2 = [1.7,-0.3]

order  = [order_x1,order_x2]
nodes  = (nodes_x1,nodes_x2)
domain = [domain_x1 domain_x2]

w = chebyshev_weights(y,nodes,order,domain)
```

would compute the weights, `w`, (a 2D array in this example) in a tensor-product polynomial.  The domain-argument is optional, needed only if one or more variable does not have domain [1.0,-1.0].  The nodes-argument can be an array-of-arrays (instead of a tuple).  Alternatively, the polynominals can be computed and entered directly into the chebyshev_weights() function:

```
p1   = chebyshev_polynomial(order_x1,nodes_x1)
p2   = chebyshev_polynomial(order_x2,nodes_x2)
poly = (p1,p2)

w = chebyshev_weights(y,poly,order)
```

The poly-argument can be an array-of-arrays (instead of a tuple).  Further, using the ChebInterp structure:

```
w = chebyshev_weights(cheb)
```

For all of these functions the `weights` are returned in a (multi-dimensional) array.

If the solution nodes are instead the Chebyshev-extrema, then the analogue to the above is the use the chebyshev_weights_extrema() function.  For example,

```
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

If the solution nodes are the Chebyshev-extended-roots, then use the chebyshev_weights_extended() function, as per:

```
order_x1  = 5
nodes_x1  = chebyshev_extended(11)
domain_x1 = [3.5,0.5]

order_x2  = 7
nodes_x2  = chebyshev_extended(15)
domain_x2 = [1.7,-0.3]

order  = [order_x1,order_x2]
nodes  = (nodes_x1,nodes_x2)
domain = [domain_x1 domain_x2]

w = chebyshev_weights_extended(y,nodes,order,domain)
```

Finally, if the solution nodes are the Vertesi-nodes, then use the chebyshev_weights_vertesi() function:

```
order_x1  = 5
nodes_x1  = vertesi_nodes(11)
domain_x1 = [3.5,0.5]

order_x2  = 7
nodes_x2  = vertesi_nodes(15)
domain_x2 = [1.7,-0.3]

order  = [order_x1,order_x2]
nodes  = (nodes_x1,nodes_x2)
domain = [domain_x1 domain_x2]

w = chebyshev_weights_vertesi(y,nodes,order,domain)
```

Function evaluation
-------------------

ChebyshevApprox uses the chebyshev_evaluate() function, which accommodates several methods, to evaluate Chebyshev polynomials.  If `x` is a 1D array representing the point at which the polynomial is to be evaluated, then:

```
yhat = chebyshev_evaluate(w,x,order,domain)
```

or

```
yhat = chebyshev_evaluate(chebpoly,x)
```

or

```
yhat = chebyshev_evaluate(cheby,x)
```

are equivalent.  For the case where a complete polynomial rather than a tensor-product polynomial is to be evaluated, the commands are the same as above, but the `order` variable is now simply an integer rather than a 1D array of integers.

ChebyshevApprox also allows the following:

```
cheb = chebyshev_evaluate(w,order,domain)
```

or

```
cheb = chebyshev_evaluate(chebpoly)
```

or

```
cheb = cheb_interp(cheby)
```

followed by

```
yhat = cheb(x)
```

allowing polynomials to be easily evaluated at point `x`.

Derivatives
-----------

The chebyshev_derivative() function can be used to approximate the partial derivative of a function with respect to a designated variable.  For example, the partial derivative with respect to the 3'rd variable evaluated at point `x` can be computed by:

```
deriv = chebyshev_derivative(w,x,3,order,domain)
```

or

```
deriv = chebyshev_derivative(chebpoly,x,3)
```

or

```
deriv = chebyshev_derivative(cheby,x,3)
```

where `deriv` that is returned is a floating point number.

Gradients
---------

Gradients are computed using the chebyshev_gradient() function.

```
grad = chebyshev_gradient(w,x,order,domain)
```

or

```
grad = chebyshev_gradient(chebpoly,x)
```

or

```
grad = chebyshev_gradient(cheby,x)
```

where `grad` that is returned is a 2D array with one row.

Multi-threading
---------------

Computing the weights in a multivariate Chebyshev polynomial can be time-consuming for functions whose dimensions are large, or where the number of nodes and/or the order of the polynomals is large.  For this reason, multi-threaded functions for computing the weights are provided.  If the nodes are Chebyshev-roots:

```
w = chebyshev_weights_threaded(y,nodes,order,domain)
```

and

```
w = chebyshev_weights_threaded(y,poly,order)
```

and

```
w = chebyshev_weights_threaded(cheby)
```

As earlier, these functions can be used to compute weights for either tensor-product polynomials or complete polynomials.  Threaded versions are also provided for the cases where the nodes are Chebyshev-extrema, Chebyshev-extended-roots, and Chebyshev-Vertesi-points, such as:

```
w = chebyshev_weights_extrema_threaded(y,nodes,order,domain)
w = chebyshev_weights_extended_threaded(y,nodes,order,domain)
w = chebyshev_weights_vertesi_threaded(y,nodes,order,domain)
```

What if I want to use my own set of nodes?
------------------------------------------

You can.  The routines used for Chebyshev-extended-nodes and Vertesi-nodes are general, so applying the functions for these node-cases to your own set of nodes should work.