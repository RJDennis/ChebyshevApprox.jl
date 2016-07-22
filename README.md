ChebyshevApprox
===============

Julia code to approximate continuous functions using Chebyshev polynomials.

In many applications it is important to be able to approximate accurately continuous functions of several variables.  This code can approximate functions that depend upon as many as six variables.

Chebyshev nodes in the [1.0, -1.0] interval are computed using the command

nodes = chebyshev_nodes(n)

where n, an integer, is the number of nodes.  Similarly, Chebyshev nodes scaled to a bounded interval are computed using the command

nodes = chebyshev_nodes(n,range)

where range is a 1-d array with two elements, the first representing the upper bound on the interval and the second the lower bound.

Approximating functions of one variable
=======================================

To approximate a 1-d function using a Chebyshev polynomial one first computes the Chebyshev weights.  There are several ways of doing this.  The central command is

w = chebyshev_weights(y,nodes,order,range)

where y is a 1-d array containing function evaluations at nodes and order is an integer (or 1-d integer array containing one element) representing the order of the polynomial.  nodes can be either a 1-d array or an array containing a 1-d array.  An alternative is the command

w = chebyshev_weights(y,poly,order,range)

where poly is an array containing a 2-d array.  Regardless of which command is used, the Chebyshev weights are computed using Chebyshev regression.

Given the Chebyshev weights, the function can be approximated at a given point, x, within range using either

y_approx = chebyshev_evaluate(w,[x],order,range)

or

y_approx = clenshaw_evaluate(w,[x],order,range)

The latter command evaluates the polynominal using Clenshaw's recursion.  For functions where the domain of x coincides with the [1.0, -1.0] interval the range argument can be omitted, both when computing the weights and when evaluating the polynomial.

Approximating functions of several variables
============================================

The codes can approximate functions of an arbitrary number of variables.  The only real difference is that in the higher dimension cases the codes allow the polynomials to be either tensor products of 1-d polynomials or complete polynomials.  When the function depends on six or fewer variables, and in the tensor-product case, to approximate a function (of three variables, say) the relevant commands are

w = chebyshev_weights(y,nodes_1,nodes_2,nodes_3,order,range)      --- (1)

where order is a 1-d integer array whose elements represents the polynomial orders in each successive dimension, or

w = chebyshev_weights(y,nodes,order,range)                        --- (2)

where nodes is an array containing 1-d arrays, or

w = chebyshev_weights(y,poly,order,range)                         --- (3)

where poly is an array containing 2-d arrays.  To go beyond functions of six variables, commands (2) or (3) can be used.

Then to evaluate the polynomial, the relevant commands are either

y_approx = chebyshev_evaluate(w,x,order,range)

or

y_approx = clenshaw_evaluate(w,x,order,range)

In each of these commands, the variable range is a 2-d array (2*3 matrix for the 3-variable case) whose columns summarize the upper- and lower-bound on each variable.

For the case where a complete polynomial rather than a tensor-product polynomial is desired, the commands are the same as above, but the order variable is now simply an integer rather than a 1-d array of integers.

As with the 1-d case, if the domain of each variable is [1.0, -1.0], then the range variable can be omitted.

Working with complete polynomials rather than tensor-product polynomials often leads to a considerable decrease in computation time with little loss of accuracy.

Summary
=======

Have fun.
