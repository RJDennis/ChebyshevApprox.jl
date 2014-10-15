ChebyshevApprox
===============

Julia code to approximate continuous functions using Chebyshev polynomials.

In many applications it is important to be able to approximate accurately continuous functions of several variables.  This code can approximate functions that depend upon as many as five variables.  Chebyshev nodes in the [1.0, -1.0] interval are computed using the command

nodes = chebyshev_nodes(n)

where n is the number of nodes.  Similarly, Chebyshev nodes scaled to a bounded interval are computed using the command

nodes = chebyshev_nodes(n,range)

where range is a 1-d array with two elements, the first representing the upper bound on the interval and the second the lower bound.

Approximating functions of one variable
=======================================

To approximate a 1-d function using a Chebyshev polynomial one first computes the Chebyshev weights using the command

w = chebyshev_weights(y,nodes,order,range) 

where y is a 1-d array containing function evaluations at nodes and order is an integer (or 1-d integer array containing one element) representing the order of the polynomial.  The Chebyshev weights are computed using Chebyshev regression.

Given the Chebyshev weights, the function can be approximated at a given point, x, within range using either

y_approx = chebyshev_evaluate(w,[x],order,range)

or

y_approx = clenshaw_evaluate(w,[x],order,range)

The latter command evaluates the polynominal using Clenshaw's recursion.  For functions where the domain of x coincides with the [1.0, -1.0] interval the range argument can be omitted. 

Approximating functions of several variables
============================================

The codes can approximate functions of as many as five variables and the extension from the 1-d case to the higher-dimension cases is straight-forward.  The only real difference is that in the higher dimension cases the codes allow the polynomials to be either tensor products of 1-d polynomials or complete polynomials.  In the tensor-product case, to approximate a function of three variables (say) the relevant commands are

w = chebyshev_weights(y,nodes_1,nodes_2,nodes_3,order,range)

where order is a 1-d integer array whose elements represents the polynomial orders in each successive dimension, and either

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