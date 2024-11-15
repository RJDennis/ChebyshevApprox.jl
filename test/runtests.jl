using ChebyshevApprox
using LinearAlgebra
using GenericLinearAlgebra

#######################
### Testing on Float64
#######################

##################################################
# Testing the construction of approximating points
##################################################

# Chebyshev nodes

r1 = chebyshev_nodes(3)
r2 = chebyshev_nodes(3,[2.0,0.0])
r3 = nodes(3,:chebyshev_nodes)
r4 = nodes(3,:chebyshev_nodes,[2.0,0.0])

r5 = chebyshev_nodes(BigFloat,3)
r6 = chebyshev_nodes(BigFloat,3,[2.0,0.0])
r7 = nodes(BigFloat,3,:chebyshev_nodes)
r8 = nodes(BigFloat,3,:chebyshev_nodes,[2.0,0.0])

# Chebyshev extrema

e1 = chebyshev_extrema(3)
e2 = chebyshev_extrema(3,[2.0,0.0])
e3 = nodes(3,:chebyshev_extrema)
e4 = nodes(3,:chebyshev_extrema,[2.0,0.0])

e5 = chebyshev_extrema(BigFloat,3)
e6 = chebyshev_extrema(BigFloat,3,[2.0,0.0])
e7 = nodes(BigFloat,3,:chebyshev_extrema)
e8 = nodes(BigFloat,3,:chebyshev_extrema,[2.0,0.0])

# Chebyshev extended points

ex1 = chebyshev_extended(3)
ex2 = chebyshev_extended(3,[2.0,0.0])
ex3 = nodes(3,:chebyshev_extended)
ex4 = nodes(3,:chebyshev_extended,[2.0,0.0])

ex5 = chebyshev_extended(BigFloat,3)
ex6 = chebyshev_extended(BigFloat,3,[2.0,0.0])
ex7 = nodes(BigFloat,3,:chebyshev_extended)
ex8 = nodes(BigFloat,3,:chebyshev_extended,[2.0,0.0])

########################################################
# Test the construction of various Chebyshev polynomials
########################################################

r1 = chebyshev_nodes(11)
r2 = chebyshev_nodes(11,[6.5,3.4])
r3 = nodes(11,:chebyshev_nodes)

r5 = chebyshev_nodes(BigFloat,11)
r6 = chebyshev_nodes(BigFloat,11,[6.5,3.4])
r7 = nodes(BigFloat,11,:chebyshev_nodes)

# Chebyshev polynomials

chebyshev_polynomial(6, r1)
chebyshev_polynomial(6, r2, [6.5,3.4])
chebyshev_polynomial(6, r3)
chebyshev_polynomial(6, r5)
chebyshev_polynomial(6, r6, [6.5, 3.4])
chebyshev_polynomial(6, r7)

# Derivative of Chebyshev polynomials

chebyshev_polynomial_deriv(6, r1)
chebyshev_polynomial_deriv(6, r3)
chebyshev_polynomial_deriv(6, r5)
chebyshev_polynomial_deriv(6, r7)

# Second derivative of Chebyshev polynomials

chebyshev_polynomial_sec_deriv(6, r1)
chebyshev_polynomial_sec_deriv(6, r3)
chebyshev_polynomial_sec_deriv(6, r5)
chebyshev_polynomial_sec_deriv(6, r7)

#################################
# Test the construction of a Grid
#################################

grid1 = Grid((r3,))
grid2 = Grid((r3,r3))

grid3 = Grid((r7,))
grid4 = Grid((r7,r7))

##############################################
# Test the construction of CApproximationPlans
##############################################

# CApproxPlan

plan1 = CApproxPlan(grid1,6,[1.0,-1.0])
plan2 = CApproxPlan(grid2,6,[1.0 1.0; -1.0 -1.0])
plan3 = CApproxPlan(grid2,(6,6),[1.0 1.0; -1.0 -1.0])

plan4 = CApproxPlan(grid3,6,[1.0,-1.0])
plan5 = CApproxPlan(grid4,6,[1.0 1.0; -1.0 -1.0])
plan6 = CApproxPlan(grid4,(6,6),[1.0 1.0; -1.0 -1.0])

# CApproxPlanPoly

poly1 = chebyshev_polynomial(6, r3)
poly2 = chebyshev_polynomial(6, r3)

b1  = CApproxPlanPoly((poly1,), 6, [1.0, -1.0])
b2  = CApproxPlanPoly((poly1,poly2), 6, [1.0 1.0; -1.0 -1.0])
b3  = CApproxPlanPoly((poly1,poly2), (6,6), [1.0 1.0; -1.0 -1.0])

poly3 = chebyshev_polynomial(6, r7)
poly4 = chebyshev_polynomial(6, r7)

b4  = CApproxPlanPoly((poly3,), 6, [1.0, -1.0])
b5  = CApproxPlanPoly((poly3,poly4), 6, [1.0 1.0; -1.0 -1.0])
b6  = CApproxPlanPoly((poly3,poly4), (6,6), [1.0 1.0; -1.0 -1.0])

##############################################
# Testing the computation of chebyshev weights
##############################################

###################################
# Test the 1-variable case, Float64
###################################

n = 10
order = 5
dom = [2.0, -3.0]

nodesr = chebyshev_nodes(n, dom)
nodese = chebyshev_extrema(n, dom)
nodesex = chebyshev_extended(n, dom)

gr = nodes(n, :chebyshev_nodes, dom)
ge = nodes(n, :chebyshev_extrema, dom)
gex = nodes(n, :chebyshev_extended, dom)

yr  = [sqrt(nodesr[i]+4) for i in 1:n]
ye  = [sqrt(nodese[i]+4) for i in 1:n]
yex = [sqrt(nodesex[i]+4) for i in 1:n]

wr_tensor = chebyshev_weights(yr, nodesr, [order], dom)
wr_complete = chebyshev_weights(yr, nodesr, order, dom)
we_tensor = chebyshev_weights_extrema(ye, nodese, [order], dom)
we_complete = chebyshev_weights_extrema(ye, nodese, order, dom)
wex_tensor = chebyshev_weights_extended(yex, nodesex, [order], dom)
wex_complete = chebyshev_weights_extended(yex, nodesex, order, dom)

hr = Grid((gr,))
he = Grid((ge,))
hex = Grid((gex,))

ar = CApproxPlan(hr, order, dom)
ae = CApproxPlan(he, order, dom)
aex = CApproxPlan(hex, order, dom)

pr = chebyshev_polynomial(order, gr)
pe = chebyshev_polynomial(order, ge)
pex = chebyshev_polynomial(order, gex)

br = CApproxPlanPoly((pr,), order, dom)
be = CApproxPlanPoly((pe,), order, dom)
bex = CApproxPlanPoly((pex,), order, dom)

wr = chebyshev_weights(yr, ar)
we = chebyshev_weights(ye, ae)
wex = chebyshev_weights(yex, aex)

wr_t = chebyshev_weights_threaded(yr, ar)
we_t = chebyshev_weights_threaded(ye, ae)
wex_t = chebyshev_weights_threaded(yex, aex)

wr = chebyshev_weights(yr, br)
we = chebyshev_weights(ye, be)
wex = chebyshev_weights(yex, bex)

wr_t = chebyshev_weights_threaded(yr, br)
we_t = chebyshev_weights_threaded(ye, be)
wex_t = chebyshev_weights_threaded(yex, bex)

point = [1.44]

yr_hat = chebyshev_evaluate(wr_tensor, point, order, dom)
ye_hat = chebyshev_evaluate(we_tensor, point, order, dom)
yex_hat = chebyshev_evaluate(wex_tensor, point, order, dom)

yr_fn = chebyshev_interp(yr, ar)
ye_fn = chebyshev_interp(ye, ae)
yex_fn = chebyshev_interp(yex, aex)

yr_fn(point)
ye_fn(point)
yex_fn(point)

yr_fn = chebyshev_interp_threaded(yr, ar)
ye_fn = chebyshev_interp_threaded(ye, ae)
yex_fn = chebyshev_interp_threaded(yex, aex)

yr_fn(point)
ye_fn(point)
yex_fn(point)

yr_grad = chebyshev_gradient(wr_tensor, point, order, dom)
ye_grad = chebyshev_gradient(we_tensor, point, order, dom)
yex_grad = chebyshev_gradient(wex_tensor, point, order, dom)

yr_grad = chebyshev_gradient(yr, ar)
ye_grad = chebyshev_gradient(ye, ae)
yex_grad = chebyshev_gradient(yex, aex)

yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_grad = chebyshev_gradient_threaded(yr, ar)
ye_grad = chebyshev_gradient_threaded(ye, ae)
yex_grad = chebyshev_gradient_threaded(yex, aex)

yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_grad = chebyshev_gradient(yr, br)
ye_grad = chebyshev_gradient(ye, be)
yex_grad = chebyshev_gradient(yex, bex)

yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_grad = chebyshev_gradient_threaded(yr, br)
ye_grad = chebyshev_gradient_threaded(ye, be)
yex_grad = chebyshev_gradient_threaded(yex, bex)

yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_hess = chebyshev_hessian(yr, ar)
ye_hess = chebyshev_hessian(ye, ae)
yex_hess = chebyshev_hessian(yex, aex)

yr_hess(point)
ye_hess(point)
yex_hess(point)

yr_hess = chebyshev_hessian_threaded(yr, ar)
ye_hess = chebyshev_hessian_threaded(ye, ae)
yex_hess = chebyshev_hessian_threaded(yex, aex)

yr_hess(point)
ye_hess(point)
yex_hess(point)

yr_hess = chebyshev_hessian(yr, br)
ye_hess = chebyshev_hessian(ye, be)
yex_hess = chebyshev_hessian(yex, bex)

yr_hess(point)
ye_hess(point)
yex_hess(point)

yr_hess = chebyshev_hessian_threaded(yr, br)
ye_hess = chebyshev_hessian_threaded(ye, be)
yex_hess = chebyshev_hessian_threaded(yex, bex)

yr_hess(point)
ye_hess(point)
yex_hess(point)

###################################
# Test the 2-variable case, Float64
###################################

n1 = 51
n2 = 51

dom_1 = [2.0, -3.0]
dom_2 = [1.5, 0.5]
dom = [dom_1 dom_2]

nodesr_1 = chebyshev_nodes(n1, dom_1)
nodesr_2 = chebyshev_nodes(n2, dom_2)
nodese_1 = chebyshev_extrema(n1, dom_1)
nodese_2 = chebyshev_extrema(n2, dom_2)
nodesex_1 = chebyshev_extended(n1, dom_1)
nodesex_2 = chebyshev_extended(n2, dom_2)

gr_1 = nodes(n1, :chebyshev_nodes, dom_1)
gr_2 = nodes(n2, :chebyshev_nodes, dom_2)
ge_1 = nodes(n1, :chebyshev_extrema, dom_1)
ge_2 = nodes(n2, :chebyshev_extrema, dom_2)
gex_1 = nodes(n1, :chebyshev_extended, dom_1)
gex_2 = nodes(n2, :chebyshev_extended, dom_2)

yr  = [(nodesr_1[i] + 4.0)^0.5 + nodesr_1[i] * sqrt(nodesr_2[j]) for i in 1:n1, j in 1:n2]
ye  = [(nodese_1[i] + 4.0)^0.5 + nodese_1[i] * sqrt(nodese_2[j]) for i in 1:n1, j in 1:n2]
yex = [(nodesex_1[i] + 4.0)^0.5 + nodesex_1[i] * sqrt(nodesex_2[j]) for i in 1:n1, j in 1:n2]

order_tensor = [10, 10]
order_complete = 10

wr_tensor = chebyshev_weights(yr, [nodesr_1, nodesr_2], order_tensor, dom)
wr_complete = chebyshev_weights(yr, (nodesr_1, nodesr_2), order_complete, dom)

we_tensor = chebyshev_weights_extrema(ye, (nodese_1, nodese_2), order_tensor, dom)
we_complete = chebyshev_weights_extrema(ye, (nodese_1, nodese_2), order_complete, dom)

wex_tensor = chebyshev_weights_extended(yex, (nodesex_1, nodesex_2), order_tensor, dom)
wex_complete = chebyshev_weights_extended(yex, (nodesex_1, nodesex_2), order_complete, dom)

order = (10, 10)
hr = Grid((gr_1, gr_2))
he = Grid((ge_1, ge_2))
hex = Grid((gex_1, gex_2))

ar = CApproxPlan(hr, order, dom)
ae = CApproxPlan(he, order, dom)
aex = CApproxPlan(hex, order, dom)

pr_1 = chebyshev_polynomial(order[1], gr_1)
pr_2 = chebyshev_polynomial(order[2], gr_2)
pe_1 = chebyshev_polynomial(order[1], ge_1)
pe_2 = chebyshev_polynomial(order[2], ge_2)
pex_1 = chebyshev_polynomial(order[1], gex_1)
pex_2 = chebyshev_polynomial(order[2], gex_2)

br = CApproxPlanPoly((pr_1, pr_2), order, dom)
be = CApproxPlanPoly((pe_1, pe_2), order, dom)
bex = CApproxPlanPoly((pex_1, pex_2), order, dom)

wr = chebyshev_weights(yr, ar)
we = chebyshev_weights(ye, ae)
wex = chebyshev_weights(yex, aex)

wr_t = chebyshev_weights_threaded(yr, ar)
we_t = chebyshev_weights_threaded(ye, ae)
wex_t = chebyshev_weights_threaded(yex, aex)

wr = chebyshev_weights(yr, br)
we = chebyshev_weights(ye, be)
wex = chebyshev_weights(yex, bex)

wr_t = chebyshev_weights_threaded(yr, br)
we_t = chebyshev_weights_threaded(ye, be)
wex_t = chebyshev_weights_threaded(yex, bex)

point = [1.748, 0.753]

yr_hat = chebyshev_evaluate(wr_tensor, point, order, dom)
ye_hat = chebyshev_evaluate(we_tensor, point, order, dom)
yex_hat = chebyshev_evaluate(wex_tensor, point, order, dom)

yr_fn = chebyshev_interp(yr, ar)
ye_fn = chebyshev_interp(ye, ae)
yex_fn = chebyshev_interp(yex, aex)
yr_fn(point)
ye_fn(point)
yex_fn(point)

yr_fn = chebyshev_interp_threaded(yr, ar)
ye_fn = chebyshev_interp_threaded(ye, ae)
yex_fn = chebyshev_interp_threaded(yex, aex)
yr_fn(point)
ye_fn(point)
yex_fn(point)

yr_grad = chebyshev_gradient(wr_tensor, point, order, dom)
ye_grad = chebyshev_gradient(we_tensor, point, order, dom)
yex_grad = chebyshev_gradient(wex_tensor, point, order, dom)

yr_grad = chebyshev_gradient(yr, ar)
ye_grad = chebyshev_gradient(ye, ae)
yex_grad = chebyshev_gradient(yex, aex)
yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_grad = chebyshev_gradient_threaded(yr, ar)
ye_grad = chebyshev_gradient_threaded(ye, ae)
yex_grad = chebyshev_gradient_threaded(yex, aex)
yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_grad = chebyshev_gradient(yr, br)
ye_grad = chebyshev_gradient(ye, be)
yex_grad = chebyshev_gradient(yex, bex)
yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_grad = chebyshev_gradient_threaded(yr, br)
ye_grad = chebyshev_gradient_threaded(ye, be)
yex_grad = chebyshev_gradient_threaded(yex, bex)
yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_hess = chebyshev_hessian(yr, ar)
ye_hess = chebyshev_hessian(ye, ae)
yex_hess = chebyshev_hessian(yex, aex)
yr_hess(point)
ye_hess(point)
yex_hess(point)

yr_hess = chebyshev_hessian_threaded(yr, ar)
ye_hess = chebyshev_hessian_threaded(ye, ae)
yex_hess = chebyshev_hessian_threaded(yex, aex)
yr_hess(point)
ye_hess(point)
yex_hess(point)

yr_hess = chebyshev_hessian(yr, br)
ye_hess = chebyshev_hessian(ye, be)
yex_hess = chebyshev_hessian(yex, bex)
yr_hess(point)
ye_hess(point)
yex_hess(point)

yr_hess = chebyshev_hessian_threaded(yr, br)
ye_hess = chebyshev_hessian_threaded(ye, be)
yex_hess = chebyshev_hessian_threaded(yex, bex)
yr_hess(point)
ye_hess(point)
yex_hess(point)

########################
### TESTING ON BIGFLOATS
########################

# Test the construction of nodes

nodes(BigFloat, 111, :chebyshev_nodes)
nodes(BigFloat, 111, :chebyshev_extrema)
nodes(BigFloat, 111, :chebyshev_extended)

# Test the construction of polynomials

g = nodes(BigFloat, 1001, :chebyshev_nodes)

chebyshev_polynomial(6, g.points)
chebyshev_polynomial(6, g)
chebyshev_polynomial_deriv(6, g.points)
chebyshev_polynomial_deriv(6, g)
chebyshev_polynomial_sec_deriv(6, g.points)
chebyshev_polynomial_sec_deriv(6, g)

# Test the construction of a Grid

h = Grid((g, g))

# Test the construction of ApproximationPlans

ord = 5
dom = [1.0 1.0; -1.0 -1.0]
a = CApproxPlan(h, ord, dom)

p1 = chebyshev_polynomial(ord, g)
p2 = chebyshev_polynomial(ord, g)
b = CApproxPlanPoly((p1, p2), ord, dom)

####################################
# Test the 1-variable case, BigFloat
####################################

n = 10
order = 5
dom = [2.0, -3.0]

nodesr = chebyshev_nodes(BigFloat, n, dom)
nodese = chebyshev_extrema(BigFloat, n, dom)
nodesex = chebyshev_extended(BigFloat, n, dom)

gr = nodes(BigFloat, n, :chebyshev_nodes, dom)
ge = nodes(BigFloat, n, :chebyshev_extrema, dom)
gex = nodes(BigFloat, n, :chebyshev_extended, dom)

yr  = [sqrt(nodesr[i]+4) for i in 1:n]
ye  = [sqrt(nodese[i]+4) for i in 1:n]
yex = [sqrt(nodesex[i]+4) for i in 1:n]

wr_tensor = chebyshev_weights(yr, nodesr, [order], dom)
wr_complete = chebyshev_weights(yr, nodesr, order, dom)
we_tensor = chebyshev_weights_extrema(ye, nodese, [order], dom)
we_complete = chebyshev_weights_extrema(ye, nodese, order, dom)
wex_tensor = chebyshev_weights_extended(yex, nodesex, [order], dom)
wex_complete = chebyshev_weights_extended(yex, nodesex, order, dom)

hr = Grid((gr,))
he = Grid((ge,))
hex = Grid((gex,))

ar = CApproxPlan(hr, order, dom)
ae = CApproxPlan(he, order, dom)
aex = CApproxPlan(hex, order, dom)

pr = chebyshev_polynomial(order, gr)
pe = chebyshev_polynomial(order, ge)
pex = chebyshev_polynomial(order, gex)

br = CApproxPlanPoly((pr,), order, dom)
be = CApproxPlanPoly((pe,), order, dom)
bex = CApproxPlanPoly((pex,), order, dom)

wr = chebyshev_weights(yr, ar)
we = chebyshev_weights(ye, ae)
wex = chebyshev_weights(yex, aex)

wr_t = chebyshev_weights_threaded(yr, ar)
we_t = chebyshev_weights_threaded(ye, ae)
wex_t = chebyshev_weights_threaded(yex, aex)

wr = chebyshev_weights(yr, br)
we = chebyshev_weights(ye, be)
wex = chebyshev_weights(yex, bex)

wr_t = chebyshev_weights_threaded(yr, br)
we_t = chebyshev_weights_threaded(ye, be)
wex_t = chebyshev_weights_threaded(yex, bex)

point = [BigFloat(1.44)]

yr_hat = chebyshev_evaluate(wr_tensor, point, order, dom)
ye_hat = chebyshev_evaluate(we_tensor, point, order, dom)
yex_hat = chebyshev_evaluate(wex_tensor, point, order, dom)

yr_fn = chebyshev_interp(yr, ar)
ye_fn = chebyshev_interp(ye, ae)
yex_fn = chebyshev_interp(yex, aex)

yr_fn(point)
ye_fn(point)
yex_fn(point)

yr_fn = chebyshev_interp_threaded(yr, ar)
ye_fn = chebyshev_interp_threaded(ye, ae)
yex_fn = chebyshev_interp_threaded(yex, aex)

yr_fn(point)
ye_fn(point)
yex_fn(point)

yr_grad = chebyshev_gradient(wr_tensor, point, order, dom)
ye_grad = chebyshev_gradient(we_tensor, point, order, dom)
yex_grad = chebyshev_gradient(wex_tensor, point, order, dom)

yr_grad = chebyshev_gradient(yr, ar)
ye_grad = chebyshev_gradient(ye, ae)
yex_grad = chebyshev_gradient(yex, aex)

yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_grad = chebyshev_gradient_threaded(yr, ar)
ye_grad = chebyshev_gradient_threaded(ye, ae)
yex_grad = chebyshev_gradient_threaded(yex, aex)

yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_grad = chebyshev_gradient(yr, br)
ye_grad = chebyshev_gradient(ye, be)
yex_grad = chebyshev_gradient(yex, bex)

yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_grad = chebyshev_gradient_threaded(yr, br)
ye_grad = chebyshev_gradient_threaded(ye, be)
yex_grad = chebyshev_gradient_threaded(yex, bex)

yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_hess = chebyshev_hessian(yr, ar)
ye_hess = chebyshev_hessian(ye, ae)
yex_hess = chebyshev_hessian(yex, aex)

yr_hess(point)
ye_hess(point)
yex_hess(point)

yr_hess = chebyshev_hessian_threaded(yr, ar)
ye_hess = chebyshev_hessian_threaded(ye, ae)
yex_hess = chebyshev_hessian_threaded(yex, aex)

yr_hess(point)
ye_hess(point)
yex_hess(point)

yr_hess = chebyshev_hessian(yr, br)
ye_hess = chebyshev_hessian(ye, be)
yex_hess = chebyshev_hessian(yex, bex)

yr_hess(point)
ye_hess(point)
yex_hess(point)

yr_hess = chebyshev_hessian_threaded(yr, br)
ye_hess = chebyshev_hessian_threaded(ye, be)
yex_hess = chebyshev_hessian_threaded(yex, bex)

yr_hess(point)
ye_hess(point)
yex_hess(point)

####################################
# Test the 2-variable case, BigFloat
####################################

n1 = 51
n2 = 51

dom_1 = [2.0, -3.0]
dom_2 = [1.5, 0.5]
dom = [dom_1 dom_2]

nodesr_1 = chebyshev_nodes(BigFloat, n1, dom_1)
nodesr_2 = chebyshev_nodes(BigFloat, n2, dom_2)
nodese_1 = chebyshev_extrema(BigFloat, n1, dom_1)
nodese_2 = chebyshev_extrema(BigFloat, n2, dom_2)
nodesex_1 = chebyshev_extended(BigFloat, n1, dom_1)
nodesex_2 = chebyshev_extended(BigFloat, n2, dom_2)

gr_1 = nodes(BigFloat, n1, :chebyshev_nodes, dom_1)
gr_2 = nodes(BigFloat, n2, :chebyshev_nodes, dom_2)
ge_1 = nodes(BigFloat, n1, :chebyshev_extrema, dom_1)
ge_2 = nodes(BigFloat, n2, :chebyshev_extrema, dom_2)
gex_1 = nodes(BigFloat, n1, :chebyshev_extended, dom_1)
gex_2 = nodes(BigFloat, n2, :chebyshev_extended, dom_2)

yr  = [(nodesr_1[i] + 4.0)^0.5 + nodesr_1[i] * sqrt(nodesr_2[j]) for i in 1:n1, j in 1:n2]
ye  = [(nodese_1[i] + 4.0)^0.5 + nodese_1[i] * sqrt(nodese_2[j]) for i in 1:n1, j in 1:n2]
yex = [(nodesex_1[i] + 4.0)^0.5 + nodesex_1[i] * sqrt(nodesex_2[j]) for i in 1:n1, j in 1:n2]

order_tensor = [10, 10]
order_complete = 10

wr_tensor = chebyshev_weights(yr, [nodesr_1, nodesr_2], order_tensor, dom)
wr_complete = chebyshev_weights(yr, (nodesr_1, nodesr_2), order_complete, dom)

we_tensor = chebyshev_weights_extrema(ye, (nodese_1, nodese_2), order_tensor, dom)
we_complete = chebyshev_weights_extrema(ye, (nodese_1, nodese_2), order_complete, dom)

wex_tensor = chebyshev_weights_extended(yex, (nodesex_1, nodesex_2), order_tensor, dom)
wex_complete = chebyshev_weights_extended(yex, (nodesex_1, nodesex_2), order_complete, dom)

order = (10, 10)
hr = Grid((gr_1, gr_2))
he = Grid((ge_1, ge_2))
hex = Grid((gex_1, gex_2))

ar = CApproxPlan(hr, order, dom)
ae = CApproxPlan(he, order, dom)
aex = CApproxPlan(hex, order, dom)

pr_1 = chebyshev_polynomial(order[1], gr_1)
pr_2 = chebyshev_polynomial(order[2], gr_2)
pe_1 = chebyshev_polynomial(order[1], ge_1)
pe_2 = chebyshev_polynomial(order[2], ge_2)
pex_1 = chebyshev_polynomial(order[1], gex_1)
pex_2 = chebyshev_polynomial(order[2], gex_2)

br = CApproxPlanPoly((pr_1, pr_2), order, dom)
be = CApproxPlanPoly((pe_1, pe_2), order, dom)
bex = CApproxPlanPoly((pex_1, pex_2), order, dom)

wr = chebyshev_weights(yr, ar)
we = chebyshev_weights(ye, ae)
wex = chebyshev_weights(yex, aex)

wr_t = chebyshev_weights_threaded(yr, ar)
we_t = chebyshev_weights_threaded(ye, ae)
wex_t = chebyshev_weights_threaded(yex, aex)

wr = chebyshev_weights(yr, br)
we = chebyshev_weights(ye, be)
wex = chebyshev_weights(yex, bex)

wr_t = chebyshev_weights_threaded(yr, br)
we_t = chebyshev_weights_threaded(ye, be)
wex_t = chebyshev_weights_threaded(yex, bex)

point = BigFloat.([1.748, 0.753])

yr_hat = chebyshev_evaluate(wr_tensor, point, order, dom)
ye_hat = chebyshev_evaluate(we_tensor, point, order, dom)
yex_hat = chebyshev_evaluate(wex_tensor, point, order, dom)

yr_fn = chebyshev_interp(yr, ar)
ye_fn = chebyshev_interp(ye, ae)
yex_fn = chebyshev_interp(yex, aex)
yr_fn(point)
ye_fn(point)
yex_fn(point)

yr_fn = chebyshev_interp_threaded(yr, ar)
ye_fn = chebyshev_interp_threaded(ye, ae)
yex_fn = chebyshev_interp_threaded(yex, aex)
yr_fn(point)
ye_fn(point)
yex_fn(point)

yr_grad = chebyshev_gradient(wr_tensor, point, order, dom)
ye_grad = chebyshev_gradient(we_tensor, point, order, dom)
yex_grad = chebyshev_gradient(wex_tensor, point, order, dom)

yr_grad = chebyshev_gradient(yr, ar)
ye_grad = chebyshev_gradient(ye, ae)
yex_grad = chebyshev_gradient(yex, aex)
yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_grad = chebyshev_gradient_threaded(yr, ar)
ye_grad = chebyshev_gradient_threaded(ye, ae)
yex_grad = chebyshev_gradient_threaded(yex, aex)
yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_grad = chebyshev_gradient(yr, br)
ye_grad = chebyshev_gradient(ye, be)
yex_grad = chebyshev_gradient(yex, bex)
yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_grad = chebyshev_gradient_threaded(yr, br)
ye_grad = chebyshev_gradient_threaded(ye, be)
yex_grad = chebyshev_gradient_threaded(yex, bex)
yr_grad(point)
ye_grad(point)
yex_grad(point)

yr_hess = chebyshev_hessian(yr, ar)
ye_hess = chebyshev_hessian(ye, ae)
yex_hess = chebyshev_hessian(yex, aex)
yr_hess(point)
ye_hess(point)
yex_hess(point)

yr_hess = chebyshev_hessian_threaded(yr, ar)
ye_hess = chebyshev_hessian_threaded(ye, ae)
yex_hess = chebyshev_hessian_threaded(yex, aex)
yr_hess(point)
ye_hess(point)
yex_hess(point)

yr_hess = chebyshev_hessian(yr, br)
ye_hess = chebyshev_hessian(ye, be)
yex_hess = chebyshev_hessian(yex, bex)
yr_hess(point)
ye_hess(point)
yex_hess(point)

yr_hess = chebyshev_hessian_threaded(yr, br)
ye_hess = chebyshev_hessian_threaded(ye, be)
yex_hess = chebyshev_hessian_threaded(yex, bex)
yr_hess(point)
ye_hess(point)
yex_hess(point)