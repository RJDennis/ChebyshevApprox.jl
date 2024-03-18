##################################

### Testing

##################################

using ChebyshevApprox

# Test the construction of nodes

nodes(111, :chebyshev_nodes)
nodes(111, :chebyshev_extrema)
nodes(111, :chebyshev_extended)

# Test the construction of polynomials

g = nodes(1001, :chebyshev_nodes)

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

# Test the 1-variable case

n = 10
order = 5
dom = [2.0, -3.0]

nodesr = chebyshev_nodes(n, dom)
nodese = chebyshev_extrema(n, dom)
nodesex = chebyshev_extended(n, dom)

gr = nodes(n, :chebyshev_nodes, dom)
ge = nodes(n, :chebyshev_extrema, dom)
gex = nodes(n, :chebyshev_extended, dom)

yr = zeros(n)
ye = zeros(n)
yex = zeros(n)

for i = 1:n

  yr[i] = (nodesr[i] + 4.0)^0.5
  ye[i] = (nodese[i] + 4.0)^0.5
  yex[i] = (nodesex[i] + 4.0)^0.5
  
end

wr_tensor = chebyshev_weights(yr, nodesr, order, dom)
wr_complete = chebyshev_weights(yr, nodesr, order, dom)
we_tensor = chebyshev_weights_extrema(ye, nodese, order, dom)
we_complete = chebyshev_weights_extrema(ye, nodese, [order], dom)
wex_tensor = chebyshev_weights_extended(yex, nodesex, order, dom)
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

# Test the 2-variable case

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

yr = zeros(n1, n2)
ye = zeros(n1, n2)
yex = zeros(n1, n2)

for i = 1:n1
  for j = 1:n2

    yr[i, j] = (nodesr_1[i] + 4.0)^0.5 + nodesr_1[i] * sqrt(nodesr_2[j])
    ye[i, j] = (nodese_1[i] + 4.0)^0.5 + nodese_1[i] * sqrt(nodese_2[j])
    yex[i, j] = (nodesex_1[i] + 4.0)^0.5 + nodesex_1[i] * sqrt(nodesex_2[j])

  end
end

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

br = CApproxPlanPoly((pr_1, pr_2), ord, dom)
be = CApproxPlanPoly((pe_1, pe_2), ord, dom)
bex = CApproxPlanPoly((pex_1, pex_2), ord, dom)

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
