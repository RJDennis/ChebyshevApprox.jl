##################################

### Testing

##################################

using ChebyshevApprox

# Test the construction of nodes

nodes(111, :chebyshev_nodes)
nodes(111, :chebyshev_extrema)
nodes(111, :chebyshev_extended)
nodes(111, :vertesi_nodes)
nodes(111, :legendre_nodes)

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
a = ApproxPlan(h, ord, dom)

p1 = chebyshev_polynomial(ord, g)
p2 = chebyshev_polynomial(ord, g)
b = ApproxPlanPoly((p1, p2), ord, dom)

# Test the 1-variable case

n = 10
order = 5
dom = [2.0, -3.0]

nodesr = chebyshev_nodes(n, dom)
nodese = chebyshev_extrema(n, dom)
nodesex = chebyshev_extended(n, dom)
nodesv = vertesi_nodes(n, dom)
nodesl = legendre_nodes(n, dom)

gr = nodes(n, :chebyshev_nodes, dom)
ge = nodes(n, :chebyshev_extrema, dom)
gex = nodes(n, :chebyshev_extended, dom)
gv = nodes(n, :vertesi_nodes, dom)
gl = nodes(n, :legendre_nodes, dom)

yr = zeros(n)
ye = zeros(n)
yex = zeros(n)
yv = zeros(n)
yl = zeros(n)

for i = 1:n

  yr[i] = (nodesr[i] + 4.0)^0.5
  ye[i] = (nodese[i] + 4.0)^0.5
  yex[i] = (nodesex[i] + 4.0)^0.5
  yv[i] = (nodesv[i] + 4.0)^0.5
  yl[i] = (nodesl[i] + 4.0)^0.5

end

wr_tensor = chebyshev_weights(yr, nodesr, order, dom)
wr_complete = chebyshev_weights(yr, nodesr, order, dom)
we_tensor = chebyshev_weights_extrema(ye, nodese, order, dom)
we_complete = chebyshev_weights_extrema(ye, nodese, [order], dom)
wex_tensor = chebyshev_weights_extended(yex, nodesex, order, dom)
wex_complete = chebyshev_weights_extended(yex, nodesex, order, dom)
wv_tensor = chebyshev_weights_vertesi(yv, nodesv, order, dom)
wv_complete = chebyshev_weights_vertesi(yv, nodesv, order, dom)
wl_tensor = chebyshev_weights_legendre(yl, nodesl, order, dom)
wl_complete = chebyshev_weights_legendre(yl, nodesl, order, dom)

hr = Grid((gr,))
he = Grid((ge,))
hex = Grid((gex,))
hv = Grid((gv,))
hl = Grid((gl,))

ar = ApproxPlan(hr, order, dom)
ae = ApproxPlan(he, order, dom)
aex = ApproxPlan(hex, order, dom)
av = ApproxPlan(hv, order, dom)
al = ApproxPlan(hl, order, dom)

pr = chebyshev_polynomial(order, gr)
pe = chebyshev_polynomial(order, ge)
pex = chebyshev_polynomial(order, gex)
pv = chebyshev_polynomial(order, gv)
pl = chebyshev_polynomial(order, gl)

br = ApproxPlanPoly((pr,), order, dom)
be = ApproxPlanPoly((pe,), order, dom)
bex = ApproxPlanPoly((pex,), order, dom)
bv = ApproxPlanPoly((pv,), order, dom)
bl = ApproxPlanPoly((pl,), order, dom)

wr = chebyshev_weights(yr, ar)
we = chebyshev_weights(ye, ae)
wex = chebyshev_weights(yex, aex)
wv = chebyshev_weights(yv, av)
wl = chebyshev_weights(yl, al)

wr_t = chebyshev_weights_threaded(yr, ar)
we_t = chebyshev_weights_threaded(ye, ae)
wex_t = chebyshev_weights_threaded(yex, aex)
wv_t = chebyshev_weights_threaded(yv, av)
wl_t = chebyshev_weights_threaded(yl, al)

wr = chebyshev_weights(yr, br)
we = chebyshev_weights(ye, be)
wex = chebyshev_weights(yex, bex)
wv = chebyshev_weights(yv, bv)
wl = chebyshev_weights(yl, bl)

wr_t = chebyshev_weights_threaded(yr, br)
we_t = chebyshev_weights_threaded(ye, be)
wex_t = chebyshev_weights_threaded(yex, bex)
wv_t = chebyshev_weights_threaded(yv, bv)
wl_t = chebyshev_weights_threaded(yl, bl)

point = [1.44]

yr_hat = chebyshev_evaluate(wr_tensor, point, order, dom)
ye_hat = chebyshev_evaluate(we_tensor, point, order, dom)
yex_hat = chebyshev_evaluate(wex_tensor, point, order, dom)
yv_hat = chebyshev_evaluate(wv_tensor, point, order, dom)
yl_hat = chebyshev_evaluate(wl_tensor, point, order, dom)

yr_fn = chebyshev_interp(yr, ar)
ye_fn = chebyshev_interp(ye, ae)
yex_fn = chebyshev_interp(yex, aex)
yv_fn = chebyshev_interp(yv, av)
yl_fn = chebyshev_interp(yl, al)
yr_fn(point)
ye_fn(point)
yex_fn(point)
yv_fn(point)
yl_fn(point)

yr_fn = chebyshev_interp_threaded(yr, ar)
ye_fn = chebyshev_interp_threaded(ye, ae)
yex_fn = chebyshev_interp_threaded(yex, aex)
yv_fn = chebyshev_interp_threaded(yv, av)
yl_fn = chebyshev_interp_threaded(yl, al)
yr_fn(point)
ye_fn(point)
yex_fn(point)
yv_fn(point)
yl_fn(point)

yr_grad = chebyshev_gradient(wr_tensor, point, order, dom)
ye_grad = chebyshev_gradient(we_tensor, point, order, dom)
yex_grad = chebyshev_gradient(wex_tensor, point, order, dom)
yv_grad = chebyshev_gradient(wv_tensor, point, order, dom)
yl_grad = chebyshev_gradient(wl_tensor, point, order, dom)

yr_grad = chebyshev_gradient(yr, ar)
ye_grad = chebyshev_gradient(ye, ae)
yex_grad = chebyshev_gradient(yex, aex)
yv_grad = chebyshev_gradient(yv, av)
yl_grad = chebyshev_gradient(yl, al)
yr_grad(point)
ye_grad(point)
yex_grad(point)
yv_grad(point)
yl_grad(point)

yr_grad = chebyshev_gradient_threaded(yr, ar)
ye_grad = chebyshev_gradient_threaded(ye, ae)
yex_grad = chebyshev_gradient_threaded(yex, aex)
yv_grad = chebyshev_gradient_threaded(yv, av)
yl_grad = chebyshev_gradient_threaded(yl, al)
yr_grad(point)
ye_grad(point)
yex_grad(point)
yv_grad(point)
yl_grad(point)

yr_grad = chebyshev_gradient(yr, br)
ye_grad = chebyshev_gradient(ye, be)
yex_grad = chebyshev_gradient(yex, bex)
yv_grad = chebyshev_gradient(yv, bv)
yl_grad = chebyshev_gradient(yl, bl)
yr_grad(point)
ye_grad(point)
yex_grad(point)
yv_grad(point)
yl_grad(point)

yr_grad = chebyshev_gradient_threaded(yr, br)
ye_grad = chebyshev_gradient_threaded(ye, be)
yex_grad = chebyshev_gradient_threaded(yex, bex)
yv_grad = chebyshev_gradient_threaded(yv, bv)
yl_grad = chebyshev_gradient_threaded(yl, bl)
yr_grad(point)
ye_grad(point)
yex_grad(point)
yv_grad(point)
yl_grad(point)

yr_hess = chebyshev_hessian(yr, ar)
ye_hess = chebyshev_hessian(ye, ae)
yex_hess = chebyshev_hessian(yex, aex)
yv_hess = chebyshev_hessian(yv, av)
yl_hess = chebyshev_hessian(yl, al)
yr_hess(point)
ye_hess(point)
yex_hess(point)
yv_hess(point)
yl_hess(point)

yr_hess = chebyshev_hessian_threaded(yr, ar)
ye_hess = chebyshev_hessian_threaded(ye, ae)
yex_hess = chebyshev_hessian_threaded(yex, aex)
yv_hess = chebyshev_hessian_threaded(yv, av)
yl_hess = chebyshev_hessian_threaded(yl, al)
yr_hess(point)
ye_hess(point)
yex_hess(point)
yv_hess(point)
yl_hess(point)

yr_hess = chebyshev_hessian(yr, br)
ye_hess = chebyshev_hessian(ye, be)
yex_hess = chebyshev_hessian(yex, bex)
yv_hess = chebyshev_hessian(yv, bv)
yl_hess = chebyshev_hessian(yl, bl)
yr_hess(point)
ye_hess(point)
yex_hess(point)
yv_hess(point)
yl_hess(point)

yr_hess = chebyshev_hessian_threaded(yr, br)
ye_hess = chebyshev_hessian_threaded(ye, be)
yex_hess = chebyshev_hessian_threaded(yex, bex)
yv_hess = chebyshev_hessian_threaded(yv, bv)
yl_hess = chebyshev_hessian_threaded(yl, bl)
yr_hess(point)
ye_hess(point)
yex_hess(point)
yv_hess(point)
yl_hess(point)

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
nodesv_1 = vertesi_nodes(n1, dom_1)
nodesv_2 = vertesi_nodes(n2, dom_2)
nodesl_1 = legendre_nodes(n1, dom_1)
nodesl_2 = legendre_nodes(n2, dom_2)

gr_1 = nodes(n1, :chebyshev_nodes, dom_1)
gr_2 = nodes(n2, :chebyshev_nodes, dom_2)
ge_1 = nodes(n1, :chebyshev_extrema, dom_1)
ge_2 = nodes(n2, :chebyshev_extrema, dom_2)
gex_1 = nodes(n1, :chebyshev_extended, dom_1)
gex_2 = nodes(n2, :chebyshev_extended, dom_2)
gv_1 = nodes(n1, :vertesi_nodes, dom_1)
gv_2 = nodes(n2, :vertesi_nodes, dom_2)
gl_1 = nodes(n1, :legendre_nodes, dom_1)
gl_2 = nodes(n2, :legendre_nodes, dom_2)

yr = zeros(n1, n2)
ye = zeros(n1, n2)
yex = zeros(n1, n2)
yv = zeros(n1, n2)
yl = zeros(n1, n2)

for i = 1:n1
  for j = 1:n2

    yr[i, j] = (nodesr_1[i] + 4.0)^0.5 + nodesr_1[i] * sqrt(nodesr_2[j])
    ye[i, j] = (nodese_1[i] + 4.0)^0.5 + nodese_1[i] * sqrt(nodese_2[j])
    yex[i, j] = (nodesex_1[i] + 4.0)^0.5 + nodesex_1[i] * sqrt(nodesex_2[j])
    yv[i, j] = (nodesv_1[i] + 4.0)^0.5 + nodesv_1[i] * sqrt(nodesv_2[j])
    yl[i, j] = (nodesl_1[i] + 4.0)^0.5 + nodesl_1[i] * sqrt(nodesl_2[j])

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
wv_tensor = chebyshev_weights_vertesi(yv, (nodesv_1, nodesv_2), order_tensor, dom)
wv_complete = chebyshev_weights_vertesi(yv, (nodesv_1, nodesv_2), order_complete, dom)
wl_tensor = chebyshev_weights_legendre(yl, (nodesl_1, nodesl_2), order_tensor, dom)
wl_complete = chebyshev_weights_legendre(yl, (nodesl_1, nodesl_2), order_complete, dom)

order = (10, 10)
hr = Grid((gr_1, gr_2))
he = Grid((ge_1, ge_2))
hex = Grid((gex_1, gex_2))
hv = Grid((gv_1, gv_2))
hl = Grid((gl_1, gl_2))

ar = ApproxPlan(hr, order, dom)
ae = ApproxPlan(he, order, dom)
aex = ApproxPlan(hex, order, dom)
av = ApproxPlan(hv, order, dom)
al = ApproxPlan(hl, order, dom)

pr_1 = chebyshev_polynomial(order[1], gr_1)
pr_2 = chebyshev_polynomial(order[2], gr_2)
pe_1 = chebyshev_polynomial(order[1], ge_1)
pe_2 = chebyshev_polynomial(order[2], ge_2)
pex_1 = chebyshev_polynomial(order[1], gex_1)
pex_2 = chebyshev_polynomial(order[2], gex_2)
pv_1 = chebyshev_polynomial(order[1], gv_1)
pv_2 = chebyshev_polynomial(order[2], gv_2)
pl_1 = chebyshev_polynomial(order[1], gl_1)
pl_2 = chebyshev_polynomial(order[2], gl_2)

br = ApproxPlanPoly((pr_1, pr_2), ord, dom)
be = ApproxPlanPoly((pe_1, pe_2), ord, dom)
bex = ApproxPlanPoly((pex_1, pex_2), ord, dom)
bv = ApproxPlanPoly((pv_1, pv_2), ord, dom)
bl = ApproxPlanPoly((pl_1, pl_2), ord, dom)

wr = chebyshev_weights(yr, ar)
we = chebyshev_weights(ye, ae)
wex = chebyshev_weights(yex, aex)
wv = chebyshev_weights(yv, av)
wl = chebyshev_weights(yl, al)

wr_t = chebyshev_weights_threaded(yr, ar)
we_t = chebyshev_weights_threaded(ye, ae)
wex_t = chebyshev_weights_threaded(yex, aex)
wv_t = chebyshev_weights_threaded(yv, av)
wl_t = chebyshev_weights_threaded(yl, al)

wr = chebyshev_weights(yr, br)
we = chebyshev_weights(ye, be)
wex = chebyshev_weights(yex, bex)
wv = chebyshev_weights(yv, bv)
wl = chebyshev_weights(yl, bl)

wr_t = chebyshev_weights_threaded(yr, br)
we_t = chebyshev_weights_threaded(ye, be)
wex_t = chebyshev_weights_threaded(yex, bex)
wv_t = chebyshev_weights_threaded(yv, bv)
wl_t = chebyshev_weights_threaded(yl, bl)

point = [1.748, 0.753]

yr_hat = chebyshev_evaluate(wr_tensor, point, order, dom)
ye_hat = chebyshev_evaluate(we_tensor, point, order, dom)
yex_hat = chebyshev_evaluate(wex_tensor, point, order, dom)
yv_hat = chebyshev_evaluate(wv_tensor, point, order, dom)
yl_hat = chebyshev_evaluate(wl_tensor, point, order, dom)

yr_fn = chebyshev_interp(yr, ar)
ye_fn = chebyshev_interp(ye, ae)
yex_fn = chebyshev_interp(yex, aex)
yv_fn = chebyshev_interp(yv, av)
yl_fn = chebyshev_interp(yl, al)
yr_fn(point)
ye_fn(point)
yex_fn(point)
yv_fn(point)
yl_fn(point)

yr_fn = chebyshev_interp_threaded(yr, ar)
ye_fn = chebyshev_interp_threaded(ye, ae)
yex_fn = chebyshev_interp_threaded(yex, aex)
yv_fn = chebyshev_interp_threaded(yv, av)
yl_fn = chebyshev_interp_threaded(yl, al)
yr_fn(point)
ye_fn(point)
yex_fn(point)
yv_fn(point)
yl_fn(point)

yr_grad = chebyshev_gradient(wr_tensor, point, order, dom)
ye_grad = chebyshev_gradient(we_tensor, point, order, dom)
yex_grad = chebyshev_gradient(wex_tensor, point, order, dom)
yv_grad = chebyshev_gradient(wv_tensor, point, order, dom)
yl_grad = chebyshev_gradient(wl_tensor, point, order, dom)

yr_grad = chebyshev_gradient(yr, ar)
ye_grad = chebyshev_gradient(ye, ae)
yex_grad = chebyshev_gradient(yex, aex)
yv_grad = chebyshev_gradient(yv, av)
yl_grad = chebyshev_gradient(yl, al)
yr_grad(point)
ye_grad(point)
yex_grad(point)
yv_grad(point)
yl_grad(point)

yr_grad = chebyshev_gradient_threaded(yr, ar)
ye_grad = chebyshev_gradient_threaded(ye, ae)
yex_grad = chebyshev_gradient_threaded(yex, aex)
yv_grad = chebyshev_gradient_threaded(yv, av)
yl_grad = chebyshev_gradient_threaded(yl, al)
yr_grad(point)
ye_grad(point)
yex_grad(point)
yv_grad(point)
yl_grad(point)

yr_grad = chebyshev_gradient(yr, br)
ye_grad = chebyshev_gradient(ye, be)
yex_grad = chebyshev_gradient(yex, bex)
yv_grad = chebyshev_gradient(yv, bv)
yl_grad = chebyshev_gradient(yl, bl)
yr_grad(point)
ye_grad(point)
yex_grad(point)
yv_grad(point)
yl_grad(point)

yr_grad = chebyshev_gradient_threaded(yr, br)
ye_grad = chebyshev_gradient_threaded(ye, be)
yex_grad = chebyshev_gradient_threaded(yex, bex)
yv_grad = chebyshev_gradient_threaded(yv, bv)
yl_grad = chebyshev_gradient_threaded(yl, bl)
yr_grad(point)
ye_grad(point)
yex_grad(point)
yv_grad(point)
yl_grad(point)

yr_hess = chebyshev_hessian(yr, ar)
ye_hess = chebyshev_hessian(ye, ae)
yex_hess = chebyshev_hessian(yex, aex)
yv_hess = chebyshev_hessian(yv, av)
yl_hess = chebyshev_hessian(yl, al)
yr_hess(point)
ye_hess(point)
yex_hess(point)
yv_hess(point)
yl_hess(point)

yr_hess = chebyshev_hessian_threaded(yr, ar)
ye_hess = chebyshev_hessian_threaded(ye, ae)
yex_hess = chebyshev_hessian_threaded(yex, aex)
yv_hess = chebyshev_hessian_threaded(yv, av)
yl_hess = chebyshev_hessian_threaded(yl, al)
yr_hess(point)
ye_hess(point)
yex_hess(point)
yv_hess(point)
yl_hess(point)

yr_hess = chebyshev_hessian(yr, br)
ye_hess = chebyshev_hessian(ye, be)
yex_hess = chebyshev_hessian(yex, bex)
yv_hess = chebyshev_hessian(yv, bv)
yl_hess = chebyshev_hessian(yl, bl)
yr_hess(point)
ye_hess(point)
yex_hess(point)
yv_hess(point)
yl_hess(point)

yr_hess = chebyshev_hessian_threaded(yr, br)
ye_hess = chebyshev_hessian_threaded(ye, be)
yex_hess = chebyshev_hessian_threaded(yex, bex)
yv_hess = chebyshev_hessian_threaded(yv, bv)
yl_hess = chebyshev_hessian_threaded(yl, bl)
yr_hess(point)
ye_hess(point)
yex_hess(point)
yv_hess(point)
yl_hess(point)