# Test file for 1-d Chebyshev polynomials

using ChebyshevApprox

n = 10
dom = [2.0,-3.0]

nodesr  = chebyshev_nodes(n,dom)
nodese  = chebyshev_extrema(n,dom)
nodesex = chebyshev_extended(n,dom)
nodesv  = chebyshev_vertesi(n,dom)

yr  = zeros(n)
ye  = zeros(n)
yex = zeros(n)
yv  = zeros(n)

for i = 1:n

  yr[i]  = (nodesr[i]+4.0)^0.5
  ye[i]  = (nodese[i]+4.0)^0.5
  yex[i] = (nodesex[i]+4.0)^0.5
  yv[i]  = (nodesv[i]+4.0)^0.5

end

order = 5

wr_tensor    = chebyshev_weights(yr,nodesr,[order],dom)
wr_complete  = chebyshev_weights(yr,nodesr,order,dom)
we_tensor    = chebyshev_weights_extrema(ye,nodese,[order],dom)
we_complete  = chebyshev_weights_extrema(ye,nodese,order,dom)
wex_tensor   = chebyshev_weights_extended(yex,nodesex,[order],dom)
wex_complete = chebyshev_weights_extended(yex,nodesex,order,dom)
wv_tensor    = chebyshev_weights_vertesi(yv,nodesv,[order],dom)
wv_complete  = chebyshev_weights_vertesi(yv,nodesv,order,dom)

println(maximum(abs,wr_tensor  - wr_complete))   # We run this becuase in the 1-d case the weights hsould be the same.
println(maximum(abs,we_tensor  - we_complete))   # We run this becuase the weights should be the same.
println(maximum(abs,wex_tensor - wex_complete))  # We run this becuase the weights should be the same.
println(maximum(abs,wv_tensor  - wv_complete))   # We run this becuase the weights should be the same.

point = [1.748]

yr_tensor    = chebyshev_evaluate(wr_tensor,point,[order],reshape(dom,2,1))
yr_complete  = chebyshev_evaluate(wr_complete,point,order,reshape(dom,2,1))
ye_tensor    = chebyshev_evaluate(we_tensor,point,[order],reshape(dom,2,1))
ye_complete  = chebyshev_evaluate(we_complete,point,order,reshape(dom,2,1))
yex_tensor   = chebyshev_evaluate(wex_tensor,point,[order],reshape(dom,2,1))
yex_complete = chebyshev_evaluate(wex_complete,point,order,reshape(dom,2,1))
yv_tensor    = chebyshev_evaluate(wv_tensor,point,[order],reshape(dom,2,1))
yv_complete  = chebyshev_evaluate(wv_complete,point,order,reshape(dom,2,1))

y_actual = (point.+4.0).^0.5

println([y_actual yr_tensor yr_complete ye_tensor ye_complete yex_tensor yex_complete yv_tensor yv_complete])

yr_grad_tensor      = chebyshev_gradient(wr_tensor,point,[order],dom)
yr_grad_complete    = chebyshev_gradient(wr_complete,point,order,dom)
yr_deriv_tensor_1   = chebyshev_derivative(wr_tensor,point,1,[order],dom)
yr_deriv_complete_1 = chebyshev_derivative(wr_complete,point,1,order,dom)

ye_grad_tensor      = chebyshev_gradient(we_tensor,point,[order],dom)
ye_grad_complete    = chebyshev_gradient(we_complete,point,order,dom)
ye_deriv_tensor_1   = chebyshev_derivative(we_tensor,point,1,[order],dom)
ye_deriv_complete_1 = chebyshev_derivative(we_complete,point,1,order,dom)

yex_grad_tensor      = chebyshev_gradient(wex_tensor,point,[order],dom)
yex_grad_complete    = chebyshev_gradient(wex_complete,point,order,dom)
yex_deriv_tensor_1   = chebyshev_derivative(wex_tensor,point,1,[order],dom)
yex_deriv_complete_1 = chebyshev_derivative(wex_complete,point,1,order,dom)

yv_grad_tensor      = chebyshev_gradient(wv_tensor,point,[order],dom)
yv_grad_complete    = chebyshev_gradient(wv_complete,point,order,dom)
yv_deriv_tensor_1   = chebyshev_derivative(wv_tensor,point,1,[order],dom)
yv_deriv_complete_1 = chebyshev_derivative(wv_complete,point,1,order,dom)

println(yr_grad_tensor, yr_grad_complete, yr_deriv_tensor_1, yr_deriv_complete_1)
println(ye_grad_tensor, ye_grad_complete, ye_deriv_tensor_1, ye_deriv_complete_1)
println(yex_grad_tensor, yex_grad_complete, yex_deriv_tensor_1, yex_deriv_complete_1)
println(yv_grad_tensor, yv_grad_complete, yv_deriv_tensor_1, yv_deriv_complete_1)
