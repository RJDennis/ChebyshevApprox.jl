# Test file for 1-d Chebyshev polynomials

using ChebyshevApprox

n = 10
dom = [2.0,-3.0]

nodesr = chebyshev_nodes(n,dom)
nodese = chebyshev_extrema(n,dom)

yr = zeros(n)
ye = zeros(n)

for i = 1:n

  yr[i] = (nodesr[i]+4.0)^0.5
  ye[i] = (nodese[i]+4.0)^0.5

end

order = 5

wr_tensor       = chebyshev_weights(yr,nodesr,[order],dom)
wr_complete     = chebyshev_weights(yr,nodesr,order,dom)
we_tensor       = chebyshev_weights_extrema(ye,nodese,[order],dom)
we_complete     = chebyshev_weights_extrema(ye,nodese,order,dom)

println(maximum(abs,wr_tensor - wr_complete)) # We run this becuase in the 1-d case the weights hsould be the same.
println(maximum(abs,wr_tensor - we_tensor))   # We run this becuase the weights should be the same.

point = [1.748]

yr_tensor   = chebyshev_evaluate(wr_tensor,point,[order],reshape(dom,2,1))
yr_complete = chebyshev_evaluate(wr_complete,point,order,reshape(dom,2,1))
ye_tensor   = chebyshev_evaluate(we_tensor,point,[order],reshape(dom,2,1))
ye_complete = chebyshev_evaluate(we_complete,point,order,reshape(dom,2,1))

y_actual = (point.+4.0).^0.5

println([y_actual yr_tensor yr_complete ye_tensor ye_complete])

yr_deriv_tensor     = chebyshev_gradient(wr_tensor,point,[order],dom)
yr_deriv_complete   = chebyshev_gradient(wr_tensor,point,order,dom)
yr_deriv_tensor_1   = chebyshev_derivative(wr_tensor,point,1,[order],dom)
yr_deriv_complete_1 = chebyshev_derivative(wr_tensor,point,1,order,dom)
ye_deriv_tensor     = chebyshev_gradient(we_tensor,point,[order],dom)
ye_deriv_complete   = chebyshev_gradient(we_tensor,point,order,dom)
ye_deriv_tensor_1   = chebyshev_derivative(we_tensor,point,1,[order],dom)
ye_deriv_complete_1 = chebyshev_derivative(we_tensor,point,1,order,dom)

println(yr_deriv_tensor, yr_deriv_complete, yr_deriv_tensor_1, yr_deriv_complete_1)
println(ye_deriv_tensor, ye_deriv_complete, ye_deriv_tensor_1, ye_deriv_complete_1)
