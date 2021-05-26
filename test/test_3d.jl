# Test file for 3-d Chebyshev polynomials

using ChebyshevApprox

n1 = 10
n2 = 15
n3 = 7

dom_1 = [2.0, -3.0]
dom_2 = [1.5,  0.5]
dom_3 = [0.5, -0.5]
dom = [dom_1 dom_2 dom_3]

nodesr_1  = chebyshev_nodes(n1,dom_1)
nodesr_2  = chebyshev_nodes(n2,dom_2)
nodesr_3  = chebyshev_nodes(n3,dom_3)
nodesex_1 = chebyshev_extended(n1,dom_1)
nodesex_2 = chebyshev_extended(n2,dom_2)
nodesex_3 = chebyshev_extended(n3,dom_3)

yr  = zeros(n1,n2,n3)
yex = zeros(n1,n2,n3)

for i = 1:n1
  for j = 1:n2
    for k = 1:n3

      yr[i,j,k]  = (nodesr_1[i]+4.0)^0.5+nodesr_1[i]*sqrt(nodesr_2[j])+exp(nodesr_3[k])
      yex[i,j,k] = (nodesex_1[i]+4.0)^0.5+nodesex_1[i]*sqrt(nodesex_2[j])+exp(nodesex_3[k])

    end
  end
end

order_tensor = [7, 6, 6]
order_complete = 6

wr_tensor    = chebyshev_weights(yr,[nodesr_1,nodesr_2,nodesr_3],order_tensor,dom)
wr_complete  = chebyshev_weights(yr,(nodesr_1,nodesr_2,nodesr_3),order_complete,dom)
wex_tensor   = chebyshev_weights_extended(yex,(nodesex_1,nodesex_2,nodesex_3),order_tensor,dom)
wex_complete = chebyshev_weights_extended(yex,(nodesex_1,nodesex_2,nodesex_3),order_complete,dom)

point = [1.748, 0.753, 0.119]

yr_tensor     = chebyshev_evaluate(wr_tensor,point,order_tensor,dom)
yr_complete   = chebyshev_evaluate(wr_complete,point,order_complete,dom)
yex_tensor    = chebyshev_evaluate(wex_tensor,point,order_tensor,dom)
yex_complete  = chebyshev_evaluate(wex_complete,point,order_tensor,dom)

y_actual = (point[1]+4.0)^0.5+point[1]*sqrt(point[2])+exp(point[3])

println([y_actual yr_tensor yr_complete yex_tensor yex_complete])

yr_grad_tensor       = chebyshev_gradient(wr_tensor,point,order_tensor,dom)
yr_grad_complete     = chebyshev_gradient(wr_complete,point,order_complete,dom)
yr_deriv_tensor_1    = chebyshev_derivative(wr_tensor,point,1,order_tensor,dom)
yr_deriv_complete_1  = chebyshev_derivative(wr_tensor,point,1,order_complete,dom)
yex_grad_tensor      = chebyshev_gradient(wex_tensor,point,order_tensor,dom)
yex_grad_complete    = chebyshev_gradient(wex_complete,point,order_complete,dom)
yex_deriv_tensor_1   = chebyshev_derivative(wex_tensor,point,1,order_tensor,dom)
yex_deriv_complete_1 = chebyshev_derivative(wex_complete,point,1,order_complete,dom)

println(yr_grad_tensor, yr_grad_complete, yr_deriv_tensor_1, yr_deriv_complete_1)
println(yex_grad_tensor, yex_grad_complete, yex_deriv_tensor_1, yex_deriv_complete_1)
