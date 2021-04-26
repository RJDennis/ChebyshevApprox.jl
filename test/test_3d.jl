# Test file for 3-d Chebyshev polynomials

using ChebyshevApprox

n1 = 10
n2 = 15
n3 = 7

dom_1 = [2.0, -3.0]
dom_2 = [1.5,  0.5]
dom_3 = [0.5, -0.5]
dom = [dom_1 dom_2 dom_3]

nodesr_1 = chebyshev_nodes(n1,dom_1)
nodesr_2 = chebyshev_nodes(n2,dom_2)
nodesr_3 = chebyshev_nodes(n3,dom_3)
nodese_1 = chebyshev_extrema(n1,dom_1)
nodese_2 = chebyshev_extrema(n2,dom_2)
nodese_3 = chebyshev_extrema(n3,dom_3)

yr = zeros(n1,n2,n3)
ye = zeros(n1,n2,n3)

for i = 1:n1
  for j = 1:n2
    for k = 1:n3

      yr[i,j,k] = (nodesr_1[i]+4.0)^0.5+nodesr_1[i]*sqrt(nodesr_2[j])+exp(nodesr_3[k])
      ye[i,j,k] = (nodese_1[i]+4.0)^0.5+nodese_1[i]*sqrt(nodese_2[j])+exp(nodese_3[k])

    end
  end
end

order_tensor = [7, 6, 6]
order_complete = 6

wr_tensor   = chebyshev_weights(yr,[nodesr_1,nodesr_2,nodesr_3],order_tensor,dom)
wr_complete = chebyshev_weights(yr,(nodesr_1,nodesr_2,nodesr_3),order_complete,dom)
we_tensor   = chebyshev_weights_extrema(ye,(nodese_1,nodese_2,nodese_3),order_tensor,dom)
we_complete = chebyshev_weights_extrema(ye,(nodese_1,nodese_2,nodese_3),order_complete,dom)

point = [1.748, 0.753, 0.119]

yr_tensor    = chebyshev_evaluate(wr_tensor,point,order_tensor,dom)
yr_complete  = chebyshev_evaluate(wr_complete,point,order_complete,dom)
ye_tensor    = chebyshev_evaluate(we_tensor,point,order_tensor,dom)
ye_complete  = chebyshev_evaluate(we_complete,point,order_tensor,dom)

y_actual = (point[1]+4.0)^0.5+point[1]*sqrt(point[2])+exp(point[3])

println([y_actual yr_tensor yr_complete ye_tensor ye_complete])

yr_deriv_tensor     = chebyshev_gradient(wr_tensor,point,order_tensor,dom)
yr_deriv_complete   = chebyshev_gradient(wr_complete,point,order_complete,dom)
yr_deriv_tensor_1   = chebyshev_derivative(wr_tensor,point,1,order_tensor,dom)
yr_deriv_complete_1 = chebyshev_derivative(wr_tensor,point,1,order_complete,dom)
ye_deriv_tensor     = chebyshev_gradient(we_tensor,point,order_tensor,dom)
ye_deriv_complete   = chebyshev_gradient(we_complete,point,order_complete,dom)
ye_deriv_tensor_1   = chebyshev_derivative(we_tensor,point,1,order_tensor,dom)
ye_deriv_complete_1 = chebyshev_derivative(we_complete,point,1,order_complete,dom)

println(yr_deriv_tensor, yr_deriv_complete, yr_deriv_tensor_1, yr_deriv_complete_1)
println(ye_deriv_tensor, ye_deriv_complete, ye_deriv_tensor_1, ye_deriv_complete_1)
