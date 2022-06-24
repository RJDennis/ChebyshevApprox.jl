# Test file for 5-d Chebyshev polynomials

using ChebyshevApprox

n1 = 10
n2 = 15
n3 = 7
n4 = 16
n5 = 31

dom_1 = [2.0, -3.0]
dom_2 = [1.5,  0.5]
dom_3 = [0.5, -0.5]
dom_4 = [-0.66, -1.79]
dom_5 = [1.0, -1.0]
dom = [dom_1 dom_2 dom_3 dom_4 dom_5]

nodesr_1  = chebyshev_nodes(n1,dom_1)
nodesr_2  = chebyshev_nodes(n2,dom_2)
nodesr_3  = chebyshev_nodes(n3,dom_3)
nodesr_4  = chebyshev_nodes(n4,dom_4)
nodesr_5  = chebyshev_nodes(n5,dom_5)
nodese_1  = chebyshev_extrema(n1,dom_1)
nodese_2  = chebyshev_extrema(n2,dom_2)
nodese_3  = chebyshev_extrema(n3,dom_3)
nodese_4  = chebyshev_extrema(n4,dom_4)
nodese_5  = chebyshev_extrema(n5,dom_5)
nodesex_1 = chebyshev_extended(n1,dom_1)
nodesex_2 = chebyshev_extended(n2,dom_2)
nodesex_3 = chebyshev_extended(n3,dom_3)
nodesex_4 = chebyshev_extended(n4,dom_4)
nodesex_5 = chebyshev_extended(n5,dom_5)
nodesv_1  = vertesi_nodes(n1,dom_1)
nodesv_2  = vertesi_nodes(n2,dom_2)
nodesv_3  = vertesi_nodes(n3,dom_3)
nodesv_4  = vertesi_nodes(n4,dom_4)
nodesv_5  = vertesi_nodes(n5,dom_5)

yr  = zeros(n1,n2,n3,n4,n5)
ye  = zeros(n1,n2,n3,n4,n5)
yex = zeros(n1,n2,n3,n4,n5)
yv  = zeros(n1,n2,n3,n4,n5)

for i = 1:n1
  for j = 1:n2
    for k = 1:n3
      for l = 1:n4
        for m = 1:n5

          yr[i,j,k,l,m]  = (nodesr_1[i]+4.0)^0.5+nodesr_1[i]*sqrt(nodesr_2[j])+exp(nodesr_3[k])*nodesr_4[l]-(1.0+nodesr_5[m])^2
          ye[i,j,k,l,m]  = (nodese_1[i]+4.0)^0.5+nodese_1[i]*sqrt(nodese_2[j])+exp(nodese_3[k])*nodese_4[l]-(1.0+nodese_5[m])^2
          yex[i,j,k,l,m] = (nodesex_1[i]+4.0)^0.5+nodesex_1[i]*sqrt(nodesex_2[j])+exp(nodesex_3[k])*nodesex_4[l]-(1.0+nodesex_5[m])^2
          yv[i,j,k,l,m]  = (nodesv_1[i]+4.0)^0.5+nodesv_1[i]*sqrt(nodesv_2[j])+exp(nodesv_3[k])*nodesv_4[l]-(1.0+nodesv_5[m])^2

        end
      end
    end
  end
end

order_tensor = [6, 6, 6, 6, 6]
order_complete = 6

wr_tensor    = chebyshev_weights_threaded(yr,[nodesr_1,nodesr_2,nodesr_3,nodesr_4,nodesr_5],order_tensor,dom)
wr_complete  = chebyshev_weights_threaded(yr,(nodesr_1,nodesr_2,nodesr_3,nodesr_4,nodesr_5),order_complete,dom)
we_tensor    = chebyshev_weights_extrema_threaded(ye,(nodese_1,nodese_2,nodese_3,nodese_4,nodese_5),order_tensor,dom)
we_complete  = chebyshev_weights_extrema_threaded(ye,(nodese_1,nodese_2,nodese_3,nodese_4,nodese_5),order_complete,dom)
wex_tensor   = chebyshev_weights_extended_threaded(yex,(nodesex_1,nodesex_2,nodesex_3,nodesex_4,nodesex_5),order_tensor,dom)
wex_complete = chebyshev_weights_extended_threaded(yex,(nodesex_1,nodesex_2,nodesex_3,nodesex_4,nodesex_5),order_complete,dom)
wv_tensor    = chebyshev_weights_vertesi_threaded(yv,(nodesv_1,nodesv_2,nodesv_3,nodesv_4,nodesv_5),order_tensor,dom)
wv_complete  = chebyshev_weights_vertesi_threaded(yv,(nodesv_1,nodesv_2,nodesv_3,nodesv_4,nodesv_5),order_complete,dom)

point = [1.748, 0.753, 0.119, -0.947, -0.23]

yr_tensor    = chebyshev_evaluate(wr_tensor,point,order_tensor,dom)
yr_complete  = chebyshev_evaluate(wr_complete,point,order_complete,dom)
ye_tensor    = chebyshev_evaluate(we_tensor,point,order_tensor,dom)
ye_complete  = chebyshev_evaluate(we_complete,point,order_complete,dom)
yex_tensor   = chebyshev_evaluate(wex_tensor,point,order_tensor,dom)
yex_complete = chebyshev_evaluate(wex_complete,point,order_complete,dom)
yv_tensor    = chebyshev_evaluate(wv_tensor,point,order_tensor,dom)
yv_complete  = chebyshev_evaluate(wv_complete,point,order_complete,dom)

y_actual = (point[1]+4.0)^0.5+point[1]*sqrt(point[2])+exp(point[3])*point[4]-(1.0+point[5])^2

println([y_actual yr_tensor yr_complete ye_tensor ye_complete yex_tensor yex_complete yv_tensor yv_complete])

yr_grad_tensor       = chebyshev_gradient(wr_tensor,point,order_tensor,dom)
yr_grad_complete     = chebyshev_gradient(wr_complete,point,order_complete,dom)
yr_deriv_tensor_1    = chebyshev_derivative(wr_tensor,point,1,order_tensor,dom)
yr_deriv_complete_1  = chebyshev_derivative(wr_tensor,point,1,order_complete,dom)
ye_grad_tensor       = chebyshev_gradient(we_tensor,point,order_tensor,dom)
ye_grad_complete     = chebyshev_gradient(we_complete,point,order_complete,dom)
ye_deriv_tensor_1    = chebyshev_derivative(we_tensor,point,1,order_tensor,dom)
ye_deriv_complete_1  = chebyshev_derivative(we_complete,point,1,order_complete,dom)
yex_grad_tensor      = chebyshev_gradient(wex_tensor,point,order_tensor,dom)
yex_grad_complete    = chebyshev_gradient(wex_complete,point,order_complete,dom)
yex_deriv_tensor_1   = chebyshev_derivative(wex_tensor,point,1,order_tensor,dom)
yex_deriv_complete_1 = chebyshev_derivative(wex_complete,point,1,order_complete,dom)
yv_grad_tensor       = chebyshev_gradient(wv_tensor,point,order_tensor,dom)
yv_grad_complete     = chebyshev_gradient(wv_complete,point,order_complete,dom)
yv_deriv_tensor_1    = chebyshev_derivative(wv_tensor,point,1,order_tensor,dom)
yv_deriv_complete_1  = chebyshev_derivative(wv_complete,point,1,order_complete,dom)

println(yr_grad_tensor, yr_grad_complete, yr_deriv_tensor_1, yr_deriv_complete_1)
println(ye_grad_tensor, ye_grad_complete, ye_deriv_tensor_1, ye_deriv_complete_1)
println(yex_grad_tensor, yex_grad_complete, yex_deriv_tensor_1, yex_deriv_complete_1)
println(yv_grad_tensor, yv_grad_complete, yv_deriv_tensor_1, yv_deriv_complete_1)

# Test the structures

chebpoly = ChebPoly(wr_complete,order_complete,dom)
chebyshev_evaluate(chebpoly,point)

chebr  = ChebInterpRoots(yr,(nodesr_1,nodesr_2,nodesr_3,nodesr_4,nodesr_5),order_complete,dom)
chebe  = ChebInterpExtrema(ye,(nodese_1,nodese_2,nodese_3,nodese_4,nodese_5),order_complete,dom)
chebex = ChebInterpExtended(yex,(nodesex_1,nodesex_2,nodesex_3,nodesex_4,nodesex_5),order_complete,dom)
chebv  = ChebInterpVertesi(yv,(nodesv_1,nodesv_2,nodesv_3,nodesv_4,nodesv_5),order_complete,dom)
cheb_interp(chebr,point)
cheb_interp(chebe,point)
cheb_interp(chebex,point)
cheb_interp(chebv,point)

chebr_fn  = cheb_interp(chebr)
chebe_fn  = cheb_interp(chebe)
chebex_fn = cheb_interp(chebex)
chebv_fn  = cheb_interp(chebv)

chebr_fn(point)
chebe_fn(point)
chebex_fn(point)
chebv_fn(point)

