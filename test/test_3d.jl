# Test file for 3-d Chebyshev polynomials

using ChebyshevApprox

n1 = 10
n2 = 15
n3 = 7

range_1 = [2.0, -3.0]
range_2 = [1.5,  0.5]
range_3 = [0.5, -0.5]
range = [range_1 range_2 range_3]

nodes_1 = chebyshev_nodes(n1,range_1)
nodes_2 = chebyshev_nodes(n2,range_2)
nodes_3 = chebyshev_nodes(n3,range_3)

y = zeros(n1,n2,n3)

for i = 1:n1
  for j = 1:n2
    for k = 1:n3

      y[i,j,k] = (nodes_1[i]+4.0)^0.5+nodes_1[i]*sqrt(nodes_2[j])+exp(nodes_3[k])

    end
  end
end

order_tensor = [7, 6, 6]
order_complete = 6

w_tensor   = chebyshev_weights(y,nodes_1,nodes_2,nodes_3,order_tensor,range)
w_complete = chebyshev_weights(y,nodes_1,nodes_2,nodes_3,order_complete,range)
#w_tensor_gen   = chebyshev_weights(y,(nodes_1,nodes_2,nodes_3),order_tensor,range)
#w_complete_gen = chebyshev_weights(y,(nodes_1,nodes_2,nodes_3),order_complete,range)

point = [1.748, 0.753, 0.119]

y_chebyshev_tensor   = chebyshev_evaluate(w_tensor,point,order_tensor,range)
y_chebyshev_complete = chebyshev_evaluate(w_complete,point,order_complete,range)
y_clenshaw_tensor    = clenshaw_evaluate(w_tensor,point,order_tensor,range)
y_clenshaw_complete  = clenshaw_evaluate(w_complete,point,order_complete,range)

y_actual = (point[1]+4.0)^0.5+point[1]*sqrt(point[2])+exp(point[3])

println([y_actual y_chebyshev_tensor y_chebyshev_complete y_clenshaw_tensor y_clenshaw_complete])
