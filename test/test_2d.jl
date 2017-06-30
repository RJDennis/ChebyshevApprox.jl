# Test file for 2-d Chebyshev polynomials

using ChebyshevApprox

n1 = 10
n2 = 15

range_1 = [2.0,-3.0]
range_2 = [1.5, 0.5]
range = [range_1 range_2]

nodes_1 = chebyshev_nodes(n1,range_1)
nodes_2 = chebyshev_nodes(n2,range_2)

y = zeros(n1,n2)

for i = 1:n1
  for j = 1:n2

    y[i,j] = (nodes_1[i]+4.0)^0.5+nodes_1[i]*sqrt(nodes_2[j])

  end
end

order_tensor = [5, 6]
order_complete = 6

w_tensor   = chebyshev_weights(y,nodes_1,nodes_2,order_tensor,range)
w_complete = chebyshev_weights(y,nodes_1,nodes_2,order_complete,range)
#w_tensor_gen   = chebyshev_weights(y,(nodes_1,nodes_2),order_tensor,range)
#w_complete_gen = chebyshev_weights(y,(nodes_1,nodes_2),order_complete,range)

point = [1.748, 0.753]

y_chebyshev_tensor   = chebyshev_evaluate(w_tensor,point,order_tensor,range)
y_chebyshev_complete = chebyshev_evaluate(w_complete,point,order_complete,range)
y_clenshaw_tensor    = clenshaw_evaluate(w_tensor,point,order_tensor,range)
y_clenshaw_complete  = clenshaw_evaluate(w_complete,point,order_complete,range)

y_actual = (point[1]+4.0)^0.5+point[1]*sqrt(point[2])

println([y_actual y_chebyshev_tensor y_chebyshev_complete y_clenshaw_tensor y_clenshaw_complete])
