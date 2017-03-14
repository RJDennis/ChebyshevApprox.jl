# Test file for 4-d Chebyshev polynomials

using ChebyshevApprox

n1 = 10
n2 = 15
n3 = 7
n4 = 16

range_1 = [2.0, -3.0]
range_2 = [1.5,  0.5]
range_3 = [0.5, -0.5]
range_4 = [-0.66, -1.79]
range = [range_1 range_2 range_3 range_4]

nodes_1 = chebyshev_nodes(n1,range_1)
nodes_2 = chebyshev_nodes(n2,range_2)
nodes_3 = chebyshev_nodes(n3,range_3)
nodes_4 = chebyshev_nodes(n4,range_4)

y = zeros(n1,n2,n3,n4)

for i = 1:n1
  for j = 1:n2
    for k = 1:n3
      for l = 1:n4

        y[i,j,k,l] = (nodes_1[i]+4.0)^0.5+nodes_1[i]*sqrt(nodes_2[j])+exp(nodes_3[k])*nodes_4[l]

      end
    end
  end
end

order_tensor = [7, 6, 6, 6]
order_complete = 6

w_tensor   = chebyshev_weights(y,nodes_1,nodes_2,nodes_3,nodes_4,order_tensor,range)
w_complete = chebyshev_weights(y,nodes_1,nodes_2,nodes_3,nodes_4,order_complete,range)
w_tensor_gen   = chebyshev_weights(y,(nodes_1,nodes_2,nodes_3,nodes_4),order_tensor,range)
w_complete_gen = chebyshev_weights(y,(nodes_1,nodes_2,nodes_3,nodes_4),order_complete,range)

point = [1.748, 0.753, 0.119, -0.947]

y_chebyshev_tensor   = chebyshev_evaluate(w_tensor,point,order_tensor,range)
y_chebyshev_complete = chebyshev_evaluate(w_complete,point,order_complete,range)
y_clenshaw_tensor    = clenshaw_evaluate(w_tensor,point,order_tensor,range)
y_clenshaw_complete  = clenshaw_evaluate(w_complete,point,order_complete,range)

y_actual = (point[1]+4.0)^0.5+point[1]*sqrt(point[2])+exp(point[3])*point[4]

println([y_actual y_chebyshev_tensor y_chebyshev_complete y_clenshaw_tensor y_clenshaw_complete])
