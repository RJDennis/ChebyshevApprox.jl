# Test file for 1-d Chebyshev polynomials

using ChebyshevApprox

n = 10
dom = [2.0,-3.0]

nodes = chebyshev_nodes(n,dom)

y = zeros(n)

for i = 1:n

  y[i] = (nodes[i]+4.0)^0.5

end

order = 5

w_tensor       = chebyshev_weights(y,nodes,[order],dom)
w_tensor_gen   = chebyshev_weights(y,(nodes,),[order],reshape(dom,2,1))
w_complete     = chebyshev_weights(y,nodes,order,dom)
w_complete_gen = chebyshev_weights(y,(nodes,),order,reshape(dom,2,1))

println(maximum(abs,w_tensor_gen - w_complete_gen)) # We run this in the 1-d case because the tensor polynomial and the complete polynomial should coincide

point = [1.748]

y_chebyshev_tensor   = chebyshev_evaluate(w_tensor_gen,point,[order],reshape(dom,2,1))
y_chebyshev_complete = chebyshev_evaluate(w_complete_gen,point,order,reshape(dom,2,1))
y_clenshaw_tensor    = clenshaw_evaluate(w_tensor_gen,point,dom)
y_clenshaw_complete  = clenshaw_evaluate(w_complete_gen,point,dom)

y_actual = (point+4.0).^0.5

println([y_actual y_chebyshev_tensor y_chebyshev_complete y_clenshaw_tensor y_clenshaw_complete])

y_deriv_tensor   = chebyshev_derivative(w_tensor_gen,point,[order],dom)
y_deriv_complete = chebyshev_derivative(w_tensor_gen,point,order,dom)
y_deriv_tensor_1   = chebyshev_derivative(w_tensor_gen,point,[order],[1],dom)
y_deriv_complete_1 = chebyshev_derivative(w_tensor_gen,point,order,[1],dom)

println(y_deriv_tensor, y_deriv_complete, y_deriv_tensor_1, y_deriv_complete_1)
