function vertesi_nodes(n::S,domain = [1.0,-1.0]) where {S <: Integer}

  if n <= 0
    error("The number of nodes must be positive.")
  end
        
  nodes = zeros(n)
        
  if isodd(n)               
    nodes[Int((n-1)/2)+1] = (domain[1]+domain[2])/2.0
  end
  
  if n == 1

    return nodes
  
  elseif n == 2
  
    nodes[1] = domain[2]
    nodes[n] = domain[1]
  
    return nodes
  
  else
        
    nodes[1] = domain[2]
    nodes[n] = domain[1]
  
    for i = 2:div(n,2)
      x = -cos((pi/2)*(2(i-1)+1.0)/(n))/cos((pi/(2*(n)))*(1.0+1/(4*log(n))))*(domain[1]-domain[2])/2.0
      nodes[i]     = (domain[1]+domain[2])/2.0 + x
      nodes[n-i+1] = (domain[1]+domain[2])/2.0 - x
    end
        
    return nodes
      
  end
      
end