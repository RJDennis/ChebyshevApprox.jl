function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,1},x::Array{T,1},order::Array{S,1},range::Array{T,1})

  x1 = normalize_node(x[1],range)

  z = Array(T,(order[1]+1)+2)

  for i in (order[1]+1):-1:1

	  z[i] = weights[i]+2*x1*z[i+1]-z[i+2]

  end

  y = z[1]-x1*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,2},x::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])

  p = Array(T,order[2]+1)

	for j = 1:(order[2]+1)

    z = Array(T,(order[1]+1)+2)

    for i in (order[1]+1):-1:1

      z[i] = weights[i,j]+2*x1*z[i+1]-z[i+2]

    end

    p[j] = z[1]-x1*z[2]

  end

  z = Array(T,(order[2]+1)+2)

  for j in (order[2]+1):-1:1

    z[j] = p[j]+2*x2*z[j+1]-z[j+2]

  end

  y = z[1]-x2*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,3},x::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])

	pp = Array(T,(order[2]+1),(order[3]+1))

	for k  = 1:(order[3]+1)
		for j = 1:(order[2]+1)

      z = Array(T,(order[1]+1)+2)

      for i in (order[1]+1):-1:1

        z[i] = weights[i,j,k]+2*x1*z[i+1]-z[i+2]

      end

      pp[j,k] = z[1]-x1*z[2]

    end
	end

	p = Array(T,(order[3]+1))

	for k = 1:(order[3]+1)

    z = Array(T,(order[2]+1)+2)

    for j in (order[2]+1):-1:1

      z[j] = pp[j,k]+2*x2*z[j+1]-z[j+2]

    end

    p[k] = z[1]-x2*z[2]

  end

  z = Array(T,(order[3]+1)+2)

  for k in (order[3]+1):-1:1

    z[k] = p[k]+2*x3*z[k+1]-z[k+2]

  end

  y = z[1]-x3*z[2]

	return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,4},x::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])
  x4 = normalize_node(x[4],range[:,4])

	ppp = Array(T,(order[2]+1),(order[3]+1),(order[4]+1))

	for l = 1:(order[4]+1)
	  for k = 1:(order[3]+1)
		  for j = 1:(order[2]+1)

        z = Array(T,(order[1]+1)+2)

        for i in (order[1]+1):-1:1

          z[i] = weights[i,j,k,l]+2*x1*z[i+1]-z[i+2]

        end

        ppp[j,k,l] = z[1]-x1*z[2]

      end
	  end
	end

	pp = Array(T,(order[3]+1),(order[4]+1))

	for l = 1:(order[4]+1)
		for k = 1:(order[3]+1)

      z = Array(T,(order[2]+1)+2)

      for j in (order[2]+1):-1:1

        z[j] = ppp[j,k,l]+2*x2*z[j+1]-z[j+2]

      end

      pp[k,l] = z[1]-x2*z[2]

    end
	end

	p = Array(T,(order[4]+1))

	for l = 1:(order[4]+1)

    z = Array(T,(order[3]+1)+2)

    for k in (order[3]+1):-1:1

      z[k] = pp[k,l]+2*x3*z[k+1]-z[k+2]

    end

    p[l] = z[1]-x3*z[2]

  end

  z = Array(T,(order[4]+1)+2)

  for l in (order[4]+1):-1:1

    z[l] = p[l]+2*x4*z[l+1]-z[l+2]

  end

  y = z[1]-x4*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,5},x::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])
  x4 = normalize_node(x[4],range[:,4])
  x5 = normalize_node(x[5],range[:,5])

  pppp = Array(T,(order[2]+1),(order[3]+1),(order[4]+1),(order[5]+1))

  for m = 1:(order[5]+1)
	  for l = 1:(order[4]+1)
	    for k = 1:(order[3]+1)
		    for j = 1:(order[2]+1)

          z = Array(T,(order[1]+1)+2)

          for i in (order[1]+1):-1:1

            z[i] = weights[i,j,k,l,m]+2*x1*z[i+1]-z[i+2]

          end

          pppp[j,k,l,m] = z[1]-x1*z[2]

        end
	    end
	  end
  end

	ppp = Array(T,(order[3]+1),(order[4]+1),(order[5]+1))

	for m = 1:(order[5]+1)
	  for l = 1:(order[4]+1)
		  for k = 1:(order[3]+1)

        z = Array(T,(order[2]+1)+2)

        for j in (order[2]+1):-1:1

          z[j] = pppp[j,k,l,m]+2*x2*z[j+1]-z[j+2]

        end

        ppp[k,l,m] = z[1]-x2*z[2]

      end
	  end
	end

	pp = Array(T,(order[4]+1),(order[5]+1))

	for m = 1:(order[5]+1)
		for l = 1:(order[4]+1)

      z = Array(T,(order[3]+1)+2)

      for k in (order[3]+1):-1:1

        z[k] = ppp[k,l,m]+2*x3*z[k+1]-z[k+2]

      end

      pp[l,m] = z[1]-x3*z[2]

    end
	end

	p = Array(T,(order[5]+1))

	for m = 1:(order[5]+1)

    z = Array(T,(order[4]+1)+2)

    for l in (order[4]+1):-1:1

      z[l] = pp[l,m]+2*x4*z[l+1]-z[l+2]

    end

    p[m] = z[1]-x4*z[2]

  end

  z = Array(T,(order[5]+1)+2)

  for m in (order[5]+1):-1:1

    z[m] = p[m]+2*x5*z[m+1]-z[m+2]

  end

  y = z[1]-x5*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,6},x::Array{T,1},order::Array{S,1},range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])
  x4 = normalize_node(x[4],range[:,4])
  x5 = normalize_node(x[5],range[:,5])
  x6 = normalize_node(x[6],range[:,6])

  ppppp = Array(T,(order[2]+1),(order[3]+1),(order[4]+1),(order[5]+1),(order[6]+1))

  for n = 1:(order[6]+1)
    for m = 1:(order[5]+1)
	    for l = 1:(order[4]+1)
	      for k = 1:(order[3]+1)
		      for j = 1:(order[2]+1)

            z = Array(T,(order[1]+1)+2)

            for i in (order[1]+1):-1:1

              z[i] = weights[i,j,k,l,m,n]+2*x1*z[i+1]-z[i+2]

            end

            ppppp[j,k,l,m,n] = z[1]-x1*z[2]

          end
	      end
	    end
    end
  end

  pppp = Array(T,(order[3]+1),(order[4]+1),(order[5]+1),(order[6]+1))

  for m = 1:(order[6]+1)
	  for l = 1:(order[5]+1)
	    for k = 1:(order[4]+1)
		    for j = 1:(order[3]+1)

          z = Array(T,(order[2]+1)+2)

          for i in (order[2]+1):-1:1

            z[i] = ppppp[i,j,k,l,m]+2*x1*z[i+1]-z[i+2]

          end

          pppp[j,k,l,m] = z[1]-x1*z[2]

        end
	    end
	  end
  end

	ppp = Array(T,(order[4]+1),(order[5]+1),(order[6]+1))

	for m = 1:(order[6]+1)
	  for l = 1:(order[5]+1)
		  for k = 1:(order[4]+1)

        z = Array(T,(order[3]+1)+2)

        for j in (order[3]+1):-1:1

          z[j] = pppp[j,k,l,m]+2*x2*z[j+1]-z[j+2]

        end

        ppp[k,l,m] = z[1]-x2*z[2]

      end
	  end
	end

	pp = Array(T,(order[5]+1),(order[6]+1))

	for m = 1:(order[6]+1)
		for l = 1:(order[5]+1)

      z = Array(T,(order[4]+1)+2)

      for k in (order[4]+1):-1:1

        z[k] = ppp[k,l,m]+2*x3*z[k+1]-z[k+2]

      end

      pp[l,m] = z[1]-x3*z[2]

    end
	end

	p = Array(T,(order[6]+1))

	for m = 1:(order[6]+1)

    z = Array(T,(order[5]+1)+2)

    for l in (order[5]+1):-1:1

      z[l] = pp[l,m]+2*x4*z[l+1]-z[l+2]

    end

    p[m] = z[1]-x4*z[2]

  end

  z = Array(T,(order[6]+1)+2)

  for m in (order[6]+1):-1:1

    z[m] = p[m]+2*x5*z[m+1]-z[m+2]

  end

  y = z[1]-x5*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,1},x::Array{T,1},order::S,range::Array{T,1})

  x1 = normalize_node(x[1],range)

  z = Array(T,(order+1)+2)

  for i in (order+1):-1:1

	  z[i] = weights[i]+2*x1*z[i+1]-z[i+2]

  end

  y = z[1]-x1*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,2},x::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])

  p = Array(T,order+1)

	for j = 1:(order+1)

    z = Array(T,(order+1)+2)

    for i in (order+1):-1:1

      z[i] = weights[i,j]+2*x1*z[i+1]-z[i+2]

    end

    p[j] = z[1]-x1*z[2]

  end

  z = Array(T,(order+1)+2)

  for j in (order+1):-1:1

    z[j] = p[j]+2*x2*z[j+1]-z[j+2]

  end

  y = z[1]-x2*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,3},x::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])

	pp = Array(T,(order+1),(order+1))

	for k  = 1:(order+1)
		for j = 1:(order+1)

      z = Array(T,(order+1)+2)

      for i in (order+1):-1:1

        z[i] = weights[i,j,k]+2*x1*z[i+1]-z[i+2]

      end

      pp[j,k] = z[1]-x1*z[2]

    end
	end

	p = Array(T,(order+1))

	for k = 1:(order+1)

    z = Array(T,(order+1)+2)

    for j in (order+1):-1:1

      z[j] = pp[j,k]+2*x2*z[j+1]-z[j+2]

    end

    p[k] = z[1]-x2*z[2]

  end

  z = Array(T,(order+1)+2)

  for k in (order+1):-1:1

    z[k] = p[k]+2*x3*z[k+1]-z[k+2]

  end

  y = z[1]-x3*z[2]

	return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,4},x::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])
  x4 = normalize_node(x[4],range[:,4])

	ppp = Array(T,(order+1),(order+1),(order+1))

	for l = 1:(order+1)
	  for k = 1:(order+1)
		  for j = 1:(order+1)

        z = Array(T,(order+1)+2)

        for i in (order+1):-1:1

          z[i] = weights[i,j,k,l]+2*x1*z[i+1]-z[i+2]

        end

        ppp[j,k,l] = z[1]-x1*z[2]

      end
	  end
	end

	pp = Array(T,(order+1),(order+1))

	for l = 1:(order+1)
		for k = 1:(order+1)

      z = Array(T,(order+1)+2)

      for j in (order+1):-1:1

        z[j] = ppp[j,k,l]+2*x2*z[j+1]-z[j+2]

      end

      pp[k,l] = z[1]-x2*z[2]

    end
	end

	p = Array(T,(order+1))

	for l = 1:(order+1)

    z = Array(T,(order+1)+2)

    for k in (order+1):-1:1

      z[k] = pp[k,l]+2*x3*z[k+1]-z[k+2]

    end

    p[l] = z[1]-x3*z[2]

  end

  z = Array(T,(order+1)+2)

  for l in (order+1):-1:1

    z[l] = p[l]+2*x4*z[l+1]-z[l+2]

  end

  y = z[1]-x4*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,5},x::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])
  x4 = normalize_node(x[4],range[:,4])
  x5 = normalize_node(x[5],range[:,5])

  pppp = Array(T,(order+1),(order+1),(order+1),(order+1))

	for m = 1:(order+1)
	  for l = 1:(order+1)
	    for k = 1:(order+1)
		    for j = 1:(order+1)

          z = Array(T,(order+1)+2)

          for i in (order+1):-1:1

            z[i] = weights[i,j,k,l,m]+2*x1*z[i+1]-z[i+2]

          end

          pppp[j,k,l,m] = z[1]-x1*z[2]

        end
	    end
	  end
  end

	ppp = Array(T,(order+1),(order+1),(order+1))

	for m = 1:(order+1)
	  for l = 1:(order+1)
		  for k = 1:(order+1)

        z = Array(T,(order+1)+2)

        for j in (order+1):-1:1

          z[j] = pppp[j,k,l,m]+2*x2*z[j+1]-z[j+2]

        end

        ppp[k,l,m] = z[1]-x2*z[2]

      end
	  end
	end

	pp = Array(T,(order+1),(order+1))

	for m = 1:(order+1)
		for l = 1:(order+1)

      z = Array(T,(order+1)+2)

      for k in (order+1):-1:1

        z[k] = ppp[k,l,m]+2*x3*z[k+1]-z[k+2]

      end

      pp[l,m] = z[1]-x3*z[2]

    end
	end

	p = Array(T,(order+1))

	for m = 1:(order+1)

    z = Array(T,(order+1)+2)

    for l in (order+1):-1:1

      z[l] = pp[l,m]+2*x4*z[l+1]-z[l+2]

    end

    p[m] = z[1]-x4*z[2]

  end

  z = Array(T,(order+1)+2)

  for m in (order+1):-1:1

    z[m] = p[m]+2*x5*z[m+1]-z[m+2]

  end

  y = z[1]-x5*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,6},x::Array{T,1},order::S,range::Array{T,2})

  x1 = normalize_node(x[1],range[:,1])
  x2 = normalize_node(x[2],range[:,2])
  x3 = normalize_node(x[3],range[:,3])
  x4 = normalize_node(x[4],range[:,4])
  x5 = normalize_node(x[5],range[:,5])
  x6 = normalize_node(x[6],range[:,6])

  ppppp = Array(T,(order+1),(order+1),(order+1),(order+1),(order+1))

	for n = 1:(order+1)
	  for m = 1:(order+1)
	    for l = 1:(order+1)
	      for k = 1:(order+1)
		      for j = 1:(order+1)

            z = Array(T,(order+1)+2)

            for i in (order+1):-1:1

              z[i] = weights[i,j,k,l,m,n]+2*x1*z[i+1]-z[i+2]

            end

            ppppp[j,k,l,m,n] = z[1]-x1*z[2]

          end
	      end
	    end
    end
  end

  pppp = Array(T,(order+1),(order+1),(order+1),(order+1))

	for m = 1:(order+1)
	  for l = 1:(order+1)
	    for k = 1:(order+1)
		    for j = 1:(order+1)

          z = Array(T,(order+1)+2)

          for i in (order+1):-1:1

            z[i] = ppppp[i,j,k,l,m]+2*x1*z[i+1]-z[i+2]

          end

          pppp[j,k,l,m] = z[1]-x1*z[2]

        end
	    end
	  end
  end

	ppp = Array(T,(order+1),(order+1),(order+1))

	for m = 1:(order+1)
	  for l = 1:(order+1)
		  for k = 1:(order+1)

        z = Array(T,(order+1)+2)

        for j in (order+1):-1:1

          z[j] = pppp[j,k,l,m]+2*x2*z[j+1]-z[j+2]

        end

        ppp[k,l,m] = z[1]-x2*z[2]

      end
	  end
	end

	pp = Array(T,(order+1),(order+1))

	for m = 1:(order+1)
		for l = 1:(order+1)

      z = Array(T,(order+1)+2)

      for k in (order+1):-1:1

        z[k] = ppp[k,l,m]+2*x3*z[k+1]-z[k+2]

      end

      pp[l,m] = z[1]-x3*z[2]

    end
	end

	p = Array(T,(order+1))

	for m = 1:(order+1)

    z = Array(T,(order+1)+2)

    for l in (order+1):-1:1

      z[l] = pp[l,m]+2*x4*z[l+1]-z[l+2]

    end

    p[m] = z[1]-x4*z[2]

  end

  z = Array(T,(order+1)+2)

  for m in (order+1):-1:1

    z[m] = p[m]+2*x5*z[m+1]-z[m+2]

  end

  y = z[1]-x5*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,1},x::Array{T,1},order::Array{S,1})

  z = Array(T,(order[1]+1)+2)

  for i in (order[1]+1):-1:1

	  z[i] = weights[i]+2*x[1]*z[i+1]-z[i+2]

  end

  y = z[1]-x[1]*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,2},x::Array{T,1},order::Array{S,1})

  p = Array(T,order[2]+1)

	for j = 1:(order[2]+1)

    z = Array(T,(order[1]+1)+2)

    for i in (order[1]+1):-1:1

      z[i] = weights[i,j]+2*x[1]*z[i+1]-z[i+2]

    end

    p[j] = z[1]-x[1]*z[2]

  end

  z = Array(T,(order[2]+1)+2)

  for j in (order[2]+1):-1:1

    z[j] = p[j]+2*x[2]*z[j+1]-z[j+2]

  end

  y = z[1]-x[2]*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,3},x::Array{T,1},order::Array{S,1})

	pp = Array(T,(order[2]+1),(order[3]+1))

	for k  = 1:(order[3]+1)
		for j = 1:(order[2]+1)

      z = Array(T,(order[1]+1)+2)

      for i in (order[1]+1):-1:1

        z[i] = weights[i,j,k]+2*x[1]*z[i+1]-z[i+2]

      end

      pp[j,k] = z[1]-x[1]*z[2]

    end
	end

	p = Array(T,(order[3]+1))

	for k = 1:(order[3]+1)

    z = Array(T,(order[2]+1)+2)

    for j in (order[2]+1):-1:1

      z[j] = pp[j,k]+2*x[2]*z[j+1]-z[j+2]

    end

    p[k] = z[1]-x[2]*z[2]

  end

  z = Array(T,(order[3]+1)+2)

  for k in (order[3]+1):-1:1

    z[k] = p[k]+2*x[3]*z[k+1]-z[k+2]

  end

  y = z[1]-x[3]*z[2]

	return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,4},x::Array{T,1},order::Array{S,1})

	ppp = Array(T,(order[2]+1),(order[3]+1),(order[4]+1))

	for l = 1:(order[4]+1)
	  for k = 1:(order[3]+1)
		  for j = 1:(order[2]+1)

        z = Array(T,(order[1]+1)+2)

        for i in (order[1]+1):-1:1

          z[i] = weights[i,j,k,l]+2*x[1]*z[i+1]-z[i+2]

        end

        ppp[j,k,l] = z[1]-x[1]*z[2]

      end
	  end
	end

	pp = Array(T,(order[3]+1),(order[4]+1))

	for l = 1:(order[4]+1)
		for k = 1:(order[3]+1)

      z = Array(T,(order[2]+1)+2)

      for j in (order[2]+1):-1:1

        z[j] = ppp[j,k,l]+2*x[2]*z[j+1]-z[j+2]

      end

      pp[k,l] = z[1]-x[2]*z[2]

    end
	end

	p = Array(T,(order[4]+1))

	for l = 1:(order[4]+1)

    z = Array(T,(order[3]+1)+2)

    for k in (order[3]+1):-1:1

      z[k] = pp[k,l]+2*x[3]*z[k+1]-z[k+2]

    end

    p[l] = z[1]-x[3]*z[2]

  end

  z = Array(T,(order[4]+1)+2)

  for l in (order[4]+1):-1:1

    z[l] = p[l]+2*x[4]*z[l+1]-z[l+2]

  end

  y = z[1]-x[4]*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,5},x::Array{T,1},order::Array{S,1})

  pppp = Array(T,(order[2]+1),(order[3]+1),(order[4]+1),(order[5]+1))

  for m = 1:(order[5]+1)
	  for l = 1:(order[4]+1)
	    for k = 1:(order[3]+1)
		    for j = 1:(order[2]+1)

          z = Array(T,(order[1]+1)+2)

          for i in (order[1]+1):-1:1

            z[i] = weights[i,j,k,l,m]+2*x[1]*z[i+1]-z[i+2]

          end

          pppp[j,k,l,m] = z[1]-x[1]*z[2]

        end
	    end
	  end
  end

	ppp = Array(T,(order[3]+1),(order[4]+1),(order[5]+1))

	for m = 1:(order[5]+1)
	  for l = 1:(order[4]+1)
		  for k = 1:(order[3]+1)

        z = Array(T,(order[2]+1)+2)

        for j in (order[2]+1):-1:1

          z[j] = pppp[j,k,l,m]+2*x[2]*z[j+1]-z[j+2]

        end

        ppp[k,l,m] = z[1]-x[2]*z[2]

      end
	  end
	end

	pp = Array(T,(order[4]+1),(order[5]+1))

	for m = 1:(order[5]+1)
		for l = 1:(order[4]+1)

      z = Array(T,(order[3]+1)+2)

      for k in (order[3]+1):-1:1

        z[k] = ppp[k,l,m]+2*x[3]*z[k+1]-z[k+2]

      end

      pp[l,m] = z[1]-x[3]*z[2]

    end
	end

	p = Array(T,(order[5]+1))

	for m = 1:(order[5]+1)

    z = Array(T,(order[4]+1)+2)

    for l in (order[4]+1):-1:1

      z[l] = pp[l,m]+2*x[4]*z[l+1]-z[l+2]

    end

    p[m] = z[1]-x[4]*z[2]

  end

  z = Array(T,(order[5]+1)+2)

  for m in (order[5]+1):-1:1

    z[m] = p[m]+2*x[5]*z[m+1]-z[m+2]

  end

  y = z[1]-x[5]*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,6},x::Array{T,1},order::Array{S,1})

  ppppp = Array(T,(order[2]+1),(order[3]+1),(order[4]+1),(order[5]+1),(order[6]+1))

  for n = 1:(order[6]+1)
    for m = 1:(order[5]+1)
	    for l = 1:(order[4]+1)
	      for k = 1:(order[3]+1)
		      for j = 1:(order[2]+1)

            z = Array(T,(order[1]+1)+2)

            for i in (order[1]+1):-1:1

              z[i] = weights[i,j,k,l,m,n]+2*x1*z[i+1]-z[i+2]

            end

            ppppp[j,k,l,m,n] = z[1]-x1*z[2]

          end
	      end
	    end
    end
  end

  pppp = Array(T,(order[3]+1),(order[4]+1),(order[5]+1),(order[6]+1))

  for m = 1:(order[6]+1)
	  for l = 1:(order[5]+1)
	    for k = 1:(order[4]+1)
		    for j = 1:(order[3]+1)

          z = Array(T,(order[2]+1)+2)

          for i in (order[2]+1):-1:1

            z[i] = ppppp[i,j,k,l,m]+2*x1*z[i+1]-z[i+2]

          end

          pppp[j,k,l,m] = z[1]-x1*z[2]

        end
	    end
	  end
  end

	ppp = Array(T,(order[4]+1),(order[5]+1),(order[6]+1))

	for m = 1:(order[6]+1)
	  for l = 1:(order[5]+1)
		  for k = 1:(order[4]+1)

        z = Array(T,(order[3]+1)+2)

        for j in (order[3]+1):-1:1

          z[j] = pppp[j,k,l,m]+2*x2*z[j+1]-z[j+2]

        end

        ppp[k,l,m] = z[1]-x2*z[2]

      end
	  end
	end

	pp = Array(T,(order[5]+1),(order[6]+1))

	for m = 1:(order[6]+1)
		for l = 1:(order[5]+1)

      z = Array(T,(order[4]+1)+2)

      for k in (order[4]+1):-1:1

        z[k] = ppp[k,l,m]+2*x3*z[k+1]-z[k+2]

      end

      pp[l,m] = z[1]-x3*z[2]

    end
	end

	p = Array(T,(order[6]+1))

	for m = 1:(order[6]+1)

    z = Array(T,(order[5]+1)+2)

    for l in (order[5]+1):-1:1

      z[l] = pp[l,m]+2*x4*z[l+1]-z[l+2]

    end

    p[m] = z[1]-x4*z[2]

  end

  z = Array(T,(order[6]+1)+2)

  for m in (order[6]+1):-1:1

    z[m] = p[m]+2*x5*z[m+1]-z[m+2]

  end

  y = z[1]-x5*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,1},x::Array{T,1},order::S)

  z = Array(T,(order+1)+2)

  for i in (order+1):-1:1

	  z[i] = weights[i]+2*x[1]*z[i+1]-z[i+2]

  end

  y = z[1]-x[1]*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,2},x::Array{T,1},order::S)

  p = Array(T,order+1)

	for j = 1:(order+1)

    z = Array(T,(order+1)+2)

    for i in (order+1):-1:1

      z[i] = weights[i,j]+2*x[1]*z[i+1]-z[i+2]

    end

    p[j] = z[1]-x[1]*z[2]

  end

  z = Array(T,(order+1)+2)

  for j in (order+1):-1:1

    z[j] = p[j]+2*x[2]*z[j+1]-z[j+2]

  end

  y = z[1]-x[2]*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,3},x::Array{T,1},order::S)

	pp = Array(T,(order+1),(order+1))

	for k  = 1:(order+1)
		for j = 1:(order+1)

      z = Array(T,(order+1)+2)

      for i in (order+1):-1:1

        z[i] = weights[i,j,k]+2*x[1]*z[i+1]-z[i+2]

      end

      pp[j,k] = z[1]-x[1]*z[2]

    end
	end

	p = Array(T,(order+1))

	for k = 1:(order+1)

    z = Array(T,(order+1)+2)

    for j in (order+1):-1:1

      z[j] = pp[j,k]+2*x[2]*z[j+1]-z[j+2]

    end

    p[k] = z[1]-x[2]*z[2]

  end

  z = Array(T,(order+1)+2)

  for k in (order+1):-1:1

    z[k] = p[k]+2*x[3]*z[k+1]-z[k+2]

  end

  y = z[1]-x[3]*z[2]

	return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,4},x::Array{T,1},order::S)

	ppp = Array(T,(order+1),(order+1),(order+1))

	for l = 1:(order+1)
	  for k = 1:(order+1)
		  for j = 1:(order+1)

        z = Array(T,(order+1)+2)

        for i in (order+1):-1:1

          z[i] = weights[i,j,k,l]+2*x[1]*z[i+1]-z[i+2]

        end

        ppp[j,k,l] = z[1]-x[1]*z[2]

      end
	  end
	end

	pp = Array(T,(order+1),(order+1))

	for l = 1:(order+1)
		for k = 1:(order+1)

      z = Array(T,(order+1)+2)

      for j in (order+1):-1:1

        z[j] = ppp[j,k,l]+2*x[2]*z[j+1]-z[j+2]

      end

      pp[k,l] = z[1]-x[2]*z[2]

    end
	end

	p = Array(T,(order+1))

	for l = 1:(order+1)

    z = Array(T,(order+1)+2)

    for k in (order+1):-1:1

      z[k] = pp[k,l]+2*x[3]*z[k+1]-z[k+2]

    end

    p[l] = z[1]-x[3]*z[2]

  end

  z = Array(T,(order+1)+2)

  for l in (order+1):-1:1

    z[l] = p[l]+2*x[4]*z[l+1]-z[l+2]

  end

  y = z[1]-x[4]*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,5},x::Array{T,1},order::S)

  pppp = Array(T,(order+1),(order+1),(order+1),(order+1))

	for m = 1:(order+1)
	  for l = 1:(order+1)
	    for k = 1:(order+1)
		    for j = 1:(order+1)

          z = Array(T,(order+1)+2)

          for i in (order+1):-1:1

            z[i] = weights[i,j,k,l,m]+2*x[1]*z[i+1]-z[i+2]

          end

          pppp[j,k,l,m] = z[1]-x[1]*z[2]

        end
	    end
	  end
  end

	ppp = Array(T,(order+1),(order+1),(order+1))

	for m = 1:(order+1)
	  for l = 1:(order+1)
		  for k = 1:(order+1)

        z = Array(T,(order+1)+2)

        for j in (order+1):-1:1

          z[j] = pppp[j,k,l,m]+2*x[2]*z[j+1]-z[j+2]

        end

        ppp[k,l,m] = z[1]-x[2]*z[2]

      end
	  end
	end

	pp = Array(T,(order+1),(order+1))

	for m = 1:(order+1)
		for l = 1:(order+1)

      z = Array(T,(order+1)+2)

      for k in (order+1):-1:1

        z[k] = ppp[k,l,m]+2*x[3]*z[k+1]-z[k+2]

      end

      pp[l,m] = z[1]-x[3]*z[2]

    end
	end

	p = Array(T,(order+1))

	for m = 1:(order+1)

    z = Array(T,(order+1)+2)

    for l in (order+1):-1:1

      z[l] = pp[l,m]+2*x[4]*z[l+1]-z[l+2]

    end

    p[m] = z[1]-x[4]*z[2]

  end

  z = Array(T,(order+1)+2)

  for m in (order+1):-1:1

    z[m] = p[m]+2*x[5]*z[m+1]-z[m+2]

  end

  y = z[1]-x[5]*z[2]

  return y

end

function clenshaw_evaluate{T<:AbstractFloat,S<:Integer}(weights::Array{T,6},x::Array{T,1},order::S)

  ppppp = Array(T,(order+1),(order+1),(order+1),(order+1),(order+1))

	for n = 1:(order+1)
	  for m = 1:(order+1)
	    for l = 1:(order+1)
	      for k = 1:(order+1)
		      for j = 1:(order+1)

            z = Array(T,(order+1)+2)

            for i in (order+1):-1:1

              z[i] = weights[i,j,k,l,m,n]+2*x1*z[i+1]-z[i+2]

            end

            ppppp[j,k,l,m,n] = z[1]-x1*z[2]

          end
	      end
	    end
    end
  end

  pppp = Array(T,(order+1),(order+1),(order+1),(order+1))

	for m = 1:(order+1)
	  for l = 1:(order+1)
	    for k = 1:(order+1)
		    for j = 1:(order+1)

          z = Array(T,(order+1)+2)

          for i in (order+1):-1:1

            z[i] = ppppp[i,j,k,l,m]+2*x1*z[i+1]-z[i+2]

          end

          pppp[j,k,l,m] = z[1]-x1*z[2]

        end
	    end
	  end
  end

	ppp = Array(T,(order+1),(order+1),(order+1))

	for m = 1:(order+1)
	  for l = 1:(order+1)
		  for k = 1:(order+1)

        z = Array(T,(order+1)+2)

        for j in (order+1):-1:1

          z[j] = pppp[j,k,l,m]+2*x2*z[j+1]-z[j+2]

        end

        ppp[k,l,m] = z[1]-x2*z[2]

      end
	  end
	end

	pp = Array(T,(order+1),(order+1))

	for m = 1:(order+1)
		for l = 1:(order+1)

      z = Array(T,(order+1)+2)

      for k in (order+1):-1:1

        z[k] = ppp[k,l,m]+2*x3*z[k+1]-z[k+2]

      end

      pp[l,m] = z[1]-x3*z[2]

    end
	end

	p = Array(T,(order+1))

	for m = 1:(order+1)

    z = Array(T,(order+1)+2)

    for l in (order+1):-1:1

      z[l] = pp[l,m]+2*x4*z[l+1]-z[l+2]

    end

    p[m] = z[1]-x4*z[2]

  end

  z = Array(T,(order+1)+2)

  for m in (order+1):-1:1

    z[m] = p[m]+2*x5*z[m+1]-z[m+2]

  end

  y = z[1]-x5*z[2]

  return y

end
