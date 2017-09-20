function clenshaw_evaluate(weights::Array{T,1},x::Array{T,1},order::Array{S,1},domain=[ones(T,1);-ones(T,1)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain)

  z = zeros(T,(order[1]+1)+2)

  for i in (order[1]+1):-1:1

	  z[i] = weights[i]+2*x1*z[i+1]-z[i+2]

  end

  y = z[1]-x1*z[2]

  return y

end

function clenshaw_evaluate(weights::Array{T,2},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,2);-ones(T,1,2)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])

  p = zeros(T,order[2]+1)

	for j = 1:(order[2]+1)

    z = zeros(T,(order[1]+1)+2)

    for i in (order[1]+1):-1:1

      z[i] = weights[i,j]+2*x1*z[i+1]-z[i+2]

    end

    p[j] = z[1]-x1*z[2]

  end

  z = zeros(T,(order[2]+1)+2)

  for j in (order[2]+1):-1:1

    z[j] = p[j]+2*x2*z[j+1]-z[j+2]

  end

  y = z[1]-x2*z[2]

  return y

end

function clenshaw_evaluate(weights::Array{T,3},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,3);-ones(T,1,3)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])
  x3 = normalize_node(x[3],domain[:,3])

	pp = zeros(T,(order[2]+1),(order[3]+1))

	for k  = 1:(order[3]+1)
		for j = 1:(order[2]+1)

      z = zeros(T,(order[1]+1)+2)

      for i in (order[1]+1):-1:1

        z[i] = weights[i,j,k]+2*x1*z[i+1]-z[i+2]

      end

      pp[j,k] = z[1]-x1*z[2]

    end
	end

	p = zeros(T,(order[3]+1))

	for k = 1:(order[3]+1)

    z = zeros(T,(order[2]+1)+2)

    for j in (order[2]+1):-1:1

      z[j] = pp[j,k]+2*x2*z[j+1]-z[j+2]

    end

    p[k] = z[1]-x2*z[2]

  end

  z = zeros(T,(order[3]+1)+2)

  for k in (order[3]+1):-1:1

    z[k] = p[k]+2*x3*z[k+1]-z[k+2]

  end

  y = z[1]-x3*z[2]

	return y

end

function clenshaw_evaluate(weights::Array{T,4},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,4);-ones(T,1,4)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])
  x3 = normalize_node(x[3],domain[:,3])
  x4 = normalize_node(x[4],domain[:,4])

	ppp = zeros(T,(order[2]+1),(order[3]+1),(order[4]+1))

	for l = 1:(order[4]+1)
	  for k = 1:(order[3]+1)
		  for j = 1:(order[2]+1)

        z = zeros(T,(order[1]+1)+2)

        for i in (order[1]+1):-1:1

          z[i] = weights[i,j,k,l]+2*x1*z[i+1]-z[i+2]

        end

        ppp[j,k,l] = z[1]-x1*z[2]

      end
	  end
	end

	pp = zeros(T,(order[3]+1),(order[4]+1))

	for l = 1:(order[4]+1)
		for k = 1:(order[3]+1)

      z = zeros(T,(order[2]+1)+2)

      for j in (order[2]+1):-1:1

        z[j] = ppp[j,k,l]+2*x2*z[j+1]-z[j+2]

      end

      pp[k,l] = z[1]-x2*z[2]

    end
	end

	p = zeros(T,(order[4]+1))

	for l = 1:(order[4]+1)

    z = zeros(T,(order[3]+1)+2)

    for k in (order[3]+1):-1:1

      z[k] = pp[k,l]+2*x3*z[k+1]-z[k+2]

    end

    p[l] = z[1]-x3*z[2]

  end

  z = zeros(T,(order[4]+1)+2)

  for l in (order[4]+1):-1:1

    z[l] = p[l]+2*x4*z[l+1]-z[l+2]

  end

  y = z[1]-x4*z[2]

  return y

end

function clenshaw_evaluate(weights::Array{T,5},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,5);-ones(T,1,5)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])
  x3 = normalize_node(x[3],domain[:,3])
  x4 = normalize_node(x[4],domain[:,4])
  x5 = normalize_node(x[5],domain[:,5])

  pppp = zeros(T,(order[2]+1),(order[3]+1),(order[4]+1),(order[5]+1))

  for m = 1:(order[5]+1)
	  for l = 1:(order[4]+1)
	    for k = 1:(order[3]+1)
		    for j = 1:(order[2]+1)

          z = zeros(T,(order[1]+1)+2)

          for i in (order[1]+1):-1:1

            z[i] = weights[i,j,k,l,m]+2*x1*z[i+1]-z[i+2]

          end

          pppp[j,k,l,m] = z[1]-x1*z[2]

        end
	    end
	  end
  end

	ppp = zeros(T,(order[3]+1),(order[4]+1),(order[5]+1))

	for m = 1:(order[5]+1)
	  for l = 1:(order[4]+1)
		  for k = 1:(order[3]+1)

        z = zeros(T,(order[2]+1)+2)

        for j in (order[2]+1):-1:1

          z[j] = pppp[j,k,l,m]+2*x2*z[j+1]-z[j+2]

        end

        ppp[k,l,m] = z[1]-x2*z[2]

      end
	  end
	end

	pp = zeros(T,(order[4]+1),(order[5]+1))

	for m = 1:(order[5]+1)
		for l = 1:(order[4]+1)

      z = zeros(T,(order[3]+1)+2)

      for k in (order[3]+1):-1:1

        z[k] = ppp[k,l,m]+2*x3*z[k+1]-z[k+2]

      end

      pp[l,m] = z[1]-x3*z[2]

    end
	end

	p = zeros(T,(order[5]+1))

	for m = 1:(order[5]+1)

    z = zeros(T,(order[4]+1)+2)

    for l in (order[4]+1):-1:1

      z[l] = pp[l,m]+2*x4*z[l+1]-z[l+2]

    end

    p[m] = z[1]-x4*z[2]

  end

  z = zeros(T,(order[5]+1)+2)

  for m in (order[5]+1):-1:1

    z[m] = p[m]+2*x5*z[m+1]-z[m+2]

  end

  y = z[1]-x5*z[2]

  return y

end

function clenshaw_evaluate(weights::Array{T,6},x::Array{T,1},order::Array{S,1},domain=[ones(T,1,6);-ones(T,1,6)]) where {T<:AbstractFloat,S<:Integer}

  x1 = normalize_node(x[1],domain[:,1])
  x2 = normalize_node(x[2],domain[:,2])
  x3 = normalize_node(x[3],domain[:,3])
  x4 = normalize_node(x[4],domain[:,4])
  x5 = normalize_node(x[5],domain[:,5])
  x6 = normalize_node(x[6],domain[:,6])

  ppppp = zeros(T,(order[2]+1),(order[3]+1),(order[4]+1),(order[5]+1),(order[6]+1))

  for n = 1:(order[6]+1)
    for m = 1:(order[5]+1)
	    for l = 1:(order[4]+1)
	      for k = 1:(order[3]+1)
		      for j = 1:(order[2]+1)

            z = zeros(T,(order[1]+1)+2)

            for i in (order[1]+1):-1:1

              z[i] = weights[i,j,k,l,m,n]+2*x1*z[i+1]-z[i+2]

            end

            ppppp[j,k,l,m,n] = z[1]-x1*z[2]

          end
	      end
	    end
    end
  end

  pppp = zeros(T,(order[3]+1),(order[4]+1),(order[5]+1),(order[6]+1))

  for m = 1:(order[6]+1)
	  for l = 1:(order[5]+1)
	    for k = 1:(order[4]+1)
		    for j = 1:(order[3]+1)

          z = zeros(T,(order[2]+1)+2)

          for i in (order[2]+1):-1:1

            z[i] = ppppp[i,j,k,l,m]+2*x2*z[i+1]-z[i+2]

          end

          pppp[j,k,l,m] = z[1]-x2*z[2]

        end
	    end
	  end
  end

	ppp = zeros(T,(order[4]+1),(order[5]+1),(order[6]+1))

	for m = 1:(order[6]+1)
	  for l = 1:(order[5]+1)
		  for k = 1:(order[4]+1)

        z = zeros(T,(order[3]+1)+2)

        for j in (order[3]+1):-1:1

          z[j] = pppp[j,k,l,m]+2*x3*z[j+1]-z[j+2]

        end

        ppp[k,l,m] = z[1]-x3*z[2]

      end
	  end
	end

	pp = zeros(T,(order[5]+1),(order[6]+1))

	for m = 1:(order[6]+1)
		for l = 1:(order[5]+1)

      z = zeros(T,(order[4]+1)+2)

      for k in (order[4]+1):-1:1

        z[k] = ppp[k,l,m]+2*x4*z[k+1]-z[k+2]

      end

      pp[l,m] = z[1]-x4*z[2]

    end
	end

	p = zeros(T,(order[6]+1))

	for m = 1:(order[6]+1)

    z = zeros(T,(order[5]+1)+2)

    for l in (order[5]+1):-1:1

      z[l] = pp[l,m]+2*x5*z[l+1]-z[l+2]

    end

    p[m] = z[1]-x5*z[2]

  end

  z = zeros(T,(order[6]+1)+2)

  for m in (order[6]+1):-1:1

    z[m] = p[m]+2*x6*z[m+1]-z[m+2]

  end

  y = z[1]-x6*z[2]

  return y

end
