abstract type ChebyshevPolynomial end

struct ChebyshevPolyTensor{T<:AbstractFloat, S<:Integer, N} <: ChebyshevPolynomial

    weights::Array{T,N}
    order::Array{S,1}
    domain::Union{Array{T,1},Array{T,2}}

end

struct ChebyshevPolyComplete{T<:AbstractFloat, S<:Integer, N} <: ChebyshevPolynomial

    weights::Array{T,N}
    order::S
    domain::Union{Array{T,1},Array{T,2}}

end
