abstract type ChebPolynomial end
abstract type ChebInterp end

struct ChebPolyTensor{T<:AbstractFloat, S<:Integer, N} <: ChebPolynomial

    weights::Array{T,N}
    order::Array{S,1}
    domain::Union{Array{T,1},Array{T,2}}

end

struct ChebPolyComplete{T<:AbstractFloat, S<:Integer, N} <: ChebPolynomial

    weights::Array{T,N}
    order::S
    domain::Union{Array{T,1},Array{T,2}}

end

struct ChebInterpTensor{T<:AbstractFloat, S<:Integer, N} <: ChebInterp

    data::Array{T,N}
    nodes::NTuple{N,Array{T,1}}
    order::Array{S,1}
    domain::Union{Array{T,1},Array{T,2}}

end

struct ChebInterpComplete{T<:AbstractFloat, S<:Integer, N} <: ChebInterp

    data::Array{T,N}
    nodes::NTuple{N,Array{T,1}}
    order::S
    domain::Union{Array{T,1},Array{T,2}}

end
