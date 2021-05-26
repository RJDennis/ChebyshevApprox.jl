abstract type ChebInterp end

struct ChebPoly{T<:AbstractFloat, S<:Integer, N}

  weights::Array{T,N}
  order::Union{S,Array{S,1}}
  domain::Union{Array{T,1},Array{T,2}}

end

struct ChebInterpRoots{T<:AbstractFloat, S<:Integer, N} <: ChebInterp

  data::Array{T,N}
  nodes::NTuple{N,Array{T,1}}
  order::Union{S,Array{S,1}}
  domain::Union{Array{T,1},Array{T,2}}

end

struct ChebInterpExtrema{T<:AbstractFloat, S<:Integer, N} <: ChebInterp

  data::Array{T,N}
  nodes::NTuple{N,Array{T,1}}
  order::Union{S,Array{S,1}}
  domain::Union{Array{T,1},Array{T,2}}

end

struct ChebInterpExtended{T<:AbstractFloat, S<:Integer, N} <: ChebInterp

  data::Array{T,N}
  nodes::NTuple{N,Array{T,1}}
  order::Union{S,Array{S,1}}
  domain::Union{Array{T,1},Array{T,2}}

end

struct ChebInterpVertesi{T<:AbstractFloat, S<:Integer, N} <: ChebInterp

  data::Array{T,N}
  nodes::NTuple{N,Array{T,1}}
  order::Union{S,Array{S,1}}
  domain::Union{Array{T,1},Array{T,2}}

end
