




import Base: +, -, *, /
import Base.LinAlg.norm


typealias Point{N,T} NTuple{N,T}

typealias Point1{T} Tuple{T}
typealias Point2{T} Tuple{T,T}
typealias Point3{T} Tuple{T,T,T}
typealias Point4{T} Tuple{T,T,T,T}

typealias Point1f Point1{Float64}
typealias Point2f Point2{Float64}
typealias Point3f Point3{Float64}
typealias Point4f Point4{Float64}
typealias Point1i Point1{Int}
typealias Point2i Point2{Int}
typealias Point3i Point3{Int}
typealias Point4i Point4{Int}



# Point{T}(x::T, y::T, z::T) = (x,y,z)

precompile(Point, (Float64,Float64,Float64))


for op in (:+, :-, :*, :/, :.+, :.-, :.*, :./)
    @eval begin
        Base.$op{T}(x::Point2{T}, y::Point2{T}) = $op(x[1], y[1]), $op(x[2], y[2])
        Base.$op{T}(x::Point3{T}, y::Point3{T}) = $op(x[1], y[1]), $op(x[2], y[2]), $op(x[3], y[3])
        Base.$op{T}(x::Point2{T}, a::Number) = $op(x[1], a), $op(x[2], a)
        Base.$op{T}(x::Point3{T}, a::Number) = $op(x[1], a), $op(x[2], a), $op(x[3], a)
        Base.$op{T}(a::Number, x::Point2{T}) = $op(x,a)
        Base.$op{T}(a::Number, x::Point3{T}) = $op(x,a)
        precompile($op, (Point2f, Point2f))
        precompile($op, (Point3f, Point3f))
        precompile($op, (Point2f, Float64))
        precompile($op, (Point3f, Float64))
    end
end


norm(x::Point2f) = sqrt(x[1]*x[1] + x[2]*x[2])
norm(x::Point3f) = sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3])

precompile(norm, (Point2f,))
precompile(norm, (Point3f,))

function Base.(:*){T}(A::Array{T,2}, x::Point3{T})
    A[1,1]*x[1] + A[1,2]*x[2] + A[1,3]*x[3], A[2,1]*x[1] + A[2,2]*x[2] + A[2,3]*x[3], A[3,1]*x[1] + A[3,2]*x[2] + A[3,3]*x[3]
end


precompile(*, (Array{Float64,2}, Point3f))




@doc """
`dist(u, v)` computes the Euclidean distance betweens points `u` and `v`.
""" ->
function dist{T}(u::Point3{T}, v::Point3{T})
    dx = u[1]-v[1]
    dy = u[2]-v[2]
    dz = u[3]-v[3]
    sqrt(dx*dx + dy*dy + dz*dz)
end

dist{T}(u::AbstractVector{Point3{T}}, v::Point3{T}) = T[dist(u[i],v) for i in 1:length(u)]
dist{T}(u::Point3{T}, v::AbstractVector{Point3{T}}) = T[dist(u,v[i]) for i in 1:length(v)]

function dist!{T}(dst::AbstractMatrix{T}, u::AbstractVector{Point3{T}}, v::AbstractVector{Point3{T}})
    @assert size(dst,1) == length(u) && size(dst,2) == length(v)
    for i in 1:length(v)
        for j in 1:length(u)
            dst[j,i] = dist(u[j],v[i])
        end
    end
end

function dist{T}(u::AbstractVector{Point3{T}}, v::AbstractVector{Point3{T}})
    dst = Array{T}(length(u),length(v))
    dist!(dst,u,v)
    dst
end

function dist!{T}(dst::AbstractMatrix{T}, u::AbstractVector{Point3{T}})
    n = length(u)
    @assert size(dst,1) == size(dst,2) == n
    for i in 1:n
        for j in 1:i
            dst[j,i] = dist(u[j],u[i])
        end
    end
end

function dist{T}(u::AbstractVector{Point3{T}})
    dst = Array{T}(length(u),length(u))
    dist!(dst,u)
    Symmetric(dst)
end


function randp(d::Int, dims...)
    dst = Array{NTuple{d,Float64}}(dims...)
    for i in eachindex(dst)
        dst[i] = ntuple(i->rand(),d)
    end
    dst
end


precompile(dist , (Point3f, Point3f))
precompile(dist , (Vector{Point3f}, Point3f))
precompile(dist , (Point3f, Vector{Point3f}))
precompile(dist!, (Matrix{Point3f}, Vector{Point3f}, Vector{Point3f}))
precompile(dist , (Vector{Point3f}, Vector{Point3f}))
precompile(dist!, (Matrix{Point3f}, Vector{Point3f}))
precompile(dist , (Vector{Point3f},))



