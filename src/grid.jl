
abstract AbstractGrid{T,N} <: AbstractArray{T,N}


abstract AbstractGrid1d{T<:Number,N} <: AbstractGrid{T,N}
abstract AbstractGrid2d{T<:Number,N} <: AbstractGrid{Point2{T},N}
abstract AbstractGrid3d{T<:Number,N} <: AbstractGrid{Point3{T},N}


Base.length(g::AbstractGrid) = length(g.c)
Base.eltype{T}(g::AbstractGrid{T}) = T
Base.ndims(g::AbstractGrid) = ndims(g.c)
Base.size(g::AbstractGrid) = size(g.c)
Base.size(g::AbstractGrid, i::Integer) = size(g.c, i)
Base.getindex(g::AbstractGrid, i...) = g.c[i...]
function Base.setindex!(g::AbstractGrid, value, i...)
    g.c[i...] = value
end
Base.transpose(g::AbstractGrid) = transpose(g.c)
Base.ctranspose(g::AbstractGrid) = ctranspose(g.c)
Base.eachindex(g::AbstractGrid) = eachindex(g.c)
Base.linearindexing{T<:AbstractGrid}(::Type{T}) = Base.LinearFast()
# Base.convert(T::Type{Array},g::AbstractGrid) = convert(T,g.c)


integrate(g::AbstractGrid, f) = integrate(g, f(g))
integrate(g::AbstractGrid, f::AbstractVector) = vecdot(g.w, f)

function integrate(g::AbstractGrid, f::AbstractArray)
    n,m = length(f),length(g)
    shp = (m,Int(n/m))
    At_mul_B(reshape(f, shp), vec(g.w))
end

# integrate(g::AbstractGrid, f::AbstractArray) = sum(g.w .* f, 1:ndims(g.w)) |> squeezeall

# integrate{T}(g::AbstractGrid{T,1}, f::AbstractMatrix) = At_mul_B(f, g.w)


∫ = integrate


function Base.writemime{T}(io::IO, mime::MIME"text/plain", g::AbstractGrid{T,1})
    write(io, summary(g), "\n")
    write(io, "[coords weights] =\n")
    writemime(io, mime, [g.c g.w])
end


function Base.writemime{T}(io::IO, mime::MIME"text/plain", g::AbstractGrid{T,2})
    write(io, summary(g), "\n")
    write(io, "coords =\n")
    writemime(io, mime, g.c)
end



function setposition(g::AbstractGrid3d, R::Point3{Float64})
    for i in 1:length(g)
        g[i] += R - g.R
    end
end

type LogarithmicGrid{T,A<:AbstractVector} <: AbstractGrid1d{T,1}
    c::A  # coords
    w::A  # weights
    dx::T  # spacing in log space
    function LogarithmicGrid(c::AbstractVector{T},w::AbstractVector{T},dx::T)
        @assert size(c) == size(w)
        new(c,w,dx)
    end
end

LogarithmicGrid{T}(c::AbstractVector{T},w::AbstractVector{T},dx::T) = 
    LogarithmicGrid{eltype(c),typeof(c)}(c,w,dx)



function loggrid(start::Real, stop::Real, n::Int)
    x = collect(linspace(start, stop, n))
    dx = x[2] - x[1]
    r = exp(x)
    w = r.*dx
    LogarithmicGrid(r, w, dx)
end

# returns a logarithmic grid over same range but with different number of points
loggrid(r::LogarithmicGrid, n::Int) = loggrid(log(r[1]),log(r[end]),n)







type UniformGrid{T,A<:AbstractVector} <: AbstractGrid1d{T,1}
    c::A  # coordinates
    w::A  # weights
    dx::T  # spacing
    function UniformGrid(c::AbstractVector{T},w::AbstractVector{T},dx::T)
        @assert size(c) == size(w)
        new(c,w,dx)
    end
end

UniformGrid{T}(c::AbstractVector{T},w::AbstractVector{T},dx::T) = 
    UniformGrid{eltype(c),typeof(c)}(c,w,dx)



function lingrid(start::Real, stop::Real, n::Int)
    x = collect(linspace(start, stop, n))
    dx = x[2]-x[1]
    UniformGrid(x, fill(dx, length(x)), dx)
end


function simpsgrid(start::Real, stop::Real, n::Int)
    @assert isodd(n) "Simpson's rule requires an odd number of grid points"
    x = collect(linspace(start, stop, n))
    dx = x[2] - x[1]
    w = fill(dx/3, n)
    w[2:2:end] .*= 4
    w[3:2:end-1] .*= 2
    UniformGrid(x, w, dx)
end




type AtomicRadialGrid{T,A<:AbstractVector} <: AbstractGrid1d{T,1}
    c::A  # coords
    w::A  # weights
    rmax::T
    function AtomicRadialGrid(c::AbstractVector{T},w::AbstractVector{T},rmax::T)
        @assert size(c) == size(w)
        new(c,w,rmax)
    end
end


AtomicRadialGrid{T}(c::AbstractVector{T},w::AbstractVector{T},rmax::T) = 
    AtomicRadialGrid{eltype(c),typeof(c)}(c,w,rmax)



function atomicradgrid(n::Int, rmax::Real)
    s = collect(1:n)
    r = rmax * log(1 - (s/(n+1.)).^2) ./ log(1 - (n/(n+1))^2)
    w = -2rmax*s ./ ( (1+n-s) .* (1+n+s) .* log(1 - (n/(n+1.))^2) )
    AtomicRadialGrid(r, w, rmax)
end


loggrid(r::AtomicRadialGrid, n::Int) = loggrid(log(r[1]),log(r[end]),n)


type AngularGrid{T,A<:AbstractVector,B<:AbstractVector} <: AbstractGrid3d{T,1}
    c::A  # coords
    w::B  # weights
    function AngularGrid(c::AbstractVector{Point3{T}},w::AbstractVector{T})
        @assert size(c) == size(w)
        new(c,w)
    end
end

AngularGrid{T}(c::AbstractVector{Point3{T}},w::AbstractVector{T}) = 
    AngularGrid{eltype(w),typeof(c),typeof(w)}(c,w)



typealias AngularGridf AngularGrid{Float64,Vector{Tuple{Float64,Float64,Float64}},Vector{Float64}}



lebgrid(n::Int) = lebedev[n]::AngularGridf




# function expand(grid, basis, f)
#     S = basis.data' * (basis .* x.w)
#     (integrate(x, basis.*h)' / S)'
# end


# function expand(grid, basis, f)
#     S = Symmetric(LinAlg.BLAS.syrk('U', 'T', basis.data .* sqrt(grid.w)))
#     (integrate(grid, basis.*f)' / S)'
# end

# LinAlg.BLAS.herk('U', 'C', 1.0, mat(basis.data .* sqrt(g.w)))

function qrgramschmidt(grid, basis)
    Q,R = qr(basis .* sqrt(grid.w));
    Q ./= sqrt(grid.w)
    Q
end








