


doc"""
Spherically symmetric function 

$f(r,\theta,\phi) = u(r)$

where $u(r)$ is a user-defined radial function.

The function vanishes for $r>r_{max}$.

f = RadialFunc(u; R=(0.,0.,0.), rmax=Inf)

`u` - a user-defined radial function

`R` - origin of spherical coordinate system
"""
type RadialFunc{T} <: ScalarFunc3d{Float64}
    u::T
    R::Point3{Float64}
    rmax::Float64
end

RadialFunc(u; R=(0.,0.,0.), rmax=Inf) = RadialFunc(u,R,rmax)

RadialFunc(r::AbstractVector, u::AbstractVector; R=(0.,0.,0.),
        rmax=Inf, npts::Int=100) =
    RadialFunc(Spline1D(r,u,npts),R,rmax)


typealias RadialFuncSpl RadialFunc{Spline1D{Float64,Vector{Float64}}}

function Base.writemime{T}(io::IO, mime::MIME"text/plain", f::RadialFunc{T})
    write(io, summary(f), "\n")
    write(io, "u(r) = ")
    writemime(io, mime, f.u)
    write(io, "\nR = ")
    writemime(io, mime, f.R)
    write(io, "\nrmax = ")
    writemime(io, mime, f.rmax)
end


function Base.call(f::RadialFunc, x::Float64, y::Float64, z::Float64)
    x0,y0,z0 = f.R
    r = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)
    r > f.rmax ? 0. : f.u(r)
end

Base.call(f::RadialFunc, p::Point3{Float64}) = f(p[1],p[2],p[3])
@vectorize2 Base.call RadialFunc Point3{Float64} Float64



function grad(f::RadialFunc{Spline1D}, x::Real, y::Real, z::Real)
    x0,y0,z0 = f.R
    r = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)
    dfdr = r > f.rmax ? 0. : deriv1(f.u,r)
    (dfdr,0.,0.)
end

grad(f::RadialFunc, p::Point3{Float64}) = grad(f,p[1],p[2],p[3])
@vectorize2 grad RadialFunc Point3{Float64} Point3{Float64}


∇ = grad


function laplacian(f::RadialFunc{Spline1D}, x::Real, y::Real, z::Real)
    x0,y0,z0 = f.R
    r = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)
    _,df,df2 = eval_deriv12(f.u,r)
    r > f.rmax ? 0. : 2*df/r + df2
end

laplacian(f::RadialFunc, p::Point3{Float64}) = laplacian(f,p[1],p[2],p[3])
@vectorize2 laplacian RadialFunc Point3{Float64} Float64



function setposition(f::RadialFunc, R::Point3{Float64})
    f.R = R
end


