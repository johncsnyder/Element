



doc"""
Multipole function 

$f_{lm}(r,\theta,\phi) = u(r) Y_{lm}(\theta,\phi)$,

where $u(r)$ is a user-defined radial function.

If $r_{cut}<r<r_{max}$, then the long-range analytical
form $q Y_{lm}(\theta,\phi)/r^{l+1}$ is used, where
`q` is the multipole moment.

The function vanishes for $r>r_{max}$.

f = MultipoleFunc(u; l=0, m=0, R=(0.,0.,0.), rcut=Inf, rmax=Inf, q=0.)

`u` - a user-defined radial function

e.g.
```
u = r -> exp(-r^2)
u = Spline1D(r, exp(-r.^2))
```

`l`,`m` - angular quantum numbers

`R` - origin of spherical coordinate system
"""
type MultipoleFunc{T} <: ScalarFunc3d{Float64}
    u::T
    l::Int
    m::Int
    R::Point3{Float64}
    rcut::Float64
    rmax::Float64
    q::Float64
end

MultipoleFunc(u; l=0, m=0, R=(0.,0.,0.), rcut=Inf, rmax=Inf, q=0.) = MultipoleFunc{typeof(u)}(u,l,m,R,rcut,rmax,q)

MultipoleFunc(r::AbstractVector, u::AbstractVector; l=0, m=0, R=(0.,0.,0.), rcut=Inf, 
        rmax=Inf, q=0., npts::Int=100) =
    MultipoleFunc(Spline1D(r,u,npts),l,m,R,rcut,rmax,q)



function Base.call(f::MultipoleFunc, x::Float64, y::Float64, z::Float64)
    r,θ,ϕ = cart2sph(x,y,z,f.R)
    r > f.rmax ? 0.0 : ( (r > f.rcut ? f.q/r^(f.l+1) : f.u(r)) * realsphharm(f.l,f.m,θ,ϕ))
end

Base.call(f::MultipoleFunc, p::Point3{Float64}) = f(p[1],p[2],p[3])

@vectorize2 Base.call MultipoleFunc Point3{Float64} Float64




doc"""
`enumeratelm(lmax)` lists all pairs `(l,m)` for `l=0,...,lmax` and `m=-l,...,l`.
"""
function enumeratelm(lmax)
    q = Point2{Int}[]
    for l in 0:lmax
        for m in -l:l
            push!(q, (l,m))
        end
    end
    q
end


typealias MultipoleFuncSpl MultipoleFunc{Spline1D{Float64,Vector{Float64}}}
typealias MultipoleBasis Vector{MultipoleFuncSpl}



function eval_basis_sum(b::MultipoleBasis, g::AbstractGrid3d, sphharmtab::AbstractArray, disttab::AbstractArray)
    n,m = length(g),length(b)
    φ = zeros(n)

    i = 1; R = b[1].R
    ind = Int[(f.R != R ? begin R = f.R; i += 1 end : i) for f in b]  # map basis function index -> center index
    
    for j in 1:m
        f = b[j]
        k = iq(f.l,f.m)
        R = f.R
        for i in 1:n
            r = disttab[i,ind[j]]
            if r <= f.rmax
                φ[i] += (r > f.rcut ? f.q/r^(f.l+1) : f.u(r)) * sphharmtab[i,ind[j],k]
            end
        end
    end
    
    φ
end



