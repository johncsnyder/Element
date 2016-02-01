


aE(E::Float64,l::Int,r::Float64,v::Float64) = 2*(E - v - l*(l+1)/(2r*r))

aE(E::Float64,l::Int,r::AbstractVector{Float64},v::AbstractVector{Float64}) =
    Float64[aE(E,l,r[i],v[i]) for i in 1:length(r)]


aElog(E::Float64,l::Int,r::Float64,v::Float64) = 2r*r*(E-v) - (l+0.5)^2

aElog(E::Float64,l::Int,r::AbstractVector{Float64},v::AbstractVector{Float64}) =
    Float64[aElog(E,l,r[i],v[i]) for i in 1:length(r)]




doc"""
`radschrod1(r,v,n,l;reverse=false,Emin=-1.0,Emax=0.0,Elim=1e4)`

solves the radial Schrödinger equation

$\left\{ -{1\over 2}{d^2\over dr^2} + {l(l+1)\over 2r^2} + v(r) \right\} u(r) = \epsilon u(r)$

where $\varphi(r,\Omega) = {u(r)\over r} Y_l^m(\Omega)$.

for a single eigenvalue/eigenfunction, specified by `n,l`. 

the eigenfunction is normalized such that $\int dr\, u(r)^2 = 1$

`r` - uniform grid (`lingrid`, `simpsgrid`) or logarithmic grid (`loggrid`)

`v` - the radial potential $v(r)$ evaluated on the grid `r`

`reverse` - if true, integrate in reverse. 

`Emin,Emax` - initial `E` interval to search. this will be expanded automatically
 (up to `±Elim`) until the range includes all requested eigenvalues.

returns `ev,ef`, the requested eigenvalue $\epsilon$ and the eigenfunction $u(r)$

e.g. *hydrogen atom 1s orbital*
```julia
r = loggrid(-12,4,500)
v = -1./r
ev,ef = radschrod1(r,v,1,0,reverse=true)
```
"""
function radschrod1(r::UniformGrid{Float64}, v::AbstractVector{Float64},
        n::Int, l::Int; reverse=false, Emin=-1., Emax=0., Elim=1e4)
    r_ = reverse ? Base.reverse(r) : r
    v_ = reverse ? Base.reverse(v) : v

    evs, efs, q = eig(E -> aE(E,l,r_,v_), r.dx, n-l-1, n-l-1, 
        Emin=Emin, Emax=Emax, Elim=Elim)

    efs = reverse ? Base.reverse(efs[:,1]) : efs[:,1]  # recover wavefunction
    efs ./= sqrt(integrate(r, efs.^2))  # normalization

    evs[], efs
end


function radschrod1(r::LogarithmicGrid{Float64}, v::AbstractVector{Float64},
        n::Int, l::Int; reverse=false, Emin=-1., Emax=0., Elim=1e4)
    r_ = reverse ? Base.reverse(r) : r
    v_ = reverse ? Base.reverse(v) : v

    evs, efs, q = eig(E -> aElog(E,l,r_,v_), r.dx, n-l-1, n-l-1, Emin=Emin, Emax=Emax, Elim=Elim)

    efs = sqrt(r) .* (reverse ? Base.reverse(efs[:,1]) : efs[:,1])  # recover wavefunction
    efs ./= sqrt(integrate(r, efs.^2))  # normalization

    evs[], efs
end




doc"""
`radschrod(r,v,nmax,lmax;reverse=false,Emin=-1.0,Emax=0.0,Elim=1e4)`

solves the radial Schrödinger equation

$\left\{ -{1\over 2}{d^2\over dr^2} + {l(l+1)\over 2r^2} + v(r) \right\} u(r) = \epsilon u(r)$

where $\varphi(r,\Omega) = {u(r)\over r} Y_l^m(\Omega)$.

for all eigenvalues/eigenfunctions up to `nmax,lmax` (`nmax > 0`, `lmax >= 0` and `lmax < nmax`).

the eigenfunctions are normalized such that $\int dr\, u(r)^2 = 1$

`r` - logarithmic grid (see `loggrid`)

`v` - the radial potential $v(r)$ evaluated on the grid `r`

`reverse` - if true, integrate in reverse. 

`Emin,Emax` - initial `E` interval to search. this will be expanded automatically
 (up to `±Elim`) until the range includes all requested eigenvalues.

returns `evs,efs,q`, the requested eigenvalues, eigenfunctions and a list of 
the corresponding quantum numbers `n,l`

e.g. *hydrogenic atom, Z=1*
```julia
r = loggrid(-12,5,1000)
v = -1./r
nmax,lmax = 3,2
evs,efs,q = radschrod(r,v,nmax,lmax,reverse=true)
```
"""
function radschrod(r::LogarithmicGrid{Float64}, v::AbstractVector{Float64},
        nmax::Int, lmax::Int; reverse=false, Emin=-1., Emax=0., Elim=1e4)
    r_ = reverse ? Base.reverse(r) : r
    v_ = reverse ? Base.reverse(v) : v
    
    evs = Float64[]  # eigenvalues
    efs = Vector{Float64}[]  # eigenfunctions
    q = Array{Int}[]  # (n,l) quantum numbers

    for l in 0:lmax
        q_ = eig!(evs, efs, E -> aElog(E,l,r_,v_), r.dx, 0, nmax-l-1, 
            Emin=Emin, Emax=Emax, Elim=Elim)
        append!(q, [[k+l+1 l] for k in q_])
    end
    
    ind = sortperm(evs)
    evs = evs[ind]
    efs = efs[ind]
    q = vcat(q[ind]...)

    for i in 1:length(efs)
        l = q[i,2]
        efs[i] = sqrt(r) .* (reverse ? Base.reverse(efs[i]) : efs[i])  # recover wavefunction
        efs[i] ./= sqrt(integrate(r, efs[i].^2))  # normalization
    end

    efs = hcat(efs...)
    
    evs, efs, q
end




doc"""
`analytic_continuation!(r,f,g)`

same as `analytic_continuation` but mutates `f`.
"""
function analytic_continuation!(r::AbstractVector, f::AbstractVector, g::Function)
    q = f./g(r)
    i = findfirst(diff(abs(diff(q))) .> 0.0, true)
    c = q[i]
    f[1:i] = c*g(r[1:i])
    return
end


doc"""
`analytic_continuation(r,f,g)`

smoothly replaces $f(r)$ with the analytic form given by $g(r)$ as $r\to 0$

the transition point is determined as the point where 
the derivatives of `f` and `g` match.

analytic continuation as $r\to\infty$ can be acheived by reversing
`r` and `f`.

e.g.
```julia
r = loggrid(-12,4.,1000)
n,l = 3,1
ev,ef = radschrod1(r,-1./r,n,l,reverse=true)
efc = analytic_continuation(r,ef,r -> r.^(l+1))
```
"""
function analytic_continuation(r,f,g)
    fc = copy(f)
    analytic_continuation!(r,fc,g)
    return fc
end



doc"""
`analytic_continuation_coulomb_nuclei!(r,f,l)`

analytically continues multiple functions $f(r)$  (columns of `f`) with
the analytic form $r^{l+1}$ as $r\to 0$.

`f` - a matrix whose columns are functions of `r`

`l` - corresponding `l`-quantum numbers
"""
function analytic_continuation_coulomb_nuclei!(r::AbstractVector{Float64}, 
        f::AbstractMatrix{Float64}, l::AbstractVector{Int}) 
    for i in 1:size(f,2)
        analytic_continuation!(r, slice(f,:,i), r -> r.^(l[i]+1))
    end
end

function analytic_continuation_coulomb_nuclei!(r::AbstractVector{Float64}, 
        f::AbstractVector{Float64}, l::Int) 
    analytic_continuation!(r, f, r -> r.^(l+1))
end


# import GSL.sf_laguerre_n

# function R_ha(Z::Real,n::Integer,l::Integer,r::Real)
#     x = 2Z*r/n
#     N = factorial(2l+1+n-l-1)*sqrt((2Z/n)^3*factorial(n-l-1)/(2n)/factorial(n+l)^3)
#     N*x^l*exp(-x/2)*sf_laguerre_n(n-l-1,2l+1,x)
# end

# function R_ha(Z::Real,n::Integer,l::Integer,r::AbstractArray)
#     reshape([R_ha(Z,n,l,r_) for r_ in r], size(r))
# end

# u_ha(Z::Real,n::Integer,l::Integer,r::AbstractArray) = r.*R_ha(Z,n,l,r)


