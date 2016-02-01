


doc"""
`multistep_rule(a_jm1,a_j,a_jp1,y_jm1,y_j,h)`

computes the multistep rule in Numerov's method:

$y_{j+1} = { \left(2-5h^2 a_j/6\right)y_j - \left(1+h^2 a_{j-1}/12\right)y_{j-1} 
    \over 1+h^2 a_{j+1}/12 }$
"""
@fastmath function multistep_rule(a_jm1::Float64, a_j::Float64, a_jp1::Float64, 
                                    y_jm1::Float64, y_j::Float64, h::Float64)
    ((2 - 5*h*h*a_j/6)*y_j - (1 + h*h*a_jm1/12)*y_jm1) / (1 + h*h*a_jp1/12)
end


doc"""
`multistep_integrator!(y,a,h,y1,y2)` is the same as `multistep_integrator`,
but overwrites `y`.
"""
function multistep_integrator!(y::AbstractVector{Float64}, a::AbstractVector{Float64},
                                    h::Float64, y1::Float64, y2::Float64)
    y[1] = y1
    y[2] = y2
    @inbounds for j in 2:length(a)-1
        y[j+1] = multistep_rule(a[j-1], a[j], a[j+1], y[j-1], y[j], h)
    end
end


doc"""
`multistep_integrator(a,h,y1,y2)` solves 

$y''(x) + a(x) y(x) = 0$ 

on a uniform grid via the 4th order 
[Numerov's method](https://en.wikipedia.org/wiki/Numerov%27s_method).

Assuming a uniform grid $[x_1, \dots, x_n]$ with spacing `h`

`a` - discretization of $a(x)$, $[a(x_1), \dots, a(x_n)]$

`y1,y2` - initial boundary condition, specifying $y(x_1)$ and $y(x_2)$

Returns the solution $[y(x_1), \dots, y(x_n)]$.
"""
function multistep_integrator(a::AbstractVector{Float64}, h::Float64, y1::Float64, y2::Float64)
    y = similar(a)
    multistep_integrator!(y,a,h,y1,y2)
    y
end


doc"""
`multistep_integrator_endpoint(y,a,h,y1,y2)` 

is the same as `multistep_integrator`,
but returns only $y(x_n)$. 

`y` is used as temporary storage.
"""
function multistep_integrator_endpoint(y::AbstractVector{Float64}, a::AbstractVector{Float64},
                                    h::Float64, y1::Float64, y2::Float64)
    multistep_integrator!(y,a,h,y1,y2)
    y[end]
end


doc"""
`multistep_integrator_node_count(y,a,h,y1,y2)` 

is the same as `multistep_integrator`,
but returns the number of nodes in $y(x)$ (i.e. the number of sign changes). 

`y` is used as temporary storage.
"""
function multistep_integrator_node_count(y::AbstractVector{Float64}, a::AbstractVector{Float64},
                                    h::Float64, y1::Float64, y2::Float64)
    y[1] = y1
    y[2] = y2
    p = 0
    @inbounds for j in 2:length(a)-1
        y[j+1] = multistep_rule(a[j-1], a[j], a[j+1], y[j-1], y[j], h)
        if sign(y[j+1]) != sign(y[j])
            p += 1
        end
    end
    p
end





function expandrange(q::Function, qmin::Int, qmax::Int, Emin::Float64, Emax::Float64, Elim::Float64)
    while q(Emin)::Int > qmin
        Emin -= Emax - Emin
        # @assert q(Emin) <= q(Emax) "node count is not monotonically increasing"
        @assert abs(Emin) < Elim "Emin exceeded Elim. decrease grid spacing, or increase Elim"
    end

    while q(Emax)::Int <= qmax
        Emax += Emax - Emin
        @assert abs(Emax) < Elim "Emax exceeded Elim. decrease grid spacing, or increase Elim"
    end

    Emin,Emax
end



doc"""
`eig(aE,dx,qmin,qmax;Emin=0.0,Emax=1.0,Elim=10000.)`

solves the eigenvalue equation 

$y''(x) + a_E(x) y(x) = 0$

via the shooting method and Numerov's multistep integrator (see `multistep_integrator`),
on a uniform grid $[x_1, \dots, x_n]$ with spacing `dx`
and boundary conditions $y(x_1) = 0$ and $y(x_n) = 0$,

`aE` - the function $a_E(x)$, which returns a vector $[a_E(x_1), \dots, a_E(x_n)]$
for a given value of `E`.

`qmin,qmax` - solve for eigenfunctions with `qmin < q < qmax` nodes.

`Emin,Emax` - initial `E` interval to search. this will be expanded automatically
 (up to `±Elim`) until the range includes all requested eigenvalues.

returns `(u,v,q)`

where `u` are the eigenvalues, `v` are the eigenfunctions (the jth
eigenfunction given by `v[:,j]`), and `q` are the respective number of nodes in
the eigenfunctions.

e.g. *particle in a box*
```julia
x = lingrid(0,1,100)
dx = x[2]-x[1]
v = zeros(x)
evs,efs,q = eig(E -> 2*(E-v),dx,0,4)
```
"""
function eig(aE::Function, dx::Float64, qmin::Integer, qmax::Integer; 
        Emin=0., Emax=1., Elim=10000.)
    y = similar(aE(Emin))::Vector{Float64}
    
    integrator = E -> multistep_integrator(aE(E), dx, 0., 1.)
    F = E -> multistep_integrator_endpoint(y, aE(E), dx, 0., 1.)  # matching function
    q = E -> multistep_integrator_node_count(y, aE(E), dx, 0., 1.)  # node count
    
    Emin,Emax = expandrange(q,qmin,qmax,Emin,Emax,Elim)
    
    I = findsteps(q, Emin, Emax, qmin, qmax)  # intervals
    u = Float64[fzero(F, a, b) for (a,b) in I]  # compute eigenvalues
    v = hcat([integrator(E) for E in u]...)  # compute eigenvectors
    q = Int64[q(a) for (a,b) in I]  # node counts
    
    u, v, q
end



doc"""
`eig!(u,v,aE,dx,qmin,qmax;Emin=0.0,Emax=1.0,Elim=10000.)`

same as `eig` except that the eigenvalues and eigenfunctions are appended to `u` 
and `v`, respectively.

returns `q`, the number of nodes in the eigenfunctions.
"""
function eig!(u::Vector{Float64}, v::Vector{Vector{Float64}},
        aE::Function, dx::Float64, qmin::Integer, qmax::Integer; 
        Emin=0., Emax=1., Elim=10000.)
    y = similar(aE(Emin))::Vector{Float64}
    
    integrator = E -> multistep_integrator(aE(E), dx, 0., 1.)
    F = E -> multistep_integrator_endpoint(y, aE(E), dx, 0., 1.)  # matching function
    q = E -> multistep_integrator_node_count(y, aE(E), dx, 0., 1.)  # node count
    
    Emin,Emax = expandrange(q,qmin,qmax,Emin,Emax,Elim)
    
    I = findsteps(q, Emin, Emax, qmin, qmax)  # intervals
    u_ = Float64[fzero(F, a, b) for (a,b) in I]  # compute eigenvalues
    q = Int64[q(a) for (a,b) in I]  # node counts
    
    append!(u, u_)
    append!(v, [integrator(E) for E in u_])  # compute eigenvectors

    q
end
