# API-INDEX


## MODULE: Element

[setposition(at::Element.Atom,  R::Tuple{Float64, Float64, Float64})](Element.md#method__setposition.1)  `setposition(at,R)`

[Element.Atom](Element.md#type__atom.1)  `Atom(Z; R=(0.,0.,0.), params...)`

[Element.Molecule](Element.md#type__molecule.1)  `Molecule(; atoms=[], params...)`

[Element.MultipoleFunc{T}](Element.md#type__multipolefunc.1)  Multipole function 

[Element.RadialFunc{T}](Element.md#type__radialfunc.1)  Spherically symmetric function 

[analytic_continuation!(r::AbstractArray{T, 1},  f::AbstractArray{T, 1},  g::Function)](Element.md#method__analytic_continuation.1)  `analytic_continuation!(r,f,g)`

[analytic_continuation(r,  f,  g)](Element.md#method__analytic_continuation.2)  `analytic_continuation(r,f,g)`

[analytic_continuation_coulomb_nuclei!(r::AbstractArray{Float64, 1},  f::AbstractArray{Float64, 2},  l::AbstractArray{Int64, 1})](Element.md#method__analytic_continuation_coulomb_nuclei.1)  `analytic_continuation_coulomb_nuclei!(r,f,l)`

[deepcopy(at::Element.Atom,  R::Tuple{Float64, Float64, Float64})](Element.md#method__deepcopy.1)  `deepcopy(at,R)`

[eig!(u::Array{Float64, 1},  v::Array{Array{Float64, 1}, 1},  aE::Function,  dx::Float64,  qmin::Integer,  qmax::Integer)](Element.md#method__eig.1)  `eig!(u,v,aE,dx,qmin,qmax;Emin=0.0,Emax=1.0,Elim=10000.)`

[eig(aE::Function,  dx::Float64,  qmin::Integer,  qmax::Integer)](Element.md#method__eig.2)  `eig(aE,dx,qmin,qmax;Emin=0.0,Emax=1.0,Elim=10000.)`

[enumeratelm(lmax)](Element.md#method__enumeratelm.1)  `enumeratelm(lmax)` lists all pairs `(l,m)` for `l=0,...,lmax` and `m=-l,...,l`.

[findsteps!(I::Array{Tuple{Float64, Float64}, 1},  f::Function,  xl::Float64,  xr::Float64,  fmin::Int64,  fmax::Int64)](Element.md#method__findsteps.1)  `findsteps!(I,f,xl,xr,[fmin,fmax];args=())` 

[findsteps(f::Function,  xl::Float64,  xr::Float64)](Element.md#method__findsteps.2)  `findsteps(f,xl,xr,[fmin,fmax];args=())` 

[fzero(f,  x1::Float64,  x2::Float64,  args...)](Element.md#method__fzero.1)  `fzero(f,x1,x2,args...;tol=1e-12,maxsteps=100)`

[genOh_a00(v)](Element.md#method__genoh_a00.1)   C version: Dmitri Laikov

[multipolemoment(r,  rho,  l)](Element.md#method__multipolemoment.1)  `multipolemoment(r,rho,l)`

[multistep_integrator!(y::AbstractArray{Float64, 1},  a::AbstractArray{Float64, 1},  h::Float64,  y1::Float64,  y2::Float64)](Element.md#method__multistep_integrator.1)  `multistep_integrator!(y,a,h,y1,y2)` is the same as `multistep_integrator`,

[multistep_integrator(a::AbstractArray{Float64, 1},  h::Float64,  y1::Float64,  y2::Float64)](Element.md#method__multistep_integrator.2)  `multistep_integrator(a,h,y1,y2)` solves 

[multistep_integrator_endpoint(y::AbstractArray{Float64, 1},  a::AbstractArray{Float64, 1},  h::Float64,  y1::Float64,  y2::Float64)](Element.md#method__multistep_integrator_endpoint.1)  `multistep_integrator_endpoint(y,a,h,y1,y2)` 

[multistep_integrator_node_count(y::AbstractArray{Float64, 1},  a::AbstractArray{Float64, 1},  h::Float64,  y1::Float64,  y2::Float64)](Element.md#method__multistep_integrator_node_count.1)  `multistep_integrator_node_count(y,a,h,y1,y2)` 

[multistep_rule(a_jm1::Float64,  a_j::Float64,  a_jp1::Float64,  y_jm1::Float64,  y_j::Float64,  h::Float64)](Element.md#method__multistep_rule.1)  `multistep_rule(a_jm1,a_j,a_jp1,y_jm1,y_j,h)`

[radpoisson(r::Element.LogarithmicGrid{T, A<:AbstractArray{T, 1}},  rho::AbstractArray{T, 1},  l::Int64)](Element.md#method__radpoisson.1)  `radpoisson(r,rho,l)`

[radschrod(r::Element.LogarithmicGrid{Float64, A<:AbstractArray{T, 1}},  v::AbstractArray{Float64, 1},  nmax::Int64,  lmax::Int64)](Element.md#method__radschrod.1)  `radschrod(r,v,nmax,lmax;reverse=false,Emin=-1.0,Emax=0.0,Elim=1e4)`

[radschrod1(r::Element.UniformGrid{Float64, A<:AbstractArray{T, 1}},  v::AbstractArray{Float64, 1},  n::Int64,  l::Int64)](Element.md#method__radschrod1.1)  `radschrod1(r,v,n,l;reverse=false,Emin=-1.0,Emax=0.0,Elim=1e4)`

[realsphharm(l::Int64,  m::Int64,  θ::Real,  ϕ::Real)](Element.md#method__realsphharm.1)  `realsphharm(l,m,θ,ϕ)`

[sphharm(l::Int64,  m::Int64,  θ::Real,  ϕ::Real)](Element.md#method__sphharm.1)  `sphharm(l,m,θ,ϕ)`

[splineradial!(at::Element.Atom)](Element.md#method__splineradial.1)  `splineradial!(at; npts=100)`

[Element.CoulombPotential](Element.md#type__coulombpotential.1)  `CoulombPotential(Z;R=(0.,0.,0.))`

[Element.ldaxc](Element.md#type__ldaxc.1)  `ldaxc` is the pw-lda exchange-correlation functional [1,2].

