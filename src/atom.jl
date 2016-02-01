


doc"""
`Atom(Z; R=(0.,0.,0.), params...)`

creates an atom at `R` with nuclear charge `Z`,
with any additional parameters `params`.
"""
type Atom <: Associative{Symbol,Any}
    Z::Float64              # atomic number
    R::Point3f              # position
    d::Dict{Symbol,Any}     # additional properties
end

Base.delete!(at::Atom, key)                                    = delete!(at.d, key)
Base.done(at::Atom, i)                                         = done(at.d, i)
Base.empty!(at::Atom)                                          = empty!(at.d)
Base.copy(at::Atom)                                            = Atom(at.Z,at.R,copy(at.d))
Base.get(at::Atom, key, default)                               = get(at.d, key, default)
Base.get(default::Union{DataType,Function}, at::Atom, key)     = get(default, at.d, key)
Base.get!(at::Atom, key0, default)                             = get!(at.d, key0, default)
Base.get!(default::Union{DataType,Function}, at::Atom, key0)   = get!(default, at.d, key0)
Base.getindex(at::Atom, key)                                   = getindex(at.d, key)
Base.getkey(at::Atom, key, default)                            = getkey(at.d, key, default)
Base.haskey(at::Atom, key)                                     = haskey(at.d, key)
Base.isempty(at::Atom)                                         = isempty(at.d)
Base.length(at::Atom)                                          = length(at.d)
Base.next(at::Atom, i)                                         = next(at.d, i)
Base.pop!(at::Atom, key)                                       = pop!(at.d, key)
Base.pop!(at::Atom, key, default)                              = pop!(at.d, key, default)
Base.setindex!(at::Atom, v0, key0)                             = setindex!(at.d, v0, key0)
Base.similar(at::Atom)                                         = similar(at.d)
Base.sizehint!(at::Atom, newsz)                                = sizehint!(at.d, newsz)
Base.start(at::Atom)                                           = start(at.d)
Base.merge!(at::Atom, d::Dict)                                 = begin merge!(at.d, d); d end
Base.merge(at::Atom, d::Dict)                                  = Atom(merge(at.d, d))
Base.writemime(io::IO, mime::MIME"text/plain", at::Atom)       = show_symbol_dict(io,at)
Base.summary(at::Atom)                                         = string(typeof(at),"(Z=$(at.Z),R=$(at.R))")
Base.show(io::IO, at::Atom)                                    = print(io,summary(at))



function show_symbol_dict(io::IO, d::Associative{Symbol,Any}; cols=100)
    shown_set = get(task_local_storage(), :SHOWNSET, nothing)
    if shown_set === nothing
        shown_set = ObjectIdDict()
        task_local_storage(:SHOWNSET, shown_set)
    end
    d in keys(shown_set) && (print(io, "#= circular reference =#"); return)

    try
        println(io, summary(d))

        shown_set[d] = true
        ks = map(string,keys(d))
        keylen = length(ks) > 0 ? maximum(map(length,ks)) : 0

        for (i,k) in enumerate(keys(d))
            print(io, Base.rpad(ks[i], keylen))
            print(io, " = ")
            val = Base.with_output_limit(()->sprint(show, d[k]))
            val = Base._truncate_at_width_or_chars(val, cols - keylen, "\r\n")
            print(io, val)
            print(io, '\n')
        end
    finally
        delete!(shown_set, d)
    end
end

function Base.call(::Type{Atom};
        Z         ::Real         = 1.,
        R         ::Point3{Real} = (0.,0.,0.),
        makebasis ::Bool         = false,
        makegrid  ::Bool         = false,
        etacut    ::Float64      = 1e-4,
        a         ::Float64      = 10.,
        w         ::Float64      = 2.,
        r0max     ::Float64      = 20.,
        npts      ::Int          = 100,
        nr        ::Int          = 80,
        na        ::Int          = 38,
        params...)
    
    at = Atom(convert(Float64,Z), convert(Point3f,R), Dict{Symbol,Any}(params))

    if makebasis
        dft!(at)
        splineradial!(at;npts=npts)
        naobasis!(at;etacut=etacut,a=a,w=w,r0max=r0max,npts=npts)
    end

    if makegrid
        makegrid!(at;nr=nr,na=na)
    end

    at
end

typealias Atoms Vector{Atom}


doc"""
`setposition(at,R)`

sets the position of atom `at` to `R`
"""
function setposition(at::Atom, R::Point3f)
    at.R = R
    for k in keys(at)
        setposition(at[k],R)
    end
end

function setposition(a::AbstractArray, R::Point3f)
    for i in 1:length(a)
        setposition(a[i],R)
    end
end

setposition(a, R::Point3f) = return



doc"""
`deepcopy(at,R)`

returns a copy of `at` at the position `R`
"""
function Base.deepcopy(at::Atom, R::Point3f)
    at_ = deepcopy(at)
    setposition(at_,R)
    at_
end






doc"""
`splineradial!(at; npts=100)`

splines the electronic density `rho`, hartree potential `vh`,
and Kohn-Sham potential `vs` and creates radial functions.
Additionally, the Coulomb potential `v` is created.

assumes the `rho`, `vh` and `vs` have been calculated and are
given in `at` discretized on the radial grid `r`.

*sets*

`at[:v]` (see `CoulombPotential`)

`at[:rho]`, `at[:vh]`, `at[:vs]` (see `RadialFunc`).

`npts` - the number of evenly distributed spline points to use

e.g.
```julia
at = dft(Z=1.)
splineradial!(at)
plot(at[:rho].u,0,5,label="rho")
plot!(at[:vh].u,0,5,label="vh")
plot!(at[:vs].u,0,5,label="vs")
ylims!(-2,1)
```
"""
function splineradial!(at::Atom; npts::Int=100)
    r = at[:r]
    Y00 = 1/sqrt(4π)
    q = 4π*multipolemoment(r,at[:rho],0)/Y00
    at[:v]   = CoulombPotential(at.Z,R=at.R)
    at[:rho] = RadialFunc(r,at[:rho],R=at.R,rmax=r[end],npts=npts)
    at[:vh]  = MultipoleFunc(r,at[:vh]./Y00,l=0,m=0,R=at.R,q=q,rcut=r[end],npts=npts)
    at[:vs]  = RadialFunc(r,at[:vs],R=at.R,rmax=r[end],npts=npts)
    return
end

function splineradial(at::Atom; npts::Int=100)
    let at = copy(at)
        splineradial!(at;npts=npts)
        return at
    end
end


function makegrid!(at::Atom; nr::Int=80, na::Int=38, rmax::Float64=Inf)
    if isinf(rmax)
        @assert :basis in keys(at) "create an atomic basis before making a grid"
        rmax = maximum([φ.rmax for φ in at[:basis]])
    end
    r = atomicradgrid(nr,rmax)
    Ω = lebgrid(na)
    at[:g] = atomicgrid(r,Ω,R=at.R)
    return
end

function makegrid(at::Atom; nr::Int=80, na::Int=38, rmax::Float64=Inf)
    let at = copy(at)
        makegrid!(at;nr=nr,na=na,rmax=rmax)
        return at
    end
end


function makegrid!(atoms::Atoms; nr::Int=80, na::Int=38, rmax::Float64=Inf)
    for at in atoms
        makegrid!(at,nr=nr,na=na,rmax=rmax)
    end
end

function makegrid(atoms::Atoms; nr::Int=80, na::Int=38, rmax::Float64=Inf)
    let atoms = copy(atoms)
        makegrid!(atoms;nr=nr,na=na,rmax=rmax)
        return atoms
    end
end


hasbasis(at::Atom) = :basis in keys(at)
hasgrid(at::Atom) = :g in keys(at)

print_stats(at::Atom) = print_stats(at[:stats])




