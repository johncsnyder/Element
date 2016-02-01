


function vcut(r::Float64, a::Float64, w::Float64, r0::Float64)
    rcut = r0 + w
    r < r0 ? 0. : a*exp(-w/(r - r0))/(r - rcut)^2
end

function vcut(r::AbstractVector{Float64}, a::Float64, w::Float64, r0::Float64)
    Float64[vcut(r[i],a,w,r0) for i in 1:length(r)]
end

function findr0{T}(u      ::AbstractVector{Float64},
                   r      ::AbstractGrid1d{T,1};
                   etacut ::Float64 = 1e-4,
                   r0max  ::Float64 = 12.)
    I = cumsum(r.w.*u.^2)
    r0 = r[searchsortedfirst(I, 1-etacut)-1]
    min(r0, r0max)
end




function naohydro!(basis  ::Basis,
                   Z      ::Float64,
                   n      ::Int,
                   l      ::Int,
                   r      ::LogarithmicGrid;
                   R      ::Point3f = (0.,0.,0.),
                   etacut ::Float64 = 1e-4,
                   a      ::Float64 = 10.,
                   w      ::Float64 = 2.,
                   r0max  ::Float64 = 20.,
                   npts   ::Int     = 100,
                   Emin   ::Float64 = -1.,
                   Emax   ::Float64 = 0.,
                   Elim   ::Float64 = 1e4)
    v = -Z./r                                           # bare Coulomb potential
    ev,ef = radschrod1(r,v,n,l,reverse=true,            # compute hydrogenic basis function
        Emin=Emin,Emax=Emax,Elim=Elim)
    r0 = findr0(ef,r,etacut=etacut,r0max=r0max)         # find start of cutoff potential
    rmax = r0 + w                                       # wavefunction vanishes outside rmax
    rconf = loggrid(log(r[1]), log(rmax), length(r))    # logarithmic grid for confining potential
    vconf = -Z./rconf + vcut(rconf,a,w,r0)              # confining potential
    vconf[end] = 0                                      # eliminate singularity at endpoint
    ev,ef = radschrod1(rconf,vconf,n,l,reverse=true,    # recompute with confining potential
        Emin=Emin,Emax=Emax,Elim=Elim)
    u = ef./rconf
    analytic_continuation!(rconf,u,r->r.^l)             # fix eigenfunction near the nucleus
    append!(basis, [NAO(rconf,u,l=l,m=m,R=R,rmax=rmax,npts=npts) for m in -l:l])
end




azimuthal_quantum_number_to_orbital_letter = Dict(
    0=>'s',
    1=>'p',
    2=>'d',
    3=>'f',
    4=>'g',
    5=>'h',
    6=>'i'
)

orbital_letter_to_azimuthal_quantum_number = Dict(
    's'=>0,
    'p'=>1,
    'd'=>2,
    'f'=>3,
    'g'=>4,
    'h'=>5,
    'i'=>6
)



# convert e.g. "1s 2s 2p" -> [(1,0),(2,0),(2,1)]
string_to_state(s::AbstractString) = (parse(Int,s[1]), orbital_letter_to_azimuthal_quantum_number[s[2]])
string_to_states(s::AbstractString) = Vector{Point2{Int}}(map(string_to_state, split(s)))

# convert e.g. [(1,0),(2,0),(2,1)] -> "1s 2s 2p"
state_to_string(state::Point2{Int}) = string(state[1], azimuthal_quantum_number_to_orbital_letter[state[2]])
states_to_string(states::Vector{Point2{Int}}) = join(map(state_to_string, states), " ")




function naominimal!(basis::Basis, at::Atom, states::AbstractString = "";
        etacut::Real = 1e-4, a::Real = 10., w::Real = 2., r0max::Real = 20., npts::Int = 100)
    naominimal!(basis,at,string_to_states(states),etacut=etacut,a=a,w=w,r0max=r0max,npts=npts)
end

function naominimal(at::Atom, states::AbstractString = "";
        etacut::Real = 1e-4, a::Real = 10., w::Real = 2., r0max::Real = 20., npts::Int = 100)
    naominimal(at,string_to_states(states),etacut=etacut,a=a,w=w,r0max=r0max,npts=npts)
end




function naominimal!(basis  ::Basis,
                     r      ::LogarithmicGrid,
                     vs     ::Spline1D,                             # spline radial KS potential
                     q      ::Matrix{Int},                          # orbital quantum numbers
                     f      ::Vector{Float64},                      # occupation numbers
                     efs    ::Matrix{Float64},                      # eigenfunctions
                     R      ::Point3f;                              # basis center
                     states ::Vector{Point2i} = Point2i[],          # states in q to select
                     etacut ::Float64         = 1e-4,
                     a      ::Float64         = 10.,
                     w      ::Float64         = 2.,
                     r0max  ::Float64         = 20.,
                     npts   ::Int             = 100,
                     Emin   ::Float64         = -1.,
                     Emax   ::Float64         = 0.,
                     Elim   ::Float64         = 1e4)
    all_states = isempty(states)  # take all states if true

    for i in 1:size(efs,2)
        n,l = q[i,:]

        if (all_states && f[i] > 1e-14) || (n,l) in states
            deleteat!(states, findin(states, [(n,l)]))          # remove (n,l) state
            ef = efs[:,i]
            r0 = findr0(ef,r,etacut=etacut,r0max=r0max)         # find start of cutoff potential
            rmax = r0 + w                                       # wavefunction vanishes outside rmax
            rconf = loggrid(log(r[1]), log(rmax), length(r))    # logarithmic grid for confining potential
            vconf = vs(rconf) + vcut(rconf,a,w,r0)              # confining potential
            vconf[end] = 0                                      # eliminate singularity at endpoint
            ev,ef = radschrod1(rconf,vconf,n,l,reverse=true,    # recompute with confining potential
                Emin=Emin,Emax=Emax,Elim=Elim)
            u = ef./rconf
            analytic_continuation!(rconf,u,r->r.^l)             # fix eigenfunction near the nucleus
            append!(basis, [NAO(rconf,u,l=l,m=m,R=R,rmax=rmax,npts=npts) for m in -l:l])
        end
    end

    @assert isempty(states) "could not find the specified states $states"
end



function naominimal!(basis  ::Basis,
                     at     ::Atom;
                     states ::Vector{Point2i} = Point2i[],
                     etacut ::Float64         = 1e-4,
                     a      ::Float64         = 10.,
                     w      ::Float64         = 2.,
                     r0max  ::Float64         = 20.,
                     npts   ::Int             = 100,
                     Emin   ::Float64         = -1.,
                     Emax   ::Float64         = 0.,
                     Elim   ::Float64         = 1e4)
    naominimal!(basis,at[:r],at[:vs].u,at[:q],at[:f],at[:efs],at.R,states=states,etacut=etacut,
        a=a,w=w,r0max=r0max,npts=npts,Emin=Emin,Emax=Emax,Elim=Elim)
end



function naominimal(at     ::Atom;
                    states ::Vector{Point2i} = Point2i[],
                    etacut ::Float64         = 1e-4,
                    a      ::Float64         = 10.,
                    w      ::Float64         = 2.,
                    r0max  ::Float64         = 20.,
                    npts   ::Int             = 100,
                    Emin   ::Float64         = -1.,
                    Emax   ::Float64         = 0.,
                    Elim   ::Float64         = 1e4)
    basis = NAOBasis()
    naominimal!(basis,at,states=states,etacut=etacut,a=a,w=w,
        r0max=r0max,npts=npts,Emin=Emin,Emax=Emax,Elim=Elim)
    basis
end




function naobasis!(at     ::Atom;
                   etacut ::Float64 = 1e-4,
                   a      ::Float64 = 10.,
                   w      ::Float64 = 2.,
                   r0max  ::Float64 = 20.,
                   npts   ::Int     = 100,
                   Emin   ::Float64 = -1.,
                   Emax   ::Float64 = 0.,
                   Elim   ::Float64 = 1e4)
    basis = naominimal(at,etacut=etacut,a=a,w=w,r0max=r0max,npts=npts)

    if at.Z in keys(nao_ionic)
        ion = atomicdft(at,f=nao_ionic_config[at.Z])
        states = nao_ionic[at.Z][1]
        naominimal!(basis,ion,states,etacut=etacut,a=a,w=w,r0max=r0max,
            npts=npts,Emin=Emin,Emax=Emax,Elim=Elim)
    end
    
    for (n,l,Z) in nao_basis[at.Z][1]  # first tier
        naohydro!(basis,Z,n,l,at[:r],R=at.R,etacut=etacut,a=a,w=w,
            r0max=r0max,npts=npts,Emin=Emin,Emax=Emax,Elim=Elim)
    end

    at[:basis] = basis
    return
end


# function naobasis(at::Atom;
#         etacut  ::Real = 1e-4, 
#         a       ::Real = 10., 
#         w       ::Real = 2., 
#         r0max   ::Real = 20., 
#         npts    ::Int  = 100)

#     at_ = copy(at)
#     naobasis!(at_;etacut=etacut,a=a,w=w,r0max=r0max,npts=npts)
#     at_
# end


