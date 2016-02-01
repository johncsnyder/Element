



# evs = eigenvalues
# k = orbital degeneracies
# N = num of particles
# T = temperature
function compute_occupations(evs::Vector{Float64}, k::Vector{Int}, N::Float64, T::Float64)
    @assert sum(k) >= N "$(sum(k)) states < $N electrons. increase nmax,lmax."
    mu = length(evs) == 1 ? evs[1] :
        fzero(mu -> N - dot(k,fermidirac(evs,mu,T)), evs[1], evs[end]+0.1)  # chemical potential
    f = k.*fermidirac(evs,mu,T)  # occupation numbers
    f,mu
end

function compute_occupations(evs::Vector{Float64}, k::Int, N::Float64, T::Float64)
    mu = length(evs) == 1 ? evs[1] :
        fzero(mu -> N - sum(k*fermidirac(evs,mu,T)), evs[1], evs[end]+0.1)  # chemical potential
    f = k*fermidirac(evs,mu,T)  # occupation numbers
    f,mu
end



init_mixer(mixer::AbstractMixer, grid::AbstractGrid) = return

function init_mixer(mixer::PulayMixer, grid::AbstractGrid)
    mixer.grid = grid
end

function init_scf(scf::SCF, grid::AbstractGrid)
    scf.grid = grid
    scf.converged = false
end


using Base: Callable

function dft(Z         ::Float64,
             N         ::Float64,
             f         ::Vector{Float64},
             fixed_occ ::Bool,
             nmax      ::Int,
             lmax      ::Int,
             r         ::LogarithmicGrid,
             xc        ::ldaxc,
             scf       ::SCF,
             mixer     ::AbstractMixer,
             T         ::Float64,
             Emin      ::Float64,
             Emax      ::Float64,
             Elim      ::Float64,
             logger    ::Callable)
    local E,Etot,vs,vh,evs,efs,q,mu

    stats = Stats()
    init_mixer(mixer,r)
    init_scf(scf,r)

    @stats "* total" begin

        # initialization
        if fixed_occ; N = sum(f); mu = NaN end
        rho_in = isempty(scf.rho) ? N^4*exp(-N*r/2)/64pi : scf.rho[end]         # initial density
        vext = -Z./r                                                            # external potential
        dr = 4pi*r.^2                                                           # radial integration factor

        # storage
        vx = similar(r); vc = similar(r); ex = similar(r); ec = similar(r)

        while !scf.converged
            
            @stats "hartree"                    vh = 4π*radpoisson(r,rho_in,0)                     # Hartree potential
            @stats "exchange-correlation"       xc(rho_in,ex,vx,ec,vc)                              # exchange-correlation
                                                vs = vext + vx + vc + vh                            # Kohn-Sham potential
            @stats "KS equations"               evs,efs,q = radschrod(r,vs,nmax,lmax,reverse=true,
                                                    Emin=Emin,Emax=Emax,Elim=Elim)
                                                l = q[:,2]                                          # l-quantum numbers
            @stats "analytic continuation"      analytic_continuation_coulomb_nuclei!(r,efs,l)

            if !fixed_occ; f,mu = compute_occupations(evs,2(2l+1),N,T) end
            @assert length(evs) >= length(f) "$(length(evs)) states < $(length(f)) required. increase nmax,lmax."

            @stats "compute density"            rho_out = (efs.^2 * f)./dr                                      # radial density
                                                E = dot(f,evs)                                                  # sum of eigenvalues
            @stats "harris energy"              Etot = E + integrate(r, dr.*rho_in.*(ex+ec - (vx+vc) - vh/2))   # Harris functional
            @stats "check convergence"          scf(Etot, rho_out)                                              # check convergence
            @stats "density mixing"             rho_in = mixer(rho_in, rho_out)                                 # mix densities

            if !is(logger,Void) logger(scf) end

        end  # scf cycle

        Ex   = integrate(r, rho_in.*ex)         # exchange energy
        Ec   = integrate(r, rho_in.*ec)         # correlation energy
        Exc  = Ex + Ec                          # exchange-correlation energy
        Eh   = integrate(r, rho_in.*vh/2)       # Hartree energy
        Epot = integrate(r, rho_in.*vext)       # potential energy (external)
        Ekin = Etot - Epot - Eh - Ex - Ec       # kinetic energy

    end  # total

    struct(Etot=Etot,Ekin=Ekin,Epot=Epot,Eh=Eh,Ex=Ex,Ec=Ec,Exc=Exc,rho=rho_in,
        vs=vs,vh=vh,evs=evs,efs=efs,q=q,N=N,f=f,mu=mu,stats=stats)
end




function dft!(at::Atom; fixed_occ=false, Z=at.Z, R=at.R, params...)
    at.Z = Z = convert(Float64,Z)
    at.R = R = convert(Point3f,R)
    merge!(at, Dict{Symbol,Any}(params))

    res = dft(
        Z,                                                      # atomic number
        get!(at, :N, Z),                                        # particle number
        fixed_occ ? at[:f] : Float64[],                         # occupations
        fixed_occ,                                              # fixed occupations
        get!(at, :nmax,   nlmax[Z][1]),                         # max n-quantum number
        get!(at, :lmax,   nlmax[Z][2]),                         # max l-quantum number
        get!(at, :r,      loggrid(atomicloggrid[Z]...)),        # radial grid
        get!(at, :xc,     ldaxcunpol()),                        # xc functional
        get!(at, :scf,    SCF()),
        get!(at, :mixer,  PulayMixer()),                        # density mixer
        fixed_occ ? NaN : get!(at, :T, 1e-6),                   # temperature
        get!(at, :Emin,   -1.),                                 # see radschrod
        get!(at, :Emax,   0.),                                  # see radschrod
        get!(at, :Elim,   1e4),                                 # see radschrod
        get!(at, :logger, Void)
    )
    
    merge!(at, res)
    return
end


function dft(at::Atom=Atom(); fixed_occ=false, Z=at.Z, R=at.R, params...)
    let at = deepcopy(at)
        dft!(at,fixed_occ=fixed_occ,Z=Z,R=R,params...)
        return at
    end
end

function dft(; fixed_occ=false, Z=1.0, R=(0.,0.,0.), params...)
    at = Atom(Z=Z)
    let at = deepcopy(at)
        dft!(at,fixed_occ=fixed_occ,Z=Z,R=R,params...)
        return at
    end
end



for mixer in (LinearMixer,PulayMixer)
    eval(quote
        precompile(dft, (Float64,Float64,Vector{Float64},Bool,Int,Int,LogarithmicGrid,ldaxc,SCF,$mixer,Float64,Float64,Float64,Float64,Callable))
    end)
end

precompile(dft!, (Atom,))
precompile(dft, (Atom,))




