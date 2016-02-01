


function dft(atoms  ::Atoms,                # atoms
             basis  ::Basis,                # molecular basis
             g      ::MultiGrid,            # integration grid
             N      ::Float64,              # number of particles
             vembed ::AbstractArray,        # embedding potential 
             xc     ::ldaxc,                # exchange-correlation functional
             scf    ::SCF,
             mixer  ::AbstractMixer,        # density mixing
             lmax   ::Int,                  # max l-quantum number in hartree multipole expansion
             nh     ::Int,                  # number of hartree radial points
             nspl   ::Int,                  # number of spline points
             T      ::Float64,              # temperature
             logger ::Callable)
    local E,Etot,Enuc,vs,vh,δrho_mp,δvh_mp,evs,efs,f,mu,H,S,rho_out

    stats = Stats()
    init_mixer(mixer,g)
    init_scf(scf,g)

    @stats "* total" begin

        # initialize
        @stats "* evaluate vext"                vext     = sum([at[:v](g) for at in atoms])                 # external potential
                                                vext    += vembed                                           # embedding potential
        @stats "* evaluate rho_free"            rho_free = sum([at[:rho](g) for at in atoms])               # free atom density
        @stats "* evaluate vh_free"             vh_free  = sum([at[:vh](g) for at in atoms])                # free atom hartree potential
                                                rho_in   = isempty(scf.rho) ? rho_free : scf.rho[end]       # initial guess for density

                                                sort!(basis, by=φ->φ.R)                                     # sort basis by center
        @stats "* evaluate basis"               φ,∆φ = eval_basis_lapl(basis,g)                             # evaluate basis on grid
        @stats "* overlap"                      wφT  = (g.w.*φ)'                                            # grid weights times basis
        @stats "* overlap"                      S    = Symmetric(wφT*φ)                                     # overlap matrix
                                                Enuc = sum([μ != ν ? μ.Z*ν.Z/dist(μ.R,ν.R)/2 : 0.           # nuclear attraction
                                                    for μ in atoms, ν in atoms])

                                                lmax = max(lmax,maximum(Int[f.l for f in basis]))
        @stats "* sphharm basis"                realsphharmbasis!(atoms,lmax)                               # atomic spherical harmonics basis
        @stats "* sphharm tables"               sphharmtab,disttab = compute_sphharm_dist_tab(basis,g,lmax) # spherical harmonics, distances

        # storage
        vx = similar(rho_in); vc = similar(rho_in)
        ex = similar(rho_in); ec = similar(rho_in)
        vh = similar(rho_in); vs = similar(rho_in)
        hφ = similar(φ)

        while !scf.converged
                                                δrho = rho_in - rho_free
            @stats "poisson multipole"          δrho_mp,δvh_mp = poissonmultipole(g,δrho,atoms,         # solve hartree potential
                                                    nh=nh,nspl=nspl)
            @stats "eval hartree multipole"     vh[:] = isempty(δvh_mp) ? vh_free : vh_free + eval_basis_sum(δvh_mp,g,sphharmtab,disttab)

            @stats "exchange-correlation"       xc(rho_in,ex,vx,ec,vc)                                  # exchange-correlation
                                                vs[:] = vext + vx + vc + vh                             # Kohn-Sham potential

            @stats "hamiltonian"                h!(hφ,φ,∆φ,vs)                                          # apply Hamiltonian to basis
            @stats "hamiltonian"                H = Symmetric(wφT*hφ)                                   # Hamiltonian matrix
            @stats "KS equations"               evs,efs = LinAlg.eig(H,S)                               # compute eigenspectrum

                                                f,mu = compute_occupations(evs,2,N,T)
                                                ζ = f .> 1e-12                                          # nonzero occupations
                                                efs = efs[:,ζ]; evs = evs[ζ]; f = f[ζ]                  # discard unused eigenspectrum

            # @stats "normalize efs"              efs ./= sqrt(integrate(g,abs2(φ*efs)))'                 # normalize eigenfunctions
            
            @stats "compute density"            rho_out = abs2(φ*efs)*f                                 # compute density
                                                E = dot(f,evs)                                          # sum of eigenvalues
            @stats "Harris functional"          Etot = E + integrate(g, rho_in.*(ex+ec - (vx+vc) - vh/2)) + Enuc    # Harris functional

            @stats "check convergence"          scf(Etot,rho_out)                                       # check convergence
            @stats "density mixing"             rho_in = mixer(rho_in,rho_out)                          # mix densities

            if !is(logger,Void) logger(scf) end

        end  # scf cycle
        
        Ex   = integrate(g, rho_in.*ex)                         # exchange energy
        Ec   = integrate(g, rho_in.*ec)                         # correlation energy
        Exc  = Ex + Ec                                          # exchange-correlation energy
        Eh   = integrate(g, rho_in.*vh/2)                       # Hartree energy
        Epot = integrate(g, rho_in.*vext)                       # potential energy (external)
        Ekin = Etot - Epot - Eh - Ex - Ec                       # kinetic energy
        
        # functions for visualization
        rho_free_ = r -> sum([at[:rho](r) for at in atoms])     # free atom density
        rho_      = r -> (abs2(basis(r)'*efs)*f)[]              # total density
        vh_free_  = r -> sum([at[:vh](r) for at in atoms])      # free atom hartree potential
        vh_       = r -> vh_free_(r) + sum(δvh_mp(r))           # hartree potential
        vx_       = r -> xc.x.v([rho_(r)])[]                    # exchange potential
        vc_       = r -> xc.c.v([rho_(r)])[]                    # correlation potential
        vxc_      = r -> vx_(r) + vc_(r)                        # XC potential
        vhxc_     = r -> vh_(r) + vxc_(r)                       # hartree + XC potential
        ex_       = r -> xc.x.e([rho(r)])[]                     # exchange energy density
        ec_       = r -> xc.c.e([rho(r)])[]                     # correlation energy density
        exc_      = r -> ex_(r) + ec_(r)                        # XC energy density
        
    end  # total

    struct(Etot=Etot,Ekin=Ekin,Epot=Epot,Eh=Eh,Ex=Ex,Ec=Ec,Exc=Exc,Enuc=Enuc,rho=rho_in,
        rho_=rho_,rho_free_=rho_free_,vs=vs,vh=vh,vh_=vh_,vh_free_=vh_free_,vx_=vx_,vc_=vc_,
        vxc_=vxc_,vhxc_=vhxc_,ex_=ex_,ec_=ec_,exc_=exc_,drho_mp=δrho_mp,dvh_mp=δvh_mp,
        evs=evs,efs=efs,f=f,mu=mu,stats=stats)
end


function dft!(m::Molecule; params...)
    merge!(m,Dict{Symbol,Any}(params))

    res = dft(
        m.atoms,                                                # atoms
        m.basis,                                                # molecular basis
        m.g,                                                    # integration grid
        get!(m, :N,      sum([at[:N] for at in m.atoms])),      # particle number
        get!(m, :vembed, zeros(length(m.g))),                   # embedding potential
        get!(m, :xc,     ldaxcunpol()),                         # xc functional
        get!(m, :scf,    SCF()),
        get!(m, :mixer,  PulayMixer()),                         # density mixer
        get!(m, :lmax,   4),                                    # max l-quantum number in hartree multipole expansion
        get!(m, :nh,     300),                                  # number of hartree radial points
        get!(m, :nspl,   100),                                  # number of spline points
        get!(m, :T,      1e-6),                                 # temperature
        get!(m, :logger, Void)
    )

    merge!(m,res)
    return
end


function dft(m::Molecule; params...)
    let m = deepcopy(m)
        dft!(m,params...)
        return m
    end
end

function dft(; name::AbstractString="", atoms::Atoms=Atom[], params...)
    m = Molecule(;name=name,atoms=atoms,params...)
    dft!(m)
    m
end






# sparse version
# function dft!(sys::Molecule; params...)
#     local E,Etot,Enuc,vs,vh,δρ_mp,δvh_mp,evs,efs,f,mu

#     merge!(sys, Dict{Symbol,Any}(params))

#     atoms = sys[:atoms]     ::Atoms
#     basis = sys[:basis]     ::Basis             # molecular basis
#     g     = sys[:g]         ::MultiGrid         # integration grid
#     stats = Stats()

#     @defaults sys begin
#         N::Float64           = :N in keys(sys) ? at[:N] : sum([at[:N] for at in atoms])  # number of particles
#         xc::ldaxc            = ldaxcunpol()                     # exchange-correlation functional
#         scf::SCF             = SCF()
#         mixer::AbstractMixer = PulayMixer()                     # density mixing
#         lmax::Int            = 4                                # max l-quantum number in hartree multipole expansion
#         num::Int             = 300                              # num of radial grid points to solve hartree potential
#         T::Float64           = 1e-6                             # temperature
#         logger               = nothing
#     end

#     if :grid in fieldnames(mixer) mixer.grid = g end
#     scf.grid = g
#     scf.converged = false

#     @stats "* total" begin

#         # initialize
#         @stats "* evaluate vext"                vext     = sum([at[:v](g) for at in atoms])                 # external potential
#         @stats "* evaluate rho_free"            rho_free = sum([at[:rho](g) for at in atoms])               # free atom density
#         @stats "* evaluate vh_free"             vh_free  = sum([at[:vh](g) for at in atoms])                # free atom hartree potential
#                                                 rho_in   = isempty(scf.rho) ? rho_free : scf.rho[end]       # initial guess for density

#         @stats "* evaluate basis"               φ,∆φ = eval_basis_v1(basis,g)                               # evaluate basis on grid
#         @stats "* overlap"                      ζ     = abs(φ) .< 1e-12
#         @stats "* overlap"                      φ[ζ]  = 0.
#         @stats "* overlap"                      ∆φ[ζ] = 0.
#         @stats "* overlap"                      φ     = sparse(φ)
#         @stats "* overlap"                      ∆φ    = sparse(∆φ,φ)
#         @stats "* overlap"                      wφT   = sparse((g.w.*φ)')
#         @stats "* overlap"                      S     = Symmetric(full(wφT*φ))

#         # storage
#         vx = similar(rho_in); vc = similar(rho_in)
#         ex = similar(rho_in); ec = similar(rho_in)
#         vh = similar(rho_in); vs = similar(rho_in)

#         while !scf.converged
#                                                 vh[:] = vh_free
#                                                 # δrho = rho_in - rho_free
#             # @stats "hartree"                    δρ_mp,δvh_mp = poissonmultipole!(vh,g,δrho,atoms;lmax=lmax,num=num)  # Hartree potential

#             @stats "exchange-correlation"       xc(rho_in,ex,vx,ec,vc)                              # exchange-correlation
#                                                 vs[:] = vext + vx + vc + vh                         # Kohn-Sham potential

#             @stats "hamiltonian"                hφ = h(φ,∆φ,vs)                                     # apply Hamiltonian to basis
#             @stats "hamiltonian"                H = Symmetric(full(wφT*hφ))                         # Hamiltonian matrix
#             @stats "KS equations"               evs,efs = LinAlg.eig(H,S)                           # compute eigenspectrum

#             @stats "KS equations"               efs ./= sqrt(integrate(g,abs2(φ*efs)))              # normalize eigenfunctions

#                                                 f,mu = compute_occupations(evs,2,N,T)

#             @stats "compute density"            rho_out = abs2(φ*efs)*f                             # output density
#                                                 E = dot(f,evs)                                      # sum of eigenvalues
#                                                 Enuc = sum([μ != ν ? μ[:Z]*ν[:Z]/dist(μ[:R],ν[:R])/2 : 0.           # nuclear attraction
#                                                     for μ in atoms, ν in atoms])
#             @stats "harris energy"              Etot = E + integrate(g, rho_in.*(ex+ec - (vx+vc) - vh/2)) + Enuc    # Harris functional

#             @stats "check convergence"          scf(Etot, rho_out)                                  # check convergence
#             @stats "density mixing"             rho_in = mixer(rho_in, rho_out)                     # mix densities

#             if !is(logger, nothing) logger(scf) end

#         end  # scf cycle

#         Ex = integrate(g, rho_in.*ex)                           # exchange energy
#         Ec = integrate(g, rho_in.*ec)                           # correlation energy
#         Exc = Ex + Ec                                           # exchange-correlation energy
#         Eh = integrate(g, rho_in.*vh/2)                         # Hartree energy
#         Epot = integrate(g, rho_in.*vext)                       # potential energy (external)
#         Ekin = Etot - Epot - Eh - Ex - Ec                       # kinetic energy

#         # functions for visualization
#         rho_      = r -> (abs2(basis(r)'*efs)*f)[]              # total density
#         rho_free_ = r -> sum([at[:rho](r) for at in atoms])     # free atom density

#     end  # total

#     res = struct(Etot=Etot,Ekin=Ekin,Epot=Epot,Eh=Eh,Ex=Ex,Ec=Ec,Exc=Exc,Enuc=Enuc,rho=rho_in,
#         rho_=rho_,rho_free_=rho_free_,vs=vs,vh=vh,evs=evs,efs=efs,f=f,mu=mu,stats=stats)

#     merge!(sys, res)
#     return
# end




