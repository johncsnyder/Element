


function realsphharmbasis(Ω::AngularGrid, q::Vector{Point2{Int}})
    Float64[realsphharm(l,m,θ,ϕ) for (r,θ,ϕ) in cart2sph(Ω), (l,m) in q]
end

# precompute spherical harmonics basis to use in multipole expansion
function realsphharmbasis!(at::Atom, q::Vector{Point2{Int}})
    Ω = at[:g].Ω
    at[:qξ] = q
    at[:ξ]  = reshape(realsphharmbasis(Ω,q), (length(Ω),1,length(q)))
    return
end

realsphharmbasis!(at::Atom, lmax::Int) = realsphharmbasis!(at,enumeratelm(lmax))

function realsphharmbasis!(atoms::Atoms, q::Vector{Point2{Int}})
    for at in atoms
        realsphharmbasis!(at,q)
    end
end

realsphharmbasis!(atoms::Atoms, lmax::Int) = realsphharmbasis!(atoms,enumeratelm(lmax))


function expandmultipole(g    ::AtomicGrid,                         # atomic grid
                         f    ::AbstractArray,                      # function to expand
                         ξ    ::AbstractArray,                      # spherical harmonics basis
                         q    ::Vector{Point2{Int}};                # (l,m) quantum numbers
                         R    ::Point3{Float64} = (0.,0.,0.),       # multipole center
                         nspl ::Int             = 100,              # number of spline points
                         trim ::Bool            = true)             # trim negligible components
    n,m = length(g.r),length(q)
    f_mp = reshape(integrate(g.Ω, f.*ξ), (n,m))
    basis = MultipoleBasis()
    for (i,(l,m)) in enumerate(q)
        if !trim || sumabs2(f_mp[:,i]) > 1e-14
            push!(basis, MultipoleFunc(g.r,f_mp[:,i],l=l,m=m,R=R,npts=nspl))
        end
    end
    basis
end

function expandmultipole(at   ::Atom,                               # atomic grid
                         f    ::AbstractArray;                      # function to expand
                         nspl ::Int  = 100,                         # number of spline points
                         trim ::Bool = true)                        # trim negligible components
    expandmultipole(at[:g],f,at[:ξ],at[:qξ],R=at.R,nspl=nspl,trim=trim)
end

function expandmultipole(at   ::Atom,                               # atomic grid
                         f    ::AbstractVector;                     # function to expand
                         nspl ::Int  = 100,                         # number of spline points
                         trim ::Bool = true)                        # trim negligible components
    expandmultipole(at[:g],reshape(f,size(at[:g])),at[:ξ],at[:qξ],R=at.R,nspl=nspl,trim=trim)
end


function expandmultipole(g     ::MultiGrid,                    # molecular grid
                         rho   ::AbstractArray,                # electronic density
                         atoms ::Atoms;                        # atoms
                         nspl  ::Int = 100,                    # number of spline points
                         trim ::Bool = true)
    rho_mp = MultipoleBasis[]                                  # rho multipole expansion
    for (i,j) in enumerate(eachatom(g))
        rho_at = expandmultipole(atoms[i],g.p[j].*rho[j],nspl=nspl,trim=trim)
        push!(rho_mp, rho_at)
    end
    rho_mp
end


# function expandmultipole(m     ::Molecule;
#                          nspl  ::Int = 100,                    # number of spline points
#                          trim ::Bool = true)
#     expandmultipole(m.g,m[:rho],m.atoms,nspl=nspl,trim=trim)
# end


function poissonmultipole(rh   ::LogarithmicGrid,               # radial grid to solve poisson equation
                          rho  ::MultipoleFuncSpl;              # multipole component of density
                          R    ::Point3f = (0.,0.,0.),          # multipole center
                          nspl ::Int     = 100)                 # number of spline points
    l,m  = rho.l,rho.m
    rcut = rho.u.x[end]                                         # radial extent of splined density
    f  = rho.u(rh)                                              # evaluate radial part of multipole component
    vh = 4π/(2l+1)*radpoisson(rh,f,l)                           # multipole componenet of hartree potential
    q  = 4π/(2l+1)*multipolemoment(rh,f,l)                      # multipole moment
    MultipoleFunc(rh,vh,l=l,m=m,R=R,rcut=rcut,q=q,npts=nspl)
end


function poissonmultipole(rho  ::AbstractArray,
                          at   ::Atom;
                          nh   ::Int = 300,
                          nspl ::Int = 100)
    rho_mp = expandmultipole(at,rho,nspl=nspl)                  # multipole expansion of density
    rh     = loggrid(at[:g].r,nh)                               # radial grid to solve poisson equation
    vh_mp  = MultipoleFuncSpl[poissonmultipole(rh,rho_mp_lm,R=at.R,nspl=nspl) for rho_mp_lm in rho_mp]      # multipole expansion of hartree potential
    rho_mp, vh_mp
end


function poissonmultipole(g     ::MultiGrid,                    # molecular grid
                          rho   ::AbstractArray,                # electronic density
                          atoms ::Atoms;                        # atoms
                          nh    ::Int = 300,                    # number of hartree radial points
                          nspl  ::Int = 100)                    # number of spline points
    rho_mp = MultipoleBasis()                                   # rho multipole expansion for each atom
    vh_mp  = MultipoleBasis()                                   # vh multipole expansion for each atom
    for (i,j) in enumerate(eachatom(g))
        rho_at, vh_at = poissonmultipole(g.p[j].*rho[j],atoms[i],nh=nh,nspl=nspl)
        append!(rho_mp, rho_at)
        append!(vh_mp, vh_at)
    end
    rho_mp, vh_mp
end




