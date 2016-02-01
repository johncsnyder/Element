


type NAO{T} <: ScalarFunc3d{Float64}
    u::T
    l::Int
    m::Int
    R::Point3{Float64}
    rmax::Float64
end

NAO(u; l=0, m=0, R=(0.,0.,0.), rmax=Inf) = NAO(u,l,m,R,rmax)

NAO(r::AbstractVector, u::AbstractVector; l=0, m=0, R=(0.,0.,0.),
        rmax=Inf, npts::Int=100) =
    NAO(Spline1D(r,u,npts),l,m,R,rmax)


function Base.call(f::NAO, x::Real, y::Real, z::Real)
    r,θ,ϕ = cart2sph(x,y,z,f.R)
    r > f.rmax ? 0. : f.u(r)*realsphharm(f.l,f.m,θ,ϕ)
end



Base.call(f::NAO, p::Point3{Float64}) = f(p[1],p[2],p[3])
@vectorize2 Base.call NAO Point3{Float64} Float64



function laplacian(f::NAO{Spline1D{Float64,Vector{Float64}}}, x::Float64, y::Float64, z::Float64)
    r,θ,ϕ = cart2sph(x,y,z,f.R)
    r > f.rmax && return 0.
    u,ud,ud2 = eval_deriv12(f.u,r)
    ( r*(2*ud + r*ud2) - u*f.l*(f.l+1) ) * realsphharm(f.l,f.m,θ,ϕ) / (r*r)
end


laplacian(f::NAO, p::Point3{Float64}) = laplacian(f,p[1],p[2],p[3])
@vectorize2 laplacian NAO Point3{Float64} Float64


function setposition(f::NAO, R::Point3{Float64})
    f.R = R
end



typealias NAOSpl NAO{Spline1D{Float64,Vector{Float64}}}
typealias NAOBasis Vector{NAOSpl}


# function eval_basis_v1(b,g)
#     n,m = length(g),length(b)
#     φ = zeros(n,m)
#     ∆φ = zeros(n,m)
    
#     for j in 1:m
#         f = b[j]        ::NAO{Spline1D{Float64,Vector{Float64}}}

#         for i in 1:n
#             r,θ,ϕ = cart2sph(g[i],f.R)
#             if r <= f.rmax
#                 u,ud,ud2 = eval_deriv12(f.u,r)
#                 ylm      = realsphharm(f.l,f.m,θ,ϕ)
#                 φ[i,j]   = u*ylm
#                 ∆φ[i,j]  = ( r*(2*ud + r*ud2) - u*f.l*(f.l+1) ) * ylm / (r*r)
#             end
#         end
#     end
    
#     φ,∆φ
# end



# function eval_basis_v2(atoms,g)
#     n  = length(g)
#     m  = sum(Int[length(at[:basis]) for at in atoms])  # total number of basis functions
#     φ  = zeros(m,n)
#     ∆φ = zeros(m,n)

#     lmax = maximum(Int[maximum(Int[f.l for f in at[:basis]]) for at in atoms])
#     # lmax = maximum(Int[f.l for f in b])
#     Y = SphericalHarmonicTable(lmax)
    
#     for j in 1:n
#         p = g[j]    # grid point
#         i = 1       # overall basis function index

#         for at in atoms
#             R     = at[:R]       ::Point3{Float64}                                      # atom position
#             b     = at[:basis]   ::Vector{NAO{Spline1D{Float64,Vector{Float64}}}}       # atomic basis
            
#             q     = p-R         # grid point relative to atom position
#             r     = norm(q)     # distance from grid point to atom position
#             invr2 = 1/(r*r)
#             Y(q/r)              # precompute spherical harmonics
            
#             for k in 1:length(b)
#                 f = b[k]        ::NAO{Spline1D{Float64,Vector{Float64}}}        # basis function

#                 if r <= f.rmax
#                     u,ud,ud2 = eval_deriv12(f.u,r)
#                     Ylm      = Y[f.l,f.m]
#                      φ[i,j]  = u*Ylm
#                     ∆φ[i,j]  = ( r*(2*ud + r*ud2) - u*f.l*(f.l+1) ) * invr2 * Ylm
#                 end
#                 i += 1
#             end
#         end
#     end
    
#     φ',∆φ'
# end


using Iterators: groupby



function eval_basis_sum(b::MultipoleBasis, g::AbstractGrid3d)
    n_g  = length(g)
    φ    = zeros(n_g)                           # basis on grid
    lmax = maximum(Int[f.l for f in b])         # max l-quantum number
    Y    = SphericalHarmonicTable(lmax)         # spherical harmonics table
    ξ    = collect(groupby(f->f.R, b))          # partition basis grouped by center
    χ    = Point3f[first(b).R for b in ξ]       # list of centers
    n_χ  = length(χ)
    
    for i in 1:n_g
        p = g[i]                              # grid point
        for j in 1:n_χ
            R = χ[j]
            q = p-R                             # grid point relative to atom position
            r = norm(q)                         # distance from grid point to atom position
            Y(q/r)                              # precompute spherical harmonics
            
            for f in ξ[j]
                if r <= f.rmax
                    φ[i] += f.u(r)*Y[f.l,f.m]
                end
            end
        end
    end
    
    φ
end


function eval_basis(b::Union{NAOBasis,MultipoleBasis}, g::AbstractGrid3d)
    n_g,n_b = length(g),length(b)
    φ    = zeros(n_b,n_g)                       # basis on grid
    lmax = maximum(Int[f.l for f in b])         # max l-quantum number
    Y    = SphericalHarmonicTable(lmax)         # spherical harmonics table
    ξ    = collect(groupby(f->f.R, b))          # partition basis grouped by center
    χ    = Point3f[first(b).R for b in ξ]       # list of centers
    n_χ  = length(χ)
    
    for i_g in 1:n_g
        p = g[i_g]                              # grid point
        i_b = 1                                 # overall basis function index
        for i_χ in 1:n_χ
            R = χ[i_χ]
            q = p-R                             # grid point relative to atom position
            r = norm(q)                         # distance from grid point to atom position
            Y(q/r)                              # precompute spherical harmonics
            
            for f in ξ[i_χ]
                if r <= f.rmax
                    Ylm = Y[f.l,f.m]
                    φ[i_b,i_g] = f.u(r)*Ylm
                end
                i_b += 1
            end
        end
    end
    
    φ'
end

function eval_basis_lapl(b::NAOBasis, g::AbstractGrid3d)
    n_g,n_b = length(g),length(b)
    φ    = zeros(n_b,n_g)                       # basis on grid
    ∆φ   = zeros(n_b,n_g)                       # basis laplacian on grid
    lmax = maximum(Int[f.l for f in b])         # max l-quantum number
    Y    = SphericalHarmonicTable(lmax)         # spherical harmonics table
    ξ    = collect(groupby(f->f.R, b))          # partition basis grouped by center
    χ    = Point3f[first(b).R for b in ξ]       # list of centers
    n_χ  = length(χ)
    
    for i_g in 1:n_g
        p = g[i_g]                              # grid point
        i_b = 1                                 # overall basis function index
        for i_χ in 1:n_χ
            R = χ[i_χ]
            q = p-R                             # grid point relative to atom position
            r = norm(q)                         # distance from grid point to atom position
            invr2 = 1/(r*r)
            Y(q/r)                              # precompute spherical harmonics
            
            for f in ξ[i_χ]
                if r <= f.rmax
                    u,ud,ud2 = eval_deriv12(f.u,r)
                    Ylm = Y[f.l,f.m]
                     φ[i_b,i_g] = u*Ylm
                    ∆φ[i_b,i_g] = ( r*(2*ud + r*ud2) - u*f.l*(f.l+1) ) * invr2 * Ylm
                end
                i_b += 1
            end
        end
    end
    
    φ',∆φ'
end





# function eval_basis(b::Union{NAOBasis,MultipoleBasis}, g::AbstractGrid3d, sphharmtab::AbstractArray)
#     n_g,n_b = length(g),length(b)
#     φ    = zeros(n_b,n_g)                       # basis on grid
#     ξ    = collect(groupby(f->f.R, b))          # partition basis grouped by center
#     χ    = Point3f[first(b).R for b in ξ]       # list of centers
#     n_χ  = length(χ)
    
#     for i_g in 1:n_g
#         p = g[i_g]                              # grid point
#         i_b = 1                                 # overall basis function index
#         for i_χ in 1:n_χ
#             R = χ[i_χ]
#             r = norm(p-R)                       # grid point relative to atom position
#             for f in ξ[i_χ]
#                 if r <= f.rmax
#                     Ylm = sphharmtab[i_g,i_χ,iq(f.l,f.m)]
#                     φ[i_b,i_g] = f.u(r)*Ylm
#                 end
#                 i_b += 1
#             end
#         end
#     end
    
#     φ'
# end



function compute_sphharm_dist_tab(b::Union{NAOBasis,MultipoleBasis}, g::AbstractGrid3d, lmax::Int)
    Y = SphericalHarmonicTable(lmax)
    ξ = unique(Point3f[f.R for f in b])  # assume basis is sorted
    sphharmtab = zeros(length(g),length(ξ),length(Y.q));
    disttab = zeros(length(g),length(ξ))

    for i in 1:length(g)
        p = g[i]
        for j in 1:length(ξ)
            R = ξ[j]
            d = p-R
            r = norm(d)
            Y(d/r)
            sphharmtab[i,j,:] = Y.q
            disttab[i,j] = r
        end
    end

    sphharmtab, disttab
end


# function eval_basis(b::Union{NAOBasis,MultipoleBasis}, g::AbstractGrid3d, sphharmtab::AbstractArray, disttab::AbstractArray)
#     n,m = length(g),length(b)
#     φ = zeros(n,m)
#     i = 1; R = b[1].R
#     ind = Int[(f.R != R ? begin R = f.R; i += 1 end : i) for f in b]
    
#     for j in 1:m
#         f = b[j]
#         k = iq(f.l,f.m)
#         R = f.R
#         for i in 1:n
#             r = disttab[i,ind[j]]
#             if r <= f.rmax
#                 Ylm = sphharmtab[i,ind[j],k]
#                 φ[i,j] = f.u(r)*Ylm
#             end
#         end
#     end
    
#     φ
# end


function eval_basis(b::NAOBasis, g::AbstractGrid3d, sphharmtab::AbstractArray, disttab::AbstractArray)
    n,m = length(g),length(b)
    φ = zeros(n)
    i = 1; R = b[1].R
    ind = Int[(f.R != R ? begin R = f.R; i += 1 end : i) for f in b]
    
    for j in 1:m
        f = b[j]
        k = iq(f.l,f.m)
        R = f.R
        for i in 1:n
            r = disttab[i,ind[j]]
            if r <= f.rmax
                φ[i] += f.u(r)*sphharmtab[i,ind[j],k]
            end
        end
    end
    
    φ
end



# function eval_basis_v2(b,g,atoms)
#     n,m  = length(g),length(b)
#     φ  = zeros(m,n)
#     ∆φ = zeros(m,n)
#     lmax = maximum(Int[f.l for f in b])
#     Y = SphericalHarmonicTable(lmax)
    
#     for j in 1:n
#         p = g[j]    # grid point
#         i = 1       # overall basis function index

#         for at in atoms
#             R     = at[:R]       ::Point3{Float64}                                      # atom position
#             b     = at[:basis]   ::Vector{NAO{Spline1D{Float64,Vector{Float64}}}}       # atomic basis
#             w
#             q     = p-R         # grid point relative to atom position
#             r     = norm(q)     # distance from grid point to atom position
#             invr2 = 1/(r*r)
#             Y(q/r)              # precompute spherical harmonics
            
#             for k in 1:length(b)
#                 f = b[k]        ::NAO{Spline1D{Float64,Vector{Float64}}}        # basis function

#                 if r <= f.rmax
#                     u,ud,ud2 = eval_deriv12(f.u,r)
#                     Ylm      = Y[f.l,f.m]
#                      φ[i,j]  = u*Ylm
#                     ∆φ[i,j]  = ( r*(2*ud + r*ud2) - u*f.l*(f.l+1) ) * invr2 * Ylm
#                 end
#                 i += 1
#             end
#         end
#     end
    
#     φ',∆φ'
# end



# function eval_basis_v3(b,g,atoms)
#     n,m,l = length(g),length(b),length(atoms)
#     φ = zeros(m,n)
#     ∆φ = zeros(m,n)
#     lmax = maximum(Int[f.l for f in b])
#     Y,p,c,s = ylm_table_init(lmax)
    
#     for j in 1:n
#         q = g[j]
#         i = 1
#         for at in atoms
#             R = at[:R]::Point3{Float64}
#             d = q-R
#             r = norm(d)
#             invr2 = 1/(r*r)
#             x,y,z = d/r
#             b_at = at[:basis]::Vector{ScalarFunc3d{Float64}}
#             ylm_table!(Y,p,c,s,x,y,z,lmax)
#             for k in 1:length(b_at)
#                 f = b_at[k]::NAO{Spline1D{Float64,Vector{Float64}}}
#                 if r <= f.rmax
#                     u,ud,ud2 = eval_deriv12(f.u,r)
#                     ylm = Y[i_ylm(f.l,f.m)]
#                     φ[i,j] = u*ylm
#                     ∆φ[i,j] = ( r*(2*ud + r*ud2) - u*f.l*(f.l+1) ) * invr2 * ylm
#                 end
#                 i += 1
#             end
#         end
#     end
    
#     φ',∆φ'
# end
