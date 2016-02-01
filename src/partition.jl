


# type PartitionFunc <: ScalarFunc3d{Float64}
#     at::Atom
#     atoms::Vector{Atom}
# end

# function Base.call(f::PartitionFunc, x::Real, y::Real, z::Real)
#     f.at.h(x,y,z)/sum([μ.h(x,y,z) for μ in f.atoms])
# end


# Base.call(f::PartitionFunc, p::Point3{Float64}) = f(p[1],p[2],p[3])

# @vectorize2 Base.call PartitionFunc Point3{Float64} Float64


function delleypartition!(atoms::Atoms)
    for at in atoms
        r = at[:g].r
        at[:h] = RadialFunc(r,at[:rho].u(r)./r.^2,R=at.R,rmax=at[:rho].rmax,npts=length(r))
    end
    for at in atoms
        at[:p] = r -> at[:h](r)./sum([μ[:h](r) for μ in atoms])
    end
    return
end


# function delleypartition!(atoms::Atoms)
#     for at in atoms
#         r,rho = at[:rho].u.x,at[:rho].u.y
#         at[:h] = RadialFunc(r,rho./r.^2,R=at[:R],rmax=at[:rho].rmax,npts=length(r))
#     end
#     for at in atoms
#         at[:p] = r -> at[:h](r)./sum([μ[:h](r) for μ in atoms])
#     end
#     return
# end

function delleypartition(atoms::Atoms)
    let atoms = copy(atoms)
        delleypartition!(atoms)
        return atoms
    end
end

