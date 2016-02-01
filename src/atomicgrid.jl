


type AtomicGrid{T,A<:AbstractMatrix,B<:AbstractMatrix,C<:AbstractGrid1d} <: AbstractGrid3d{T,2}
    c::A  # coords
    w::B  # weights
    r::C  # radial grid
    Ω::AngularGrid{T}  # angular grid
    R::Point3{T}  # origin of spherical coord system
    function AtomicGrid(c::AbstractMatrix{Point3{T}},w::AbstractMatrix{T},
            r::AbstractGrid1d{T,1},Ω::AngularGrid{T},R::Point3{T})
        @assert size(c) == size(w)
        new(c,w,r,Ω,R)
    end
end

AtomicGrid{T}(c::AbstractMatrix{Point3{T}},w::AbstractMatrix{T},
        r::AbstractGrid1d{T,1},Ω::AngularGrid{T},R::Point3{T}) = 
    AtomicGrid{eltype(w),typeof(c),typeof(w),typeof(r)}(c,w,r,Ω,R)




function atomicgrid{T}(r::AbstractGrid1d{T,1}, Ω::AngularGrid{T}; R::Point3{Float64}=(0.,0.,0.))
    shp = (length(Ω),length(r))
    c = Array{Point3{Float64}}(shp)
    w = Array{Float64}(shp)
    for j in 1:length(r)
        for i in 1:length(Ω)
            c[i,j] = r[j]*Ω[i] + R
            w[i,j] = r[j]^2*r.w[j]*Ω.w[i]
        end
    end
    AtomicGrid(c,w,r,Ω,R)
end






# integrate_angular{T}(g::AtomicGrid{T}, f::AbstractMatrix{T}) = At_mul_B(f,g.Ω.w)


# integrate_angular{T}(g::AtomicGrid{T}, f::AbstractVector{T}) = reshape(f,size(g.w))'*g.Ω.w


# function integrate_angular{T}(g::AtomicGrid{T}, f::AbstractArray{T})
#     shp = (size(g.w,1),prod(size(f)[2:end]))
#     reshape(f,shp)'*g.Ω.w
# end







# type AtomicGrid{T,N,Tc<:AbstractArray,Tw<:AbstractArray,Tr<:AbstractGrid1d} <: AbstractGrid3d{T,N}
#     c::Tc  # coords
#     w::Tw  # weights
#     r::Tc  # radial grid
#     Ω::AngularGrid{T}  # angular grid
#     R::Point3{T}  # origin of spherical coord system
#     function AtomicGrid(c::AbstractArray{Point3{T},N},
#                         w::AbstractArray{T,N},
#                         r::AbstractGrid1d{T,1},
#                         Ω::AngularGrid{T},
#                         R::Point3{T})
#         @assert size(c) == size(w)
#         new(c,w,r,Ω,R)
#     end
# end

# function AtomicGrid{T,N}(c::AbstractArray{Point3{T},N},
#                          w::AbstractArray{T,N},
#                          r::AbstractGrid1d{T,1},
#                          Ω::AngularGrid{T},
#                          R::Point3{T})
#     AtomicGrid{eltype(w),N,typeof(c),typeof(w),typeof(r)}(c,w,r,Ω,R)
# end







