

# conversion between cartesian and spherical coordinate systems



# spherical coordinates relative to origin

@fastmath function cart2sph(x::Float64, y::Float64, z::Float64)
    xy2 = x*x + y*y
    r = sqrt(xy2 + z*z)  # r
    theta = atan2(sqrt(xy2), z)  # theta (polar angle)
    phi = atan2(y,x)  # phi (azimuthal angle)
    r,theta,phi
end

cart2sph(p::Point3{Float64}) = cart2sph(p[1],p[2],p[3])


# spherical coordinates reference to origin p0=(x0,y0,z0)

cart2sph(x::Float64, y::Float64, z::Float64, x0::Float64, y0::Float64, z0::Float64) = cart2sph(x-x0,y-y0,z-z0)
cart2sph(x::Float64, y::Float64, z::Float64, p0::Point3{Float64}) = cart2sph(x-p0[1],y-p0[2],z-p0[3])
cart2sph(p::Point3{Float64}, p0::Point3{Float64}) = cart2sph(p[1]-p0[1],p[2]-p0[2],p[3]-p0[3])


function cart2sph(dst::AbstractArray{Point3{Float64}}, c::AbstractArray{Point3{Float64}}, p0::Point3{Float64})
    @assert size(dst) == size(c)
    for i in eachindex(c)
        dst[i] = cart2sph(c[i],p0)
    end
end

function cart2sph(c::AbstractArray{Point3{Float64}}, p0::Point3{Float64})
    dst = Array{Point3{Float64}}(size(c))
    cart2sph(dst,c,p0)
    dst
end

@vectorize cart2sph Point3{Float64} Point3{Float64}




@fastmath function sph2cart(r::Float64, θ::Float64, ϕ::Float64)
    x = r*sin(θ)*cos(ϕ)
    y = r*sin(θ)*sin(ϕ)
    z = r*cos(θ)
    x,y,z
end

sph2cart(p::Point3{Float64}) = sph2cart(p[1],p[2],p[3])



@vectorize sph2cart Point3{Float64} Point3{Float64}


