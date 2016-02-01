

doc"""
`CoulombPotential(Z;R=(0.,0.,0.))`

represents the Coulomb potential $v(r) = -Z/|r-R|$
"""
type CoulombPotential <: ScalarFunc3d{Float64}
    Z::Float64
    R::Point3{Float64}
end

CoulombPotential(Z;R=(0.,0.,0.)) = CoulombPotential(Z,R)

function Base.call(f::CoulombPotential, x::Float64, y::Float64, z::Float64)
    x0,y0,z0 = f.R
    r = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)
    -f.Z/r
end

Base.call(f::CoulombPotential, p::Point3{Float64}) = f(p[1],p[2],p[3])
@vectorize2 Base.call CoulombPotential Point3{Float64} Float64



function setposition(f::CoulombPotential, R::Point3{Float64})
    f.R = R
end
