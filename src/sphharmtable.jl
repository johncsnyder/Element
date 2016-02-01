

# l,m -> plm array index
ip(l,m) = ( l*(l+1) >> 1 )+1+m  # bitshift >> 1 divides by 2

# l,m -> ylm array index
iq(l,m) = l*(l+1)+1+m



type SphericalHarmonicTable
    lmax::Int
    p::Array{Float64}  # associated legendre polynomial P_l^m(cos(θ)) / sin(θ)^m
    q::Array{Float64}  # spherical harmonic Y_{lm}(θ,ϕ)
    c::Array{Float64}  # zenith angular factor cos(mϕ) * sin(θ)^m
    s::Array{Float64}  # zenith angular factor sin(mϕ) * sin(θ)^m
end


function SphericalHarmonicTable(lmax)
    p = Array{Float64}(ip(lmax,lmax))
    q = Array{Float64}(iq(lmax,lmax))
    c = Array{Float64}(lmax+1)
    s = Array{Float64}(lmax+1)
    c[1] = 1.
    s[1] = 0.
    p[1] = 1.
    q[1] = 1/sqrt(4π)
    SphericalHarmonicTable(lmax,p,q,c,s)
end





function computesphharm(l,m,p,q,c,s)
    if m == 0
        q[iq(l,m)] = N(l,m)*p[ip(l,m)]
    else  # m > 0
        z = sqrt(2)*N(l,m)*p[ip(l,m)]
        q[iq(l,m)] = (-1)^m*c[m+1]*z
        q[iq(l,-m)] = -s[m+1]*z
    end
end




# compute sqrt((2l+1)*(l-m)! / (4π*(l+m)!))
function N(l::Int, m::Int)  # m > 0
    f = 1.
    for k in l-m+1:l+m; f *= k end  # f = factorial(l-m) / factorial(l+m)
    sqrt( (2l+1) / (4π*f) )
end



function Base.call(Y::SphericalHarmonicTable, x::Float64, y::Float64, z::Float64)
    lmax,p,q,c,s = Y.lmax,Y.p,Y.q,Y.c,Y.s

    for m in 0:lmax-1
        p[ip(m+1,m+1)] = -(2m+1)*p[ip(m,m)]  # P_{m+1,m+1} = -(2m+1)*P_{mm}
        s[m+2] = x*s[m+1] + y*c[m+1]
        c[m+2] = x*c[m+1] - y*s[m+1]
        computesphharm(m+1,m+1,p,q,c,s)
        p[ip(m+1,m)] = (2m+1)*z*p[ip(m,m)]  # P_{m+1,m} = (2m+1)*z*P_{mm}
        computesphharm(m+1,m,p,q,c,s)
        for l in m+1:lmax-1
            p[ip(l+1,m)] = ( (2l+1)*z*p[ip(l,m)] - (l+m)*p[ip(l-1,m)] ) / (l-m+1)
            computesphharm(l+1,m,p,q,c,s)
        end
    end
end


Base.call(Y::SphericalHarmonicTable, p::Point3{Float64}) = Y(p[1],p[2],p[3])


Base.getindex(Y::SphericalHarmonicTable, l::Int, m::Int) = Y.q[iq(l,m)]





