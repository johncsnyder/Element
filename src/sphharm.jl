



# using GSL.sf_legendre_sphPlm


# function sphharm(l::Integer, m::Integer, theta::Real, phi::Real)
# 	@assert l >= 0 && abs(m) <= l
# 	if m >= 0
# 		return exp(m*phi*im)*sf_legendre_sphPlm(l,m,cos(theta))
# 	else
#     	return (-1.)^m*exp(m*phi*im)*sf_legendre_sphPlm(l,-m,cos(theta))
#     end
# end




doc"""
`sphharm(l,m,θ,ϕ)`

computes the spherical harmonic function (standard definition used in quantum mechanics)

$ Y_l^m(\theta,\phi) = (-1)^m \sqrt{{2l+1 \over 4\pi} {(l-m)!\over (l+m)!}} P_l^m(cos \theta) e^{im\phi} $

which are normalized such that 

$ \int Y_l^m(\theta,\phi)^* Y_{l'}^{m'}(\theta,\phi) d\Omega = \delta_{ll'} \delta_{mm'} $
"""
# sphharm(l::Int, m::Int, θ::Real, ϕ::Real) = sphharmtab[l,m](θ,ϕ)::Complex{Float64}
sphharm(l::Int, m::Int, θ::Real, ϕ::Real) = sphharmtab[i_lm(l,m)](θ,ϕ)::Complex{Float64}


doc"""
`realsphharm(l,m,θ,ϕ)`

computes the real spherical harmonic function

$ Y_{lm} = {i \over \sqrt{2}} \left(Y_l^m - (-1)^m Y_l^{-m}\right), \quad m<0 $

$ Y_{lm} = Y_l^m, \quad m=0 $

$ Y_{lm} = {1 \over \sqrt{2}} \left(Y_l^{-m} + (-1)^m Y_l^m\right), \quad m>0 $

"""
# realsphharm(l::Int, m::Int, θ::Real, ϕ::Real) = realsphharmtab[l,m](θ,ϕ)::Float64
realsphharm(l::Int, m::Int, θ::Real, ϕ::Real) = realsphharmtab[i_lm(l,m)](θ,ϕ)::Float64


sphharm00(θ::Float64, ϕ::Float64) = complex(1/(2*sqrt(pi)))
sphharm1m1(θ::Float64, ϕ::Float64) = ((1/2)*sqrt(3/(2*pi))*sin(θ))/exp(ϕ*im)
sphharm10(θ::Float64, ϕ::Float64) = complex((1/2)*sqrt(3/pi)*cos(θ))
sphharm11(θ::Float64, ϕ::Float64) = (-(1/2))*exp(ϕ*im)*sqrt(3/(2*pi))*sin(θ)
sphharm2m2(θ::Float64, ϕ::Float64) = ((1/4)*sqrt(15/(2*pi))*sin(θ)^2)/exp(2*ϕ*im)
sphharm2m1(θ::Float64, ϕ::Float64) = ((1/2)*sqrt(15/(2*pi))*cos(θ)*sin(θ))/exp(ϕ*im)
sphharm20(θ::Float64, ϕ::Float64) = complex((1/8)*sqrt(5/pi)*(1 + 3*cos(2*θ)))
sphharm21(θ::Float64, ϕ::Float64) = (-(1/2))*exp(ϕ*im)*sqrt(15/(2*pi))*cos(θ)*sin(θ)
sphharm22(θ::Float64, ϕ::Float64) = (1/4)*exp(2*ϕ*im)*sqrt(15/(2*pi))*sin(θ)^2
sphharm3m3(θ::Float64, ϕ::Float64) = ((1/8)*sqrt(35/pi)*sin(θ)^3)/exp(3*ϕ*im)
sphharm3m2(θ::Float64, ϕ::Float64) = ((1/4)*sqrt(105/(2*pi))*cos(θ)*sin(θ)^2)/exp(2*ϕ*im)
sphharm3m1(θ::Float64, ϕ::Float64) = ((1/16)*sqrt(21/pi)*(3 + 5*cos(2*θ))*sin(θ))/exp(ϕ*im)
sphharm30(θ::Float64, ϕ::Float64) = complex((1/16)*sqrt(7/pi)*(3*cos(θ) + 5*cos(3*θ)))
sphharm31(θ::Float64, ϕ::Float64) = (-(1/16))*exp(ϕ*im)*sqrt(21/pi)*(3 + 5*cos(2*θ))*sin(θ)
sphharm32(θ::Float64, ϕ::Float64) = (1/4)*exp(2*ϕ*im)*sqrt(105/(2*pi))*cos(θ)*sin(θ)^2
sphharm33(θ::Float64, ϕ::Float64) = (-(1/8))*exp(3*ϕ*im)*sqrt(35/pi)*sin(θ)^3
sphharm4m4(θ::Float64, ϕ::Float64) = ((3/16)*sqrt(35/(2*pi))*sin(θ)^4)/exp(4*ϕ*im)
sphharm4m3(θ::Float64, ϕ::Float64) = ((3/8)*sqrt(35/pi)*cos(θ)*sin(θ)^3)/exp(3*ϕ*im)
sphharm4m2(θ::Float64, ϕ::Float64) = ((3/16)*sqrt(5/(2*pi))*(5 + 7*cos(2*θ))*sin(θ)^2)/exp(2*ϕ*im)
sphharm4m1(θ::Float64, ϕ::Float64) = ((3/32)*sqrt(5/pi)*(1 + 7*cos(2*θ))*sin(2*θ))/exp(ϕ*im)
sphharm40(θ::Float64, ϕ::Float64) = complex((3*(9 + 20*cos(2*θ) + 35*cos(4*θ)))/(128*sqrt(pi)))
sphharm41(θ::Float64, ϕ::Float64) = (-(3/32))*exp(ϕ*im)*sqrt(5/pi)*(1 + 7*cos(2*θ))*sin(2*θ)
sphharm42(θ::Float64, ϕ::Float64) = (3/16)*exp(2*ϕ*im)*sqrt(5/(2*pi))*(5 + 7*cos(2*θ))*sin(θ)^2
sphharm43(θ::Float64, ϕ::Float64) = (-(3/8))*exp(3*ϕ*im)*sqrt(35/pi)*cos(θ)*sin(θ)^3
sphharm44(θ::Float64, ϕ::Float64) = (3/16)*exp(4*ϕ*im)*sqrt(35/(2*pi))*sin(θ)^4
sphharm5m5(θ::Float64, ϕ::Float64) = ((3/32)*sqrt(77/pi)*sin(θ)^5)/exp(5*ϕ*im)
sphharm5m4(θ::Float64, ϕ::Float64) = ((3/16)*sqrt(385/(2*pi))*cos(θ)*sin(θ)^4)/exp(4*ϕ*im)
sphharm5m3(θ::Float64, ϕ::Float64) = ((1/64)*sqrt(385/pi)*(7 + 9*cos(2*θ))*sin(θ)^3)/exp(3*ϕ*im)
sphharm5m2(θ::Float64, ϕ::Float64) = ((1/16)*sqrt(1155/(2*pi))*cos(θ)*(1 + 3*cos(2*θ))*sin(θ)^2)/exp(2*ϕ*im)
sphharm5m1(θ::Float64, ϕ::Float64) = ((1/16)*sqrt(165/(2*pi))*(1 - 14*cos(θ)^2 +21*cos(θ)^4)*sin(θ))/exp(ϕ*im)
sphharm50(θ::Float64, ϕ::Float64) = complex((1/256)*sqrt(11/pi)*(30*cos(θ) + 35*cos(3*θ) +63*cos(5*θ)))
sphharm51(θ::Float64, ϕ::Float64) = (-(1/16))*exp(ϕ*im)*sqrt(165/(2*pi))*(1 - 14*cos(θ)^2 + 21*cos(θ)^4)*sin(θ)
sphharm52(θ::Float64, ϕ::Float64) = (1/16)*exp(2*ϕ*im)*sqrt(1155/(2*pi))*cos(θ)*(1 + 3*cos(2*θ))*sin(θ)^2
sphharm53(θ::Float64, ϕ::Float64) = (-(1/64))*exp(3*ϕ*im)*sqrt(385/pi)*(7 + 9*cos(2*θ))*sin(θ)^3
sphharm54(θ::Float64, ϕ::Float64) = (3/16)*exp(4*ϕ*im)*sqrt(385/(2*pi))*cos(θ)*sin(θ)^4
sphharm55(θ::Float64, ϕ::Float64) = (-(3/32))*exp(5*ϕ*im)*sqrt(77/pi)*sin(θ)^5
sphharm6m6(θ::Float64, ϕ::Float64) = ((1/64)*sqrt(3003/pi)*sin(θ)^6)/exp(6*ϕ*im)
sphharm6m5(θ::Float64, ϕ::Float64) = ((3/32)*sqrt(1001/pi)*cos(θ)*sin(θ)^5)/exp(5*ϕ*im)
sphharm6m4(θ::Float64, ϕ::Float64) = ((3/64)*sqrt(91/(2*pi))*(9 + 11*cos(2*θ))*sin(θ)^4)/exp(4*ϕ*im)
sphharm6m3(θ::Float64, ϕ::Float64) = ((1/64)*sqrt(1365/pi)*cos(θ)*(5 + 11*cos(2*θ))*sin(θ)^3)/exp(3*ϕ*im)
sphharm6m2(θ::Float64, ϕ::Float64) = ((1/64)*sqrt(1365/pi)*(1 - 18*cos(θ)^2 +33*cos(θ)^4)*sin(θ)^2)/exp(2*ϕ*im)
sphharm6m1(θ::Float64, ϕ::Float64) = ((1/16)*sqrt(273/(2*pi))*cos(θ)*(5 - 30*cos(θ)^2 +33*cos(θ)^4)*sin(θ))/exp(ϕ*im)
sphharm60(θ::Float64, ϕ::Float64) = complex((1/32)*sqrt(13/pi)*(-5 + 105*cos(θ)^2 -315*cos(θ)^4 + 231*cos(θ)^6))
sphharm61(θ::Float64, ϕ::Float64) = (-(1/16))*exp(ϕ*im)*sqrt(273/(2*pi))*cos(θ)*(5 - 30*cos(θ)^2 + 33*cos(θ)^4)*sin(θ)
sphharm62(θ::Float64, ϕ::Float64) = (1/64)*exp(2*ϕ*im)*sqrt(1365/pi)*(1 - 18*cos(θ)^2 +33*cos(θ)^4)*sin(θ)^2
sphharm63(θ::Float64, ϕ::Float64) = (-(1/64))*exp(3*ϕ*im)*sqrt(1365/pi)*cos(θ)*(5 + 11*cos(2*θ))*sin(θ)^3
sphharm64(θ::Float64, ϕ::Float64) = (3/64)*exp(4*ϕ*im)*sqrt(91/(2*pi))*(9 + 11*cos(2*θ))*sin(θ)^4
sphharm65(θ::Float64, ϕ::Float64) = (-(3/32))*exp(5*ϕ*im)*sqrt(1001/pi)*cos(θ)*sin(θ)^5
sphharm66(θ::Float64, ϕ::Float64) = (1/64)*exp(6*ϕ*im)*sqrt(3003/pi)*sin(θ)^6
sphharm7m7(θ::Float64, ϕ::Float64) = ((3/64)*sqrt(715/(2*pi))*sin(θ)^7)/exp(7*ϕ*im)
sphharm7m6(θ::Float64, ϕ::Float64) = ((3/64)*sqrt(5005/pi)*cos(θ)*sin(θ)^6)/exp(6*ϕ*im)
sphharm7m5(θ::Float64, ϕ::Float64) = ((3/64)*sqrt(385/(2*pi))*(-1 + 13*cos(θ)^2)*sin(θ)^5)/exp(5*ϕ*im)
sphharm7m4(θ::Float64, ϕ::Float64) = ((3/64)*sqrt(385/(2*pi))*cos(θ)*(7 + 13*cos(2*θ))*sin(θ)^4)/exp(4*ϕ*im)
sphharm7m3(θ::Float64, ϕ::Float64) = ((3/64)*sqrt(35/(2*pi))*(3 - 66*cos(θ)^2 +143*cos(θ)^4)*sin(θ)^3)/exp(3*ϕ*im)
sphharm7m2(θ::Float64, ϕ::Float64) = ((3/64)*sqrt(35/pi)*cos(θ)*(15 - 110*cos(θ)^2 +143*cos(θ)^4)*sin(θ)^2)/exp(2*ϕ*im)
sphharm7m1(θ::Float64, ϕ::Float64) = ((1/64)*sqrt(105/(2*pi))*(-5 + 135*cos(θ)^2 -495*cos(θ)^4 + 429*cos(θ)^6)*sin(θ))/exp(ϕ*im)
sphharm70(θ::Float64, ϕ::Float64) = complex((1/32)*sqrt(15/pi)*cos(θ)*(-35 + 315*cos(θ)^2 -693*cos(θ)^4 + 429*cos(θ)^6))
sphharm71(θ::Float64, ϕ::Float64) = (-(1/64))*exp(ϕ*im)*sqrt(105/(2*pi))*(-5 + 135*cos(θ)^2 - 495*cos(θ)^4 + 429*cos(θ)^6)*sin(θ)
sphharm72(θ::Float64, ϕ::Float64) = (3/64)*exp(2*ϕ*im)*sqrt(35/pi)*cos(θ)*(15 - 110*cos(θ)^2 + 143*cos(θ)^4)*sin(θ)^2
sphharm73(θ::Float64, ϕ::Float64) = (-(3/64))*exp(3*ϕ*im)*sqrt(35/(2*pi))*(3 - 66*cos(θ)^2 + 143*cos(θ)^4)*sin(θ)^3
sphharm74(θ::Float64, ϕ::Float64) = (3/64)*exp(4*ϕ*im)*sqrt(385/(2*pi))*cos(θ)*(7 + 13*cos(2*θ))*sin(θ)^4
sphharm75(θ::Float64, ϕ::Float64) = (-(3/64))*exp(5*ϕ*im)*sqrt(385/(2*pi))*(-1 + 13*cos(θ)^2)*sin(θ)^5
sphharm76(θ::Float64, ϕ::Float64) = (3/64)*exp(6*ϕ*im)*sqrt(5005/pi)*cos(θ)*sin(θ)^6
sphharm77(θ::Float64, ϕ::Float64) = (-(3/64))*exp(7*ϕ*im)*sqrt(715/(2*pi))*sin(θ)^7
sphharm8m8(θ::Float64, ϕ::Float64) = ((3/256)*sqrt(12155/(2*pi))*sin(θ)^8)/exp(8*ϕ*im)
sphharm8m7(θ::Float64, ϕ::Float64) = ((3/64)*sqrt(12155/(2*pi))*cos(θ)*sin(θ)^7)/exp(7*ϕ*im)
sphharm8m6(θ::Float64, ϕ::Float64) = ((1/256)*sqrt(7293/pi)*(13 + 15*cos(2*θ))*sin(θ)^6)/exp(6*ϕ*im)
sphharm8m5(θ::Float64, ϕ::Float64) = ((3/128)*sqrt(17017/(2*pi))*cos(θ)*(3 + 5*cos(2*θ))*sin(θ)^5)/exp(5*ϕ*im)
sphharm8m4(θ::Float64, ϕ::Float64) = ((3/128)*sqrt(1309/(2*pi))*(1 - 26*cos(θ)^2 +65*cos(θ)^4)*sin(θ)^4)/exp(4*ϕ*im)
sphharm8m3(θ::Float64, ϕ::Float64) = ((1/64)*sqrt(19635/(2*pi))*cos(θ)*(3 - 26*cos(θ)^2 + 39*cos(θ)^4)*sin(θ)^3)/exp(3*ϕ*im)
sphharm8m2(θ::Float64, ϕ::Float64) = ((3/128)*sqrt(595/pi)*(-1 + 33*cos(θ)^2 - 143*cos(θ)^4 + 143*cos(θ)^6)*sin(θ)^2)/exp(2*ϕ*im)
sphharm8m1(θ::Float64, ϕ::Float64) = (1/4096)*((3*sqrt(17/(2*pi))*(178 + 869*cos(2*θ) +286*cos(4*θ) + 715*cos(6*θ))*sin(2*θ))/exp(ϕ*im))
sphharm80(θ::Float64, ϕ::Float64) = complex((1/256)*sqrt(17/pi)*(35 - 1260*cos(θ)^2 + 6930*cos(θ)^4 -12012*cos(θ)^6 + 6435*cos(θ)^8))
sphharm81(θ::Float64, ϕ::Float64) = -((1/4096)*(3*exp(ϕ*im)*sqrt(17/(2*pi))*(178 + 869*cos(2*θ) + 286*cos(4*θ) +715*cos(6*θ))*sin(2*θ)))
sphharm82(θ::Float64, ϕ::Float64) = (3/128)*exp(2*ϕ*im)*sqrt(595/pi)*(-1 + 33*cos(θ)^2 -143*cos(θ)^4 + 143*cos(θ)^6)*sin(θ)^2
sphharm83(θ::Float64, ϕ::Float64) = (-(1/64))*exp(3*ϕ*im)*sqrt(19635/(2*pi))*cos(θ)*(3 - 26*cos(θ)^2 + 39*cos(θ)^4)*sin(θ)^3
sphharm84(θ::Float64, ϕ::Float64) = (3/128)*exp(4*ϕ*im)*sqrt(1309/(2*pi))*(1 - 26*cos(θ)^2 + 65*cos(θ)^4)*sin(θ)^4
sphharm85(θ::Float64, ϕ::Float64) = (-(3/128))*exp(5*ϕ*im)*sqrt(17017/(2*pi))*cos(θ)*(3 + 5*cos(2*θ))*sin(θ)^5
sphharm86(θ::Float64, ϕ::Float64) = (1/256)*exp(6*ϕ*im)*sqrt(7293/pi)*(13 + 15*cos(2*θ))*sin(θ)^6
sphharm87(θ::Float64, ϕ::Float64) = (-(3/64))*exp(7*ϕ*im)*sqrt(12155/(2*pi))*cos(θ)*sin(θ)^7
sphharm88(θ::Float64, ϕ::Float64) = (3/256)*exp(8*ϕ*im)*sqrt(12155/(2*pi))*sin(θ)^8





realsphharm00(θ::Float64, ϕ::Float64) = 1/(2*sqrt(pi)) 
realsphharm1m1(θ::Float64, ϕ::Float64) = (1/2)*sqrt(3/pi)*sin(ϕ)*sin(θ) 
realsphharm10(θ::Float64, ϕ::Float64) = (1/2)*sqrt(3/pi)*cos(θ) 
realsphharm11(θ::Float64, ϕ::Float64) = (1/2)*sqrt(3/pi)*cos(ϕ)*sin(θ) 
realsphharm2m2(θ::Float64, ϕ::Float64) = (-(1/4))*sqrt(15/pi)*sin(2*ϕ)*sin(θ)^2 
realsphharm2m1(θ::Float64, ϕ::Float64) = (1/2)*sqrt(15/pi)*cos(θ)*sin(ϕ)*sin(θ) 
realsphharm20(θ::Float64, ϕ::Float64) = (1/8)*sqrt(5/pi)*(1 + 3*cos(2*θ)) 
realsphharm21(θ::Float64, ϕ::Float64) = (1/2)*sqrt(15/pi)*cos(ϕ)*cos(θ)*sin(θ) 
realsphharm22(θ::Float64, ϕ::Float64) = (1/4)*sqrt(15/pi)*cos(2*ϕ)*sin(θ)^2 
realsphharm3m3(θ::Float64, ϕ::Float64) = (1/4)*sqrt(35/(2*pi))*sin(3*ϕ)*sin(θ)^3 
realsphharm3m2(θ::Float64, ϕ::Float64) = (-(1/4))*sqrt(105/pi)*cos(θ)*sin(2*ϕ)*sin(θ)^2 
realsphharm3m1(θ::Float64, ϕ::Float64) = (1/16)*sqrt(21/(2*pi))*sin(ϕ)*(sin(θ) + 5*sin(3*θ)) 
realsphharm30(θ::Float64, ϕ::Float64) = (1/8)*sqrt(7/pi)*cos(θ)*(-1 + 5*cos(2*θ)) 
realsphharm31(θ::Float64, ϕ::Float64) = (1/16)*sqrt(21/(2*pi))*cos(ϕ)*(sin(θ) + 5*sin(3*θ)) 
realsphharm32(θ::Float64, ϕ::Float64) = (1/4)*sqrt(105/pi)*cos(2*ϕ)*cos(θ)*sin(θ)^2 
realsphharm33(θ::Float64, ϕ::Float64) = (1/4)*sqrt(35/(2*pi))*cos(3*ϕ)*sin(θ)^3 
realsphharm4m4(θ::Float64, ϕ::Float64) = (-(3/16))*sqrt(35/pi)*sin(4*ϕ)*sin(θ)^4 
realsphharm4m3(θ::Float64, ϕ::Float64) = (3/4)*sqrt(35/(2*pi))*cos(θ)*sin(3*ϕ)*sin(θ)^3 
realsphharm4m2(θ::Float64, ϕ::Float64) = (3/8)*sqrt(5/pi)*(1 - 7*cos(θ)^2)*sin(2*ϕ)*sin(θ)^2 
realsphharm4m1(θ::Float64, ϕ::Float64) = (3/32)*sqrt(5/(2*pi))*sin(ϕ)*(2*sin(2*θ) + 7*sin(4*θ)) 
realsphharm40(θ::Float64, ϕ::Float64) = (3*(9 + 20*cos(2*θ) + 35*cos(4*θ)))/(128*sqrt(pi)) 
realsphharm41(θ::Float64, ϕ::Float64) = (3/32)*sqrt(5/(2*pi))*cos(ϕ)*(2*sin(2*θ) + 7*sin(4*θ)) 
realsphharm42(θ::Float64, ϕ::Float64) = (3/16)*sqrt(5/pi)*cos(2*ϕ)*(5 + 7*cos(2*θ))*sin(θ)^2 
realsphharm43(θ::Float64, ϕ::Float64) = (3/4)*sqrt(35/(2*pi))*cos(3*ϕ)*cos(θ)*sin(θ)^3 
realsphharm44(θ::Float64, ϕ::Float64) = (3/16)*sqrt(35/pi)*cos(4*ϕ)*sin(θ)^4 
realsphharm5m5(θ::Float64, ϕ::Float64) = (3/16)*sqrt(77/(2*pi))*sin(5*ϕ)*sin(θ)^5 
realsphharm5m4(θ::Float64, ϕ::Float64) = (-(3/16))*sqrt(385/pi)*cos(θ)*sin(4*ϕ)*sin(θ)^4 
realsphharm5m3(θ::Float64, ϕ::Float64) = (1/32)*sqrt(385/(2*pi))*(7 + 9*cos(2*θ))*sin(3*ϕ)*sin(θ)^3 
realsphharm5m2(θ::Float64, ϕ::Float64) = (-(1/32))*sqrt(1155/pi)*(5*cos(θ) + 3*cos(3*θ))*sin(2*ϕ)*sin(θ)^2 
realsphharm5m1(θ::Float64, ϕ::Float64) = (1/256)*sqrt(165/pi)*sin(ϕ)*(2*sin(θ) + 7*(sin(3*θ) + 3*sin(5*θ))) 
realsphharm50(θ::Float64, ϕ::Float64) = (1/256)*sqrt(11/pi)*(30*cos(θ) + 35*cos(3*θ) + 63*cos(5*θ)) 
realsphharm51(θ::Float64, ϕ::Float64) = (1/256)*sqrt(165/pi)*cos(ϕ)*(2*sin(θ) + 7*(sin(3*θ) + 3*sin(5*θ))) 
realsphharm52(θ::Float64, ϕ::Float64) = (1/16)*sqrt(1155/pi)*cos(2*ϕ)*cos(θ)*(1 + 3*cos(2*θ))*sin(θ)^2 
realsphharm53(θ::Float64, ϕ::Float64) = (1/32)*sqrt(385/(2*pi))*cos(3*ϕ)*(7 + 9*cos(2*θ))*sin(θ)^3 
realsphharm54(θ::Float64, ϕ::Float64) = (3/16)*sqrt(385/pi)*cos(4*ϕ)*cos(θ)*sin(θ)^4 
realsphharm55(θ::Float64, ϕ::Float64) = (3/16)*sqrt(77/(2*pi))*cos(5*ϕ)*sin(θ)^5 
realsphharm6m6(θ::Float64, ϕ::Float64) = (-(1/32))*sqrt(3003/(2*pi))*sin(6*ϕ)*sin(θ)^6 
realsphharm6m5(θ::Float64, ϕ::Float64) = (3/16)*sqrt(1001/(2*pi))*cos(θ)*sin(5*ϕ)*sin(θ)^5 
realsphharm6m4(θ::Float64, ϕ::Float64) = (-(3/64))*sqrt(91/pi)*(9 + 11*cos(2*θ))*sin(4*ϕ)*sin(θ)^4 
realsphharm6m3(θ::Float64, ϕ::Float64) = (1/64)*sqrt(1365/(2*pi))*(21*cos(θ) + 11*cos(3*θ))*sin(3*ϕ)*sin(θ)^3 
realsphharm6m2(θ::Float64, ϕ::Float64) = (-(1/256))*sqrt(1365/(2*pi))*(35 + 60*cos(2*θ) + 33*cos(4*θ))*sin(2*ϕ)*sin(θ)^2 
realsphharm6m1(θ::Float64, ϕ::Float64) = (1/512)*sqrt(273/pi)*sin(ϕ)*(5*sin(2*θ) + 12*sin(4*θ) + 33*sin(6*θ)) 
realsphharm60(θ::Float64, ϕ::Float64) = (sqrt(13/pi)*(50 + 105*cos(2*θ) + 126*cos(4*θ) + 231*cos(6*θ)))/1024 
realsphharm61(θ::Float64, ϕ::Float64) = (1/512)*sqrt(273/pi)*cos(ϕ)*(5*sin(2*θ) + 12*sin(4*θ) + 33*sin(6*θ)) 
realsphharm62(θ::Float64, ϕ::Float64) = (1/256)*sqrt(1365/(2*pi))*cos(2*ϕ)*(35 + 60*cos(2*θ) + 33*cos(4*θ))*sin(θ)^2 
realsphharm63(θ::Float64, ϕ::Float64) = (1/64)*sqrt(1365/(2*pi))*cos(3*ϕ)*(21*cos(θ) + 11*cos(3*θ))*sin(θ)^3 
realsphharm64(θ::Float64, ϕ::Float64) = (3/64)*sqrt(91/pi)*cos(4*ϕ)*(9 + 11*cos(2*θ))*sin(θ)^4 
realsphharm65(θ::Float64, ϕ::Float64) = (3/16)*sqrt(1001/(2*pi))*cos(5*ϕ)*cos(θ)*sin(θ)^5 
realsphharm66(θ::Float64, ϕ::Float64) = (1/32)*sqrt(3003/(2*pi))*cos(6*ϕ)*sin(θ)^6 
realsphharm7m7(θ::Float64, ϕ::Float64) = (3/64)*sqrt(715/pi)*sin(7*ϕ)*sin(θ)^7 
realsphharm7m6(θ::Float64, ϕ::Float64) = (-(3/32))*sqrt(5005/(2*pi))*cos(θ)*sin(6*ϕ)*sin(θ)^6 
realsphharm7m5(θ::Float64, ϕ::Float64) = (3/128)*sqrt(385/pi)*(11 + 13*cos(2*θ))*sin(5*ϕ)*sin(θ)^5 
realsphharm7m4(θ::Float64, ϕ::Float64) = (-(3/64))*sqrt(385/pi)*cos(θ)*(7 + 13*cos(2*θ))*sin(4*ϕ)*sin(θ)^4 
realsphharm7m3(θ::Float64, ϕ::Float64) = (3/512)*sqrt(35/pi)*(189 + 308*cos(2*θ) + 143*cos(4*θ))*sin(3*ϕ)*sin(θ)^3 
realsphharm7m2(θ::Float64, ϕ::Float64) = (-(3/512))*sqrt(35/(2*pi))*(350*cos(θ) + 275*cos(3*θ) + 143*cos(5*θ))*sin(2*ϕ)*sin(θ)^2 
realsphharm7m1(θ::Float64, ϕ::Float64) = (sqrt(105/pi)*sin(ϕ)*(25*sin(θ) + 81*sin(3*θ) + 165*sin(5*θ) + 429*sin(7*θ)))/4096 
realsphharm70(θ::Float64, ϕ::Float64) = (sqrt(15/pi)*(175*cos(θ) + 189*cos(3*θ) + 231*cos(5*θ) + 429*cos(7*θ)))/2048 
realsphharm71(θ::Float64, ϕ::Float64) = (sqrt(105/pi)*cos(ϕ)*(25*sin(θ) + 81*sin(3*θ) + 165*sin(5*θ) + 429*sin(7*θ)))/4096 
realsphharm72(θ::Float64, ϕ::Float64) = (3/512)*sqrt(35/(2*pi))*cos(2*ϕ)*(350*cos(θ) + 275*cos(3*θ) + 143*cos(5*θ))*sin(θ)^2 
realsphharm73(θ::Float64, ϕ::Float64) = (3/512)*sqrt(35/pi)*cos(3*ϕ)*(189 + 308*cos(2*θ) + 143*cos(4*θ))*sin(θ)^3 
realsphharm74(θ::Float64, ϕ::Float64) = (3/128)*sqrt(385/pi)*cos(4*ϕ)*(27*cos(θ) + 13*cos(3*θ))*sin(θ)^4 
realsphharm75(θ::Float64, ϕ::Float64) = (3/128)*sqrt(385/pi)*cos(5*ϕ)*(11 + 13*cos(2*θ))*sin(θ)^5 
realsphharm76(θ::Float64, ϕ::Float64) = (3/32)*sqrt(5005/(2*pi))*cos(6*ϕ)*cos(θ)*sin(θ)^6 
realsphharm77(θ::Float64, ϕ::Float64) = (3/64)*sqrt(715/pi)*cos(7*ϕ)*sin(θ)^7 
realsphharm8m8(θ::Float64, ϕ::Float64) = (-(3/256))*sqrt(12155/pi)*sin(8*ϕ)*sin(θ)^8 
realsphharm8m7(θ::Float64, ϕ::Float64) = (3/64)*sqrt(12155/pi)*cos(θ)*sin(7*ϕ)*sin(θ)^7 
realsphharm8m6(θ::Float64, ϕ::Float64) = (-(1/128))*sqrt(7293/(2*pi))*(13 + 15*cos(2*θ))*sin(6*ϕ)*sin(θ)^6 
realsphharm8m5(θ::Float64, ϕ::Float64) = (3/256)*sqrt(17017/pi)*(11*cos(θ) + 5*cos(3*θ))*sin(5*ϕ)*sin(θ)^5 
realsphharm8m4(θ::Float64, ϕ::Float64) = -((3*sqrt(1309/pi)*(99 + 156*cos(2*θ) + 65*cos(4*θ))*sin(4*ϕ)*sin(θ)^4)/1024) 
realsphharm8m3(θ::Float64, ϕ::Float64) = (sqrt(19635/pi)*(126*cos(θ) + 91*cos(3*θ) + 39*cos(5*θ))*sin(3*ϕ)*sin(θ)^3)/1024 
realsphharm8m2(θ::Float64, ϕ::Float64) = -((3*sqrt(595/(2*pi))*(210 + 385*cos(2*θ) + 286*cos(4*θ) + 143*cos(6*θ))*sin(2*ϕ)*sin(θ)^2)/2048) 
realsphharm8m1(θ::Float64, ϕ::Float64) = (3*sqrt(17/pi)*sin(ϕ)*(70*sin(2*θ) + 154*sin(4*θ) + 286*sin(6*θ) + 715*sin(8*θ)))/8192 
realsphharm80(θ::Float64, ϕ::Float64) = (sqrt(17/pi)*(1225 + 2520*cos(2*θ) + 2772*cos(4*θ) + 3432*cos(6*θ) + 6435*cos(8*θ)))/32768 
realsphharm81(θ::Float64, ϕ::Float64) = (3*sqrt(17/pi)*cos(ϕ)*(70*sin(2*θ) + 154*sin(4*θ) + 286*sin(6*θ) + 715*sin(8*θ)))/8192 
realsphharm82(θ::Float64, ϕ::Float64) = (3*sqrt(595/(2*pi))*cos(2*ϕ)*(210 + 385*cos(2*θ) + 286*cos(4*θ) + 143*cos(6*θ))*sin(θ)^2)/2048 
realsphharm83(θ::Float64, ϕ::Float64) = (sqrt(19635/pi)*cos(3*ϕ)*(126*cos(θ) + 91*cos(3*θ) + 39*cos(5*θ))*sin(θ)^3)/1024 
realsphharm84(θ::Float64, ϕ::Float64) = (3*sqrt(1309/pi)*cos(4*ϕ)*(99 + 156*cos(2*θ) + 65*cos(4*θ))*sin(θ)^4)/1024 
realsphharm85(θ::Float64, ϕ::Float64) = (3/256)*sqrt(17017/pi)*cos(5*ϕ)*(11*cos(θ) + 5*cos(3*θ))*sin(θ)^5 
realsphharm86(θ::Float64, ϕ::Float64) = (1/128)*sqrt(7293/(2*pi))*cos(6*ϕ)*(13 + 15*cos(2*θ))*sin(θ)^6 
realsphharm87(θ::Float64, ϕ::Float64) = (3/64)*sqrt(12155/pi)*cos(7*ϕ)*cos(θ)*sin(θ)^7 
realsphharm88(θ::Float64, ϕ::Float64) = (3/256)*sqrt(12155/pi)*cos(8*ϕ)*sin(θ)^8




i_lm(l::Int, m::Int) = l*(l+1)+1+m


sphharmtab = [
    sphharm00,
    sphharm1m1,
    sphharm10,
    sphharm11,
    sphharm2m2,
    sphharm2m1,
    sphharm20,
    sphharm21,
    sphharm22,
    sphharm3m3,
    sphharm3m2,
    sphharm3m1,
    sphharm30,
    sphharm31,
    sphharm32,
    sphharm33,
    sphharm4m4,
    sphharm4m3,
    sphharm4m2,
    sphharm4m1,
    sphharm40,
    sphharm41,
    sphharm42,
    sphharm43,
    sphharm44,
    sphharm5m5,
    sphharm5m4,
    sphharm5m3,
    sphharm5m2,
    sphharm5m1,
    sphharm50,
    sphharm51,
    sphharm52,
    sphharm53,
    sphharm54,
    sphharm55,
    sphharm6m6,
    sphharm6m5,
    sphharm6m4,
    sphharm6m3,
    sphharm6m2,
    sphharm6m1,
    sphharm60,
    sphharm61,
    sphharm62,
    sphharm63,
    sphharm64,
    sphharm65,
    sphharm66,
    sphharm7m7,
    sphharm7m6,
    sphharm7m5,
    sphharm7m4,
    sphharm7m3,
    sphharm7m2,
    sphharm7m1,
    sphharm70,
    sphharm71,
    sphharm72,
    sphharm73,
    sphharm74,
    sphharm75,
    sphharm76,
    sphharm77,
    sphharm8m8,
    sphharm8m7,
    sphharm8m6,
    sphharm8m5,
    sphharm8m4,
    sphharm8m3,
    sphharm8m2,
    sphharm8m1,
    sphharm80,
    sphharm81,
    sphharm82,
    sphharm83,
    sphharm84,
    sphharm85,
    sphharm86,
    sphharm87,
    sphharm88
]


realsphharmtab = [
    realsphharm00,
    realsphharm1m1,
    realsphharm10,
    realsphharm11,
    realsphharm2m2,
    realsphharm2m1,
    realsphharm20,
    realsphharm21,
    realsphharm22,
    realsphharm3m3,
    realsphharm3m2,
    realsphharm3m1,
    realsphharm30,
    realsphharm31,
    realsphharm32,
    realsphharm33,
    realsphharm4m4,
    realsphharm4m3,
    realsphharm4m2,
    realsphharm4m1,
    realsphharm40,
    realsphharm41,
    realsphharm42,
    realsphharm43,
    realsphharm44,
    realsphharm5m5,
    realsphharm5m4,
    realsphharm5m3,
    realsphharm5m2,
    realsphharm5m1,
    realsphharm50,
    realsphharm51,
    realsphharm52,
    realsphharm53,
    realsphharm54,
    realsphharm55,
    realsphharm6m6,
    realsphharm6m5,
    realsphharm6m4,
    realsphharm6m3,
    realsphharm6m2,
    realsphharm6m1,
    realsphharm60,
    realsphharm61,
    realsphharm62,
    realsphharm63,
    realsphharm64,
    realsphharm65,
    realsphharm66,
    realsphharm7m7,
    realsphharm7m6,
    realsphharm7m5,
    realsphharm7m4,
    realsphharm7m3,
    realsphharm7m2,
    realsphharm7m1,
    realsphharm70,
    realsphharm71,
    realsphharm72,
    realsphharm73,
    realsphharm74,
    realsphharm75,
    realsphharm76,
    realsphharm77,
    realsphharm8m8,
    realsphharm8m7,
    realsphharm8m6,
    realsphharm8m5,
    realsphharm8m4,
    realsphharm8m3,
    realsphharm8m2,
    realsphharm8m1,
    realsphharm80,
    realsphharm81,
    realsphharm82,
    realsphharm83,
    realsphharm84,
    realsphharm85,
    realsphharm86,
    realsphharm87,
    realsphharm88
]


