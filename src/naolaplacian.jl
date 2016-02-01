
naolaplacian00(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (2*ud + r*ud2)/(2*sqrt(π)*r)
naolaplacian1m1(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(3/π)*sin(θ)*sin(ϕ)*(-2*u + r*(2*ud + r*ud2)))/(2*r^2)
naolaplacian10(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(3/π)*cos(θ)*(-2*u + r*(2*ud + r*ud2)))/(2*r^2)
naolaplacian11(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(3/π)*cos(ϕ)*sin(θ)*(-2*u + r*(2*ud + r*ud2)))/(2*r^2)
naolaplacian2m2(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((sqrt(15/π)*sin(θ)^2*sin(2*ϕ)*(-6*u + r*(2*ud + r*ud2)))/(4*r^2))
naolaplacian2m1(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(15/π)*sin(2*θ)*sin(ϕ)*(-6*u + r*(2*ud + r*ud2)))/(4*r^2)
naolaplacian20(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(5/π)*(1 + 3*cos(2*θ))*(-6*u + r*(2*ud + r*ud2)))/(8*r^2)
naolaplacian21(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(15/π)*cos(ϕ)*sin(2*θ)*(-6*u + r*(2*ud + r*ud2)))/(4*r^2)
naolaplacian22(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(15/π)*cos(2*ϕ)*sin(θ)^2*(-6*u + r*(2*ud + r*ud2)))/(4*r^2)
naolaplacian3m3(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(35/(2*π))*sin(θ)^3*sin(3*ϕ)*(-12*u + r*(2*ud + r*ud2)))/(4*r^2)
naolaplacian3m2(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((sqrt(105/π)*cos(θ)*sin(θ)^2*sin(2*ϕ)*(-12*u + r*(2*ud + r*ud2)))/(4*r^2))
naolaplacian3m1(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(21/(2*π))*(3 + 5*cos(2*θ))*sin(θ)*sin(ϕ)*(-12*u + r*(2*ud + r*ud2)))/(8*r^2)
naolaplacian30(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(7/π)*(3*cos(θ) + 5*cos(3*θ))*(-12*u + r*(2*ud + r*ud2)))/(16*r^2)
naolaplacian31(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(21/(2*π))*(3 + 5*cos(2*θ))*cos(ϕ)*sin(θ)*(-12*u + r*(2*ud + r*ud2)))/(8*r^2)
naolaplacian32(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(105/π)*cos(θ)*cos(2*ϕ)*sin(θ)^2*(-12*u + r*(2*ud + r*ud2)))/(4*r^2)
naolaplacian33(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(35/(2*π))*cos(3*ϕ)*sin(θ)^3*(-12*u + r*(2*ud + r*ud2)))/(4*r^2)
naolaplacian4m4(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((3*sqrt(35/π)*sin(θ)^4*sin(4*ϕ)*(-20*u + r*(2*ud + r*ud2)))/(16*r^2))
naolaplacian4m3(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(35/(2*π))*cos(θ)*sin(θ)^3*sin(3*ϕ)*(-20*u + r*(2*ud + r*ud2)))/(4*r^2)
naolaplacian4m2(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((3*sqrt(5/π)*(5 + 7*cos(2*θ))*sin(θ)^2*sin(2*ϕ)*(-20*u + r*(2*ud + r*ud2)))/(16*r^2))
naolaplacian4m1(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(5/(2*π))*(2*sin(2*θ) + 7*sin(4*θ))*sin(ϕ)*(-20*u + r*(2*ud + r*ud2)))/(32*r^2)
naolaplacian40(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*(9 + 20*cos(2*θ) + 35*cos(4*θ))*(-20*u + r*(2*ud + r*ud2)))/(128*sqrt(π)*r^2)
naolaplacian41(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(5/(2*π))*cos(ϕ)*(2*sin(2*θ) + 7*sin(4*θ))*(-20*u + r*(2*ud + r*ud2)))/(32*r^2)
naolaplacian42(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(5/π)*(5 + 7*cos(2*θ))*cos(2*ϕ)*sin(θ)^2*(-20*u + r*(2*ud + r*ud2)))/(16*r^2)
naolaplacian43(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(35/(2*π))*cos(θ)*cos(3*ϕ)*sin(θ)^3*(-20*u + r*(2*ud + r*ud2)))/(4*r^2)
naolaplacian44(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(35/π)*cos(4*ϕ)*sin(θ)^4*(-20*u + r*(2*ud + r*ud2)))/(16*r^2)
naolaplacian5m5(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(77/(2*π))*sin(θ)^5*sin(5*ϕ)*(-30*u + r*(2*ud + r*ud2)))/(16*r^2)
naolaplacian5m4(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((3*sqrt(385/π)*cos(θ)*sin(θ)^4*sin(4*ϕ)*(-30*u + r*(2*ud + r*ud2)))/(16*r^2))
naolaplacian5m3(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(385/(2*π))*(7 + 9*cos(2*θ))*sin(θ)^3*sin(3*ϕ)*(-30*u + r*(2*ud + r*ud2)))/(32*r^2)
naolaplacian5m2(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((sqrt(1155/π)*(5*cos(θ) + 3*cos(3*θ))*sin(θ)^2*sin(2*ϕ)*(-30*u + r*(2*ud + r*ud2)))/(32*r^2))
naolaplacian5m1(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(165/π)*(15 + 28*cos(2*θ) + 21*cos(4*θ))*sin(θ)*sin(ϕ)*(-30*u + r*(2*ud + r*ud2)))/(128*r^2)
naolaplacian50(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(11/π)*(30*cos(θ) + 35*cos(3*θ) + 63*cos(5*θ))*(-30*u + r*(2*ud + r*ud2)))/(256*r^2)
naolaplacian51(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(165/π)*(15 + 28*cos(2*θ) + 21*cos(4*θ))*cos(ϕ)*sin(θ)*(-30*u + r*(2*ud + r*ud2)))/(128*r^2)
naolaplacian52(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(1155/π)*(5*cos(θ) + 3*cos(3*θ))*cos(2*ϕ)*sin(θ)^2*(-30*u + r*(2*ud + r*ud2)))/(32*r^2)
naolaplacian53(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(385/(2*π))*(7 + 9*cos(2*θ))*cos(3*ϕ)*sin(θ)^3*(-30*u + r*(2*ud + r*ud2)))/(32*r^2)
naolaplacian54(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(385/π)*cos(θ)*cos(4*ϕ)*sin(θ)^4*(-30*u + r*(2*ud + r*ud2)))/(16*r^2)
naolaplacian55(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(77/(2*π))*cos(5*ϕ)*sin(θ)^5*(-30*u + r*(2*ud + r*ud2)))/(16*r^2)
naolaplacian6m6(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((sqrt(3003/(2*π))*sin(θ)^6*sin(6*ϕ)*(-42*u + r*(2*ud + r*ud2)))/(32*r^2))
naolaplacian6m5(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(1001/(2*π))*cos(θ)*sin(θ)^5*sin(5*ϕ)*(-42*u + r*(2*ud + r*ud2)))/(16*r^2)
naolaplacian6m4(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((3*sqrt(91/π)*(9 + 11*cos(2*θ))*sin(θ)^4*sin(4*ϕ)*(-42*u + r*(2*ud + r*ud2)))/(64*r^2))
naolaplacian6m3(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(1365/(2*π))*(21*cos(θ) + 11*cos(3*θ))*sin(θ)^3*sin(3*ϕ)*(-42*u + r*(2*ud + r*ud2)))/(64*r^2)
naolaplacian6m2(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((sqrt(1365/(2*π))*(35 + 60*cos(2*θ) + 33*cos(4*θ))*sin(θ)^2*sin(2*ϕ)*(-42*u + r*(2*ud + r*ud2)))/(256*r^2))
naolaplacian6m1(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(273/π)*(50*cos(θ) + 45*cos(3*θ) + 33*cos(5*θ))*sin(θ)*sin(ϕ)*(-42*u + r*(2*ud + r*ud2)))/(256*r^2)
naolaplacian60(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(13/π)*(50 + 105*cos(2*θ) + 126*cos(4*θ) + 231*cos(6*θ))*(-42*u + r*(2*ud + r*ud2)))/(1024*r^2)
naolaplacian61(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(273/π)*(50*cos(θ) + 45*cos(3*θ) + 33*cos(5*θ))*cos(ϕ)*sin(θ)*(-42*u + r*(2*ud + r*ud2)))/(256*r^2)
naolaplacian62(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(1365/(2*π))*(35 + 60*cos(2*θ) + 33*cos(4*θ))*cos(2*ϕ)*sin(θ)^2*(-42*u + r*(2*ud + r*ud2)))/(256*r^2)
naolaplacian63(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(1365/(2*π))*(21*cos(θ) + 11*cos(3*θ))*cos(3*ϕ)*sin(θ)^3*(-42*u + r*(2*ud + r*ud2)))/(64*r^2)
naolaplacian64(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(91/π)*(9 + 11*cos(2*θ))*cos(4*ϕ)*sin(θ)^4*(-42*u + r*(2*ud + r*ud2)))/(64*r^2)
naolaplacian65(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(1001/(2*π))*cos(θ)*cos(5*ϕ)*sin(θ)^5*(-42*u + r*(2*ud + r*ud2)))/(16*r^2)
naolaplacian66(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(3003/(2*π))*cos(6*ϕ)*sin(θ)^6*(-42*u + r*(2*ud + r*ud2)))/(32*r^2)
naolaplacian7m7(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(715/π)*sin(θ)^7*sin(7*ϕ)*(-56*u + r*(2*ud + r*ud2)))/(64*r^2)
naolaplacian7m6(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((3*sqrt(5005/(2*π))*cos(θ)*sin(θ)^6*sin(6*ϕ)*(-56*u + r*(2*ud + r*ud2)))/(32*r^2))
naolaplacian7m5(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(385/π)*(11 + 13*cos(2*θ))*sin(θ)^5*sin(5*ϕ)*(-56*u + r*(2*ud + r*ud2)))/(128*r^2)
naolaplacian7m4(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((3*sqrt(385/π)*(27*cos(θ) + 13*cos(3*θ))*sin(θ)^4*sin(4*ϕ)*(-56*u + r*(2*ud + r*ud2)))/(128*r^2))
naolaplacian7m3(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(35/π)*(189 + 308*cos(2*θ) + 143*cos(4*θ))*sin(θ)^3*sin(3*ϕ)*(-56*u + r*(2*ud + r*ud2)))/(512*r^2)
naolaplacian7m2(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((3*sqrt(35/(2*π))*(350*cos(θ) + 275*cos(3*θ) + 143*cos(5*θ))*sin(θ)^2*sin(2*ϕ)*(-56*u + r*(2*ud + r*ud2)))/(512*r^2))
naolaplacian7m1(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(105/π)*(350 + 675*cos(2*θ) + 594*cos(4*θ) + 429*cos(6*θ))*sin(θ)*sin(ϕ)*(-56*u + r*(2*ud + r*ud2)))/(2048*r^2)
naolaplacian70(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(15/π)*(175*cos(θ) + 189*cos(3*θ) + 231*cos(5*θ) + 429*cos(7*θ))*(-56*u + r*(2*ud + r*ud2)))/(2048*r^2)
naolaplacian71(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(105/π)*(350 + 675*cos(2*θ) + 594*cos(4*θ) + 429*cos(6*θ))*cos(ϕ)*sin(θ)*(-56*u + r*(2*ud + r*ud2)))/(2048*r^2)
naolaplacian72(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(35/(2*π))*(350*cos(θ) + 275*cos(3*θ) + 143*cos(5*θ))*cos(2*ϕ)*sin(θ)^2*(-56*u + r*(2*ud + r*ud2)))/(512*r^2)
naolaplacian73(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(35/π)*(189 + 308*cos(2*θ) + 143*cos(4*θ))*cos(3*ϕ)*sin(θ)^3*(-56*u + r*(2*ud + r*ud2)))/(512*r^2)
naolaplacian74(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(385/π)*(27*cos(θ) + 13*cos(3*θ))*cos(4*ϕ)*sin(θ)^4*(-56*u + r*(2*ud + r*ud2)))/(128*r^2)
naolaplacian75(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(385/π)*(11 + 13*cos(2*θ))*cos(5*ϕ)*sin(θ)^5*(-56*u + r*(2*ud + r*ud2)))/(128*r^2)
naolaplacian76(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(5005/(2*π))*cos(θ)*cos(6*ϕ)*sin(θ)^6*(-56*u + r*(2*ud + r*ud2)))/(32*r^2)
naolaplacian77(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(715/π)*cos(7*ϕ)*sin(θ)^7*(-56*u + r*(2*ud + r*ud2)))/(64*r^2)
naolaplacian8m8(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((3*sqrt(12155/π)*sin(θ)^8*sin(8*ϕ)*(-72*u + r*(2*ud + r*ud2)))/(256*r^2))
naolaplacian8m7(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(12155/π)*cos(θ)*sin(θ)^7*sin(7*ϕ)*(-72*u + r*(2*ud + r*ud2)))/(64*r^2)
naolaplacian8m6(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((sqrt(7293/(2*π))*(13 + 15*cos(2*θ))*sin(θ)^6*sin(6*ϕ)*(-72*u + r*(2*ud + r*ud2)))/(128*r^2))
naolaplacian8m5(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(17017/π)*(11*cos(θ) + 5*cos(3*θ))*sin(θ)^5*sin(5*ϕ)*(-72*u + r*(2*ud + r*ud2)))/(256*r^2)
naolaplacian8m4(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((3*sqrt(1309/π)*(99 + 156*cos(2*θ) + 65*cos(4*θ))*sin(θ)^4*sin(4*ϕ)*(-72*u + r*(2*ud + r*ud2)))/
      (1024*r^2))
naolaplacian8m3(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(19635/π)*(126*cos(θ) + 91*cos(3*θ) + 39*cos(5*θ))*sin(θ)^3*sin(3*ϕ)*(-72*u + r*(2*ud + r*ud2)))/(1024*r^2)
naolaplacian8m2(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = -((3*sqrt(595/(2*π))*(210 + 385*cos(2*θ) + 286*cos(4*θ) + 143*cos(6*θ))*sin(θ)^2*sin(2*ϕ)*(-72*u + r*(2*ud + r*ud2)))/(2048*r^2))
naolaplacian8m1(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(17/π)*(70*sin(2*θ) + 154*sin(4*θ) + 286*sin(6*θ) + 715*sin(8*θ))*sin(ϕ)*(-72*u + r*(2*ud + r*ud2)))/(8192*r^2)
naolaplacian80(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(17/π)*(1225 + 2520*cos(2*θ) + 2772*cos(4*θ) + 3432*cos(6*θ) + 6435*cos(8*θ))*(-72*u + r*(2*ud + r*ud2)))/(32768*r^2)
naolaplacian81(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(17/π)*cos(ϕ)*(70*sin(2*θ) + 154*sin(4*θ) + 286*sin(6*θ) + 715*sin(8*θ))*(-72*u + r*(2*ud + r*ud2)))/(8192*r^2)
naolaplacian82(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(595/(2*π))*(210 + 385*cos(2*θ) + 286*cos(4*θ) + 143*cos(6*θ))*cos(2*ϕ)*sin(θ)^2*(-72*u + r*(2*ud + r*ud2)))/(2048*r^2)
naolaplacian83(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(19635/π)*(126*cos(θ) + 91*cos(3*θ) + 39*cos(5*θ))*cos(3*ϕ)*sin(θ)^3*(-72*u + r*(2*ud + r*ud2)))/(1024*r^2)
naolaplacian84(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(1309/π)*(99 + 156*cos(2*θ) + 65*cos(4*θ))*cos(4*ϕ)*sin(θ)^4*(-72*u + r*(2*ud + r*ud2)))/(1024*r^2)
naolaplacian85(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(17017/π)*(11*cos(θ) + 5*cos(3*θ))*cos(5*ϕ)*sin(θ)^5*(-72*u + r*(2*ud + r*ud2)))/(256*r^2)
naolaplacian86(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (sqrt(7293/(2*π))*(13 + 15*cos(2*θ))*cos(6*ϕ)*sin(θ)^6*(-72*u + r*(2*ud + r*ud2)))/(128*r^2)
naolaplacian87(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(12155/π)*cos(θ)*cos(7*ϕ)*sin(θ)^7*(-72*u + r*(2*ud + r*ud2)))/(64*r^2)
naolaplacian88(u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = (3*sqrt(12155/π)*cos(8*ϕ)*sin(θ)^8*(-72*u + r*(2*ud + r*ud2)))/(256*r^2)


# naolaplacian(l::Int, m::Int, u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = naolaplaciantab[l,m](u,ud,ud2,r,θ,ϕ)

naolaplacian(l::Int, m::Int, u::Float64, ud::Float64, ud2::Float64, r::Float64, θ::Float64, ϕ::Float64) = 
    naolaplaciantab[i_lm(l,m)](u,ud,ud2,r,θ,ϕ)::Float64



naolaplaciantab = [
    naolaplacian00,
    naolaplacian1m1,
    naolaplacian10,
    naolaplacian11,
    naolaplacian2m2,
    naolaplacian2m1,
    naolaplacian20,
    naolaplacian21,
    naolaplacian22,
    naolaplacian3m3,
    naolaplacian3m2,
    naolaplacian3m1,
    naolaplacian30,
    naolaplacian31,
    naolaplacian32,
    naolaplacian33,
    naolaplacian4m4,
    naolaplacian4m3,
    naolaplacian4m2,
    naolaplacian4m1,
    naolaplacian40,
    naolaplacian41,
    naolaplacian42,
    naolaplacian43,
    naolaplacian44,
    naolaplacian5m5,
    naolaplacian5m4,
    naolaplacian5m3,
    naolaplacian5m2,
    naolaplacian5m1,
    naolaplacian50,
    naolaplacian51,
    naolaplacian52,
    naolaplacian53,
    naolaplacian54,
    naolaplacian55,
    naolaplacian6m6,
    naolaplacian6m5,
    naolaplacian6m4,
    naolaplacian6m3,
    naolaplacian6m2,
    naolaplacian6m1,
    naolaplacian60,
    naolaplacian61,
    naolaplacian62,
    naolaplacian63,
    naolaplacian64,
    naolaplacian65,
    naolaplacian66,
    naolaplacian7m7,
    naolaplacian7m6,
    naolaplacian7m5,
    naolaplacian7m4,
    naolaplacian7m3,
    naolaplacian7m2,
    naolaplacian7m1,
    naolaplacian70,
    naolaplacian71,
    naolaplacian72,
    naolaplacian73,
    naolaplacian74,
    naolaplacian75,
    naolaplacian76,
    naolaplacian77,
    naolaplacian8m8,
    naolaplacian8m7,
    naolaplacian8m6,
    naolaplacian8m5,
    naolaplacian8m4,
    naolaplacian8m3,
    naolaplacian8m2,
    naolaplacian8m1,
    naolaplacian80,
    naolaplacian81,
    naolaplacian82,
    naolaplacian83,
    naolaplacian84,
    naolaplacian85,
    naolaplacian86,
    naolaplacian87,
    naolaplacian88
]



