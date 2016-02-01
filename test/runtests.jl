using Element
using Base.Test


using Element: lingrid, integrate, findsteps, eig, loggrid
using Element: SphericalHarmonicTable, radschrod1, cart2sph, realsphharm, enumeratelm, randp, radpoisson




# findsteps

f = (x,c) -> round(Int,x^2+c)
intervals = findsteps(f,0.0,2.0,args=1.0)
@test intervals == [(0.0,1.0),(1.0,1.5),(1.5,1.75),(1.75,2.0)]





# eig

x = lingrid(0,1,1000)
dx = x[2]-x[1]
v = 0*x
aE = E -> 2*(E-v)
qmin, qmax = 0,10

evs,efs,q = eig(aE,dx,qmin,qmax,Emin=0.0,Emax=20.0)
n = length(evs)
efs ./= sqrt(integrate(x, efs.^2))'  # normalization

exact_evs = [j^2*pi^2/2 for j in 1:n]
exact_efs = hcat([sqrt(2)*sin(j*pi*x) for j in 1:length(evs)]...)

@test q == collect(qmin:qmax)
@test_approx_eq integrate(x,efs.^2) fill(1.0,n)
@test_approx_eq_eps integrate(x, efs.^2 - exact_efs.^2) fill(0.0,n) 1e-12
@test_approx_eq_eps evs exact_evs 1e-5



# radschrod1


r = loggrid(-22,4,1000)
v = -1./r

n,l = 1,0
ev,ef = radschrod1(r,v,n,l)
@test_approx_eq_eps ev -1/(2*n^2) 1e-8
@test_approx_eq_eps integrate(r, (ef.^2 - (2r.*exp(-1*r)).^2).^2) 0.0 1e-12

n,l = 2,0
ev,ef = radschrod1(r,v,n,l)
@test_approx_eq_eps ev -1/(2*n^2) 1e-8
@test_approx_eq_eps integrate(r, (ef.^2 - (-1*r/sqrt(2).*(1-0.5*r).*exp(-0.5*r)).^2).^2) 0.0 1e-12

n,l = 2,1
ev,ef = radschrod1(r,v,n,l)
@test_approx_eq_eps ev -1/(2*n^2) 1e-8
@test_approx_eq_eps integrate(r, (ef.^2 - (r/3^0.5/2^(3./2).*r.*exp(-0.5*r)).^2).^2) 0.0 1e-12

n,l = 3,1
ev,ef = radschrod1(r,v,n,l)
@test_approx_eq_eps ev -1/(2*n^2) 1e-8
@test_approx_eq_eps integrate(r, (ef.^2 - (r*8./27/sqrt(6).*(1-r/6).*r.*exp(-1.*r/3)).^2).^2) 0.0 1e-12






# radpoisson
r = loggrid(-13,5,500)
rho = exp(-2r)

l = 0
vh_exact = 1/4*exp(-2r).*(1 + 2r) + (1 + exp(-2r).*(-1 - 2r.*(1 + r)))./(4r)  # l = 0
vh = radpoisson(r,rho,l)

@test_approx_eq_eps vh vh_exact 1e-5
@test_approx_eq_eps integrate(r, rho.*vh_exact/2) integrate(r, rho.*vh/2) 1e-6

l = 1
vh_exact = (3*exp(-2r).*(-1 + exp(2r) - 2r - 2r.^2))./(8r.^2)  # l = 1
vh = radpoisson(r,rho,l)

@test_approx_eq_eps vh vh_exact 1e-5
@test_approx_eq_eps integrate(r, rho.*vh_exact/2) integrate(r, rho.*vh/2) 1e-6







# spherical harmonics table

p = randp(3)[]
p = p/norm(p)
r,θ,ϕ = cart2sph(p)

lmax = 8
Y = SphericalHarmonicTable(lmax)
Y(p)  # precompute spherical harmonics

a = [realsphharm(l,m,θ,ϕ) for (l,m) in enumeratelm(lmax)]
b = [Y[l,m] for (l,m) in enumeratelm(lmax)]

@test_approx_eq a b





