


doc"""
`radpoisson(r,rho,l)`

solves the poisson equation for a radial density `rho` in
the `l`-angular momentum channel

computes the hartree potential

$v_{h,lm}(r) = \int_0^r dr_<\, r_<^2 g_l(r_<,r) \rho_{lm}(r_<) + 
	\int_r^\infty dr_>\, r_>^2 g_l(r,r_>) \rho_{lm}(r_>)$,

where $g_l(r_<,r_>) = r_<^l/r_>^{l+1}$

uses a 3rd order adams-moulton multistep integrator to compute
each integral

`r` (LogarithmicGrid) - logarithmic grid (`loggrid`)

`rho` (AbstractVector) - radial density $\rho_{lm}(r)$ evaluated on `r` 

`l` (Int) - l-quantum number

e.g.
```julia
l = 0
r = loggrid(-13,5,100)
rho = exp(-2r)
vh = radpoisson(r,rho,l)
q = multipolemoment(r,rho,l)
plot(r,vh)
plot!(r,q./r.^(l+1))  # asymptotic form of hartree potential
xlims!(0,10)
ylims!(0,0.3)
```

[`V. Blum et al., CPC 180 (2009)` p. 2185-2186]
"""
function radpoisson(r::LogarithmicGrid, rho::AbstractVector, l::Int)
    a = similar(r)
    b = similar(r)
    poisson_multistep_a(a,r,rho,r.dx,l)
    poisson_multistep_b(b,reverse(r),reverse(rho),-r.dx,l)
    a + reverse(b)
end



@adamsmoulton3_multistep 	poisson_multistep_a 	-(l+1)*y[j] + r[j]^2*rho[j]
@adamsmoulton3_multistep 	poisson_multistep_b 	     l*y[j] - r[j]^2*rho[j]


doc"""
`multipolemoment(r,rho,l)`

computes the multipole moment of radial density `rho` in
the `l`-angular momentum channel

$q = \int dr\, r^{l+2} \rho_{lm}(r)$.

`r` (LogarithmicGrid) - logarithmic grid (`loggrid`)

`rho` (AbstractVector) - radial density $\rho_{lm}(r)$ evaluated on `r` 

`l` (Int) - l-quantum number

"""
multipolemoment(r,rho,l) = integrate(r,r.^(l+2).*rho)










# @fastmath function poisson_multistep_rule(f_jm1::Float64, f_j::Float64, f_jp1::Float64, y_jm1::Float64, y_j::Float64, h::Float64, c::Float64)
#     (h*h*(f_jp1 + 10f_j + f_jm1 + c*(10y_j + y_jm1)) + 24y_j - 12y_jm1)/(12 - c*h*h)
# end

# @fastmath @inbounds function poisson_multistep_integrator!(y::Vector{Float64}, f::Vector{Float64}, h::Float64, 
#                                                      y1::Float64, y2::Float64, c::Float64)
#     y[1] = y1
#     y[2] = y2
#     for j in 2:length(f)-1
#         y[j+1] = poisson_multistep_rule(f[j-1], f[j], f[j+1], y[j-1], y[j], h, c)
#     end
# end

# @fastmath @inbounds function poisson_multistep_integrator(f::Vector{Float64}, h::Float64, y1::Float64, y2::Float64, c::Float64)
#     y = similar(f)
#     poisson_multistep_integrator!(y,f,h,y1,y2,c)
#     y
# end

# @multistep poisson_multistep_integrator y[j+1] = (h*h*(f[j+1] + 10f[j] + f[j-1] + c*(10y[j] + y[j-1])) + 24y[j] - 12y[j-1])/(12 - c*h*h)


# @multistep numerov_multistep y[j+1] = ( (2-5h*h*g[j]/6)*y[j] - (1+h*h*g[j-1]/12)*y[j-1] + h*h/12*(s[j+1]+10s[j]+s[j-1]) ) / (1+h*h*g[j+1]/12)



# function radpoisson(r::LogarithmicGrid, rho::AbstractVector)
#     const c = 0.25
#     N = integrate(r, 4pi*r.^2.*rho)
#     f = -4pi.*r.^2.5.*rho
#     U = similar(r)
#     poisson_multistep_integrator!(U, 0., 1e-16, f, c, r.dx)
#     U .*= sqrt(r)
#     R = r[end]
#     U += (N - U[end])/R.*r
#     vh = U./r
#     vh/4pi
# end







# @fastmath function multistep_rule1(f_jm1::Float64, f_j::Float64, f_jp1::Float64, y_jm1::Float64, y_j::Float64, h::Float64)
#     2y_j - y_jm1 + h*h*(f_jp1 + 10f_j + f_jm1)/12
# end

# @fastmath @inbounds function multistep_integrator1!(y::Vector{Float64}, f::Vector{Float64},
#                                                       h::Float64, y1::Float64, y2::Float64)
#     y[1] = y1
#     y[2] = y2
#     for j in 2:length(f)-1
#         y[j+1] = multistep_rule1(f[j-1], f[j], f[j+1], y[j-1], y[j], h)
#     end
# end

# @fastmath @inbounds function multistep_integrator1(f::Vector{Float64}, h::Float64, y1::Float64, y2::Float64)
#     y = similar(f)
#     multistep_integrator1!(y,f,h,y1,y2)
#     y
# end



# uniform radial grid
# function poissonrad(rho,r)
#     N = integrate(r, 4pi*r.^2.*rho)
#     dr = r[1] - r[0]
#     f = -4pi.*r.*rho
#     U = similar(r)
#     multistep_integrator1!(U, f, dr, 0., 1e-16)
#     R = r[end]
#     U += (N - U[end])/R.*r
#     vh = U./r
#     vh
# end


# function hartreeatomicrad(r::Grid, rho::AbstractArray, l::Int)
#     cumsum(r.w.*r.^(2+l).*rho)./r.^(l+1) + r.^l.*flipdim(cumsum(flipdim((r.w.*r.^(2-l-1).*rho),1)),1)
# end

