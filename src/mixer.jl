
abstract AbstractMixer

type LinearMixer <: AbstractMixer
    c::Float64  # mixing weight
end

LinearMixer(;c=0.25) = LinearMixer(c)

Base.call(mixer::LinearMixer, rho_in::AbstractArray, rho_out::AbstractArray) = (1 - mixer.c)*rho_in + mixer.c*rho_out



type PulayMixer <: AbstractMixer
    n::Int  # number of densities to mix
    beta::Float64  # damping coefficient
    grid  # integration grid
    rho_in::Array{AbstractArray}
    R::Array{AbstractArray}  # residuals
end

PulayMixer(;n=3, beta=0.25, grid=nothing) = PulayMixer(n,beta,grid,[],[])





function Base.call(mixer::PulayMixer, rho_in::AbstractArray, rho_out::AbstractArray)
    push!(mixer.rho_in, rho_in)
    push!(mixer.R, rho_out - rho_in)
    
    if length(mixer.rho_in) > mixer.n
        shift!(mixer.rho_in)
        shift!(mixer.R)
    end

    n = length(mixer.rho_in)
    
    B = Array{Float64}(n+1,n+1)
    for i in 1:n
        for j in i:n
            B[i,j] = integrate(mixer.grid, mixer.R[i].*mixer.R[j])
        end
    end
    
    B[:,end] = -1
    B[end,end] = 0
    B = Symmetric(B)
    
    c = zeros(1,n+1)
    c[end] = -1
    
    a = squeeze(c/B,1)  # solve
    
    rho = zeros(rho_in)
    for i in 1:n
        rho += a[i].*(mixer.rho_in[i] + mixer.beta.*mixer.R[i])
    end
    
    rho
end