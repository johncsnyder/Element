




type RandomFourierKernelMachine
    d::Int
    dout
    σ::Float64
    w::Array{Float64,2}
    α::Union{Float64,Array{Float64}}
    z0::Union{Float64,Array{Float64}}
    y0::Union{Float64,Array{Float64}}
end


using Distributions: Normal

function fit(X::AbstractArray,Y::AbstractArray;λ::Float64=1e-6,d::Int=100,σ::Float64=1.0,centered=true)
    # x => n x k, y => m x k
    # n = number of input dimensions
    # m = number of output dimensions
    # d = number of random features
    # k = number of samples
    x = mat(X)
    y = mat(Y)
    dout = size(Y)[1:end-1]
    dist = Normal(0,1/σ)
    w = rand(dist, (d,size(x,1)))
    z = real(exp(w*x*im)/sqrt(d))
    z0 = squeeze(mean(z,ndims(z)),ndims(z))
    y0 = squeeze(mean(y,ndims(y)),ndims(y))
    if !centered
        z0 .*= 0.
        y0 .*= 0.
    end
    zc = z .- z0
    yc = y .- y0
    α = (λ*eye(size(z,1)) + zc*zc')\(zc*yc')
    RandomFourierKernelMachine(d,dout,σ,w,α,z0,y0)
end


function Base.call(Φ::RandomFourierKernelMachine,X::AbstractArray) 
    x = mat(X)
    k = size(x,2)  # number of samples
    z = real(exp(Φ.w*x*im)/sqrt(Φ.d))
    y = Φ.α'*(z .- Φ.z0) .+ Φ.y0
    reshape(y, (Φ.dout...,k))
end



