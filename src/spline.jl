


abstract Spline{T} <: ScalarFunc{T}


# cubic 1D spline
type Spline1D{T,V} <: Spline{T}
    x::V
    y::V
    y2::V
    Spline1D(x::AbstractVector{T}, y::AbstractVector{T}, y2::AbstractVector{T}) = new(x,y,y2)
end


Spline1D{T}(x::AbstractVector{T}, y::AbstractVector{T}, y2::AbstractVector{T}) = 
    Spline1D{eltype(x),typeof(x)}(x,y,y2)



typealias Spline1Df Spline1D{Float64,Vector{Float64}}


@fastmath function Spline1D(x::AbstractVector{Float64}, y::AbstractVector{Float64})
    @assert length(x) == length(y)
    local sig,p
    y2 = similar(y)
    n = length(x)
    u = Array{Float64}(n-1)
    y2[1] = u[1] = 0.0
    @inbounds for i in 2:n-1
        sig = (x[i]-x[i-1])/(x[i+1]-x[i-1])
        p = sig*y2[i-1]+2.0
        y2[i]=(sig-1.0)/p
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1])
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p
    end
    qn=un=0.0
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0)
    @inbounds for k in n-1:-1:1  # This is the backsubstitution loop of the tridiagonal
        y2[k]=y2[k]*y2[k+1]+u[k];
    end
    Spline1D(x,y,y2)
end


function Spline1D(x::AbstractVector{Float64}, y::AbstractVector{Float64}, num::Int)
    n = length(x)
    ind = num < n ? round(Int,linspace(1,n,num)) : Colon()
    Spline1D(x[ind],y[ind])
end


Base.show(io::IO, f::Spline1D) = write(io, "Spline1D $(length(f.x))-points on interval [$(f.x[1]),$(f.x[end])]")



@fastmath function locate(xx::AbstractVector{Float64}, x::Float64)
    local k
    klo,khi = 1,length(xx)
    while khi-klo > 1
        k = (khi+klo) >> 1  # mean of khi and klo
        if xx[k] > x
            khi = k
        else
            klo = k
        end
    end
    klo,khi
end



@fastmath function Base.call(f::Spline1D, x::Float64)
    klo,khi = locate(f.x,x)
    h = f.x[khi] - f.x[klo]  # @assert h != 0.0
    a = (f.x[khi]-x)/h
    b = (x-f.x[klo])/h
    a*f.y[klo] + b*f.y[khi] + ((a*a*a-a)*f.y2[klo]+(b*b*b-b)*f.y2[khi])*(h*h)/6.
end


@vectorize2 Base.call Spline1D Float64 Float64


# function Base.call(f::Spline1D, dst::AbstractVector{Float64}, x::AbstractVector{Float64})
#     @assert length(dst) == length(x)
#     for i in 1:length(x)
#         dst[i] = f(x[i])
#     end
# end

# function Base.call(f::Spline1D, x::AbstractVector{Float64})
#     dst = similar(x)
#     f(dst,x)
#     dst
# end




@fastmath function deriv1(f::Spline1D, x::Float64)
    klo,khi = locate(f.x,x)
    h = f.x[khi] - f.x[klo]  # @assert h != 0.0
    a = (f.x[khi]-x)/h
    b = (x-f.x[klo])/h
    (f.y[khi] - f.y[klo])/h + ((1.-3.*a*a)*f.y2[klo] + (3*b*b-1.)*f.y2[khi])*h/6.
end

@vectorize2 deriv1 Spline1D Float64 Float64

# function deriv1(f::Spline1D, dst::AbstractVector{Float64}, x::AbstractVector{Float64})
#     @assert length(dst) == length(x)
#     @inbounds for i in 1:length(x)
#         dst[i] = deriv1(f,x[i])
#     end
# end

# function deriv1(f::Spline1D, x::AbstractVector{Float64})
#     dst = similar(x)
#     deriv1(f,dst,x)
#     dst
# end


@fastmath function deriv2(f::Spline1D, x::Float64)
    klo,khi = locate(f.x,x)
    h = f.x[khi] - f.x[klo]  # @assert h != 0.0
    a = (f.x[khi]-x)/h
    b = (x-f.x[klo])/h
    a*f.y2[klo] + b*f.y2[khi]
end

@vectorize2 deriv2 Spline1D Float64 Float64


# function deriv2(f::Spline1D, dst::AbstractVector{Float64}, x::AbstractVector{Float64})
#     @assert length(dst) == length(x)
#     @inbounds for i in 1:length(x)
#         dst[i] = deriv2(f,x[i])
#     end
# end

# function deriv2(f::Spline1D, x::AbstractVector{Float64})
#     dst = similar(x)
#     deriv2(f,dst,x)
#     dst
# end




@fastmath function eval_deriv1(f::Spline1D, x::Float64)
    klo,khi = locate(f.x,x)
    h = f.x[khi] - f.x[klo]  # @assert h != 0.0
    a = (f.x[khi]-x)/h
    b = (x-f.x[klo])/h
    a*f.y[klo] + b*f.y[khi] + ((a*a*a-a)*f.y2[klo]+(b*b*b-b)*f.y2[khi])*(h*h)/6., 
        (f.y[khi] - f.y[klo])/h + ((1.-3.*a*a)*f.y2[klo] + (3*b*b-1.)*f.y2[khi])*h/6.
end


@fastmath function eval_deriv2(f::Spline1D, x::Float64)
    klo,khi = locate(f.x,x)
    h = f.x[khi] - f.x[klo]  # @assert h != 0.0
    a = (f.x[khi]-x)/h
    b = (x-f.x[klo])/h
    a*f.y[klo] + b*f.y[khi] + ((a*a*a-a)*f.y2[klo]+(b*b*b-b)*f.y2[khi])*(h*h)/6., 
        a*f.y2[klo] + b*f.y2[khi]
end


@fastmath function eval_deriv12(f::Spline1D, x::Float64)
    klo,khi = locate(f.x,x)
    h = f.x[khi] - f.x[klo]  # @assert h != 0.0
    a = (f.x[khi]-x)/h
    b = (x-f.x[klo])/h
    a*f.y[klo] + b*f.y[khi] + ((a*a*a-a)*f.y2[klo]+(b*b*b-b)*f.y2[khi])*(h*h)/6., 
        (f.y[khi] - f.y[klo])/h + ((1.-3.*a*a)*f.y2[klo] + (3*b*b-1.)*f.y2[khi])*h/6.,
            a*f.y2[klo] + b*f.y2[khi]
end





precompile(Spline1D, (Vector{Float64}, Vector{Float64}))
precompile(locate, (Vector{Float64}, Float64))
precompile(Base.call, (Spline1D, Float64))
precompile(Base.call, (Spline1D, Vector{Float64}, Float64))
precompile(Base.call, (Spline1D, Vector{Float64}))
precompile(deriv1, (Spline1D, Float64))
precompile(deriv1, (Spline1D, Vector{Float64}, Float64))
precompile(deriv1, (Spline1D, Vector{Float64}))
precompile(deriv2, (Spline1D, Float64))
precompile(deriv2, (Spline1D, Vector{Float64}, Float64))
precompile(deriv2, (Spline1D, Vector{Float64}))
precompile(eval_deriv1, (Spline1D, Float64))
precompile(eval_deriv2, (Spline1D, Float64))
precompile(eval_deriv12, (Spline1D, Float64))




