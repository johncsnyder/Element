



@generated function findsteps!(I::Array{Tuple{Float64,Float64},1}, f::Function, xl::Float64, xr::Float64,
            fl::Int, fr::Int, fmin::Int, fmax::Int, args...)
    quote 
        if fl == fr || fr <= fmin || fl > fmax
            return
        end
        
        if fr - fl == 1
            push!(I, (xl,xr))  # found an interval (xl,xr)
            return
        end
        
        # subdivide the interval
        xmid = (xl + xr)/2
        fmid = f(xmid,args...)::Int
        
        # recursive search
        findsteps!(I, f, xl, xmid, fl, fmid, fmin, fmax, args...)
        findsteps!(I, f, xmid, xr, fmid, fr, fmin, fmax, args...)
    end
end




doc"""
`findsteps!(I,f,xl,xr,[fmin,fmax];args=())` 

like `findsteps` except that the intervals $(x_l,x_r)$ are appended to `I`.
"""
function findsteps!(I::Array{Tuple{Float64,Float64},1}, f::Function,
        xl::Float64, xr::Float64, fmin::Int, fmax::Int; args=())
    findsteps!(I, f, xl, xr, f(xl,args...)::Int, f(xr,args...)::Int, fmin, fmax, args...)
end

function findsteps!(I::Array{Tuple{Float64,Float64},1}, f::Function, xl::Float64, xr::Float64; args=())
    fl = f(xl,args...)::Int
    fr = f(xr,args...)::Int
    findsteps!(I, f, xl, xr, fl, fr, fl, fr, args...)
end





doc"""
`findsteps(f,xl,xr,[fmin,fmax];args=())` 

takes a function `f(x,args...)`, assumed to be a
monotonically increasing integer function of `x`.

optionally one can specify the range of `f` to search 
with `fmin`, `fmax`

performs a binary search and returns an ordered list of 
intervals $[(x_{l,1},x_{r,1}),(x_{l,2},x_{r,2}),\dots]$
within the specificed range `(xl,xr)` which contains exactly one step.

e.g.
```julia
f = (x,c) -> round(Int,x^2+c)
findsteps(f,0.0,2.0,args=1.0)
```
"""
function findsteps(f::Function, xl::Float64, xr::Float64; args=())
    I = Array{Tuple{Float64,Float64}}(0)
    findsteps!(I, f, xl, xr, args=args)
    I
end

function findsteps(f::Function, xl::Float64, xr::Float64, fmin::Int, fmax::Int; args=())
    I = Array{Tuple{Float64,Float64}}(0)
    findsteps!(I, f, xl, xr, fmin, fmax, args=args)
    I
end




# function findsteps!(I::Array{Tuple{Float64,Float64},1}, f::Function, xl::Float64, xr::Float64,
#             fl::Int, fr::Int, fmin::Int, fmax::Int)
#     if fl == fr || fr <= fmin || fl > fmax
#         return
#     end
    
#     if fr - fl == 1
#         push!(I, (xl,xr))  # found an interval (xl,xr)
#         return
#     end
    
#     # subdivide the interval
#     xmid = (xl + xr)/2
#     fmid = f(xmid)::Int
    
#     # recursive search
#     findsteps!(I, f, xl, xmid, fl, fmid, fmin, fmax)
#     findsteps!(I, f, xmid, xr, fmid, fr, fmin, fmax)
# end


# function findsteps!(I::Array{Tuple{Float64,Float64},1}, f::Function,
#         xl::Float64, xr::Float64, fmin::Int, fmax::Int)
#     findsteps!(I, f, xl, xr, f(xl)::Int, f(xr)::Int, fmin, fmax)
# end

# function findsteps!(I::Array{Tuple{Float64,Float64},1}, f::Function, xl::Float64, xr::Float64)
#     fl = f(xl)::Int
#     fr = f(xr)::Int
#     findsteps!(I, f, xl, xr, fl, fr, fl, fr)
# end

# function findsteps(f::Function, xl::Float64, xr::Float64, fmin::Int, fmax::Int)
#     I = Array{Tuple{Float64,Float64}}(0)
#     findsteps!(I, f, xl, xr, fmin, fmax)
#     I
# end



# doc"""
# `findsteps(f,xl,xr)` 

# takes a function $f(x)$, assumed to be a
# monotonically increasing integer function of $x$.

# performs a binary search and returns an ordered list of 
# intervals $[(x_{l,1},x_{r,1}),(x_{l,2},x_{r,2}),\dots]$
# within the specificed range `(xl,xr)` which contains exactly one step.

# """
# function findsteps(f::Function, xl::Float64, xr::Float64)
#     I = Array{Tuple{Float64,Float64}}(0)
#     findsteps!(I, f, xl, xr)
#     I
# end

