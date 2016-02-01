


# extending base functionality and basic utilities



import Base: +, *, .+, .*, ./, .^, .-



typealias RealArray{T<:Real,N} AbstractArray{T,N}
typealias RealVector{T<:Real} AbstractArray{T,1}
typealias RealMatrix{T<:Real} AbstractArray{T,2}

typealias ComplexArray{T<:Complex,N} AbstractArray{T,N}
typealias ComplexVector{T<:Complex} AbstractArray{T,1}
typealias ComplexMatrix{T<:Complex} AbstractArray{T,2}

typealias IntegerArray{T<:Integer,N} AbstractArray{T,N}
typealias IntegerVector{T<:Integer} AbstractArray{T,1}
typealias IntegerMatrix{T<:Integer} AbstractArray{T,2}



# if catdims are not specificed, concatenate along new axis on end
Base.cat(A::AbstractArray...) = cat(maximum(map(ndims,A))+1, A...)
Base.cat{T}(catdims, A::AbstractArray{T}...) = catdims < 0 ? Base.cat_t(maximum(map(ndims,A))+catdims+1, T, A...) : Base.cat_t(catdims, T, A...)




function memsize(a)
    f = fieldnames(a)
    length(f) > 0 ? sum([memsize(getfield(a,k)) for k in f]) : sizeof(a)
end

memsize(a::Function) = sizeof(a)

function memsize(a::AbstractArray)
    f = fieldnames(a)
    length(f) > 0 ? sum([memsize(getfield(a,k)) for k in f]) : 
        (length(a) > 0 ? sum([isdefined(a,i) ? memsize(a[i]) : 0 for i in 1:length(a)]) : 0)
end

memsize(a::Associative) = sum([memsize(a[k]) for k in keys(a)])






# function flatten(x::AbstractArray, dim::Int64)
#     if dim > 0
#         n = size(x)[1:dim]
#         d = prod(size(x)[dim+1:end])
#         return reshape(x,tuple(n...,d))
#     else
#         n = prod(size(x)[1:end+dim])
#         d = size(x)[end+dim+1:end]
#         return reshape(x,tuple(n,d...))
#     end
# end


# mat(x::AbstractArray) = flatten(x,-1)





@inline mat(x::AbstractMatrix) = x

@generated function mat(x::AbstractArray)
    d = ndims(x)
    si = [_ref(:s,i) for i in 1:d-1]             # [s[1], ..., s[d-1]]
    n = Expr(:call,:*,si...)                    # s[1] * s[2] * ... * s[d-1]
    m = _ref(:s,d)                               # s[d]
    quote
        s = size(x)
        reshape(x,($n,$m))
    end
end







macro defaults(d, expr)
    d = esc(d)
    out = Expr(:block)
    for line in expr.args
        if line.head == :(=)
            lhs = line.args[1]
            s = (:head in fieldnames(lhs) && lhs.head == :(::)) ? lhs.args[1] : lhs
            t = (:head in fieldnames(lhs) && lhs.head == :(::)) ? esc(lhs.args[2]) : :Any
            k = QuoteNode(s)
            s = esc(s)
            v = esc(line.args[2])
            push!(out.args,
                :(if $k in keys($d); $s = $d[$k]::$t; else $s = $d[$k] = ($v)::$t end)
            )
        end
    end
    out 
end




struct(; kv...) = Dict{Symbol,Any}(kv)
struct(s::Dict{Symbol,Any}; kv...) = merge(s, Dict{Symbol,Any}(kv))





# function expandshape(As::Union{AbstractArray,Number}...)
#     d = maximum(map(ndims,As))
#     shps = [map(size,As)...]
#     d = maximum(map(length,shps))
#     lshps = filter(shp->length(shp)==d, shps) |> collect
#     tshp = [maximum([shp[i] for shp in lshps]) for i in 1:d]
#     nshps = Array{Array{Int64}}(length(As))
    
#     for (i,s) in enumerate(shps)
#         shp = nshps[i] = fill(1, max(length(s),d))
#         ind = 1
#         for (j,l) in filter(a->a[2]>1, enumerate(s))
#             k = findnext(tshp,l,max(ind,j))
#             if k == 0
#                 k = findnext(tshp,1,max(ind,j))
#             end
#             if k < j || k == 0
#                 throw(DimensionMismatch("array could not be broadcast to match destination"))
#             end
#             shp[k] = l
#             ind = k+1
#         end
#     end

#     tuple([nshps[i] != [] ? reshape(A, tuple(nshps[i]...)) : A for (i,A) in enumerate(As)]...)
# end




# (.*)(As::AbstractArray...) = broadcast(*, expandshape(As...)...)
# (.*)(A::AbstractArray, B::AbstractArray) = broadcast(*, expandshape(A,B)...)

# function .+(As::AbstractArray...)
#     As_ = expandshape(As...)
#     broadcast!(+, Array(Base.Broadcast.eltype_plus(As_...), 
#         Base.Broadcast.broadcast_shape(As_...)), As_...)
# end

# function ./(A::AbstractArray, B::AbstractArray)
#     A_, B_ = expandshape(A,B)
#     broadcast!(/, Array(Base.Broadcast.type_div(eltype(A_), eltype(B_)), Base.Broadcast.broadcast_shape(A_, B_)), A_, B_)
# end

# function .^(A::AbstractArray, B::AbstractArray)
#     A_, B_ = expandshape(A,B)
#     broadcast!(^, Array(Base.Broadcast.type_pow(eltype(A_), eltype(B_)), Base.Broadcast.broadcast_shape(A_, B_)), A_, B_)
# end

# function .-(A::AbstractArray, B::AbstractArray)
#     A_, B_ = expandshape(A,B)
#     broadcast!(-, Array(Base.Broadcast.type_minus(eltype(A_), eltype(B_)), Base.Broadcast.broadcast_shape(A_,B_)), A_, B_)
# end




# function Base.(:*){T}(A::AbstractArray{T,2}, B::AbstractArray{T,3})
#     shp = size(B)
#     reshape(A*reshape(B, (shp[1],prod(shp[2:end]))), (size(A,1),shp[2:end]...))
# end


# function Base.(:*){T}(A::AbstractArray{T,4}, B::AbstractArray{T,2})
#     shp = size(A)
#     reshape(reshape(A, (prod(shp[1:end-1]),shp[end]))*B, (shp[1:end-1]...,size(B,2)))
# end


# function Base.(:*){T}(A::AbstractArray{T,3}, B::AbstractArray{T,2})
#     shp = size(A)
#     reshape(reshape(A, (prod(shp[1:end-1]),shp[end]))*B, (shp[1:end-1]...,size(B,2)))
# end



# function Base.(:*){T}(A::AbstractArray{T,4}, B::AbstractArray{T,1})
#     shp = size(A)
#     reshape(reshape(A, (prod(shp[1:end-1]),shp[end]))*B, (shp[1:end-1]...))
# end

_ref(sym,i) = Expr(:ref,sym,i)


@generated function Base.(:*)(x::AbstractArray, y::AbstractArray)
    dx = ndims(x)
    dy = ndims(y)
    sxi = [_ref(:sx,i) for i in 1:dx-1]
    syi = [_ref(:sy,i) for i in 2:dy]
    n = Expr(:call,:*,sxi...)
    m = _ref(:sx,dx)
    k = _ref(:sy,1)
    l = Expr(:call,:*,syi...)
    s = Expr(:tuple, sxi..., syi...)
    xm = dx > 2    ?    :(reshape(x,($n,$m)))    :    :x
    ym = dy > 2    ?    :(reshape(y,($k,$l)))    :    :y
    quote
        sx = size(x)
        sy = size(y)
        reshape($xm * $ym, $s)
    end
end


function Base.SparseMatrix.sparse{Tv}(A::AbstractArray{Tv,2}, B::SparseMatrixCSC)
    I,J = findn(B)
    ind = I+(J-1)*size(A,1)
    SparseMatrixCSC(B.m,B.n,B.colptr,B.rowval,A[ind])
end



# reverse slice
function rslice{T,N}(A::AbstractArray{T,N}, I::Base.ViewIndex...)
    n = ndims(A)
    m = length(I)
    if m < n
        I = tuple(I..., fill(Colon(),n-m)...)
    end
    slice(A, reverse(I)...)
end

# reverse sub
function rsub{T,N}(A::AbstractArray{T,N}, I::Base.ViewIndex...)
    n = ndims(A)
    m = length(I)
    if m < n
        I = tuple(I..., fill(Colon(),n-m)...)
    end
    sub(A, reverse(I)...)
end







function squeezeall(A::AbstractArray)
    s = tuple(filter(d -> d>1, size(A))...)  # select all dimensions greater than 1
    s == () ? reshape(A, s)[] : reshape(A, s)  # reshape to remove those dimensions
end


fermidirac(E::Union{Real,RealVector},mu::Real,T::Real) = 1./(exp((E-mu)./T) + 1)



# # converge a function f(x) wrt some parameter x
# # assuming an error approx c0 + c1 x + c2 x^2 + ... as x->0

# function converge(f_,x0,xiter;ftol=1e-6,xmin=1e-6,maxsteps=100,limit=:zero)
#     local p
#     x = Float64[x0]
#     f = Float64[]
#     ftol = 1e-6
#     i = 0
#     while true
#         push!(f,f_(x[end]))
#         if i > maxsteps
#             warn("did not converge. stopping")
#             break
#         end
#         if length(x) > 2
#             if limit == :zero
#                 p = poly_fit(x[end-2:end],f[end-2:end],2)
#             elseif limit == :infinity
#                 p = poly_fit(1./x[end-2:end],f[end-2:end],2)
#             end
#             df = p[1] - f[end]
#             # println("x=$(x[end]), f=$(f[end]), df=$(df)")
#             if abs(df) < ftol
#                 break
#             end
#         end
#         push!(x,xiter(x[end]))
#         i += 1
#     end
#     x[end], p[1], p[1] - f[end]
# end




