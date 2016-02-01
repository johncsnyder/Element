


typealias VecVec{T} Vector{Vector{T}}
typealias VecVecVec{T} Vector{Vector{Vector{T}}}


nested_array_size{T<:AbstractVector}(x::Vector{T}) = [nested_array_size(v) for v in x]
nested_array_size(x::AbstractVector) = length(x)

nested_array_length{T<:AbstractVector}(x::Vector{T}) = sum([nested_array_length(x) for x in x])
nested_array_length(x::AbstractVector) = length(x)

nested_array_type(T::DataType, d::Int) = d == 0 ? T : Vector{nested_array_type(T,d-1)}

function nested_array_depth(T::DataType)
    d = 0
    while issubtype(T,AbstractArray)
        T = T.parameters[1]
        d += 1
    end
    d
end

nested_array_depth(x::AbstractArray) = nested_array_depth(typeof(x))




function nested_array_eltype(T::DataType)
    while issubtype(T,AbstractArray)
        T = eltype(T)
    end
    T::DataType
end

nested_array_eltype(x::AbstractArray) = nested_array_eltype(typeof(x))



function rand_nested_array(d,len)
    if d == 1
        return rand(rand(len))
    else
        T = nested_array_type(Float64,d-1)
        T[rand_nested_array(d-1,len) for _ in 1:rand(len)]
    end
end



reduce_size{T<:Integer}(shp::Vector{T}) = sum(shp)
reduce_size{T}(shp::VecVec{T}) = [reduce_size(s) for s in shp]




function reduce_ptr{T<:Integer}(ptr::VecVec{T})
    p = Array{Int}(length(ptr)+1)
    for i in 1:length(ptr)
        p[i] = ptr[i][1]
    end
    p[end] = ptr[end][end]
    p
end

reduce_ptr{T<:Integer}(ptr::VecVecVec{T}) = [reduce_ptr(p) for p in ptr]
reduce_ptr{T<:Integer}(ptr::Vector{T}) = 1





@generated function expand_size(shp)
    d = nested_array_depth(shp)
    quote
        shps = Any[shp]
        for _ in 1:$d
            shp = nested_array_size(shp)
            insert!(shps,1,shp)
        end
        tuple(shps...)
    end
end

expand_size(shp::Integer) = (shp,)



function index_nested_array!{T<:AbstractVector}(x::Vector{T}, val, k)
    ptr = Array{Int}(length(x)+1)
    ptr[1] = k[1]
    for (i,v) in enumerate(x)
        ptr[i+1] = (k[1] += length(v))
        val[ptr[i]:(ptr[i+1]-1)] = v
    end
    ptr
end

index_nested_array!{T<:AbstractVector}(x::VecVec{T},val,k) = [index_nested_array!(v,val,k) for v in x]
index_nested_array!(x,val) = index_nested_array!(x,val,[1])

@generated function index_nested_array(x)
    Tv = nested_array_eltype(x)
    quote
        val = Array{$Tv}(nested_array_length(x))
        ptr = index_nested_array!(x,val)
        val,ptr
    end
end


function combine_size(shp)
    n = length(shp)
    d = length(shp[1])
    T = typeof(shp[1]).types
    cshp = Any[length(shp)]
    for j in 1:d
        push!(cshp,T[j][s[j] for s in shp])
    end
    tuple(cshp...)
end



add_offset{T<:Integer}(x::Vector{T},i,k) = begin k[1] = x[end]; x+i-1 end
add_offset{T<:Integer}(x::VecVec{T},i,k) = [add_offset(v,i,k) for v in x]

function combine_ptr(ptr)
    k = [1]
    T = typeof(ptr[1])
    T[add_offset(p,k[1],k) for p in ptr]
end





type RaggedArray{Tv,N,Tp,Ts} <: DenseVector{Tv}
    val::Vector{Tv}
    ptr::Tp
    shp::Ts
end



@generated function RaggedArray(x::AbstractVector)
    Tv = nested_array_eltype(x)
    N = nested_array_depth(x)
    quote
        val,ptr = index_nested_array(x)
        shp = nested_array_size(x)
        ptrs = Any[ptr]
        shps = Any[expand_size(shp)]
        for _ in 1:$(N-1)
            ptr = reduce_ptr(ptr)
            shp = reduce_size(shp)
            insert!(ptrs,1,ptr)
            insert!(shps,1,expand_size(shp))
        end
        ptr = tuple(ptrs...)
        shp = tuple(shps...)
        RaggedArray{$Tv,$N,typeof(ptr),typeof(shp)}(val,ptr,shp)
    end
end



function RaggedArray{Tv,N}(x::RaggedArray{Tv,N}...)
    val = vcat([v.val for v in x]...)

    ptr = Any[combine_ptr([v.ptr[i] for v in x]) for i in 2:N]
    insert!(ptr,1,reduce_ptr(ptr[1]))
    insert!(ptr,1,1)
    ptr = tuple(ptr...)

    shp = Any[combine_size([v.shp[i] for v in x]) for i in 1:N]
    insert!(shp,1,(reduce_size(shp[1][end]),))
    shp = tuple(shp...)

    RaggedArray{Tv,N,typeof(ptr),typeof(shp)}(val,ptr,shp)
end




Base.eltype{Tv}(a::RaggedArray{Tv}) = Tv
Base.length(ra::RaggedArray) = length(ra.val)
Base.ndims(::RaggedArray) = 1
Base.eachindex(ra::RaggedArray) = eachindex(ra.val)
Base.linearindexing{T<:RaggedArray}(::Type{T}) = Base.LinearFast()
Base.size(ra::RaggedArray) = size(ra.val)
Base.size(ra::RaggedArray, i::Integer) = size(ra.val, i)




Base.similar{Tv,N,Tp,Ts}(ra::RaggedArray{Tv,N,Tp,Ts}) = RaggedArray{Tv,N,Tp,Ts}(similar(val),copy(ra.ptr),copy(ra.shp))




_ref(sym, i) = Expr(:ref, sym, i)  # e.g. A[i_1,...,i_N]
_nestedref(sym, i) = _ref(sym, i)
_nestedref(sym, i...) = _nestedref(_nestedref(sym, i[1]), i[2:end]...)  # A[i_1][i_2]...[i_n]
_nindices(sym, n) = [symbol(sym, :_, i) for i = 1:n]  # [i_1,...,i_n]



@generated function len{Tv,N}(ra::RaggedArray{Tv,N}, d::Integer, i::Integer...)
    n = length(i)                                               # number of indices
    ind = [_ref(:i,j) for j in 1:n]                             # [i_1,...,i_n]
    len = _nestedref(:(ra.shp), :d, n+1, ind...)                # ra.shp[d][2][i[1]], ra.shp[d][3][i[1]][i[2]], etc...
    quote
        $len::Int
    end
end




Base.getindex(ra::RaggedArray, i...) = error("not implemented")
Base.getindex(ra::RaggedArray, i::Integer) = ra.val[i]



@generated function Base.getindex{Tv,N}(ra::RaggedArray{Tv,N}, i::Integer...)
    d = length(i)                                           # number of indices
    d > N && throw(BoundsError)
    ind = [_ref(:i,j) for j in 1:d]                         # [i_1,...,i_d]
    offset = _nestedref(:(ra.ptr), d, ind[end:-1:2]...)     # ra.ptr[1][i[2]], ra.ptr[2][i[3]][i[2]], ra.ptr[3][i[4]][i[3]][i[2]], etc...
    quote
      ra.val[$offset + i[1] - 1]                            # no bounds check ?
    end
end

@generated function Base.getindex{Tv,N}(ra::RaggedArray{Tv,N}, ::Colon, i::Integer...)
    d = length(i)+1
    d > N && throw(BoundsError)
    ind = reverse([_ref(:i,j) for j in 1:length(i)])
    a = _nestedref(:(ra.ptr), d, ind...)
    ind[end] = :($(ind[end])+1)
    b = _nestedref(:(ra.ptr), d, ind...)
    quote
      ra.val[$a:$b-1]
    end
end




# # multiply by weights w and then sum over last dimension
# function sum{Tv,N,Tp,Ts}(ra::RaggedArray{Tv,N,Tp,Ts}, w::AbstractVector)
#     raw = RaggedArray{Tv,N,Tp,Ts}(ra.val.*w,ra.ptr,ra.shp)
#     [sum(raw[:,i]) for i in 1:len(raw,N)]
# end

