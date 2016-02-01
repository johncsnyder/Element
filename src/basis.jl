



# typealias Basis{T<:ScalarFunc3d{Float64},N} Array{T,N}

typealias Basis{T<:Func,N} AbstractArray{T,N}


_nindices(sym, n) = [symbol(sym, :_, i) for i = 1:n]  # e.g. i_1,...,i_N
_nref(sym, ind...) = Expr(:ref, sym, ind...)  # e.g. A[i_1,...,i_N]



function Base.super(T,S)
    @assert issubtype(T,S)
    while issubtype(super(T),S) && !is(T,Any)
        T = super(T)
    end
    T
end


for f in (:(Base.call), :laplacian)

    eval(quote

        @generated function $f{Ti,To,T<:Func,N,M}(b::Basis{T,N}, dst::AbstractArray{To}, c::AbstractArray{Ti,M})
            Ti_,To_ = super(T,Func).parameters
            @assert Ti == Ti_ && To == To_ "invalid types for $(super(T,Func))"
            i      = _nindices(:i,M)
            j      = _nindices(:j,N)
            c_i    = _nref(:c,i...)
            b_j    = _nref(:b,j...)
            dst_ij = _nref(:dst,i...,j...)
            quote
                Base.Cartesian.@nloops $N j b begin
                    Base.Cartesian.@nloops $M i c begin
                        $dst_ij = $($f)($b_j,$c_i)
                    end
                end
            end
        end

        @generated function $f{Ti,T<:Func,N,M}(b::Basis{T,N}, c::AbstractArray{Ti,M})
            Ti_,To = super(T,Func).parameters
            @assert Ti == Ti_ "invalid input type for $(super(T,Func))"
            quote
                dst = Array{$To}((size(c)...,size(b)...))
                $($f)(b,dst,c)
                dst
            end
        end

        @generated function $f{Ti,T<:Func,N}(b::Basis{T,N}, c::Ti)
            Ti_,To = super(T,Func).parameters
            @assert Ti == Ti_ "invalid input type for $(super(T,Func))"
            quote
                dst = Array{$To}(size(b))
                for i in eachindex(b)
                    dst[i] = $($f)(b[i],c)
                end
                dst
            end
        end

        @generated function $f{T<:Func,N,X}(b::Basis{T,N}, c::X...)
            Ti,To = super(T,Func).parameters
            @assert issubtype(Ti,Point) && X == first(Ti.parameters) "invalid input type for $(super(T,Func))"
            quote
                dst = Array{$To}(size(b))
                for i in eachindex(b)
                    dst[i] = $($f)(b[i],c...)
                end
                dst
            end
        end

    end)

end

# @generated function Base.call{T,N,M}(b::Basis{T,M}, dst::AbstractArray{Float64}, c::AbstractArray{Point3{Float64},N})
#     i      = _nindices(:i,N)
#     j      = _nindices(:j,M)
#     c_i    = _nref(:c,i...)
#     b_j    = _nref(:b,j...)
#     dst_ij = _nref(:dst,i...,j...)
#     quote
#         Base.Cartesian.@nloops $M j b begin
#             Base.Cartesian.@nloops $N i c begin
#                 $dst_ij = $b_j($c_i)
#             end
#         end
#     end
# end


# function Base.call(b::Basis, c::AbstractArray{Point3{Float64}})
#     dst = Array{Float64}((size(c)...,size(b)...))
#     b(dst,c)
#     dst
# end







# Base.call(b::Basis, x::Float64, y::Float64, z::Float64) = reshape(Float64[b[i](x,y,z) for i in eachindex(b)], size(b))
# Base.call(b::Basis, p::Point3{Float64}) = reshape(Float64[b[i](p) for i in eachindex(b)], size(b))



# @generated function laplacian{T,N,M}(b::Basis{T,M}, dst::AbstractArray{Float64}, c::AbstractArray{Point3{Float64},N})
#     i      = _nindices(:i,N)
#     j      = _nindices(:j,M)
#     c_i    = _nref(:c,i...)
#     b_j    = _nref(:b,j...)
#     dst_ij = _nref(:dst,i...,j...)
#     quote
#         Base.Cartesian.@nloops $M j b begin
#             Base.Cartesian.@nloops $N i c begin
#                 $dst_ij = laplacian($b_j,$c_i)
#             end
#         end
#     end
# end


# function laplacian(b::Basis, c::AbstractArray{Point3{Float64}})
#     dst = Array{Float64}((size(c)...,size(b)...))
#     laplacian(b,dst,c)
#     dst
# end




# @generated function add!{T,N,M}(b::Basis{T,M}, dst::AbstractArray{Float64}, c::AbstractArray{Point3{Float64},N})
#     i     = _nindices(:i,N)
#     j     = _nindices(:j,M)
#     c_i   = _nref(:c,i...)
#     b_j   = _nref(:b,j...)
#     dst_i = _nref(:dst,i...)
#     quote
#         @assert size(dst) == size(c)
#         Base.Cartesian.@nloops $N i c begin
#             Base.Cartesian.@nloops $M j b begin
#                 $dst_i += $b_j($c_i)
#             end
#         end
#     end
# end


# @generated function sum!{N,M}(b::Basis{M}, dst::AbstractArray{Float64}, c::AbstractArray{Point3{Float64},N})
#     i = _nindices(:i,N)
#     j = _nindices(:j,M)
#     c_i = _nref(:c,i...)
#     b_j = _nref(:b,j...)
#     dst_i = _nref(:dst,i...)
#     quote
#         @assert size(dst) == size(c)
#         Base.Cartesian.@nloops $N i c begin
#             $dst_i = 0.0
#             Base.Cartesian.@nloops $M j b begin
#                 $dst_i += laplacian($b_j,$c_i)
#             end
#         end
#     end
# end

# function sum(b::Basis, c::AbstractArray{Point3{Float64}})
#     dst = zeros(c)
#     add!(b,dst,c)
#     dst
# end






