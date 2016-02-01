


abstract Func{Ti,To}

typealias ScalarFunc{T} Func{T,T}
typealias ScalarFunc2d{T} Func{Point2{T},T}
typealias ScalarFunc3d{T} Func{Point3{T},T}


typealias FuncArray{T<:Func} AbstractArray{T}
typealias ScalarFuncArray{T<:ScalarFunc} AbstractArray{T}
typealias ScalarFunc2dArray{T<:ScalarFunc2d} AbstractArray{T}
typealias ScalarFunc3dArray{T<:ScalarFunc3d} AbstractArray{T}




macro vectorize(f,Ti,To)
    f = esc(f)
    Ti = esc(Ti)
    To = esc(To)

    quote
        function $f(dst::AbstractArray{$To}, c::AbstractArray{$Ti})
            @assert size(dst) == size(c)
            for i in eachindex(c)
                dst[i] = $f(c[i])
            end
        end

        function $f(c::AbstractArray{$Ti})
            dst = Array{$To}(size(c))
            $f(dst,c)
            dst
        end
    end
end



macro vectorize2(f,T,Ti,To)
    f = esc(f)
    T = esc(T)
    Ti = esc(Ti)
    To = esc(To)

    quote
        function $f(f::$T, dst::AbstractArray{$To}, c::AbstractArray{$Ti})
            @assert size(dst) == size(c)
            for i in eachindex(c)
                dst[i] = $f(f,c[i])
            end
        end

        function $f(f::$T, c::AbstractArray{$Ti})
            dst = Array{$To}(size(c))
            $f(f,dst,c)
            dst
        end
    end
end



