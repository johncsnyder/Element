


abstract Struct

# Base.keys(t::Struct)       = fieldnames(t)
# Base.values(t::Struct)     = [getfield(t,i) for i in fieldnames(t)]
# Base.length(t::Struct)     = length(fieldnames(t))
# Base.start(t::Struct)      = 1
# Base.done(t::Struct, iter) = iter == length(fieldnames(t))+1
# Base.next(t::Struct, iter) = ((fieldnames(t)[iter], getfield(t,iter)), iter += 1)
# Base.endof(t::Struct)      = length(t)
# Base.last(t::Struct)       = (fieldnames(t)[end], t[end])

# Base.getindex(t::Struct, i::Int)              = getfield(t,i)
# Base.getindex(t::Struct, i::UnitRange{Int64}) = slice(t,i)
# Base.getindex(t::Struct, i::Symbol)           = getfield(t,i)
# Base.getindex(t::Struct, i::Symbol, default)  = get(t,i,default)
# Base.get(t::Struct, i::Symbol, default)       = i in keys(t) ? t[i] : default

Base.writemime(io::IO, mime::MIME"text/plain", t::Struct) = show(io,t)

function Base.show(io::IO, t::Struct; cols=100)
    shown_set = get(task_local_storage(), :SHOWNSET, nothing)
    if shown_set === nothing
        shown_set = ObjectIdDict()
        task_local_storage(:SHOWNSET, shown_set)
    end
    t in keys(shown_set) && (print(io, "#= circular reference =#"); return)

    try
        shown_set[t] = true
        fields = fieldnames(t)
        ks = map(string,fields)
        keylen = maximum(map(length,ks))

        for (i,k) in enumerate(fields)
            print(io, Base.rpad(ks[i], keylen))
            print(io, " = ")
            val = Base.with_output_limit(()->sprint(show, getfield(t,k)))
            val = Base._truncate_at_width_or_chars(val, cols - keylen, "\r\n")
            print(io, val)
            print(io, '\n')
        end
    finally
        delete!(shown_set, t)
    end
end

# function show(io::IO, t::Struct)
#     shown_set = get(task_local_storage(), :SHOWNSET, nothing)
#     if shown_set === nothing
#         shown_set = ObjectIdDict()
#         task_local_storage(:SHOWNSET, shown_set)
#     end
#     t in keys(shown_set) && (print(io, "#= circular reference =#"); return)

#     try
#         shown_set[t] = true
#         # if compact
#         #     # show in a Julia-syntax-like form: Dict(k=>v, ...)
#         #     if isempty(t)
#         #         print(io, typeof(t), "()")
#         #     else
#         #         if isleaftype(K) && isleaftype(V)
#         #             print(io, typeof(t).name)
#         #         else
#         #             print(io, typeof(t))
#         #         end
#         #         print(io, '(')
#         #         first = true
#         #         n = 0
#         #         for (k, v) in t
#         #             first || print(io, ',')
#         #             first = false
#         #             show(io, k)
#         #             print(io, "=>")
#         #             show(io, v)
#         #             n+=1
#         #             limit && n >= 10 && (print(io, "…"); break)
#         #         end
#         #         print(io, ')')
#         #     end
#         #     return
#         # end

#         # Otherwise show more descriptively, with one line per key/value pair
#         # rows, cols = sz
#         print(io, summary(t))
#         # isempty(t) && return
#         print(io, ":")
#         # if limit
#         #     rows < 2   && (print(io, " …"); return)
#         #     cols < 12  && (cols = 12) # Minimum widths of 2 for key, 4 for value
#         #     cols -= 6 # Subtract the widths of prefix "  " separator " => "
#         #     rows -= 2 # Subtract the summary and final ⋮ continuation lines

#         #     # determine max key width to align the output, caching the strings
#         #     ks = Array(AbstractString, min(rows, length(t)))
#         #     keylen = 0
#         #     for (i, k) in enumerate(keys(t))
#         #         i > rows && break
#         #         ks[i] = sprint(show, k)
#         #         keylen = clamp(length(ks[i]), keylen, div(cols, 3))
#         #     end
#         # end

#         for (i, (k, v)) in enumerate(t)
#             print(io, "\n  ")
#             # limit && i > rows && (print(io, rpad("⋮", keylen), " => ⋮"); break)

#             # if limit
#                 # key = rpad(_truncate_at_width_or_chars(ks[i], keylen, "\r\n"), keylen)
#             # else
#             key = sprint(show, k)
#             # end
#             print(io, key)
#             print(io, " = ")

#             # if limit
#             #     val = with_output_limit(()->sprint(show, v))
#             #     val = _truncate_at_width_or_chars(val, cols - keylen, "\r\n")
#             #     print(io, val)
#             # else
#             show(io, v)
#             # end
#         end
#     finally
#         delete!(shown_set, t)
#     end
# end



function Base.(:(==))(a::Struct, b::Struct)
    a === b                && return true
    typeof(a) != typeof(b) && return false
    for k in fieldnames(a)
        if getfield(a,k) != getfield(b,k)
            return false
        end
    end
    true
end

function Base.hash(t::Struct, h::UInt)
    h_ = 17
    for k in fieldnames(t)
        h_ = 23*h_ + hash(getfield(t,k), h)
    end
    h_
end


trans(expr::Expr)                     = trans(Val{expr.head}, expr)                 #
trans(sym::Symbol)                    = (sym, :Any, nothing)                        # x
trans(::Type{Val{:(::)}}, expr::Expr) = (expr.args[1], expr.args[2], nothing)       # x::T

function trans(::Type{Val{:kw}}, expr::Expr)                                        # x::T = v
    sym, typ = trans(expr.args[1])
    (sym, typ, expr.args[2])
end




function make_type(base::Symbol, fields::Vector{Symbol})
    name = symbol(base,'_',join(fields,'_'))
    if !isdefined(Element,name)
        types = [symbol('T',i) for i in 1:length(fields)]
        tfields = [Expr(:(::), fields[i], types[i]) for i in 1:length(fields)]
        def = Expr(:type, true, Expr(:(<:), Expr(:curly, name, types...), base), Expr(:block, tfields...))
        eval(def)
    end
    return name
end



function make_tuple(base::Symbol, exprs::Union{Symbol,Expr}...)
    local flag
    n = length(exprs)

    fields = Array{Symbol}(n)
    values = Array{Any}(n)
    types  = Array{Any}(n)

    for i in 1:n
        fields[i],T,v = trans(exprs[i])
        if i == 1; flag = is(v,nothing) end
        is(v,nothing) != flag && error("all values must be specified during construction")
        types[i] = esc(T)
        values[i] = is(T,:Any) ? esc(v) : esc(Expr(:call, :convert, T, v))
    end

    # sort fields alphabetically
    ind = sortperm(fields)
    fields = fields[ind]
    types = types[ind]
    values = values[ind]

    name = make_type(base, fields)  # create type if nonexistent

    if flag
        return Expr(:curly, name, types...)
    else
        return Expr(:call, name, values...)
    end
end




abstract Atom <: Struct

macro atom(exprs::Expr...)
    make_tuple(:Atom, exprs...)
end

macro struct(exprs::Expr...)
    make_tuple(:Struct, exprs...)
end

@generated function merge(a::Struct, b::Struct)
    merge_impl(a,b)
end

function merge_impl{T1<:Struct,T2<:Struct}(a::Type{T1}, b::Type{T2})
    fields = sort(unique([fieldnames(a);fieldnames(b)]))
    T = symbol(last(split(string(super(a)),'.')))
    name = make_type(T, fields)
    vals = [k in fieldnames(b) ? Expr(:call,:getfield,:b,QuoteNode(k)) : Expr(:call,:getfield,:a,QuoteNode(k)) for k in fields]
    Expr(:call, name, vals...)
end


# @doc doc"""
# Create a new Struct with the new value set on it, either overwriting
# the old value or appending a new value.
# This copies the underlying data.
# """ ->
# function setindex{V}( t::Struct, key::Symbol, val::V)
#     nt = eval( make_tuple( super(typeof(t)), [ Expr( :(=>), key, val  )]))
#     return merge( t, nt )
# end

# @doc doc"""
# Create a new Struct with the specifed element removed.
# This copies the underlying data.
# """ ->
# function delete( t::Struct, key::Symbol )
#     nms = filter( x->x!=key, fieldnames( t ) )
#     name = create_tuple( super(typeof(t)), nms )
#     vals = [ getindex( t, nm ) for nm in nms ]
#     return eval( Expr( :call, name, vals... ) )
# end



