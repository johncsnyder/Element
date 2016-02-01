

using Base.Cartesian: @nexprs






macro multistep(fname,rule)
	@assert rule.head == :(=)
	y = rule.args[1].args[1]  # output variable
    vecs = Set{Symbol}()
    consts = Set{Symbol}()
    inds = Set{Symbol}()
    offsets = Set{Int}()
    extractvars!(rule,vecs,consts,inds,offsets)
    delete!(vecs,y)
    vecs = sort([vecs...])
    consts = sort([consts...])
    @assert length(inds) == 1 "only 1 index variable allowed"
    ind = pop!(inds)
    maxoffset = maximum(offsets)
    minoffset = minimum(offsets)
    nconsts = maxoffset - minoffset
    fname! = symbol(string(fname) * "!")
    
    yarg = :($y::AbstractVector{T})
    yconsts = [symbol(string(y) * "_$i") for i in 1:nconsts]
    yconstargs = [:($c::T) for c in yconsts]
    refargs = [Expr(:(::), v, :(AbstractVector{T})) for v in vecs]
    constargs = [Expr(:(::), c, :T) for c in consts]
    args = [yconstargs;refargs;constargs]
    fsig = esc(Expr(:call, Expr(:curly, fname, :T), args...))
    fsig! = esc(Expr(:call, Expr(:curly, fname!, :T), yarg, args...))
    fcall = esc(Expr(:call, Expr(:curly, fname!, :T), y, yconsts..., vecs..., consts...))
    ind = esc(ind)
    
    quote
        $fsig! = begin
            @nexprs $nconsts j -> $y[j] = $yconsts[j]
            @inbounds for $ind in $(1-minoffset):length($y)-$maxoffset
                $rule
            end
        end

        $fsig = begin
            $(esc(y)) = similar($(esc(vecs[1])))
            $(Expr(:escape, Expr(:call, fname!, y, yconsts..., vecs..., consts...)))
            $(esc(y))
        end
    end
end




function extractvars!(e, vecs, consts, inds, offsets)
    if isa(e, Symbol) && isalpha(string(e))
        push!(consts, e)
        return
    end
    
    if isa(e, Expr)
        if e.head == :ref
            push!(vecs, e.args[1])
            ind = e.args[2]
            
            if isa(ind, Symbol)
                push!(inds, ind)
            elseif isa(ind, Expr)
                op = ind.args[1]
                @assert ind.head == :call && (op == :+ || op == :-)
                offset = op == :+ ? ind.args[3] : -ind.args[3]
                push!(inds, ind.args[2])
                push!(offsets, offset)
            end
            
            return
        end
        
        for i in 1:length(e.args)
            extractvars!(e.args[i], vecs, consts, inds, offsets)
        end
    end
end





macro adamsmoulton3_multistep(name,kernel)
    vecs = Set()
    consts = Set()
    inds = Set()
    offsets = Set()
    extractvars!(kernel,vecs,consts,inds,offsets)
    delete!(vecs,:y)
    delete!(vecs,:f)
    push!(consts,:h)
    vecs = sort([vecs...])
    consts = sort([consts...])
    @assert length(inds) == 1 "only 1 index variable allowed"
    ind = pop!(inds)
    yarg = :(y::AbstractVector)
    vecargs = [Expr(:(::), v, :(AbstractVector)) for v in vecs]
    fsig = esc(Expr(:call, :($name), yarg, vecargs..., consts...))

	y = esc(:y)
	i = esc(ind)
	f = esc(:f)
    
    quote
        $fsig = begin
            $f = similar($y)
            
            $y[1] = 0.
            $f[1] = $(replace(kernel,ind,1))
            $y[2] = $y[1] + h*$f[1]
            $f[2] = $(replace(kernel,ind,2))
            $y[3] = $y[2] + h*($f[2]*3./2. - $f[1]/2.)
            $f[3] = $(replace(kernel,ind,3))
            for $i in 4:length($y)
                $y[$i] = $y[$i-1] + h*($f[$i-1]*23./12. - $f[$i-2]*4./3. + $f[$i-3]*5./12.)
                $f[$i] = $(replace(kernel,ind,i))
            end
        end
    end
end




function replace!(e,var,val)
    !isa(e,Expr) && return
    for i in 1:length(e.args)
        if isa(e.args[i],Expr)
            replace!(e.args[i],var,val)
        elseif e.args[i] == var
            e.args[i] = val
        end
    end
end

function replace(e,var,val)
    ec = copy(e)
    replace!(ec,var,val)
    ec
end


function isvariable(s)
    if isa(s, Symbol) && isalpha(string(s))
        return true
    end
    
    if isa(s, Expr)
        return reduce(|, map(isvariable, s.args))
    end
    
    return false
end


