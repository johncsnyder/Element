
typealias Stats Dict{AbstractString,Vector{UInt64}}



macro stats(name,expr)
    name = esc(name)
    expr = esc(expr)
    stats = esc(:stats)
    quote
        local gc_stats = Base.gc_num()
        local elapsedtime = Base.time_ns()
        $expr
        elapsedtime = Base.time_ns() - elapsedtime
        local diff = Base.GC_Diff(Base.gc_num(),gc_stats)
        if $name in keys($stats)
            $stats[$name] += [elapsedtime,diff.allocd,diff.total_time,Base.gc_alloc_count(diff)]
        else
            $stats[$name] = [elapsedtime,diff.allocd,diff.total_time,Base.gc_alloc_count(diff)]
        end
    end
end


function print_stats(stats::Stats)
    for k in sort(collect(keys(stats)), by=k->stats[k][1], rev=true)  # sort by elapsed time
        @printf "%-30s" k
        Base.time_print(stats[k]...)
    end
end







