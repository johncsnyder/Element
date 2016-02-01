


type MultiGrid{T,A<:AbstractVector,B<:AbstractVector,C<:AbstractVector} <: AbstractGrid3d{T,1}
    c::A  # coords
    w::B  # weights
    p::C  # partition function
    atomptr::Vector{Int}  # atom i contained in atomptr[i]:atomptr[i+1]-1
    function MultiGrid(c::AbstractVector{Point3{T}},w::AbstractVector{T},
            p::AbstractVector{T},atomptr::Vector{Int})
        @assert size(c) == size(w)
        new(c,w,p,atomptr)
    end
end

MultiGrid{T}(c::AbstractVector{Point3{T}},w::AbstractVector{T},
        p::AbstractVector{T},atomptr::Vector{Int}) =
    MultiGrid{eltype(w),typeof(c),typeof(w),typeof(p)}(c,w,p,atomptr)

# empty grid
MultiGrid() = MultiGrid(Point3f[],Float64[],Float64[],Int[])

eachatom(g::MultiGrid) = [g.atomptr[i]:g.atomptr[i+1]-1 for i in 1:length(g.atomptr)-1]




function multigrid(atoms::Atoms)
    c = vcat([vec(at[:g]) for at in atoms]...)
    p = vcat([vec(at[:p](at[:g])) for at in atoms]...)
    w = p.*vcat([vec(at[:g].w) for at in atoms]...)
    atomptr = [1; cumsum([length(at[:g]) for at in atoms])+1...]
    MultiGrid(c,w,p,atomptr)
end








