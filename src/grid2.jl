

type Grid{T,N,d,P<:AbstractArray,Q<:AbstractArray} <: AbstractArray{T,N}
    c::P
    w::Q
    function Grid{V}(c::AbstractArray{T,N},w::AbstractArray{V,N})
        @assert size(c) == size(w)
        new(c,w)
    end
end

Grid{T<:Tuple}(c::AbstractArray{T},w::AbstractArray) = Grid{eltype(c),ndims(c),length(T.types),typeof(c),typeof(w)}(c,w)
Grid{T<:Number}(c::AbstractArray{T},w::AbstractArray) = Grid{eltype(c),ndims(c),1,typeof(c),typeof(w)}(c,w)








Base.length(g::Grid) = length(g.c)
Base.eltype{T,N}(g::Grid{T,N}) = T
Base.ndims{T,N}(g::Grid{T,N}) = N
Base.size(g::Grid) = size(g.c)
Base.size(g::Grid, i::Integer) = size(g.c, i)
Base.getindex(g::Grid, i...) = g.c[i...]
function Base.setindex!(g::Grid, value, i...)
    g.c[i...] = value
end
Base.show(io::IO, g::Grid) = print(io, summary(g))
Base.writemime(io::IO, mime::MIME"text/plain", g::Grid) = writemime(io, mime, g.c)
Base.transpose(g::Grid) = transpose(g.c)
Base.ctranspose(g::Grid) = ctranspose(g.c)
Base.convert(T::Type{Array{Float64,1}},g::Grid) = convert(T,g.c)
Base.eachindex(g::Grid) = eachindex(g.c)



integrate(g::Grid, f::AbstractVector) = vecdot(g.w, f)
integrate(g::Grid, f::AbstractArray) = sum(g.w .* f, 1:ndims(g.w)) |> squeezeall
integrate(g::Grid, f) = integrate(g, f(g))



function integrate(g::Grid, f::AbstractArray, dims::Integer)
    @assert length(dims) == ndims(g.w)
    shp = fill(1,ndims(f))
    for (i,d) in enumerate(dims)
        shp[d] = size(g.w,i)
    end
    sum(reshape(g.w,tuple(shp...)) .* f, dims) |> squeezeall
end

Base.(:*){T}(A::Array{T,2}, g::Grid{Point3{T}}) = Grid(reshape(Point3{T}[A*x for x in g.c], size(g.c)), g.w)


# function expand(grid, basis, f)
#     S = basis.data' * (basis .* x.w)
#     (integrate(x, basis.*h)' / S)'
# end


# function expand(grid, basis, f)
#     S = Symmetric(LinAlg.BLAS.syrk('U', 'T', basis.data .* sqrt(grid.w)))
#     (integrate(grid, basis.*f)' / S)'
# end

# LinAlg.BLAS.herk('U', 'C', 1.0, mat(basis.data .* sqrt(g.w)))

# function qrgramschmidt(grid, basis)
#     Q,R = qr(basis .* sqrt(grid.w));
#     Q ./= sqrt(grid.w)
#     Q
# end



function lingrid(start::Real, stop::Real, n::Integer)
    x = collect(linspace(start, stop, n))
    dx = x[2] - x[1]
    Grid(x, fill(dx, length(x)))
end

function loggrid(start::Real, stop::Real, n::Integer)
    x = collect(linspace(start, stop, n))
    dx = x[2] - x[1]
    r = exp(x)
    Grid(r, r.*dx)
end


function simpsgrid(start::Real, stop::Real, n::Integer)
    @assert isodd(n) "Simpsons Rule requires an odd number of grid points"
    x = collect(linspace(start, stop, n))
    dx = x[2] - x[1]
    w = fill(dx/3, n)
    w[2:2:end] .*= 4
    w[3:2:end-1] .*= 2
    Grid(x, w)
end


function atomicradgrid(Nr::Int, rmax::Real)
    s = collect(1:Nr)
    r = rmax * log(1 - (s/(Nr+1.)).^2) ./ log(1 - (Nr/(Nr+1))^2)
    w = -2rmax*s ./ ( (1+Nr-s) .* (1+Nr+s) .* log(1 - (Nr/(Nr+1.))^2) )
    Grid(r, w)
end



lebgrid(n::Int) = lebedev[n]


function atomicgrid(r::Grid, Ω::Grid; R::Point3{Float64}=(0.,0.,0.))
    shp = (length(r),length(Ω))
    g = Grid(Array{Point3{Float64}}(shp), Array{Float64}(shp))
    for i in 1:length(Ω)
        θ,ϕ = Ω[i]
        for j in 1:length(r)
            g.w[j,i] = r[j]^2*r.w[j]*Ω.w[i]
            g.c[j,i] = sph2cart(r[j],θ,ϕ) + R
        end
    end
    g
end


# function multigrid(atoms)
#     for at in atoms
#         at.p = at.h(at.g)./sum([μ.h(at.g) for μ in atoms])
#     end
#     p = cat([ν.p for ν in atoms]...)
#     g = Grid(cat([μ.g for μ in atoms]...), p.*cat([μ.g.w for μ in atoms]...))
#     for at in atoms
#         at.rho = at.rho(g)
#         at.vh = at.vh(g)
#         at.v = at.v(g)
#     end
#     g
# end






