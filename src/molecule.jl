


doc"""
`Molecule(; atoms=[], params...)`

creates an molecules with `atoms`.
"""
type Molecule <: Associative{Symbol,Any}
    name  ::AbstractString
    atoms ::Atoms                # atoms
    basis ::Basis                # basis
    g     ::MultiGrid            # integration grid
    d     ::Dict{Symbol,Any}     # additional properties
end


Base.delete!(m::Molecule, key)                                    = delete!(m.d, key)
Base.done(m::Molecule, i)                                         = done(m.d, i)
Base.empty!(m::Molecule)                                          = empty!(m.d)
Base.copy(m::Molecule)                                            = Atom(at.Z,at.R,copy(m.d))
Base.get(m::Molecule, key, default)                               = get(m.d, key, default)
Base.get(default::Union{DataType,Function}, m::Molecule, key)     = get(default, m.d, key)
Base.get!(m::Molecule, key0, default)                             = get!(m.d, key0, default)
Base.get!(default::Union{DataType,Function}, m::Molecule, key0)   = get!(default, m.d, key0)
Base.getindex(m::Molecule, key)                                   = getindex(m.d, key)
Base.getkey(m::Molecule, key, default)                            = getkey(m.d, key, default)
Base.haskey(m::Molecule, key)                                     = haskey(m.d, key)
Base.isempty(m::Molecule)                                         = isempty(m.d)
Base.length(m::Molecule)                                          = length(m.d)
Base.next(m::Molecule, i)                                         = next(m.d, i)
Base.pop!(m::Molecule, key)                                       = pop!(m.d, key)
Base.pop!(m::Molecule, key, default)                              = pop!(m.d, key, default)
Base.setindex!(m::Molecule, v0, key0)                             = setindex!(m.d, v0, key0)
Base.similar(m::Molecule)                                         = similar(m.d)
Base.sizehint!(m::Molecule, newsz)                                = sizehint!(m.d, newsz)
Base.start(m::Molecule)                                           = start(m.d)
Base.merge!(m::Molecule, d::Dict)                                 = begin merge!(m.d, d); d end
Base.merge(m::Molecule, d::Dict)                                  = Atom(merge(m.d, d))
Base.writemime(io::IO, mime::MIME"text/plain", m::Molecule)       = show_symbol_dict(io,m)
Base.summary(m::Molecule)                                         = string(typeof(m),"($(m.name))")
Base.show(io::IO, m::Molecule)                                    = print(io,summary(at))


function Base.call(::Type{Molecule};
                    name  ::AbstractString = "",
                    atoms ::Atoms          = Atom[],
                    params...)
    m = Molecule(name,atoms,NAOBasis(),MultiGrid(),Dict{Symbol,Any}(params))

    if !isempty(atoms)
        for at in atoms
            if !hasbasis(at)
                dft!(at)
                splineradial!(at)
                naobasis!(at)
            end

            if !hasgrid(at)
                makegrid!(at)
            end
        end

        delleypartition!(atoms)
        makegrid!(m)
        makebasis!(m)
    end

    m
end

function makegrid!(m::Molecule)
    m.g = multigrid(m.atoms)
end


function makebasis!(m::Molecule)
    m.basis = vcat([at[:basis] for at in m.atoms]...)
end

print_stats(m::Molecule) = print_stats(m[:stats])




function bondslice(atoms, ind)
    u = [atoms[ind[i+1]].R - atoms[ind[i]].R for i in 1:length(ind)-1]
    d = [norm(r) for r in u]
    u ./= d
    d = cumsum(d)
    function r(λ)
        i = findfirst(d .> λ)
        if i == 0; i = length(u) end
        atoms[ind[i+1]].R + (λ - d[i])*u[i]
    end
    d, f -> (λ -> f(r(λ))[])
end
        
bondslice(atoms) = bondslice(atoms, 1:length(atoms))



