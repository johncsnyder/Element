__precompile__()

module Element

export Atom, Molecule, dft, dft!, lingrid, loggrid, lebgrid, atomicradgrid, atomicgrid,
	AtomicGrid, AtomicRadialGrid, AngularGrid, makebasis!, makegrid!, eachatom, NAO, 
	MultipoleFunc, MultipoleFuncSpl, RadialFunc, MultipoleBasis, integrate, ∫, grad, ∇, 
	mat, expandmultipole, Point3, Point3f, Point2, dist, randp, rotmat, randrotmat, 
	rotate!, setposition, simple_logger, bondslice, print_stats

include("base.jl")
include("multistep.jl")
include("fzero.jl")
include("findsteps.jl")
include("point.jl")
include("func.jl")
include("spline.jl")
include("sphharm.jl")
include("sphharmtable.jl")
include("stats.jl")
# include("raggedarray.jl")
include("lebedev.jl")
include("coord.jl")  		# point, func
include("grid.jl")  		# point
include("atomicgrid.jl")
include("radialfunc.jl")  	# point, func
include("multipole.jl")  	# point, func, grid
include("coulomb.jl")
include("basis.jl")  		# point, func
include("eig.jl")
include("radschrod.jl")
include("radpoisson.jl")
include("xc.jl")
include("mixer.jl")
include("scf.jl")
include("atomicdata.jl")
include("atom.jl")
include("atomicdft.jl")
include("nao.jl")
include("naobasis.jl")
include("partition.jl")
include("multigrid.jl")
include("molecule.jl")
include("hamiltonian.jl")
include("poissonmultipole.jl")
include("dft.jl")
include("rotation.jl")
# include("ml.jl")


function __init__()
	global lebedev = open(joinpath(Pkg.dir("Element"), "deps/lebedev.jls")) do io deserialize(io) end
end

end  # module Element
