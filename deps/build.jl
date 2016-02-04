


using Element: generate_lebedev

# generate lebedev angular grids
generate_lebedev()






using BinDeps
using Compat

@BinDeps.setup

libxc = library_dependency("libxc")

@osx_only begin
	using Homebrew
	provides(Homebrew.HB, "libxc", libxc, os = :Darwin)
end

provides(AptGet, "libxc", libxc)
provides(Yum, "libxc-devel", libxc)

julia_usrdir = normpath(JULIA_HOME*"/../") # This is a stopgap, we need a better built-in solution to get the included libraries
libdirs = AbstractString["$(julia_usrdir)/lib"]
includedirs = AbstractString["$(julia_usrdir)/include"]

provides(
	Sources,
	URI("http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.2.2.tar.gz"),
	libxc
)

provides(
	BuildProcess,
	Autotools(lib_dirs = libdirs, include_dirs = includedirs),
	libxc
)



@compat @BinDeps.install Dict(:libxc => :libxc)









