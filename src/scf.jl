

type SCF
    Etol::Float64
    rhotol::Float64
    maxsteps::Int  # max num of iterations
    n::Int  # how many previous iterations to store
    grid  # integration grid
    i::Int  # iteration count
    converged::Bool
    E::Array{Float64}  #
    rho::Array{AbstractArray}  #
    dE::Float64
    drho::Float64
end

SCF(;Etol=1e-6, rhotol=1e-8, maxsteps=100, n::Int=1, grid=nothing) = SCF(Etol,rhotol,maxsteps,n,grid,1,false,[],[],NaN,NaN)



function Base.call(scf::SCF, E::Float64, rho::AbstractArray)
    if length(scf.E) > 0
        scf.dE = abs(E - scf.E[end])
        scf.drho = integrate(scf.grid, (rho - scf.rho[end]).^2)
    end

    if scf.dE < scf.Etol && scf.drho < scf.rhotol
        scf.converged = true
    end

    scf.i += 1

    if scf.i > scf.maxsteps
        warn("max # of iterations exceeded")
        scf.converged = true
    end
    
    push!(scf.rho, rho)
    push!(scf.E, E)
    
    if scf.n > 0 && length(scf.rho) > scf.n
        shift!(scf.rho)
        shift!(scf.E)
    end
end




simple_logger = scf -> !isnan(scf.dE) ? @printf("|δE| = %8.2e     |δρ| = %8.2e\n", scf.dE, scf.drho) : return




