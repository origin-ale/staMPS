# Program taken from the ITensor DMRG tutorial at https://docs.itensor.org/ITensorMPS/stable/tutorials/DMRG.html
# in order to get familiar with the package

using ITensors, ITensorMPS # Import-like statement

let # Create a local scope; need this to keep Julia from freaking out about global vars\
    N = 100 # 100-site DMRG
    sites = siteinds("S=1", N) # Array of N Index objects w/ properties (dimension and tag) of S=1 spins

    os = OpSum() # Accumulate Hamiltonian terms
    for j=1:N-1 # Range statement (Julia arrays start at 1)
        os += "Sz", j, "Sz", j+1 # Each term is a comma-sep list of factors (?)
        os += 1/2, "S+", j, "S-", j+1
        os += 1/2, "S-", j, "S+", j+1
    end
    H = MPO(os, sites) # Construct Hamiltonian w/physical indices given by array "sites"

    psi0 = random_mps(sites; linkdims=10) # Random MPS w/ physical indices given by array "sites" and bond dim 10

    nsweeps = 5 # Do 5 DMRG sweeps
    maxdim = [10,20,100,100,200] # Maximum bond dims for each sweep
    cutoff = [1E-10] # Truncation error goal for each (ie. all) sweep
    
    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff) # Call ITensor DMRG starting from psi0; automatically prints at each iteration
    return
end