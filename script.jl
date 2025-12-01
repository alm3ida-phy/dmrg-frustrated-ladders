using ITensors
using Random
using LinearAlgebra
using JLD	
using MKL
using ITensorMPS
using HDF5

# --- Performance Configuration (Block-Sparse Optimization) ---
# IMPORTANT: Only use this configuration if `conserve_qns=true`.
# It optimizes block-parallelism by disabling BLAS threading to avoid oversubscription.
#
# If running Dense DMRG (`conserve_qns=false`), comment these lines out 
# to allow BLAS to use all available threads.
ITensors.Strided.disable_threads()
BLAS.set_num_threads(1)
ITensors.enable_threaded_blocksparse()

"""
    dot(num::Float64)
Helper function to convert floating point numbers to strings 
and remove the decimal point for clean filenames (e.g., 0.5 -> "05").
"""
function dot(num::Float64)
    str = replace(string(num), "." => "")
    return str
end

"""
    heisenberg(...)
Main simulation function for the coupled frustrated ladders.
"""
function heisenberg(Jx,Jy,Jcy,Jcd,Jd,Lx,spin,n_sweeps,n_maxdim)
    # Define physical sites (Spin 1/2). 
    # The total number of sites is 4 * Lx (4 sites per rung column).
    # conserve_qns=true utilizes U(1) symmetry (conservation of Sz).
    sites = siteinds("S=1/2", 4*Lx; conserve_qns=true)
    
    os = OpSum()

    # --- Hamiltonian Construction ---
    
    # 1. Jx: Interaction along the legs (Horizontal)
    # Connects site j to site j+4 (the site in the same leg in the next column).
    for j=1:4*(Lx-1)
        os += Jx,"Sz",j,"Sz",j+4
        os += Jx/2,"S+",j,"S-",j+4
        os += Jx/2,"S-",j,"S+",j+4
    end

    # 2. Jy: Rungs within each ladder (Vertical, Intra-ladder)
    # Iterates with step 2: (1-2), (3-4), etc.
    # Connects the bottom leg to the top leg within the same ladder.
    for j=1:2:4*Lx
        os += Jy,"Sz",j,"Sz",j+1
        os += Jy/2,"S+",j,"S-",j+1
        os += Jy/2,"S-",j,"S+",j+1
    end

    # 3. Jcy: Coupling between the two ladders (Vertical, Inter-ladder)
    # Iterates with step 4 starting at 2: (2-3), (6-7), etc.
    # Connects the top leg of the bottom ladder (even site) 
    # to the bottom leg of the top ladder (odd site).
    for j =2:4:4*Lx
        os += Jcy,"Sz",j,"Sz",j+1
        os += Jcy/2,"S+",j,"S-",j+1
        os += Jcy/2,"S-",j,"S+",j+1
    end

    # 4. Jcd: Diagonal frustration between ladders (Inter-ladder)
    # Part A: Connects top of bottom ladder (j) to bottom of top ladder next col (j+5).
    for j = 2:4:4*(Lx-1)
        os += Jcd,"Sz",j,"Sz",j+5
        os += Jcd/2,"S+",j,"S-",j+5
        os += Jcd/2,"S-",j,"S+",j+5
    end
    # Part B: Connects bottom of top ladder (j) to top of bottom ladder next col (j+3).
    for j=3:4:4*(Lx-1)
        os += Jcd,"Sz",j,"Sz",j+3
        os += Jcd/2,"S+",j,"S-",j+3
        os += Jcd/2,"S-",j,"S+",j+3
    end

    # 5. Jd: Diagonal frustration within ladders (Intra-ladder)
    # Part A: Diagonal 'up' within a ladder (e.g., site 1 to 6).
    for j=1:2:4*(Lx-1)
        os += Jd,"Sz",j,"Sz",j+5
        os += Jd/2,"S+",j,"S-",j+5
        os += Jd/2,"S-",j,"S+",j+5
    end
    # Part B: Diagonal 'down' within a ladder (e.g., site 2 to 5).
    for j =2:2:4*(Lx-1)
        os += Jd,"Sz",j,"Sz",j+3
        os += Jd/2,"S+",j,"S-",j+3
        os += Jd/2,"S-",j,"S+",j+3
    end
 
    # Convert the operator sum to a Matrix Product Operator (MPO)
    H = MPO(os,sites)
        
    # --- State Initialization ---
    # Create an initial product state with a specific Sz sector.
    # If spin=0, we have an equal number of Up and Dn (half-filling).
    # 'spin' acts as an offset to the number of Up spins.
    state = [n <= spin + 2*Lx ? "Up" : "Dn" for n=1:4*Lx]
    
    # Shuffle the state to randomize the distribution of Up/Dn spins.
    # This helps avoid getting stuck in local minima during DMRG.
    for i=1:10
        shuffle!(state)
    end
    
    psi0 = productMPS(sites,state)

    # --- DMRG Sweep Settings ---
    sweeps = Sweeps(n_sweeps)
    # Set the maximum bond dimension allowed for the MPS
    setmaxdim!(sweeps,maxdim)
    # Set the minimum bond dimension to ensure enough states are kept early on
    setmindim!(sweeps,16,64,128,256,300,420,500,650,700,760,830,900,1000,1000)
    # Truncation error cutoff
    setcutoff!(sweeps,1E-7)
    # Add noise to the density matrix to prevent getting stuck in local minima.
    # Noise is gradually reduced to 0 in later sweeps.
    noise!(sweeps,1E-5,1E-5,1E-5,1E-5,1E-5,1E-8,1E-8,1E-8,1E-8,1E-8,1E-10,1E-10,0)
    
    # Run the DMRG algorithm
    energy, psi = dmrg(H,psi0, sweeps)
    println(" energy = $energy")
    return energy, psi
end

# --- Main Execution Block ---

# System parameters
Lx = 32         # Length of the ladder (number of columns)
sweeps = 20     # Number of DMRG sweeps
maxdim = 3000   # Maximum bond dimension
sz = 0          # Magnetization sector (offset from half-filling)

# Coupling constants
Jx = 0.8        # Leg coupling
Jd = 0.7        # Intra-ladder diagonal
Jcy = 0.5       # Inter-ladder rung
Jy = 1          # Intra-ladder rung
Jcd = 1         # Inter-ladder diagonal

# Run the simulation
energy, psi = heisenberg(Jx,Jy,Jcy,Jcd,Jd,Lx,sz,sweeps,maxdim)

# Check quantum number conservation (Total Sz)
println(flux(psi))

# --- Data Saving ---

##### If you wish to save the wave function (HDF5 format) ######
#f = h5open("psi_$(Lx)x4_spin$(sz)_$(dot(Jx))jx_$(dot(Jd))jd_$(dot(Jcd))jcd_$(dot(Jcy))jcy.h5","w")
#write(f,"psi",psi)
#close(f)

# Compute observables
avg_sz = expect(psi,"Sz")               # Local magnetization
czz = correlation_matrix(psi,"Sz","Sz") # Spin-spin correlation (Longitudinal)
cpm = correlation_matrix(psi,"S+","S-") # Transverse correlation (+-)
cmp = correlation_matrix(psi,"S-","S+") # Transverse correlation (-+)

# Save observables to JLD files
save("energy_$(Lx)x4_spin$(sz)_$(dot(Jy))jy_$(dot(Jx))jx_$(dot(Jd))jd_$(dot(Jcy))jcy_$(dot(Jcd))jcd.jld","data",energy)
save("sz_$(Lx)x4_spin$(sz)_$(dot(Jy))jy_$(dot(Jx))jx_$(dot(Jd))jd_$(dot(Jcy))jcy_$(dot(Jcd))jcd.jld","data",avg_sz)
save("czz_$(Lx)x4_spin$(sz)_$(dot(Jy))jy_$(dot(Jx))jx_$(dot(Jd))jd_$(dot(Jcy))jcy_$(dot(Jcd))jcd.jld","data",czz)
save("cpm_$(Lx)x4_spin$(sz)_$(dot(Jy))jy_$(dot(Jx))jx_$(dot(Jd))jd_$(dot(Jcy))jcy_$(dot(Jcd))jcd.jld","data",cpm)
save("cmp_$(Lx)x4_spin$(sz)_$(dot(Jy))jy_$(dot(Jx))jx_$(dot(Jd))jd_$(dot(Jcy))jcy_$(dot(Jcd))jcd.jld","data",cmp)

