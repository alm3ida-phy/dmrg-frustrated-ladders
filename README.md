# dmrg-frustrated-ladders
DMRG simulation of coupled frustrated spin-1/2 ladders using ITensors.jl. Includes calculation of energy, magnetization, and correlation matrices.


# DMRG Simulation of Coupled Frustrated Spin-1/2 Ladders

This repository contains a Julia script to simulate the ground state properties of a quantum magnetic system composed of two coupled, frustrated spin-1/2 ladders. The simulation uses the **Density Matrix Renormalization Group (DMRG)** algorithm implemented via the [ITensors.jl](https://github.com/ITensor/ITensors.jl) library.

## Physical Model

The code models a Hamiltonian for two parallel Heisenberg spin ladders coupled together. The system includes various interaction terms allowing for the study of magnetic frustration:

 **$J_x$**: Antiferromagnetic coupling along the legs (horizontal).
 **$J_y$**: Coupling along the rungs of each individual ladder (intra-ladder).
 **$J_{cy}$**: Coupling along the rungs connecting the two ladders (inter-ladder).
 **$J_{d}$**: Diagonal frustration within each ladder.
 **$J_{cd}$**: Diagonal frustration between the two ladders.

The total system size is defined by $4 \times L_x$, where $L_x$ is the length of the ladders.

## Dependencies

To run this code, you need [Julia](https://julialang.org/) installed along with the following packages:

* `ITensors`
* `ITensorMPS`
* `HDF5`
* `JLD`
* `MKL` (Optional, but recommended for performance on Intel hardware)

You can install these dependencies in the Julia REPL:

```julia
using Pkg
Pkg.add(["ITensors", "ITensorMPS", "HDF5", "JLD", "MKL"])

## Usage

The parameters for the simulation (couplings $J$, system size $L_x$​, bond dimensions, etc.) are defined within the script variables. To change the physics of the simulation, edit the variables at the bottom of the .jl file

## Running the Simulation

To execute the code, use the Julia command line. It is highly recommended to use the -t flag to specify the number of threads. The code is optimized to use ITensors' threaded block-sparse operations.

## Command syntax:

julia -t [number_of_cores] [script_name].jl


