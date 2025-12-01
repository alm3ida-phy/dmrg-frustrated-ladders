# dmrg-frustrated-ladders
DMRG simulation of coupled frustrated spin-1/2 ladders using ITensors.jl. Includes calculation of energy, magnetization, and correlation matrices.


# DMRG Simulation of Coupled Frustrated Spin-1/2 Ladders

This repository contains a Julia script to simulate the ground state properties of a quantum magnetic system composed of two coupled, frustrated spin-1/2 ladders. The simulation uses the **Density Matrix Renormalization Group (DMRG)** algorithm implemented via the [ITensors.jl](https://github.com/ITensor/ITensors.jl) library.

## Physical Model

The code models a Hamiltonian for two parallel Heisenberg spin ladders coupled together. The system includes various interaction terms allowing for the study of magnetic frustration:

 **$Jx$**: Antiferromagnetic coupling along the legs (horizontal).
 **$Jy$**: Coupling along the rungs of each individual ladder (intra-ladder).
 **$Jcy$**: Coupling along the rungs connecting the two ladders (inter-ladder).
 **$Jd$**: Diagonal frustration within each ladder.
 **$Jcd$**: Diagonal frustration between the two ladders.

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
