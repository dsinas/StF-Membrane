# Characteristic length of membrane composition

<br>This repository contains source codes for simulation of the 2D lattice-based model of the membrane that comprises a binary lipid phase mixtures and incorporates out-of-plane fluctuations (Monge representation).<br>

## Introduction

The program implements Monte Carlo simulations of the lipid-lipid interaction in the membrane and update the height fluctuation modes probabilistically in Fourier space. The program calculates the structure factor and the characteristic length of the membrane composition. The model considers a coupling between the local membrane composition and the local curvature so that with no coupling the model shall retrieve the 2D Ising universality class.<br>

## Requirements

This program uses fftw3 library for fast Fourier transform.<br>

## Compilation

Compile the program with **make**:

```bash
make clean
make
```

Run the program for different system sizes, temperature as well as coupling to curvature strength by editing ``run.sh`` and run:

```bash
./run.sh
```

## Configuration

### Parameters:
The system size, the temperature (Ising coupling), and the composition-curvature coupling are given by default and can be hardcoded as system setting in the **main.cpp**. They may alternatively be set through input arguments by flags ``-L``, ``-J`` and ``-g``, e.g.:

```bash
./stfMem -L 100 -J 0.5 -g 1.00
```

### Hyperparameters:
The number of initial iterations to get equilibrium and the total number of sweeps for MC moves are given experimentaly. They can be set according to the desired accuracy and the available resources.<br>

## Output Data

* **st-factor\*.dat**: The structure factor as a function of mode numbers.

* **cmps\*.dat**: The snapshot of the system composition.

* **height\*.dat**: The height profile.

* **lh\*.xyz**: The composition + height data as **xyz** file for 3D plot.

* **sys_data\*.dat**: All the physical system parameters.
