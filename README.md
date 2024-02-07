# Gas Phase Molecular Dynamics Simulator using Verlet Integration

A molecular dynamics program was developed to simulate atom particles in a gas phase for a university exam. Verlet
integration (Velocity Verlet) was employed to numerically integrate the Newton’s equation of motion, allowing for the calculation
of particle trajectories. The interaction between gas particles was simulated through a Lennard-Jones potential,
where every atom experience a cumulative Lennard-Jones force generated by other particles.
The program, written in FORTRAN 90, is a combination of two distinct codes. The first code MD.f90 serves
as the main program, while the second utilities.f90 comprises a set of subroutines that are invoked within the main
program.

## Folder Structure

- **input_files**: Includes sample input files. Make sure that the input file is named `input` before running the program.

- **output**: Stores all output files generated by the program for the Hydrogen and Argon gas-phase system, including energy plots and the Gnuplot code used to generate them.

## Instructions

1. Create an input file named `input` with a structure similar to the provided samples in the `input_files` folder.
2. Compile and run the program using the following commands:

```bash
gfortran -c utilities.f90
gfortran utilities.o MD.f90 -o MD
./MD
```
3. Two outputs will be generated:
   - `md.xyz`: a file that contain the coordinates of the system at every step (trajectory). 
   - `energy.dat`: a file that contain energies data of the system. 


