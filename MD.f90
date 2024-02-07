! Author: Gianluca Regni
! Copyright (c) 2024 Gianluca Regni
! License: GPL-3.0 (https://www.gnu.org/licenses/)
! Credits: Please cite the author if you use or reference this program.
!
! Description: 
! This is a molecular dynamics program developed to simulate atom particles in a gas phase. Verlet integration
! (Velocity Verlet) was employed to numerically integrate the Newton’s equation of motion, allowing for
! the calculation of particle trajectories. The interaction between gas particles was simulated
! through a Lennard-Jones potential, where every atom experience a cumulative Lennard-Jones force generated
! by other particles.
! MD.f90 serves as the main program, while utilities.f90 comprises a set of subroutines. Ensure that they are
! located in the same folder before compilation.
!
! Run:
! 1. Create an input file named "input" with a structure similar to the provided samples in the input_files
!    folder (https://github.com/gianwski/molecular-dynamics/tree/main).
! 2. Compile and run the program using the following commands:
!    gfortran -c utilities.f90
!    gfortran utilities.o MD.f90 -o MD
!    ./MD
! 3. Two outputs will be generated:
!    - "md.xyz": a file that contain the coordinates of the system at every step (trajectory). 
!    - "energy.dat": a file that contain energies data of the system. 

  program md

  use utilities

  implicit none

! Parameter for unit conversions to atomic unit

  double precision, parameter  :: kcal_to_Hartree = 0.0015936015d0
  double precision, parameter  :: ev_to_Hartree = 0.0367492929d0
  double precision, parameter  :: ang_to_bohr = 1.8897261d0
  double precision, parameter  :: fs_to_timeau = 41.34137314d0
  double precision, parameter  :: amu_to_au = 1822.8884850d0

! double precision, parameter  :: Hartree_to_kcal = 627.50947d0
! double precision, parameter  :: Hartree_to_ev = 27.211407953d0
! double precision, parameter  :: bohr_to_ang = 0.52917721d0
  double precision, parameter  :: timeau_to_fs = 0.02418884254d0
! double precision, parameter  :: au_to_amu = 0.00054857990943d0

! Definition of variables

  integer                                        :: i, j, k, l, step, nstep, natom
  double precision                               :: timestep, mass, eps, sig, kin_ener, pot_ener, energy, time
  double precision, allocatable, dimension(:,:)  :: positions, velocities, distances, forces
  character(len=100)                             :: filename = "input" ! Input file name
  character(len=2), allocatable, dimension(:)    :: atom


! Open output file for trajectory and energies

  open(10, file="md.xyz")
  open(11, file="energy.dat")

! Read input variables

  call lecture(filename, positions, atom, velocities, nstep, mass, timestep, sig, eps, natom)

! Allocate memory

  allocate(forces(natom,3), distances(natom,natom))

! Write input parameters
  
  write(6,"(A)") "Initial parameters"
  write(6,*)
  write(6,"(A,I4)") "number of steps : ", nstep
  write(6,"(A,F6.3)") "time step (fs)  : ", timestep
  write(6,"(A,I2)") "number of atoms : ", natom
  write(6,"(A,F9.6)") "masses (amu)    : ", mass
  write(6,"(A,F9.6)") "epsilon (eV)   : ", eps 
  write(6,"(A,F9.6)") "sigma (ang)     : ", sig
  write(6,*)

! Write geometry information

  write(6,"(A)") "Geometry (ang)"
  write(6,"(8X,3(A1,10X))") "x", "y", "z"
  do i=1,natom
    write(6,"(A2,1X,3(F10.6,1X))") atom(i), positions(i,:)
  enddo
  write(6,*)

! Write velocities information

  write(6,"(A)") "Velocities (ang/fs)"
  write(6,"(8X,3(A1,10X))") "x", "y", "z"
  do i=1,natom
    write(6,"(A2,1X,3(F10.6,1X))") atom(i), velocities(i,:)
  enddo

! Conversion to au

  mass = mass * amu_to_au
  timestep = timestep * fs_to_timeau
  sig = sig * ang_to_bohr
  eps = eps * ev_to_Hartree
  positions = positions * ang_to_bohr
  velocities = velocities * ang_to_bohr / fs_to_timeau

! Write input converted

  write(6,*)
  write(6,"(A)") "-------------------------------------"
  write(6,"(A)") "Initial parameters in atomic units"
  write(6,*) " "
  write(6,"(A,F10.6)") "time step  : ", timestep
  write(6,"(A,F12.6)") "masses     : ", mass
  write(6,"(A,E13.6)") "epsilon    : ", eps 
  write(6,"(A,F9.6)") "sigma      : ", sig
  write(6,*)

! Write geometry information

  write(6,"(A)") "Geometry"
  write(6,"(8X,3(A1,10X))") "x", "y", "z"
  do i=1,natom
    write(6,"(A2,1X,3(F10.6,1X))") atom(i), positions(i,:)
  enddo
  write(6,*)

! Write velocities information

  write(6,"(A)") "Velocities"
  write(6,"(8X,3(A1,10X))") "x", "y", "z"
  do i=1,natom
    write(6,"(A2,1X,3(F10.6,1X))") atom(i), velocities(i,:)
  enddo

! Compute distances for LJ

  call get_distances(positions, distances, natom)

! Compute LJ forces

  call get_forces(forces, distances, positions, sig, eps, natom)

! Calculate total energy

  call get_kinetic(velocities, natom, kin_ener, mass)
  call get_potential(distances, natom, sig, eps, pot_ener)
  energy = kin_ener + pot_ener
  
! Initialization of step

  step = 0

! Write geometry and energy on the output

  time = step * timestep * timeau_to_fs
  write(10,"(I3)") natom
  write(10,"(A,I4,A,F6.3,3(A,E15.6))") "# step: ", step, ",  t(fs): ", time, &
                                     & ",   E =", energy, ",   Ek =", kin_ener, ",   Ep =", pot_ener
  do i = 1,natom
     write(10,*) atom(i), positions(i,:)
  enddo
  write(11,*) "# step, time, energy, kin_ener, pot_ener"
  write(11,*) step, time, energy, kin_ener, pot_ener

! Debug
! write(6,*) 
! write(6,*) "step, Total Energy, Kinetic, Potential"
! write(6,*)
! write(6,*) step, ") ", ",   E =", energy, ",   Ek =", kin_ener, ",   Ep =", pot_ener

! *****************************************************

! Loop for the trajectory

  do step = 1, nstep

! Velocity Verlet

  do i=1,natom
    do j=1,3
      velocities(i,j)=velocities(i,j)+0.5d0*timestep*forces(i,j)/mass
    enddo
  enddo

  do i=1,natom
    do j=1,3
      positions(i,j)=positions(i,j)+timestep*velocities(i,j)
    enddo
  enddo

  call get_distances(positions,distances,natom)
  call get_forces(forces,distances,positions,sig,eps,natom)

  do i=1,natom
    do j=1,3
      velocities(i,j)=velocities(i,j)+0.5d0*timestep*forces(i,j)/mass
    enddo
  enddo

! Calculate energies

  call get_kinetic(velocities, natom, kin_ener, mass)
  call get_potential(distances, natom, sig, eps, pot_ener)
  energy = kin_ener + pot_ener

! Debug
! write(6,*) step, ") ", ",   E =", energy, ",   Ek =", kin_ener, ",   Ep =", pot_ener

! Write geometry and energy on the output

  time = step * timestep * timeau_to_fs
  write(10,"(I3)") natom                                
  write(10,"(A,I4,A,F9.3,3(A,E15.6))") "# step: ", step, ",  t(fs): ", time, &
                                    & ",   E =", energy, ",   Ek =", kin_ener, ",   Ep =", pot_ener
  do i = 1, natom
     write(10,*) atom(i), positions(i,:)
  enddo

  write(11,*) step, time, energy, kin_ener, pot_ener
  
  enddo

! *****************************************************

  close(10)
  close(11)

  end program md