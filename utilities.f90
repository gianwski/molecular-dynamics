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
  
  module utilities

  public::lecture_xyz             ! Read and save positions in an .xyz file
  public::get_distances           ! Calculate square of distances between atoms
  public::get_forces              ! Calculate forces from the positions
  public::get_kinetic             ! Calculate kinetic energy
  public::get_potential           ! Calculate potential energy

  contains

! **********************************************************

  subroutine lecture (filename, positions, atom, velocities, nstep, mass, timestep, sig, eps, natom)

! Record input parameters

  implicit none

  integer                                        :: i, j, nstep, natom
  double precision                               :: timestep, mass, eps, sig
  double precision, allocatable, dimension(:,:)  :: positions, velocities
  character(len=100)                             :: filename
  character(len=2), allocatable, dimension(:)    :: atom

  open(1,file=filename,status="old")

  read(1,*) 
  read(1,*) nstep
  
  read(1,*)
  read(1,*) timestep

  read(1,*)
  read(1,*) natom

  ! Allocate arrays

  allocate(positions(natom,3))
  allocate(atom(natom))
  allocate(velocities(natom,3))

  read(1,*)
  read(1,*) mass
  
  read(1,*)
  do i=1,natom
     read(1,*) atom(i), positions(i,:)
  enddo

  read(1,*)
  read(1,*) sig

  read(1,*)
  read(1,*) eps

  read(1,*)
  do i=1,natom
     read(1,*) velocities(i,:)
  enddo

  ! Close file
  close(1)

  endsubroutine lecture

! **********************************************************

  subroutine get_distances (positions, distances, natom)

! Calculate distances between atoms

  implicit none

  integer, intent(in)           :: natom
  double precision, intent(in)  :: positions(natom,3)
  double precision, intent(out) :: distances(natom,natom)

  integer                       :: i, j

!  Initialization of the distances
  distances=0.0d0

  do i = 1, natom
      do j = 1, natom
        distances(i,j) = dsqrt((positions(i,1) - positions(j,1))**2 + &
                             & (positions(i,2) - positions(j,2))**2 + &
                             & (positions(i,3) - positions(j,3))**2)
      enddo
  enddo

  endsubroutine get_distances

! **********************************************************

  subroutine get_forces(forces, distances, positions, sig, eps, natom)

! Calculate LJ forces

  implicit none

  integer, intent(in)             :: natom
  double precision, intent(in)    :: eps, sig
  double precision, intent(inout) :: forces(natom,3)
  double precision, intent(in)    :: positions(natom,3), distances(natom,natom)

  integer                         :: i, j, k


! Initialization of the forces
  forces = 0.d0

  do i = 1, natom  ! the atom experiencing the forces
    do k = 1,3  ! xyz
        do j = 1, natom  ! the atom generating the force
          if (i == j) cycle  ! skip the interaction with itself
          forces(i,k) = forces(i,k) + 48*eps*(positions(i,k)-positions(j,k)) * &
                                    & (sig**12/distances(i,j)**14 - 0.5d0*sig**6/distances(i,j)**8) 
        enddo
    enddo
  enddo

  endsubroutine get_forces

! **********************************************************

  subroutine get_kinetic(velocities,natom,kin_ener,mass)

! Calculate Kinetic energy

  implicit none

  double precision, intent(in)     :: velocities(natom,3), mass
  double precision, intent(out)    :: kin_ener
  integer, intent(in)              :: natom
  integer                          :: i,k

! Initialization of kinetic energy
  kin_ener=0.0d0

  do i=1,natom  ! iterate over all atoms
    do k=1,3  ! iterate over x, y and z
      kin_ener = kin_ener + 0.5d0 * mass * velocities(i,k)*velocities(i,k)
    enddo
  enddo

  endsubroutine get_kinetic

! **********************************************************

  subroutine get_potential (distances,natom,sig,eps,pot_ener)

! Calculate Potential energy

  implicit none

  double precision, intent(in)  :: distances(natom,natom), sig, eps
  double precision, intent(out) :: pot_ener
  integer, intent(in)           :: natom
  integer                       :: i,j

! Initialization of potential energy
  pot_ener = 0.0d0

! Indexes are set up to avoid calculating double interactions (between the same pair of atoms)

  do i = 1,natom-1
     do j = i+1,natom
        pot_ener = pot_ener + 4.d0 * eps * ((sig/distances(i,j))**12-(sig/distances(i,j))**6)
     enddo
  enddo

  endsubroutine get_potential

end module utilities
