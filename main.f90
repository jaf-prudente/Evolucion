
  program main

  implicit none

  common / gauss   / amp1,amp2,r0,sigma1,sigma2
  common / Vparam  / lambda,gamma,epsilon,Vtype
  common / gauss3  / r3,amp3,sigma3,kr
  common / Vparam3 / lambda3,epsilon3,Vtype3

  integer Nr,Nt,Noutput0D,Noutput1D
  integer Vtype,Vtype3
  integer mattertype

  real(8) dr,dt,dtfac
  real(8) amp1,amp2,r0,sigma1,sigma2
  real(8) r3,amp3,sigma3,kr
  real(8) lambda,gamma,epsilon
  real(8) lambda3,epsilon3

  mattertype=0

! The variables introduced above are:
!
! dr:         Spatial interval.
! dt:         Time step.
!
! Nr:         Total number of grid points.
! Nt:         Total number of time steps.
!
! Noutput0D:  How often to do 0D output (in time steps).
! Noutput1D:  How often to do 1D output (in time steps).
!
! amp1:        Amplitude1 of initial gaussian profile for phi1.
! amp2:        Amplitude2 of initial gaussian profile for phi2.
! r0:          Center of initial gaussian profile for phi1 & phi2.
! sigma1:      Width1 of initial gaussian profile for phi1.
! sigma2:      Width2 of initial gaussian profile for phi2.
!
! Vtype1:     Type of potential:
!               0. No potential.
!               1. Quadratic (only mass term)
!               2. Quadratic plus quartic
!
! lambda:     Potential parameter for V.
! gamma:      Potential parameter for V.
! epsilon:    Potential parameter for V.


! --->   MESSAGES   <---

  print *
  print *, '---------------------------------------'
  print *, '>                                     <'
  print *, '>             PROGRAM Sphere          <'
  print *, '>                                     <'
  print *, '>                                     <'
  print *, '>          Evolving a spacetime       <'
  print *, '>         with spherical symmetry	   <'
  print *, '>                                     <'
  print *, '---------------------------------------'

  print *


! GET PARAMETERS

! Matter type

  print *, 'Que tipo de materia quieres evolucionar?'
  read(*,*) mattertype
  print *

! Grid spacing.

  print *, 'Size of grid spacing (dr)'
  read(*,*) dr
  print *

! Courant parameter (dtfac).

  print *, 'Size of Courant parameter (dtfac=dt/dr)'
  read(*,*) dtfac
  print *

  dt = dtfac*dr

! Number of grid points (Nr).

  print *, 'Number of grid points (Nr)'
  read(*,*) Nr
  print *

! Number of time steps (Nt).

  print *, 'Number of time steps (Nt)'
  read(*,*) Nt
  print *

! 0D output.

  print *, 'How often do you want 0D output? (Noutput0D)'
  read(*,*) Noutput0D
  print *

  if (Noutput0D>Nt) Noutput0D=Nt

! 1D output.

  print *, 'How often do you want 1D output? (Noutput1D)'
  read(*,*) Noutput1D
  print *

  if (Noutput1D>Nt) Noutput1D=Nt

! Parameters for the initial gaussian.

  r0 = 0.0D0

  print *, 'Amplitude1 of initial gaussian profile for phi1 (amp)'
  read(*,*) amp1
  print *

  print *, 'Amplitude2 of initial gaussian profile for phi2 (amp)'
  read(*,*) amp2
  print *

  print *, 'Width1 of initial gausian profile for phi1 (sigma1)'
  read(*,*) sigma1
  print *

  print *, 'Width2 of initial gausian profile for phi2 (sigma2)'
  read(*,*) sigma2
  print *

! Ask for the type of potential for phi.

  print *, 'Type of potential for phi1:'
  print *, '   0. No potential'
  print *, '   1. Quadratic (mass term)'
  print *, '   2. Quartic (lambda)'
  read(*,*) Vtype
  print *

! Ask for parameters for potential.

  if (Vtype.eq.1) then
     epsilon = 0.0D0
     print *, 'Coefficient of quadratic term'
     read(*,*) lambda
     print *
  else if (Vtype.eq.2) then
     print *, 'Coefficient of quadratic term'
     read(*,*) lambda
     print *
     print *, 'Coefficient of quartic term'
     read(*,*) epsilon
     print *
 
  end if

! ----------------------------------

! --->   EVOLUTION   <---

  call evolve(Nr,Nt,Noutput0D,Noutput1D,dr,dt,mattertype)


! --->   END   <---

  print *
  print *, 'Program has finished'
  print *

  end program main

