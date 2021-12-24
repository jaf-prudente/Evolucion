
  module arrays

! --->  MODULE FOR ARRAYS   <---

  implicit none

  real(kind=8), allocatable, dimension (:) :: r,x,ra,V,VP1,VP2,V3,DV3
  real(kind=8), allocatable, dimension (:) :: a,alpha,krr
  real(kind=8), allocatable, dimension (:) :: a_p,alpha_p,a_pp
  real(kind=8), allocatable, dimension (:) :: phi1,psi1,pi1
  real(kind=8), allocatable, dimension (:) :: phi2,psi2,pi2
  real(kind=8), allocatable, dimension (:) :: phi3,psi3,pi3
  real(kind=8), allocatable, dimension (:) :: phi1_p,psi1_p,pi1_p
  real(kind=8), allocatable, dimension (:) :: phi2_p,psi2_p,pi2_p
  real(kind=8), allocatable, dimension (:) :: phi3_p,psi3_p,pi3_p
  real(kind=8), allocatable, dimension (:) :: sphi1,spsi1,spi1
  real(kind=8), allocatable, dimension (:) :: sphi2,spsi2,spi2
  real(kind=8), allocatable, dimension (:) :: sphi3,spsi3,spi3
  real(kind=8), allocatable, dimension (:) :: spi1_1,spi1_2,spi2_1,spi2_2,spi3_1,spi3_2
  real(kind=8), allocatable, dimension (:) :: phi
  real(kind=8), allocatable, dimension (:) :: rho,rho3
  real(kind=8), allocatable, dimension (:) :: p,p3,p_t,p3_t
  real(kind=8), allocatable, dimension (:) :: mass,omega,anisotropy,omega_t
  real(kind=8), allocatable, dimension (:) :: ADMmass,NofP,j0
  real(kind=8), allocatable, dimension (:) :: mass3,omega3
  real(kind=8), allocatable, dimension (:) :: rho_total,mass_total
  real(kind=8), allocatable, dimension (:) :: adot
  real(kind=8), allocatable, dimension (:) :: trash
  real(kind=8), allocatable, dimension (:) :: RicciScalar

  real(kind=8), allocatable, dimension (:) :: F1,F2,G1,G2
  real(kind=8), allocatable, dimension (:) :: F1_p,F2_p,G1_p,G2_p
  real(kind=8), allocatable, dimension (:) :: sF1,sF2,sG1,sG2

! All these are 1D arrays.
! The variables introduced above are.
!
! r:      Array for the radial coordinate.
! V:      Array for the potential.
! VP1:    Array for the real part of the derivative of potential.
! VP2:    Array for the real part of the derivative of potential.
!
! a:      Array for the metric coefficient (dl = a dr).
! alpha:  Array for the lapse function.
!
! phi1:   Array for phi1.
! psi1:   Array for radial derivative of phi1 (psi = dphi/dr).
! pi1:    Array for rescaled time derivative of phi1 (pi = a/alpha
!         dphi/dt).
!
! phi2:   Array for phi2.
! psi2:   Array for radial derivative of phi2 (psi = dphi/dr).
! pi2:    Array for rescaled time derivative of phi2 (pi = a/alpha
!         dphi/dt).
! phi3:   Array for phi3 [the perturbations SF].
! psi3:   Array for radial derivative of phi3 (psi = dphi/dr).
! pi3:    Array for rescaled time derivative of phi3 (pi = a/alpha
!         dphi/dt).
!
! sphi1:  Array for source term for phi1 evolution.
! spsi1:  Array for source term for phi1 evolution. 
! spi1:   Array for source term for phi1 evolution.
!
! sph2:   Array for source term for phi2 evolution.
! spsi2:  Array for source term for phi2 evolution.
! spi2:   Array for source term for phi2 evolution.
! 
! sph3:   Array for source term for phi3 evolution.
! spsi3:  Array for source term for phi3 evolution.
! spi3:   Array for source term for phi3 evolution.
!
! spix_y: Arrays for the surces using strang splitting,
!         the x indicates the label of the variable and
!         the y labels the source inside strang 1 or 2.
!
! omega:  Array for the equation of state.
! omega_t:Array for the tangential equation of state.
! p:      Array for the radial pressure
! p_t:    Array for the tangential pressure
! anisotropy: Array for the anisotropy (p_r - p_t)/p_r
! rho1:   Array for energy density of phi and phi2 [the complex field].
! rho3:   Array for energy density of phi3.
! rho:    Array for total energy density.
! phi:    Array for the modulus of phi.
! mass:   Array for mass of the Boson Star.
! ADMmass:Array for mass of the Boson Star in terms of the grr metric function.
! NofP:   Array for the number of particles.
! mass3:  Array for mass of the perturbation.
! mass_t: Array for the total mass.
! 
! adot:   Time derivative of a metric coefficient minus its   
!         source term (should converge to zero).

! --->   END   <---

  end module arrays
