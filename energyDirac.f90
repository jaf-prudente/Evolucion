
  subroutine energyDirac(Nr,dr,dt)

! --->   ENERGY DENSITY, MASS, NUMBER OF PARTICLES et al.  <---

  use arrays

  implicit none

  common / numbers / zero,half,one,two,kappa0

  integer i,Nr

  real(8) dr,dt, idr
  real(8) zero,half,one,two,pi,kappa0
  real(8) radius95,max95,min95
  real(8) rAH

  real(8), allocatable, dimension (:) :: dF1,dF2,dG1,dG2
  allocate(dF1(0:Nr),dF2(0:Nr),dG1(0:Nr),dG2(0:Nr))


! --->   ENERGY DENSITY

! The energy density is found to be:
!
!             2       2       2       2       2
! rho = ( pi1  +  psi1  + pi2  +  psi2 ) / 2 a  +  V
!
! omega = p/rho
!
!                 2       2       2       2       2
! where p = ( pi1  +  psi1  + pi2  +  psi2 ) / 2 a  -  V
!
!

!  pi = 4.0*atan(1.0)


! dF and dG using finite differencing.

    idr   = 1.0D0/dr

    do i=1,Nr-1
      dF1(i) = half*(F1(i+1) - F1(i-1))*idr
      dF2(i) = half*(F2(i+1) - F2(i-1))*idr
      dG1(i) = half*(G1(i+1) - G1(i-1))*idr
      dG2(i) = half*(G2(i+1) - G2(i-1))*idr
    end do

!  One sided derivative at outer boundary.

    dF1(Nr) = half*(3.0d0*F1(Nr) - 4.0d0*F1(Nr-1) + F1(Nr-2))*idr
    dF2(Nr) = half*(3.0d0*F2(Nr) - 4.0d0*F2(Nr-1) + F2(Nr-2))*idr
    dG1(Nr) = half*(3.0d0*G1(Nr) - 4.0d0*G1(Nr-1) + G1(Nr-2))*idr
    dG2(Nr) = half*(3.0d0*G2(Nr) - 4.0d0*G2(Nr-1) + G2(Nr-2))*idr

    dF1(0) = half*(-3.0d0*F1(0) + 4.0d0*F1(1) - F1(2))*idr
    dF2(0) = half*(-3.0d0*F2(0) + 4.0d0*F2(1) - F2(2))*idr
    dG1(0) = half*(-3.0d0*G1(0) + 4.0d0*G1(1) - G1(2))*idr
    dG2(0) = half*(-3.0d0*G2(0) + 4.0d0*G2(1) - G2(2))*idr

    do i = 1, Nr
        rho(i) =  (1.0/(r(i)**2*a(i)**2))*( &
            a(i)*(F1(i)**2 + F2(i)**2 - G1(i)**2 - G2(i)**2) &
            + (2.0*a(i)/r(i))*(F1(i)*G2(i) - F2(i)*G1(i)) &
            + F1(i)*dG2(i) - F2(i)*dG1(i) + G1(i)*dF2(i) - G2(i)*dF1(i))
 
    end do

    rho(0)= rho(1)

! -> The current

     j0 = F1**2 + F2**2 + G1**2 + G2**2

!     rho = half*( pi1**2 + pi2**2 + &
!	          psi1**2 + psi2**2 )/a**2 + 0.5*V


!     p   = half*( pi1**2 + pi2**2 + &
!              psi1**2 + psi2**2 )/a**2 - 0.5*V
!!              psi1**2 + psi2**2 )/a**2 - V

!     p_t = half*( pi1**2 + pi2**2 - &
!              psi1**2 - psi2**2 )/a**2 - 0.5*V
!!              psi1**2 - psi2**2 )/a**2 - V

!     omega = p/rho
!     omega_t = p_t/rho
!     anisotropy = (p - p_t)/p

! Quantities for the extra field of the perturbation

!     rho3 = half*(pi3**2 + psi3**2)/a**2 + V3
!     p3   = half*(pi3**2 + psi3**2)/a**2 - V3
!     p3_t = half*(pi3**2 - psi3**2)/a**2 - V3

!     omega3 = p3/rho3
!     omega3 = zero

     rho_total = rho 


! --->   MASS, NofP et al.

! The mass function is just the integral of the energy density
! times r**2.  For this integral I use the trapezoidal rule.

! The mass at negative r makes no sense.

  mass(0)   = zero
!  mass3(0)  = zero
  mass_total(0)  = zero
  ADMmass(0) = zero
  NodP(0)    = zero

! The first point is just dr/2 away from r=0.  And since the
! integrand vanishes at the origin, we get another factor of 1/2
! from the trapezoidal rule.

  mass(1)  = 0.25D0*dr*r(1)**2*rho(1)
!  mass3(1) = 0.25D0*dr*r(1)**2*rho3(1)
  mass_total(1) = 0.25D0*dr*r(1)**2*rho_total(1)
  ADMmass(1) = 0.5*r(1)*(1.0 - 1.0/a(1)**2)
!!  NofP(1)  = 0.25D0*dr*r(1)**2*(phi1(1)*pi2(1) - phi2(1)*pi1(1))
  NodP(1)  = 0.25D0*dr*r(1)**2*j0(1)

! Integrate using trapezoidal rule for the rest of the points.
! Except the ADM mass, which is calculated using the formula
!
!      r   /     1   \
! M = --- | 1 - ----  |
!      2   \     grr /
!


  do i=2,Nr
     mass(i) = mass(i-1) &
             + half*dr*(r(i-1)**2*rho(i-1) + r(i)**2*rho(i)) 

 !    mass3(i) = mass3(i-1) &
 !             + half*dr*(r(i-1)**2*rho3(i-1) + r(i)**2*rho3(i))

     mass_total(i) = mass_total(i-1) &
             + half*dr*(r(i-1)**2*(rho(i-1)) &
             + r(i)**2*(rho(i)))

     ADMmass(i) = 0.5*r(i)*(1.0 - 1.0/a(Nr)**2)

!     NofP(i) = NofP(i-1) &
!             + half*dr*(r(i-1)**2*(phi1(i-1)*pi2(i-1) - phi2(i-1)*pi1(i-1)) + &
!                        r(i)**2*(phi1(i)*pi2(i) - phi2(i)*pi1(i)))

     NodP(i) = NodP(i-1) &
             + half*dr*(r(i-1)**2*j0(i-1) + r(i)**2*j0(i) )
  end do

!       do i=2,Nr
!          if ((two*ADMmass(Nr)/r(i).ge.0.99).and.(two*ADMmass(Nr)/r(i).le.1.0D0)) then
!              rAH = r(i)
!              print *, 'rAH =',rAH, 'M_AH=',rAH/2.
!          end if
!       end do

! The RicciScalar

  do i=1,Nr-1
    RicciScalar(i) = two*half*(alpha(i+1) - alpha(i-1))*half*(a(i+1) - a(i-1))/dr**2/alpha(i)/a(i)**3 &
                   - two*(alpha(i+1) - two*alpha(i) + alpha(i-1))/dr**2/alpha(i)/a(i)**2 &
                   + two*(a(i) - two*a_p(i) + a_pp(i))/dt*2/alpha(i)**2/a(i) &
                   - two*(a(i)-a_p(i))*(alpha(i) - alpha_p(i))/dt**2/a(i)/alpha(i)**3 &
                   - 4.0D0*half*(alpha(i+1) - alpha(i-1))/dr/alpha(i)/a(i)**2/r(i) &
                   + 4.0D0*half*(a(i+1) - a(i-1))/dr/a(i)**3/r(i) &
                   + two/r(i)**2 - two/a(i)**2/r(i)**2
  end do

  RicciScalar(0)   = 3.0D0*RicciScalar(1) - 3.0D0*RicciScalar(2) + RicciScalar(3)
  RicciScalar(Nr)  = RicciScalar(Nr-1)

! --->   END   <---

  end subroutine energyDirac

