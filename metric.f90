
  subroutine metric(Nr,dr)

! METRIC COMPONENTS  
! Integrate the equations for the metric

  use arrays

  implicit none

  common / numbers / zero,half,one,two,kappa0

  integer i,Nr

  real(8) dr,idr
  real(8) zero,half,one,two,kappa0
  real(8) aux, rescale,ssa,rk1,rhoa, rhoa1, rhoa2, rrhoa
  real(8) ssalpha, aux1, ssa2, aux2


! --->   NUMBERS

  idr = 1.0d0/dr
  V3=0.0

! --->   METRIC {a}

! For this we solve the hamiltonian constraint, which
! in this case takes the form:
!
!                   2                       2      2     2      2
! d  (ln a) = (1 - a )/2r +  kappa/4 r ( pi1 + psi1 + pi2 + psi2 + 2a V )
!  r


! Boundary values at origin. Regularity implies that a(r=0)=1.
! We need a(0)=a(r=-dr/2) and a(1)=a(r=dr/2).  Since we also know that
! da/dr=0 at the origin, with a Taylor expansion we get
! a(0) = a(1) = 1 + O(dr^3)


  a(0) = one
  a(1) = one

  do i=2,Nr

!    Advance half a grid spacing.
     rhoa = half*( (pi1(i-1)**2  + pi2(i-1)**2 &
                   + psi1(i-1)**2 + psi2(i-1)**2)/a(i-1)**2 &
                   + V(i-1) )

     ssa = a(i-1)*( half*(one - a(i-1)**2)/r(i-1) &
            + half*kappa0*r(i-1)*a(i-1)**2*rhoa )

     rk1 = a(i-1) + half*dr*ssa

!    Now evaluate the source at the mid point and
!    leapfrog over it.

     rhoa1 = half*( (pi1(i-1)**2  + pi2(i-1)**2 &
                   + psi1(i-1)**2 + psi2(i-1)**2)/rk1**2 &
                   + V(i-1) )


     rhoa2 = half*( (pi1(i)**2  + pi2(i)**2 &
                   + psi1(i)**2 + psi2(i)**2)/rk1**2 &
                   + V(i) )
 

     rrhoa = 0.5* (r(i-1)*rhoa1 + r(i)*rhoa2)


     ssa =rk1* ( (one - rk1**2)/(r(i-1) + r(i)) + kappa0*half*rk1**2*( rrhoa ))

     a(i) = a(i-1) + dr*(ssa)

  end do


! --->   LAPSE {alpha}

! This must be done after finding the metric
! coefficient "a", as it depends on it.

! The equation for the lapse (the "slicing condition")
! comes from the condition that the angular metric must
! always be r**2 ("polar-areal" slicing):
!
!                      2                             2
! d alpha = alpha ( ( a - 1 ) / r  +  d ln(a)  -  r a  V )
!  r                                   r
!
! We solve this equation integrating from the outside
! in order to impose the correct asymptotic
! boundary conditions.
!
! The outer boundary condition we use is simply:
!
! alpha = 1/a
!
! This is true in Schwarzschild, and also has to be
! true in vacuum. 
   alpha(Nr) = -a(Nr)*(phi1(Nr)/r(Nr) + psi1(Nr))/pi1(Nr)

! Integrate from the outside using second order Runge-Kutta.

!
!                      2                             2
! d alpha = alpha ( ( a - 1 ) / r  +  d ln(a)  -  r a  V )
!  r                                   r


  do i=Nr-1,1,-1

!    Move half a grid spacing.
     rhoa = half*( (pi1(i+1)**2  + pi2(i+1)**2 &
                   + psi1(i+1)**2 + psi2(i+1)**2)/a(i+1)**2 &
                   + V(i+1) )

     ssa =  half*(one - a(i+1)**2)/r(i+1) &
            + half*kappa0*r(i+1)*a(i+1)**2*rhoa 

     aux1 = ( a(i+1)**2 - 1.0d0 )/r(i+1) +  ssa  & 
            - kappa0*half*r(i+1)*a(i+1)**2*V(i+1) 

     ssalpha = alpha(i+1)*aux1 

     aux = alpha(i+1) - half*dr*ssalpha

!    Now evaluate the source at the mid point and
!    leapfrog over it.

     rhoa = half*( (pi1(i)**2  + pi2(i)**2 &
                   + psi1(i)**2 + psi2(i)**2)/a(i)**2 &
                   + V(i) )

     ssa2 =  half*(one - a(i)**2)/r(i) &
            + half*kappa0*r(i)*a(i)**2*rhoa 

     aux2 = ( a(i)**2 - 1.0d0 )/r(i) +  ssa2  & 
            - kappa0*half*r(i)*a(i)**2*V(i) 

     ssalpha = half* (aux1+aux2)
  
     alpha(i) = alpha(i+1) - dr*aux*(ssalpha)


  end do

! And at the origin.

  alpha(0) = alpha(1)


! --->   END   <---

  end subroutine metric
