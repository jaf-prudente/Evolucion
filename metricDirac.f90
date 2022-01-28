
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

! RESOLVER (36) DE DAKA
! aquí van las ecuaciones (36) de Daka


  a(0) = one
  a(1) = one

  do i=2,Nr

    ! ESTE ES EL MÉTODO EN SÍ PARA RESOLVER (36)
!    rhoa es rho de (37).
     rhoa = half*( (pi1(i-1)**2  + pi2(i-1)**2 &
                   + psi1(i-1)**2 + psi2(i-1)**2)/a(i-1)**2 &
                   + V(i-1) )
!   ssa es (36) en sí
     ssa = a(i-1)*( half*(one - a(i-1)**2)/r(i-1) &
            + half*kappa0*r(i-1)*a(i-1)**2*rhoa )
!   este es el método de rk1 con un paso intermedio
     rk1 = a(i-1) + half*dr*ssa

!    Método de la rana
     ! rhoa1 es rho evaluada en i-1
     rhoa1 = half*( (pi1(i-1)**2  + pi2(i-1)**2 &
                   + psi1(i-1)**2 + psi2(i-1)**2)/rk1**2 &
                   + V(i-1) )

    ! rhoa2 es rho evaluada en i
     rhoa2 = half*( (pi1(i)**2  + pi2(i)**2 &
                   + psi1(i)**2 + psi2(i)**2)/rk1**2 &
                   + V(i) )
 
    ! rrhoa es el promedio
     rrhoa = 0.5* (r(i-1)*rhoa1 + r(i)*rhoa2)

    ! ssa es (36)con el promedio de rho y con rk1 en lugar de a
     ssa =rk1* ( (one - rk1**2)/(r(i-1) + r(i)) + kappa0*half*rk1**2*( rrhoa ))
     
    ! El método de Euler
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
