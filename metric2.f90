
  subroutine metric(Nr,dr)

! METRIC COMPONENTS  
! Integrate the equations for the metric

  use arrays

  implicit none

  common / numbers / zero,half,one,two,kappa0

  integer i,Nr,l

  real(8) dr,idr
  real(8) zero,half,one,two,kappa0
  real(8) aux, rescale,ssa,rk1,rhoa, rhoa1, rhoa2, rrhoa
  real(8) ssalpha, aux1, ssa2, aux2
  real(8) A_rk4,phi1_rk4,pi1_rk4,psi1_rk4, alpha_rk4,rm
  real(8) V_rk4,phi2_rk4,pi2_rk4,psi2_rk4
  real(8) k11,k12,k13,k14,J1c


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

  do l=2,Nr-1

!    First step.

     rm        =  r(l)
     a_rk4     =  a(l)

     alpha_rk4 =  alpha(l)
     phi1_rk4   =  phi1(l)
     pi1_rk4    =  pi1(l)
     psi1_rk4 =    psi1(l)
     phi2_rk4   =  phi2(l)
     pi2_rk4    =  pi2(l)
     psi2_rk4 =    psi2(l)
     V_rk4 =    V(l)

    k11   =  dr*J1c(a_rk4,alpha_rk4,phi1_rk4,pi1_rk4,&
                    psi1_rk4,phi2_rk4,pi2_rk4,psi2_rk4,V_rk4,rm)

!    Second step.  We use forth order interpolation for all other variables.

     rm    = r(l) + 0.5*dr
     a_rk4 = a(l) + 0.5*k11


     if (l==Nr-1) then
        alpha_rk4  = 0.5*(alpha(l)      +  alpha(l+1))
        phi1_rk4   = 0.5*(phi1(l) + phi1(l+1))
        psi1_rk4   = 0.5*(psi1(l) + psi1(l+1))
        pi1_rk4    = 0.5*(pi1(l)  + pi1(l+1))
        phi2_rk4   = 0.5*(phi2(l) + phi2(l+1))
        psi2_rk4   = 0.5*(psi2(l) + psi2(l+1))
        pi2_rk4    = 0.5*(pi2(l)  + pi2(l+1))
        V_rk4    = 0.5*(V(l)  + V(l+1))
     else

        alpha_rk4   = (9.0*(alpha(l ) + alpha(l+1)) &
                       -  (alpha(l-1) + alpha(l+2)))/16.0

        phi1_rk4   = (9.0*(phi1(l  ) + phi1(l+1)) &
                       - (phi1(l-1) + phi1(l+2)))/16.0

        psi1_rk4    = (9.0*(psi1(l  ) + psi1(l+1)) &
                       - (psi1(l-1) + psi1(l+2)))/16.0

        pi1_rk4    = (9.0*(pi1(l  ) + pi1(l+1)) &
                       - (pi1(l-1) + pi1(l+2)))/16.0

        phi2_rk4   = (9.0*(phi2(l  ) + phi2(l+1)) &
                       - (phi2(l-1) + phi2(l+2)))/16.0

        psi2_rk4    = (9.0*(psi2(l  ) + psi2(l+1)) &
                       - (psi2(l-1) + psi2(l+2)))/16.0

        pi2_rk4    = (9.0*(pi2(l  ) + pi2(l+1)) &
                       - (pi2(l-1) + pi2(l+2)))/16.0

        V_rk4    = (9.0*(V(l  ) + V(l+1)) &
                       - (V(l-1) + V(l+2)))/16.0

     end if

     k12 = dr*J1c(a_rk4,alpha_rk4,phi1_rk4,pi1_rk4,psi1_rk4,phi2_rk4,pi2_rk4,psi2_rk4,V_rk4,rm)

!    Third step.

     a_rk4 = a(l) + 0.5*k12

     k13 = dr*J1c(a_rk4,alpha_rk4,phi1_rk4,pi1_rk4,psi1_rk4,phi2_rk4,pi2_rk4,psi2_rk4,V_rk4,rm)

!    Fourth step.

     rm    = r(l) + dr
     A_rk4 = A(l) + k13

     alpha_rk4 = alpha(l+1)
     phi1_rk4   = phi1(l+1)
     pi1_rk4    = pi1(l+1)
     psi1_rk4    = psi1(l+1)

     V_rk4 = V(l+1)
     phi2_rk4   = phi2(l+1)
     pi2_rk4    = pi2(l+1)
     psi2_rk4    = psi2(l+1)


     k14 = dr*J1c(a_rk4,alpha_rk4,phi1_rk4,pi1_rk4,psi1_rk4,phi2_rk4,pi2_rk4,psi2_rk4,V_rk4,rm)

!    Advance variable.

     a(l+1) = a(l) + (k11 + 2.0*(k12 + k13) + k14)/6.0

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

! **********************
! *** End subroutine ***
! **********************

  end subroutine metric

! *****************
! *** Function  ***
! *****************

! Sources for the metric function 
   
  function J1c(a_rk4,alpha_rk4,phi1_rk4,pi1_rk4,psi1_rk4,phi2_rk4,pi2_rk4,psi2_rk4,V_rk4,rm)
    implicit none
    common / numbers / zero,half,one,two,kappa0    

    real(8) J1c,A_rk4,alpha_rk4,rm
    real(8) pi1_rk4,phi1_rk4,psi1_rk4
    real(8) pi2_rk4,phi2_rk4,psi2_rk4, V_rk4
    real(8) zero,half,one,two,kappa0
    real(8) rhoa


     rhoa = 0.5*( (pi1_rk4**2  + pi2_rk4**2 &
                   + psi1_rk4**2 + psi2_rk4**2)/a_rk4**2 &
                   + V_rk4 )

       J1C = a_rk4*( half*(one - a_rk4**2)/rm &
            + half*kappa0*rm*a_rk4**2*rhoa )

 end function J1c




