
  subroutine sources(Nr,dr)

! --->  SOURCES FOR EVOLUTION EQUATIONS   <---

! This subrotuine calculates the RHSs of the
! evolution equations for (phi,psi,pi).

  use arrays

  implicit none

  common / numbers / zero,half,one,two,kappa0

  integer i,Nr

  real(kind=8) dr,idrh
  real(kind=8) zero,half,one,two,kappa0


! --->   NUMBERS

  idrh = half/dr


!  -------------------------
! |   SOURCES FOR PHI 1     |
! |           &             |
! |   SOURCES FOR PHI 2     |
!  -------------------------

! The evolution equation for phi is:
!
! d phi1  =  alpha pi1 / a
!  t
!
! d phi2  =  alpha pi2 / a
!  t
!
! Notice that the source term in these equations
! involves no spatial derivatives, so it can be
! defined all the way to the boundaries.

  do i=0,Nr
     sphi1(i) = alpha(i)*pi1(i)/a(i)
     sphi2(i) = alpha(i)*pi2(i)/a(i)
     sphi3(i) = alpha(i)*pi3(i)/a(i)
  end do


!  --------------------------
! |   SOURCES FOR PSI's      |
!  --------------------------

! The evolution equation for psi is:
!
! d psi1  =  d  ( alpha pi1 / a )
!  t           r
!
! d psi2  =  d  ( alpha pi2 / a )
!  t           r
!
! Since the source term now involves spatial
! derivatives, it can not be calculated at
! the boundaries.

  do i=1,Nr-1
     spsi1(i) = idrh*(alpha(i+1)*pi1(i+1)/a(i+1) &
             - alpha(i-1)*pi1(i-1)/a(i-1))
     spsi2(i) = idrh*(alpha(i+1)*pi2(i+1)/a(i+1) &
             - alpha(i-1)*pi2(i-1)/a(i-1))
     spsi3(i) = idrh*(alpha(i+1)*pi3(i+1)/a(i+1) &
             - alpha(i-1)*pi3(i-1)/a(i-1))
  end do


!  ---------------------------
! |     SOURCES FOR PI's      |
!  ---------------------------

! The evolution equation for pi1 is:
!
!           -2       2
! d pi1  =  r   d  ( r  alpha phi1 / a )  -  a alpha VP1
!  t             r
!
! This can be rewritten as:
!
!                   2                      3
! d pi1  =  3  d ( r  alpha phi1 / a ) / dr  -  a alpha VP1
!  t            
!
! The evolution equation for pi1 is:
!
!           -2       2
! d pi2  =  r   d  ( r  alpha phi2 / a )  -  a alpha VP2
!  t             r
!
! This can be rewritten as:
!
!                   2                      3
! d pii  =  3  d ( r  alpha phii / a ) / dr  -  a alpha VPi
!  t            
!
! All we have done here is rewrite the first term as a
! derivative with respect to r^3 (this is not a third
! derivative!).  The reason for this is to ensure that
! our finite difference approximation is consistent at
! the origin.
!
! Again, since the source term now involves spatial
! derivatives, it can not be calculated up to the boundaries.

  do i=1,10
     spi1(i) = 3.0D0*((alpha(i+1)*r(i+1)**2*psi1(i+1)/a(i+1) &
            - alpha(i-1)*r(i-1)**2*psi1(i-1)/a(i-1)) &
            /(r(i+1)**3 - r(i-1)**3)) &
            - a(i)*alpha(i)*VP1(i)

     spi2(i) = 3.0D0*((alpha(i+1)*r(i+1)**2*psi2(i+1)/a(i+1) &
            - alpha(i-1)*r(i-1)**2*psi2(i-1)/a(i-1)) &
            /(r(i+1)**3 - r(i-1)**3)) &
            - a(i)*alpha(i)*VP2(i)
    
  end do




  do i=10,Nr-1
     spi1(i) = idrh*((alpha(i+1)*r(i+1)**2*psi1(i+1)/a(i+1) &
            - alpha(i-1)*r(i-1)**2*psi1(i-1)/a(i-1)) &
            )/r(i)**2 &
            - a(i)*alpha(i)*VP1(i)

     spi2(i) = idrh*((alpha(i+1)*r(i+1)**2*psi2(i+1)/a(i+1) &
            - alpha(i-1)*r(i-1)**2*psi2(i-1)/a(i-1)) &
            )/r(i)**2 &
            - a(i)*alpha(i)*VP2(i)
  end do







! If strang is used, then the following sources apply.

  do i=1,Nr-1
     spi1_1(i) = 3.0D0*((alpha(i+1)*r(i+1)**2*psi1(i+1)/a(i+1) &
               - alpha(i-1)*r(i-1)**2*psi1(i-1)/a(i-1)) &
                 /(r(i+1)**3 - r(i-1)**3))
     spi1_2(i) = - a(i)*alpha(i)*VP1(i)
     spi2_1(i) = 3.0D0*((alpha(i+1)*r(i+1)**2*psi2(i+1)/a(i+1) &
               - alpha(i-1)*r(i-1)**2*psi2(i-1)/a(i-1)) &
                 /(r(i+1)**3 - r(i-1)**3))
     spi2_2(i) = - a(i)*alpha(i)*VP2(i)
     spi3_1(i) = 3.0D0*((alpha(i+1)*r(i+1)**2*psi3(i+1)/a(i+1) &
               - alpha(i-1)*r(i-1)**2*psi3(i-1)/a(i-1)) &
                 /(r(i+1)**3 - r(i-1)**3))
     spi3_2(i) = - a(i)*alpha(i)*DV3(i)
  end do

! And the sources labeled with 2 can be known up to 0 and Nr

  spi1_2(0 ) = - a(0 )*alpha(0 )*VP1(0 )
  spi1_2(Nr) = - a(Nr)*alpha(Nr)*VP1(Nr)
  spi2_2(0 ) = - a(0 )*alpha(0 )*VP2(0 )
  spi2_2(Nr) = - a(Nr)*alpha(Nr)*VP2(Nr)
  spi3_2(0 ) = - a(0 )*alpha(0 )*DV3(0 )
  spi3_2(Nr) = - a(Nr)*alpha(Nr)*DV3(Nr)

! --->   END   <---

  end subroutine sources

