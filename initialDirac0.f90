
  subroutine initialDirac0(Nr,dr)

! GAUSSIAN INITIAL DATA 

! This routine sets an initial gaussian for Re(psi)

  use arrays

  implicit none

  common / numbers / zero,half,one,two,kappa0
  common / gauss   / amp1,amp2,r0,sigma1,sigma2
  common / gauss3  / r3,amp3,sigma3

  integer i,Nr

  real(8) dr,idr
  real(8) amp1,r0,sigma1
  real(8) amp2,sigma2
  real(8) r3,amp3,sigma3
  real(8) zero,half,one,two,kappa0,third
  real(8) aux
  real(8) smallpi
  real(8) ro, Srr, ro_aux, Srr_aux

  real(8), allocatable, dimension (:) :: dF1,dF2,dG1,dG2
  allocate(dF1(0:Nr),dF2(0:Nr),dG1(0:Nr),dG2(0:Nr))

! --->   NUMBERS

  idr = 1.0D0/dr
  smallpi = 4.0*atan(1.0)
  zero = 0.0D0
  third = 1.0D0/3.0D0
  half = 0.5D0
  one = 1.0D0
  two = 2.0D0
  kappa0 = 2.0D0

! Initial gaussian profile.  

  do i=0,Nr
     F1(i) = amp1*dexp(-(r(i)-r0)**2/sigma1**2) + 1.0d-20
     F2(i) = zero
     G1(i) = zero
     G2(i) = zero
  end do

! F and G using finite differencing.

  do i=1,Nr-1
    dF1(i) = half*(F1(i+1) - F1(i-1))*idr
    dF2(i) = half*(F2(i+1) - F2(i-1))*idr
    dG1(i) = half*(G1(i+1) - G1(i-1))*idr
    dG2(i) = half*(G2(i+1) - G2(i-1))*idr
  end do

  ! One sided derivative at outer boundary.

  dF1(Nr) = half*(3.0d0*F1(Nr) - 4.0d0*F1(Nr-1) + F1(Nr-2))*idr
  dF2(Nr) = half*(3.0d0*F2(Nr) - 4.0d0*F2(Nr-1) + F2(Nr-2))*idr
  dG1(Nr) = half*(3.0d0*G1(Nr) - 4.0d0*G1(Nr-1) + G1(Nr-2))*idr
  dG2(Nr) = half*(3.0d0*G2(Nr) - 4.0d0*G2(Nr-1) + G2(Nr-2))*idr

  dF1(Nr) = half*(-3.0d0*F1(0) + 4.0d0*F1(1) - F1(2))*idr
  dF2(Nr) = half*(-3.0d0*F2(0) + 4.0d0*F2(1) - F2(2))*idr
  dG1(Nr) = half*(-3.0d0*G1(0) + 4.0d0*G1(1) - G1(2))*idr
  dG2(Nr) = half*(-3.0d0*G2(0) + 4.0d0*G2(1) - G2(2))*idr

! Regularity at the origin.

  F1(0) = - F1(1)
  F2(0) = -F2(1)
  G1(0) = G1(1)
  G2(0) = G2(1)

! Find the metric coefficients.

       a(0) = one                           
       a(1) = one

       do i=2,Nr

          ro = (1/(r(i-1)**2*a(i-1)**2))*( &
            a(i-1)*(F1(i-1)**2 + F2(i-1)**2 - G1(i-1)**2 - G2(i-1)**2) &
            + (2*a(i-1)/r(i-1))*(F1(i-1)*G2(i-1) - F2(i-1)*G1(i-1)) &
            + F1(i-1)*dG2(i-1) - F2(i-1)*dG1(i-1) + G1(i-1)*dF2(i-1) - G2(i-1)*dF1(i-1))

          aux = a(i-1) + &
            (a(i-1)*(half*(one - a(i-1)**2)/r(i-1) &
            + r(i-1)*a(i-1)**2*ro))*dr*half

          ro_aux = (1/(r(i-1)**2*aux**2))*( &
          aux*(F1(i-1)**2 + F2(i-1)**2 - G1(i-1)**2 - G2(i-1)**2) &
          + (2*aux/r(i-1))*(F1(i-1)*G2(i-1) - F2(i-1)*G1(i-1)) &
          + F1(i-1)*dG2(i-1) - F2(i-1)*dG1(i-1) + G1(i-1)*dF2(i-1) - G2(i-1)*dF1(i-1))

          a(i) = a(i-1) + &
            (aux*(half*(one - aux**2)/r(i-1) &
            + r(i-1)*aux**2*ro_aux))*dr

       end do

!      alpha(Nr) = alpha(Nr)
       alpha(Nr) = one/a(Nr)

       do i=Nr-1,1,-1
          
          Srr = (1/(r(i+1)**2*a(i+1)**2))*(F1(i+1)*dG2(i+1) &
            - F2(i+1)*dG1(i+1) + G1(i+1)*dF2(i+1) - G2(i+1)*dF1(i+1))

          aux = alpha(i+1) - &
            (alpha(i+1)*(half*(one - a(i+1)**2)/r(i+1) &
            - r(i+1)*a(i+1)**2*Srr))*dr*half
          
          Srr_aux = (1/(r(i+1)**2*aux**2))*(F1(i+1)*dG2(i+1) &
          - F2(i+1)*dG1(i+1) + G1(i+1)*dF2(i+1) - G2(i+1)*dF1(i+1))

          alpha(i) = alpha(i+1) - &
            (aux*(half*(one - a(i+1)**2)/r(i+1) &
            - r(i+1)*a(i+1)**2*Srr_aux))*dr
                                
       end do

       alpha(0) = alpha(1)

! -------------------------
! -------------------------

       print *,'Recalculating the metric functions'
       print *
       print *, alpha(Nr-1)
       print *, alpha(Nr)
       print *



  end subroutine initialDirac0

