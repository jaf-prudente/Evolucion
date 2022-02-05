
  subroutine initial(Nr,dr)

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

! --->   NUMBERS

  idr   = 1.0D0/dr
  smallpi    = 4.0*atan(1.0)
  zero  = 0.0D0
  third = 1.0D0/3.0D0
  half  = 0.5D0
  one   = 1.0D0
  two   = 2.0D0
  kappa0 = 2.0D0

! Initial gaussian profile.  

  do i=0,Nr
     phi1(i) = amp1*dexp(-(r(i)-r0)**2/sigma1**2) + 1.0d-20
     phi2(i) = zero
  end do

! psi using finite differencing.

  do i=1,Nr-1
     psi1(i) = half*(phi1(i+1) - phi1(i-1))*idr
     psi2(i) = half*(phi2(i+1) - phi2(i-1))*idr
  end do

! Regularity at the origin.

  psi1(0)  = - psi1(1)
  psi2(0)  = - psi2(1)

! One sided derivative at outer boundary.

  psi1(Nr) = half*(3.0d0*phi1(Nr) - 4.0d0*phi1(Nr-1) + phi1(Nr-2))
  psi2(Nr) = half*(3.0d0*phi2(Nr) - 4.0d0*phi2(Nr-1) + phi2(Nr-2))

  do i=1,Nr-1
     pi1(i) = zero
     pi2(i) = amp2*dexp(-(r(i)-r0)**2/sigma2**2) + 1.0D-20
  end do

! Find potential and its derivative.

  call potential(Nr,dr)


! Find the metric coefficients.

 
       a(0) = one                           
       a(1) = one

       do i=2,Nr

          aux = a(i-1) + (a(i-1)*(half*(one - a(i-1)**2)/r(i-1) &
             + 0.25D0*kappa0*r(i-1)*( &
             pi1(i-1)**2  + pi2(i-1)**2  + pi3(i-1)**2 &
           + psi1(i-1)**2 + psi2(i-1)**2 + psi3(i-1)**2 &
           + a(i-1)**2*(V(i-1) + V3(i-1)))))*half*dr


          a(i) = a(i-1) + &
                 dr*aux*((one - aux**2)/(r(i-1) + r(i)) + &
                  kappa0*0.25*(0.5*(r(i-1)*(pi1(i-1)**2  + pi2(i-1)**2  + pi3(i-1)**2 + &
                                            psi1(i-1)**2 + psi2(i-1)**2 + psi3(i-1)**2) + &
                                    r(i) * (pi1(i)**2    + pi2(i)**2    + pi3(i)**2 + &
                                            psi1(i)**2   + psi2(i)**2   + psi3(i)**2)) + &
                               0.5*aux**2*(r(i-1)*(V(i-1) + V3(i-1)) + &
                                           r(i)  *(V(i)   + V3(i)))))

       end do

!       alpha(Nr) = alpha(Nr)
       alpha(Nr) = one/a(Nr)

       do i=Nr-1,1,-1

           aux = alpha(i+1) - &
                 half*dr*alpha(i+1)*((a(i+1)**2-one)/r(i+1) + &
                                      ( &
                                       half*(one - a(i+1)**2)/r(i+1) + &
                                       0.25*kappa0*r(i+1)*( & 
                                                    pi1(i+1)**2 + pi2(i+1)**2 + pi3(i+1)**2 + &
                                                    psi1(i+1)**2 + psi2(i+1)**2 + psi3(i+1)**2 + &
                                                    a(i+1)**2*(V(i+1) + V3(i+1)) &
                                                    ) &
                                      ) - &
                                     kappa0*half*r(i+1)*a(i+1)**2*(V(i+1) + V3(i+1)) &
                                    )

           alpha(i) = alpha(i+1) - &
                      dr*aux*( &
                              half*((a(i+1)**2 - one)/r(i+1) + (a(i)**2 - one)/r(i)) + &
                              half*(( &
                                       half*(one - a(i+1)**2)/r(i+1) + &
                                       0.25*kappa0*r(i+1)*( &
                                                    pi1(i+1)**2 + pi2(i+1)**2 + pi3(i+1)**2 + &
                                                    psi1(i+1)**2 + psi2(i+1)**2 + psi3(i+1)**2 + &
                                                    a(i+1)**2*(V(i+1) + V3(i+1)) &
                                                    ) &
                                    ) + &
                                    ( &
                                       half*(one - a(i)**2)/r(i) + &
                                       0.25*kappa0*r(i)*( &
                                                    pi1(i)**2 + pi2(i)**2 + pi3(i)**2 + &
                                                    psi1(i)**2 + psi2(i)**2 + psi3(i)**2 + &
                                                    a(i)**2*(V(i) + V3(i)) &
                                                    ) &
                                    )) - &
                              0.25*kappa0*(r(i+1)*a(i+1)**2*V(i+1) + r(i)*a(i)**2*V(i)) &
                             )
                                
       end do

       alpha(0) = alpha(1)

! -------------------------
! -------------------------

       print *,'Recalculating the metric functions'
       print *
       print *, alpha(Nr-1)
       print *, alpha(Nr)
       print *



  end subroutine initial

