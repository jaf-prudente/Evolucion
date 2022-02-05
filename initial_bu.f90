
  subroutine initial(Nr,dr)

! --->   READING INITIAL DATA   <---

! This routine reads initial data for equilibrium
! configurations calculated elsewhere.

  use arrays

  implicit none

  common / gauss   / amp1,amp2,r0,sigma1,sigma2
  common / numbers / zero,half,one,two,kappa0
!  common / gauss   / r0,amp1,amp2,sigma1,sigma2
  common / gauss3  / r3,amp3,sigma3,kr

  integer i,Nr

  real(kind=8) dr,idr,dt
  real(kind=8) r0,amp1,sigma1
  real(kind=8) amp2,sigma2
  real(kind=8) r3,amp3,sigma3,kr
  real(kind=8) zero,half,one,two,third,kappa0
  real(kind=8) pi,aux
  real(kind=8), allocatable, dimension (:) :: auxt

  logical      perturb_complex, unperturbed
  character*20 filestatus

! --->   NUMBERS

  idr   = 1.0D0/dr
  pi    = 4.0*atan(1.0)
  zero  = 0.0D0
  third = 1.0D0/3.0D0
  half  = 0.5D0
  one   = 1.0D0
  two   = 2.0D0
  kappa0 = 2.0D0

! Variable to write in files

  filestatus = 'replace'

! Extra variable to allocate some numbers

  allocate(auxt(0:Nr))

! This array will serve to set zeros where required

! --->   INITIAL data

     phi1 = amp1*dexp(-r**2/sigma1**2) + 1.0D-20
     phi2 = zero
     phi3 = amp3*dexp(-r**2/sigma3**2) + 1.0D-20


! --> Apply a perturbation (zero mos of the times)

     phi1 = phi1 + amp3*dexp(-(r)**2/sigma3**2) + 1.0D-20
     phi2 = phi2
     phi3 = 0.0

     pi1 = zero
     pi2 = amp1*dexp(-r**2/sigma1**2) + 1.0D-20
     pi3 = zero

  print *, ' --------------------------------------------------'
  print *, ' -----  You are evolving an initial GAUSSIAN  -----'
  print *, ' -------     amp1   =', amp1
  print *, ' -------     sigma1 =', sigma1
  print *, ' --------------------------------------------------'

! --> Here it finishes

       do i=1,Nr-1
            psi1(i) = half*(phi1(i+1) - phi1(i-1))*idr
            psi2(i) = half*(phi2(i+1) - phi2(i-1))*idr
            psi3(i) = half*(phi3(i+1) - phi3(i-1))*idr
       end do
  
! Regularity at origin.

            psi1(0)  = - psi1(1)
            psi2(0)  = - psi2(1)
            psi3(0)  = - psi3(1)
        
! One sided derivative at outer boundary.

            psi1(Nr) = half*(3.0D0*phi1(Nr) - 4.0D0*phi1(Nr-1) + phi1(Nr-2))
            psi2(Nr) = half*(3.0D0*phi2(Nr) - 4.0D0*phi2(Nr-1) + phi2(Nr-2))
            psi3(Nr) = half*(3.0D0*phi3(Nr) - 4.0D0*phi3(Nr-1) + phi3(Nr-2))

       call potential(Nr,dr)

! ---------------------------------------------------------------
!      CALCULATING THE METRIC FUNCTIONS AT THE INITIAL SLICE
! ---------------------------------------------------------------

! >-<       call metric(Nr,dr)  

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

!      Set up alpha at the outermost becasuse it has been forgotten.

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


!   WRITING ONCE AGAIN TO THE FILE IDfile.dat
! ------>

      print *,'Writing once again to IDfile.dat the perturbed data'
      print *

      auxt = 0.0
      do i=0,Nr
!!            r(i) = r(i) + half*dr
            r(i) = x(i)
      end do
      
     do i=0,Nr
       r(i) = dble(i)*dr
     end do

      if (filestatus == 'replace') then
         open(1,file='IDfile.dat',form='formatted',status=filestatus)
      else  
         open(1,file='IDfile.dat',form='formatted',status=filestatus, &
         position='append')
      end if

      do i=0,Nr
         write(1,"(8d24.16)") r(i),phi1(i),auxt(i),auxt(i),phi1(i)*a(i)/alpha(i),a(i),alpha(i),rho(i)
      end do

      close(1)

!     REREADING THE DATA FROM IDfile.dat

      r     = zero
      phi1  = zero
      phi2  = zero
      pi1   = zero
      pi2   = zero
      a     = zero
      alpha = zero
      rho   = zero

      print *,'Reading once again the data from IDfile.dat'

      open(unit = 78,file='IDfile.dat',status='old')
      print *, 'Reading data from ID file'
      print *
  
      do i=0,Nr
           read(78,*) x(i), phi1(i), phi2(i), pi1(i), pi2(i), a(i), alpha(i), rho(i)
      end do

      close(unit = 78)

       do i=0,Nr
          r(i) = x(i)
!          r(i) = x(i) - half*dr
       end do

!#######################################################################


! Radial derivative of phi (psi) using finite differencing.

     do i=1,Nr-1
        psi1(i) = half*(phi1(i+1) - phi1(i-1))*idr
        psi2(i) = half*(phi2(i+1) - phi2(i-1))*idr
        psi3(i) = half*(phi3(i+1) - phi3(i-1))*idr
     end do

! Regularity at origin.

     psi1(0)  = - psi1(1)
     psi2(0)  = - psi2(1)
     psi3(0)  = - psi3(1)

! One sided derivative at outer boundary.

     psi1(Nr) = half*(3.0D0*phi1(Nr) - 4.0D0*phi1(Nr-1) + phi1(Nr-2))
     psi2(Nr) = half*(3.0D0*phi2(Nr) - 4.0D0*phi2(Nr-1) + phi2(Nr-2))
     psi3(Nr) = half*(3.0D0*phi3(Nr) - 4.0D0*phi3(Nr-1) + phi3(Nr-2))


! --->   Calculates pi for the main field and the perturbation

     do i=1,Nr-1
        pi1(i) = zero
        pi3(i) = zero
     end do

! --->   Calculates the potential

     call potential(Nr,dr)

! --->   END   <---

  end subroutine initial
