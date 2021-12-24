
  subroutine reader(Nr,dr)

! --->   READING INITIAL DATA   <---

! This routine reads initial data for equilibrium
! configurations calculated elsewhere.

  use arrays

  implicit none

  common / numbers / zero,half,one,two,kappa0
  common / gauss   / r0,amp1,amp2,sigma1,sigma2
  common / gauss3  / r3,amp3,sigma3,kr

  integer i,Nr

  real(kind=8) dr,idr,dt
  real(kind=8) r0,amp1,sigma1
  real(kind=8) amp2,sigma2
  real(kind=8) r3,amp3,sigma3, aleatorio, kr
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

!-----------------------------
! --> THREE POSSIBLE USES <--
!-----------------------------

! ----------------
! README OF THIS:
! ----------------
!                 perturbed_complex will allow only to perturb the
!                 initial Re and Im of the scalar field. Moreover
!                 as I have detected addition of noise when applying
!                 perturbations of zero amplitude, I decided to add 
!                 a new variable that will permit to add or not the 
!                 perturbation above. This will make ths initial 
!                 configuration hopefully more stable.
! ----------------
! Just notice that one has to recompile if depending on the use
! one is making of the code.
! ----------------

! The following flags serve to choose a perturbation of the
! complex field (true) or a perturbation due to an extra
! field (false). If true, then it is possible to select the
! option unperturbed, for which no matters if values for
! the 3gaussian appear in the input file.

! Changing the commented lines will allow perturbed or
! unperturbed configurations....

!  perturb_complex = .false.
  perturb_complex = .true.

!  unperturbed     = .true.
  unperturbed     = .false.

! Variable to write in files

  filestatus = 'replace'

! Extra variable to allocate some numbers

  allocate(auxt(0:Nr))

! This array will serve to set zeros where required

! ---------------------------------------------------
! ----->  Reading data FROM the temporal file  <-----
! ---------------------------------------------------

   if (perturb_complex) then

! -----> If SS (1990) <-------

!    PERTURBING THE REAL AND IMAGINARY PARTS OF THE COMPLEX FIELD
! ------>

! --->   INITIAL data

     open(unit = 78,file='IDfilet.dat',status='old')
     print *, 'Reading data from ID file'
     print *

     do i=0,Nr
          read(78,*) x(i), phi1(i), phi2(i), pi1(i), pi2(i), a(i), alpha(i), rho(i)
     end do

     close(unit = 78)

! ->  Applying the perturbation

       print *, 'Applying the perturbations to the Re(phi) and Im(phi)'
       print *

       do i=0,Nr
            r(i) = x(i)
!            r(i) = x(i) - half*dr
       end do

! --> Here enters the second boolean flag

       if (unperturbed) then

          do i=0,Nr
               phi1(i) = phi1(i)
               phi2(i) = phi2(i)
               phi3(i) = 0.0D0
          end do

       else

! -----> applying the RANDOM NOISE perturbation of amplitude amp3 <-----


          do i=0,Nr
               call random_number(aleatorio)
!               phi1(i) = phi1(i) + amp3*(aleatorio-0.5D0) + 1.0D-20
!++ the kr stuff added               phi1(i) = phi1(i) + amp3*dexp(-(r(i)-r3)**2/sigma3**2) + 1.0D-20
               phi1(i) = phi1(i) + amp3*dexp(-(r(i)-r3)**2/sigma3**2)*cos(kr*r(i)) + 1.0D-20
               phi2(i) = phi2(i) - amp3*dexp(-(r(i)-r3)**2/sigma3**2)*sin(kr*r(i))
               phi3(i) = 0.0
          end do

       end if

! --> Here it finishes

! --> Note here that somethimes the perturbation has plus and other minus signs,
! --> this should simulate a perturbation by adding or removing perticles
! --> from the boson star configuration respectively. Notice this also in pi2.

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

! --> Here enter the second boolean flag

       if (unperturbed) then

          do i=1,Nr-1
               pi1(i) = zero
               pi2(i) = pi2(i)
               pi3(i) = zero
          end do

       else     

          do i=1,Nr-1
               pi1(i) = pi1(i) + a(i)*(-amp3*dexp(-(r(i)-r3)**2/sigma3**2)*sin(kr*r(i)) + 1.0D-20)/alpha(i)
!               pi1(i) = pi1(i) 
!++               pi1(i) = zero
               pi2(i) = pi2(i) + a(i)*(amp3*dexp(-(r(i)-r3)**2/sigma3**2)*cos(kr*r(i)) + 1.0D-20)/alpha(i)
!               pi2(i) = pi2(i) 
               pi3(i) = zero
          end do

       end if

       call potential(Nr,dr)

! ---------------------------------------------------------------
!      CALCULATING THE METRIC FUNCTIONS AT THE INITIAL SLICE
! ---------------------------------------------------------------

! -------------------------
! --------> AQUI <---------
! -------------------------
!       call metric(Nr,dr)  

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

       alpha(Nr) = alpha(Nr)

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

! ---------------------------------------------------
!   Initial gaussian profile for the perturbation.
! ---------------------------------------------------

   else
! -----> If Hawley and Choptuik paper Phys Rev D something (2000) <-------

!      PERTURBING WITH AN EXTERNAL SCALAR FIELD

       open(unit = 78,file='IDfilet.dat',status='old')
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

       do i=0,Nr
          phi3(i) = amp3*dexp(-(r(i)-r3)**2/sigma3**2) + 1.0D-20
       end do

   end if

! -----> Ends this if <-------


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

  end subroutine reader
