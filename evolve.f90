
  subroutine evolve(Nr,Nt,Noutput0D,Noutput1D,dr,dt,mattertype)

! *************************
! *** Evolving in time  ***
! *************************

! This subroutine evolves the scalar field and the
! metric variables using a three step Iterative 
! Cranck-Nicholson method.

  use arrays

  implicit none

  common / numbers / zero,half,one,two,kappa0
  common / initialD / IData  

  logical icnthird, rk3

  integer i,k,l
  integer Nr,Nt,Noutput0D,Noutput1D,mattertype
  integer IData
  integer num_steps

  real(8) dr,dt,t
  real(8) dtfac,dtw
  real(8) zero,third,half,one,two,kappa0
  real(8) aux,r95
  real(8), allocatable, dimension (:) :: temp,m95
  real(8) quarter

! Numbers

  zero  = 0.0D0
  third = 1.0D0/3.0D0
  half  = 0.5D0
  one   = 1.0D0
  two   = 2.0D0
  kappa0 = 2.0D0
  quarter = 0.25D0


! Allocate arrays 

  allocate(r(0:Nr),x(0:Nr),ra(0:Nr),V(0:Nr),VP1(0:Nr),VP2(0:Nr))
  allocate(V3(0:Nr),DV3(0:Nr))  

  allocate(a(0:Nr),alpha(0:Nr),krr(0:Nr))
  allocate(a_p(0:Nr),alpha_p(0:Nr),a_pp(0:Nr))

  allocate(phi1(0:Nr),psi1(0:Nr),pi1(0:Nr))
  allocate(phi1_p(0:Nr),psi1_p(0:Nr),pi1_p(0:Nr))
  allocate(sphi1(0:Nr),spsi1(0:Nr),spi1(0:Nr))

  allocate(phi2(0:Nr),psi2(0:Nr),pi2(0:Nr))
  allocate(phi2_p(0:Nr),psi2_p(0:Nr),pi2_p(0:Nr))
  allocate(sphi2(0:Nr),spsi2(0:Nr),spi2(0:Nr))

  allocate(phi3(0:Nr),psi3(0:Nr),pi3(0:Nr))
  allocate(phi3_p(0:Nr),psi3_p(0:Nr),pi3_p(0:Nr))
  allocate(sphi3(0:Nr),spsi3(0:Nr),spi3(0:Nr))

  allocate(F1(0:Nr),F2(0:Nr),G1(0:Nr),G2(0:Nr))
  allocate(F1_p(0:Nr),F2_p(0:Nr),G1_p(0:Nr),G2_p(0:Nr))
  allocate(sF1(0:Nr),sF2(0:Nr),sG1(0:Nr),sG2(0:Nr))

  allocate(phi(0:Nr))
  allocate(rho(0:Nr),rho3(0:Nr))
  allocate(p(0:Nr),p3(0:Nr),p_t(0:Nr),p3_t(0:Nr))
  allocate(mass(0:Nr),omega(0:Nr),omega_t(0:Nr),ADMmass(0:Nr),NofP(0:Nr))
  allocate(mass3(0:Nr),omega3(0:Nr),anisotropy(0:Nr),j0(0:Nr))
  allocate(rho_total(0:Nr),mass_total(0:Nr),NodP(0:Nr))

  allocate(adot(0:Nr))
  allocate(trash(0:Nr))
  allocate(m95(0:Nr),temp(0:Nr))
  allocate(RicciScalar(0:Nr))

  allocate(spi1_1(0:Nr),spi1_2(0:Nr),spi2_1(0:Nr),spi2_2(0:Nr),spi3_1(0:Nr),spi3_2(0:Nr))

! Initialise time

  t = zero

! Grid
!  Notice that the origin has been staggered,
!  the grid point i=0 is at r=-dr/2.

  do i=0,Nr
     r(i) = dble(i)*dr - half*dr
  end do

! Info to screen

  print *,'----------------------------'
  print *,'|  Time step  |    Time    |'
  print *,'----------------------------'

  write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',0,'    | ',t,'  |'

! clean all the variables

  a_p = 0.0d0
  a_pp = 0.0d0
  alpha_p = 0.0d0

  phi1_p = 0.0d0
  psi1_p = 0.0d0
  pi1_p  = 0.0d0

  phi2_p = 0.0d0
  psi2_p = 0.0d0
  pi2_p  = 0.0d0

  adot = 0.0d0

  F1_p = 0.0d0
  F2_p = 0.0d0
  G1_p  = 0.0d0
  G2_p  = 0.0d0
  

! Call Initial data
  if(mattertype==0) then
   call initial(Nr,dr)

  else if(mattertype==1) then
   call initialDirac0(Nr,dr)

  end if


! Compute Energy 

  if(mattertype==0) then

    call energy(Nr,dr,dt)

  else if(mattertype==1) then

    call energyDirac(Nr,dr,dt)
    print*, "energy dirac"
   
  end if



       temp(0) = zero
       temp(1) = 0.25D0*dr*r(1)**2*rho(1)
       m95(0) = zero
       m95(1) = 0.25D0*dr*r(1)**2*rho(1)   
     
       do i=2,Nr
          temp(i) = temp(i-1) &
                  + half*dr*(r(i-1)**2*rho(i-1) + r(i)**2*rho(i))
          ADMmass(i) = 0.5*r(i)*(1.0 - 1.0/a(Nr)**2)
       end do

       do i=2,Nr
          m95(i) = m95(i-1) &
                  + half*dr*(r(i-1)**2*rho(i-1) + r(i)**2*rho(i))
          if ((m95(i)/temp(Nr).ge.0.95).and.(m95(i)/temp(Nr).le.0.96)) then
              r95 = r(i)
              print *, 'r95 =',r95
          end if
       end do
              print *, 'The initial ADM mass is     =',ADMmass(Nr)
              print *, 'The size of the grid is =',dr*Nr


! Save initial data

  call save0Ddata(Nr,t,dr)
  call save1Ddata(Nr,t)

! Set number of iterations
! this is done thinking in other evolving
! algorithms 
  num_steps = 3

! Main iteration loop.

  do l=1,Nt

!    Time.

     t = t + dt

!    Time step informtaion to screen.

     if (mod(l,Noutput1D).eq.0) then
        write(*,"(A5,I6,A6,ES9.2,A3)") ' |   ',l,'    | ',t,'  |'
     end if

!    Save old time step.

     a_pp = a_p
     a_p = a
     alpha_p = alpha

! ---> matter

!    scalar 
     phi1_p = phi1
     psi1_p = psi1
     pi1_p  = pi1

     phi2_p = phi2
     psi2_p = psi2
     pi2_p  = pi2

!    Dirac
     F1_p = F1
     F2_p = F2
     G1_p = G1
     G2_p = G2

!    Use Iterative Crank-Nicholson (ICN) with 3 iterations.

     do k=1,num_steps

         if (k == 1) then
            dtw = half*dt
         else if (k == 2) then
            dtw = half*dt
         else
            dtw = dt
         end if

       if(mattertype==0) then

!       Calculate RHSs.

        call sources(Nr,dr)
 
        do i=0,Nr
           phi1(i) = phi1_p(i) + dtw*sphi1(i)
           phi2(i) = phi2_p(i) + dtw*sphi2(i)
        end do

        do i=1,Nr-1
           psi1(i) = psi1_p(i) + dtw*spsi1(i)
           psi2(i) = psi2_p(i) + dtw*spsi2(i)
        end do

        do i=1,Nr-1
           pi1(i)  = pi1_p(i)  + dtw*spi1(i)
           pi2(i)  = pi2_p(i)  + dtw*spi2(i)
        end do

!     Regularity boundary conditions at origin for (psi,pi):
!     1) psi must be odd.
!     2) pi must be even.

        psi1(0) = - psi1(1)
        pi1(0)  = + pi1(1)

        psi2(0) = - psi2(1)
        pi2(0)  = + pi2(1)

!       Outer boundary condition for pi.  For pi I use an
!       outgoing wave boundary condition:
!
!       f  =  u(r - t)/r
!
!       For some unknown function u.  This condition implies:
!
!       d f  +  d f  +  f / r  =  0
!        r       t
!
!       For the boundary condition I use a second order accurate
!       finite differencing of this last expression, and solve for
!       the boundary value.

        dtfac = dtw/dr
        aux = dtw/(r(Nr) + r(Nr-1))

        pi1(Nr) = (pi1(Nr-1)*(dtfac - one - aux) &
                  + pi1_p(Nr  )*(one - dtfac - aux) &
                  + pi1_p(Nr-1)*(one + dtfac - aux)) &
                   /(one + dtfac + aux)

        pi2(Nr) = (pi2(Nr-1)*(dtfac - one - aux) &
                  + pi2_p(Nr  )*(one - dtfac - aux) &
                  + pi2_p(Nr-1)*(one + dtfac - aux)) &
                   /(one + dtfac + aux)

!       Outer boundary condition for psi.  Psi does not obey a simple
!       outgoing wave boundary condition.  But both phi and pi do.
!       The condition applied to phi in particular implies that:
!
!       psi1 = - pi1 - phi1/r
!       psi2 = - pi2 - phi2/r
!
!       where we assumed that alpha=a=1 at the bounadry.

!       psi1(Nr) = - pi1(Nr) - phi1(Nr)/r(Nr)
        psi1(Nr) = - alpha(Nr)*pi1(Nr)/a(Nr) - phi1(Nr)/r(Nr)

!       psi2(Nr) = - pi2(Nr) - phi2(Nr)/r(Nr)
        psi2(Nr) = - alpha(Nr)*pi2(Nr)/a(Nr) - phi2(Nr)/r(Nr)


!       Get new potential (must be done before integrating
!       the metric components because it is needed then).

        call potential(Nr,dr)

!       Integrate metric coefficients (a,alpha).

        call metric(Nr,dr)

       else if(mattertype==1) then

        call sourcesDirac(Nr,dr)
 
        do i=1,Nr-1
           F1(i) = F1_p(i) + dtw*sF1(i)
           F2(i) = F2_p(i) + dtw*sF2(i)
           G1(i) = G1_p(i) + dtw*sG1(i)
           G2(i) = G2_p(i) + dtw*sG2(i)
        end do

!     Regularity boundary conditions at origin for (f's and g's):
!     1) f1,f2 must be odd.
!     2) g1,g2 must be even.

        F1(0) = - F1(1)
        G1(0)  = + G1(1)

        F2(0) = - F2(1)
        G2(0)  = + G2(1)

!       Outer boundary condition 

        F1(Nr) = F1(Nr-1)
        G1(Nr) = G1(Nr-1)

        F2(Nr) = F2(Nr-1)
        G2(Nr) = G2(Nr-1)

!
!       Integrate metric coefficients (a,alpha).

        call metricDirac(Nr,dr)

       end if

     end do


!  Compute Energy 

    if(mattertype==0) then

      call energy(Nr,dr,dt)

    else if(mattertype==1) then

      call energyDirac(Nr,dr,dt)
    
    end if



!    Calculate the Krr component of the extrinsic curvature

!    krr = half*d(grr)/dt/alpha

     krr = 0.5D0*(a**2-a_p**2)/dt/alpha

!    Calculate the constraint "adot"

!    "adot" is the time derivative of a,
!     minus its source term:
!
!    adot  :=  da/dt  -  r alpha ( pi1 psi1 + pi2 psi2 )  =  0
!
!    This quantity should be zero analytically, but not numerically.
!    Numerically, it should converge to zero to second order.

   if(mattertype==0) then

     do i=0,Nr
        adot(i) = (a(i) - a_p(i))/dt - 0.25D0*r(i) &
           *(alpha(i)*pi1(i)*psi1(i) &
           + alpha_p(i)*pi1_p(i)*psi1_p(i) &  
         + alpha(i)*pi2(i)*psi2(i) &
           + alpha_p(i)*pi2_p(i)*psi2_p(i)  )
     end do

    else if(mattertype==1) then

    
    
    end if





!    Save 0D data every Noutput0D time steps.

     if (mod(l,Noutput0D).eq.0) then
        call save0Ddata(Nr,t,dr)
     end if

!    Save 1D data every Noutput1D time steps.

     if (mod(l,Noutput1D).eq.0) then
        call save1Ddata(Nr,t)
     end if

! En main evolution loop

  end do

  print *,'-----------------------------'


 end subroutine 
