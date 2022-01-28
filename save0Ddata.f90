
  subroutine save0Ddata(Nr,t,dr)

! --->   SAVE 0D DATA TO FILE   <---

! This routine saves "0D" data to files.
! By 0D data I mean reduced quantities as
! fucntions of time (minimum, maximum, norm).

  use arrays

  implicit none

  logical firstcall

  integer i,Nr

  real*8 t,dr
  real*8 max,min,nm1,nm2
  real*8 aux1,rAH

  character*20 filestatus

  data firstcall / .true. /

  save firstcall


! --->   IS THIS THE FIRST CALL?

! On first call, replace file.

  if (firstcall) then

     firstcall = .false.
     filestatus = 'replace'

! Otherwise, just append to it.

  else

     filestatus = 'old'

  end if


! --->   SAVE DATA   <---

! Radial metric function a.

  max = a(0)
  min = a(0)

  do i=1,Nr
     if (a(i)>max) then
        max  = a(i)
        aux1 = r(i)
     end if
     if (a(i)<min) min = a(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='a.tl',form='formatted',status=filestatus)
     write(1,*) '"t      a_max     a(0)'
  else
     open(1,file='a.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,*) t,max,a(0)
!  write(1,"(4ES14.6)") t,max,min,a(0)
  close(1)

  if (filestatus == 'replace') then
     open(11,file='ra_max.tl',form='formatted',status=filestatus)
     write(11,*) '"ra_max.tl'
  else
     open(11,file='ra_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(11,"(2ES14.6)") t,aux1
  close(11)

! Lapse function alpha.

  max = alpha(0)
  min = alpha(0)

  do i=1,Nr
     if (alpha(i)>max) max = alpha(i)
     if (alpha(i)<min) min = alpha(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='alpha.tl',form='formatted',status=filestatus)
     write(1,*) '"t       alpha_min      alpha_max      alpha(0)'
  else
     open(1,file='alpha.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(4ES14.6)") t,min,max,alpha(0)
  close(1)

! TENEMOS QUE HACER ESTO MISMO PARA F_1, F_2, G_1 Y G_2.
! Scalar field phi [the modulus of the wave function].

  max = phi(0)
  min = phi(0)

  do i=1,Nr
     if (phi(i)>max) max = phi(i)
     if (phi(i)<min) min = phi(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='phi_max.tl',form='formatted',status=filestatus)
     write(1,*) '"phi_max.tl'
  else
     open(1,file='phi_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

! Real part of the scalar field (phi1).

  max = phi1(0)
  min = phi1(0)

  do i=1,Nr
     if (phi1(i)>max) max = phi1(i)
     if (phi1(i)<min) min = phi1(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='phi1_max.tl',form='formatted',status=filestatus)
     write(1,*) '"phi1_max.tl'
  else
     open(1,file='phi1_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

  if (filestatus == 'replace') then
     open(1,file='phi1_0.tl',form='formatted',status=filestatus)
     write(1,*) '"phi1_0.tl'
  else
     open(1,file='phi1_0.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,phi1(0)
  close(1)

  if (filestatus == 'replace') then
     open(1,file='phi2_0.tl',form='formatted',status=filestatus)
     write(1,*) '"phi2_0.tl'
  else
     open(1,file='phi2_0.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,phi2(0)
  close(1)

! Imaginari part of the scalar field (phi2).

  max = phi2(0)
  min = phi2(0)

  do i=1,Nr
     if (phi2(i)>max) max = phi2(i)
     if (phi2(i)<min) min = phi2(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='phi2_max.tl',form='formatted',status=filestatus)
     write(1,*) '"phi2_max.tl'
  else
     open(1,file='phi2_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

!  if (filestatus == 'replace') then
!     open(1,file='phi_min.tl',form='formatted',status=filestatus)
!     write(1,*) '"phi_min.tl'
!  else
!     open(1,file='phi_min.tl',form='formatted',status=filestatus, &
!     position='append')
!  end if
!  write(1,"(2ES14.6)") t,min
!  close(1)

! Radial derivative of phi (psi1).

  max = psi1(0)
  min = psi1(0)

  do i=1,Nr
     if (psi1(i)>max) max = psi1(i)
     if (psi1(i)<min) min = psi1(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='psi1_max.tl',form='formatted',status=filestatus)
     write(1,*) '"psi1_max.tl'
  else
     open(1,file='psi1_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

!  if (filestatus == 'replace') then
!     open(1,file='psi_min.tl',form='formatted',status=filestatus)
!     write(1,*) '"psi_min.tl'
!  else
!     open(1,file='psi_min.tl',form='formatted',status=filestatus, &
!     position='append')
!  end if
!  write(1,"(2ES14.6)") t,min
!  close(1)

! Radial derivative of phi (psi2).

  max = psi2(0)
  min = psi2(0)

  do i=1,Nr
     if (psi2(i)>max) max = psi2(i)
     if (psi2(i)<min) min = psi2(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='psi2_max.tl',form='formatted',status=filestatus)
     write(1,*) '"psi2_max.tl'
  else
     open(1,file='psi2_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

!  if (filestatus == 'replace') then
!     open(1,file='psi_min.tl',form='formatted',status=filestatus)
!     write(1,*) '"psi_min.tl'
!  else
!     open(1,file='psi_min.tl',form='formatted',status=filestatus, &
!     position='append')
!  end if
!  write(1,"(2ES14.6)") t,min
!  close(1)

! Radial derivative of phi (psi3).

  max = psi3(0)
  min = psi3(0)

  do i=1,Nr
     if (psi3(i)>max) max = psi3(i)
     if (psi3(i)<min) min = psi3(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='psi3_max.tl',form='formatted',status=filestatus)
     write(1,*) '"psi3_max.tl'
  else
     open(1,file='psi3_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

!  if (filestatus == 'replace') then
!     open(1,file='psi_min.tl',form='formatted',status=filestatus)
!     write(1,*) '"psi_min.tl'
!  else
!     open(1,file='psi_min.tl',form='formatted',status=filestatus, &
!     position='append')
!  end if
!  write(1,"(2ES14.6)") t,min
!  close(1)


!--------------------------------------
! Time derivative of phi (pi).

  max = pi1(0)
  min = pi1(0)

  do i=1,Nr
     if (pi1(i)>max) max = pi1(i)
     if (pi1(i)<min) min = pi1(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='pi1_max.tl',form='formatted',status=filestatus)
     write(1,*) '"pi1_max.tl'
  else
     open(1,file='pi1_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

  max = pi2(0)
  min = pi2(0)

  do i=1,Nr
     if (pi2(i)>max) max = pi2(i)
     if (pi2(i)<min) min = pi2(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='pi2_max.tl',form='formatted',status=filestatus)
     write(1,*) '"pi2_max.tl'
  else
     open(1,file='pi2_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

!  if (filestatus == 'replace') then
!     open(1,file='pi_min.tl',form='formatted',status=filestatus)
!     write(1,*) '"pi_min.tl'
!  else
!     open(1,file='pi_min.tl',form='formatted',status=filestatus, &
!     position='append')
!  end if
!  write(1,"(2ES14.6)") t,min
!  close(1)

! ---------------------------------------------------------

! Energy density (for the BOSON STAR).

  max = rho(0)
  min = rho(0)

  do i=1,Nr
     if (rho(i)>max) max = rho(i)
     if (rho(i)<min) min = rho(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='rho_max.tl',form='formatted',status=filestatus)
     write(1,*) '"rho_max.tl'
  else
     open(1,file='rho_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

! Mass.

  max = mass(0)
  min = mass(0)

  do i=1,Nr
     if (mass(i)>max) max = mass(i)
     if (mass(i)<min) min = mass(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='mass_max.tl',form='formatted',status=filestatus)
     write(1,*) '"mass_max.tl'
  else
     open(1,file='mass_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

! ADM Mass.

  max = ADMmass(0)
  min = ADMmass(0)

  do i=1,Nr
     if (ADMmass(i)>max) max = ADMmass(i)
     if (ADMmass(i)<min) min = ADMmass(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='ADMmass_max.tl',form='formatted',status=filestatus)
     write(1,*) '"ADMmass_max.tl'
  else
     open(1,file='ADMmass_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

! Number of particles.

  max = NofP(0)
  min = NofP(0)

  do i=1,Nr
     if (NofP(i)>max) max = NofP(i)
     if (NofP(i)<min) min = NofP(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='NofP_max.tl',form='formatted',status=filestatus)
     write(1,*) '"NofP_max.tl'
  else
     open(1,file='NofP_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)


!  if (filestatus == 'replace') then
!     open(1,file='mass_min.tl',form='formatted',status=filestatus)
!     write(1,*) '"mass_min.tl'
!  else
!     open(1,file='mass_min.tl',form='formatted',status=filestatus, &
!     position='append')
!  end if
!  write(1,"(2ES14.6)") t,min
!  close(1)

! -----
! Energy density (for the PERTURBATION).

  max = rho3(0)
  min = rho3(0)

  do i=1,Nr
     if (rho3(i)>max) max = rho3(i)
     if (rho3(i)<min) min = rho3(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='rho3_max.tl',form='formatted',status=filestatus)
     write(1,*) '"rho3_max.tl'
  else
     open(1,file='rho3_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

! Mass.

  max = mass3(0)
  min = mass3(0)

  do i=1,Nr
     if (mass3(i)>max) max = mass3(i)
     if (mass3(i)<min) min = mass3(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='mass3_max.tl',form='formatted',status=filestatus)
     write(1,*) '"mass3_max.tl'
  else
     open(1,file='mass3_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

!  if (filestatus == 'replace') then
!     open(1,file='mass_min.tl',form='formatted',status=filestatus)
!     write(1,*) '"mass_min.tl'
!  else
!     open(1,file='mass_min.tl',form='formatted',status=filestatus, &
!     position='append')
!  end if
!  write(1,"(2ES14.6)") t,min
!  close(1)

! Energy density (TOTAL).

  max = rho_total(0)
  min = rho_total(0)

  do i=1,Nr
     if (rho_total(i)>max) max = rho_total(i)
     if (rho_total(i)<min) min = rho_total(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='rho_total_max.tl',form='formatted',status=filestatus)
     write(1,*) '"rho_total_max.tl'
  else
     open(1,file='rho_total_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

! Mass.

  max = mass_total(0)
  min = mass_total(0)

  do i=1,Nr
     if (mass_total(i)>max) max = mass_total(i)
     if (mass_total(i)<min) min = mass_total(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='mass_total_max.tl',form='formatted',status=filestatus)
     write(1,*) '"mass_total_max.tl'
  else
     open(1,file='mass_total_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

! Anisotropy.

  max = anisotropy(0)
  min = anisotropy(0)

  do i=1,Nr
     if (anisotropy(i)>max) max = anisotropy(i)
     if (anisotropy(i)<min) min = anisotropy(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='anisotropy_max.tl',form='formatted',status=filestatus)
     write(1,*) '"anisotropy_max.tl'
  else
     open(1,file='anisotropy_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

! Pressure

  max = p(0)
  min = p(0)

  do i=1,Nr
     if (p(i)>max) max = p(i)
     if (p(i)<min) min = p(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='p_max.tl',form='formatted',status=filestatus)
     write(1,*) '"p_max.tl'
  else
     open(1,file='p_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)


! Tangential pressure

  max = p_t(0)
  min = p_t(0)

  do i=1,Nr
     if (p_t(i)>max) max = p_t(i)
     if (p_t(i)<min) min = p_t(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='p_t_max.tl',form='formatted',status=filestatus)
     write(1,*) '"p_t_max.tl'
  else
     open(1,file='p_t_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)


!  if (filestatus == 'replace') then
!     open(1,file='mass_min.tl',form='formatted',status=filestatus)
!     write(1,*) '"mass_min.tl'
!  else
!     open(1,file='mass_min.tl',form='formatted',status=filestatus, &
!     position='append')
!  end if
!  write(1,"(2ES14.6)") t,min
!  close(1)

! ---------------------------------------------------------

! Potential.

!  max = V(0)
!  min = V(0)
!
!  do i=1,Nr
!     if (V(i)>max) max = V(i)
!     if (V(i)<min) min = V(i)
!  end do
!
!  if (filestatus == 'replace') then
!     open(1,file='V_max.tl',form='formatted',status=filestatus)
!     write(1,*) '"V_max.tl'
!  else
!     open(1,file='V_max.tl',form='formatted',status=filestatus, &
!     position='append')
!  end if
!  write(1,"(2ES14.6)") t,max
!  close(1)

!  if (filestatus == 'replace') then
!     open(1,file='V_min.tl',form='formatted',status=filestatus)
!     write(1,*) '"V_min.tl'
!  else
!     open(1,file='V_min.tl',form='formatted',status=filestatus, &
!     position='append')
!  end if
!  write(1,"(2ES14.6)") t,min
!  close(1)

! Extrinsic curvature (krr)

  max = krr(0)
  min = krr(0)

  do i=1,Nr
     nm1 = nm1 + 0.5D0*(dabs(krr(i-1)) + dabs(krr(i)))*dr
     nm2 = nm2 + 0.5D0*(krr(i-1)**2 + krr(i)**2)*dr
     if (krr(i)>max) max = krr(i)
     if (krr(i)<min) min = krr(i)
  end do

  nm2 = dsqrt(nm2)

  if (filestatus == 'replace') then
     open(1,file='krr_max.tl',form='formatted',status=filestatus)
     write(1,*) '"krr_max.tl'
  else
     open(1,file='krr_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(3ES14.6)") t,krr(0),max
  close(1)

! Constraint adot.

  nm1 = 0.0D0
  nm2 = 0.0D0

  max = adot(0)
  min = adot(0)

  do i=1,Nr
     nm1 = nm1 + 0.5D0*(dabs(adot(i-1)) + dabs(adot(i)))*dr
     nm2 = nm2 + 0.5D0*(adot(i-1)**2 + adot(i)**2)*dr
     if (adot(i)>max) max = adot(i)
     if (adot(i)<min) min = adot(i)
  end do

  nm2 = dsqrt(nm2)

  if (filestatus == 'replace') then
     open(1,file='adot_max.tl',form='formatted',status=filestatus)
     write(1,*) '"adot_max.tl'
  else
     open(1,file='adot_max.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,max
  close(1)

!  if (filestatus == 'replace') then
!     open(1,file='adot_min.tl',form='formatted',status=filestatus)
!     write(1,*) '"adot_min.tl'
!  else
!     open(1,file='adot_min.tl',form='formatted',status=filestatus, &
!     position='append')
!  end if
!  write(1,"(2ES14.6)") t,min
!  close(1)

  if (filestatus == 'replace') then
     open(1,file='adot_nm1.tl',form='formatted',status=filestatus)
     write(1,*) '"adot_nm1.tl'
  else
     open(1,file='adot_nm1.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,nm1
  close(1)

  if (filestatus == 'replace') then
     open(1,file='adot_nm2.tl',form='formatted',status=filestatus)
     write(1,*) '"adot_nm2.tl'
  else
     open(1,file='adot_nm2.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,nm2
  close(1)


! Apparent horizon mass

    do i=2,Nr
       if ((2.0D0*ADMmass(Nr)/r(i).ge.0.99).and.(2.0D0*ADMmass(Nr)/r(i).le.1.0D0)) then
           rAH = r(i)
!           print *, 'rAH =',rAH, 'M_AH=',rAH/2.
       end if
    end do


  if (filestatus == 'replace') then
     open(1,file='Mass_AH.tl',form='formatted',status=filestatus)
     write(1,*) '"Radius_AH   Mass_AH.tl'
  else
     open(1,file='Mass_AH.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(3ES14.6)") t,rAH/2.0D0,rAH
  close(1)

! RicciScalar.

  max = RicciScalar(0)
  min = RicciScalar(0)

  do i=1,Nr
     if (RicciScalar(i)>max) max = RicciScalar(i)
!     if (ADMmass(i)<min) min = ADMmass(i)
  end do

  if (filestatus == 'replace') then
     open(1,file='RicciScalar.tl',form='formatted',status=filestatus)
     write(1,*) '#t Central RicciScalar'
  else
     open(1,file='RicciScalar.tl',form='formatted',status=filestatus, &
     position='append')
  end if
  write(1,"(2ES14.6)") t,RicciScalar(1)
  close(1)


! --->   END   <---

  end subroutine save0Ddata
