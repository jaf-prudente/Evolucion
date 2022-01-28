
  subroutine save1Ddata(Nr,t)

! --->   SAVE 1D DATA TO FILE   <---

! This routine saves 1D data to the files,
! that is, the whole spatial arrays at a
! given time.

  use arrays

  implicit none

  logical firstcall

  integer i,Nr

  real*8 t

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

  if (filestatus=='replace') then
     open(1,file='a.rl',form='formatted',status=filestatus)
  else
     open(1,file='a.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),a(i)
  end do

  write(1,*)

  close(1)

! TENEMOS QUE HACER LO CORRESPONDIENTE A ESTO PERO CON F.
! Lapse function alpha.

  if (filestatus=='replace') then
     open(1,file='alpha.rl',form='formatted',status=filestatus)
  else
     open(1,file='alpha.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),alpha(i)
  end do

  write(1,*)

  close(1)

! Modulus of the wave function phi.

  if (filestatus=='replace') then
     open(1,file='phi.rl',form='formatted',status=filestatus)
  else
     open(1,file='phi.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),phi(i)
  end do

  write(1,*)

! Sclar Field phi1.

  if (filestatus=='replace') then
     open(1,file='phi1.rl',form='formatted',status=filestatus)
  else
     open(1,file='phi1.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),phi1(i)
  end do

  write(1,*)

! Scalar Field phi2.

  if (filestatus=='replace') then
     open(1,file='phi2.rl',form='formatted',status=filestatus)
  else
     open(1,file='phi2.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),phi2(i)
  end do

  write(1,*)

  close(1)

! Scalar Field phi3.

  if (filestatus=='replace') then
     open(1,file='phi3.rl',form='formatted',status=filestatus)
  else
     open(1,file='phi3.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),phi3(i)
  end do

  write(1,*)

  close(1)

! Radial derivative of phi (psi1).

  if (filestatus=='replace') then
     open(1,file='psi1.rl',form='formatted',status=filestatus)
  else
     open(1,file='psi1.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),psi1(i)
  end do

  write(1,*)

  close(1)

! Radial derivative of phi (psi2).

  if (filestatus=='replace') then
     open(1,file='psi2.rl',form='formatted',status=filestatus)
  else
     open(1,file='psi2.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),psi2(i)
  end do

  write(1,*)

  close(1)

! Radial derivative of phi (psi3).

  if (filestatus=='replace') then
     open(1,file='psi3.rl',form='formatted',status=filestatus)
  else
     open(1,file='psi3.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),psi3(i)
  end do

  write(1,*)

  close(1)

! Time derivative of phi (pi1).

  if (filestatus=='replace') then
     open(1,file='pi1.rl',form='formatted',status=filestatus)
  else
     open(1,file='pi1.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),pi1(i)
  end do

  write(1,*)

  close(1)

! Time derivative of phi2 (pi2).

  if (filestatus=='replace') then
     open(1,file='pi2.rl',form='formatted',status=filestatus)
  else
     open(1,file='pi2.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),pi2(i)
  end do

  write(1,*)

  close(1)

! Time derivative of phi (pi3).

  if (filestatus=='replace') then
     open(1,file='pi3.rl',form='formatted',status=filestatus)
  else
     open(1,file='pi3.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),pi3(i)
  end do

  write(1,*)

  close(1)

! -------------------------------------------
! Energy density fo the BOSON STAR (rho).

  if (filestatus=='replace') then
     open(1,file='rho.rl',form='formatted',status=filestatus)
  else
     open(1,file='rho.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),rho(i)
  end do

  write(1,*)

  close(1)

! Energy density fo the PERTURBATION (rho3).

  if (filestatus=='replace') then
     open(1,file='rho3.rl',form='formatted',status=filestatus)
  else
     open(1,file='rho3.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),rho3(i)
  end do

  write(1,*)

  close(1)

! Energy density fo the WHOLE SYSTEM (rho_total).

  if (filestatus=='replace') then
     open(1,file='rho_total.rl',form='formatted',status=filestatus)
  else
     open(1,file='rho_total.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),rho_total(i)
  end do

  write(1,*)

  close(1)

! -----------------------------------------

! Radial pressure for the BOSON STAR (p).

  if (filestatus=='replace') then
     open(1,file='p.rl',form='formatted',status=filestatus)
  else
     open(1,file='p.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),p(i)
  end do

  write(1,*)

  close(1)

! Tangential pressure for the BOSON STAR (p_t).

  if (filestatus=='replace') then
     open(1,file='p_t.rl',form='formatted',status=filestatus)
  else
     open(1,file='p_t.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),p_t(i)
  end do

  write(1,*)

  close(1)

! Radial pressure for the PERTURBATION (p3).

  if (filestatus=='replace') then
     open(1,file='p3.rl',form='formatted',status=filestatus)
  else
     open(1,file='p3.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),p3(i)
  end do

  write(1,*)

  close(1)

! Tangential pressure for the BOSON STAR (p3_t).

  if (filestatus=='replace') then
     open(1,file='p3_t.rl',form='formatted',status=filestatus)
  else
     open(1,file='p3_t.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),p3_t(i)
  end do

  write(1,*)

  close(1)


! --------------------------------------------

! Mass of the BOSON STAR.

  if (filestatus=='replace') then
     open(1,file='mass.rl',form='formatted',status=filestatus)
  else
     open(1,file='mass.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),mass(i)
  end do

  write(1,*)

  close(1)

! Number of particles in the system.

  if (filestatus=='replace') then
     open(1,file='NofP.rl',form='formatted',status=filestatus)
  else
     open(1,file='NofP.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),nofp(i)
  end do

  write(1,*)

  close(1)

! ADM Mass of the BOSON STAR.

  if (filestatus=='replace') then
     open(1,file='ADMmass.rl',form='formatted',status=filestatus)
  else
     open(1,file='ADMmass.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),ADMmass(i)
  end do

  write(1,*)

  close(1)

! Mass of the PERTURBATION.

  if (filestatus=='replace') then
     open(1,file='mass3.rl',form='formatted',status=filestatus)
  else
     open(1,file='mass3.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),mass3(i)
  end do

  write(1,*)

  close(1)

! Mass of the WHOLE SYSTEM.

  if (filestatus=='replace') then
     open(1,file='mass_total.rl',form='formatted',status=filestatus)
  else
     open(1,file='mass_total.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),mass_total(i)
  end do

  write(1,*)

  close(1)

! Omega of the BOSON STAR.

  if (filestatus=='replace') then
     open(1,file='omega.rl',form='formatted',status=filestatus)
  else
     open(1,file='omega.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),omega(i)
  end do

  write(1,*)

  close(1)

! Tangential omega of the BOSON STAR.

  if (filestatus=='replace') then
     open(1,file='omega_t.rl',form='formatted',status=filestatus)
  else
     open(1,file='omega_t.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),omega_t(i)
  end do

  write(1,*)

  close(1)

! Fractional anisotropy of the BOSON STAR

  if (filestatus=='replace') then
     open(1,file='anisotropy.rl',form='formatted',status=filestatus)
  else
     open(1,file='anisotropy.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),anisotropy(i)
  end do

  write(1,*)

  close(1)

! Omega of the PERTURBATION

  if (filestatus=='replace') then
     open(1,file='omega3.rl',form='formatted',status=filestatus)
  else
     open(1,file='omega3.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),omega3(i)
  end do

  write(1,*)

  close(1)

! Potential of the BOSON STAR.

  if (filestatus=='replace') then
     open(1,file='V.rl',form='formatted',status=filestatus)
  else
     open(1,file='V.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),V(i)
  end do

  write(1,*)

  close(1)

! Potential of the PERTURBATION

  if (filestatus=='replace') then
     open(1,file='V3.rl',form='formatted',status=filestatus)
  else
     open(1,file='V3.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),V3(i)
  end do

  write(1,*)

  close(1)

! Extrinsic curvature.

  if (filestatus=='replace') then
     open(1,file='krr.rl',form='formatted',status=filestatus)
  else
     open(1,file='krr.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),krr(i)
  end do

  write(1,*)

  close(1)

! Constraint adot.

  if (filestatus=='replace') then
     open(1,file='adot.rl',form='formatted',status=filestatus)
  else
     open(1,file='adot.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),adot(i)
  end do

  write(1,*)

  close(1)

! RicciScalar.

  if (filestatus=='replace') then
     open(1,file='RicciScalar.rl',form='formatted',status=filestatus)
  else
     open(1,file='RicciScalar.rl',form='formatted',status=filestatus, &
     position='append')
  end if

  write(1,"(A9,ES14.6)") ''
  write(1,"(A9,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(1,"(2ES14.6)") r(i),RicciScalar(i)
  end do

  write(1,*)

  close(1)

! --->   END   <---

  end subroutine save1Ddata

