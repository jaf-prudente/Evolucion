
  subroutine potential(Nr,dr)

! --->   POTENTIALS FOR SCALAR FIELDS   <---

! This routine calculates the potential and its
! derivative with respect to the field at every point
! along the radial coordinate.

  use arrays

  implicit none

  common / numbers / zero,half,one,two,kappa0
  common / Vparam  / lambda,gamma,epsilon,Vtype
  common / Vparam3 / lambda3,epsilon3,Vtype3

  integer i,Nr
  integer Vtype,Vtype3

  real(kind=8) dr
  real(kind=8) lambda,gamma,epsilon
  real(kind=8) lambda3,epsilon3
  real(kind=8) zero,half,one,two,kappa0


! --->  POTENTIAL for PHI

  do i=0,Nr

! The modulus of the wave function

     phi(i) = sqrt(phi1(i)**2 + phi2(i)**2)

!    No potential.

     if (Vtype.eq.0) then

        V(i)   = zero
        VP1(i) = zero
        VP2(i) = zero

!    Quadratic.

     else if (Vtype.eq.1) then

        V(i)   = lambda*(phi1(i)**2 + phi2(i)**2)
        VP1(i) = lambda*phi1(i)
        VP2(i) = lambda*phi2(i)

!    Quartic.

     else if (Vtype.eq.2) then

        V(i)   = lambda*(phi1(i)**2 + phi2(i)**2) + &
                 epsilon*(phi1(i)**2 + phi2(i)**2)**2
        VP1(i) = lambda*phi1(i) + two*epsilon*(phi1(i)**2 &
&+ phi2(i)**2)*phi1(i)
        VP2(i) = lambda*phi2(i) + two*epsilon*(phi1(i)**2 &
&+ phi2(i)**2)*phi2(i)


     end if

  end do



! --->   END   <---

  end subroutine potential


