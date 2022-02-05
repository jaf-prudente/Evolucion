subroutine metricDirac(Nr,dr)

! METRIC COMPONENTS  
! Integrate the equations for the metric

  use arrays

  implicit none

  common / numbers / zero,half,one,two,kappa0

  integer i,Nr

  real(8) dr,idr
  real(8) zero,half,one,two,kappa0
!  real(8) aux, rescale,ssa,rk1,rhoa, rhoa1, rhoa2, rrhoa
!  real(8) ssalpha, aux1, ssa2, aux2

    real(8) aux, ro, Srr, ro_aux, Srr_aux, ro_aux2, Srr_aux2
    real(8) ro_a1, ro_a2, roa, ssa, rk1, aux1, aux2
    real(8) ssalpha

    real(8), allocatable, dimension (:) :: dF1,dF2,dG1,dG2
    allocate(dF1(0:Nr),dF2(0:Nr),dG1(0:Nr),dG2(0:Nr))


! --->   NUMBERS

  idr = 1.0d0/dr


! --->   METRIC {a}


! dF and dG using finite differencing.

    do i=1,Nr-1
      dF1(i) = half*(F1(i+1) - F1(i-1))*idr
      dF2(i) = half*(F2(i+1) - F2(i-1))*idr
      dG1(i) = half*(G1(i+1) - G1(i-1))*idr
      dG2(i) = half*(G2(i+1) - G2(i-1))*idr
    end do

!  One sided derivative at outer boundary.

    dF1(Nr) = half*(3.0d0*F1(Nr) - 4.0d0*F1(Nr-1) + F1(Nr-2))*idr
    dF2(Nr) = half*(3.0d0*F2(Nr) - 4.0d0*F2(Nr-1) + F2(Nr-2))*idr
    dG1(Nr) = half*(3.0d0*G1(Nr) - 4.0d0*G1(Nr-1) + G1(Nr-2))*idr
    dG2(Nr) = half*(3.0d0*G2(Nr) - 4.0d0*G2(Nr-1) + G2(Nr-2))*idr

    dF1(0) = half*(-3.0d0*F1(0) + 4.0d0*F1(1) - F1(2))*idr
    dF2(0) = half*(-3.0d0*F2(0) + 4.0d0*F2(1) - F2(2))*idr
    dG1(0) = half*(-3.0d0*G1(0) + 4.0d0*G1(1) - G1(2))*idr
    dG2(0) = half*(-3.0d0*G2(0) + 4.0d0*G2(1) - G2(2))*idr

!  Regularity at the origin.

    F1(0) = -F1(1)
    F2(0) = -F2(1)
    G1(0) = G1(1)
    G2(0) = G2(1)

! Find the metric coefficients.

    a(0) = one                           
    a(1) = one

    do i=2,Nr

       ro = (1.0/(r(i-1)**2*a(i-1)**2))*( &
            a(i-1)*(F1(i-1)**2 + F2(i-1)**2 - G1(i-1)**2 - G2(i-1)**2) &
            + (2.0*a(i-1)/r(i-1))*(F1(i-1)*G2(i-1) - F2(i-1)*G1(i-1)) &
            + F1(i-1)*dG2(i-1) - F2(i-1)*dG1(i-1) + G1(i-1)*dF2(i-1) - G2(i-1)*dF1(i-1))

       ssa = a(i-1)*(half*(one - a(i-1)**2)/r(i-1) &
                + r(i-1)*a(i-1)**2*ro)

       rk1 = a(i-1) + dr*half*ssa

       ro_a1 = (1.0/(r(i-1)**2*rk1**2))*( &
          rk1*(F1(i-1)**2 + F2(i-1)**2 - G1(i-1)**2 - G2(i-1)**2) &
          + (2.0*rk1/r(i-1))*(F1(i-1)*G2(i-1) - F2(i-1)*G1(i-1)) &
          + F1(i-1)*dG2(i-1) - F2(i-1)*dG1(i-1) + G1(i-1)*dF2(i-1) - G2(i-1)*dF1(i-1))

          ro_a2 = (1.0/(r(i)**2*rk1**2))*( &
          rk1*(F1(i)**2 + F2(i)**2 - G1(i)**2 - G2(i)**2) &
          + (2.0*rk1/r(i))*(F1(i)*G2(i) - F2(i)*G1(i)) &
          + F1(i)*dG2(i) - F2(i)*dG1(i) + G1(i)*dF2(i) - G2(i)*dF1(i))

          roa = half*(ro_a1*r(i-1) + ro_a2*r(i))   
          
          ssa = rk1* ( (one - rk1**2)/(r(i-1) + r(i)) + half*(r(i-1) + r(i))*rk1**2*(roa) )  

          a(i) = a(i-1) + dr*ssa
    end do


!   Compute alpha from outside

    alpha(Nr) = one/a(Nr)

    do i=Nr-1,1,-1
          
        Srr = (1.0/(r(i+1)**2*a(i+1)**2))*(F1(i+1)*dG2(i+1) &
               - F2(i+1)*dG1(i+1) + G1(i+1)*dF2(i+1) - G2(i+1)*dF1(i+1))
               
        aux1 = -half*(one - a(i+1)**2)/r(i+1) &
                + r(i+1)*a(i+1)**2*Srr
           
        ssalpha = alpha(i+1)*aux1       
 
        aux = alpha(i+1) - half*dr*ssalpha

        Srr_aux2 = (1.0/(r(i)**2*aux**2))*(F1(i)*dG2(i) &
                    - F2(i)*dG1(i) + G1(i)*dF2(i) - G2(i)*dF1(i))

        aux2 = -half*(one - a(i)**2)/r(i) &
                + r(i)*a(i)**2*Srr_aux2

        ssalpha = half* (aux1+aux2)

        alpha(i) = alpha(i+1) - aux*(ssalpha)*dr
                                
    end do

    alpha(0) = alpha(1)

! -------------------------
! -------------------------


! --->   END   <---

  end subroutine 
