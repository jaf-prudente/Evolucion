
  subroutine sourcesDirac(Nr,dr)

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
    ! |   SOURCES FOR F 1     |
    ! |           &           |
    ! |   SOURCES FOR F 2     |
    !  -------------------------
    
    ! Son las ecuaciones de evolución de Daka, en resumen las derivadas temporales de F


    ! AQUÍ VA 35 DE DAKA

    
    !  --------------------------------
    ! |   SOURCES FOR F's and G's      |
    !  --------------------------------
    
    ! Since the source term now involves spatial
    ! derivatives, it can not be calculated at
    ! the boundaries.
    
      do i=1,Nr-1
         sF1(i) = -alpha(i)/r(i)*(G1(i)-r(i)*F2(i))& 
                  -sqrt(alpha(i)/a(i))*idrh*(sqrt(alpha(i+1)/a(i+1))*G1(i+1)&
                  -sqrt(alpha(i-1)/a(i-1))*G1(i-1))

         sF2(i) = -alpha(i)/r(i)*(G2(i)+r(i)*F1(i))& 
                  -sqrt(alpha(i)/a(i))*idrh*(sqrt(alpha(i+1)/a(i+1))*G2(i+1)&
                  -sqrt(alpha(i-1)/a(i-1))*G2(i-1))


         sG1(i) = alpha(i)/r(i)*(F1(i)-r(i)*G2(i))& 
                  -sqrt(alpha(i)/a(i))*idrh*(sqrt(alpha(i+1)/a(i+1))*F1(i+1)&
                  -sqrt(alpha(i-1)/a(i-1))*F1(i-1))

         sG2(i) = alpha(i)/r(i)*(F2(i)+r(i)*G1(i))& 
                  -sqrt(alpha(i)/a(i))*idrh*(sqrt(alpha(i+1)/a(i+1))*F2(i+1)&
                  -sqrt(alpha(i-1)/a(i-1))*F2(i-1))

      end do
    
!  =============== 
!  ***   END   ***
!  ===============
  
  end subroutine
    
    
