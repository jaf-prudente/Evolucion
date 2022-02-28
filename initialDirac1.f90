subroutine initialDirac1(Nr,dr,F0,w,Nr_rk4)

!                                              _______________________
!    _______________________-------------------                       `\
!  /:--_                                                                |
! ||[ ] |                                   ___________________________/
! | \__/_________________-------------------                         |
!  |                                                                  |
!  |                              INITIALDIRAC1                       |
!   |                                                                  |
!    |            Este programa soluciona el sistema de ecuaciones      |
!     |            diferenciales para la ecuación de Dirac usando        |
!    |           el método de Runge-Kutta a 4to orden, así mismo        |
!   |            calcula los coeficientes métricos.                    |
!   |                                                                  |
!    |                                                                  |
!     |                                                                  |
!    |                                                                  |
!   |                                              ____________________|_
!   |  ___________________-------------------------                      `\
!   |/`--_                                                                 |
!   ||[ ]||                                            ___________________/
!    \===/___________________--------------------------

    use arrays

    implicit none

    !---------------------------------------------------------------
    ! Declaramos variables útiles
    integer Nr, Nr_rk4, i

    real(8) sexto, tercio, cuarto, medio, cero, uno, dos, ocho, diezseis
    real(8) dr

    real(8) F0, w

    ! Declaramos variables del método
    real(8) k1m, k1s, k1F, k1G
    real(8) k2m, k2s, k2F, k2G
    real(8) k3m, k3s, k3F, k3G
    real(8) k4m, k4s, k4F, k4G

    ! Declaramos variables que guarden a las funciones derivadas
    real(8) dm, ds, dF, dG

    ! Declaramos los arreglos que usaremos
    real(8), allocatable, dimension (:) :: m, s, F, G
    allocate( m(0:Nr_rk4), s(0:Nr_rk4), F(0:Nr_rk4), G(0:Nr_rk4) )

    ! Iniciamos variables
    sexto = 1.0D0/6.0D0
    tercio = 1.0D0/3.0D0
    cuarto = 1.0D0/4.0D0
    medio = 1.0D0/2.0D0
    cero = 0.0D0
    uno = 1.0D0
    dos = 2.0D0
    ocho = 8.0D0
    diezseis = 16.0D0

    !---------------------------------------------------------------
    ! Definimos las condiciones iniciales (con simetrías)
    m(0) = -(dos*tercio)*(F0**2)*(dr**3)*w
    s(0) = uno + tercio*(F0**2)*(dos*dos*w - uno)*dr**2
    F(0) = -F0*dr
    G(0) = tercio*F0*(w - uno)*(dr**2)

    m(1) = -m(0)
    s(1) = s(0)
    F(1) = -F(0)
    G(1) = G(0)

    !---------------------------------------------------------------
    ! Hacemos el método de Runge-Kutta en sí
    do i=1, Nr_rk4-1

        !-------------------------------------------
        k1m = dm(w, r(i), m(i), s(i), F(i), G(i))
        k1s = ds(w, r(i), m(i), s(i), F(i), G(i))
        k1F = dF(w, r(i), m(i), s(i), F(i), G(i))
        k1G = dG(w, r(i), m(i), s(i), F(i), G(i))

        !-------------------------------------------
        k2m = dm(w, r(i) + medio*dr, m(i) + medio*k1m*dr, s(i) + medio*k1s*dr, F(i) + medio*k1F*dr, G(i) + medio*k1G*dr)
        k2s = ds(w, r(i) + medio*dr, m(i) + medio*k1m*dr, s(i) + medio*k1s*dr, F(i) + medio*k1F*dr, G(i) + medio*k1G*dr)
        k2F = dF(w, r(i) + medio*dr, m(i) + medio*k1m*dr, s(i) + medio*k1s*dr, F(i) + medio*k1F*dr, G(i) + medio*k1G*dr)
        k2G = dG(w, r(i) + medio*dr, m(i) + medio*k1m*dr, s(i) + medio*k1s*dr, F(i) + medio*k1F*dr, G(i) + medio*k1G*dr)

        !-------------------------------------------
        k3m = dm(w, r(i) + medio*dr, m(i) + medio*k2m*dr, s(i) + medio*k2s*dr, F(i) + medio*k2F*dr, G(i) + medio*k2G*dr)
        k3s = ds(w, r(i) + medio*dr, m(i) + medio*k2m*dr, s(i) + medio*k2s*dr, F(i) + medio*k2F*dr, G(i) + medio*k2G*dr)
        k3F = dF(w, r(i) + medio*dr, m(i) + medio*k2m*dr, s(i) + medio*k2s*dr, F(i) + medio*k2F*dr, G(i) + medio*k2G*dr)
        k3G = dG(w, r(i) + medio*dr, m(i) + medio*k2m*dr, s(i) + medio*k2s*dr, F(i) + medio*k2F*dr, G(i) + medio*k2G*dr)

        !-------------------------------------------
        k4m = dm(w, r(i) + dr, m(i) + k3m*dr, s(i) + k3s*dr, F(i) + k3F*dr, G(i) + k3G*dr)
        k4s = ds(w, r(i) + dr, m(i) + k3m*dr, s(i) + k3s*dr, F(i) + k3F*dr, G(i) + k3G*dr)
        k4F = dF(w, r(i) + dr, m(i) + k3m*dr, s(i) + k3s*dr, F(i) + k3F*dr, G(i) + k3G*dr)
        k4G = dG(w, r(i) + dr, m(i) + k3m*dr, s(i) + k3s*dr, F(i) + k3F*dr, G(i) + k3G*dr)

        !-------------------------------------------
        m(i+1) = m(i) + (k1m + dos*k2m + dos*k3m + k4m)*(dr*sexto)
        s(i+1) = s(i) + (k1s + dos*k2s + dos*k3s + k4s)*(dr*sexto)
        F(i+1) = F(i) + (k1F + dos*k2F + dos*k3F + k4F)*(dr*sexto)
        G(i+1) = G(i) + (k1G + dos*k2G + dos*k3G + k4G)*(dr*sexto)

    end do

    ! Hacemos los arreglos globlales y extrapolamos
    do i=0, Nr_rk4
        F1(i) = F(i)
        F2(i) = 0
        G1(i) = 0
        G2(i) = G(i)
    end do

    do i=Nr_rk4+1, Nr 
        m(i) = m(Nr_rk4)
        s(i) = s(Nr_rk4)
        F1(i) = F(Nr_rk4)
        F2(i) = 0
        G1(i) = G(Nr_rk4)
        G2(i) = 0
    end do

    !---------------------------------------------------------------
    ! Encontremos los coeficientes métricos, primero a.

    do i=0, Nr
        a(i) = uno/sqrt(uno - (2*m(i)/r(i)) )
        alpha(i) = s(i)*sqrt(uno - (2*m(i)/r(i)) )
    end do

end subroutine

!---------------------------------------------------------------
! Definimos funciones para las derivadas.

!-------------------------------------------
function dm(w, r, m, s, F, G)

    implicit none
    real(8) w, r, m, s, F, G, N
    real(8) uno, dos, ocho
    real(8) dm

    uno = 1.0D0
    dos = 2.0D0
    ocho = 8.0D0

    N = uno - (dos*m)/(r)

    dm = (dos*w*( F**2 + G**2 ))/(s)

    return

end function

!-------------------------------------------
function ds(w, r, m, s, F, G)

    implicit none
    real(8) w, r, m, s, F, G, N
    real(8) uno, dos, diezseis
    real(8) ds

    uno = 1.0D0
    dos = 2.0D0
    diezseis = 16.0D0

    N = uno - (dos*m)/(r)

    ds = ( -dos*dos*F*G*sqrt(N)*s + dos*r*G**2*( dos*w + sqrt(N)*s ) + F**2*( dos*dos*r*w - dos*r*sqrt(N)*s ) )/( r*(r-dos*m) )

    return

end function

!-------------------------------------------
function dF(w, r, m, s, F, G)

    implicit none
    real(8) w, r, m, s, F, G, N
    real(8) uno, dos, ocho, diezseis
    real(8) dF

    uno = 1.0D0
    dos = 2.0D0
    ocho = 8.0D0
    diezseis = 16.0D0

    N = uno - (dos*m)/(r)

    dF = ( ( -m + (r + r*F**2 + dos*F*G - r*G**2)*sqrt(N) )/( r*(r-2*m) ) )*F - ( (1)/(sqrt(N)) + (r*w)/(r*s - dos*m*s) )*G

    return

end function

!-------------------------------------------
function dG(w, r, m, s, F, G)

    implicit none
    real(8) w, r, m, s, F, G, N
    real(8) uno, dos, ocho, diezseis
    real(8) dG

    uno = 1.0D0
    dos = 2.0D0
    ocho = 8.0D0
    diezseis = 16.0D0

    N = uno - (dos*m)/(r)

    dG = -( ( m + (r - r*F**2 - dos*F*G + r*G**2)*sqrt(N) )/( r*(r-2*m) ) )*G + ( (r*w - r*sqrt(N)*s)/( s*(r - dos*m) ) )*F

    return

end function