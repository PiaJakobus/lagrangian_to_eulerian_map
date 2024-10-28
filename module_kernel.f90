MODULE KERNEL
  USE CONSTANTS
  IMPLICIT NONE
  CONTAINS

  SUBROUTINE kernel_WC6(w, h, r2, N)

    !***************************************
    !                                      * 
    ! WC6 Kernel                           * 
    !                                      * 
    !***************************************

    INTEGER, INTENT(IN) :: N
    INTEGER :: p, zaehler
    DOUBLE PRECISION, INTENT(IN) :: h, r2(N)
    DOUBLE PRECISION, INTENT(INOUT) :: w(N)
    DOUBLE PRECISION :: dr, q
    zaehler = 0
    DO p = 1, N
     dr = SQRT(r2(p))
     q = dr / h
     IF ((q < 1.0d0) .AND. (q >= 0.0d0))  THEN
        zaehler = zaehler + 1 
        w(p) = (1365d0/(64d0*pi))*(1d0-q)**8d0*(1d0+8d0*q+25d0*q**2d0+32d0*q**3d0)
      ELSE IF (q >= 1d0) THEN
        w(p) = 0d0
      ELSE
        print*, "kernel error"
      END IF
      !print*, "========", dr, h, q, w(p), zaehler, "========"
    END DO 
    print*, "particles: ", zaehler
    print*, "====================="
  END SUBROUTINE kernel_WC6

  SUBROUTINE  kernel_M4(w, h, r2, N)
    !***************************************
    !                                      * 
    ! M4 Kernel                            * 
    !                                      * 
    !***************************************
    INTEGER, INTENT(IN) :: N
    INTEGER :: p
    DOUBLE PRECISION, INTENT(IN) :: h, r2(N)
    DOUBLE PRECISION, INTENT(INOUT) :: w(N)
    DOUBLE PRECISION :: dr,q, const
    const = 1365d0 / (512d0 * pi * h**3d0)
    DO p = 1, N
     dr = SQRT(r2(p))
      q = dr / h
      IF ((q > 0d0) .AND. (q <= 2d0)) THEN
        w(p) = const * (1d0-0.5d0*q)**8d0*(4d0*q**3d0+6.25d0*q**2d0+3d0*q+1d0)
      ELSE IF (q > 2d0) THEN 
        w(p) = 0d0
      ELSE 
        print*, "ELSE KERNEL"
        stop
      END IF
    END DO
  END SUBROUTINE kernel_M4


END MODULE KERNEL
