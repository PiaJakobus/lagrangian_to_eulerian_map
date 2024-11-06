MODULE KERNEL

  !**********************************************
  !                                             *
  ! Module KERNEL                               *
  !                                             *
  ! This module contains several kernel         *
  ! functions                                   *
  !                                             *
  ! Included subroutines:                       *
  !   - get_h: Calculates the smoothing length  *
  !     to match a given number of nearest      *
  !     particles as neighbors.                 *
  !   - Wendland_C6_3D: Implements the Wendland *
  !     kernel "Phi_3,3" with 3D normalization, *
  !     often used for smoothing in SPH methods.*
  !   - kernel_WC6: Computes the WC6 kernel.    *
  !   - kernel_M4: Computes the M4 kernel.      *
  !                                             * 
  !**********************************************
  USE constants, ONLY: pi
  USE NR, ONLY : indexx 
  IMPLICIT NONE
  CONTAINS



SUBROUTINE get_h(r2, r_particles, r, p_max, h, rdim, N)

    !**********************************************
    !                                             *
    ! Subroutine to calculate smoothing length h  *  
    ! to match pmax partices as neighbours        *
    !                                             *
    ! Inputs:                                     *
    !   - r_particles: Array of particle          *
    !     positions, dimensions rdim x N          *
    !   - r: Current position vector with         *
    !     dimension rdim                          *
    !   - rdim: Number of dimensions              *
    !   - N: Number of particles                  *
    !   - p_max: Index for the maximum distance   *
    !     to be considered                        *
    ! Outputs:                                    *
    !   - r2: Array of distances between          *
    !     particles, length N                     *
    !   - h: Calculated distance value            *
    !                                             *
    !**********************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: r_particles(rdim, N), r(rdim)
    INTEGER, INTENT(IN) :: rdim, N, p_max
    DOUBLE PRECISION, INTENT(OUT) :: r2(N), h
    INTEGER :: index_p_max, j, index_arr(N)
    DOUBLE PRECISION :: delta_r(rdim)

    !**********************************************
    ! Loop to calculate distances between         *
    ! particles                                   *
    !**********************************************
    DO j = 1, N
        delta_r = r_particles(:, j) - r
        ! Calculate the difference between the j-th particle position and the current position
        r2(j) = delta_r(1)**2d0 + delta_r(2)**2d0 + delta_r(3)**2d0
        ! Compute the squared distance and store in the array r2
    END DO

    !**********************************************
    ! Sorting distances and selecting the index   *
    ! with p_max                                  *
    !**********************************************
    CALL indexx(N, r2, index_arr)
    ! Call a sorting routine (indexx) to sort the indices based on values in r2
    index_p_max = index_arr(p_max)
    ! Select the index corresponding to the p_max element

    !**********************************************
    ! Calculating h based on the maximum distance *
    !**********************************************
    h = SQRT(r2(index_p_max))
    ! Set h as the square root of the distance of the selected particle

  END SUBROUTINE get_h
  

  PURE SUBROUTINE Wendland_C6_3D(q,W,zaehler)
    
    !**********************************************
    !                                             *
    ! Wendland kernel "Phi_3,3" from Table 5.1 in *
    ! Schaback and Wendland "Kernel techniques "  *
    ! Acta Numerica 15, 543 (2006); with a        *
    ! support of 2*h;                             *
    !                                             *
    ! ATTENTION: this version assumes a           *
    !            3D-normalization !!!             *
    ! SKR 25.08.2014                              *
    !                                             *
    !**********************************************
    
    IMPLICIT NONE
    
    DOUBLE PRECISION,INTENT(IN)  :: q
    DOUBLE PRECISION,INTENT(OUT) :: W
    DOUBLE PRECISION, PARAMETER  :: norm_3D= 1365.D0/(512.D0*Pi)
    INTEGER,INTENT(INOUT) :: zaehler 
    DOUBLE PRECISION                q2,q3,norm,pref,brac

    ! so that support is 2.*h
    !r=    0.5D0*v
    IF(q < 1.D0)THEN
       zaehler = zaehler + 1 
       q2=       q * q 
       q3=       q2*q
       pref=     (1.0D0 - q)**8
       brac=     32.D0*q3 + 25.D0*q2 +  8.D0*q + 1.D0
       W=        pref*brac
    ELSE
       W=        0.D0
    ENDIF
    norm= norm_3D
    W=    norm*W
  END SUBROUTINE Wendland_C6_3D

  SUBROUTINE kernel_WC6(v,w)

    !***************************************
    !                                      * 
    ! WC6 Kernel                           * 
    !                                      * 
    !***************************************

    DOUBLE PRECISION, INTENT(IN) :: v
    DOUBLE PRECISION, INTENT(OUT) :: w
    DOUBLE PRECISION :: q
    q = 0.5d0*v 
    !q = dr / h
    IF (q < 1.0d0)  THEN
      w = (1365d0/(64d0*pi))*(1d0-q)**8d0*(1d0+8d0*q+25d0*q**2d0+32d0*q**3d0)
    ELSE IF (q >= 1d0) THEN
      w = 0d0
    END IF
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
    DOUBLE PRECISION, INTENT(OUT) :: w(N)
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
