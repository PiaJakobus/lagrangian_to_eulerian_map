MODULE project_2D

  !**********************************************
  !                                             *
  ! Module project_2D                           *
  !                                             *
  ! This module contains subroutines to project *
  ! 3D particle data onto a 2D plane            *
  !                                             *
  ! The code is adapted from the visualisation  *
  ! technique described in:                     *
  ! Daniel Price 2007, "Splash: An Interactive  *
  ! Visualisation Tool for Smoothed Particle    *
  ! Hydrodynamics Simulations."                 *
  !                                             *
  ! Uses:                                       *
  !   - KERNEL module, specifically the         *
  !     Wendland_C6_3D subroutine for kernel    *
  !     calculations.                           *
  !                                             *
  !**********************************************

  USE KERNEL, ONLY : Wendland_C6_3D
  IMPLICIT NONE
  CONTAINS

  SUBROUTINE integrate_kernel(r_G, rp, h, rdim, N, f_grid, fp,m,rho)
    !**********************************************
    !                                             *
    ! Subroutine to integrate the SPH kernel over *
    ! particle data in order to compute the       *
    ! contribution to a 2D grid point.            *
    !                                             *
    ! Inputs:                                     *
    !   - r_G: Coordinates of the 2D grid point   *
    !     (rdim dimensions)                       *
    !   - rp: Particle positions (rdim x N)       *
    !   - h: Smoothing lengths for each particle  *
    !   - rdim: Number of dimensions (usually 3)  *
    !   - N: Number of particles                  *
    !   - fp: Scalar field values for particles   *
    ! Outputs:                                    *
    !   - f_grid: Computed field value for the    *
    !     grid point                              *
    !**********************************************

    INTEGER, INTENT(IN) :: N, rdim
    DOUBLE PRECISION, INTENT(IN) :: r_G(rdim), rp(rdim, N), h(N),fp(N),m(N),rho(N)
    DOUBLE PRECISION, INTENT(OUT) :: f_grid
    DOUBLE PRECISION :: w_tilde, w, pos_ab(rdim), dz, omega
    DOUBLE PRECISION :: r2_ab, norm, q_xy2,r2, h2, r_z, q
    DOUBLE PRECISION, PARAMETER :: R2_kernel = 4.D0  ! Squared kernel radius 
    INTEGER :: p, i, zaehler, num_steps

    num_steps = 100  ! Number of integration steps along the z-axis
    f_grid = 0.D0   ! Initialize the field value for the grid point
    norm = 0.D0

    DO p = 1, N
      pos_ab = ABS(rp(1:2, p) - r_G(1:2))  ! Calculate distance in x-y plane
      r2_ab = pos_ab(1)**2d0 + pos_ab(2)**2d0   
      h2 = h(p) * h(p) 
      w_tilde = 0.D0  
      omega = 0.D0
      q_xy2 = r2_ab / h2  
      IF (q_xy2 <= R2_kernel) THEN 
        ! integration from -SQRT(R^2 - q_xy^2) to +SQRT(R^2 - q_xy^2)
        ! but we integrate from 0 to SQRT(R2 - q_xy^2) times a factor 2
        dz = SQRT(R2_kernel - q_xy2) / num_steps 

        DO i = 0, num_steps
          r_z = i * dz  
          r2 = r2_ab + r_z * r_z   
          q = 0.5D0 * SQRT(r2) / h(p)  

          CALL Wendland_C6_3D(q, w, zaehler)
          w_tilde = w_tilde + w * dz  
        END DO
        w_tilde = 2.D0 * w_tilde ! since we used that integral is symmetric
        omega = w_tilde * m(p) / (rho(p) * h(p)*h2) 
        norm = norm + omega
      END IF 
      !ELSE IF: No contribution to summation in f_grid
      f_grid = f_grid + omega * fp(p)  ! Add weighted particle value to the grid point
    END DO
    f_grid = f_grid / norm 
  END SUBROUTINE

END MODULE project_2D
