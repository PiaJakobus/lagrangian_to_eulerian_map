MODULE SET_UP

  !**********************************************
  !                                             *
  ! Module SET_UP                               *
  !                                             *
  ! This module contains routines to set up     *
  ! particle distributions and grid generation  *
  ! for use in SPH 
  !                                             *
  !   - particle_distr: Sets up a spherical     *
  !     particle distribution.                  *
  !   - generate_particles: Places particles    *
  !     on a grid with random noise.            *
  !   - generate_grid: Generates a Cartesian    *
  !     grid for calculations.                  *
  !                                             *
  !**********************************************

  USE NR
  IMPLICIT NONE 
  CONTAINS

  PURE SUBROUTINE particle_distr(rp, fp, rdim, N)
    !****************************************
    !                                       * 
    ! Setup spherical particle distribution * 
    ! This subroutine calculates a property *
    ! for each particle based on its        *
    ! distance from the origin.             *
    !                                       * 
    ! Inputs:                               *
    !   - rp: Positions of particles (rdim x N)
    !   - rdim: Number of dimensions        *
    !   - N: Number of particles            *
    ! Output:                               *
    !   - fp: Calculated property for each  *
    !     particle (N)                      *
    !****************************************

    INTEGER, INTENT(IN) :: rdim, N
    DOUBLE PRECISION, INTENT(IN) :: rp(rdim, N)  ! Particle positions
    DOUBLE PRECISION, INTENT(OUT) :: fp(N)       ! Property array to be filled
    DOUBLE PRECISION :: r2                       ! Squared distance from the origin
    INTEGER :: i, j                              ! Loop indices

    DO j = 1, N
      r2 = 0d0
      DO i = 1, rdim
        r2 = r2 + rp(i, j)**2d0
      END DO
      ! Calculate property based on distance (example function used here)
      fp(j) = 10d0 - SQRT(r2)
    END DO
  END SUBROUTINE particle_distr

  SUBROUTINE generate_particles(rp, fp, rdim, N, border, ni, fac)
    !***************************************
    !                                      * 
    ! Generate particle positions on a     * 
    ! regular grid with added noise.       *
    !                                      *
    ! Inputs:                              *
    !   - rdim: Number of dimensions       *
    !   - N: Number of particles           *
    !   - border: Boundary for grid        *
    !   - ni: Number of intervals along    *
    !     each axis                        *
    !   - fac: Factor to scale random noise*
    ! Outputs:                             *
    !   - rp: Positions of particles (rdim x N)
    !   - fp: Property values for particles*
    !***************************************

    INTEGER, INTENT(IN) :: rdim, N, ni
    INTEGER :: i, j, k, p                      ! Loop indices and particle counter
    DOUBLE PRECISION, INTENT(OUT) :: rp(rdim, N), fp(N)  ! Particle positions and properties
    DOUBLE PRECISION :: range_(ni), rand(rdim)            ! Grid range and random values
    DOUBLE PRECISION, INTENT(IN) :: border, fac          ! Grid boundary and noise factor

    p = 0
    range_ = linspace(-border, border, ni)  ! Create evenly spaced points within the boundary

    DO i = 1, ni
      DO j = 1, ni
        DO k = 1, ni
          CALL RANDOM_NUMBER(rand)  ! Generate random numbers
          rand = (rand - 0.5d0 * border)  ! Center random noise around zero
          p = p + 1
          ! Assign particle positions with added noise
          rp(1, p) = range_(i) + fac * rand(1)
          rp(2, p) = range_(j) + fac * rand(2)
          rp(3, p) = range_(k) + fac * rand(3)
        END DO
      END DO
    END DO

    ! Call particle distribution subroutine to calculate properties
    CALL particle_distr(rp, fp, rdim, N)
  END SUBROUTINE generate_particles

  SUBROUTINE generate_grid(x_grid, rdim, N_grid, border, ni_grid)
    !***************************************
    !                                      * 
    ! Setup Cartesian grid for simulations * 
    ! Generates grid points within a       *
    ! defined boundary.                    *
    !                                      *
    ! Inputs:                              *
    !   - rdim: Number of dimensions       *
    !   - N_grid: Number of grid points    *
    !   - border: Boundary for grid        *
    !   - ni_grid: Number of intervals     *
    !     along each axis                  *
    ! Outputs:                             *
    !   - x_grid: Grid points (rdim x N_grid) *
    !***************************************

    INTEGER, INTENT(IN) :: rdim, N_grid, ni_grid
    DOUBLE PRECISION, INTENT(IN) :: border
    DOUBLE PRECISION, INTENT(OUT) :: x_grid(rdim,N_grid)
    DOUBLE PRECISION :: range_(ni_grid)  ! Grid points and range array
    INTEGER :: i, j, k, ig  ! Loop indices and grid point counter

    range_ = linspace(-border, border, ni_grid)  ! Create evenly spaced grid points within boundary

    ! If ni_grid is 1, set range to zero (single point)
    IF (ni_grid == 1) THEN
      range_ = (/0/)
    END IF 

    ig = 0  ! Initialize grid point counter
    !DO i = 1, ni_grid
      DO j = 1, ni_grid
        DO k = 1, ni_grid
          ig = ig + 1
          ! Assign grid points based on range
          x_grid(1, ig) = range_(k)
          x_grid(2, ig) = range_(j)
          x_grid(3, ig) = 0.D0!range_(i)
        END DO
      END DO
    !END DO
  END SUBROUTINE generate_grid

END MODULE SET_UP

