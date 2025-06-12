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
  USE CONSTANTS
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
    DOUBLE PRECISION :: range_x(ni),range_y(ni),range_z(ni)

    p = 0
    !range_ = linspace(-border,border, ni)  ! Create evenly spaced points within the boundary
    range_x = linspace(-250.d0,10.d0,ni)
    range_z = linspace(-5.d0,5.d0,ni)
    range_y = linspace(-170.d0,60.d0,ni)

    DO i = 1, ni
      DO j = 1, ni
        DO k = 1, ni
          CALL RANDOM_NUMBER(rand)  ! Generate random numbers
          rand = (rand - 0.5d0 * border)  ! Center random noise around zero
          p = p + 1
          ! Assign particle positions with added noise
          rp(1, p) = range_x(i) + fac * rand(1)
          rp(2, p) = range_y(j) + fac * rand(2)
          rp(3, p) = range_z(k) + fac * rand(3)
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
    DOUBLE PRECISION :: range_x(ni_grid), range_y(ni_grid), range_z(ni_grid)
    INTEGER :: i, j, k, ig  ! Loop indices and grid point counter

    range_ = linspace(-border, border, ni_grid)  ! Create evenly spaced grid points within boundary
    !range_x = linspace(0.d0,30.d0,ni_grid)
    !range_z = linspace(0.d0,30.d0,ni_grid)
    !range_y = linspace(0.d0,30.d0,ni_grid)


    ! If ni_grid is 1, set range to zero (single point)
    IF (ni_grid == 1) THEN
      range_ = (/0/)
    END IF 

    ig = 0  ! Initialize grid point counter
    DO i = 1, ni_grid
      DO j = 1, ni_grid
        DO k = 1, ni_grid
          ig = ig + 1
          ! Assign grid points based on range
          x_grid(1, ig) = range_(k)
          x_grid(2, ig) = range_(j)
          x_grid(3, ig) = range_(i)
        END DO
      END DO
    END DO
  END SUBROUTINE generate_grid

  SUBROUTINE generate_polar(x_grid, rdim, N_grid, border_r, ni_grid_r,ni_grid_th,&
                  &ni_grid_ph, range_r,range_th,range_ph)

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

    INTEGER, INTENT(IN) :: rdim, N_grid, ni_grid_r,ni_grid_ph,ni_grid_th
    DOUBLE PRECISION, INTENT(IN) :: border_r
    DOUBLE PRECISION, INTENT(OUT) :: x_grid(rdim,N_grid)
    DOUBLE PRECISION,INTENT(OUT) :: range_r(ni_grid_r) 
    DOUBLE PRECISION,INTENT(OUT) :: range_th(ni_grid_th) 
    DOUBLE PRECISION,INTENT(OUT) :: range_ph(ni_grid_ph) 
    DOUBLE PRECISION, PARAMETER :: border_th = pi
    DOUBLE PRECISION, PARAMETER :: border_ph = 2.d0 * pi 
    INTEGER :: i, j, k, ig  ! Loop indices and grid point counter
    DOUBLE PRECISION :: dr, dtheta, dphi
    dr = 0.05 !border_r / (ni_grid_r - 1.D0)
    dr = 500.D0 !border_r / (ni_grid_r - 1.D0)
    dphi = 2.D0 * pi / (ni_grid_ph - 1.D0)
    dtheta = pi / (ni_grid_th)

    range_r = logspace(log10(dr), log10(border_r), ni_grid_r)  ! Create evenly spaced grid points within boundary
    range_ph = linspace(-pi+0.017453D0, pi-0.017453D0, ni_grid_ph)  ! Create evenly spaced grid points within boundaARY
    range_th = linspace(0.017453D0, pi-0.017453D0, ni_grid_th)  ! Create evenly spaced grid points within boundary

    ig = 0  ! Initialize grid point counter
    DO i = 1, ni_grid_ph
      DO j = 1, ni_grid_th
        DO k = 1, ni_grid_r
          ig = ig + 1
          ! Assign grid points based on range
          x_grid(1, ig) = range_r(k)
          x_grid(2, ig) = range_th(j)
          x_grid(3, ig) = range_ph(i)
        END DO
      END DO
    END DO
  END SUBROUTINE generate_polar

  SUBROUTINE xy_plane(x_grid, N_grid, border, ni_grid)
    !***************************************
    !                                      * 
    ! Setup Cartesian plane for simulations * 
    ! Inputs:                              *
    !   - rdim: Number of dimensions       *
    !   - N_grid: Number of grid points    *
    !   - border: Boundary for grid        *
    !   - ni_grid: Number of intervals     *
    !     along each axis                  *
    ! Outputs:                             *
    !   - x_grid: Grid points (rdim x N_grid) *
    !***************************************

    INTEGER, INTENT(IN) :: N_grid, ni_grid
    DOUBLE PRECISION, INTENT(IN) :: border
    DOUBLE PRECISION, INTENT(OUT) :: x_grid(2,N_grid)
    DOUBLE PRECISION :: range_(ni_grid)  ! Grid points and range array
    INTEGER :: j, k, ig  ! Loop indices and grid point counter
    range_ = linspace(-border, border, ni_grid)  ! Create evenly spaced grid points within boundary
    ! If ni_grid is 1, set range to zero (single point)
    IF (ni_grid == 1) THEN
      range_ = (/0/)
    END IF 
    ig = 0  ! Initialize grid point counter
    DO j = 1, ni_grid
      DO k = 1, ni_grid
        ig = ig + 1
        ! Assign grid points based on range
        x_grid(1, ig) = range_(k)
        x_grid(2, ig) = range_(j)
      END DO
    END DO
  END SUBROUTINE xy_plane


END MODULE SET_UP

