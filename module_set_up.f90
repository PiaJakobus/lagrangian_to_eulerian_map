MODULE SET_UP
  USE NR
  IMPLICIT NONE 
  CONTAINS

  SUBROUTINE particle_distr(rp,fp,rdim,N)
    !****************************************
    !                                       * 
    ! setup spherical particle distribution * 
    !                                       * 
    !****************************************

    INTEGER, INTENT(IN) :: rdim, N
    DOUBLE PRECISION, INTENT(IN) :: rp(rdim,N)
    DOUBLE PRECISION, INTENT(INOUT) :: fp(N)
    DOUBLE PRECISION :: r2
    INTEGER :: i, j
    DO j = 1, N
      r2 = 0d0
      DO i = 1, rdim
        r2 = r2 + rp(i,j)**2d0
      END DO
      !fp(j) = EXP(-r2)
      fp(j) = 10d0 - SQRT(r2)
    END DO
  END SUBROUTINE particle_distr
  

  SUBROUTINE get_h(r2,r_particles,r, p_max,h,rdim, N)
    DOUBLE PRECISION, INTENT(IN) :: r_particles(rdim,N), r(rdim)
    INTEGER, INTENT(IN) :: rdim, N, p_max
    DOUBLE PRECISION, INTENT(INOUT) :: r2(N), h
    INTEGER :: index_p_max, j, index_arr(N)
    DOUBLE PRECISION :: delta_r(rdim)
    DO j = 1, N
       delta_r = ABS(r_particles(:,j)-r)
       r2(j) = delta_r(1)**2d0 + delta_r(2)**2d0 + delta_r(3)**2d0
    END DO 
    CALL indexx(N, r2, index_arr)
    index_p_max = index_arr(p_max)
    h = SQRT(r2(index_p_max))
  END SUBROUTINE get_h



   SUBROUTINE generate_particles(rp,fp,rdim,N,border,ni,fac)

    !***************************************
    !                                      * 
    ! Put particles on a grid with         * 
    ! some noise in [0,1] fac              * 
    !                                      *
    !***************************************
 
    INTEGER, INTENT(IN) :: rdim, N,ni
    INTEGER :: i,j,k, p
    DOUBLE PRECISION, INTENT(INOUT) :: rp(rdim,N), fp(N)
    DOUBLE PRECISION :: range_(ni), rand(rdim)
    DOUBLE PRECISION, INTENT(IN) :: border, fac
    p = 0
    range_ = linspace(-border,border,ni)
    CALL RANDOM_NUMBER(rand)
    rand = (rand - 0.5d0*border) 
    DO i = 1, ni
      DO j = 1, ni
        DO k = 1, ni
        !IF (SQRT(range_(i)**2d0 + range_(j)**2d0 + range_(k)**2d0) <= border) THEN 
          p = p + 1
          rp(1,p) = range_(i) + fac * rand(1)
          rp(2,p) = range_(j) + fac * rand(2)
          rp(3,p) = range_(k) + fac * rand(3)
          !print*, range_
        !END IF  
        END DO
      END DO
    END DO
    !CALL RANDOM_NUMBER(rp)
    !rp = (rp - 0.5d0) * 2d0
    CALL particle_distr(rp,fp,rdim,N)
  END SUBROUTINE generate_particles

  SUBROUTINE generate_grid(x_grid,rdim,N_grid,border,ni_grid)
    
    !***************************************
    !                                      * 
    ! Setup cartesian grid                 * 
    !                                      * 
    !***************************************
    
    INTEGER, INTENT(IN) :: rdim, N_grid, ni_grid
    DOUBLE PRECISION, INTENT(IN) :: border
    DOUBLE PRECISION :: x_grid(rdim,N_grid), range_(ni_grid)
    INTEGER :: i, j, k, ig
    range_ = linspace(-border,border,ni_grid)
    IF (ni_grid == 1) THEN
      range_ = (/0/)
    END IF 
    ig = 0 
    DO i = 1, ni_grid
      DO j = 1,ni_grid
        DO k = 1,ni_grid
          ig = ig + 1
          x_grid(1,ig) = range_(i)!i * dx
          x_grid(2,ig) = range_(j)!j * dx
          x_grid(3,ig) = range_(k)!k * dx
        END DO
      END DO
    END DO
  END SUBROUTINE generate_grid

END MODULE SET_UP
