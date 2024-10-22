MODULE set_up
  USE NR
  USE MLS
  ! schreibe w um!!
  ! schreibe header
  ! Schreibe 05d0, etc
  !******************************
  ! Generate particles in a box
  ! Generate grid
  !******************************

  IMPLICIT NONE
  CONTAINS

  SUBROUTINE particle_distr(rp,fp,rdim,N)
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
      fp(j) = EXP(-r2)
    END DO
  END SUBROUTINE particle_distr

  SUBROUTINE kernel(w,rdim,h,r,rp)
    INTEGER, INTENT(IN) :: rdim
    INTEGER :: i
    DOUBLE PRECISION, INTENT(IN) :: h, r(rdim), rp(rdim)
    DOUBLE PRECISION, INTENT(INOUT) :: w
    DOUBLE PRECISION :: dr, dr3(rdim),q
    dr3 = r - rp
    dr = 0d0
    DO i = 1, rdim
      dr = dr + dr3(i)**2d0
    END DO
    dr = SQRT(dr)
    q = dr / h
    IF (q < 1.0) THEN
      w = 1d0 - 2.5d0 * q**2d0 + (3d0/2d0) * q**3d0
    ELSE IF ((q .GE. 1d0) .AND. (q < 2d0)) THEN
      w = 0.5d0 * (2d0 - q)**2d0 * (1d0 - q)
    ELSE
      w = 0d0
    END IF
  END SUBROUTINE kernel


  SUBROUTINE generate_particles(rp,fp,rdim,N,border,ni)
    INTEGER, INTENT(IN) :: rdim, N,ni
    INTEGER :: i,j,k, p
    DOUBLE PRECISION, INTENT(INOUT) :: rp(rdim,N), fp(N)
    DOUBLE PRECISION :: range_(ni)
    DOUBLE PRECISION, INTENT(IN) :: border
    p = 0
    range_ = linspace(-border,border,ni)
      DO i = 1, ni
        DO j = 1,ni
          DO k = 1,ni
            p = p + 1
            rp(1,p) = range_(i)!i * dx
            rp(2,p) = range_(j)!j * dx
            rp(3,p) = range_(k)!k * dx
          END DO
        END DO
      END DO
    !CALL RANDOM_NUMBER(rp)
    !rp = (rp - 0.5d0) * 2d0
    CALL particle_distr(rp,fp,rdim,N)
  END SUBROUTINE generate_particles

  SUBROUTINE generate_grid(x_grid,rdim,N_grid,border,ni_grid)
    INTEGER, INTENT(IN) :: rdim, N_grid, ni_grid
    DOUBLE PRECISION, INTENT(IN) :: border
    DOUBLE PRECISION :: x_grid(rdim,N_grid), range_(ni_grid)
    INTEGER :: i, j, k, ig
    range_ = linspace(-border,border,ni_grid)
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

   SUBROUTINE assemble_moment_matrix(N,NN3D_lin,moment_matrix,P_i,w,rp,r,rdim,h)
    INTEGER, INTENT(IN) :: NN3D_lin, N, rdim
    DOUBLE PRECISION, INTENT(IN) :: h, rp(rdim, N), r(rdim)
    DOUBLE PRECISION, INTENT(INOUT) :: moment_matrix(NN3D_lin, NN3D_lin), w, P_i(NN3D_lin)
    INTEGER :: i,p,k
    moment_matrix = 0d0
    DO p = 1, N
      CALL get_base3D_lin(rp(:,p),r,P_i)
      DO k = 1, NN3D_lin
        DO i = 1, NN3D_lin
          CALL kernel(w,rdim,h,r,rp(:,p))
          moment_matrix(i,k) = moment_matrix(i,k) + P_i(i) * P_i(k) * w
        END DO
      END DO
    END DO
  END SUBROUTINE assemble_moment_matrix

  SUBROUTINE get_vector_b(B,fp,P_i,N,w,rp,r,rdim,h)
    INTEGER :: k,p
    INTEGER, intent(IN) :: N, rdim
    DOUBLE PRECISION, INTENT(IN) :: rp(rdim,N), r(rdim), fp(N), h
    DOUBLE PRECISION, INTENT(INOUT) :: w, B(NN3D_lin),P_i(NN3D_lin)
    B = 0.
    DO p = 1, N
      CALL get_base3D_lin(r,rp(:,p),P_i)
      CALL kernel(w,rdim,h,r,rp(:,p))
      !print*, w
      !print*, "w", w
      DO k = 1, NN3D_lin
          B(k) = B(k) + fp(p) * P_i(k) * w
        END DO
      END DO
  END SUBROUTINE


  SUBROUTINE algebra(r,r_G,rdim,NN3D_lin,A_inv,moment_matrix,B,beta,f_grid)
    INTEGER, INTENT(IN) :: NN3D_lin, rdim
    DOUBLE PRECISION, INTENT(IN) :: r(rdim), r_G(rdim),moment_matrix(NN3D_lin,NN3D_lin), B(NN3D_lin)
    DOUBLE PRECISION :: pl(NN3D_lin)
    DOUBLE PRECISION, INTENT(INOUT) :: A_inv(NN3D_lin,NN3D_lin), f_grid, beta(NN3D_lin)
    f_grid = 0.
    IF (NN3D_lin == 1) THEN
      A_inv = 1./moment_matrix
    ELSE
      A_inv = matinv4(moment_matrix)
    END IF
    beta = matmul(A_inv,B)
    call get_base3D_lin(r,r_G,pl)
    f_grid = dot_product(beta,pl)
  END SUBROUTINE algebra

  SUBROUTINE file_out(rp,fp,grid,N,rdim,N_grid)
    INTEGER, INTENT(IN) :: rdim, N, N_grid
    INTEGER :: i
    DOUBLE PRECISION, INTENT(IN) :: rp(rdim,N), fp(N), grid(rdim,N_grid)
    CHARACTER(LEN=20) :: file_name1 = "particles.dat", file_name2 = 'grid.dat'
    OPEN(unit=99, file = file_name1, status = 'unknown', action = 'write')
    OPEN(unit=98, file = file_name2, status = 'unknown', action = 'write')
    DO i = 1, N
      WRITE(99,*) i, rp(:,i), fp(i)
    END DO
    DO i = 1, N_grid
      WRITE(98,*) grid(:,i)
    END DO
    CLOSE(98)
    CLOSE(99)
  END SUBROUTINE file_out

END MODULE set_up





PROGRAM main
  USE set_up
  IMPLICIT NONE
  INTEGER, PARAMETER :: ni = 21, N = ni**3, rdim = 3, ni_grid = 5, N_grid = ni_grid**3
  DOUBLE PRECISION, PARAMETER :: border = 2d0, border_g = 1d0, dx = 2d0 * border / (REAL(ni) - 1d0)
  DOUBLE PRECISION, PARAMETER :: pi = 3.1415d0!, dx_g = 2d0
  !DOUBLE PRECISION, parameter :: dx = 0.1d0, pi = 3.1415, border = 2d0, border_g = 1d0
  !INTEGER, PARAMETER :: ni = INT(2d0 * border / dx) + 1, N= (INT(2d0 * border / dx) + 1)**3, N_grid = 2, rdim = 3 
  !INTEGER, PARAMETER :: ni_grid = 5
  DOUBLE PRECISION :: rp(rdim,N), r(rdim), r_G(rdim), r_grid(rdim,N_grid), fp(N), f_grid
  DOUBLE PRECISION :: moment_matrix(NN3D_lin,NN3D_lin), B(NN3D_lin), beta(NN3D_lin), A_inv(NN3D_lin,NN3D_lin)
  DOUBLE PRECISION :: x(rdim,N), w, f_analytic(1)
  DOUBLE PRECISION :: P_i(NN3D_lin), h
  INTEGER :: i
  print*, "Number of particles: ", N, dx
  !h = ((4d0/3d0)*pi / 40d0)**(1d0/4d0) * dx
  h = ((4d0/3d0)*pi / 40d0)**(1d0/4d0) * dx
  !h = ((4d0/3d0)*pi / 40d0)**(1d0/4d0) * dx_g
  CALL generate_particles(rp,fp,rdim,N,border,ni)
  !CALL generate_grid(r_grid,rdim, N_grid,dx_g)
  CALL generate_grid(r_grid,rdim,N_grid,border_g,ni_grid)
  x(1,1) = 0.4d0
  x(2,1) = 0.4d0
  x(3,1) = -0.4d0
  !r_G = x(:,1)
  !r = r_G
  DO i = 1, N_grid
    r = r_grid(:,i)
    !r = x(:,1) 
    !stop
    CALL get_vector_b(B,fp,P_i,N,w,rp,r,rdim,h)
    CALL assemble_moment_matrix(N,NN3D_lin,moment_matrix,P_i,w,rp,r,rdim,h)
    CALL algebra(r,r_G,rdim,NN3D_lin,A_inv,moment_matrix,B,beta,f_grid)
    CALL file_out(rp,fp,r_grid,N,rdim,N_grid)
    CALL particle_distr(r,f_analytic,rdim,1)
    print*, "radius", r, "f_grid", f_grid, "f_analytic", f_analytic, "error: ", ABS(f_grid - f_analytic) / f_analytic
  END DO 
END PROGRAM main
