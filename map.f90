MODULE set_up

! beliebige dimension (INT_param)
! Mit standard SPH ausrechnen und vergleichen 
! Quadru Order 
! Rand effekte mit delta * W^1 + (1-delta) * W^2
! Gaussian noise auf die Teilchen addieren 
! SVD()
  !**********************************************
  !                                             * 
  ! Map particles onto a grid                   * 
  ! Pia Jakobus                                 * 
  ! 10.10.2024                                  * 
  !                                             * 
  !**********************************************
 
  USE NR
  USE MLS
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: pi = 3.1415d0
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



  SUBROUTINE kernel_WC6(w, h, r2, N)

    !***************************************
    !                                      * 
    ! M4 Kernel                            * 
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


  SUBROUTINE generate_particles(rp,fp,rdim,N,border,ni,fac)

    !***************************************
    !                                      * 
    ! Put particles on a grid with         * 
    ! some noise in [0,1] fac              * 
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

   SUBROUTINE assemble_moment_matrix(N,INT_param,moment_matrix,P_i,w,rp,r_G,rdim)

    !********************************************
    !                                           * 
    !   moment matrix: after a set of basis     *
    !   functions bas and a weighting function  *
    !   W have been chosen, we can assemble     * 
    !   the moment matrix:                      *
    !                                           *
    !   M_ij= sum_b W(v_g)*P_i(x_b)*P_j(x_b)    *
    !                                           * 
    !   where v_g= |x_b - x_g|/h and P_k        *
    !   is the kth basis function; this matrix  *
    !   is actually INDEPENDENT of the          *
    !   quantities to be approximated and only  *
    !   dependent on the positions, basis       *
    !   functions and weights                   * 
    !                                           * 
    !********************************************

    INTEGER, INTENT(IN) :: INT_param, N, rdim
    DOUBLE PRECISION, INTENT(IN) :: rp(rdim, N), r_G(rdim), w(N)
    DOUBLE PRECISION, INTENT(INOUT) :: moment_matrix(INT_param, INT_param), P_i(INT_param)
    INTEGER :: i,p,k
    moment_matrix = 0d0
    DO p = 1, N
      IF (INT_param == NN3D_quad) THEN 
        CALL get_base3D_quad(rp(:,p),r_G,P_i)
      ELSE IF (INT_param == NN3D_lin) THEN 
        CALL get_base3D_lin(rp(:,p),r_G,P_i)
      ELSE IF (INT_param == NN3D_cub) THEN
        CALL get_base3D_cubic(rp(:,p),r_G,P_i)
      ELSE IF (INT_param == NN3D_1d) THEN 
        P_i = (/1d0/)
      END IF 
      DO k = 1, INT_param
        DO i = 1, INT_param
          !CALL kernel(w,rdim,h,r_G,rp(:,p))
          moment_matrix(i,k) = moment_matrix(i,k) + P_i(i) * P_i(k) * w(p)
        END DO
      END DO
    END DO
  END SUBROUTINE assemble_moment_matrix

  SUBROUTINE get_vector_b(B,fp,P_i,N,w,rp,r_G,rdim,INT_param)

    !********************************************
    !                                           *
    !    the vector contains the                *
    !    quantities that we actually want       *
    !    to approximate, say A:                 *
    !                                           *
    !    (B)_i= sum_b A_b*W(v_g)*bas_i(x_b)     *
    !                                           *
    ! So to approximate different functions     *
    ! (for the same configuration) we can       *
    ! use the SAME moment matrix, but need      *
    ! a different function vector               *
    !                                           * 
    !********************************************

    INTEGER :: k,p
    INTEGER, intent(IN) :: N, rdim, INT_param
    DOUBLE PRECISION, INTENT(IN) :: rp(rdim,N), r_G(rdim), fp(N), w(N)
    DOUBLE PRECISION, INTENT(INOUT) :: B(INT_param),P_i(INT_param)
    B = 0.
    DO p = 1, N
      IF (INT_param == NN3D_lin) THEN 
        CALL get_base3D_lin(rp(:,p),r_G,P_i)
        !CALL kernel(w,rdim,h,r_G,rp(:,p))
      ELSE IF (INT_param == NN3D_quad) THEN
        CALL get_base3D_quad(rp(:,p),r_G,P_i)
      ELSE IF (INT_param == NN3D_cub) THEN 
        CALL get_base3D_cubic(rp(:,p),r_G,P_i)
      ELSE IF (INT_param == NN3D_1d) THEN 
        P_i = (/1d0/)
      END IF 
      DO k = 1, INT_param
          B(k) = B(k) + fp(p) * P_i(k) * w(p)
      END DO
   END DO 
  END SUBROUTINE


  SUBROUTINE algebra(r,r_G,rdim,INT_param,A_inv,moment_matrix,B,beta,f_grid)

    !********************************************
    !                                           * 
    ! The solution for the coefficient          *
    ! vector is then:                           *
    !                                           *
    ! a_i= (M_ij)^(-1)*B_j                      *
    !                                           *
    ! and the function approximation is:        *
    !                                           *
    ! A(x)= sum_i a_i P_i(x)                    *
    !                                           * 
    !********************************************

    INTEGER, INTENT(IN) :: INT_param, rdim
    DOUBLE PRECISION, INTENT(IN) :: r(rdim), r_G(rdim), B(INT_param)
    DOUBLE PRECISION, INTENT(INOUT) :: moment_matrix(INT_param, INT_param)
    DOUBLE PRECISION :: pl(INT_param)
    DOUBLE PRECISION, INTENT(INOUT) :: A_inv(INT_param,INT_param), f_grid, beta(INT_param)
    f_grid = 0d0
    A_inv = 0d0
    IF (INT_param == NN3D_1d) THEN
      A_inv = 1d0 / moment_matrix
    ELSE 
      A_inv = inv(moment_matrix)
    END IF
    beta = MATMUL(A_inv,B)
    IF (INT_param == NN3D_lin) THEN 
      CALL get_base3D_lin(r,r_G,pl)
    ELSE IF (INT_param == NN3D_quad) THEN 
      CALL get_base3D_quad(r,r_G,pl)
    ELSE IF (INT_param == NN3D_cub) THEN 
      CALL get_base3D_cubic(r,r_G,pl)
    ELSE IF (INT_param == NN3D_1d) THEN 
      pl = (/1d0/)
    END IF 
    f_grid = DOT_PRODUCT(beta,pl)
    !f_grid = beta(1)
  END SUBROUTINE algebra

  SUBROUTINE file_out(rp,fp,N,rdim)
    INTEGER, INTENT(IN) :: rdim, N
    INTEGER :: i
    DOUBLE PRECISION, INTENT(IN) :: rp(rdim,N), fp(N)
    CHARACTER(LEN=20) :: file_name1 = "particles.dat"
    OPEN(unit=99, file = file_name1, status = 'unknown', action = 'write')
    DO i = 1, N
      WRITE(99,*) i, rp(1,i), rp(2,i),rp(3,i),fp(i)
    END DO
    CLOSE(99)
  END SUBROUTINE file_out

END MODULE set_up





PROGRAM main
  USE set_up
  IMPLICIT NONE
  INTEGER, PARAMETER :: ni = 18, N = ni**3, rdim = 3, ni_grid = 14, N_grid = ni_grid**3
  INTEGER, PARAMETER :: INT_param = NN3D_cub
  DOUBLE PRECISION, PARAMETER :: border = 0.9d0, border_g = 0.5d0
  DOUBLE PRECISION, PARAMETER :: dx_g = 2d0 * border_g / (REAL(ni_grid) - 1d0)
  DOUBLE PRECISION, PARAMETER :: dx = 2d0 * border / (REAL(ni) - 1d0)
  DOUBLE PRECISION, PARAMETER :: fac = 0.01d0!, dx_g = 2d0
  DOUBLE PRECISION :: r_particles(rdim,N), r2(N), r(rdim), r_G(rdim), r_grid(rdim,N_grid), fp(N), f_grid
  DOUBLE PRECISION :: moment_matrix(INT_param,INT_param), B(INT_param), beta(INT_param), A_inv(INT_param,INT_param)
  DOUBLE PRECISION :: w(rdim,N), f_analytic(1), error 
  DOUBLE PRECISION :: P_i(INT_param), h, h_g
  INTEGER, PARAMETER :: p_max = 301
  INTEGER :: i
  DOUBLE PRECISION :: T1, T2
  print*, "test"
  OPEN(unit=97, file = "error_cub.dat", status = 'unknown', action = 'write')
  !OPEN(unit=97, file = "error_lin.dat", status = 'unknown', action = 'write')
  !OPEN(unit=97, file = "error_quad.dat", status = 'unknown', action = 'write')
  OPEN(unit=99, file = "particles.dat", status = 'unknown', action = 'write')
  h = 0.02d0 !0.015d0 (100 particles (20parts), !0.03d0 (20 particles)!0.018d0 (with 100 particles in quad (30particles)) !((4d0/3d0) * pi / 30d0)**(1d0/3d0) * dx
  h_g = ((4d0/3d0) * pi / 30d0)**(1d0/3d0) * dx_g
  print*, "=================================================="
  print*, "Number of particles:  ", N
  print*, "Smoothing length:     ", h 
  print*, "Smoothing length grid:", h_g 
  print*, "Spacing dx:           ", dx
  print*, "Spacing grid:         ", dx_g
  print*, "border_part.:         ", border
  print*, "border_grid:          ", border_g
  print*, "=================================================="
  CALL generate_particles(r_particles,fp,rdim,N,border,ni,fac)
  CALL generate_grid(r_grid,rdim,N_grid,border_g,ni_grid)
  CALL cpu_time(t1)
  DO i = 1, N_grid
    r_G = r_grid(:,i)
    r = r_G
    CALL get_h(r2,r_particles, r_G, p_max,h,rdim, N)
    CALL kernel_WC6(w, h, r2, N)
    CALL get_vector_b(B,fp,P_i,N,w,r_particles,r_G,rdim,INT_param)
    CALL assemble_moment_matrix(N,INT_param,moment_matrix,P_i,w,r_particles,r_G,rdim)
    CALL algebra(r,r_G,rdim,INT_param,A_inv,moment_matrix,B,beta,f_grid)
    CALL particle_distr(r,f_analytic,rdim,1)
    error = ABS(f_grid - f_analytic(1)) / f_analytic(1) 
    WRITE(97, *) r, error 
    print*, "INDEX: ",i, "f_grid", f_grid, "f_analytic", f_analytic, "error: ", error
  END DO 
  CALL cpu_time(t2)
  print*, "======================="
  print*, "time", T2-T1
  print*, "======================="
  !CALL file_out(rp,fp,N,rdim)
  DO i = 1, N
    WRITE(99,*) i, r_particles(:,i), fp(i)
  END DO
  CLOSE(99)
  CLOSE(97)
END PROGRAM main
