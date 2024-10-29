PROGRAM main
  USE LRE
  USE KERNEL
  USE SET_UP
  USE OMP_LIB
  USE IO
  USE INPUT_OUTPUT
  USE SPH_variables

  IMPLICIT NONE
  INTEGER, PARAMETER :: ni = 18, N = ni**3, rdim = 3, ni_grid = 30, N_grid = ni_grid**3
  INTEGER, PARAMETER :: INT_param = NN3D_lin
  DOUBLE PRECISION, PARAMETER :: border = 0.9d0, border_g = 0.5d0
  DOUBLE PRECISION, PARAMETER :: dx_g = 2d0 * border_g / (REAL(ni_grid) - 1d0)
  DOUBLE PRECISION, PARAMETER :: dx = 2d0 * border / (REAL(ni) - 1d0)
  DOUBLE PRECISION, PARAMETER :: fac = 0.01d0!, dx_g = 2d0
  DOUBLE PRECISION :: r(rdim), r_G(rdim), f_grid
  DOUBLE PRECISION :: moment_matrix(INT_param,INT_param), B(INT_param), beta(INT_param), A_inv(INT_param,INT_param)
  DOUBLE PRECISION :: f_analytic(1), error 
  DOUBLE PRECISION :: P_i(INT_param), h_smooth, h_g
  INTEGER, PARAMETER :: p_max = 31
  INTEGER :: i_g
  DOUBLE PRECISION :: T1, T2
  DOUBLE PRECISION, ALLOCATABLE:: w(:), r_particles(:,:), r_grid(:,:), fp(:), r2(:)
  CHARACTER(LEN= 21) :: namefile = "data/BHWD.00268"

  ALLOCATE(r_particles(rdim,N),r_grid(rdim,N_grid),fp(N),r2(N),w(N))

  OPEN(unit=97, file = "error_quad.dat", status = 'unknown', action = 'write')
  !OPEN(unit=97, file = "error_lin.dat", status = 'unknown', action = 'write')
  !OPEN(unit=97, file = "error_quad.dat", status = 'unknown', action = 'write')
  OPEN(unit=99, file = "particles.dat", status = 'unknown', action = 'write')
  print*, "=================================================="
  print*, "Number of particles:  ", N
  print*, "Smoothing length:     ", h_smooth
  print*, "Smoothing length grid:", h_g 
  print*, "Spacing dx:           ", dx
  print*, "Spacing grid:         ", dx_g
  print*, "border_part.:         ", border
  print*, "border_grid:          ", border_g
  print*, "=================================================="
  CALL generate_particles(r_particles,fp,rdim,N,border,ni,fac)
  CALL generate_grid(r_grid,rdim,N_grid,border_g,ni_grid)
  T1 = OMP_GET_WTIME()
  CALL OMP_SET_NUM_THREADS(6)
  CALL read_SPHINCS_dump(namefile)

  stop
!$OMP PARALLEL DO DEFAULT(NONE) & 
!$OMP& PRIVATE(i,r2,r_G,r,h_smooth,w,B,P_i,A_inv, moment_matrix,beta,f_grid,f_analytic,error) &
!$OMP& SHARED(r_grid, fp, r_particles)
!OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(r_grid,fp,r_particles)
  DO i_g = 1, N_grid
    r_G = r_grid(:,i_g)
    r = r_G
    CALL get_h(r2,r_particles, r_G, p_max,h_smooth,rdim, N)
    CALL kernel_WC6(w, h_smooth, r2, N)
    CALL get_vector_b(B,fp,P_i,N,w,r_particles,r_G,rdim,INT_param)
    CALL assemble_moment_matrix(N,INT_param,moment_matrix,P_i,w,r_particles,r_G,rdim)
    CALL algebra(r,r_G,rdim,INT_param,A_inv,moment_matrix,B,beta,f_grid)
    CALL particle_distr(r,f_analytic,rdim,1)
    error = ABS(f_grid - f_analytic(1)) / f_analytic(1) 
    WRITE(97, *) r, error 
    print*, "thread: ", omp_get_thread_num(), "INDEX: ",i_g, "f_grid", f_grid, "f_analytic", f_analytic, "error: ", error
  END DO 
!$OMP END PARALLEL DO 
  T2 = OMP_GET_WTIME()
  print*, "======================="
  print*, "time", T2-T1
  print*, "======================="
  !CALL file_out(r_particles,fp,N,rdim)
  DO i = 1, N
    WRITE(99,*) i, r_particles(:,i), fp(i)
  END DO
  CLOSE(99)
  CLOSE(97)
  DEALLOCATE(r_particles,r_grid,fp,r2,w)
  CALL DEALLOCATE_SPH_memory
END PROGRAM main
