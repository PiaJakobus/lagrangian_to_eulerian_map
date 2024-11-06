PROGRAM main
  USE LRE
  USE KERNEL
  USE SET_UP
  USE OMP_LIB
  USE IO
  USE INPUT_OUTPUT, ONLY : read_SPHINCS_dump
  USE SPH_variables

  IMPLICIT NONE
  INTEGER, PARAMETER :: rdim = 3, ni_grid = 50, N_grid = ni_grid**3
  INTEGER :: N
  INTEGER, PARAMETER :: INT_param = NN3D_1d
  INTEGER, PARAMETER :: p_max = 31
  DOUBLE PRECISION, PARAMETER :: border_g = 300d0
  DOUBLE PRECISION, PARAMETER :: dx_g = 2d0 * border_g / (REAL(ni_grid) - 1d0)
  DOUBLE PRECISION :: r(rdim), r_G(rdim), f_grid
  DOUBLE PRECISION :: moment_matrix(INT_param,INT_param), B(INT_param), beta(INT_param), A_inv(INT_param,INT_param)
  DOUBLE PRECISION :: P_i(INT_param), h_smooth, h_g, f_grid_vec(N_grid)
  INTEGER :: i_g, i_p, i 
  DOUBLE PRECISION :: T1, T2
  DOUBLE PRECISION, ALLOCATABLE:: w(:), r_particles(:,:), r_grid(:,:), fp(:), r2(:), index_list(:)
  CHARACTER(LEN= 21) :: namefile = "data/BHWD.00268"


  OPEN(unit=97, file = "mapped_T5.dat", status = 'unknown', action = 'write')
  OPEN(unit=99, file = "particles_T5.dat", status = 'unknown', action = 'write')
  print*, "=================================================="
  print*, "Number of particles:  ", N
  print*, "Smoothing length:     ", h_smooth
  print*, "Smoothing length grid:", h_g 
  print*, "Spacing grid:         ", dx_g
  print*, "border_grid:          ", border_g
  print*, "=================================================="
  !CALL generate_particles(r_particles,fp,rdim,N,border,ni,fac)
  CALL read_SPHINCS_dump(namefile)
  N = npart 
  ALLOCATE(r_particles(rdim,N),r_grid(rdim,N_grid),fp(N),r2(N),w(N),index_list(N))
  CALL generate_grid(r_grid,rdim,N_grid,border_g,ni_grid)
  r_particles = pos_u
  fp = rho
  !r_grid = pos_u
  CALL OMP_SET_NUM_THREADS(24)
  T1 = OMP_GET_WTIME()
  !DO i_p = 1, N_grid
  !print*, r_grid(:,i_p)
  !END DO 
  !stop 
 
!$OMP PARALLEL DO DEFAULT(NONE) & 
!$OMP& PRIVATE(i_g,i_p,r2,r_G,r,h_smooth,w,B,P_i,A_inv, moment_matrix,beta,f_grid) &
!$OMP& SHARED(r_grid, fp,r_particles,N,f_grid_vec,p_max,index_list)
!OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(r_grid,fp,r_particles)
  DO i_g = 1, N_grid
    r_G = r_grid(:,i_g)
    r = r_G
    CALL get_h(r2,r_particles, r_G, p_max,h_smooth,rdim, N,index_list)
    !DO i_p = 1, N 
    !  r2(i_p) = (r_particles(1,i_p)-r(1))**2d0 + &
    !          & (r_particles(2,i_p)-r(2))**2d0 + & 
    !          & (r_particles(3,i_p)-r(3))**2d0
    !END DO 
    CALL kernel_WC6(w, h_smooth, r2, N)
    CALL get_vector_b(B,fp,P_i,N,w,r_particles,r_G,rdim,INT_param)
    CALL assemble_moment_matrix(N,INT_param,moment_matrix,P_i,w,r_particles,r_G,rdim)
    CALL algebra(r,r_G,rdim,INT_param,A_inv,moment_matrix,B,beta,f_grid)
    f_grid_vec(i_g) = f_grid
    !error = ABS(f_grid - fp(i_g)) / fp(i_g) 
    !WRITE(97, *) i_g, r, f_grid
    print*, "thread: ", omp_get_thread_num(), "INDEX: ",i_g, "f_grid", f_grid 
  END DO 
!$OMP END PARALLEL DO 
  T2 = OMP_GET_WTIME()
  print*, "======================="
  print*, "time", T2-T1
  print*, "======================="
  !CALL file_out(r_particles,fp,N,rdim)
  DO i_g = 1, N_grid
    WRITE(97,*) i_g, r_grid(:,i_g), f_grid_vec(i_g)
  END DO
  DO i = 1, N
    WRITE(99,*) i, r_particles(:,i), fp(i)
  END DO
  CLOSE(99)
  CLOSE(97)
  DEALLOCATE(r_particles,r_grid,fp,r2,w)
  CALL DEALLOCATE_SPH_memory
END PROGRAM main
