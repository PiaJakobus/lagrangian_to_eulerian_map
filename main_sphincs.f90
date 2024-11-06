PROGRAM main

  !**********************************************
  !                                             *
  ! Main program to perform SPH calculations    * 
  ! Reads particle data, generates a grid, and  *
  ! computes properties using different kernels *
  ! and algebraic methods.                      *
  !                                             *
  ! Uses OpenMP for parallel processing to      *
  ! speed up grid-based calculations.           *
  !  Pia Jakobus 2024                           *
  !**********************************************

  USE LRE
  USE KERNEL
  USE SET_UP
  USE OMP_LIB
  USE IO
  USE INPUT_OUTPUT, ONLY : read_SPHINCS_dump
  USE SPH_variables
  USE project_2D, ONLY : integrate_kernel

  IMPLICIT NONE

  INTEGER, PARAMETER :: rdim = 3                    ! Number of dimensions (3D)
  INTEGER, PARAMETER :: ni_grid = 100                ! Number of grid points along each axis
  INTEGER, PARAMETER :: N_grid = ni_grid**2         ! Total number of grid points
  INTEGER :: N                                      ! Number of particles
  INTEGER, PARAMETER :: INT_param = NN3D_1d         ! Number of interpolation parameters
  DOUBLE PRECISION, PARAMETER :: border_g = 300d0   ! Grid boundary size
  DOUBLE PRECISION, PARAMETER :: dx_g = 2d0 * border_g / (REAL(ni_grid) - 1d0)  ! Grid spacing
  DOUBLE PRECISION :: r(rdim), r_G(rdim), f_grid, f_grid2  ! Position and grid values
  DOUBLE PRECISION :: moment_matrix(INT_param,INT_param), B(INT_param)  ! Moment matrix and vector B
  DOUBLE PRECISION :: h_smooth, h_g, f_grid_vec(N_grid)  ! Smoothing lengths and grid results vector
  INTEGER :: i_g, i                                 ! Loop indices
  DOUBLE PRECISION :: T1, T2                        ! Timing variables for performance analysis
  DOUBLE PRECISION, ALLOCATABLE :: r_particles(:,:), r_grid(:,:), fp(:), r2(:)  ! Allocatable arrays
  CHARACTER(LEN=21) :: namefile = "data/BHWD.00268" ! Filename for input data

  ! Open files for output
  OPEN(unit=97, file = "mapped_T.dat", status = 'unknown', action = 'write')
  OPEN(unit=99, file = "particles_T.dat", status = 'unknown', action = 'write')

  ! Output simulation parameters to the console
  print*, "=================================================="
  print*, "Number of particles:  ", N
  print*, "Smoothing length:     ", h_smooth
  print*, "Smoothing length grid:", h_g 
  print*, "Spacing grid:         ", dx_g
  print*, "border_grid:          ", border_g
  print*, "=================================================="

  ! Read particle data from file
  CALL read_SPHINCS_dump(namefile)
  N = npart  ! Number of particles read from file

  ! Allocate arrays based on number of particles and grid size
  ALLOCATE(r_particles(rdim, N), r_grid(rdim, N_grid), fp(N), r2(N))

  ! Generate grid for SPH calculations
  CALL generate_grid(r_grid, rdim, N_grid, border_g, ni_grid)
  r_particles = pos_u   ! Set particle positions
  fp = rho              ! Set density values for particles

  ! Set the number of OpenMP threads for parallel processing
  CALL OMP_SET_NUM_THREADS(24)
  T1 = OMP_GET_WTIME()  ! Start timing for performance measurement

  ! Parallel loop over all grid points
!!$OMP PARALLEL DO DEFAULT(NONE) & 
!!$OMP& PRIVATE(i_g,r_G,r,r2,B,moment_matrix,f_grid,f_grid2) &
!!$OMP& SHARED(r_grid,fp,r_particles,N,f_grid_vec,h, pmass,rho)
!$OMP PARALLEL DO DEFAULT(NONE) & 
!$OMP& PRIVATE(f_grid,r,r_G,i_g) & 
!$OMP& SHARED(r_particles,r_grid,fp,N,h,f_grid_vec)
  DO i_g = 1, N_grid
    r_G = r_grid(:, i_g)  ! Get the coordinates of the current grid point
    r = r_G               ! Set r to the current grid point position
    CALL integrate_kernel(r_G, r_particles, h, rdim, N, f_grid,fp)
    print*, f_grid, r
    !print*, f_grid
    !stop
    !CALL get_h(r2,r_particles, r_G, p_max,h_smooth,rdim, N)
    ! Calculate SPH values for the grid point using standard SPH method
    !CALL standard_SPH(fp, N, r_particles, r_G, rdim, h, rho, pmass, f_grid2)

    ! Calculate vector B for the given grid point and assemble the moment matrix
    !CALL get_vector_b(B, fp, N, r_particles, r_G, rdim, INT_param, h)
    !CALL assemble_moment_matrix(N, INT_param, moment_matrix, r_particles, r_G, rdim, h)

    ! Perform algebraic operations to compute grid values
    !CALL algebra(r, r_G, rdim, INT_param, moment_matrix, B, f_grid)
    f_grid_vec(i_g) = f_grid  ! Store the computed value for the grid point

    ! Output information for debugging
    !print*, "thread: ", omp_get_thread_num(), "INDEX: ", i_g, "f_grid", f_grid, "f_grid2", f_grid2, "ratio", f_grid / f_grid2
  END DO 
!$OMP END PARALLEL DO 

  T2 = OMP_GET_WTIME()  ! End timing for performance measurement
  print*, "======================="
  print*, "time", T2-T1           ! Output the time taken for the calculations
  print*, "======================="

  ! Write the computed grid values to an output file
  DO i_g = 1, N_grid
    WRITE(97, *) i_g, r_grid(1:2, i_g), f_grid_vec(i_g)
  END DO

  ! Write particle data to an output file
  DO i = 1, N
    WRITE(99, *) i, r_particles(:, i), fp(i)
  END DO

  ! Close output files
  CLOSE(99)
  CLOSE(97)

  ! Deallocate arrays to free memory
  DEALLOCATE(r_particles, r_grid, fp, r2)
  CALL DEALLOCATE_SPH_memory

END PROGRAM main

