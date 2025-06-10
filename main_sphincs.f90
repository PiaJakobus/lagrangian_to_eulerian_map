PROGRAM main

  !****************************************************
  ! Main program to perform SPH calculations
  ! --------------------------------------------------
  ! - Reads particle data and a list of unbound indices
  ! - Generates a grid (spherical or cartesian)
  ! - Projects SPH particle data onto the grid
  ! - Outputs grid data for further analysis
  ! Uses OpenMP for parallel processing.
  ! Pia Jakobus, 2024
  !****************************************************

  ! === Module imports ===
  USE LRE
  USE KERNEL
  USE SET_UP
  USE OMP_LIB
  USE IO
  USE INPUT_OUTPUT, ONLY : read_SPHINCS_dump
  USE SPH_variables
  USE project_2D, ONLY : integrate_kernel
  USE UNITS

  IMPLICIT NONE

  ! === Parameters and constants ===
  LOGICAL, PARAMETER           :: polar        = .true.
  INTEGER, PARAMETER           :: rdim         = 3
  INTEGER, PARAMETER           :: ni_grid      = 100
  INTEGER, PARAMETER           :: ni_grid_r    = 100
  INTEGER, PARAMETER           :: ni_grid_th   = 90
  INTEGER, PARAMETER           :: ni_grid_ph   = 180
  INTEGER, PARAMETER           :: N_grid       = ni_grid_r * ni_grid_ph * ni_grid_th
  ! === If polar = true then you need to comment in this: ===
  !INTEGER, PARAMETER           :: N_grid       = ni_grid * ni_grid * ni_grid
  INTEGER, PARAMETER           :: N_variable   = 5
  INTEGER, PARAMETER           :: INT_param    = NN3D_lin
  DOUBLE PRECISION, PARAMETER  :: border_g     = 30.d0
  DOUBLE PRECISION, PARAMETER  :: border_r     = 5000.d0

  ! === Variables ===
  INTEGER                      :: stat, i, i_g, k, l
  INTEGER                      :: zaehler, zaehler0
  INTEGER                      :: N

  DOUBLE PRECISION             :: T1, T2
  DOUBLE PRECISION             :: x, y, z, vx, vy, vz, dist, xy
  DOUBLE PRECISION             :: f_grid, f_grid2
  DOUBLE PRECISION             :: rho_to_cgs, test_var
  DOUBLE PRECISION             :: moment_matrix(INT_param, INT_param), B(INT_param)
  DOUBLE PRECISION             :: f_grid_pol(N_variable)
  DOUBLE PRECISION             :: range_r(ni_grid_r), range_th(ni_grid_th), range_ph(ni_grid_ph)

  ! === Allocatable arrays ===
  DOUBLE PRECISION, ALLOCATABLE :: r_particles(:,:), r_grid(:,:), fp(:,:), f_grid_vec(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: ye_const(:), h_subst(:)
  INTEGER, ALLOCATABLE          :: unbound_list(:)

  CHARACTER(LEN=21)            :: namefile = "data/NSNS.00740"

  ! === File I/O ===
  OPEN(unit=97, file="output/mapped_grid.dat", &
       form="unformatted", status="unknown", action="write", access="stream")

  OPEN(unit=99, file="data/unbound_list_MPA1_2x1.3.dat", status="old", action="read", iostat=stat)

  ! === Count number of unbound particles (lines in list) ===
  print *, "=================================================="
  print *, "Reading unbound particle list..."
  l = 0
  DO
    READ(99, *, iostat=stat)
    IF (stat /= 0) EXIT
    l = l + 1
  END DO

  print *, "Number of unbound particles: ", l
  print *, "=================================================="

  ! === Read SPH particle data ===
  CALL read_SPHINCS_dump(namefile)
  N = l ! This is because our number of particles is reduced to the number of unbound particles

  ! === Allocate arrays ===
  ALLOCATE(r_particles(rdim, l), r_grid(rdim, N_grid), fp(l, N_variable), &
           f_grid_vec(N_grid, N_variable), ye_const(N_grid), h_subst(l), unbound_list(l))

  REWIND(99)
  DO i = 1, l
    READ(99, *) unbound_list(i)
  END DO

  ! === Generate grid ===
  IF (.NOT. polar) THEN
    CALL generate_grid(r_grid, rdim, N_grid, border_g, ni_grid)
  ELSE
    CALL generate_polar(r_grid, rdim, N_grid, border_r, ni_grid_r, ni_grid_th, &
                        ni_grid_ph, range_r, range_th, range_ph)
  END IF

  ! === Select unbound particle data ===
  r_particles      = pos_u(:, unbound_list)
  fp(:,1)          = rho(unbound_list)
  fp(:,2)          = vel_u(1, unbound_list)
  fp(:,3)          = vel_u(2, unbound_list)
  fp(:,4)          = vel_u(3, unbound_list)
  fp(:,5)          = SQRT(SUM(r_particles**2, DIM=1)) ! Test function for error

  ! === Unit conversion ===
  CALL set_units('NSM')
  rho_to_cgs = rho_cu_fm_3 * baryon_fm3_to_cgs
  ye_const   = 0.1D0

  ! === Start parallel  ===
  CALL OMP_SET_NUM_THREADS(1)
  T1 = OMP_GET_WTIME()

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP& PRIVATE(i_g, B, moment_matrix, f_grid, vx, vy, vz, xy, f_grid_pol, dist, x, y, z, zaehler, zaehler0) &
!$OMP& SHARED(r_grid, fp, r_particles, N, f_grid_vec, h, pmass, rho, h_subst, unbound_list)
  DO i_g = 1, N_grid
    CALL get_grid_value(i_g, INT_param, rdim, N_variable, zaehler, zaehler0, &
                        fp, N, r_grid, N_grid, f_grid_vec, rho(unbound_list), &
                        pmass(unbound_list), h(unbound_list), moment_matrix, r_particles, &
                        unbound_list, f_grid, B, polar)

    IF (polar) THEN
      CALL to_polar(r_grid(:, i_g), rdim, N_variable, c_light, dist, f_grid_vec(i_g,:), f_grid_pol)
      f_grid_vec(i_g,:) = f_grid_pol
    END IF
    ! === Comment in following line only for testing! ===
    print *, dist*1e-5, "rho", f_grid_vec(i_g,1), "error", f_grid_vec(i_g,5), "zaehler", zaehler
  END DO
!$OMP END PARALLEL DO

  T2 = OMP_GET_WTIME()
  print *, "======================="
  print *, "Total runtime (s):", T2 - T1
  print *, "======================="

  ! === Write output data ===
  WRITE(97) t * 1e-3
  WRITE(97) ni_grid_r
  WRITE(97) ni_grid_th
  WRITE(97) ni_grid_ph
  WRITE(97) (range_r(i) * 1.0D5, i = 1, ni_grid_r)
  WRITE(97) (range_th(i), i = 1, ni_grid_th)
  WRITE(97) (range_ph(i), i = 1, ni_grid_ph)
  WRITE(97) f_grid_vec(:,1) * rho_to_cgs  ! rho
  WRITE(97) ye_const                      ! Ye
  WRITE(97) f_grid_vec(:,2)              ! vr
  WRITE(97) f_grid_vec(:,3)              ! vth
  WRITE(97) f_grid_vec(:,4)              ! vph
  WRITE(97) f_grid_vec(:,5)              ! test function

  ! === Cleanup ===
  CLOSE(99)
  CLOSE(97)
  DEALLOCATE(r_particles, r_grid, fp, f_grid_vec, ye_const, h_subst)
  CALL DEALLOCATE_SPH_memory

END PROGRAM main

