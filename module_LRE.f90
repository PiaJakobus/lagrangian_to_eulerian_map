MODULE LRE

  !=============================================================
  ! Module: LRE (Least Regression Estimate)
  ! ------------------------------------------------------------
  ! Contains routines for assembling moment matrices, computing
  ! function vectors, evaluating function values via SPH/LRE
  ! methods, and transforming between coordinate systems.
  ! SPH = 0th order LRE
  !
  ! Author: SKR, modified by Pia Jakobus 2024
  !=============================================================

  USE KERNEL, ONLY: Wendland_C6_3D, kernel_WC6
  USE NR
  IMPLICIT NONE

  ! === Constants for basis function dimensionality ===
  INTEGER, PARAMETER :: NN3D_1d   = 1
  INTEGER, PARAMETER :: NN3D_lin  = 4
  INTEGER, PARAMETER :: NN3D_quad = 10
  INTEGER, PARAMETER :: NN3D_cub  = 20

  CONTAINS

  SUBROUTINE assemble_moment_matrix(N, INT_param, moment_matrix, rp, r_G, rdim, h, rho, m, polar, zaehler)
    !=============================================================
    ! Assemble the moment matrix M_ij for the LRE method
    ! method using a given basis and weighting function.
    !
    ! M_ij = sum_b W(q) * P_i(x_b) * P_j(x_b)
    !
    ! Inputs:
    !   N            - Number of particles
    !   INT_param    - Number of basis functions (depends on order)
    !   rp           - Particle positions (rdim x N)
    !   r_G          - Grid point position
    !   rdim         - Number of dimensions (usually 3)
    !   h            - Smoothing lengths
    !   rho          - Densities
    !   m            - Masses
    !   polar        - If true, uses spherical coordinate distances
    !
    ! Outputs:
    !   moment_matrix - The resulting moment matrix
    !   zaehler       - Counter for contributing particles (q <= 1)
    !=============================================================

    INTEGER, INTENT(IN)           :: INT_param, N, rdim
    INTEGER, INTENT(INOUT)        :: zaehler
    DOUBLE PRECISION, INTENT(IN)  :: rp(rdim, N), r_G(rdim), h(N), rho(N), m(N)
    DOUBLE PRECISION, INTENT(OUT) :: moment_matrix(INT_param, INT_param)
    LOGICAL, INTENT(IN)           :: polar

    DOUBLE PRECISION              :: P_i(INT_param)
    DOUBLE PRECISION              :: w, w_tilde, q
    DOUBLE PRECISION              :: pos_ab(rdim)
    INTEGER                       :: i, p, k

    moment_matrix = 0d0
    zaehler = 0

    DO p = 1, N
      IF (.NOT. polar) THEN
        pos_ab = ABS(rp(:,p) - r_G)
        q = 0.5D0 * SQRT(SUM(pos_ab**2)) / h(p)
      ELSE
        CALL distance_polar(rp(:,p), r_G, q)
        q = 0.5D0 * q / h(p)
      END IF

      IF (q <= 1.D0) THEN
        CALL Wendland_C6_3D(q, w, zaehler)

        SELECT CASE (INT_param)
        CASE (NN3D_quad)
          CALL get_base3D_quad(rp(:,p), r_G, P_i, polar)
        CASE (NN3D_lin)
          CALL get_base3D_lin(rp(:,p), r_G, P_i, polar)
        CASE (NN3D_cub)
          CALL get_base3D_cubic(rp(:,p), r_G, P_i, polar)
        CASE (NN3D_1d)
          P_i = (/1d0/)
        END SELECT

        w_tilde = w * m(p) / (rho(p) * h(p)**3.D0)

        DO k = 1, INT_param
          DO i = 1, INT_param
            moment_matrix(i,k) = moment_matrix(i,k) + P_i(i) * P_i(k) * w_tilde
          END DO
        END DO
      END IF
    END DO

  END SUBROUTINE assemble_moment_matrix

  SUBROUTINE get_vector_b(B, fp, N, rp, r_G, rdim, INT_param, h, rho, m, polar)
    !=============================================================
    ! Assemble the function vector B_i = sum_b A_b * W(q) * P_i(x_b)
    !
    ! Inputs:
    !   fp         - Function values at particle positions
    !   rp         - Particle positions
    !   r_G        - Grid point position
    !   h, rho, m  - Smoothing lengths, densities, masses
    !   polar      - Whether to use spherical coordinates
    !
    ! Output:
    !   B          - Resulting function vector 
    !=============================================================

    INTEGER, INTENT(IN)           :: N, rdim, INT_param
    DOUBLE PRECISION, INTENT(IN)  :: rp(rdim, N), r_G(rdim), fp(N), h(N), rho(N), m(N)
    DOUBLE PRECISION, INTENT(OUT) :: B(INT_param)
    LOGICAL, INTENT(IN)           :: polar

    DOUBLE PRECISION              :: P_i(INT_param)
    DOUBLE PRECISION              :: w, w_tilde, q
    DOUBLE PRECISION              :: pos_ab(rdim)
    INTEGER                       :: k, p, zaehler

    B = 0.D0
    zaehler = 0

    DO p = 1, N
      IF (.NOT. polar) THEN
        pos_ab = rp(:,p) - r_G
        q = 0.5D0 * SQRT(SUM(pos_ab**2)) / h(p)
      ELSE
        CALL distance_polar(rp(:,p), r_G, q)
        q = 0.5D0 * q / h(p)
      END IF

      IF (q <= 1.D0) THEN
        CALL Wendland_C6_3D(q, w, zaehler)

        SELECT CASE (INT_param)
        CASE (NN3D_lin)
          CALL get_base3D_lin(rp(:,p), r_G, P_i, polar)
        CASE (NN3D_quad)
          CALL get_base3D_quad(rp(:,p), r_G, P_i, polar)
        CASE (NN3D_cub)
          CALL get_base3D_cubic(rp(:,p), r_G, P_i, polar)
        CASE (NN3D_1d)
          P_i = (/1d0/)
        END SELECT

        w_tilde = w * m(p) / (rho(p) * h(p)**3.D0)

        DO k = 1, INT_param
          B(k) = B(k) + fp(p) * P_i(k) * w_tilde
        END DO
      END IF
    END DO

  END SUBROUTINE get_vector_b

  PURE SUBROUTINE distance_polar(xp, rg, dist)
    !=============================================================
    ! Compute Euclidean distance between a particle (xp) and a point
    ! given in spherical coordinates (rg), assuming rg = (r, theta, phi).
    !
    ! Inputs:
    !   xp   - Cartesian position of particle
    !   rg   - Spherical coordinate position (r, theta, phi)
    !
    ! Output:
    !   dist - Distance between xp and the Cartesian equivalent of rg
    !=============================================================

    DOUBLE PRECISION, INTENT(IN)    :: xp(3), rg(3)
    DOUBLE PRECISION, INTENT(INOUT) :: dist
    DOUBLE PRECISION                :: x, y, z

    x = rg(1) * SIN(rg(2)) * COS(rg(3))
    y = rg(1) * SIN(rg(2)) * SIN(rg(3))
    z = rg(1) * COS(rg(2))

    dist = SQRT((xp(1) - x)**2 + (xp(2) - y)**2 + (xp(3) - z)**2)

  END SUBROUTINE distance_polar

  SUBROUTINE standard_SPH(fp, N, rp, r_G, rdim, h, rho, m, f_grid, polar, zaehler)
    !=============================================================
    ! Compute SPH estimate of a function at a grid point r_G.
    !
    ! Inputs:
    !   fp      - Function values at particle positions
    !   rp      - Particle positions (rdim x N)
    !   r_G     - Grid point position
    !   rdim    - Number of spatial dimensions (typically 3)
    !   h       - Smoothing lengths
    !   rho     - Densities
    !   m       - Particle masses
    !   polar   - If true, use polar coordinate distance
    !
    ! Outputs:
    !   f_grid  - Resulting interpolated function value
    !   zaehler - Counter for included neighbors (q <= 1)
    !=============================================================

    INTEGER, INTENT(IN)           :: N, rdim
    LOGICAL, INTENT(IN)           :: polar
    DOUBLE PRECISION, INTENT(IN)  :: fp(N), rho(N), m(N), rp(rdim,N), r_G(rdim), h(N)
    DOUBLE PRECISION, INTENT(OUT) :: f_grid
    INTEGER, INTENT(INOUT)        :: zaehler

    DOUBLE PRECISION              :: w, q, pos_ab(rdim), norm, omega
    INTEGER                       :: p

    norm = 1e-40
    f_grid = 0.d0
    zaehler = 0

    DO p = 1, N
      IF (.NOT. polar) THEN
        pos_ab = rp(:,p) - r_G
        q = 0.5D0 * SQRT(SUM(pos_ab**2)) / h(p)
      ELSE
        CALL distance_polar(rp(:,p), r_G, q)
        q = 0.5D0 * q / h(p)
      END IF

      omega = 0.D0

      IF (q <= 1.D0) THEN
        CALL Wendland_C6_3D(q, w, zaehler)
        omega = (m(p) / rho(p)) * (w / h(p)**3.D0)
        norm = norm + omega
        f_grid = f_grid + omega * fp(p)
      END IF
    END DO

    f_grid = f_grid / norm

  END SUBROUTINE standard_SPH

  SUBROUTINE algebra(INT_param, moment_matrix, B, f_grid)
    !=============================================================
    ! Solve the system M * beta = B for the coefficients beta.
    !
    ! For INT_param = 1 (constant basis), use direct inversion.
    ! Otherwise, use general matrix inversion.
    !
    ! Inputs:
    !   INT_param     - Number of basis functions
    !   moment_matrix - Moment matrix M
    !   B             - Function vector
    !
    ! Output:
    !   f_grid        - Approximated function value at grid point
    !=============================================================

    INTEGER, INTENT(IN)           :: INT_param
    DOUBLE PRECISION, INTENT(IN)  :: B(INT_param), moment_matrix(INT_param, INT_param)
    DOUBLE PRECISION, INTENT(OUT) :: f_grid

    DOUBLE PRECISION              :: A_inv(INT_param, INT_param), beta(INT_param)

    f_grid = 0.D0
    A_inv  = 0.D0

    IF (INT_param == NN3D_1d) THEN
      IF (ABS(moment_matrix(1,1)) > 0.D0) THEN
        A_inv(1,1) = 1.D0 / moment_matrix(1,1)
      ELSE
        A_inv(1,1) = 0.D0
      END IF
    ELSE
      A_inv = inv(moment_matrix)
    END IF

    beta   = MATMUL(A_inv, B)
    f_grid = beta(1)

  END SUBROUTINE algebra

  SUBROUTINE get_grid_value(i_g, INT_param, rdim, N_variable, zaehler, zaehler0, &
                           fp, N, r_grid, N_grid, f_grid_vec, rho_red, pmass_red, h_red, &
                           moment_matrix, r_particles, unbound_list, f_grid, B, polar)
    !=============================================================
    ! Compute grid values for each variable at a specific grid point.
    !
    ! This routine uses either LRE or SPH depending on particle density
    ! to project particle values onto the grid.
    !
    ! Inputs:
    !   i_g          - Index of the grid point
    !   INT_param    - Number of basis functions
    !   rdim         - Number of dimensions
    !   N_variable   - Number of physical variables
    !   fp           - Particle function values (N x N_variable)
    !   N            - Number of particles
    !   r_grid       - Grid positions (rdim x N_grid)
    !   N_grid       - Total number of grid points
    !   rho_red      - Particle densities
    !   pmass_red    - Particle masses
    !   h_red        - Particle smoothing lengths
    !   r_particles  - Particle positions (rdim x N)
    !   unbound_list - Indices of unbound particles
    !   polar        - Use polar coordinates if .true.
    !
    ! Outputs:
    !   f_grid_vec   - Grid values for each variable (N_grid x N_variable)
    !   f_grid       - Temporary storage for individual variable value
    !   moment_matrix, B - Working arrays for LRE
    !   zaehler, zaehler0 - Diagnostic counters
    !=============================================================

    INTEGER, INTENT(IN)             :: i_g, N, N_variable, N_grid, rdim, INT_param
    INTEGER, INTENT(IN)             :: unbound_list(N)
    INTEGER, INTENT(INOUT)          :: zaehler, zaehler0
    DOUBLE PRECISION, INTENT(IN)    :: fp(N, N_variable), r_grid(rdim, N_grid)
    DOUBLE PRECISION, INTENT(IN)    :: rho_red(N), pmass_red(N), h_red(N)
    DOUBLE PRECISION, INTENT(IN)    :: r_particles(rdim, N)
    DOUBLE PRECISION, INTENT(INOUT) :: moment_matrix(INT_param, INT_param), B(INT_param)
    DOUBLE PRECISION, INTENT(OUT)   :: f_grid_vec(N_grid, N_variable)
    DOUBLE PRECISION, INTENT(OUT)   :: f_grid
    LOGICAL, INTENT(IN)             :: polar

    DOUBLE PRECISION                :: h_subst(N)
    INTEGER                         :: k

    CALL assemble_moment_matrix(N, INT_param, moment_matrix, r_particles, r_grid(:, i_g), rdim, &
                                h_red, rho_red, pmass_red, polar, zaehler)

    DO k = 1, N_variable
      IF (zaehler <= 10) THEN
        h_subst = 2.5D0 * h_red
        CALL standard_SPH(fp(:, k), N, r_particles, r_grid(:, i_g), rdim, h_subst, &
                          rho_red, pmass_red, f_grid, polar, zaehler0)
      ELSE
        CALL get_vector_b(B, fp(:, k), N, r_particles, r_grid(:, i_g), rdim, INT_param, &
                          h_red, rho_red, pmass_red, polar)
        CALL algebra(INT_param, moment_matrix, B, f_grid)
      END IF

      f_grid_vec(i_g, k) = f_grid
    END DO

  END SUBROUTINE get_grid_value

  SUBROUTINE to_polar(r, rdim, N_variable, c_light, dist, f_grid_vec_i, f_grid_pol)
    !=============================================================
    ! Transform physical vector components to spherical (polar) coordinates.
    !
    ! Inputs:
    !   r             - Grid point in spherical coordinates (r, theta, phi)
    !   rdim          - Number of dimensions (should be 3)
    !   N_variable    - Number of physical variables
    !   c_light       - Speed of light 
    !   f_grid_vec_i  - Original variable values at grid point (cartesian)
    !
    ! Outputs:
    !   dist          - Radial distance of the grid point (for diagnostics)
    !   f_grid_pol    - Transformed variables in spherical coordinates
    !=============================================================

    INTEGER, INTENT(IN)           :: rdim, N_variable
    DOUBLE PRECISION, INTENT(IN)  :: r(rdim), c_light
    DOUBLE PRECISION, INTENT(IN)  :: f_grid_vec_i(N_variable)
    DOUBLE PRECISION, INTENT(OUT) :: dist, f_grid_pol(N_variable)

    DOUBLE PRECISION              :: x, y, z, vx, vy, vz, xy

    ! Velocity components scaled
    vx = f_grid_vec_i(2) * c_light
    vy = f_grid_vec_i(3) * c_light
    vz = f_grid_vec_i(4) * c_light

    ! Convert spherical position to cartesian coordinates
    x = r(1) * 1.0D5 * SIN(r(2)) * COS(r(3))
    y = r(1) * 1.0D5 * SIN(r(2)) * SIN(r(3))
    z = r(1) * 1.0D5 * COS(r(2))
    xy = x**2 + y**2
    dist = SQRT(xy + z**2)

    ! Transform components
    f_grid_pol(1) = f_grid_vec_i(1)                                 ! density
    f_grid_pol(2) = (x*vx + y*vy + z*vz) / dist                     ! radial velocity
    f_grid_pol(4) = (z*f_grid_pol(2) - dist*vz) / (dist**2 * SIN(r(2))) ! phi velocity
    f_grid_pol(3) = (x*vy - y*vx) / xy                              ! theta velocity
    f_grid_pol(5) = ABS(1.0D0 - (f_grid_vec_i(5) / (dist * 1.0D-5))) ! test function error

  END SUBROUTINE to_polar

  PURE SUBROUTINE get_base3D_1d(bases_3D_lin)
    !=============================================================
    ! Return constant basis for 1D approximation (NN3D_1d = 1).
    !=============================================================
    DOUBLE PRECISION, INTENT(OUT) :: bases_3D_lin(NN3D_1d)
    bases_3D_lin(1) = 1.0D0
  END SUBROUTINE get_base3D_1d

  PURE SUBROUTINE get_base3D_lin(pos, pos0, bases_3D_lin, polar)
    !=============================================================
    ! Return values of linear 3D basis functions.
    !
    ! Inputs:
    !   pos, pos0     - Particle and grid point positions (cartesian or spherical)
    !   polar         - If true, convert spherical to cartesian before evaluating
    !
    ! Output:
    !   bases_3D_lin  - Linear basis function values
    !=============================================================

    DOUBLE PRECISION, INTENT(IN)  :: pos(3), pos0(3)
    LOGICAL, INTENT(IN)           :: polar
    DOUBLE PRECISION, INTENT(OUT) :: bases_3D_lin(NN3D_lin)

    DOUBLE PRECISION              :: x(3)

    bases_3D_lin(1) = 1.0D0
    x = pos0

    IF (polar) THEN
      x(1) = pos0(1) * SIN(pos0(2)) * COS(pos0(3))
      x(2) = pos0(1) * SIN(pos0(2)) * SIN(pos0(3))
      x(3) = pos0(1) * COS(pos0(2))
    END IF

    bases_3D_lin(2) = pos(1) - x(1)
    bases_3D_lin(3) = pos(2) - x(2)
    bases_3D_lin(4) = pos(3) - x(3)

  END SUBROUTINE get_base3D_lin

  SUBROUTINE get_base3D_quad(pos, pos0, bases_3D_quad, polar)
    !=============================================================
    ! Return values for 3D quadratic basis functions.
    !
    ! Inputs:
    !   pos, pos0        - Position of particle and grid point (3D)
    !   polar            - If true, transform spherical to cartesian
    !
    ! Output:
    !   bases_3D_quad    - Evaluated quadratic basis functions
    !=============================================================

    DOUBLE PRECISION, INTENT(IN)  :: pos(3), pos0(3)
    LOGICAL, INTENT(IN)           :: polar
    DOUBLE PRECISION, INTENT(OUT) :: bases_3D_quad(NN3D_quad)

    DOUBLE PRECISION              :: dx, dy, dz
    DOUBLE PRECISION              :: x(3)

    x = pos0
    IF (polar) THEN
      x(1) = pos0(1) * SIN(pos0(2)) * COS(pos0(3))
      x(2) = pos0(1) * SIN(pos0(2)) * SIN(pos0(3))
      x(3) = pos0(1) * COS(pos0(2))
    END IF

    dx = pos(1) - x(1)
    dy = pos(2) - x(2)
    dz = pos(3) - x(3)

    bases_3D_quad(1)  = 1.0D0
    bases_3D_quad(2)  = dx
    bases_3D_quad(3)  = dy
    bases_3D_quad(4)  = dz
    bases_3D_quad(5)  = dx * dx
    bases_3D_quad(6)  = dx * dy
    bases_3D_quad(7)  = dx * dz
    bases_3D_quad(8)  = dy * dy
    bases_3D_quad(9)  = dy * dz
    bases_3D_quad(10) = dz * dz

  END SUBROUTINE get_base3D_quad

  SUBROUTINE get_base3D_cubic(pos, pos0, bases_3D_cubic, polar)
    !=============================================================
    ! Return values for 3D cubic basis functions.
    !
    ! Inputs:
    !   pos, pos0        - Position of particle and grid point (3D)
    !   polar            - If true, transform spherical to cartesian
    !
    ! Output:
    !   bases_3D_cubic   - Evaluated cubic basis functions
    !=============================================================

    DOUBLE PRECISION, INTENT(IN)  :: pos(3), pos0(3)
    LOGICAL, INTENT(IN)           :: polar
    DOUBLE PRECISION, INTENT(OUT) :: bases_3D_cubic(NN3D_cub)

    DOUBLE PRECISION              :: dx, dy, dz, dxx, dyy, dzz, dxy, dxz, dyz
    DOUBLE PRECISION              :: x(3)

    x = pos0
    IF (polar) THEN
      x(1) = pos0(1) * SIN(pos0(2)) * COS(pos0(3))
      x(2) = pos0(1) * SIN(pos0(2)) * SIN(pos0(3))
      x(3) = pos0(1) * COS(pos0(2))
    END IF

    dx  = pos(1) - x(1)
    dy  = pos(2) - x(2)
    dz  = pos(3) - x(3)

    dxx = dx * dx
    dyy = dy * dy
    dzz = dz * dz
    dxy = dx * dy
    dxz = dx * dz
    dyz = dy * dz

    bases_3D_cubic(1)  = 1.0D0
    bases_3D_cubic(2)  = dx
    bases_3D_cubic(3)  = dy
    bases_3D_cubic(4)  = dz
    bases_3D_cubic(5)  = dxx
    bases_3D_cubic(6)  = dxy
    bases_3D_cubic(7)  = dxz
    bases_3D_cubic(8)  = dyy
    bases_3D_cubic(9)  = dyz
    bases_3D_cubic(10) = dzz
    bases_3D_cubic(11) = dxx * dx
    bases_3D_cubic(12) = dyy * dy
    bases_3D_cubic(13) = dzz * dz
    bases_3D_cubic(14) = dxx * dy
    bases_3D_cubic(15) = dxx * dz
    bases_3D_cubic(16) = dyy * dx
    bases_3D_cubic(17) = dyy * dz
    bases_3D_cubic(18) = dzz * dx
    bases_3D_cubic(19) = dzz * dy
    bases_3D_cubic(20) = dxy * dz

  END SUBROUTINE get_base3D_cubic

END MODULE LRE

