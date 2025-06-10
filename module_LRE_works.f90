MODULE LRE
  USE KERNEL, ONLY : Wendland_C6_3D, kernel_WC6

  !********************************************
  !                                           *
  ! collect here routines related to MLS;     *
  !                                           *
  ! main strategy:                            *
  ! =============                             *
  !                                           *
  ! a) moment matrix: after a set of basis    *
  !    functions bas and a weighting function *
  !    W have been chosen, we can assemble    * 
  !    the moment matrix:                     *
  !                                           *
  ! M_ij= sum_b W(v_g)*bas_i(x_b)*bas_j(x_b)  *
  !                                           *
  !   where v_g= |x_g - x_b|/h_g and bas_k    *
  !   is the kth basis function; this matrix  *
  !   is actually INDEPENDENT of the          *
  !   quantities to be approximated and only  *
  !   dependent on the positions, basis       *
  !   functions and weights                   *
  !                                           *
  !                                           *
  ! b) function vector: it contains the       *
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
  ! The solution for the coefficient          *
  ! vector is then:                           *
  !                                           *
  ! a_i= (M_ij)^(-1)*B_j                      *
  !                                           *
  ! and the function approximation is:        *
  !                                           *
  ! A(x)= sum_i a_i bas_i(x)                  *
  !                                           *
  ! NOTE ON MATRIX INVERSION:                 *
  ! inversion is considered "tricky" since    *
  ! matrices are often "nearly singular";     *
  ! use special matrix-inversion routines,    *
  ! e.g. "Numerical Recipes", LAPACK or       *
  ! similar...                                *
  !                                           *
  ! DEGREES OF FREEDOM (DOF):                 *
  ! we need at least as many particles as     *
  ! DOFs; there are                           *
  !                                           *
  ! DOF= (d+m)!/(d! m!),                      *
  !                                           *
  ! where d is the number of spatial          *
  ! dimensions and m the order of the basis   *
  ! functions:                                *
  ! in 3D, linear basis:    DOF=  4           *
  !        quadratic basis: DOF= 10           *
  !        cubic basis:     DOF= 20           *
  !                                           *
  ! SKR 27.01.2020                            *
  !                                           *
  !********************************************
  USE NR
  IMPLICIT NONE

  INTEGER, PARAMETER :: NN3D_1d =   1
  INTEGER, PARAMETER :: NN3D_lin=   4
  INTEGER, PARAMETER :: NN3D_quad= 10
  INTEGER, PARAMETER :: NN3D_cub=  20

  
  CONTAINS

     SUBROUTINE assemble_moment_matrix(N,INT_param,moment_matrix,rp,r_G,rdim,h,rho,m,polar,zaehler)

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
    INTEGER, INTENT(INOUT) :: zaehler
    DOUBLE PRECISION, INTENT(IN) :: rp(rdim, N), r_G(rdim), h(N),rho(N),m(N)
    DOUBLE PRECISION, INTENT(OUT) :: moment_matrix(INT_param, INT_param)
    DOUBLE PRECISION :: P_i(INT_param), w,w_tilde
    DOUBLE PRECISION :: pos_ab(rdim), q
    INTEGER :: i,p,k
    LOGICAL, INTENT(IN) :: polar
    moment_matrix = 0d0
    zaehler = 0
    DO p = 1, N
      IF (polar .eqv. .false.) THEN 
        pos_ab = ABS(rp(:,p) - r_G) 
        q = 0.5D0 * SQRT(pos_ab(1)**2d0 + pos_ab(2)**2d0 + pos_ab(3)**2d0) / h(p) 
      ELSE 
        CALL distance_polar(rp(:,p),r_G,q)
        q = 0.5D0 * q / h(p)
      END IF 
      IF (q <= 1.D0) THEN 
        CALL Wendland_C6_3D(q,w,zaehler)
 
        IF (INT_param == NN3D_quad) THEN 
          CALL get_base3D_quad(rp(:,p),r_G,P_i,polar)
        ELSE IF (INT_param == NN3D_lin) THEN 
          CALL get_base3D_lin(rp(:,p),r_G,P_i,polar)
        ELSE IF (INT_param == NN3D_cub) THEN
          CALL get_base3D_cubic(rp(:,p),r_G,P_i,polar)
        ELSE IF (INT_param == NN3D_1d) THEN 
          P_i = (/1d0/)
        END IF 
        w_tilde = w * m(p) / (rho(p) * h(p)**3.D0)
        DO k = 1, INT_param
          DO i = 1, INT_param
            !CALL kernel(w,rdim,h,r_G,rp(:,p))
            moment_matrix(i,k) = moment_matrix(i,k) + P_i(i) * P_i(k) * w_tilde
          END DO
        END DO
      END IF 
    END DO
   !print*, ":::::::::::::: particles: ", zaehler 
  END SUBROUTINE assemble_moment_matrix

  SUBROUTINE get_vector_b(B,fp,N,rp,r_G,rdim,INT_param,h,rho,m,polar)

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
    DOUBLE PRECISION, INTENT(IN) :: rp(rdim,N), r_G(rdim), fp(N), h(N),rho(N),m(N)
    DOUBLE PRECISION, INTENT(OUT) :: B(INT_param)
    DOUBLE PRECISION :: P_i(INT_param), w, w_tilde
    DOUBLE PRECISION :: q, pos_ab(rdim)
    LOGICAL, INTENT(IN) :: polar
    INTEGER :: zaehler 
    zaehler = 0 
    B = 0.
    DO p = 1, N
      IF (polar .eqv. .false.) THEN
        pos_ab = (rp(:,p) - r_G) 
        q = 0.5D0 * SQRT(pos_ab(1)**2d0 + pos_ab(2)**2d0 + pos_ab(3)**2d0) / h(p)
      ELSE 
        CALL distance_polar(rp(:,p),r_G,q)
        q = 0.5d0 * q  / h(p)
      END IF 
      IF (q <= 1.D0) THEN 
        CALL Wendland_C6_3D(q,w,zaehler)
        !CALL kernel_WC6(v,w,zaehler)
        IF (INT_param == NN3D_lin) THEN 
          CALL get_base3D_lin(rp(:,p),r_G,P_i,polar)
          !CALL kernel(w,rdim,h,r_G,rp(:,p))
        ELSE IF (INT_param == NN3D_quad) THEN
          CALL get_base3D_quad(rp(:,p),r_G,P_i,polar)
        ELSE IF (INT_param == NN3D_cub) THEN 
          CALL get_base3D_cubic(rp(:,p),r_G,P_i,polar)
        ELSE IF (INT_param == NN3D_1d) THEN 
          P_i = (/1d0/)
        END IF 
        w_tilde = w * m(p) / (rho(p) * h(p)**3.D0) 
        DO k = 1, INT_param
          B(k) = B(k) + fp(p) * P_i(k) * w_tilde
        END DO
      END IF 
   END DO 
  END SUBROUTINE get_vector_b


  PURE SUBROUTINE distance_polar(xp,rg,dist)
    DOUBLE PRECISION, INTENT(INOUT) :: dist
    DOUBLE PRECISION, INTENT(IN) :: xp(3), rg(3)
    DOUBLE PRECISION :: x, y, z
    x = rg(1) * sin(rg(2)) * cos(rg(3))
    y = rg(1) * sin(rg(2)) * sin(rg(3))
    z = rg(1) * cos(rg(2))
    dist = SQRT((xp(1)-x)**2.d0 + (xp(2)-y)**2.d0 + (xp(3)-z)**2.d0)
  END SUBROUTINE distance_polar


  SUBROUTINE standard_SPH(fp,N,rp,r_G,rdim,h,rho,m,f_grid,polar,zaehler)
    INTEGER, INTENT(IN) :: N, rdim
    LOGICAL, INTENT(IN) :: polar 
    DOUBLE PRECISION, INTENT(IN) :: fp(N), rho(N), m(N), rp(rdim,N), r_G(rdim),h(N)
    DOUBLE PRECISION, INTENT(OUT) :: f_grid
    INTEGER, INTENT(INOUT) :: zaehler
    DOUBLE PRECISION :: w, q, pos_ab(rdim), norm, omega
    INTEGER :: p
    norm = 1e-40
    f_grid = 0.d0
    zaehler = 0
    DO p = 1, N
      IF (polar .eqv. .false.) THEN
        pos_ab = (rp(:,p) - r_G)
        q = 0.5D0 * SQRT(pos_ab(1)**2d0 + pos_ab(2)**2d0 + pos_ab(3)**2d0) / h(p)
      ELSE 
        CALL distance_polar(rp(:,p),r_G,q)
        q = 0.5d0 * q  / h(p)
      END IF 
      omega = 0.D0
      !q = SQRT(pos_ab(1)**2d0 + pos_ab(2)**2d0 + pos_ab(3)**2d0) / h(p)
      IF (q <= 1.D0) THEN
        CALL Wendland_C6_3D(q,w,zaehler)
        omega = (m(p) / rho(p)) * (w / h(p)**3.D0)
        !omega = (m(p) / rho(p)) * w
        norm = norm + omega
      END IF 
      f_grid = f_grid + omega * fp(p)
    END DO 
    !print*, "Zaehler SPH", zaehler
    f_grid = f_grid / norm 
  END SUBROUTINE standard_SPH

  SUBROUTINE algebra(INT_param,moment_matrix,B,f_grid)

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

    INTEGER, INTENT(IN) :: INT_param
    DOUBLE PRECISION, INTENT(IN) :: B(INT_param)
    DOUBLE PRECISION, INTENT(IN) :: moment_matrix(INT_param, INT_param)
    DOUBLE PRECISION, INTENT(OUT) :: f_grid
    !DOUBLE PRECISION :: pl(INT_param)
    DOUBLE PRECISION :: A_inv(INT_param,INT_param), beta(INT_param)
    f_grid = 0d0
    A_inv = 0d0
    IF (INT_param == NN3D_1d) THEN
      IF (ABS(moment_matrix(1,1)) > 0d0) THEN 
        A_inv = 1d0 / moment_matrix
      ELSE
        A_inv = 0d0 
      END IF 
    ELSE 
      A_inv = inv(moment_matrix)
      !A_inv = matinv4(moment_matrix)
    END IF
    beta = MATMUL(A_inv,B)
    !IF (INT_param == NN3D_lin) THEN 
    !  CALL get_base3D_lin(r,r_G,pl,polar)
    !ELSE IF (INT_param == NN3D_quad) THEN 
    !  CALL get_base3D_quad(r,r_G,pl,polar)
    !ELSE IF (INT_param == NN3D_cub) THEN 
    !  CALL get_base3D_cubic(r,r_G,pl,polar)
    !ELSE IF (INT_param == NN3D_1d) THEN 
    !  pl = (/1d0/)
    !END IF 
    !f_grid = DOT_PRODUCT(beta,pl)
    f_grid = beta(1)
    !stop
  END SUBROUTINE algebra


  PURE SUBROUTINE get_base3D_1d(bases_3D_lin)
    !*****************************************
    !                                         *
    ! returns values for all base functions   *
    ! in standard SPH appraoch; PJ 14.10.2024 *
    !                                         *
    !*****************************************
    DOUBLE PRECISION, INTENT(OUT) :: bases_3D_lin(NN3D_1d)
    bases_3D_lin(1) = 1.0D0
  END SUBROUTINE get_base3D_1d


  PURE SUBROUTINE get_base3D_lin(pos,pos0,bases_3D_lin,polar)

    !*****************************************
    !                                        *
    ! returns values for all base functions  *
    ! in 3D and linear basis; SKR 23.11.2018 *
    !                                        *
    !*****************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: pos(3),pos0(3)
    LOGICAL, INTENT(IN) :: polar 
    DOUBLE PRECISION, INTENT(OUT):: bases_3D_lin(NN3D_lin)
    DOUBLE PRECISION :: x(3)
    bases_3D_lin(1)= 1.0D0
    x = pos0
    IF (polar .eqv. .true.) THEN 
      x(1) = pos0(1) * sin(pos0(2)) * cos(pos0(3))
      x(2) = pos0(1) * sin(pos0(2)) * sin(pos0(3))
      x(3) = pos0(1) * cos(pos0(2))
    END IF 
    bases_3D_lin(2)= pos(1) - x(1)
    bases_3D_lin(3)= pos(2) - x(2)
    bases_3D_lin(4)= pos(3) - x(3)
  END SUBROUTINE get_base3D_lin


  SUBROUTINE get_base3D_quad(pos,pos0,bases_3D_quad,polar)

    !********************************************
    !                                           *
    ! returns values for all base functions     *
    ! in 3D and quadratic basis; SKR 23.11.2018 *
    !                                           *
    !********************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: pos(3),pos0(3)
    DOUBLE PRECISION, INTENT(OUT) :: bases_3D_quad(NN3D_quad)
    LOGICAL, INTENT(IN) :: polar 
    ! auxiliary
    DOUBLE PRECISION dx,dy,dz
    DOUBLE PRECISION :: x(3)
    x = pos0
    IF (polar .eqv. .true.) THEN 
      x(1) = pos0(1) * sin(pos0(2)) * cos(pos0(3))
      x(2) = pos0(1) * sin(pos0(2)) * sin(pos0(3))
      x(3) = pos0(1) * cos(pos0(2))
    END IF 
    !print*, "x(1)-x(3) ", x(1), x(2),x(3)
 
    
    dx= pos(1) - x(1)
    dy= pos(2) - x(2)
    dz= pos(3) - x(3)
    
    bases_3D_quad(1)=  1.0D0
    bases_3D_quad(2)=  dx
    bases_3D_quad(3)=  dy
    bases_3D_quad(4)=  dz
    bases_3D_quad(5)=  dx*dx
    bases_3D_quad(6)=  dx*dy
    bases_3D_quad(7)=  dx*dz
    bases_3D_quad(8)=  dy*dy
    bases_3D_quad(9)=  dy*dz
    bases_3D_quad(10)= dz*dz
    
  END SUBROUTINE get_base3D_quad


  PURE SUBROUTINE get_base3D_cubic(pos,pos0,bases_3D_cubic,polar)

    !****************************************
    !                                       *
    ! returns values for all base functions *
    ! in 3D and cubic basis; SKR 23.11.2018 *
    !                                       *
    !****************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: pos(3),pos0(3)
    LOGICAL, INTENT(IN) :: polar 
    DOUBLE PRECISION, INTENT(OUT) :: bases_3D_cubic(NN3D_cub)           

    ! auxiliary
    DOUBLE PRECISION dx,dy,dz,dxx,dyy,dzz,dxy,dxz,dyz
    DOUBLE PRECISION :: x(3)
    x = pos0
    IF (polar .eqv. .true.) THEN 
      x(1) = pos0(1) * sin(pos0(2)) * cos(pos0(3))
      x(2) = pos0(1) * sin(pos0(2)) * sin(pos0(3))
      x(3) = pos0(1) * cos(pos0(2))
    END IF 
    
    dx=  pos(1) - x(1)
    dy=  pos(2) - x(2)
    dz=  pos(3) - x(3)

    dxx= dx*dx
    dyy= dy*dy
    dzz= dz*dz
    dxy= dx*dy
    dxz= dx*dz
    dyz= dy*dz

    ! now the function...
    bases_3D_cubic(1)=  1.0D0
    bases_3D_cubic(2)=  dx
    bases_3D_cubic(3)=  dy
    bases_3D_cubic(4)=  dz
    bases_3D_cubic(5)=  dxx
    bases_3D_cubic(6)=  dxy
    bases_3D_cubic(7)=  dxz
    bases_3D_cubic(8)=  dyy
    bases_3D_cubic(9)=  dyz
    bases_3D_cubic(10)= dzz
    bases_3D_cubic(11)= dxx*dx
    bases_3D_cubic(12)= dyy*dy
    bases_3D_cubic(13)= dzz*dz
    bases_3D_cubic(14)= dxx*dy
    bases_3D_cubic(15)= dxx*dz
    bases_3D_cubic(16)= dyy*dx
    bases_3D_cubic(17)= dyy*dz
    bases_3D_cubic(18)= dzz*dx
    bases_3D_cubic(19)= dzz*dy
    bases_3D_cubic(20)= dxy*dz
    
  END SUBROUTINE get_base3D_cubic
  
  PURE SUBROUTINE get_base3D_MLS_cubic(pos,pos0,bases_3D_cubic)

    !****************************************
    !                                       *
    ! returns values for all base functions *
    ! in 3D and cubic basis; SKR 23.11.2018 *
    !                                       *
    !****************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: pos(3),pos0(3)
    DOUBLE PRECISION, INTENT(OUT) :: bases_3D_cubic(NN3D_cub)           

    ! auxiliary
    DOUBLE PRECISION dx,dy,dz,dxx,dyy,dzz,dxy,dxz,dyz
    
    dx=  pos(1) - pos0(1)
    dy=  pos(2) - pos0(2)
    dz=  pos(3) - pos0(3)

    dxx= dx*dx
    dyy= dy*dy
    dzz= dz*dz
    dxy= dx*dy
    dxz= dx*dz
    dyz= dy*dz

    ! now the function...
    bases_3D_cubic(1)=  1.0D0
    bases_3D_cubic(2)=  dx
    bases_3D_cubic(3)=  dy
    bases_3D_cubic(4)=  dz
    bases_3D_cubic(5)=  dxx
    bases_3D_cubic(6)=  dxy
    bases_3D_cubic(7)=  dxz
    bases_3D_cubic(8)=  dyy
    bases_3D_cubic(9)=  dyz
    bases_3D_cubic(10)= dzz
    bases_3D_cubic(11)= dxx*dx
    bases_3D_cubic(12)= dyy*dy
    bases_3D_cubic(13)= dzz*dz
    bases_3D_cubic(14)= dxx*dy
    bases_3D_cubic(15)= dxx*dz
    bases_3D_cubic(16)= dyy*dx
    bases_3D_cubic(17)= dyy*dz
    bases_3D_cubic(18)= dzz*dx
    bases_3D_cubic(19)= dzz*dy
    bases_3D_cubic(20)= dxy*dz
    
  END SUBROUTINE get_base3D_MLS_cubic
    
END MODULE LRE
