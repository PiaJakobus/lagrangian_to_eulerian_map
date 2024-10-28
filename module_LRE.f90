MODULE LRE

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


  SUBROUTINE get_base3D_1d(bases_3D_lin)
    !*****************************************
    !                                         *
    ! returns values for all base functions   *
    ! in standard SPH appraoch; PJ 14.10.2024 *
    !                                         *
    !*****************************************
    DOUBLE PRECISION, INTENT(OUT) :: bases_3D_lin(NN3D_1d)
    bases_3D_lin(1) = 1.0D0
  END SUBROUTINE get_base3D_1d


  SUBROUTINE get_base3D_lin(pos,pos0,bases_3D_lin)

    !*****************************************
    !                                        *
    ! returns values for all base functions  *
    ! in 3D and linear basis; SKR 23.11.2018 *
    !                                        *
    !*****************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: pos(3),pos0(3)
    DOUBLE PRECISION, INTENT(OUT):: bases_3D_lin(NN3D_lin)
    bases_3D_lin(1)= 1.0D0
    bases_3D_lin(2)= pos(1) - pos0(1)
    bases_3D_lin(3)= pos(2) - pos0(2)
    bases_3D_lin(4)= pos(3) - pos0(3)
  END SUBROUTINE get_base3D_lin

  
  SUBROUTINE get_base3D_quad(pos,pos0,bases_3D_quad)

    !********************************************
    !                                           *
    ! returns values for all base functions     *
    ! in 3D and quadratic basis; SKR 23.11.2018 *
    !                                           *
    !********************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: pos(3),pos0(3)
    DOUBLE PRECISION, INTENT(OUT) :: bases_3D_quad(NN3D_quad)

    ! auxiliary
    DOUBLE PRECISION dx,dy,dz
    
    dx= pos(1) - pos0(1)
    dy= pos(2) - pos0(2)
    dz= pos(3) - pos0(3)
    
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


  SUBROUTINE get_base3D_cubic(pos,pos0,bases_3D_cubic)

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
    
  END SUBROUTINE get_base3D_cubic
  
  SUBROUTINE get_base3D_MLS_cubic(pos,pos0,bases_3D_cubic)

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
