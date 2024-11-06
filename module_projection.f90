MODULE project_2D
  USE KERNEL, ONLY : Wendland_C6_3D
  IMPLICIT NONE
  CONTAINS


  SUBROUTINE integrate_kernel(r_G, rp, h, rdim, N, f_grid,fp)
    INTEGER, INTENT(IN) :: N, rdim
    DOUBLE PRECISION, INTENT(IN) :: r_G(rdim), rp(rdim,N), h(N),fp(N)
    DOUBLE PRECISION, INTENT(OUT) :: f_grid
    DOUBLE PRECISION:: w_tilde, w,pos_ab(rdim),dz,omega
    DOUBLE PRECISION:: r_ab,norm,r2, r_z,q
    DOUBLE PRECISION, PARAMETER :: R = 2000
    INTEGER :: p,i,zaehler,num_steps
    num_steps = 50
    f_grid = 0.D0 
    DO p = 1, N,100
      pos_ab = ABS(rp(1:2,p) - r_G(1:2))
      r_ab = 0.5D0 * SQRT(pos_ab(1)**2d0 + pos_ab(2)**2d0) 
      dz = 2.D0 * SQRT(R*R - r_ab) / num_steps
      w_tilde = 0.D0
      DO i = 0, num_steps 
        r_z = i * dz
        r2 = r_ab + SQRT(0.5D0 * r_z * r_z)  
        q = SQRT(r2) / h(p)
        CALL Wendland_C6_3D(q,w,zaehler)
        w_tilde = w_tilde + w * dz
      !print*, SQRT(r2)/h(p)
      END DO 
      !print*, w_tilde
      omega =  w_tilde!(m(p) / rho(p)) * w
      norm = norm + omega 
      f_grid = f_grid + omega * fp(p)
    END DO
  END SUBROUTINE



END MODULE project_2D
