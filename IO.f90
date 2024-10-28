MODULE IO
  IMPLICIT NONE 
  CONTAINS
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
END MODULE IO
