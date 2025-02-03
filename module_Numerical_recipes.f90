MODULE NR
use iso_fortran_env, only: dp => real64
  !***************************************************
  !                                                  *
  ! collect here subroutines from Numerical recipes  *
  ! or slight adaptions of them;                     *
  ! some are really (ugly!) old-style....            *
  ! SKR 28.05.2017                                   *
  !                                                  *
  !***************************************************

  IMPLICIT NONE


CONTAINS
  
  function inv(A) result(Ainv)
    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
    double precision, dimension(:,:), intent(in) :: A
    double precision, dimension(size(A,1),size(A,2)) :: Ainv

    double precision :: work = -1  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: info!, n
    ! =========================================================
    integer :: i,j, m,n, lda, ldu, ldvt
    double precision :: s(size(A,1)), u(size(A,1),size(A,2)), trace(size(A,1),size(A,2))
    double precision :: err, A_save(size(A,1),size(A,2)),B(size(A,1),size(A,2)),sigma_inv(size(A,1),size(A,2))
    double precision :: vt(size(A,1),size(A,2)), tmp(size(A,1),size(A,2))
    integer :: lwork 
    ! =========================================================
    !External procedures defined in LAPACK
    !external dgetrf
    !external dgetri

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    lda = size(A,1) 
    ldu = size(A,2)
    ldvt = size(A,2)
 
    A_save = A 
    Ainv = A
    n = size(A,1)
    m = size(A,2)
    lwork = 2*MAX(3*MIN(m,n)+MAX(m,n),5*MIN(m,n))
    !lwork = MAX(3*MIN(m,n)+MAX(m,n),5*MIN(m,n))

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call dgetrf(n, n, Ainv, n, ipiv, info)
    ! DGETRI computes the inverse of a matrix using the LU factorization
    call dgetri(n, Ainv, n, ipiv, work, n, info)
    trace = MATMUL(A,Ainv)
    err = 0.D0
    DO i = 1, size(A,1)
      err = err + trace(i,i)
    END DO 
    IF ((info /= 0) .OR. (ABS(err) - size(A,1) > 0.001)) THEN
      call dgesvd('All','All',m,n,A,lda,s,u,ldu,vt,ldvt,work,lwork,info)
      print*, " ---- SVD ---- "
      sigma_inv = 0.d0
      DO i = 1, size(A,1)
        IF (s(i) >= 1e-9) THEN 
          sigma_inv(i,i) = 1.D0 / s(i)
        ELSE 
          sigma_inv(i,i) = 0.D0
        END IF 
      END DO 
      B = matmul(transpose(vt), sigma_inv)
      Ainv = matmul(B, transpose(u))
      !DO i = 1, size(A,1)
      !  DO j = 1,size(A,1)
      !    IF (s(i) >= 1e-23) THEN
      !      tmp(i,j) = vt(i,j) * (1d0/s(i))
      !    ELSE
      !      tmp(j,i) = 0d0
      !    END IF 
      !  END DO 
      !END DO 
      !Ainv = MATMUL(tmp, u)
   END IF 
   if (info /= 0) then
     print*, 'Matrix inversion failed!'
   end if 
  end function inv

DOUBLE PRECISION FUNCTION pythag(a,b)

    !******************************************
    !                                         *
    ! like in "Numerical Recipes", slightly   *
    ! modified;                               *
    ! calculates SQRT(a**2 + b**2) without    *
    ! deconstructive over- or underflow;      *
    ! SKR 13.08.2018                          *
    !                                         *
    !******************************************

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: a,b
    DOUBLE PRECISION absa,absb

    absa= ABS(a)
    absb= ABS(b)
    IF( absa > absb )THEN
       pythag= absa*SQRT(1.0D0 + (absb/absa)**2)
    ELSE
       IF( absb == 0.0D0 )THEN
          pythag= 0.0D0
       ELSE
          pythag= absb*SQRT(1.0D0 + (absa/absb)**2)
       ENDIF
    ENDIF

  END FUNCTION pythag


  SUBROUTINE indexx(n,arr,indx)

    !**********************************************************
    !                                                         *
    !     indexing routine from Numerical Recipes, p.330;     *
    !     SKR 2.2.2005                                        *
    !                                                         *
    !**********************************************************

    IMPLICIT NONE

    INTEGER n,indx(n),M,NSTACK
    DOUBLE PRECISION arr(n)
    PARAMETER (M=7,NSTACK=50)
    INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
    DOUBLE PRECISION a
    DO j=1,n
       indx(j)=j
    ENDDO
    jstack=0
    l=1
    ir=n
1   if(ir-l.lt.M)then
       DO j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          DO i=j-1,1,-1
             if(arr(indx(i)).le.a)goto 2
             indx(i+1)=indx(i)
          ENDDO
          i=0
2         indx(i+1)=indxt
       ENDDO
       if(jstack.eq.0)return
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       k=(l+ir)/2
       itemp=indx(k)
       indx(k)=indx(l+1)
       indx(l+1)=itemp
       if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
       endif
       if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
       endif
       if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
       endif
       i=l+1
       j=ir
       indxt=indx(l)
       a=arr(indxt)
3      continue
       i=i+1
       if(arr(indx(i)).lt.a)goto 3
4      continue
       j=j-1
       if(arr(indx(j)).gt.a)goto 4
       if(j.lt.i)goto 5
       itemp=indx(i)
       indx(i)=indx(j)
       indx(j)=itemp
       goto 3
5      indx(l)=indx(j)
       indx(j)=indxt
       jstack=jstack+2
       if(jstack.gt.NSTACK)THEN
          PRINT*,'NSTACK too small in indexx'
          STOP
       ENDIF
       if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
       else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
       endif
    endif
    goto 1
  END SUBROUTINE indexx


 
  FUNCTION matinv4(A) result(B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    double precision, intent(in) :: A(4,4)   !! Matrix
    double precision           :: B(4,4),unity(4,4)   !! Inverse matrix
    double precision             :: detinv, one, two, three, four 
    ! Calculate the inverse determinant of the matrix
    !detinv = &
    !  1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
    !   - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
    !   + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
    !   - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    one = A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    two =  A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    three =  A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    four =  A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    
    detinv = 1 / (one - two + three - four) 
    !print*, one, two, three, four, detinv
    ! Calculate the inverse of the matrix
    B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    unity = MATMUL(B,A)
    IF (unity(1,1) + unity(2,2) + unity(3,3) + unity(4,4) - 4.D0 > 0.01) THEN
      print*, "UNITY?", unity
      stop
    END IF 
 END FUNCTION
function linspace(start,end,num,endpoint,step) result(samples)
        
        ! PARAMETERS
        real(dp), intent(in) :: start 
            !! The starting value of the sequence.
        real(dp), intent(in) :: end
            !! The end value of the sequence, unless `endpoint` is set to `.false.`. 
            !! In that case, the sequence consists of all but the last of `num + 1` 
            !! evenly spaced samples, so that `end` is excluded. Note that the 
            !! step size changes when `endpoint` is `.false.`.
        integer, intent(in), optional :: num
            !! Number of samples to generate. Default value is 50.
        logical, intent(in), optional :: endpoint
            !! If `.true.`, `end` is the last sample. Otherwise, it is not included. Default is `.true.`.
        real(dp), intent(out), optional :: step
            !! If present, `step` is the size of spacing between samples.

        ! RETURNS
        real(dp), allocatable :: samples(:)
            !! There are `num` equally spaced samples in the closed interval `[start, stop]` or 
            !! the half-open interval `[start, stop)` (depending on whether `endpoint` is `.true.` or `.false.`).

        integer :: num_, i
        logical :: endpoint_
        real(dp) :: step_

        num_ = 50
        if (present(num)) num_ = num

        endpoint_ = .true.
        if (present(endpoint)) endpoint_ = endpoint

        ! find step size
        if (endpoint_) then
            step_ = (end - start)/real(num_-1,dp)
        else
            step_ = (end - start)/real(num_,dp)
        end if

        if (present(step)) step = step_

        allocate(samples(num_))
        do i = 1, num_
            samples(i) = start + (i-1)*step_
        end do
    end function linspace
END MODULE NR
