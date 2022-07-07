
      PROGRAM TESTSVD
!
!
! Program to test (SVD) singular value decomposition.
!
      use svd_mod

      interface
         subroutine SVD(A,U,S,V,M,N)
            DOUBLE PRECISION  :: A
            integer           :: M,N
            DOUBLE PRECISION  :: U,S,V
         end subroutine SVD
      end interface

!      use svd_mod

      DOUBLE PRECISION A(4,2),U(4,4),S(2),V(2,2)
      integer:: I,J,K
       K=1
       DO I=1,4
        DO J=1,2
           A(I,J)=K*1.0D0
           K=K+1
         END DO
       END DO
       WRITE(*,*) 'A ='
       WRITE(*,'(2F10.4)') ((A(I,J),J=1,2),I=1,4)

       CALL SVD(A,U,S,V,4,2)

       WRITE(*,*) 'U ='
       WRITE(*,'(4F10.4)') ((U(I,J),J=1,4),I=1,4)

       WRITE(*,*) 'V='
       WRITE(*,'(2F10.4)') V
       WRITE(*,*) 'S ='
       WRITE(*,'(2F10.4)')S

      END PROGRAM

