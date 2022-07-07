
      PROGRAM TESTSVD

      implicit none

      real(8):: A(4,2),U(4,4),S(2),V(2,2),VT(2,2),SS(4,2)
      real(8):: test1(4,2),test2(4,2),pA(2,4)
      integer:: I,J,K,M,N
      real(8) ,allocatable:: work(:)
      integer  :: info_svd, LWORK

       K=1
       DO I=1,4
        DO J=1,2
           A(I,J)=K*1.0D0
           K=K+1
         END DO
       END DO

       WRITE(*,*) 'A ='
       WRITE(*,'(2F10.4)') ((A(I,J),J=1,2),I=1,4)

      M=4
      N=2

      LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
      allocate(work(LWORK))

      call dgesvd('A','A',4,2,A,4,S,U,4,VT,2,WORK,LWORK,info_svd)

      V=transpose(VT)

       WRITE(*,*) 'U ='
       WRITE(*,'(4F10.4)') ((U(I,J),J=1,4),I=1,4)

       WRITE(*,*) 'V='
       WRITE(*,'(2F10.4)') ((V(I,J),J=1,2),I=1,2)
       WRITE(*,*) 'S ='
       WRITE(*,'(2F10.4)') S

      ! CALL svd_pseudo_inverse(A,pA,4,2,0.001)

       SS=0
       SS(1,1) = S(1)
       SS(2,2) = S(2)

       write(*,*) 'SS='
       WRITE(*,'(2F10.4)') ((SS(I,J),J=1,2),I=1,4)

       test1 = matmul(SS,VT)
       test2 = matmul(U,test1)

       WRITE(*,*) 'A(ret) = '
       WRITE(*,'(2F10.4)') ((test2(I,J),J=1,2),I=1,4)


      END PROGRAM
