
      PROGRAM TESTSVD

      implicit none

      real(8):: A(4,2),U(4,4),S(2),V(2,2),VT(2,2),SS(4,2)
      real(8):: test1(4,2),test2(4,2),pA(2,4)
      real(8):: cond_num=0.001
      integer:: I,J,K,M,N

       K=2
       DO I=1,4
        DO J=1,2
           A(I,J)=K*1.0D0
           K=K+1
         END DO
       END DO

       CALL svd_pseudo_inverse(A,pA,4,2,cond_num)

      ! SS=0
      ! SS(1,1) = S(1)
      ! SS(2,2) = S(2)

      ! write(*,*) 'SS='
      ! WRITE(*,'(2F10.4)') ((SS(I,J),J=1,2),I=1,4)

      ! test1 = matmul(SS,VT)
      ! test2 = matmul(U,test1)

      ! WRITE(*,*) 'A(ret) = '
!       write(*,*) ' test1='
      ! WRITE(*,'(2F10.4)') ((test2(I,J),J=1,2),I=1,4)

!       write(*,*) ' test2='
!       WRITE(*,'(2F10.4)') ((test2(I,J),J=1,2),I=1,4)

      END PROGRAM



      SUBROUTINE svd_pseudo_inverse(A,pA,m,n,cond_num)
      ! A is the input matrix, which has dimension mxn
      ! cond_num is the largest condition number that A can be 
      ! if cond(A)>cond_num, this subroutine makes cond(A)=cond_num
      ! before calculating the (pseudo)inverse
      ! pA is the pseudo-inverse of A (dimension: nxm)

      real(8), intent(in)   :: A(m,n)
      real(8), intent(in)   :: cond_num
      integer, intent(in)   :: m,n
      real(8), intent(out)  :: pA(n,m)
      real(8) :: U(m,m), UT(m,m), V(n,n), VT(n,n), S_inv(n,m) 
      real(8), allocatable:: S(:), work(:)
      integer  :: info_svd, LWORK, i, len_s
      character:: stringA*50, stringS*50
!      real(8) :: test1(m,n), test2(m,n), SS(m,n)

      WRITE(*,*) 'A ='
      WRITE(stringA, '( "(" ,I4, "f10.4)" )' ) n ! # of column A
      !write(*,*) stringA
      WRITE(*,stringA) ((A(I,J),J=1,n),I=1,m)

      len_s = min(m,n)
      allocate(S(len_s))

      ! calculate the svd of A:

      LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
      allocate(work(LWORK))

      call dgesvd('A','A',m,n,A,m,S,U,m,VT,n,WORK,LWORK,info_svd)

      ! WRITE(*,*) 'U ='
      ! WRITE(*,'(4F10.4)') ((U(I,J),J=1,4),I=1,4)

      ! WRITE(*,*) 'VT='
      ! WRITE(*,'(2F10.4)') ((VT(I,J),J=1,2),I=1,2)
      ! WRITE(*,*) 'S ='
      ! WRITE(*,'(2F10.4)') S

      ! below is trying to get A back based on U,S,VT =======
      ! SS=0
      ! SS(1,1) = S(1)
      ! SS(2,2) = S(2)

      ! write(*,*) 'SS='
      ! WRITE(*,'(2F10.4)') ((SS(I,J),J=1,2),I=1,4)

      ! test1 = matmul(SS,VT)
      ! test2 = matmul(U,test1)

      ! WRITE(*,*) 'A(ret) = '
      ! WRITE(*,'(2F10.4)') ((test2(I,J),J=1,2),I=1,4)
      ! =====================================================

     ! calculate the inverse of A (pA) below:
      V  = transpose(VT)
      UT = transpose(U)

      S_inv = 0
      do i=1,len_s
         if ( S(i).ge.cond_num*S(1) ) then
             S_inv(i,i) = 1./S(i)
         endif
      enddo

      write(*,*) 'pinv(S) = '
      WRITE(stringS, '( "(" ,I4, "f10.4)" )' ) m ! # of column S
      write(*,stringS) ((S_inv(i,j),j=1,m),i=1,n)

      pA = matmul(V,matmul(S_inv,UT))
  
      write(*,*) 'pinv(A) = '
      write(*,stringS) ((pA(i,j),j=1,m),i=1,n)

      END SUBROUTINE
