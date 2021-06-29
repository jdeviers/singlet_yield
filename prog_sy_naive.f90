PROGRAM prog_sy
	USE mod_sy_proc
	implicit none

	COMPLEX(8),ALLOCATABLE :: Sxyz1(:,:,:),Sxyz2(:,:,:)
	REAL(dp)  ,ALLOCATABLE :: lambda1(:),lambda2(:)
	REAL(dp)               :: v
	REAL(dp),PARAMETER     :: k = 1.
	INTEGER ,PARAMETER     :: N = 5
	INTEGER                :: i

!	.. Timing vars ..
	REAL(dp)               :: t0,t1


! -- Allocations
	ALLOCATE( Sxyz1(3,N,N),Sxyz2(3,N,N) )
!

! -- Init Sxyz1,2 with random complex values
	DO i=1,3
		CALL INIT_A_RND(N,Sxyz1(i,:,:))
!		WRITE(*,'(/,A,I0,A)') 'Block ',i,' of Sxyz1:' 
!		CALL PRINT_CMAT(Sxyz1(i,:,:))

		CALL INIT_A_RND(N,Sxyz2(i,:,:))
!		WRITE(*,'(/,A,I0,A)') 'Block ',i,' of Sxyz2:' 
!		CALL PRINT_CMAT(Sxyz2(i,:,:))
	END DO
!

! -- Init lambda1,2 with random complex values
	ALLOCATE( lambda1(N),lambda2(N) )
	CALL RANDOM_NUMBER(lambda1)
	CALL RANDOM_NUMBER(lambda2)
!

! -- Calc singlet yield
	CALL CPU_TIME(t0)
	v = evalYield(k,Sxyz1,lambda1,Sxyz2,lambda2)
	CALL CPU_TIME(t1)

	WRITE(*,'(/,A,E10.3)') 'v = ',v
	WRITE(*,'(A,E10.3,A,I0)') 'Timing: ',t1-t0, 's for N = ',N
!

! -- Deallocations 
	DEALLOCATE(Sxyz1,Sxyz2)
	DEALLOCATE(lambda1,lambda2)

	contains

	SUBROUTINE INIT_A_RND(N,A) ! Creates Hermitian matrix with random coeffs
		implicit none

!	.. Parameters ..
		COMPLEX(8),INTENT(OUT)  :: A(:,:)
		INTEGER,   INTENT(IN)   :: N
!	.. Local arrays ..
		REAL(dp),DIMENSION(N,N) :: Re,Im
		INTEGER                 :: i,j

		CALL RANDOM_NUMBER(Re)
		CALL RANDOM_NUMBER(Im)
		DO i=1,N
			DO j=i,N
				A(i,j) = COMPLEX( Re(i,j),Im(i,j) )
				A(j,i) = CONJG( A(i,j) )
			END DO
		END DO

	END SUBROUTINE INIT_A_RND

	SUBROUTINE PRINT_CMAT(M)
		implicit none

!	.. Parameters ..
		COMPLEX(8),INTENT(IN) :: M(:,:)
!	.. Local scalars ..
		INTEGER               :: i,j

		DO i = 1,UBOUND(M,1)
			WRITE(*,FMT_SCMPLX) (M(i,j), j=1,UBOUND(M,2))
		END DO

	END SUBROUTINE PRINT_CMAT


END PROGRAM prog_sy