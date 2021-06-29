PROGRAM prog_sy_naive
	USE mod_sy_naive_proc
	implicit none

	COMPLEX(8),ALLOCATABLE :: Sxyz1(:,:,:),Sxyz2(:,:,:)
	REAL(dp)  ,ALLOCATABLE :: lambda1(:),lambda2(:)
	REAL(dp)               :: v
	REAL(dp),PARAMETER     :: k = 1.
	INTEGER ,PARAMETER     :: N = 5
	INTEGER                :: i,j


	ALLOCATE( Sxyz1(3,N,N),Sxyz2(3,N,N))
	WRITE(11,*) Sxyz1,Sxyz2

	DO i=1,3
		CALL INIT_A_RND(N,Sxyz1(i,:,:))
		WRITE(*,'(/,A,I0,A)') 'Block ',i,' of Sxyz1:' 
		CALL PRINT_CMAT(Sxyz1(i,:,:))

		CALL INIT_A_RND(N,Sxyz2(i,:,:))
		WRITE(*,'(/,A,I0,A)') 'Block ',i,' of Sxyz2:' 
		CALL PRINT_CMAT(Sxyz2(i,:,:))
	END DO

	ALLOCATE( lambda1(N),lambda2(N) )
	CALL RANDOM_NUMBER(lambda1)
	CALL RANDOM_NUMBER(lambda2)

	v = evalYield(k,Sxyz1,lambda1,Sxyz2,lambda2)
	WRITE(*,'(/,A,E10.3)') 'v = ',v

	DEALLOCATE(Sxyz1,Sxyz2)
	DEALLOCATE(lambda1,lambda2)

	contains

	SUBROUTINE INIT_A_RND(N,A) ! Creates square complex matrix with random coeffs
		implicit none

!	.. Parameters ..
		COMPLEX(8),INTENT(OUT) :: A(:,:)
		INTEGER,   INTENT(IN)  :: N

!	.. Local scalars ..
		INTEGER,PARAMETER      :: IDIST = 1
		INTEGER                :: i,ISEED(4)

		DO i=1,N
			CALL ZLARNV(IDIST,ISEED,N,A(:,i))
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


END PROGRAM prog_sy_naive