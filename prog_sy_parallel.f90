PROGRAM prog_sy
	USE mod_sy_proc
	implicit none

	COMPLEX(8),ALLOCATABLE :: Sxyz1(:,:,:),Sxyz2(:,:,:)
	REAL(dp)  ,ALLOCATABLE :: lambda1(:),lambda2(:)
	REAL(dp)               :: v
	REAL(dp),PARAMETER     :: k = 1.   ! Initial value for k_f
	INTEGER ,PARAMETER     :: N = 100  ! NxN S_(x,y,z) operators
	INTEGER                :: i,j,io

!	.. Timing vars ..
	INTEGER                :: it0,it1,rate ! CPU_TIME() unsuitable for parallel runs

!	.. Input files
	CHARACTER(LEN=40)      :: Sxyz1_file,Sxyz2_file

! -- Allocations
	ALLOCATE( Sxyz1(3,N,N),Sxyz2(3,N,N) )
!
! -- If Sxyz1 and 2 files are provided as arguments, read them
	IF ( IARGC() .GT. 0 ) THEN
		CALL GETARG(1,Sxyz1_file)
		OPEN(10,file=Sxyz1_file,status='OLD',action='READ')
		DO i=1,3
			DO j=1,N
				READ(10,*,iostat=io) Sxyz1(i,j,:)
			END DO
		END DO
		CLOSE(10)

		CALL GETARG(2,Sxyz2_file)
		OPEN(10,file=Sxyz2_file,status='OLD',action='READ')
		DO i=1,3
			DO j=1,N
				READ(10,*,iostat=io) Sxyz2(i,j,:)
			END DO
		END DO
		CLOSE(10)
	ELSE
! -- Else, init Sxyz1,2 with random complex values
		DO i=1,3
			CALL INIT_A_RND(N,Sxyz1(i,:,:))
!			WRITE(*,'(/,A,I0,A)') 'Block ',i,' of Sxyz1:' 
!			CALL PRINT_CMAT(Sxyz1(i,:,:))

			CALL INIT_A_RND(N,Sxyz2(i,:,:))
!			WRITE(*,'(/,A,I0,A)') 'Block ',i,' of Sxyz2:' 
!			CALL PRINT_CMAT(Sxyz2(i,:,:))
		END DO
!
	END IF

! -- Init lambda1,2 with random real values
	ALLOCATE( lambda1(N),lambda2(N) )
	CALL RANDOM_NUMBER(lambda1)
	CALL RANDOM_NUMBER(lambda2)
!
! -- Calc singlet yield with parallelised second method
	CALL SYSTEM_CLOCK(count_rate = rate)
	CALL SYSTEM_CLOCK(it0)
	v = evalYield_offdiag2p(k,Sxyz1,lambda1,Sxyz2,lambda2)
	CALL SYSTEM_CLOCK(it1)

	WRITE(*,'(/,A,E10.3)') 'Parallelised offdiag method: v = ',v
	WRITE(*,'(A,E10.3,A,I0)') 'Timing: ',REAL(it1-it0)/REAL(rate), 's for N = ',N
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