PROGRAM prog_sy
	USE mod_rwfile
	USE mod_sy_proc
	implicit none

	COMPLEX(8),ALLOCATABLE :: Sxyz1(:,:,:),Sxyz2(:,:,:)
	REAL(dp)  ,ALLOCATABLE :: lambda1(:),lambda2(:)
	REAL(dp)               :: v
	REAL(dp),PARAMETER     :: k = 1.   ! Initial value for k_f
	INTEGER, PARAMETER     :: N = 250  ! NxN S_(x,y,z) operators
	INTEGER                :: i,j

!	.. Timing vars ..
	INTEGER                :: it0,it1,rate ! CPU_TIME() unsuitable for parallel runs

!	.. Switch ..
	LOGICAL,PARAMETER      :: WRITE_MAT = .TRUE.

!
! ---------- SYSTEM SETUP SECTION: EITHER CREATE RANDOM Sxyz AND lambda, OR READ FROM FILE ----------
!
! -- If Sxyz1 and 2 files are provided as arguments, read them
	IF ( IARGC().GT.0 ) THEN
        CALL READ_3xNx2N_FROM_F(N,1,Sxyz1)
        CALL READ_3xNx2N_FROM_F(N,2,Sxyz2)
	ELSE
! -- Else, init Sxyz1,2 with random complex values
	    ALLOCATE( Sxyz1(3,N,N),Sxyz2(3,N,N) )
		DO i=1,3
			CALL INIT_A_RND(N,Sxyz1(i,:,:))
			CALL INIT_A_RND(N,Sxyz2(i,:,:))
			DO j=1,N
				Sxyz1(i,j,j) = 10.d0 * Sxyz1(i,j,j)
				Sxyz2(i,j,j) = 10.d0 * Sxyz2(i,j,j)
			END DO
		END DO
	END IF
!
! -- If lambda1 and 2 files are provided as arguments, read them
	IF ( IARGC().GT.0 ) THEN
        CALL READ_N_FROM_F(N,3,lambda1)
        CALL READ_N_FROM_F(N,4,lambda2)
	ELSE
! -- Init lambda1,2 with random real values
	    ALLOCATE( lambda1(N),lambda2(N) )
		CALL RANDOM_NUMBER(lambda1)
		CALL RANDOM_NUMBER(lambda2)
	END IF
!
	IF (WRITE_MAT) THEN ! SAVING MATRICES USED TO FILE, FOR BENCHMARKING PURPOSES
        CALL WRITE_3xNxN_TO_F(N,'Sxyz1.dat',Sxyz1)
        CALL WRITE_3xNxN_TO_F(N,'Sxyz2.dat',Sxyz2)

        CALL WRITE_N_TO_F(N,'lambda1.dat',lambda1)
        CALL WRITE_N_TO_F(N,'lambda2.dat',lambda2)
	END IF
!
! ---------- SINGLET YIELD CALCULATION SECTION ----------
!
! -- Calc singlet yield with parallelised second method
	CALL SYSTEM_CLOCK(count_rate = rate)
	CALL SYSTEM_CLOCK(it0)
	v = evalYield_offdiag2p(k,Sxyz1,lambda1,Sxyz2,lambda2)
	CALL SYSTEM_CLOCK(it1)

	WRITE(*,'(/,A,E12.5)') 'Parallelised offdiag method: v = ',v
	WRITE(*,'(A,E10.3,A,I0)') 'Timing: ',REAL(it1-it0)/REAL(rate), 's for N = ',N
!
! -- Deallocations 
	DEALLOCATE(Sxyz1,Sxyz2)
	DEALLOCATE(lambda1,lambda2)

END PROGRAM prog_sy