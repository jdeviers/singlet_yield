PROGRAM prog_sy_timing
	USE mod_rwfile
	USE mod_sy_proc
	implicit none

	COMPLEX(8),ALLOCATABLE :: Sxyz1(:,:,:),Sxyz2(:,:,:)
	REAL(dp)  ,ALLOCATABLE :: lambda1(:),lambda2(:)
	REAL(dp)               :: v
	REAL(dp),PARAMETER     :: k = 1. ! Initial value for k_f
	INTEGER(8)             :: N = 50  ! NxN S_(x,y,z) operators
	INTEGER(8)             :: i

!	.. Timing vars ..
	INTEGER(8)               :: it0,it1,rate ! CPU_TIME() unsuitable for parallel runs
	REAL(dp)               :: t_1,t_2p,t_2s


	OPEN(10,file='timings.dat',status='unknown',action='write',position='append')
	DO N=5,250,5
! -- Allocations
		ALLOCATE( Sxyz1(3,N,N),Sxyz2(3,N,N) )
!

! -- Init Sxyz1,2 with random complex values
		DO i=1,3
			CALL INIT_A_RND(N,Sxyz1(i,:,:))
!			WRITE(*,'(/,A,I0,A)') 'Block ',i,' of Sxyz1:' 
!			CALL PRINT_CMAT(Sxyz1(i,:,:))

			CALL INIT_A_RND(N,Sxyz2(i,:,:))
!			WRITE(*,'(/,A,I0,A)') 'Block ',i,' of Sxyz2:' 
!			CALL PRINT_CMAT(Sxyz2(i,:,:))
		END DO
!

! -- Init lambda1,2 with random real values
		ALLOCATE( lambda1(N),lambda2(N) )
		CALL RANDOM_NUMBER(lambda1)
		CALL RANDOM_NUMBER(lambda2)
!
! -- Calc singlet yield with first method
		CALL SYSTEM_CLOCK(count_rate = rate)
		CALL SYSTEM_CLOCK(it0)
		v = evalYield(k,Sxyz1,lambda1,Sxyz2,lambda2)
		CALL SYSTEM_CLOCK(it1)

		t_1 = REAL(it1-it0)/REAL(rate)

!		WRITE(*,'(/,A,E10.3)') 'First method: v = ',v
!		WRITE(*,'(A,E10.3,A,I0)') 'Timing: ',REAL(it1-it0)/REAL(rate), 's for N = ',N
!
! -- Calc singlet yield with parallelised second method
		CALL SYSTEM_CLOCK(count_rate = rate)
		CALL SYSTEM_CLOCK(it0)
		v = evalYield_offdiag2p(k,Sxyz1,lambda1,Sxyz2,lambda2)
		CALL SYSTEM_CLOCK(it1)

		t_2p = REAL(it1-it0)/REAL(rate)

!		WRITE(*,'(/,A,E10.3)') 'Parallelised offdiag method: v = ',v
!		WRITE(*,'(A,E10.3,A,I0)') 'Timing: ',REAL(it1-it0)/REAL(rate), 's for N = ',N
!

! -- Calc singlet yield with serial second method
		CALL SYSTEM_CLOCK(count_rate = rate)
		CALL SYSTEM_CLOCK(it0)
		v = evalYield_offdiag2p_serial(k,Sxyz1,lambda1,Sxyz2,lambda2)
		CALL SYSTEM_CLOCK(it1)

		t_2s = REAL(it1-it0)/REAL(rate)

!		WRITE(*,'(/,A,E10.3)') 'Serial offdiag method: v = ',v
!		WRITE(*,'(A,E10.3,A,I0)') 'Timing: ',REAL(it1-it0)/REAL(rate), 's for N = ',N
!

! -- Deallocations 
		DEALLOCATE(Sxyz1,Sxyz2)
		DEALLOCATE(lambda1,lambda2)

		WRITE(10,'(I4,3F12.2)') N,t_1,t_2s,t_2p

	END DO
	CLOSE(10)

END PROGRAM prog_sy_timing
