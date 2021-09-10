MODULE random_sampling ! Can also be called from prog_sy_parallel
  USE mod_rwfile
  implicit none

  DOUBLE PRECISION :: last_ravg   ! ravg = running average
  INTEGER(8)       :: last_rcount ! rcount = running count

  contains

! ----------

  REAL(dp) FUNCTION R_S(threshold,k,Sxyz1,lambda1,Sxyz2,lambda2) 
    USE OMP_LIB

!	.. Arguments ..
		COMPLEX(8),ALLOCATABLE,INTENT(IN) :: Sxyz1(:,:,:),Sxyz2(:,:,:)
		REAL(dp),  ALLOCATABLE,INTENT(IN) :: lambda1(:),lambda2(:)
		REAL(dp)              ,INTENT(IN) :: k,threshold
!	.. Local scalars ..
		INTEGER(8)           :: N,Z,thread_count,current_rcount,a
    INTEGER(16)          :: N_max
    INTEGER              :: a1,a2,b1,b2,thread_id
		REAL(dp)             :: thread_sum,current_rsum, current_ravg,last_ravg
    REAL(dp)             :: dla,dlb
!	.. Local tensors ..
    REAL(dp),ALLOCATABLE :: indices(:,:)


    N = UBOUND(Sxyz1,2); Z = FLOOR( UBOUND(Sxyz1,1)*UBOUND(Sxyz2,2)/4. )
    current_rcount = 0; current_rsum = 0.d0; last_ravg = 0.d0
    thread_count = 0; thread_sum = 0.d0

    ALLOCATE( indices(N*N,4) )

    !$OMP PARALLEL PRIVATE(thread_id,thread_count,thread_sum) SHARED(indices,current_rcount,current_rsum,current_ravg)
    DO ! Threshold condition on diff b/w current and previous cumulated average
      thread_id = OMP_GET_THREAD_NUM()
      ! Thread 0 regenerates list of indices to sample 
      IF (thread_id .EQ. 0) THEN
        CALL RANDOM_NUMBER(indices)
        indices = 1_4+indices*(N-1_4) ! Prevents an index to be =0
      END IF
        
      !$OMP DO
      DO a = 1,N*N
        a1 = INT(indices(a,1)); a2 = INT(indices(a,2))
        b1 = INT(indices(a,3)); b2 = INT(indices(a,4))
        dla = lambda1(a1) - lambda1(a2); dlb = lambda2(b1) - lambda2(b2)

        ! Skip diagonal elements
        IF ( a1.EQ.a2 ) THEN
            IF ( b1.EQ.b2 ) CYCLE ! Nested IF avoid performing too many checks at every step
        END IF

        ! -- Update the running sum of Ps products
        thread_sum = thread_sum + Ps2_kernel(Sxyz1(:,a1,a2),Sxyz2(:,b1,b2),k,dla,dlb)
        ! -- Update the count of actually explored combination (thread_count .LE. N)
        thread_count = thread_count+1_8
      END DO
      !$OMP END DO

      !$OMP CRITICAL
        current_rsum   = current_rsum   + thread_sum
        current_rcount = current_rcount + thread_count

        IF (thread_id .EQ. 0) THEN
          current_ravg = current_rsum/current_rcount
          WRITE(10,*) ( (current_rsum/current_rcount)/(k*k / Z) ),current_rsum,current_rcount
        END IF

        IF ( ABS(current_ravg - last_ravg) .LE. threshold) EXIT ! All threads exit, not just 0 -- not sure it makes a difference
        
        IF (thread_id .EQ. 0) THEN
          last_ravg = current_ravg
        END IF
      !$OMP END CRITICAL

    END DO
    !$OMP END PARALLEL
    DEALLOCATE(indices) 


    ! Final normalisation and return:
    R_S = (current_rsum / current_rcount) / (k*k / Z)

  END FUNCTION R_S

! ----------

  REAL(dp) FUNCTION Ps2_kernel(S1,S2,k,dla,dlb)

    COMPLEX(8) :: S1(3),S2(3)
    REAL(dp)   :: k,dla,dlb

    Ps2_kernel = ( ABS(S1(1)*S2(1) + S1(2)*S2(2) + S1(3)*S2(3))**2.d0 / (k*k + (dla + dlb)**2.d0) )

  END FUNCTION Ps2_kernel

! ----------

END MODULE random_sampling
