MODULE random_sampling ! Can be called from prog_sy_parallel
  implicit none

	INTEGER(4),PARAMETER ::          &
		sp = SELECTED_REAL_KIND(6,37), & ! 32-bits precision
		dp = SELECTED_REAL_KIND(15,307)  ! 64-bits precision

	CHARACTER(LEN=*),PARAMETER ::                          &
		FMT_SCMPLX='(*( "(" SP F6.3 X SP F6.3 " i)" 2X ))',  &
		FMT_LCMPLX='(*( "(" SP E10.3 X SP E10.3 " i)" 2X ))'

  DOUBLE PRECISION :: last_ravg   ! ravg = running average
  INTEGER          :: last_rcount ! rcount = running count

  contains

! ----------

  REAL(dp) FUNCTION R_S(threshold,k,Sxyz1,lambda1,Sxyz2,lambda2) 
    USE OMP_LIB

!	.. Arguments ..
		COMPLEX(8),ALLOCATABLE,INTENT(IN) :: Sxyz1(:,:,:),Sxyz2(:,:,:)
		REAL(dp),  ALLOCATABLE,INTENT(IN) :: lambda1(:),lambda2(:)
		REAL(dp)              ,INTENT(IN) :: k,threshold
!	.. Local scalars ..
		INTEGER                :: N,thread_count,running_count
    INTEGER                :: a1,a2,b1,b2
		REAL(dp)               :: thread_sum,running_sum
    REAL(dp)               :: dla,dlb
!	.. Local tensors ..
    LOGICAL,ALLOCATABLE    :: lookup_1(:,:),lookup_2(:,:)
    REAL,   ALLOCATABLE    :: indices(:,:)

    N = UBOUND(Sxyz1,2)
    ALLOCATE( lookup1(N,N),lookup2(N,N) )
    ALLOCATE( indices(N,4) )
    lookup1 = .TRUE.; lookup2 = .TRUE.
    running_count = 0; running_sum = 0.d0

    DO ! Threshold condition on diff b/w current and previous cumulated average
      ! Regenerate list of indices to sample 
      CALL RANDOM_NUMBER(indices)
      indices = INT(indices*N)

      !$OMP PARALLEL PRIVATE(thread_count,thread_sum) SHARED(lookup1,lookup2,running_count,running_sum)
        thread_count = 0; thread_sum = 0.d0
        
        !$OMP DO
        DO a = 1,N

          a1 = indices(a,1); a2 = indices(a,2)
          b1 = indices(a,3); b2 = indices(a,4)
          dla = lambda(a1) - lambda(a2); dlb = lambda(b1) - lambda(b2)

          ! Skip diagonal elements
          IF ( a1.EQ.a2 ) THEN
            IF ( a2.EQ.b1 ) THEN
              IF ( b1.EQ.b2 ) CYCLE ! Nested IF avoid performing too many checks at every step
            END IF
          END IF

          ! Skip doubly .FALSE. combinations
          IF ( .NOT.(lookup1(a1,a2).AND.lookup2(b1,b2)) ) CYCLE

          ! Calc sums on unexplored indices combinations:
          ! -- Mark the index combination as explored:
          lookup1( a1,a2 ) = .NOT.lookup1( a1,a2 )
          lookup2( b1,b2 ) = .NOT.lookup2( b1,b2 )
          ! -- Update the running sum of Ps products
          thread_sum = thread_sum + Ps2_kernel(Sxyz1(:,a1,a2),Sxyz2(:,b1,b2),k,dla,dlb)
          ! -- Update the count of actually explored combination (thread_count .LE. N)
          thread_count = thread_count+1

        END DO
        !$OMP END DO

      !$OMP CRITICAL
        current_rsum   = current_rsum   + thread_sum
        current_rcount = current_rcount + thread_count
      !$OMP END CRITICAL

      !$OMP END PARALLEL

      current_ravg = current_rsum/current_rcount
      IF ( ABS(current_ravg - last_ravg) .LE. threshold) EXIT
      last_ravg = current_ravg

    END DO
    DEALLOCATE(lookup1,lookup2,indices)

    RS = current_ravg

  END FUNCTION R_S

! ----------

  REAL(dp) FUNCTION Ps2_kernel()

    COMPLEX(8) :: 

    Ps2_kernel = ( ABS(S1(1)*S2(1) + S1(2)*S2(2) + S1(3)*S2(3))**2. / (k*k + (dla + dlb)**2.) ) ! Not done, check all normalisations etc in the other kernel

  END FUNCTION Ps2_kernel

! ----------

END MODULE random_sampling
