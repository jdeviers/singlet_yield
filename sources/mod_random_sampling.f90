MODULE mod_random_sampling ! Can also be called from prog_sy_parallel
  USE mod_rwfile
  implicit none

  contains

! ----------

  REAL(dp) FUNCTION R_S(k,Sxyz1,lambda1,Sxyz2,lambda2,nr_draws) 
    USE OMP_LIB

!  .. Arguments ..
    COMPLEX(8),ALLOCATABLE,INTENT(IN) :: Sxyz1(:,:,:),Sxyz2(:,:,:)
    REAL(dp),  ALLOCATABLE,INTENT(IN) :: lambda1(:),lambda2(:)
    REAL(dp)              ,INTENT(IN) :: k
    INTEGER(4)            ,INTENT(IN) :: nr_draws
!  .. Sampling vars ..
    INTEGER(8)       :: a
!  .. Params ..
    INTEGER(8)       :: d1,d2,Z,N_max
    DOUBLE PRECISION :: k2
!  .. OMP vars ..
    INTEGER(4)       :: thread_id
    INTEGER(8)       :: thread_count
    DOUBLE PRECISION :: thread_sum
!  .. Scalars .. 
    INTEGER(8)       :: a1,a2,b1,b2
    INTEGER(8)       :: current_rcount
    DOUBLE PRECISION :: dla,dlb
    DOUBLE PRECISION :: current_rsum
!  .. Local tensors ..
    DOUBLE PRECISION,ALLOCATABLE :: indices(:,:)


    d1 = UBOUND(Sxyz1,3) ; d2 = UBOUND(Sxyz2,3)
    Z  = FLOOR( (d1*d2)/4.d0 )
    k2 = k*k

    N_max = INT( ((d1**2.d0 * d2**2.d0) - d1*d2),KIND=8 )
    current_rcount = 0; current_rsum = 0.d0


    ALLOCATE( indices(d1*d2,4) )

    DO WHILE (current_rcount .LT. nr_draws)
      CALL RANDOM_NUMBER(indices)
      indices = indices*(d1-1_4)+1_4 ! Prevents an index to be =0
  

      !$OMP PARALLEL                                    &
      !$OMP PRIVATE(thread_id,thread_count,thread_sum)  &
      !$OMP SHARED(indices,current_rcount,current_rsum)

        thread_count = 0; thread_sum = 0.d0 ! Reset the threadwise sums and counts
        
        !$OMP DO
        DO a = 1,UBOUND(indices,1)

          a1 = INT(indices(a,1),KIND=8); a2 = INT(indices(a,2),KIND=8)
          b1 = INT(indices(a,3),KIND=8); b2 = INT(indices(a,4),KIND=8)
          dla = lambda1(a1) - lambda1(a2); dlb = lambda2(b1) - lambda2(b2)

          ! Skip diagonal elements
          IF ( a1.EQ.a2 ) THEN
              IF ( b1.EQ.b2 ) CYCLE ! Nested IF avoid performing too many checks at every step
          END IF

          thread_sum = thread_sum + Ps2_kernel(Sxyz1(:,a1,a2),Sxyz2(:,b1,b2),k2,dla,dlb)
          thread_count = thread_count+1_8

        END DO
        !$OMP END DO

        !$OMP CRITICAL
          current_rsum   = current_rsum   + thread_sum
          current_rcount = current_rcount + thread_count
          WRITE(10,*) ( (current_rsum/current_rcount) * N_max * k2 / Z ), current_rsum, current_rcount
        !$OMP END CRITICAL

      !$OMP END PARALLEL

    END DO
    DEALLOCATE(indices) 

    ! Final normalisation and return:
    R_S = (current_rsum/current_rcount) * N_max * k2 / Z 

  END FUNCTION R_S

! ----------

  REAL(dp) FUNCTION Ps2_kernel(S1,S2,k2,dla,dlb)

    COMPLEX(8) :: S1(3),S2(3)
    REAL(dp)   :: k2,dla,dlb

    Ps2_kernel = ( ABS(S1(1)*S2(1) + S1(2)*S2(2) + S1(3)*S2(3))**2.d0 / (k2 + (dla + dlb)**2.d0) )

  END FUNCTION Ps2_kernel

! ----------

END MODULE mod_random_sampling
