MODULE evalYield_randomAccess
  implicit none
! -- This module is meant for having python bindings.
! -- It accepts TRANSPOSED Sxyzi, i.e of dim (di,di,3).
! -- Unlike the random_sampling program, which has (3,di,di).


! -- Defines compiler-independent data types.


  contains

  DOUBLE PRECISION FUNCTION offdiag_random(k,Sxyz1T,lambda1,Sxyz2T,lambda2,nr_draws)
!
! -- This standalone R_S corresponds to evalYield_offdiag_random in singlet_yield_noInter_eigenbasis.ipynb.
! -- Its returned value is scaled in the notebook, not in the return statement as in random_sampling.f90.
! -- It is normal that R_S makes no mention of N_max and is very small.
! -- Caution: this program will likely overshoot nr_draws, by a maximum of (d1*d2)-1 draws.
!
  	USE OMP_LIB

!	.. Arguments ..
  	COMPLEX(8),		 INTENT(IN) :: Sxyz1T(:,:,:),Sxyz2T(:,:,:) ! dim: (di,di,3)
  	DOUBLE PRECISION,INTENT(IN) :: lambda1(:),lambda2(:)       ! dim: (di)
  	DOUBLE PRECISION,INTENT(IN) :: k
  	INTEGER(8),		 INTENT(IN) :: nr_draws
!	.. Sampling vars ..
  	INTEGER(8)		 :: a
!	.. Params ..
  	INTEGER(8)		 :: d1,d2,Z
  	DOUBLE PRECISION :: k2
!	.. OMP vars ..
  	INTEGER(4)		 :: thread_id
  	INTEGER(8)		 :: thread_count
  	DOUBLE PRECISION :: thread_sum
!	.. Scalars .. 
  	INTEGER(4)		 :: a1,a2,b1,b2
  	INTEGER(8)		 :: current_rcount
  	DOUBLE PRECISION :: dla,dlb
  	DOUBLE PRECISION :: current_rsum
!	.. Local tensors ..
  	DOUBLE PRECISION,ALLOCATABLE :: indices(:,:)


  	d1 = UBOUND(Sxyz1T,1); d2 = UBOUND(Sxyz2T,1)
  	Z = FLOOR( (d1*d2)/4.d0 )
  	k2 = k*k

  	current_rcount = 0; current_rsum = 0.d0
  	thread_count = 0; thread_sum = 0.d0

  	ALLOCATE( indices(d1*d2,4) )

  	DO WHILE (current_rcount .LT. nr_draws)! This setup will most likely overshoot nr_draws
 
  	  CALL RANDOM_NUMBER(indices)
  	  indices = indices*(d1-1_4)+1_4 ! Prevents an index to be =0
				
  	  !$OMP PARALLEL                                    &
	  !$OMP PRIVATE(thread_id,thread_count,thread_sum)  &
	  !$OMP SHARED(indices,current_rcount,current_rsum)
  		
	  thread_sum = 0.d0; thread_count = 0.d0 ! Reset the threadwise sums and counts				

	  !$OMP DO
  	  DO a = 1,UBOUND(indices,1)

	  	a1 = INT(indices(a,1)); a2 = INT(indices(a,2))
	  	b1 = INT(indices(a,3)); b2 = INT(indices(a,4))
	  	dla = lambda1(a1) - lambda1(a2); dlb = lambda2(b1) - lambda2(b2)

		IF ( a1.EQ.a2 ) THEN
		  IF ( b1.EQ.b2 ) CYCLE
		END IF

  		thread_sum   = thread_sum + Ps2_kernel(Sxyz1T(a1,a2,:),Sxyz2T(b1,b2,:),k2,dla,dlb)
  		thread_count = thread_count+1_8
  	  END DO
  	  !$OMP END DO

  	  !$OMP CRITICAL
  		current_rsum   = current_rsum	+ thread_sum
  		current_rcount = current_rcount + thread_count
	  !$OMP END CRITICAL

  	  !$OMP END PARALLEL

  	END DO
  	DEALLOCATE(indices) 

  	! Final normalisation and return:
  	offdiag_random = (current_rsum/current_rcount) * k2 / Z 

  END FUNCTION offdiag_random

! ----------

  DOUBLE PRECISION FUNCTION diag_ordered(k,Sxyz1T,Sxyz2T)
	USE OMP_LIB

!	.. Arguments ..
  	COMPLEX(8),		 INTENT(IN) :: Sxyz1T(:,:,:),Sxyz2T(:,:,:) ! dim: (di,di,3)
  	DOUBLE PRECISION,INTENT(IN) :: k
!	.. Sampling vars ..
	INTEGER(8)       :: a1,b1
!	.. Params ..
  	INTEGER(8)		 :: d1,d2,Z
  	DOUBLE PRECISION :: k2
!	.. Scalars ..
	DOUBLE PRECISION :: rsum

  	d1 = UBOUND(Sxyz1T,1); d2 = UBOUND(Sxyz2T,1)
  	Z = FLOOR( (d1*d2)/4.d0 )
  	k2 = k*k

	rsum = 0.d0

	DO a1 = 1,d1
	  DO b1 = 1,d2
        rsum = rsum + Ps2_kernel( Sxyz1T(a1,a1,:),Sxyz2T(b1,b1,:),k2,0.d0,0.d0 )
	  END DO
    END DO

    diag_ordered = 0.25d0 + (rsum / Z)

  END FUNCTION diag_ordered
! ----------

  DOUBLE PRECISION FUNCTION Ps2_kernel(S1,S2,k2,dla,dlb)

  	COMPLEX(8)		 :: S1(3),S2(3)
  	DOUBLE PRECISION :: k2,dla,dlb

  	Ps2_kernel = ( ABS(S1(1)*S2(1) + S1(2)*S2(2) + S1(3)*S2(3))**2.d0 / (k2 + (dla + dlb)**2.d0) )
  END FUNCTION Ps2_kernel

! ----------

END MODULE evalYield_randomAccess