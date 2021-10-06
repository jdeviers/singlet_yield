MODULE evalYield
  implicit none
! -- This module is meant for having python bindings.
! -- evalYield_diag and evalYield_offDiag_random accept TRANSPOSED Sxyzi, i.e of dim (di,di,3).
! -- Unlike evalYield_offDiag2p, evalYield_offDiag2p_serial and evalYield_full which accept (3,di,di).

  contains


! ----------

  DOUBLE PRECISION FUNCTION Ps2_kernel(S1,S2,k2,dla,dlb)

    COMPLEX(8)       :: S1(3),S2(3)
    DOUBLE PRECISION :: k2,dla,dlb

    Ps2_kernel = ( ABS(S1(1)*S2(1) + S1(2)*S2(2) + S1(3)*S2(3))**2.d0 / (k2 + (dla + dlb)**2.d0) )
  END FUNCTION Ps2_kernel

! ----------

  DOUBLE PRECISION FUNCTION evalYield_diag(k,Sxyz1T,Sxyz2T)
    USE OMP_LIB

!  .. Arguments ..
    COMPLEX(8),      INTENT(IN) :: Sxyz1T(:,:,:),Sxyz2T(:,:,:) ! dim: (di,di,3)
    DOUBLE PRECISION,INTENT(IN) :: k
!  .. Sampling vars ..
    INTEGER(8)       :: a1,b1
!  .. Params ..
    INTEGER(8)       :: d1,d2,Z
    DOUBLE PRECISION :: k2
!  .. Scalars ..
    DOUBLE PRECISION :: rsum

    d1 = UBOUND(Sxyz1T,1); d2 = UBOUND(Sxyz2T,1)
    Z  = FLOOR( (d1*d2)/4.d0 )
    k2 = k*k

    rsum = 0.d0

    DO a1 = 1,d1
      DO b1 = 1,d2
        rsum = rsum + Ps2_kernel( Sxyz1T(a1,a1,:),Sxyz2T(b1,b1,:),k2,0.d0,0.d0 )
      END DO
    END DO

    evalYield_diag = 0.25d0 + (rsum / Z)

  END FUNCTION evalYield_diag
  
! ----------

  DOUBLE PRECISION FUNCTION evalYield_full(k,Sxyz1,lambda1,Sxyz2,lambda2)

!  .. Parameters ..
    COMPLEX(8),      INTENT(IN) :: Sxyz1(:,:,:),Sxyz2(:,:,:)
    DOUBLE PRECISION,INTENT(IN) :: lambda1(:),lambda2(:)
    DOUBLE PRECISION,INTENT(IN) :: k
!  .. Local scalars ..
    INTEGER(8)                  :: a1,a2,b1,b2
    INTEGER(8)                  :: d1,d2,z
    DOUBLE PRECISION            :: v,k2,dl1,dl2
!  .. Local tensors ..
    COMPLEX(8),ALLOCATABLE      :: sA(:),sB(:)

    ALLOCATE( sA(UBOUND(Sxyz1,1)),sB(UBOUND(Sxyz2,1)) )

    d1 = UBOUND(Sxyz1,3)    ; d2 = UBOUND(Sxyz2,3) ! NOTE: 1-indexing
    z  = FLOOR( (d1*d2)/4. ); v  = 0.d0; k2 = k*k

    DO a1 = 1,d1
      DO a2 = 1,d1
        dl1 = lambda1(a1) - lambda1(a2)
        sA = Sxyz1(:,a1,a2)
        DO b1 = 1,d2
          DO b2 = 1,d2
            dl2 = lambda2(b1) - lambda2(b2)
            sB = Sxyz2(:,b1,b2)
            v = v + Ps2_kernel(sA,sB,k2,dl1,dl2)
          END DO
        END DO
      END DO
    END DO

    DEALLOCATE(sA,sB)
    v = v * k2 / z
    evalYield_full = 0.25d0 + v

  END FUNCTION evalYield_full

! ----------

  DOUBLE PRECISION FUNCTION evalYield_offDiag_random(k,Sxyz1T,lambda1,Sxyz2T,lambda2,nr_draws)
    USE OMP_LIB
!
! -- This standalone R_S corresponds to evalYield_offdiag_random in singlet_yield_noInter_eigenbasis.ipynb.
! -- Its returned value is scaled in the notebook, not in the return statement as in random_sampling.f90.
! -- It is normal that R_S makes no mention of N_max and is very small.
! -- Caution: this program will likely overshoot nr_draws, by a maximum of (d1*d2)-1 draws.
!

!  .. Arguments ..
    COMPLEX(8),      INTENT(IN) :: Sxyz1T(:,:,:),Sxyz2T(:,:,:) ! dim: (di,di,3)
    DOUBLE PRECISION,INTENT(IN) :: lambda1(:),lambda2(:)       ! dim: (di)
    DOUBLE PRECISION,INTENT(IN) :: k
    INTEGER(8),      INTENT(IN) :: nr_draws
!  .. Sampling vars ..
    INTEGER(8)       :: a
!  .. Params ..
    INTEGER(8)       :: d1,d2,Z
    DOUBLE PRECISION :: k2
!  .. OMP vars ..
    INTEGER(4)       :: thread_id
    INTEGER(8)       :: thread_count
    DOUBLE PRECISION :: thread_sum
!  .. Scalars .. 
    INTEGER(4)       :: a1,a2,b1,b2
    INTEGER(8)       :: current_rcount
    DOUBLE PRECISION :: dla,dlb
    DOUBLE PRECISION :: current_rsum
!  .. Local tensors ..
    DOUBLE PRECISION,ALLOCATABLE :: indices(:,:)


    d1 = UBOUND(Sxyz1T,1); d2 = UBOUND(Sxyz2T,1)
    Z = FLOOR( (d1*d2)/4.d0 )
    k2 = k*k
    current_rcount = 0; current_rsum = 0.d0

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

          a1 = INT(indices(a,1),KIND=4); a2 = INT(indices(a,2),KIND=4)
          b1 = INT(indices(a,3),KIND=4); b2 = INT(indices(a,4),KIND=4)
          dla = lambda1(a1) - lambda1(a2); dlb = lambda2(b1) - lambda2(b2)

          IF ( a1.EQ.a2 ) THEN
            IF ( b1.EQ.b2 ) CYCLE
          END IF

          thread_sum   = thread_sum + Ps2_kernel(Sxyz1T(a1,a2,:),Sxyz2T(b1,b2,:),k2,dla,dlb)
          thread_count = thread_count+1_8
        END DO
        !$OMP END DO

        !$OMP CRITICAL
          current_rsum   = current_rsum  + thread_sum
          current_rcount = current_rcount + thread_count
        !$OMP END CRITICAL

      !$OMP END PARALLEL

    END DO
    DEALLOCATE(indices) 

    ! Final normalisation and return:
    evalYield_offDiag_random = (current_rsum/current_rcount) * k2 / Z 

  END FUNCTION evalYield_offDiag_random

! ----------

  DOUBLE PRECISION FUNCTION evalYield_offdiag2p(d1,d2,k,Sxyz1,lambda1,Sxyz2,lambda2)
    USE OMP_LIB

!  .. Parameters ..
    INTEGER(8),      INTENT(IN) :: d1,d2
    COMPLEX(8),      INTENT(IN) :: Sxyz1(3,d1,d1),Sxyz2(3,d2,d2)
    DOUBLE PRECISION,INTENT(IN) :: lambda1(d1),lambda2(d2)
    DOUBLE PRECISION,INTENT(IN) :: k
!  .. Local scalars ..
    INTEGER(8)                  :: a1,z
    DOUBLE PRECISION            :: v,thread_v,k2


    z = FLOOR( (d1*d2)/4. ); k2 = k*k
    v = 0.d0

    !$OMP PARALLEL PRIVATE(thread_v) SHARED(v)
      thread_v = 0.d0

      !$OMP DO
      DO a1 = 1,d1
        thread_v = thread_v + evalYield_offdiag2p_kernel_F( d1,d2,k2,INT(a1,KIND=8),Sxyz1(:,:,a1),lambda1,Sxyz2,lambda2 )
      END DO
      !$OMP END DO

      !$OMP CRITICAL
        v = v + thread_v
      !$OMP END CRITICAL

    !$OMP END PARALLEL

    v = 2 * (v * k2 / z)
    evalYield_offdiag2p = v

  END FUNCTION evalYield_offdiag2p

! ----------

  DOUBLE PRECISION FUNCTION evalYield_offdiag2p_serial(k,Sxyz1,lambda1,Sxyz2,lambda2)

!  .. Parameters ..
    COMPLEX(8),      INTENT(IN) :: Sxyz1(:,:,:),Sxyz2(:,:,:)
    DOUBLE PRECISION,INTENT(IN) :: lambda1(:),lambda2(:)
    DOUBLE PRECISION,INTENT(IN) :: k
!  .. Local scalars ..
    INTEGER(8)                  :: a1,d1,d2,z
    DOUBLE PRECISION            :: v,k2

    d1 = UBOUND(Sxyz1,3); d2 = UBOUND(Sxyz2,3) ! NOTE: 1-indexing
    z = FLOOR( (d1*d2)/4. ); k2 = k*k; v = 0.d0

    DO a1 = 1,d1
      v = v + evalYield_offdiag2p_kernel_F( d1,d2,k2,INT(a1,KIND=8),Sxyz1(:,:,a1),lambda1,Sxyz2,lambda2 )
    END DO
    v = 2 * (v * k2 / z)
    evalYield_offdiag2p_serial = v

  END FUNCTION evalYield_offdiag2p_serial

! ----------

  DOUBLE PRECISION FUNCTION evalYield_offdiag2p_kernel_F(d1,d2,k2,a1,Sxyz1_a1,lambda1,Sxyz2,lambda2)

!  .. Parameters ..
    COMPLEX(8),      INTENT(IN) :: Sxyz2(3,d2,d2)
    COMPLEX(8),      INTENT(IN) :: Sxyz1_a1(3,d1)
    DOUBLE PRECISION,INTENT(IN) :: lambda1(d1),lambda2(d2)
    DOUBLE PRECISION,INTENT(IN) :: k2
    INTEGER(8),      INTENT(IN) :: a1,d1,d2
!  .. Local arrays ..
    COMPLEX(8)       :: Sxyz2_b1(3,d1)
!  .. Local scalars ..
    INTEGER(8)       :: a2,b1,b2
    DOUBLE PRECISION :: lambda1_a1,y,dl1,dl2
    COMPLEX(8)       :: sAx,sAy,sAz,sBx,sBy,sBz

    lambda1_a1 = lambda1(a1); y = 0.d0

    DO b1 = 1,d2
      Sxyz2_b1 = Sxyz2(:,:,b1)
      a2 = a1; b2 = b1
      sAx = Sxyz1_a1(1,a2); sAy = Sxyz1_a1(2,a2); sAz = Sxyz1_a1(3,a2)
      dl1 = lambda1_a1 - lambda1(a2)
      DO
        b2 = b2 + 1
        IF (b2 .EQ. d2+1) THEN
          b2 = 1; a2 = a2 + 1
          IF (a2 .EQ. d1+1) THEN
            EXIT ! exits the innermost DO loop
          END IF
          sAx = Sxyz1_a1(1,a2); sAy = Sxyz1_a1(2,a2); sAz = Sxyz1_a1(3,a2)
          dl1 = lambda1_a1 - lambda1(a2)
        END IF
        sBx = Sxyz2_b1(1,b2); sBy = Sxyz2_b1(2,b2); sBz = Sxyz2_b1(3,b2)
        dl2 = lambda2(b1) - lambda2(b2)

        y = y + ( ABS(sAx*sBx + sAy*sBy + sAz*sBz)**2. / (k2 + (dl1 + dl2)**2.) ) 
      END DO
    END DO

    evalYield_offdiag2p_kernel_F = y

  END FUNCTION evalYield_offdiag2p_kernel_F

! ----------

END MODULE evalYield
