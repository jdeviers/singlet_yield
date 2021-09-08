DOUBLE PRECISION FUNCTION evalYield_offdiag2p(d1,d2,k,Sxyz1,lambda1,Sxyz2,lambda2)
  USE OMP_LIB
  implicit none

!  .. Parameters ..
  INTEGER,          INTENT(IN) :: d1,d2
  COMPLEX(8),       INTENT(IN) :: Sxyz1(3,d1,d1),Sxyz2(3,d2,d2)
  DOUBLE PRECISION,         INTENT(IN) :: lambda1(d1),lambda2(d2)
  DOUBLE PRECISION,         INTENT(IN) :: k
!  .. Local scalars ..
  INTEGER(8)                     :: a1,z
  DOUBLE PRECISION                     :: v,thread_v,k2


  z = FLOOR( (d1*d2)/4. ); k2 = k*k

  v        = 0.d0
  !$OMP PARALLEL PRIVATE(thread_v) SHARED(v)
    thread_v = 0.d0

    !$OMP DO
    DO a1 = 1,d1
      thread_v = thread_v + evalYield_offdiag2p_kernel_F( d1,d2,k2,INT(a1),Sxyz1(:,:,a1),lambda1,Sxyz2,lambda2 )
    END DO
    !$OMP END DO

    !$OMP CRITICAL
      v = v + thread_v
    !$OMP END CRITICAL

  !$OMP END PARALLEL
  v = 2 * (v * k2 / z)
  evalYield_offdiag2p = v

  contains

  DOUBLE PRECISION FUNCTION evalYield_offdiag2p_kernel_F(d1,d2,k2,a1,Sxyz1_a1,lambda1,Sxyz2,lambda2)
    implicit none

!  .. Parameters ..
    COMPLEX(8),INTENT(IN) :: Sxyz2(3,d2,d2)
    COMPLEX(8),INTENT(IN) :: Sxyz1_a1(3,d1)
    DOUBLE PRECISION,  INTENT(IN) :: lambda1(d1),lambda2(d2)
    DOUBLE PRECISION,  INTENT(IN) :: k2
    INTEGER,   INTENT(IN) :: a1,d1,d2
!  .. Local arrays ..
    COMPLEX(8)           :: Sxyz2_b1(3,d1)
!  .. Local scalars ..
    INTEGER(8)              :: a2,b1,b2
    DOUBLE PRECISION              :: lambda1_a1,y,dl1,dl2
    COMPLEX(8)            :: sAx,sAy,sAz,sBx,sBy,sBz

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

END FUNCTION evalYield_offdiag2p

