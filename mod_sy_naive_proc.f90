MODULE mod_sy_proc
	implicit none

	INTEGER(4),PARAMETER ::            &
		sp = SELECTED_REAL_KIND(6,37), & ! 32-bits precision
		dp = SELECTED_REAL_KIND(15,307)  ! 64-bits precision

	CHARACTER(LEN=*),PARAMETER ::                            &
		FMT_SCMPLX='(*( "(" SP F6.3 X SP F6.3 " i)" 2X ))',  &
		FMT_LCMPLX='(*( "(" SP E10.3 X SP E10.3 " i)" 2X ))'


	contains

! ----------
	REAL(dp) FUNCTION evalYield(k,Sxyz1,lambda1,Sxyz2,lambda2)
		implicit none

!	.. Parameters ..
		COMPLEX(8),ALLOCATABLE :: Sxyz1(:,:,:),Sxyz2(:,:,:)
		REAL(dp),  ALLOCATABLE :: lambda1(:),lambda2(:)
		REAL(dp)               :: k
!	.. Local scalars ..
		INTEGER                :: a1,a2,b1,b2
		INTEGER                :: d1,d2,z
		REAL(dp)               :: v,k2,dl1,dl2
!	.. Local tensors ..
		COMPLEX(8),ALLOCATABLE :: sA(:),sB(:)

		ALLOCATE( sA(UBOUND(Sxyz1,1)),sB(UBOUND(Sxyz2,1)) )

		d1 = UBOUND(Sxyz1,3); d2 = UBOUND(Sxyz2,3) ! NOTE: 1-indexing
		z = FLOOR( (d1*d2)/4. ); v=0.; k2 = k*k

		DO a1 = 1,d1
			DO a2 = 1,d1
				dl1 = lambda1(a1) - lambda1(a2)
				sA = Sxyz1(:,a1,a2)
				DO b1 = 1, d2
					DO b2 = 1,d2
						dl2 = lambda2(b1) - lambda2(b2)
						sB = Sxyz2(:,b1,b2)
						v = v + REAL(( ABS( DOT_PRODUCT(sA,sB) )**2. / (k2 + (dl1+dl2)**2.) )) ! IMAGPART = 0 because of the square
					END DO
				END DO
			END DO
		END DO

		DEALLOCATE(sA,sB)
		v = v * k2 / z
		evalYield = 0.25d0 + v

	END FUNCTION evalYield
! ----------
	REAL(dp) FUNCTION evalYield_offdiag2p(k,Sxyz1,lambda1,Sxyz2,lambda2)
		implicit none

!	.. Parameters ..
		COMPLEX(8),ALLOCATABLE :: Sxyz1(:,:,:),Sxyz2(:,:,:)
		REAL(dp),  ALLOCATABLE :: lambda1(:),lambda2(:)
		REAL(dp)               :: k
!	.. Local scalars ..
		INTEGER                :: a1,d1,d2,z
		REAL(dp)               :: v,k2


		d1 = UBOUND(Sxyz1,3); d2 = UBOUND(Sxyz2,3) ! NOTE: 1-indexing
		z = FLOOR( (d1*d2)/4. ); v=0.; k2 = k*k

		DO a1 = 1,d1
			v = v + evalYield_offdiag2p_kernel_F( k2,INT(a1),Sxyz1(:,:,a1),lambda1,Sxyz2,lambda2 )
		END DO
		v = 2 * (v * k2 / z)
		evalYield_offdiag2p = v

	END FUNCTION evalYield_offdiag2p
! ----------
	REAL(dp) FUNCTION evalYield_offdiag2p_kernel_F(k2,a1,Sxyz1_a1,lambda1,Sxyz2,lambda2)
		implicit none

!	.. Parameters ..
		COMPLEX(8),ALLOCATABLE :: Sxyz2(:,:,:)
		REAL(dp),  ALLOCATABLE :: lambda1(:),lambda2(:)
		COMPLEX(8)             :: Sxyz1_a1(:,:)
		REAL(dp)               :: k2
		INTEGER                :: a1	


	END FUNCTION evalYield_offdiag2p_kernel_F
! ----------
END MODULE mod_sy_proc