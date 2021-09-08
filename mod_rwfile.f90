MODULE mod_rwfile
    implicit none

	INTEGER(4),PARAMETER ::            &
		sp = SELECTED_REAL_KIND(6,37), & ! 32-bits precision
		dp = SELECTED_REAL_KIND(15,307)  ! 64-bits precision

	CHARACTER(LEN=*),PARAMETER ::                            &
		FMT_SCMPLX='(*( "(" SP F6.3 X SP F6.3 " i)" 2X ))',  &
		FMT_LCMPLX='(*( "(" SP E10.3 X SP E10.3 " i)" 2X ))'

    contains
!
! ----------
!
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
!
! ----------
!
    SUBROUTINE READ_3xNx2N_FROM_F(N,ARG_ID,Sxyz)
        implicit none

!   .. Arguments ..
        COMPLEX(8),ALLOCATABLE,INTENT(OUT) :: Sxyz(:,:,:)
        INTEGER,               INTENT(IN)  :: N,ARG_ID
!
!   .. Local tensors ..
        REAL(dp),ALLOCATABLE            :: tmpSxyz(:,:,:)
!
!   .. Local scalars ..
        INTEGER                         :: i,j,k,io
        CHARACTER(LEN=40)               :: infile

		ALLOCATE(tmpSxyz(3,N,2*N))
		CALL GETARG(ARG_ID,infile)
		OPEN(10,file=infile,status='OLD',action='READ')
		DO i=1,3
			DO j=1,N
				READ(10,*,iostat=io) tmpSxyz(i,j,:)
			END DO
		END DO
		CLOSE(10)
        ALLOCATE(Sxyz(3,N,N))
        DO i=1,3
            DO j=1,N
		        DO k=1,N
			        Sxyz(i,j,k) = COMPLEX(tmpSxyz(i,j,INT(2*k-1)),tmpSxyz(i,j,INT(2*k)))
                END DO
            END DO
		END DO
		DEALLOCATE(tmpSxyz)

    END SUBROUTINE READ_3xNx2N_FROM_F
!
! ----------
!
    SUBROUTINE READ_N_FROM_F(N,ARG_ID,L)
        implicit none

!   .. Arguments ..
        REAL(dp),ALLOCATABLE,INTENT(OUT) :: L(:)
        INTEGER,             INTENT(IN)  :: N,ARG_ID
!
!   .. Local scalars ..
        INTEGER                         :: i,io
        CHARACTER(LEN=40)               :: infile

        ALLOCATE(L(N))
		CALL GETARG(ARG_ID,infile)
		OPEN(10,file=infile,status='OLD',action='READ')
		DO i=1,N
			READ(10,*,iostat=io) L(i)
		END DO
		CLOSE(10)

    END SUBROUTINE READ_N_FROM_F
!
! ----------
!
    SUBROUTINE WRITE_3xNxN_TO_F(N,outfile,Sxyz)
        implicit none

!   .. Arguments ..
        COMPLEX(8),ALLOCATABLE,INTENT(IN) :: Sxyz(:,:,:)
        CHARACTER(LEN=*),      INTENT(IN)  :: outfile
        INTEGER,               INTENT(IN)  :: N
!
!   .. Local scalars ..
        INTEGER                            :: i,j
        
        OPEN(10,file=outfile,status='UNKNOWN',action='WRITE')
		DO i=1,3
			DO j=1,N
				WRITE(10,'(*(F8.4))') Sxyz(i,j,:)
			END DO
		END DO
		CLOSE(10)

    END SUBROUTINE WRITE_3xNxN_TO_F
!
! ----------
!
    SUBROUTINE WRITE_N_TO_F(N,outfile,L)
        implicit none

!   .. Arguments ..
        REAL(dp),ALLOCATABLE,INTENT(IN) :: L(:)
        CHARACTER(LEN=*),    INTENT(IN)  :: outfile
        INTEGER,             INTENT(IN)  :: N
!
!   .. Local scalars ..
        INTEGER                          :: i

		OPEN(10,file=outfile,status='UNKNOWN',action='WRITE')
		DO i=1,N
			WRITE(10,'(F8.4)') L(i)
		END DO
		CLOSE(10)

    END SUBROUTINE WRITE_N_TO_F
!
! ----------
!
END MODULE mod_rwfile