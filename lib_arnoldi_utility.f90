MODULE arnoldi_utility

  IMPLICIT NONE

  INTEGER, PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  !Storage of L, R matrices if we use the common helper function
  COMPLEX(KIND=DP), ALLOCATABLE :: Arnoldi_MPO_L(:,:), Arnoldi_MPO_R(:,:) 
  COMPLEX(KIND=DP), ALLOCATABLE :: Arnoldi_MPO_temp(:)

  !Used in the alternative approach
  COMPLEX(KIND=DP), ALLOCATABLE :: Arnoldi_MPO_theta(:,:)

  !Machine precision
  REAL(KIND=DP) :: macheps=epsilon(1.0_DP)

  !Factor above machine precision to require on eigenvalues
  REAL(KIND=DP) :: tol_factor=20.0D0

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   app_LR_mpo_to_vec
  ! VARIABLES:
  !           vi - Vector in
  !           vo - Vector out
  !           hc - Whether we want HC of matrix
  ! SYNOPSIS:
  ! Standard wrapper for apply to operator, when what we have
  ! is an MPO of predefined L and R matrices
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE app_LR_mpo_to_vec(vi, vo, hc)

    COMPLEX (KIND=DP), INTENT(IN) :: vi(:)
    COMPLEX (KIND=DP), INTENT(OUT) :: vo(:)
    LOGICAL, INTENT(IN), OPTIONAL :: hc
    LOGICAL :: local_hc

    INTEGER :: ldim, rdim, mpodim
    COMPLEX (KIND=DP) :: alpha=CMPLX(1.0D0,0.0D0,DP), &
         &               beta =CMPLX(0.0D0, 0.0D0, DP)

    EXTERNAL ZGEMV

    ! Expected size from L, R matrices
    ldim=SIZE(Arnoldi_MPO_L,1); rdim=SIZE(Arnoldi_MPO_R,2)
    mpodim=SIZE(Arnoldi_MPO_L,2)
    CALL check_dim_match(SIZE(Arnoldi_MPO_R,1), mpodim, "LR consistency")
    CALL check_dim_match(SIZE(Arnoldi_MPO_temp), mpodim, "Temp vector")

    ! Default to not being hermitian conjugate
    local_hc=.FALSE.; IF (PRESENT(hc)) local_hc=hc

    IF (local_hc) THEN
       ! Hermitian conjugate of operator applied to vin
       CALL check_dim_match(SIZE(vi), ldim, "Input vec (HC)")
       CALL check_dim_match(SIZE(vo), rdim, "Output vec (HC)")

       ! Use BLAS to do hermitian conjugate multiplication,
       ! A'* x = R'* L'* x (prime is transpose)
       CALL zgemv('C', ldim, mpodim, alpha, &
            & Arnoldi_MPO_L, ldim, vi, 1, beta, Arnoldi_MPO_temp, 1)
       CALL zgemv('C', mpodim, rdim, alpha, &
            & Arnoldi_MPO_R, mpodim, Arnoldi_MPO_temp, 1, beta, vo, 1)

    ELSE
       ! Operator to be applied to vin.
       CALL check_dim_match(SIZE(vi), rdim, "Input vec")
       CALL check_dim_match(SIZE(vo), ldim, "Output vec")

       ! Project onto and out of vectors, using BLAS
       CALL zgemv('N', mpodim ,rdim,  alpha, &
            & Arnoldi_MPO_R, mpodim, vi, 1, beta, Arnoldi_MPO_temp, 1)
       CALL zgemv('N', ldim, mpodim, alpha, &
            & Arnoldi_MPO_L, ldim, Arnoldi_MPO_temp, 1, beta, vo, 1)
    end IF

  end SUBROUTINE app_LR_mpo_to_vec



  SUBROUTINE check_dim_match(v1,v2,text)
    INTEGER, INTENT(IN) :: v1, v2
    CHARACTER(LEN=*), INTENT(IN) :: text
    IF (v1 .NE. v2) THEN
       WRITE(*,*) "ARNOLDI: Mismatch of sizes: ", text, v1, v2
       STOP
    end IF
  end SUBROUTINE check_dim_match



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   update_nev(min, cur, max, direction)
  ! VARIABLES:
  !         la, lb - Range
  !         cur      - Value to be updated
  !         direction- +/-1 to choose which way to go
  ! SYNOPSIS:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE update_nev(la, cur, lb, direction, unchanged)

    INTEGER, INTENT(INOUT) :: cur
    INTEGER, INTENT(IN) :: la, lb, direction
    LOGICAL, INTENT(OUT) :: unchanged

    REAL (KIND=DP), PARAMETER :: frac=0.9
    INTEGER :: old

    ! Record old value
    old=cur

    ! Update toward appropriate limit.  Make sure we move by at least
    ! one unit whatever happens.
    IF (direction .GT. 0) THEN
       cur= INT(frac*cur + (1-frac)*lb)
       IF ((cur .EQ. old) .AND. (old .LT. lb)) cur=cur+1
    ELSE
       cur= INT(frac*cur + (1-frac)*la)
       IF ((cur .EQ. old) .AND. (old .GT. la)) cur=cur-1
    end IF

    ! Ensure sensible range
    cur=MAX(la,cur)
    cur=MIN(lb,cur)

    ! Write what we did
    !WRITE(*,'(" ARNOLDI: Adjusting nev from ",I4," to ",I4)')  old, cur

    ! Return whether it changed, to avoid infinite loops.
    unchanged=(cur .EQ. old)

  END SUBROUTINE update_nev



 !Alternative approach
 SUBROUTINE apply_theta_to_vec(vi, vo, hc)

    COMPLEX (KIND=DP), INTENT(IN) :: vi(:)
    COMPLEX (KIND=DP), INTENT(OUT) :: vo(:)
    LOGICAL, INTENT(IN), OPTIONAL :: hc
    LOGICAL :: local_hc

    INTEGER :: ldim, rdim
    COMPLEX (KIND=DP) :: alpha=CMPLX(1.0D0,0.0D0,DP), &
         &               beta=CMPLX(0.0D0, 0.0D0, DP)

    EXTERNAL ZGEMV

    ! Expected size from L, R matrices
    ldim=SIZE(Arnoldi_MPO_theta,1); rdim=SIZE(Arnoldi_MPO_theta,2)

    ! Default to not being hermitian conjugate
    local_hc=.FALSE.; IF (PRESENT(hc)) local_hc=hc

    IF(local_hc) THEN 
       ! Hermitian conjugate of operator applied to vin
       CALL check_dim_match(SIZE(vi), ldim, "Input vec (HC)")
       CALL check_dim_match(SIZE(vo), rdim, "Output vec (HC)")

       ! Use BLAS to do hermitian conjugate multiplication,
       ! A'* x = R'* L'* x (prime is transpose)
       CALL zgemv('C', ldim, rdim, alpha, &
            & Arnoldi_MPO_theta, ldim, vi, 1, beta, vo, 1)
    ELSE
       ! Operator to be applied to vin.
       CALL check_dim_match(SIZE(vi), rdim, "Input vec")
       CALL check_dim_match(SIZE(vo), ldim, "Output vec")

       ! Project onto and out of vectors, using BLAS
       CALL zgemv('N', ldim ,rdim,  alpha, &
            & Arnoldi_MPO_theta, ldim, vi, 1, beta, vo, 1)
    end IF

  end SUBROUTINE apply_theta_to_vec

END MODULE arnoldi_utility
