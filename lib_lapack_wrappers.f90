MODULE lapack_wrappers

  IMPLICIT NONE

  INTEGER, PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  ! This module contains wrappers around lapack routines, to cope with
  ! workspace etc, and implicit sizes, as a replacement for the NAG
  ! wrappers.
  ! Note that these are defined as recursive in order to ensure thread
  ! safe variables --- they obviously are never called recursively, but
  ! they are called in omp parallel sections.
  
CONTAINS
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   check_error
  ! VARIABLES:
  !           condition
  !           errormsg
  ! SYNOPSIS:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE check_error(prefix, info, errormsg)
    INTEGER, INTENT(IN) :: info
    CHARACTER(LEN=*), INTENT(IN) :: errormsg, prefix

    IF (info .NE. 0) THEN
       WRITE(*,*) prefix
       WRITE(*,*) errormsg
       WRITE(*,*) "Error code", info
       STOP
    end IF

  END SUBROUTINE check_error




  SUBROUTINE assert(prefix, condition, errormsg)
    LOGICAL, INTENT(IN) ::condition
    CHARACTER(LEN=*), INTENT(IN) :: errormsg,prefix

    IF (.NOT. condition) THEN
       WRITE(*,*) prefix
       WRITE(*,*) errormsg
       STOP
    end IF
    
  END SUBROUTINE assert


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   lapack_nsym_eig_all
  ! VARIABLES:
  !           matrix  - To get evals
  !           evals   - Evals returned
  ! SYNOPSIS:
  ! Wrapper for Lapack
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE lapack_nsym_eig_all(matrix,evals)
    COMPLEX (KIND=DP), INTENT(INOUT) :: matrix(:,:)
    COMPLEX (KIND=DP), INTENT(OUT) :: evals(:)

    COMPLEX (KIND=DP), ALLOCATABLE :: WORK(:)
    REAL (KIND=DP), ALLOCATABLE :: RWORK(:)
    COMPLEX (KIND=DP) :: dummy(1)
    INTEGER :: dim, lwork, info
    CHARACTER :: prefix*80

    EXTERNAL ZGEEV

    prefix="lapack_nsym_eigv_all"

    dim=SIZE(matrix,1)
    CALL assert(prefix,SIZE(matrix,2) .EQ. dim, "Non square matrix")
    CALL assert(prefix,SIZE(evals) .EQ. dim, "Wrong size evals")

    lwork=2*dim
    ALLOCATE(work(lwork), rwork(2*dim))

    CALL ZGEEV('N', 'N', dim, matrix, dim, evals, &
         & dummy, 1, dummy, 1, &
         & work, lwork, rwork, info)

    CALL check_error(prefix,info, "ZGEEV Failed")
    
    DEALLOCATE(work, rwork)

    
  end SUBROUTINE lapack_nsym_eig_all
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   lapack_nsym_eigv_all
  ! VARIABLES:
  !           matrix
  !           evals
  !           evecsl
  !           evecsr
  ! SYNOPSIS:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE lapack_nsym_eigv_all(matrix,evals,evecsl,evecsr)
    COMPLEX (KIND=DP), INTENT(INOUT) :: matrix(:,:)
    COMPLEX (KIND=DP), INTENT(OUT) :: evals(:),evecsl(:,:),evecsr(:,:)

    COMPLEX (KIND=DP), ALLOCATABLE :: WORK(:)
    REAL (KIND=DP), ALLOCATABLE :: RWORK(:)
    INTEGER :: dim, lwork, info
    CHARACTER :: prefix*80

    EXTERNAL ZGEEV

    prefix="lapack_nsym_eigv_all"

    dim=SIZE(matrix,1)
    CALL assert(prefix,SIZE(matrix,2) .EQ. dim, "Non square matrix")
    CALL assert(prefix,SIZE(evals) .EQ. dim, "Wrong size evals")
    CALL assert(prefix,(SIZE(evecsl,1) .EQ. dim) .AND.&
         &      (SIZE(evecsl,2) .EQ. dim), "Wrong size evecs matrix")
    CALL assert(prefix,(SIZE(evecsr,1) .EQ. dim) .AND.&
         &      (SIZE(evecsr,2) .EQ. dim), "Wrong size evecs matrix")

    lwork=2*dim
    ALLOCATE(work(lwork), rwork(2*dim))

    CALL ZGEEV('V', 'V', dim, matrix, dim, evals, &
         & evecsl, dim, evecsr, dim, &
         & work, lwork, rwork, info)

    CALL check_error(prefix,info, "ZGEEV Failed")
    
    DEALLOCATE(work, rwork)
    
  END SUBROUTINE lapack_nsym_eigv_all


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   lapack_sym_eigv_all
  ! VARIABLES:
  !           matrix
  !           evals
  !           evecs
  ! SYNOPSIS:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE lapack_sym_eigv_all(matrix,evals,evecs)
    COMPLEX (KIND=DP), INTENT(INOUT) :: matrix(:,:)
    REAL (KIND=DP), INTENT(OUT) :: evals(:)
    COMPLEX (KIND=DP), INTENT(OUT) :: evecs(:,:)

    COMPLEX (KIND=DP), ALLOCATABLE :: WORK(:)
    REAL (KIND=DP), ALLOCATABLE :: RWORK(:)
    INTEGER :: dim, lwork, info
    CHARACTER :: prefix*80

    EXTERNAL ZHEEV

    prefix="lapack_nsym_eigv_all"

    dim=SIZE(matrix,1)
    CALL assert(prefix,SIZE(matrix,2) .EQ. dim, "Non square matrix")
    CALL assert(prefix,SIZE(evals) .EQ. dim, "Wrong size evals")
    CALL assert(prefix,(SIZE(evecs,1) .EQ. dim) .AND.&
         &      (SIZE(evecs,2) .EQ. dim), "Wrong size evecs matrix")

    lwork=5*dim
    ALLOCATE(work(lwork), rwork(3*dim))

    CALL ZHEEV('V', 'U', dim, matrix, dim, evals, &
         & work, lwork, rwork, info)
    evecs=matrix

    CALL check_error(prefix,info, "ZHEEV Failed")
    
    DEALLOCATE(work, rwork)
    
  END SUBROUTINE lapack_sym_eigv_all


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   lapack_gen_svd
  ! VARIABLES:
  !           a
  !           sigma
  !           u
  !           vh
  ! SYNOPSIS:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE lapack_gen_svd(a,sigma,u,vh)
    COMPLEX (KIND=DP), INTENT(INOUT) ::a(:,:)
    REAL (KIND=DP), INTENT(OUT) :: sigma(:)
    COMPLEX (KIND=DP), INTENT(OUT) ::u(:,:),vh(:,:)
    
    COMPLEX (KIND=DP), ALLOCATABLE :: WORK(:)
    REAL (KIND=DP), ALLOCATABLE :: RWORK(:)
    INTEGER (KIND=DP), ALLOCATABLE :: iwork(:)
    INTEGER :: diml, dimr, mindim, lwork, lrwork, info
    CHARACTER :: prefix*80

    EXTERNAL ZGESVD

    prefix="lapack_gen_svd"


    diml=SIZE(a,1); dimr=SIZE(a,2)
    mindim=MIN(diml, dimr)
    CALL assert(prefix,SIZE(sigma,1) .EQ. mindim, "Sigma wrong size")
    CALL assert(prefix,(SIZE(u,1) .EQ. diml) .AND. &
         &      (SIZE(u,2) .EQ. diml), "U wrong size")
    CALL assert(prefix,(SIZE(vh,1) .EQ. dimr) .AND. &
         &      (SIZE(vh,2) .EQ. dimr), "VH wrong size")	!this was previously diml, causing an error
    
    lwork = mindim*mindim + 2*mindim+MAX(diml,dimr)
    lrwork = mindim*MAX(5*mindim+7,2*MAX(diml,dimr)+2*mindim+1)
    
    ALLOCATE(WORK(lwork), rwork(lrwork), iwork(8*mindim))

    CALL ZGESVD('A','A', diml, dimr,  a, diml,  sigma, &
         & u, diml, vh, dimr,  &
         & work, lwork, rwork, info)

    CALL check_error(prefix,info, "ZGESVD Failed")

    DEALLOCATE(work, rwork, iwork)

    
  END SUBROUTINE lapack_gen_svd


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   lapack_gen_mat_inv
  ! VARIABLES:
  !           matrix
  ! SYNOPSIS:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE lapack_gen_mat_inv(matrix)
    COMPLEX (KIND=DP), INTENT(INOUT) ::matrix(:,:)

    INTEGER, ALLOCATABLE :: IPIV(:)
    COMPLEX (KIND=DP), ALLOCATABLE :: work(:)
    INTEGER :: dim, info, lwork
    CHARACTER :: prefix*80

    prefix="lapack_gen_mat_inv"

    dim=SIZE(matrix,1)
    CALL assert(prefix,SIZE(matrix,2) .EQ. dim, "Non square matrix")

    lwork = 5*dim
    ALLOCATE(IPIV(dim), work(lwork))
   
    CALL ZGETRF(dim, dim, matrix, dim, ipiv, info)
    CALL check_error(prefix,info, "ZGETRF Failed")

    CALL ZGETRI(dim, matrix, dim, ipiv, work, lwork, info)
    CALL check_error(prefix,info, "ZGETRI Failed")

    DEALLOCATE(ipiv,work)

  END SUBROUTINE lapack_gen_mat_inv




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   lapack_gen_linsol
  ! VARIABLES:
  !           X
  !           A
  !           B
  ! SYNOPSIS:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RECURSIVE SUBROUTINE lapack_gen_linsol(X, A, B)

    COMPLEX(KIND=DP), INTENT(INOUT) :: X(:,:)
    COMPLEX(KIND=DP), INTENT(IN)    :: A(:,:), B(:,:)
    
    !! Local arrays
    COMPLEX(KIND=DP), ALLOCATABLE :: Amat(:,:), Bmat(:,:)
    INTEGER,          ALLOCATABLE :: iPIV(:)    

    !! Local dims && indices
    INTEGER   :: N, NRHS, ldimA, ldimB, info
    CHARACTER :: prefix*80

    EXTERNAL ZGESV
    prefix="lapack_gen_linsol"

    !! Set dims
    ldimA = SIZE(A, 1);  ldimB = SIZE(B, 1)  
    N     = SIZE(A, 2);  NRHS  = SIZE(B, 2)

    !! Check X dims
    CALL assert(prefix, (SIZE(X,1) .EQ. N) .AND. (SIZE(X,2) .EQ. NRHS), "X wrong size")

    !! Allocate && initialize arrays for ZGESV
    ALLOCATE(Amat(ldimA, N), Bmat(ldimB, NRHS), iPIV(N))
    Amat = A
    Bmat = B
    
    !! Call general linear solver ZGESV
    CALL ZGESV(N, NRHS, Amat, ldimA, iPIV, Bmat, ldimB, info)

    !! Check for errors
    CALL check_error(prefix, info, "ZGESV Failed")

    !! Output Bmat contains the solution matrix X
    X = Bmat

  END SUBROUTINE lapack_gen_linsol


END MODULE lapack_wrappers
