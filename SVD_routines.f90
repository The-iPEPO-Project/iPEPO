MODULE svd_routines

 USE basic_functions
 USE error_handling
 USE array_utility
 USE TENS_reshape
 USE TENS_transpose
 USE TENS_mult
 USE TENS_mult_extension

 USE matrix_functions, ONLY: DiagonalAsMat !! FIXME

 IMPLICIT NONE

 !Module level definitions of matrices/vectors needed for apply_theta_to_vector routine called inside arnoldi SVD
 COMPLEX(KIND=DP), ALLOCATABLE :: arnoldi_theta_op(:,:), arnoldi_theta_op_hc(:,:)

 interface compute_lapack_svd
    module procedure compute_lapack_svd_U_Theta
    module procedure compute_lapack_svd_U_VH_Sigma
 end interface compute_lapack_svd

 interface compute_arnoldi_svd
    module procedure compute_arnoldi_svd_U_Theta
    module procedure compute_arnoldi_svd_U_VH_Sigma
 end interface compute_arnoldi_svd

 interface truncate_svd_matrices
    module procedure truncate_svd_matrices_U_Sigma
    module procedure truncate_svd_matrices_U_VH_Sigma
 end interface truncate_svd_matrices

 interface filter_svd_matrices
    module procedure filter_svd_matrices_U_Sigma
    module procedure filter_svd_matrices_U_VH_Sigma
 end interface filter_svd_matrices

 private
 public :: compute_lapack_svd, compute_arnoldi_svd, setChi, is_lapack
 public :: filter_svd_matrices, rescale_svd_matrices

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPUTE LAPACK SVD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute lapack SVD = (U, Theta) --- (input chi=-1 means 'do not truncate') !!!
 SUBROUTINE compute_lapack_svd_U_Theta(theta, U, chi, eps, chiMax, errTrunc)

  USE lapack_wrappers

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(IN)              :: theta(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT)             :: U(:,:)
  INTEGER,                       INTENT(INOUT)           :: chi
  REAL(KIND=DP),                 INTENT(IN)              :: eps
  INTEGER,                       INTENT(INOUT), OPTIONAL :: chiMax
  REAL(KIND=DP),                 INTENT(INOUT), OPTIONAL :: errTrunc

  !! Temp vars needed for svd
  COMPLEX(KIND=DP), ALLOCATABLE :: VH(:,:), theta_tmp(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)

  !! Dims & indices
  INTEGER :: dimT(2), sigma_dim
  INTEGER :: i

  !! (0A) Get dims & allocate untruncated SVD matrices
  dimT = shape(theta); sigma_dim = MIN(dimT(1), dimT(2))
  ALLOCATE(U(dimT(1), dimT(1)), VH(dimT(2), dimT(2)), Sigma(sigma_dim))

  !! (0B) Create a copy of theta (since theta input is modified by lapack_gen_svd)
  ALLOCATE(theta_tmp(dimT(1), dimT(2))); theta_tmp(:,:) = theta(:,:)

  !! (1) Lapack svd
  CALL lapack_gen_svd(theta_tmp, Sigma, U, VH)
  DEALLOCATE(theta_tmp, VH)

  !! (2) Adjust value of chi for truncation (important in fixed-chi mode)
  chi = setChi(chi, sigma_dim)

  !! --- Filter out small evals
  CALL filter_svd_matrices(U, Sigma, 1.0D-08)

  !! (3) Truncate svd matrices 
  CALL truncate_svd_matrices(U, Sigma, chi, eps, errTrunc)

  !! (4) Impose upper limit on bond dimension
  IF(PRESENT(chiMax)) THEN
     IF((chi .GT. chiMax) .AND. (chiMax .GT. 0)) CALL truncate_svd_matrices(U, Sigma, chiMax, 2.0d0, errTrunc)
  end IF

 END SUBROUTINE compute_lapack_svd_U_Theta





 !!! Compute lapack SVD = (U, VH, Sigma) --- (input chi=-1 means 'do not truncate') !!!
 SUBROUTINE compute_lapack_svd_U_VH_Sigma(thetaIn, U, VH, Sigma, chi, eps, chiMax, errTrunc)

  USE lapack_wrappers

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(IN)              :: thetaIn(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT)             :: U(:,:), VH(:,:)
  REAL(KIND=DP),    ALLOCATABLE, INTENT(OUT)             :: Sigma(:)
  INTEGER,                       INTENT(INOUT)           :: chi
  REAL(KIND=DP),                 INTENT(IN)              :: eps
  INTEGER,                       INTENT(INOUT), OPTIONAL :: chiMax
  REAL(KIND=DP),                 INTENT(INOUT), OPTIONAL :: errTrunc

  !! Local copy of thetaIn
  COMPLEX(KIND=DP), ALLOCATABLE :: theta(:,:)

  !! Dims & indices
  INTEGER :: dimT(2), sigma_dim
  INTEGER :: i, j

  !! (0A) Create a local copy of thetaIn (since theta will be modified by lapack_gen_svd)
  ALLOCATE(theta(SIZE(thetaIn,1), SIZE(thetaIn,2))); theta(:,:) = thetaIn(:,:)

  !! (0B) Set dims && allocate untruncated svd matrices
  dimT = shape(theta); sigma_dim = MIN(dimT(1), dimT(2))
  ALLOCATE(U(dimT(1), dimT(1)), VH(dimT(2), dimT(2)), Sigma(sigma_dim))

  !! (1) Compute lapack svd
  CALL lapack_gen_svd(theta, Sigma, U, VH)
  !WRITE(*,*) "Maxvals: ", MAXVAL(ABS(theta)), MAXVAL(ABS(U)), MAXVAL(ABS(VH)), MAXVAL(ABS(Sigma))
  DEALLOCATE(theta)

  !! (2) Adjust value of chi for truncation
  !!     -- chi = dim of fixed-chi truncation
  !!     -- chi = lower-bound of fixed-eps truncation
  chi = setChi(chi, sigma_dim)

  !! (3) Truncate svd matrices (if eps < 1.0d0, use fixed-eps mode -- else, use fixed-chi mode)
  CALL truncate_svd_matrices(U, VH, Sigma, chi, eps, errTrunc)

  !! (4) Impose upper limit on bond dimension (if upper dim chiMaxP is present)
  IF(PRESENT(chiMax)) THEN
     IF((chi .GT. chiMax) .AND. (chiMax .GT. 0)) CALL truncate_svd_matrices(U, VH, Sigma, chiMax, 2.0d0, errTrunc)
  end IF

  !! (5) Rescale svd matrices so that SUM(Sigma**2) = 1
  CALL rescale_svd_matrices(U, VH, Sigma)

  !WRITE(*,*) "Sigmas: "
  !CALL printVector(complx(Sigma))

 END SUBROUTINE compute_lapack_svd_U_VH_Sigma

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPUTE ARNOLDI SVD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute Arnoldi SVD = (U, Theta) -- (using Arnoldi algorithm & by constructing a linear operator) !!!
 SUBROUTINE compute_arnoldi_svd_U_Theta(theta, U, chi, eps, accuracy_not_achieved)

  USE arnoldi_wrappers

  COMPLEX(KIND=DP),  ALLOCATABLE, INTENT(IN)    :: theta(:,:)
  COMPLEX(KIND=DP),  ALLOCATABLE, INTENT(OUT)   :: U(:,:)
  INTEGER,                        INTENT(INOUT) :: chi
  REAL(KIND=DP),                  INTENT(IN)    :: eps
  LOGICAL,                        INTENT(INOUT) :: accuracy_not_achieved

  COMPLEX(KIND=DP), ALLOCATABLE :: VH(:,:)   
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)

  !dims & indices
  INTEGER :: dimT(2), sigma_dim
  INTEGER :: delta_chi

  !Determine dims of SVD matrices
  dimT = shape(theta); sigma_dim = MIN(dimT(1), dimT(2))

  !Set the value of chi
  chi=NINT(0.11*sigma_dim)

  !Allocate SVD matrices
  ALLOCATE(U(dimT(1), chi), VH(chi, dimT(2)), Sigma(chi))

  !Allocate other objects needed for apply_op_to_vec routine called in arnoldi SVD
  ALLOCATE(arnoldi_theta_op(dimT(1), dimT(2)), arnoldi_theta_op_hc(dimT(2), dimT(1)))

  !Arnoldi SVD (& module-level copy of theta, needed for apply_op_to_vec)
  arnoldi_theta_op = theta
  arnoldi_theta_op_hc = CONJG(TRANSPOSE(theta))

  CALL mpo_arnoldi_svd(apply_arnoldi_op_to_vec, U, Sigma, VH)

  !clear memory
  DEALLOCATE(VH, arnoldi_theta_op, arnoldi_theta_op_hc)

  !Have we achieved the required accuracy in this SVD?
  accuracy_not_achieved = .NOT. is_optimal_eps(Sigma, eps)

  !Increase chi (acc not achieved), or proceed to trunc SVD matrices (acc achieved)
  IF(accuracy_not_achieved) THEN

      !determine the increment of chi
      delta_chi = MAX(2, FLOOR(ABS(0.11*sigma_dim - chi)))

      !increase chi
      chi = chi + delta_chi 
      WRITE(*,*) "Arnoldi SVD: target accuracy not achieved, increasing chi to", chi
  ELSE

      !truncate & return SVD matrices
      WRITE(*,*) "Arnoldi SVD: target accuracy achieved"
      CALL truncate_svd_matrices(U, Sigma, chi, eps)
  end IF

 END SUBROUTINE compute_arnoldi_svd_U_Theta




 !!! Compute Arnoldi SVD = (U, VH, Sigma) -- (using Arnoldi algorithm & linear operator) !!!
 SUBROUTINE compute_arnoldi_svd_U_VH_Sigma(theta, U, VH, Sigma, chi, eps, accuracy_not_achieved)

  USE arnoldi_wrappers

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(IN)    :: theta(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT)   :: U(:,:), VH(:,:)
  REAL(KIND=DP),    ALLOCATABLE, INTENT(OUT)   :: Sigma(:)
  INTEGER,                       INTENT(INOUT) :: chi
  REAL(KIND=DP),                 INTENT(IN)    :: eps
  LOGICAL,                       INTENT(INOUT) :: accuracy_not_achieved

  !dims & indices
  INTEGER :: dimT(2), sigma_dim
  INTEGER :: delta_chi

  !Determine dims of SVD matrices
  dimT = shape(theta); sigma_dim = MIN(dimT(1), dimT(2))

  !Set the value of chi
  chi=NINT(0.11*sigma_dim)

  !Allocate SVD matrices
  ALLOCATE(U(dimT(1), chi), VH(chi, dimT(2)), Sigma(chi))

  !Allocate other objects needed for apply_op_to_vec routine called in arnoldi SVD
  ALLOCATE(arnoldi_theta_op(dimT(1), dimT(2)), arnoldi_theta_op_hc(dimT(2), dimT(1)))

  !Arnoldi SVD (& module-level copy of theta, needed for apply_op_to_vec)
  arnoldi_theta_op = theta
  arnoldi_theta_op_hc = CONJG(TRANSPOSE(theta))

  CALL mpo_arnoldi_svd(apply_arnoldi_op_to_vec, U, Sigma, VH)

  !clear memory
  DEALLOCATE(arnoldi_theta_op, arnoldi_theta_op_hc)

  !Have we achieved the required accuracy in this SVD?
  accuracy_not_achieved = .NOT. is_optimal_eps(Sigma, eps)

  !Increase chi (acc not achieved), or proceed to trunc SVD matrices (acc achieved)
  IF(accuracy_not_achieved) THEN

      !determine the increment of chi
      delta_chi = MAX(2, FLOOR(ABS(0.11*sigma_dim - chi)))

      !increase chi
      chi = chi + delta_chi 
      WRITE(*,*) "Arnoldi SVD: target accuracy not achieved, increasing chi to", chi
  ELSE

      !truncate & return SVD matrices
      WRITE(*,*) "Arnoldi SVD: target accuracy achieved"
      CALL truncate_svd_matrices(U, VH, Sigma, chi, eps)
  end IF

 END SUBROUTINE compute_arnoldi_svd_U_VH_Sigma



 !!! Apply Arnoldi Operator to Vector (vi = input vec, vo = output vec, hc = TRUE if we apply HC of theta matrix to vi) !!!
 SUBROUTINE apply_arnoldi_op_to_vec(vi, vo, hc)

  COMPLEX(KIND=DP), INTENT(IN)           :: vi(:)
  COMPLEX(KIND=DP), INTENT(OUT)          :: vo(:)
  LOGICAL,          INTENT(IN), OPTIONAL :: hc

  LOGICAL :: local_hc

  !Default to not being hermitian conjugate
  local_hc=.FALSE.; IF(PRESENT(hc)) local_hc=hc

  !Initialize the output to zero
  vo=0.0D0

  !apply op to input vec ---> get output vec
  IF(local_hc) THEN
     !apply HC theta op to |vi>
     vo = MATMUL(arnoldi_theta_op_hc, vi)
  ELSE
     !apply theta op to |vi>
     vo = MATMUL(arnoldi_theta_op, vi)
  end IF

 END SUBROUTINE apply_arnoldi_op_to_vec

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TRUNCATE SVD OUTPUT MATRICES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Truncate SVD matrices (U, Sigma) -- NB. chi on exit = number of Sigmas we keep !!!
 SUBROUTINE truncate_svd_matrices_U_Sigma(U, Sigma, chi, eps, errTrunc)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT)           :: U(:,:)
  REAL(KIND=DP),    ALLOCATABLE, INTENT(INOUT)           :: Sigma(:)
  INTEGER,                       INTENT(INOUT)           :: chi
  REAL(KIND=DP),                 INTENT(IN)              :: eps
  REAL(KIND=DP),                 INTENT(INOUT), OPTIONAL :: errTrunc

  !! Temp copies of svd matrix U
  COMPLEX(KIND=DP), ALLOCATABLE :: U_tmp(:,:)

  !! Dims & indices 
  INTEGER :: dimU(2), dimS
  INTEGER :: i

  !! Find dims
  dimU = shape(U); dimS = SIZE(Sigma)

  !! Create temp copy of pre-trunc SVD matrix U
  ALLOCATE(U_tmp(dimU(1), dimU(2))); U_tmp = U
  
  !! DETERMINE truncation chi
  CALL find_svd_truncation_dim(chi, eps, Sigma)

  !! Calc SVD truncation error
  IF(PRESENT(errTrunc)) CALL calc_svd_trunc_error(errTrunc, Sigma, chi)

  !! Allocate new post-trunc SVD matrix U
  DEALLOCATE(U); ALLOCATE(U(dimU(1), chi))

  !! Perform truncation using chi obtained above
  U(:,:) = U_tmp(:, 1:chi)

  !! Clear memory
  DEALLOCATE(Sigma, U_tmp)

 END SUBROUTINE truncate_svd_matrices_U_Sigma







 !!! Truncate SVD matrices (U, VH, Sigma) -- NB. chi on exit = number of Sigmas we keep !!!
 SUBROUTINE truncate_svd_matrices_U_VH_Sigma(U, VH, Sigma, chi, eps, errTrunc)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT)           :: U(:,:), VH(:,:)
  REAL(KIND=DP),    ALLOCATABLE, INTENT(INOUT)           :: Sigma(:)
  INTEGER,                       INTENT(INOUT)           :: chi
  REAL(KIND=DP),                 INTENT(IN)              :: eps
  REAL(KIND=DP),                 INTENT(INOUT), OPTIONAL :: errTrunc

  !! Temp copies of svd matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: U_tmp(:,:), VH_tmp(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma_tmp(:) 

  !! Dims & indices 
  INTEGER :: dimU(2), dimV(2), dimS
  INTEGER :: i

  !! Find dims
  dimU = shape(U); dimV = shape(VH); dimS = SIZE(Sigma)

  !! Create temp copies of pre-trunc SVD matrices
  ALLOCATE(U_tmp(dimU(1), dimU(2)), VH_tmp(dimV(1), dimV(2)), Sigma_tmp(dimS))
  U_tmp     = U
  VH_tmp    = VH
  Sigma_tmp = Sigma
  
  !! DETERMINE truncation chi
  CALL find_svd_truncation_dim(chi, eps, Sigma)

  !! Calc SVD truncation error
  IF(PRESENT(errTrunc)) CALL calc_svd_trunc_error(errTrunc, Sigma, chi)

  !! Allocate new post-trunc SVD matrices
  DEALLOCATE(U,VH,Sigma); ALLOCATE(U(dimU(1), chi), VH(chi, dimV(2)), Sigma(chi))

  !! Perform truncation using chi obtained above
  U(:,:)   = U_tmp(:, 1:chi)
  VH(:,:)  = VH_tmp(1:chi, :)
  Sigma(:) = Sigma_tmp(1:chi)

  !! Clear memory
  DEALLOCATE(U_tmp, VH_tmp, Sigma_tmp)

 END SUBROUTINE truncate_svd_matrices_U_VH_Sigma





 !!! Find SVD truncation dim given Sigma = sval matrix, eps = desired precision, chi = lower bound on bond dim !!!
 SUBROUTINE find_svd_truncation_dim(chi, eps, Sigma)

  INTEGER,                       INTENT(INOUT)        :: chi
  REAL(KIND=DP),                 INTENT(IN)           :: eps
  REAL(KIND=DP),    ALLOCATABLE, INTENT(IN)           :: Sigma(:)

  !! Indices && dims
  INTEGER :: dimS, i

  !! Find dims of Sigma
  dimS = SIZE(Sigma)

  !! If (eps<1.0) --> find what chi is needed to achieved desired precision
  !! Else         --> chi remains unchanged
  IF(eps .LT. 1.0) THEN

      !! Run precision loop
      DO i=1,dimS

         !! Find what chi is needed to achieve the required trunc err
         IF(is_optimal_eps(Sigma(1:i), eps)) THEN

             chi = MAX(i-1, chi)
             EXIT

         ELSEIF(i .EQ. dimS) THEN

             chi = dimS
         end IF
      END DO
  end IF

 END SUBROUTINE find_svd_truncation_dim




 !!! Calculate SVD truncation error !!!
 SUBROUTINE calc_svd_trunc_error(errTrunc, Sigma, chi)

  REAL(KIND=DP), INTENT(INOUT) :: errTrunc
  REAL(KIND=DP), INTENT(IN)    :: Sigma(:)
  INTEGER,       INTENT(IN)    :: chi
 
  INTEGER :: chiTOT

  !! Size of Sigma before truncation
  chiTOT = SIZE(Sigma)
  
  !! Calculate truncation error
  errTrunc = SUM(Sigma(chi+1 : chiTOT)**2) / SUM(Sigma(1 : chiTOT)**2)
  
 END SUBROUTINE calc_svd_trunc_error





 !!! Rescale SVD matrices U, VH, Sigma so that SUM(Sigma**2) = 1
 SUBROUTINE rescale_svd_matrices(U, VH, Sigma)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: U(:,:), VH(:,:)
  REAL(KIND=DP),    ALLOCATABLE, INTENT(INOUT) :: Sigma(:)

  !! Rescaling factor 
  REAL(KIND=DP) :: fac

  !! Calculate fac = sqrt(SUM(Sigma**2))
  fac = sqrt(SUM(Sigma**2)) 

  !! Rescale SVD matrices
  Sigma = Sigma/fac
  U  = U  * sqrt(fac)
  VH = VH * sqrt(fac)

 END SUBROUTINE rescale_svd_matrices





 !!! Filter out small evecs && evals (U, VH, Sigma) !!!
 SUBROUTINE filter_svd_matrices_U_VH_Sigma(U, VH, Sigma, eps)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT)  :: U(:,:), VH(:,:)
  REAL(KIND=DP),    ALLOCATABLE, INTENT(INOUT)  :: Sigma(:)
  REAL(KIND=DP),                 INTENT(IN)     :: eps

  INTEGER       :: i
  REAL(KIND=DP) :: fac

  !! Calc SUM(Sigma**2) for norm
  fac = SUM(Sigma**2)

  !! Nullify evecs whose sigmas are below threshold
  DO i=1,SIZE(Sigma)
     IF(Sigma(i)/fac .LT. eps) THEN
         U(:,i) = (0.0D0, 0.0D0)
        VH(i,:) = (0.0D0, 0.0D0)
       Sigma(i) =  0.0D0
     end IF
  end DO
  
 END SUBROUTINE filter_svd_matrices_U_VH_Sigma




 !!! Filter out small evecs && evals (U, Sigma) !!! 
 SUBROUTINE filter_svd_matrices_U_Sigma(U, Sigma, eps)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT)  :: U(:,:)
  REAL(KIND=DP),    ALLOCATABLE, INTENT(INOUT)  :: Sigma(:)
  REAL(KIND=DP),                 INTENT(IN)     :: eps

  INTEGER       :: i
  REAL(KIND=DP) :: fac

  !! Calc SUM(Sigma**2) for norm
  fac = SUM(Sigma**2)

  !! Nullify evecs whose sigmas are below threshold
  DO i=1,SIZE(Sigma)
     IF(Sigma(i)/fac .LT. eps) THEN
        U(:,i)   = (0.0D0, 0.0D0)
        Sigma(i) =  0.0D0
     end IF
  end DO

 END SUBROUTINE filter_svd_matrices_U_Sigma

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!! SETTING SVD PARAMS (Chi, Eps, Lapack vs Arnoldi, etc) !!!!!!!!!!!!!!!!!!!!!!!!!!!

 !set chi_input=-1 to get chi=sigma_dim, chi_input=0 to keep mps dim unchanged
 !input sigma_dim = MIN(dimT(1), dimT(2), chi_current)
 FUNCTION setChi(chi_in, sigma_dim) result(chi)

  INTEGER, INTENT(IN) :: chi_in, sigma_dim
  INTEGER :: chi !output

  !find the value of chi (truncation param), sigma_dim = max possible dim     
  IF(chi_in .GT. 0) THEN
     chi = MIN(chi_in, sigma_dim)
  ELSE
     chi = sigma_dim
  end IF

 END FUNCTION setChi



 !eps_optimized, optimal_eps
 FUNCTION is_optimal_eps(Sigma, eps0) 

  REAL(KIND=DP), INTENT(IN) :: Sigma(:)
  REAL(KIND=DP), INTENT(IN) :: eps0
  LOGICAL :: is_optimal_eps

  REAL(KIND=DP) :: eps

  eps = MINVAL(Sigma(:))/MAXVAL(Sigma(:))
  is_optimal_eps = (eps .LE. eps0)

 end FUNCTION is_optimal_eps



 !make decision whether to use lapack or arnoldi at a given site of mps/mpo network
 FUNCTION is_lapack(chi, sigma_dim, Edim_is_one)

  INTEGER, INTENT(IN) :: chi, sigma_dim
  LOGICAL, INTENT(IN) :: Edim_is_one
  LOGICAL :: is_lapack

  !Use lapack if the following conditions hold:
  is_lapack = ((sigma_dim .LT. 1000) .OR. (chi .GT. 0.11*sigma_dim) .OR. Edim_is_one)

 END FUNCTION is_lapack

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE svd_routines
