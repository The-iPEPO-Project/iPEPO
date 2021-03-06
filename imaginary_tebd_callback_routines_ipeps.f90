MODULE imaginary_tebd_callback_routines_ipeps

 USE utility
 USE simulation_parameters
 USE imaginary_observables_ipeps

 IMPLICIT NONE

 interface create_initial_state_vec
    module procedure create_psi_vec_from_input_vec
    module procedure create_rho_vec_from_input_vec
    module procedure create_rho_vec_from_rho_mat
 end interface create_initial_state_vec

 INTEGER :: DUMP_NUM, G2_SUMMARY

 private create_psi_vec_from_input_vec, create_rho_vec_from_input_vec, create_rho_vec_from_rho_mat

CONTAINS

 !!! Test if both 1P && 2P observables have converged !!!
 SUBROUTINE test_tebd_convergence(is_converged, new_obs_2P, old_obs_2P, new_obs_1P, old_obs_1P, & 
                                                & rho_onesite, rho_twosite, TEBD, iter, time, use_last_env)

  LOGICAL,                       INTENT(INOUT) :: is_converged                        !Whether observables have converged
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: new_obs_2P(:,:), old_obs_2P(:,:)    !New && old values of observables (to be compared)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: new_obs_1P(:),   old_obs_1P(:)      !New && old values of observables (to be compared)
  COMPLEX(KIND=DP),              INTENT(IN)    :: rho_onesite(:,:)                    !One-site reduced density matrix
  COMPLEX(KIND=DP),              INTENT(IN)    :: rho_twosite(:,:,:,:,:)              !Two-site reduced density matrix
  TYPE(tebd_params),             INTENT(IN)    :: TEBD                                !TEBD params
  INTEGER,                       INTENT(IN)    :: iter                                !Time evolution -- iter we are at
  REAL(KIND=DP),                 INTENT(IN)    :: time                                !Time evolution -- time we are at (corresponds to iter)
  LOGICAL,                       INTENT(IN)    :: use_last_env                        !Whether we've used previous environment as initialization 

  LOGICAL :: is_converged_1P, is_converged_2P

  !! Initialize
  is_converged = .FALSE.

  !! Check 1P observables
  CALL test_tebd_conv_1P(is_converged_1P, new_obs_1P, old_obs_1P, & 
                          & rho_onesite,               tebd, iter, time, use_last_env)

  !! Check 2P observables
  !! CALL test_tebd_conv_2P(is_converged_2P, new_obs_2P, old_obs_2P, &
                          !& rho_onesite, rho_twosite, tebd, iter, time, use_last_env)

  !! Both 1P && 2P must be converged
  is_converged = is_converged_1P !.AND. is_converged_2P

 END SUBROUTINE test_tebd_convergence





 !!! Test if 1P observables have converged !!!
 SUBROUTINE test_tebd_conv_1P(is_converged, new_obs, old_obs, rho_onesite, TEBD, iter, time, use_last_env)

  USE datafile_utility

  LOGICAL,                       INTENT(INOUT) :: is_converged            !whether observables have converged
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: new_obs(:), old_obs(:)  !new && old values of observables (to be compared)
  COMPLEX(KIND=DP),              INTENT(IN)    :: rho_onesite(:,:)        !one-site reduced density matrix for computing 1P observables
  TYPE(tebd_params),             INTENT(IN)    :: TEBD                    !TEBD params
  INTEGER,                       INTENT(IN)    :: iter                    !Time evolution -- iter we are on
  REAL(KIND=DP),                 INTENT(IN)    :: time                    !Time evolution -- time we are at (corresponds to iter)
  LOGICAL,                       INTENT(IN)    :: use_last_env            !Whether we've used previous environment as initialization 

  INTEGER              :: Nop, i           !operator indexing
  LOGICAL, ALLOCATABLE :: obs_converged(:) !convergence tracking for individual obs
  REAL(KIND=DP)        :: etaR             !rescaled eta
  COMPLEX(KIND=DP)     :: norm             !Norm of TN

  !number of operators
  Nop = SIZE(TRACE_OP_1P, 1) 

  !Convergence tracking for diff obs
  ALLOCATE(obs_converged(Nop)); obs_converged(:) = .FALSE.

  !Allocate storage for transient observables
  IF(.NOT. ALLOCATED(new_obs)) ALLOCATE(new_obs(Nop))
  IF(.NOT. ALLOCATED(old_obs)) ALLOCATE(old_obs(Nop))

  !Rescaled convergence cutoff
  etaR = (tebd % eta) * (tebd % test_interval) * ABS(tebd % dt)

  !Get norm of TN
  norm = TRACE(rho_onesite) 

  !Compute all observables on the list i={1,nop}
  DO i=1,Nop
     CALL compute_obs_1P(new_obs(i), rho_onesite, TRACE_OP_1P(i,:,:))
     new_obs(i) = new_obs(i)/norm
     CALL test_convergence(obs_converged(i), new_obs(i), old_obs(i), & 
                               & etaR, TRACE_OUT_1P(i), iter, tebd % test_interval, time, norm, use_last_env)  !1, time, norm, use_last_env) 
  end DO

  !Check convergence
  is_converged = isTrueVec(obs_converged, 'AND') 

  !Set old --> old = new for next iter
  old_obs = new_obs 

 END SUBROUTINE test_tebd_conv_1P




 !!! Test if 2P observables have converged !!!
 SUBROUTINE test_tebd_conv_2P(is_converged, new_obs, old_obs, rho_onesite, rho_twosite, TEBD, iter, time, use_last_env)

  USE datafile_utility

  LOGICAL,                       INTENT(INOUT) :: is_converged                 !Whether observables have converged
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: new_obs(:,:), old_obs(:,:)   !New && old values of observables (to be compared)
  COMPLEX(KIND=DP),              INTENT(IN)    :: rho_onesite(:,:)             !One-site reduced density matrix 
  COMPLEX(KIND=DP),              INTENT(IN)    :: rho_twosite(:,:,:,:,:)       !Two-site reduced density matrix 
  TYPE(tebd_params),             INTENT(IN)    :: TEBD                         !TEBD params
  INTEGER,                       INTENT(IN)    :: iter                         !Time evolution -- iter we are at
  REAL(KIND=DP),                 INTENT(IN)    :: time                         !Time evolution -- time we are at (corresponds to iter)
  LOGICAL,                       INTENT(IN)    :: use_last_env                 !Whether we've used previous environment as initialization 

  INTEGER              :: Nop, i                       !operator indexing
  INTEGER              :: min_sep, max_sep, Nsep, sep  !separations
  LOGICAL, ALLOCATABLE :: obs_converged(:,:)           !convergence tracking for individual obs
  REAL(KIND=DP)        :: etaR                         !rescaled eta
  COMPLEX(KIND=DP)     :: norm                         !Norm of TN

  !number of operators
  Nop = SIZE(TRACE_OP_2P, 1)

  !separations
  min_sep = tebd % min_sep
  max_sep = tebd % max_sep
  Nsep = max_sep - min_sep + 1

  !Convergence tracking for diff obs
  ALLOCATE(obs_converged(Nop, Nsep)); obs_converged(:,:) = .FALSE.

  !Allocate storage for transient observables
  IF(.NOT. ALLOCATED(new_obs)) ALLOCATE(new_obs(Nop, Nsep))
  IF(.NOT. ALLOCATED(old_obs)) ALLOCATE(old_obs(Nop, Nsep))

  !Rescaled convergence cutoff
  etaR = (tebd % eta) * (tebd % test_interval) * ABS(tebd % dt)

  !Get norm of TN
  norm = TRACE(rho_onesite) 

  !Compute all 2P observables on the list i={1,nop}
  DO i=1,Nop

     !compute obs
     CALL compute_obs_2P(new_obs(i,:), rho_twosite, rho_onesite, TRACE_OP_2P(i,:,:,:), min_sep, max_sep)

     !normalize obs, and test convergence
     DO sep = min_sep, max_sep
        new_obs(i,sep+1-min_sep) = new_obs(i,sep+1-min_sep)/norm**(sep/2 + 1)
        CALL test_convergence(obs_converged(i,sep), new_obs(i,sep), old_obs(i,sep), & 
                                  & etaR, TRACE_OUT_2P(i,sep), iter, tebd % test_interval, time, norm, use_last_env)
     end DO
  end DO

  !Check convergence
  is_converged = isTrueMat(obs_converged, 'AND')

  !Set old --> old = new for next iter
  old_obs = new_obs

 END SUBROUTINE test_tebd_conv_2P

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Lambda Convergence !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Test if Lambdas have converged !!!
 SUBROUTINE test_lambda_convergence(is_converged, newLambda, oldLambda, tebd, datafile, iter, time)

  LOGICAL,                         INTENT(INOUT) :: is_converged            !whether observables have converged
  TYPE(block_lambda),              INTENT(IN)    :: newLambda(:)            !new Lambdas
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: oldLambda(:)            !old Lambdas
  TYPE(tebd_params),               INTENT(IN)    :: tebd                    !TEBD params
  INTEGER,                         INTENT(INOUT) :: datafile                !datafile for convergence data
  INTEGER,                         INTENT(IN)    :: iter                    !Time evolution -- iter we are on
  REAL(KIND=DP),                   INTENT(IN)    :: time                    !Time evolution -- time we are at (corresponds to iter)

  INTEGER              :: N_sites, site       !! lambda sites
  LOGICAL, ALLOCATABLE :: lambda_converged(:) !! convergence tracking for individual lambdas
  REAL(KIND=DP)        :: etaR                !! Rescaled eta

  !! Number of lambdas
  N_sites = SIZE(newLambda) 

  !! Convergence tracking for diff lambdas
  ALLOCATE(lambda_converged(N_sites)); lambda_converged(:) = .FALSE.  

  !! Rescaled convergence cutoff
  etaR = (tebd % eta) * (tebd % test_interval) * ABS(tebd % dt)

  !! Test Lambda convergence
  DO site=1,N_sites
     CALL test_convergence_vec(lambda_converged(site), newLambda(site) % m, oldLambda(site) % m, etaR, &
                                   & datafile, iter, TEBD % test_interval, xval=TEBD % dt, PRINT_ON=.TRUE.)
  end DO

  !! Check for convergence
  is_converged = isTrueVec(lambda_converged, 'AND')

  !! Set old --> old = new for next iter
  CALL copy_lambda_block(oldLambda, newLambda)

 END SUBROUTINE test_lambda_convergence





 !!! Get lambda MIN/MAX ratio to estimate trunc error !!!
 !!! (especially when the actual SVD info is not available, e.g. after reloading a state) !!!
 FUNCTION lambda_min_max(iLambda) result(eps)

  TYPE(block_lambda), INTENT(IN) :: iLambda(:)
  REAL(KIND=DP)                  :: eps 

  !! Trunc error on each bond
  REAL(KIND=DP), ALLOCATABLE :: eps_lambda(:)

  !! Indices && dims
  INTEGER :: bond, N_bonds

  !! Get num of bonds/lambdas
  N_bonds = SIZE(iLambda)

  !! (1) Find eps on each bond
  ALLOCATE(eps_lambda(N_bonds))

  DO bond=1,N_bonds
     eps_lambda(bond) = MINVAL(ABS(iLambda(bond) % m)) / MAXVAL(ABS(iLambda(bond) % m))
  end DO

  !! (2) Take the worst trunc prec out of all lambdas 
  !!     -- i.e. find the dominant truncation error
  eps = MAXVAL(eps_lambda)

 END FUNCTION lambda_min_max

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Decide which param to use in the param-loop, based on "xname" !!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE read_xname(x, xname, v1, v1_str, v2, v2_str, v3, v3_str, v4, v4_str)
  
  REAL(KIND=DP),    INTENT(IN)              :: x
  CHARACTER(LEN=*), INTENT(IN)              :: xname

  REAL(KIND=DP),    INTENT(INOUT)           :: v1
  CHARACTER(LEN=*), INTENT(IN)              :: v1_str

  REAL(KIND=DP),    INTENT(INOUT), OPTIONAL :: v2
  CHARACTER(LEN=*), INTENT(IN),    OPTIONAL :: v2_str

  REAL(KIND=DP),    INTENT(INOUT), OPTIONAL :: v3
  CHARACTER(LEN=*), INTENT(IN),    OPTIONAL :: v3_str

  REAL(KIND=DP),    INTENT(INOUT), OPTIONAL :: v4
  CHARACTER(LEN=*), INTENT(IN),    OPTIONAL :: v4_str

  !Local copies of optional string vars
  CHARACTER(LEN=16) :: v2_strX, v3_strX, v4_strX

  !Set default values (If var not present -- set its string to INVALID_STRING, to guarantee it's not selected)
  v2_strX = OptArg(v2_str, "INVALID_STRING")
  v3_strX = OptArg(v3_str, "INVALID_STRING")
  v4_strX = OptArg(v4_str, "INVALID_STRING")
  
  !Decide which param to use
  IF(xname .EQ. v1_str) THEN

           WRITE(*,*) "Setting ", v1_str, " = ", x
           v1 = x

  ELSEIF(xname .EQ. v2_strX) THEN

           WRITE(*,*) "Setting ", v2_str, " = ", x
           CALL setOptVar(v2, x, "set_xvar: argument must be present")

  ELSEIF(xname .EQ. v3_strX) THEN

           WRITE(*,*) "Setting ", v3_str, " = ", x
           CALL setOptVar(v3, x, "set_xvar: argument must be present")

  ELSEIF(xname .EQ. v4_strX) THEN

           WRITE(*,*) "Setting ", v4_str, " = ", x
           CALL setOptVar(v4, x, "set_xvar: argument must be present")

  ELSE
           CALL invalid_flag("set_xvar: invalid xname ", xname)
  end IF

 END SUBROUTINE read_xname

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Create initial state vector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !Create initial state vector directly by reading in the input vector
 SUBROUTINE create_psi_vec_from_input_vec(psi_vec, s0, psi_flag) 

  USE simulation_parameters, ONLY: local_dim

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: psi_vec(:)
  REAL(KIND=DP),                 INTENT(INOUT) :: s0(:)
  CHARACTER(LEN=*),              INTENT(IN)    :: psi_flag

  !Consistency checks
  IF(local_dim .LT. 2)     CALL invalid_value("create_psi_vec_from_input_vec: invalid local dim ",   local_dim)
  CALL check_sizes_equal(SIZE(s0), local_dim, "create_psi_vec_from_input_vec: s0 must have a size of local_dim ")

  ALLOCATE(psi_vec(local_dim))

  !Set psi_vec = sqrt of probability amplitude coefficients
  s0 = s0/sqrt(sum(s0*s0))
  psi_vec = sqrt(s0)

 END SUBROUTINE create_psi_vec_from_input_vec





 !Create initial rho vector directly by reading in the input vector
 SUBROUTINE create_rho_vec_from_input_vec(rho_vec, s0) 

  USE simulation_parameters, ONLY: local_dim

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: rho_vec(:)
  REAL(KIND=DP),                 INTENT(INOUT) :: s0(:)

  !Consistency checks
  IF(local_dim .LT. 2) CALL invalid_value("create_rho_vec_from_input_vec: invalid local dim ", local_dim)
  CALL check_sizes_equal(SIZE(s0), local_dim - 1, "create_rho_vec_from_input_vec: s0 must have a size of local_dim - 1 ")

  ALLOCATE(rho_vec(local_dim))

  !Read rho_vec (NB. to ensure positivity, we must make sure s0 is not too outside the Bloch sphere)
  IF(SUM(s0*s0) .GT. 1.0D0) s0 = s0/sqrt(sum(s0*s0))
  rho_vec(1) = 1.0D0
  rho_vec(2:local_dim) = s0

 END SUBROUTINE create_rho_vec_from_input_vec






 !Create initial rho vector from input density matrix
 SUBROUTINE create_rho_vec_from_rho_mat(rho_vec, rho_mat, represent) 

  USE simulation_parameters, ONLY: local_dim

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: rho_vec(:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: rho_mat(:,:)

  INTERFACE 
    FUNCTION represent(rho_mat)
      USE utility
      USE simulation_parameters, ONLY: local_dim
      COMPLEX(KIND=DP), INTENT(IN) :: rho_mat(INT(sqrt(local_dim*1.0d0)), INT(sqrt(local_dim*1.0d0)))
      COMPLEX(KIND=DP)             :: represent(local_dim)
    end FUNCTION represent
  end INTERFACE

  COMPLEX(KIND=DP) :: norm

  !Consistency checks
  IF(local_dim .LT. 2) CALL invalid_value("create_rho_vec_from_rho_mat: invalid local dim ", local_dim)

  !Rho normalization (trace)
  norm = TRACE(rho_mat)

  !Represent rho matrix as a state vector using Gell-Mann basis matrices
  ALLOCATE(rho_vec(local_dim))
  rho_vec = represent(rho_mat / norm)

 END SUBROUTINE create_rho_vec_from_rho_mat

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 !Test if 1P observables have converged
 SUBROUTINE simplified_test_tebd_conv(is_converged, new_obs, old_obs, ipeps, Cmat, Tmat, TEBD, CTMRG, iter, time, use_last_env)

  USE datafile_utility

  LOGICAL,                              INTENT(INOUT) :: is_converged            !whether observables have converged
  COMPLEX(KIND=DP),        ALLOCATABLE, INTENT(INOUT) :: new_obs(:), old_obs(:)  !new && old values of observables (to be compared)
  TYPE(block_peps),                     INTENT(IN)    :: ipeps(:)                !iPEPS 
  TYPE(ctm_corner_type),                INTENT(IN)    :: Cmat(:)                 !CTM environment -- Cmat
  TYPE(ctm_transfer_type),              INTENT(IN)    :: Tmat(:)                 !CTM environment -- Tmat
  TYPE(tebd_params),                    INTENT(IN)    :: TEBD                    !TEBD params
  TYPE(ctm_params),                     INTENT(IN)    :: CTMRG                   !CTM params
  INTEGER,                              INTENT(IN)    :: iter                    !Time evolution -- iter we are on
  REAL(KIND=DP),                        INTENT(IN)    :: time                    !Time evolution -- time we are at (corresponds to iter)
  LOGICAL,                              INTENT(IN)    :: use_last_env            !Whether we've used previous environment as initialization 

  INTEGER              :: Nop, i           !operator indexing
  LOGICAL, ALLOCATABLE :: obs_converged(:) !convergence tracking for individual obs
  REAL(KIND=DP)        :: etaR             !rescaled eta
  COMPLEX(KIND=DP)     :: norm             !Norm of TN

  !number of operators
  Nop = SIZE(TRACE_OP_1P, 1) 

  !Convergence tracking for diff obs
  ALLOCATE(obs_converged(Nop)); obs_converged(:) = .FALSE.

  !Allocate storage for transient observables
  IF(.NOT. ALLOCATED(new_obs)) ALLOCATE(new_obs(Nop))
  IF(.NOT. ALLOCATED(old_obs)) ALLOCATE(old_obs(Nop))

  !Rescaled convergence cutoff
  etaR = (tebd % eta) * (tebd % test_interval) * ABS(tebd % dt)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Compute all observables on the list i={1,nop}
  DO i=1,Nop
     CALL simplified_compute_obs_1P(new_obs(i), norm, ipeps, Cmat, Tmat, CTMRG, TRACE_OP_1P(i,:,:))
     CALL test_convergence(obs_converged(i), new_obs(i), old_obs(i), & 
                                                         & etaR, TRACE_OUT_1P(i), iter, tebd % test_interval, time, norm, use_last_env) 
  end DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Check convergence
  is_converged = isTrueVec(obs_converged, 'AND') 

  !Set old --> old = new for next iter
  old_obs = new_obs 

 END SUBROUTINE simplified_test_tebd_conv


END MODULE imaginary_tebd_callback_routines_ipeps
