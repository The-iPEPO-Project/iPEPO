MODULE tebd_callback_routines_impo

 USE utility
 USE iteration_helper
 USE simulation_parameters
 USE observables_impo

 IMPLICIT NONE

 INTEGER :: DUMP_NUM, G2_SUMMARY

CONTAINS

 !!! Test if both 1P && 2P observables have converged !!!
 SUBROUTINE test_tebd_convergence(is_converged, new_obs_2P, old_obs_2P, new_obs_1P, old_obs_1P, rho_onesite, rho_twosite, tebd, iter, time)

  LOGICAL,                       INTENT(INOUT) :: is_converged                        !Whether observables have converged
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: new_obs_2P(:,:), old_obs_2P(:,:)    !New && old values of observables (to be compared)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: new_obs_1P(:),   old_obs_1P(:)      !New && old values of observables (to be compared)
  COMPLEX(KIND=DP),              INTENT(IN)    :: rho_onesite(:,:)                    !One-site reduced density matrix
  COMPLEX(KIND=DP),              INTENT(IN)    :: rho_twosite(:,:,:,:,:)              !Two-site reduced density matrix
  TYPE(tebd_params),             INTENT(IN)    :: tebd                                !Time evolution -- propagation params (dt, N_steps, test_interval)
  INTEGER,                       INTENT(IN)    :: iter                                !Time evolution -- iter we are at
  REAL(KIND=DP),                 INTENT(IN)    :: time                                !Time evolution -- time we are at (corresponds to iter)

  LOGICAL :: is_converged_1P, is_converged_2P

  !! Initialize
  is_converged = .FALSE.

  !! Check 1P observables
  CALL test_tebd_conv_1P(is_converged_1P, new_obs_1P, old_obs_1P, & 
                          & rho_onesite,              tebd, iter, time)

  !! Check 2P observables
  CALL test_tebd_conv_2P(is_converged_2P, new_obs_2P, old_obs_2P, &
                          & rho_onesite, rho_twosite, tebd, iter, time)

  !! Both 1P && 2P must be converged
  is_converged = is_converged_1P .AND. is_converged_2P

 END SUBROUTINE test_tebd_convergence



 !!! Test if 1P observables have converged !!!
 SUBROUTINE test_tebd_conv_1P(is_converged, new_obs, old_obs, rho_onesite, tebd, iter, time)

  USE datafile_utility

  LOGICAL,                       INTENT(INOUT) :: is_converged            !whether observables have converged
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: new_obs(:), old_obs(:)  !new && old values of observables (to be compared)
  COMPLEX(KIND=DP),              INTENT(IN)    :: rho_onesite(:,:)        !One-site reduced density matrix for computing 1P observables
  TYPE(tebd_params),             INTENT(IN)    :: tebd                    !Time evolution -- propagation params (dt, N_steps, test_interval)
  INTEGER,                       INTENT(IN)    :: iter                    !Time evolution -- iter we are on
  REAL(KIND=DP),                 INTENT(IN)    :: time                    !Time evolution -- time we are at (corresponds to iter)

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

  !Get norm of rho TN (i.e. trace)
  norm = TRACE(rho_onesite) 

  !Compute all observables on the list i={1,nop}
  DO i=1,Nop
     CALL compute_obs_1P(new_obs(i), rho_onesite, TRACE_OP_1P(i,:,:))
     new_obs(i) = new_obs(i)/norm
     CALL test_convergence(obs_converged(i), new_obs(i), old_obs(i), &
                                                         & etaR, TRACE_OUT_1P(i), iter, tebd % test_interval, time, norm)
  end DO

  !Check convergence
  is_converged = isTrueVec(obs_converged, 'AND')

  !Set old --> old = new for next iter
  old_obs = new_obs

 END SUBROUTINE test_tebd_conv_1P




 !!! Test if 2P observables have converged !!!
 SUBROUTINE test_tebd_conv_2P(is_converged, new_obs, old_obs, rho_onesite, rho_twosite, tebd, iter, time)

  USE datafile_utility

  LOGICAL,                       INTENT(INOUT) :: is_converged                 !Whether observables have converged
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: new_obs(:,:), old_obs(:,:)   !New && old values of observables (to be compared)
  COMPLEX(KIND=DP),              INTENT(IN)    :: rho_onesite(:,:)             !One-site reduced density matrix 
  COMPLEX(KIND=DP),              INTENT(IN)    :: rho_twosite(:,:,:,:,:)       !Two-site reduced density matrix 
  TYPE(tebd_params),             INTENT(IN)    :: tebd                         !Time evolution -- propagation params (dt, N_steps, test_interval)
  INTEGER,                       INTENT(IN)    :: iter                         !Time evolution -- iter we are at
  REAL(KIND=DP),                 INTENT(IN)    :: time                         !Time evolution -- time we are at (corresponds to iter)

  INTEGER              :: Nop, i                          !operator indexing
  INTEGER              :: min_sep, max_sep, Nsep, sep     !separations
  LOGICAL, ALLOCATABLE :: obs_converged(:,:)              !convergence tracking for individual obs
  REAL(KIND=DP)        :: etaR                            !rescaled eta
  COMPLEX(KIND=DP)     :: norm                            !Norm of TN

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

  !Get norm of rho TN (i.e. trace)
  norm = TRACE(rho_onesite) 

  !Compute all 2P observables on the list i={1,nop}
  DO i=1,Nop

     !compute obs
     CALL compute_obs_2P(new_obs(i,:), rho_twosite, rho_onesite, TRACE_OP_2P(i,:,:,:), min_sep, max_sep)

     !normalize obs, and test convergence
     DO sep = min_sep, max_sep
        new_obs(i, sep+1-min_sep) = new_obs(i, sep+1-min_sep) / norm**(sep/2 + 1)
        CALL test_convergence(obs_converged(i,sep), new_obs(i,sep), old_obs(i,sep), &
                                                                    & etaR, TRACE_OUT_2P(i,sep), iter, tebd % test_interval, time, norm)
     end DO
  end DO

  !Check convergence
  is_converged = isTrueMat(obs_converged, 'AND')

  !Set old --> old = new for next iter
  old_obs = new_obs

 END SUBROUTINE test_tebd_conv_2P




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





 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Decide which param to use in the param-loop, based on "xname" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

  !! Local copies of optional string vars
  CHARACTER(LEN=16) :: v2_strX, v3_strX, v4_strX

  !! Set default values (If var not present -- set its string to INVALID_STRING, to guarantee it's not selected)
  v2_strX = OptArg(v2_str, "INVALID_STRING")
  v3_strX = OptArg(v3_str, "INVALID_STRING")
  v4_strX = OptArg(v4_str, "INVALID_STRING")
  
  !! Decide which param to use
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

END MODULE tebd_callback_routines_impo
