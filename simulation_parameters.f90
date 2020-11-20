MODULE simulation_parameters

 USE utility, ONLY: OptArg, DP

 IMPLICIT NONE

 !Local (physical) dim of Quantum MPS, MPO -- Liouville/Hilbert space dim
 !trace_basis -- needed by pauli and gell mann modules
 INTEGER       :: local_dim=-1
 INTEGER       :: hs_dim=-1
 REAL(KIND=DP) :: trace_basis 

 !MPS chain: infinite/finite & pbc/obc
 LOGICAL :: pbc=.TRUE.
 LOGICAL :: infinite=.TRUE.

 !Version of MPS files
 LOGICAL :: use_old_formatting = .FALSE.

 !TEBD params
 TYPE tebd_params
   REAL(KIND=DP)    :: eps                          !SVD prec
   INTEGER          :: chi                          !SVD chi
   REAL(KIND=DP)    :: epsK                         !SVD prec for extra dim K
   INTEGER          :: chiK                         !SVD chi for extra dim K
   REAL(KIND=DP)    :: eta                          !TEBD convergence prec
   COMPLEX(KIND=DP) :: dt                           !Timestep
   INTEGER          :: N_steps                      !Max number of timesteps
   INTEGER          :: test_interval                !Testing interval
   INTEGER          :: min_sep, max_sep             !Min/max sep for corrs during convergence testing
   INTEGER          :: exit_level                   !Soft exit (exit_level = 1) vs Hard exit (exit_level = 2)
 END TYPE tebd_params

 !Precision params (for an iterative method with SVD truncation)
 TYPE prec_params
   REAL(KIND=DP) :: eps
   INTEGER       :: chi
   REAL(KIND=DP) :: eta
 END TYPE prec_params

CONTAINS

 !!! Setup the type that stores all TEBD params !!!
 SUBROUTINE create_tebd_params(TEBD, eps, chi, eta, dt, N_steps, test_interval, min_sep, max_sep, exit_level)

  TYPE(tebd_params), INTENT(INOUT)        :: TEBD                        !TEBD params
  REAL(KIND=DP),     INTENT(IN)           :: eps                         !SVD prec
  INTEGER,           INTENT(IN)           :: chi                         !SVD chi
  REAL(KIND=DP),     INTENT(IN)           :: eta                         !TEBD convergence prec
  COMPLEX(KIND=DP),  INTENT(IN)           :: dt                          !timestep
  INTEGER,           INTENT(IN)           :: N_steps                     !num of timesteps, 
  INTEGER,           INTENT(IN)           :: test_interval               !interval between convergence tests
  INTEGER,           INTENT(IN)           :: min_sep, max_sep            !min && max sep for checking convergence
  INTEGER,           INTENT(IN), OPTIONAL :: exit_level                  !specify how strictly we should check convergence

  !! Set params
  tebd % eps = eps
  tebd % chi = chi
  tebd % eta = eta

  tebd % epsK = 0.0D0
  tebd % chiK = -1

  tebd % dt            = dt
  tebd % N_steps       = N_steps
  tebd % test_interval = test_interval

  tebd % min_sep = min_sep
  tebd % max_sep = max_sep

  tebd % exit_level  = OptArg(exit_level, 2)  !! Set the level of convergence checks 
                                              !! (default = 2, i.e. the full convergence check)
 END SUBROUTINE create_tebd_params




 !!! Create a variable of type svd_params -- packs chi, eps together !!!
 SUBROUTINE create_prec_params(prec, eps, chi, eta)

  TYPE(prec_params), INTENT(INOUT) :: prec      !Precision params
  INTEGER,           INTENT(IN)    :: chi       !bond dim
  REAL(KIND=DP),     INTENT(IN)    :: eps, eta  !precision

  !! Set params
  prec % eps = eps
  prec % chi = chi
  prec % eta = eta

 END SUBROUTINE create_prec_params


END MODULE simulation_parameters
