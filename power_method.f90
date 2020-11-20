MODULE power_method

 USE utility
 USE definitions_mps_mpo
 USE iteration_helper
 USE mps_mpo_utility
 USE mps_mpo_algebra_finite

 IMPLICIT NONE

 !! Precision params (for an iterative method with SVD truncation)
 TYPE power_params
   REAL(KIND=DP) :: eps(2)
   INTEGER       :: chi(2)
   REAL(KIND=DP) :: eta
 END TYPE power_params

CONTAINS


 !!! Create a variable of type svd_params -- packs chi, eps together !!!
 SUBROUTINE create_power_params(POWER, eps, chi, eta)

  TYPE(power_params), INTENT(INOUT) :: POWER     !Power method params
  REAL(KIND=DP),      INTENT(IN)    :: eps(2)    !SVD precision
  INTEGER,            INTENT(IN)    :: chi(2)    !SVD bond dim
  REAL(KIND=DP),      INTENT(IN)    :: eta       !Convergence precision

  !! Set params
  POWER % eps = eps(2)
  POWER % chi = chi(2)
  POWER % eta = eta

 END SUBROUTINE create_power_params




 !!! Compute evecs of finite MPO using power method !!!
 SUBROUTINE calc_mpsevecs_power_meth(Rvec, Lvec, Tmat, POWER, use_last_evec, N_sites, dN, evecFlag) 

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: Rvec(:), Lvec(:)    !MPS evecs
  TYPE(block_mpo),              INTENT(IN)    :: Tmat(:,:)           !MPO Tmat
  TYPE(power_params),           INTENT(IN)    :: POWER               !Power method params
  LOGICAL,                      INTENT(IN)    :: use_last_evec       !Whether to use last evec as initial pt
  INTEGER,                      INTENT(IN)    :: N_sites             !Size of evec to calculate                
  INTEGER,                      INTENT(IN)    :: dN                  !Diff b/n [last evec size] and [new evec size N_sites] 
  CHARACTER(LEN=*),             INTENT(IN)    :: evecFlag            !Specify evec to calculate

  !! Transposed Tmat
  TYPE(block_mpo), ALLOCATABLE :: TmatTransp(:,:)

  !! Array of SN dims for initializing new evec from scratch
  INTEGER, ALLOCATABLE :: SNdim(:)
  INTEGER              :: site

  !! Set SNdim of evec MPS
  ALLOCATE(SNdim(N_sites))
  SNdim(1)         =    Tmat(2,1) % Ndim
  SNdim(2:N_sites) = (/(Tmat(2,2) % Ndim, site=2, N_sites)/)

  SELECT CASE(evecFlag)
  CASE('Rvec') 

     !! Find Rvec
     CALL initialize_mps_block(Rvec, SNdim=SNdim, chi=2, use_last_mps=use_last_evec, N_sites=N_sites, dN=dN)
     CALL calc_single_mpsevec_power_meth(Rvec, Tmat, POWER, use_last_evec, 'Rvec') 

  CASE('Lvec')

     !! Find Lvec
     CALL initialize_mps_block(Lvec, SNdim=SNdim, chi=2, use_last_mps=use_last_evec, N_sites=N_sites, dN=dN)
     CALL transpose_bimpo_block(TmatTransp, Tmat)

     CALL calc_single_mpsevec_power_meth(Lvec, TmatTransp, POWER, use_last_evec, 'Lvec') 

  CASE DEFAULT

     !! Find both evecs
     CALL initialize_mps_block(Rvec, SNdim=SNdim, chi=2, use_last_mps=use_last_evec, N_sites=N_sites, dN=dN)
     CALL initialize_mps_block(Lvec, SNdim=SNdim, chi=2, use_last_mps=use_last_evec, N_sites=N_sites, dN=dN)

     CALL calc_single_mpsevec_power_meth(Rvec, Tmat, POWER, use_last_evec, 'Rvec') 
     CALL calc_single_mpsevec_power_meth(Lvec, Tmat, POWER, use_last_evec, 'Lvec') 
  end SELECT

 END SUBROUTINE calc_mpsevecs_power_meth







 !!! Compute single evec of finite MPO using power method !!!
 SUBROUTINE calc_single_mpsevec_power_meth(evec, Tmat, POWER, use_last_evec, evecFlag) 

  USE datafile_utility,     ONLY: setupConvDatafile
  USE mps_peps_INOUT_files, ONLY: mps_print_function

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: evec(:)             !MPS evec
  TYPE(block_mpo),              INTENT(IN)    :: Tmat(:,:)           !MPO Tmat
  TYPE(power_params),           INTENT(IN)    :: POWER               !Power method params
  LOGICAL,                      INTENT(IN)    :: use_last_evec       !Whether to use last evec as initial pt
  CHARACTER(LEN=*),             INTENT(IN)    :: evecFlag            !Specify evec to calculate

  !! Evals && norm
  COMPLEX(KIND=DP) :: eval, evalOLD
  COMPLEX(KIND=DP) :: C_norm

  !! Evec from previous iter
  TYPE(block_mps), ALLOCATABLE :: evecOLD(:), evecOLD_HC(:)
  INTEGER                      :: N_sites

  !! Power method settings & params
  INTEGER, PARAMETER :: N_iters = 300
  INTEGER            :: datafile, iter 
  LOGICAL            :: is_converged

  !! Get evec size
  N_sites = SIZE(evec)

  !! Setup file for convergence data
  CALL setupConvDatafile(datafile=datafile, descr=TRIM(ADJUSTL(evecFlag))//"_power_meth", &
       & errstr="Err_eval", valstr="eval", xstr="C_norm", xlog1="reuse_evec", chi=POWER % chi(2), N_sites=N_sites)

  !! Init convergence params
  evalOLD = (0.0D0, 0.0D0); is_converged = .FALSE.

  !! Starting power method
  WRITE(*,*) "Starting power method for N_sites = ", N_sites

  powerloop:DO iter=1,N_iters

     !! Create copies of evec
     CALL copy_mps_block(evecOLD,    evec)
     CALL copy_mps_block(evecOLD_HC, evecOLD, 'HC')

     !! Multiply Evec by Tmat
     CALL multiply_mps_bimpo(evec, Tmat, POWER % chi, POWER % eps)

     !! Calc evals
     CALL contract_mps_mps(eval,   evec,    evecOLD_HC)
     CALL contract_mps_mps(C_norm, evecOLD, evecOLD_HC)
     eval = eval / C_norm
     
     !! Test convergence
     CALL test_convergence(is_converged, new_val=eval, old_val=evalOLD, eta=POWER % eta, &
                           & datafile=datafile, iter=iter, xval=C_norm, xlog=use_last_evec)

     !! Normalize evec
     CALL normalize_mps(evec)

     !! If converged: exit loop && output evec
     IF(is_converged) THEN                      
         CALL mps_print_function(evec, evecFlag)
         EXIT powerloop
     ELSEIF(iter .EQ. N_iters) THEN   
         WRITE(*,*) "POWER METHOD has failed to converge after N_iter = ", iter
         STOP
     end IF

     !! Pass eval from [iter-1] to [iter]
    evalOLD = eval

  END DO powerloop

  !! Close files
  CLOSE(datafile) 

 END SUBROUTINE calc_single_mpsevec_power_meth


END MODULE power_method
