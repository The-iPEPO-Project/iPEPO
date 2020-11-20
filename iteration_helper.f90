MODULE iteration_helper

 USE utility
 USE definitions_mps_mpo
 USE definitions_peps_pepo

 IMPLICIT NONE

 !initialize using old mps
 interface initialize_using_old_mps
    module procedure extend_mps
    module procedure initialize_recycled_mps
    module procedure initialize_to_last_imps
 end interface initialize_using_old_mps

 !initialize random mps block
 interface initialize_rand_mps
    module procedure init_rand_mps_finite
    module procedure init_rand_imps_gammas
 end interface initialize_rand_mps

 !initialize factorized mps block
 interface initialize_factorized_mps
    module procedure init_factorized_mps_non_uniform
    module procedure init_factorized_mps_uniform
 end interface initialize_factorized_mps

 !initialize random tensor
 interface initialize_rand_tens
    module procedure initialize_rand_tens_1D
    module procedure initialize_rand_tens_2D
    module procedure initialize_rand_tens_3D
    module procedure initialize_rand_tens_4D
    module procedure initialize_rand_tens_5D
 end interface initialize_rand_tens


 private  !! hides all items not listed in public statement 
 public test_convergence, test_convergence_vec, vec_error, test_cost_function

 !! MPS 
 public initialize_mps_block, initialize_imps
 public initialize_using_old_mps, initialize_rand_mps, initialize_factorized_mps
 public initialize_rand_full_imps

 !! MPO
 public initialize_factorized_mpo
 public initialize_rand_full_impo

 !! Lambda
 public initialize_factorized_lambda, initialize_rand_lambda

 !! PEPS
 public initialize_rand_peps, initialize_factorized_peps
 public initialize_rand_full_ipeps

 !! Boundary vectors
 public initialize_bvec, bvec_dims
 public initialize_rand_tens

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Single-value Convergence !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Test convergence of a value !!!
 SUBROUTINE test_convergence(is_converged, new_val, old_val, eta, datafile, iter, test_interval, iter_param, xval, xlog)

  USE datafile_utility, ONLY: write_conv_data

  LOGICAL,          INTENT(INOUT)        :: is_converged       !whether or not val has converged
  COMPLEX(KIND=DP), INTENT(IN)           :: new_val, old_val   !new vs old vals
  REAL(KIND=DP),    INTENT(IN)           :: eta                !convergence cutoff
  INTEGER,          INTENT(INOUT)        :: datafile           !datafile for writing convergence data
  INTEGER,          INTENT(IN)           :: iter               !iteration we are at
  INTEGER,          INTENT(IN), OPTIONAL :: test_interval      !test interval
  REAL(KIND=DP),                OPTIONAL :: iter_param         !optional args
  COMPLEX(KIND=DP), INTENT(IN), OPTIONAL :: xval               !optional args
  LOGICAL,          INTENT(IN), OPTIONAL :: xlog               !optional args

  !! Convergence errors && test interval
  REAL(KIND=DP) :: Err_val

  !! No convergence check when iter=1 
  IF(iter .LE. OptArg(test_interval, 1)) RETURN 

  !! (1) Find relative error in eval (diff between evals at current & previous steps) 
  IF(ABS(new_val) .GT. eta) THEN
     Err_val = ABS(new_val - old_val)/ABS(new_val)
  ELSE
     Err_val = ABS(new_val - old_val)
  end IF

  !! (2) Check if eval & evec convergence criteria have been fulfilled 
  is_converged = (Err_val .LT. eta)

  !! (3) record convergence data in a file
  CALL write_conv_data(datafile=datafile, iter=iter, ipar=iter_param, errval=Err_val, val=new_val, xval=xval, xlog1=xlog)

 END SUBROUTINE test_convergence


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Vector Convergence !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Test convergence of a vector of values !!!
 SUBROUTINE test_convergence_vec(is_converged, new_val, old_val, eta, datafile, iter, test_interval, xval, xlog, PRINT_ON)

  USE datafile_utility, ONLY: write_conv_data

  LOGICAL,            INTENT(INOUT)        :: is_converged             !whether or not values have converged           
  COMPLEX(KIND=DP),   INTENT(IN)           :: new_val(:), old_val(:)   !new vs old values (vectors)
  REAL(KIND=DP),      INTENT(IN)           :: eta                      !convergence cutoff
  INTEGER,            INTENT(INOUT)        :: datafile                 !file handle to output conv data
  INTEGER,            INTENT(IN)           :: iter                     !num of iterations
  INTEGER,            INTENT(IN), OPTIONAL :: test_interval            !test interval        
  COMPLEX(KIND=DP),   INTENT(IN), OPTIONAL :: xval                     !optional args
  LOGICAL,            INTENT(IN), OPTIONAL :: xlog                     !optional args
  LOGICAL,            INTENT(IN), OPTIONAL :: PRINT_ON

  !! Convergence errors, overlapping length of new_val && old_val arrays
  REAL(KIND=DP) :: Err_val
  INTEGER       :: nval

  !! No convergence tests on first iter
  IF(iter .LE. OptArg(test_interval, 1)) RETURN

  !! Get the overlapping length of new_val && old_val arrays
  nval = MIN(SIZE(new_val), SIZE(old_val))

  !! Find [new error] = [MAX relative error in the array of values] 
  !! -- i.e. compare vecs at current && previous steps
  Err_val      = vec_error(new_val(1:nval), old_val(1:nval)) 
  is_converged = Err_val .LT. eta

  !! Record convergence data in a file
  IF(OptArg(PRINT_ON, .FALSE.)) CALL write_conv_data(datafile=datafile, iter=iter, errval=Err_val, &
                                                     & val=complx(MAXVAL(ABS(new_val(1:nval)))), xval=xval, xlog1=xlog)

 END SUBROUTINE test_convergence_vec





 !!! Vector convergence, relative error: vecA vs vecB !!!
 FUNCTION vec_error(vecA, vecB)

  COMPLEX(KIND=DP), INTENT(IN) :: vecA(:), vecB(:)
  REAL(KIND=DP)                :: vec_error 

  INTEGER :: sizeA, sizeB, minDim

  sizeA  = SIZE(vecA)
  sizeB  = SIZE(vecB)
  minDim = MIN(sizeA, sizeB)

  vec_error = MAXVAL(ABS(vecA(1:minDim) - vecB(1:minDim))) / MAXVAL(ABS(vecA(1:minDim))) 

 END FUNCTION vec_error


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Convergence using Cost Function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Test convergence using cost function !!!
 SUBROUTINE test_cost_function(is_converged, costF, eta, datafile, iter, iter_param, xval, xlog)

  USE datafile_utility, ONLY: write_conv_data

  LOGICAL,          INTENT(INOUT)        :: is_converged  !whether or not val has converged
  COMPLEX(KIND=DP), INTENT(IN)           :: costF         !new vs old vals
  REAL(KIND=DP),    INTENT(IN)           :: eta           !convergence cutoff
  INTEGER,          INTENT(INOUT)        :: datafile      !datafile for writing convergence data
  INTEGER,          INTENT(IN)           :: iter          !iteration we are at
  REAL(KIND=DP),                OPTIONAL :: iter_param    !optional args
  COMPLEX(KIND=DP), INTENT(IN), OPTIONAL :: xval          !optional args (should contain a rescaling factor for costF, such as eval)
  LOGICAL,          INTENT(IN), OPTIONAL :: xlog          !optional args

  !Relative error
  REAL(KIND=DP) :: Err_val

  !No convergence check when iter=1
  IF(iter .EQ. 1) RETURN

  !Set relative error (rescale by xval if provided -- xval input should contain a rescaling factor, such as eval)
  Err_val = ABS(costF)
  IF(PRESENT(xval)) Err_val = Err_val/ABS(xval)**2

  !Check if costF has converged
  is_converged = Err_val .LT. eta

  !Record convergence data in a file
  CALL write_conv_data(datafile=datafile, iter=iter, ipar=iter_param, errval=Err_val, val=costF, xval=xval, xlog1=xlog)

 END SUBROUTINE test_cost_function


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














 !!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize MPS BLOCK or iMPS -- Top-level routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Initialize mps block !!!
 SUBROUTINE initialize_mps_block(mps, SNdim, chi, use_last_mps, N_sites, dN)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps(:)
  INTEGER,                      INTENT(IN)    :: SNdim(:), chi !mps site dims of evec (if initialized from scratch)
  LOGICAL,                      INTENT(IN)    :: use_last_mps  !whether to reuse last MPS
  INTEGER,                      INTENT(IN)    :: N_sites, dN   !N_sites = size of new evec we're initializing, 
                                                               !dN = size diff from old evec

  !!!!! setup new evec from scratch (IF clause) or insert extra sites into existing evec (ELSE clause) !!!!!!!
  IF(.NOT. use_last_mps) THEN 

     !Initialize MPS to a rand state
     CALL initialize_rand_mps(mps, SNdim, chi, N_sites)

  ELSEIF(dN .EQ. 0) THEN

     !! Recycle evec to ensure SNdim remains compatible
     CALL initialize_using_old_mps(mps, SNdim)
  ELSE

     !! Use an existing evec as a basis for new evec: take last converged evec & insert [dN] copies  
     !! of its mid-site into its mid part to extend it to the size of new evec = [N_new] = [N_sites]
     CALL initialize_using_old_mps(mps, N_sites, dN)                                          
  end IF

 END SUBROUTINE initialize_mps_block







 !!! Initialize (infinite) iMPS !!!
 SUBROUTINE initialize_imps(iGamma, iLambda, SNdim, chi, use_last_mps)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  INTEGER,                         INTENT(IN)    :: SNdim, chi !mps site dims of evec (if initialized from scratch)
  LOGICAL,                         INTENT(IN)    :: use_last_mps

  !Local iMPS copy
  TYPE(block_mps),    ALLOCATABLE :: iGammaCpy(:)
  TYPE(block_lambda), ALLOCATABLE :: iLambdaCpy(:)
  INTEGER :: EWdim, site

  !Vector along SNdim at each site of iMPS
  COMPLEX(KIND=DP), ALLOCATABLE :: vecSN(:) 

  IF(.NOT. use_last_mps) THEN

     !Initialize random Gammas && Lambdas
     ALLOCATE(vecSN(SNdim)); vecSN(:) = (1.0d0,0.0d0)

     CALL initialize_factorized_mps(iGamma, SNdim, 2, vecSN) 
     CALL initialize_factorized_lambda(iLambda, 2)

     !CALL initialize_rand_mps(iGamma, SNdim, chi)
     !CALL initialize_rand_lambda(iLambda, chi)
  ELSE  
     !Initialize to the last converged evec
     CALL initialize_using_old_mps(iGamma, iLambda)
  end IF

 END SUBROUTINE initialize_imps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize MPS block or iMPS using an old MPS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Insert sites into an existing MPS BLOCK !!!
 SUBROUTINE extend_mps(mps_block, N_new, dN)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:)
  INTEGER,                      INTENT(IN)    :: N_new, dN !size of new mps, number of sites to insert

  TYPE(block_mps), ALLOCATABLE :: mps_old(:)   !copy of old mps
  TYPE(block_mps)              :: ExtraSite    !copy of site we intend to insert

  !! Dims & indices
  INTEGER :: N_old
  INTEGER :: x_ins, x_start, x_end
  INTEGER :: SNdim, min_Wdim, Edim
  INTEGER :: site

  !! Create a copy of old mps:
  CALL copy_mps_block(mps_old, mps_block)

  !! Size of old mps; check that pre- && post- insertion sizes are consistent
  N_old = SIZE(mps_old)
  CALL check_sizes_equal(N_new, N_old + dN, "extend_mps_block: N_new and N_old + dN must be equal ")

  !! Pos where we do insertion. Insert dN copies of [ExtraSite]: start && end of insertion:
  x_ins   = N_old
  x_start = x_ins + 1
  x_end   = x_ins + dN

  !! Take mps_site at [site=x_ins] of [mps_block], create its copy
  CALL allocate_mps_site(ExtraSite, mps_old(x_ins) % m)

  !! Dims of insertion
  SNdim    = ExtraSite % SNdim
  min_Wdim = MIN(ExtraSite % Wdim, ExtraSite % Edim)
  Edim     = ExtraSite % Edim

  !! Allocate new empty mps block
  CALL allocate_empty_mps_block(mps_block, N_new)

  !! Insert [block_to_insert] into [mps_block]
  DO site=1,N_new

      !! Do the insertion:
      IF(site .LT. x_start) THEN

           !! Copy the first half of mps_old to mps_block
           CALL allocate_mps_site(mps_block(site), mps_old(site) % m)

      ELSEIF(site .GT. x_end) THEN

           !! Copy the last half of mps_old to mps_block
           CALL allocate_mps_site(mps_block(site), mps_old(site - dN) % m)
      ELSE

           !! Insert copy of the middle site to the mid-part of mps_block
           CALL allocate_mps_site(mps_block(site), SNdim, Edim, Edim)
           mps_block(site) % m(:, 1:min_Wdim, :) = ExtraSite % m(:, 1:min_Wdim, :)
      end IF 
  end DO

 END SUBROUTINE extend_mps






 !!! Initialize recycled mps_block 
 !!! (i.e. with new SNdim -- MPS recycled by mults along SN axis)
 SUBROUTINE initialize_recycled_mps(mps_block, SNdim)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:)
  INTEGER,                      INTENT(IN)    :: SNdim(:)

  !! Dims & indices
  INTEGER :: N_sites, site
  
  !! Temp tensor for recycled sites of mps_old
  COMPLEX(KIND=DP), ALLOCATABLE :: tempTens(:,:,:)

  !! Get size of mps block
  N_sites = SIZE(mps_block)

  !! Recycle all sites of mps_block
  DO site=1,N_sites 
     CALL copyTens(tempTens, mps_block(site) % m,  (/SNdim(site), mps_block(site) % Wdim, mps_block(site) % Edim/))
     CALL allocate_mps_site(mps_block(site), tempTens)
  END DO

 END SUBROUTINE initialize_recycled_mps








 !!! Initialize iMPS evec to the last converged evec 
 !!! (with rounded-off values)
 SUBROUTINE initialize_to_last_imps(iGamma, iLambda)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)

  !! Local iMPS copy
  TYPE(block_mps),    ALLOCATABLE :: iGammaCpy(:)
  TYPE(block_lambda), ALLOCATABLE :: iLambdaCpy(:)

  COMPLEX(KIND=DP), ALLOCATABLE :: tempGamma(:,:,:), tempLambda(:)
  INTEGER :: EWdim(2), minEWdim, site

  !! Rounding factor
  REAL(KIND=DP), PARAMETER :: RoundFac = 1.0d-04

  !! Round off the values in iMPS
  DO site=1,2

     !! Round off iGamma, iLambda
     CALL roundTens(iGamma(site) % m,  RoundFac)
     CALL roundTens(iLambda(site) % m, RoundFac)

     !! Get the new bond dim (number of nonzero lambdas)
     EWdim(site) = COUNT(ABS(iLambda(site) % m) .GT. RoundFac)
  end DO

  !! Use the smaller one of two bond dims at sites=1,2
  !! Must guarantee that sites=1,2 have equal chi, and no zero lambdas
  minEWdim = MIN(EWdim(1), EWdim(2)) 

  !! Copy rounded tensors into iMPS, ensuring that bond dims match
  DO site=1,2

     !! Cut the dims of G,L to exclude any zero lambdas (i.e. to avoid non-invertible boundary matrices)
     CALL copyTens(tempGamma,  iGamma(site) % m,   (/ iGamma(site) % SNdim, minEWdim, minEWdim /))
     CALL copyTens(tempLambda, iLambda(site) % m,  (/ minEWdim /))

     !! Allocate new iMPS
     CALL allocate_mps_site(iGamma(site), tempGamma)
     CALL allocate_lambda_site(iLambda(site), tempLambda)
  END DO

 END SUBROUTINE initialize_to_last_imps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize random MPS block or iMPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Generate random finite mps_block !!!
 SUBROUTINE init_rand_mps_finite(mps_block, SNdim, chi, N_sites)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:)
  INTEGER,                      INTENT(IN)    :: SNdim(:), chi, N_sites

  !! Dims & indices
  INTEGER :: dimA(3)
  INTEGER :: site, i, j, k

  !! Rand seed
  LOGICAL, PARAMETER :: use_const_seed = .TRUE.
  REAL               :: cpu_times 
  INTEGER            :: seed 

  !! Make sure SNdim input array has the same size as mps_block
  CALL check_sizes_equal(SIZE(SNdim), N_sites, "init_rand_mps_finite -- SNdim must have length = N_sites")

  !! Generate rand nums using a const seed or a seed from cpu_time
  IF(use_const_seed) THEN
     seed = 86456
  ELSE
     CALL cpu_time(cpu_times)
     seed = NINT(cpu_times)
  end IF
  CALL srand(seed)

  !! Allocate empty evec
  CALL allocate_empty_mps_block(mps_block, N_sites)

  !! Loop over all sites && do initialization
  DO site=1,N_sites 

     !! Get dims of mps_block
     IF(site .EQ. 1) THEN
            dimA = (/SNdim(site), 1, chi/)
     ELSEIF(site .EQ. N_sites) THEN
            dimA = (/SNdim(site), chi, 1/)
     ELSE
            dimA = (/SNdim(site), chi, chi/)
     end IF

     !! Initialize as a random number state
     CALL allocate_mps_site(mps_block(site), dimA) 
        
     Do i=1,dimA(1)
       Do j=1,dimA(2)
         Do k=1,dimA(3)
            mps_block(site) % m(i,j,k) = CMPLX(rand(), 0d0)
         end Do
       end Do
     end Do
  END DO

 END SUBROUTINE init_rand_mps_finite




 !!! Generate random MPS gamma block !!!
 SUBROUTINE init_rand_imps_gammas(mps_block, SNdim, chi, N_sites_in)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT)        :: mps_block(:)
  INTEGER,                      INTENT(IN)           :: SNdim, chi
  INTEGER,                      INTENT(IN), OPTIONAL :: N_sites_in

  !! Temp site tensor
  COMPLEX(KIND=DP), ALLOCATABLE :: tempSite(:,:,:)

  !! Dims & indices
  INTEGER :: i, j, k
  INTEGER :: N_sites, site

  !! Rand seed
  LOGICAL, PARAMETER :: use_const_seed = .TRUE.
  REAL :: cpu_times 
  INTEGER :: seed 

  !! Set optional args
  N_sites = OptArg(N_sites_in, 2)

  !! Generate rand nums using a const seed or a seed from cpu_time
  IF(use_const_seed) THEN
     seed = 86456
  ELSE
     CALL cpu_time(cpu_times)
     seed = NINT(cpu_times)
  end IF
  CALL srand(seed)

  !! Init random site
  ALLOCATE(tempSite(SNdim, chi, chi)) 

  Do i=1,SNdim
    Do j=1,chi  
      Do k=1,chi  
         tempSite(i,j,k) = CMPLX(rand(), 0d0)
      end Do
    end Do
  end Do

  !! Allocate empty imps
  CALL allocate_empty_mps_block(mps_block, N_sites)

  !! Generate random sites
  DO site=1,N_sites
     CALL allocate_mps_site(mps_block(site), tempSite)
  END DO

 END SUBROUTINE init_rand_imps_gammas

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize factorized MPS block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Initialize factorized mps_block (with chi=1) !!!
 SUBROUTINE init_factorized_mps_non_uniform(mps_block, SNdim, N_sites)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:)
  INTEGER,                      INTENT(IN)    :: SNdim(:), N_sites

  !vecA && bond dim for initializing a product state
  INTEGER :: chi, site

  !Bond dim of a factorized (product) MPS
  chi = 1

  !Allocate factorized mps
  CALL allocate_empty_mps_block(mps_block, N_sites)

  DO site=1,N_sites
     CALL allocate_mps_site(mps_block(site),  (/SNdim(site), chi, chi/))
     mps_block(site) % m(:,1,1) = (1.0d0,0.0d0)
  end DO

 END SUBROUTINE init_factorized_mps_non_uniform




 !!! Initialize factorized MPS gamma block !!!
 SUBROUTINE init_factorized_mps_uniform(mps_block, SNdim, N_sites, vecA_In)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT)        :: mps_block(:)
  INTEGER,                      INTENT(IN)           :: SNdim, N_sites
  COMPLEX(KIND=DP),             INTENT(IN), OPTIONAL :: vecA_In(:)

  !! vecA && bond dim for initializing a product state
  COMPLEX(KIND=DP), ALLOCATABLE :: vecA(:)
  INTEGER :: chi, site

  !! Bond dim of a factorized (product) MPS
  chi = 1

  !! Find vecA
  IF(PRESENT(vecA_In)) THEN
     CALL copyTens(vecA, vecA_In)
  ELSE
     ALLOCATE(vecA(SNdim))
     vecA(:) = (0.0D0, 0.0D0)
     vecA(1) = (1.0D0, 0.0D0)
  end IF

  !! Allocate product state gammas
  CALL allocate_empty_mps_block(mps_block, N_sites)

  DO site=1,N_sites
     CALL allocate_mps_site(mps_block(site),  (/SNdim, chi, chi/))
     mps_block(site) % m(:,1,1) = vecA
  end DO

 END SUBROUTINE init_factorized_mps_uniform


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize complete iMPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE initialize_rand_full_imps(iGamma, iLambda, localDim, chi)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  INTEGER,                         INTENT(IN)    :: localDim, chi

  !! Rand tensors
  COMPLEX(KIND=DP), ALLOCATABLE :: tmpLAM(:), tmpGAM(:,:,:)

  !! Num of sites
  INTEGER, PARAMETER :: N_sites = 2
  INTEGER            :: site

  !! Norm fac
  COMPLEX(KIND=DP) :: fac

  !! Init random site tensors
  ALLOCATE(tmpLAM(chi), tmpGAM(localDim, chi, chi))
  CALL initialize_rand_tens(tmpGAM)
  CALL initialize_rand_tens(tmpLAM)

  !! Normalize lambda
  fac = SQRT(SUM(tmpLAM**2))
  tmpLAM = tmpLAM / fac
  tmpGAM = tmpGAM * SQRT(fac)**2

  !! Sort lambdas
  CALL sort_array(tmpLAM)

  !! Setup random iGamma 
  CALL allocate_empty_mps_block(iGamma, N_sites)
  DO site=1,N_sites
     CALL allocate_mps_site(iGamma(site), tmpGAM)
  END DO

  !! Setup random iLambda
  CALL allocate_empty_lambda_block(iLambda, N_sites)
  DO site=1,N_sites
     CALL allocate_lambda_site(iLambda(site), tmpLAM)
  END DO

 END SUBROUTINE initialize_rand_full_imps 


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize MPO block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Initialize factorized MPO block !!!
 SUBROUTINE initialize_factorized_mpo(mpo_block, Sdim, Ndim, N_sites, matA_In)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT)        :: mpo_block(:)
  INTEGER,                      INTENT(IN)           :: Sdim, Ndim, N_sites
  COMPLEX(KIND=DP),             INTENT(IN), OPTIONAL :: matA_In(:,:)

  !! matA && bond dim for initializing a product state
  COMPLEX(KIND=DP), ALLOCATABLE :: matA(:,:)
  INTEGER :: chi, site

  !! Bond dim of a factorized (product) MPO
  chi = 1

  !! Find matA
  IF(PRESENT(matA_In)) THEN
     CALL copyTens(matA, matA_In)
  ELSE
     ALLOCATE(matA(Sdim, Ndim))
     matA(:,:) = (0.0D0, 0.0D0)
     matA(1,1) = (1.0D0, 0.0D0)
  end IF

  !! Allocate product state gammas
  CALL allocate_empty_mpo_block(mpo_block, N_sites)

  DO site=1,N_sites
     CALL allocate_mpo_site(mpo_block(site),  (/Sdim, Ndim, chi, chi/))
     mpo_block(site) % m(:,:,1,1) = matA
  end DO

 END SUBROUTINE initialize_factorized_mpo




 !!! Initialize random iMPO !!!
 SUBROUTINE initialize_rand_full_impo(iGamma, iLambda, Sdim, Ndim, chi) 

  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  INTEGER,                         INTENT(IN)    :: Sdim, Ndim
  INTEGER,                         INTENT(IN)    :: chi

  !! Rand tensors
  COMPLEX(KIND=DP), ALLOCATABLE :: tmpLAM(:), tmpGAM(:,:,:,:)

  !! Num of sites
  INTEGER, PARAMETER :: N_sites = 2
  INTEGER            :: site

  !! Norm fac
  COMPLEX(KIND=DP) :: fac

  !! Init random site tensors
  ALLOCATE(tmpLAM(chi), tmpGAM(Sdim, Ndim, chi, chi))
  CALL initialize_rand_tens(tmpGAM)
  CALL initialize_rand_tens(tmpLAM)

  !! Normalize lambda
  fac = SQRT(SUM(tmpLAM**2))
  tmpLAM = tmpLAM / fac
  tmpGAM = tmpGAM * SQRT(fac)**2

  !! Sort lambdas
  CALL sort_array(tmpLAM)

  !! Setup random iGamma 
  CALL allocate_empty_mpo_block(iGamma, N_sites)
  DO site=1,N_sites
     CALL allocate_mpo_site(iGamma(site), tmpGAM)
  END DO

  !! Setup random iLambda
  CALL allocate_empty_lambda_block(iLambda, N_sites)
  DO site=1,N_sites
     CALL allocate_lambda_site(iLambda(site), tmpLAM)
  END DO

 END SUBROUTINE initialize_rand_full_impo 


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize Lambda block (random, factorized) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !Generate random infinite Lambda block
 SUBROUTINE initialize_rand_lambda(lambda_block, chi, N_sites_in)

  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)        :: lambda_block(:)
  INTEGER,                         INTENT(IN)           :: chi
  INTEGER,                         INTENT(IN), OPTIONAL :: N_sites_in

  !temp site tensor
  COMPLEX(KIND=DP), ALLOCATABLE :: tempSite(:)

  !dims & indices
  INTEGER :: i
  INTEGER :: N_sites, site

  !rand seed
  LOGICAL, PARAMETER :: use_const_seed = .TRUE. 
  REAL               :: cpu_times 
  INTEGER            :: seed 

  !set optional args
  N_sites = OptArg(N_sites_in, 2)

  !Generate rand nums using a const seed or a seed from cpu_time
  IF(use_const_seed) THEN
     seed = 86456
  ELSE
     CALL cpu_time(cpu_times)
     seed = NINT(cpu_times)
  end IF
  CALL srand(seed)

  !Init random site
  ALLOCATE(tempSite(chi)) 
  Do i=1,chi
     tempSite(i) = CMPLX(rand(), 0d0)
  end Do

  !allocate empty imps
  CALL allocate_empty_lambda_block(lambda_block, N_sites)

  !generate random sites
  DO site=1,N_sites
     CALL allocate_lambda_site(lambda_block(site), tempSite)
  END DO

 END SUBROUTINE initialize_rand_lambda






 !Initialize factorized lambda_block
 SUBROUTINE initialize_factorized_lambda(lambda_block, N_sites, chi_In)

  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)        :: lambda_block(:)
  INTEGER,                         INTENT(IN)           :: N_sites
  INTEGER,                         INTENT(IN), OPTIONAL :: chi_In

  !dims && indices
  INTEGER :: chi, site

  !bond dim of a factorized (product) peps
  chi = OptArg(chi_In, 1)

  !Allocate product state gammas
  CALL allocate_empty_lambda_block(lambda_block, N_sites)

  DO site=1,N_sites
     CALL allocate_lambda_site(lambda_block(site), chi)
     lambda_block(site) % m = (1.0d0, 0.0d0)
  end DO

 END SUBROUTINE initialize_factorized_lambda

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! initialize iPEPS block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !Initialize factorized peps_block
 SUBROUTINE initialize_factorized_peps(peps_block, LocalDim, N_sites, vecA_In, chi_In)

  TYPE(block_peps), ALLOCATABLE, INTENT(INOUT)        :: peps_block(:)
  INTEGER,                       INTENT(IN)           :: LocalDim, N_sites
  COMPLEX(KIND=DP),              INTENT(IN), OPTIONAL :: vecA_In(:)
  INTEGER,                       INTENT(IN), OPTIONAL :: chi_In

  !vecA && bond dim for initializing a product state
  COMPLEX(KIND=DP), ALLOCATABLE :: vecA(:)
  INTEGER :: chi, site

  !bond dim of a factorized (product) peps
  chi = OptArg(chi_In, 1)

  !find vecA
  IF(PRESENT(vecA_in)) THEN
     CALL copyTens(vecA, vecA_In)
  ELSE
     vecA(:) = (0.0d0, 0.0d0)
     vecA(1) = (1.0d0, 0.0d0)
  end IF

  !Allocate product state gammas
  CALL allocate_empty_peps_block(peps_block, N_sites)

  DO site=1,N_sites
     CALL allocate_peps_site(peps_block(site),  (/LocalDim, chi, chi, chi, chi/))
     peps_block(site) % m(:,1,1,1,1) = vecA
  end DO

 END SUBROUTINE initialize_factorized_peps





 !Initialize random peps_block
 SUBROUTINE initialize_rand_peps(peps_block, LocalDim, SNdim, EWdim, N_sites_in)

  TYPE(block_peps), ALLOCATABLE, INTENT(INOUT)        :: peps_block(:)
  INTEGER,                       INTENT(IN)           :: LocalDim, SNdim, EWdim
  INTEGER,                       INTENT(IN), OPTIONAL :: N_sites_in

  !temp site tensor
  COMPLEX(KIND=DP), ALLOCATABLE :: tempSite(:,:,:,:,:)

  !dims & indices
  INTEGER :: i, j, k, l, r
  INTEGER :: N_sites, site

  !rand seed
  LOGICAL, PARAMETER :: use_const_seed = .TRUE. 
  REAL               :: cpu_times 
  INTEGER            :: seed 

  !set optional args
  N_sites = OptArg(N_sites_in, 2)

  !Generate rand nums using a const seed or a seed from cpu_time
  IF(use_const_seed) THEN
     seed = 86456
  ELSE
     CALL cpu_time(cpu_times)
     seed = NINT(cpu_times)
  end IF
  CALL srand(seed)

  !Init random site
  ALLOCATE(tempSite(LocalDim, SNdim, SNdim, EWdim, EWdim))

  Do i=1,LocalDim
    Do j=1,SNdim
      Do k=1,SNdim
        Do l=1,EWdim
          Do r=1,EWdim
             tempSite(i,j,k,l,r) = CMPLX(rand(), 0d0)
          end Do
        end Do
      end Do
    end Do
  end Do

  !Allocate empty iPEPS
  CALL allocate_empty_peps_block(peps_block, N_sites)

  !Generate random iPEPS sites
  DO site=1,N_sites
     CALL allocate_peps_site(peps_block(site), tempSite)
  END DO

 END SUBROUTINE initialize_rand_peps








 !!! Initialize iPEPS evec to the last converged evec 
 !!! (with rounded-off values) 
 SUBROUTINE roundoff_ipeps(iGamma, iLambda, eps)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  REAL(KIND=DP),                   INTENT(IN)    :: eps

  !! Local iMPS copy
  TYPE(block_peps),   ALLOCATABLE :: iGammaCpy(:)
  TYPE(block_lambda), ALLOCATABLE :: iLambdaCpy(:)

  COMPLEX(KIND=DP), ALLOCATABLE :: tempGamma(:,:,:,:,:), tempLambda(:)
  INTEGER :: chi(4), minChi
  INTEGER :: site, bond

  !! (1A) Round off iGamma
  DO site=1,2
     CALL roundTens(iGamma(site) % m, eps)
  end DO

  !! (1B) Round off iLambda, get the new bond dim (number of nonzero lambdas)
  DO bond=1,4
     CALL roundTens(iLambda(bond) % m, eps)
     chi(bond) = COUNT(ABS(iLambda(bond) % m) .GT. eps)
  end DO

  !! (2) Use the smaller one of two bond dims at sites=1,2 
  !!     (Must guarantee that sites=1,2 have equal chi, and no zero lambdas)
  minChi = MIN(chi(1), chi(2), chi(3), chi(4)) 

  !!!!!! (3) Copy rounded tensors to iPEPS, ensuring that bond dims match !!!!!!

  !!  (3A) Cut dims of G,L to exclude any zero lambdas 
  !!       (i.e. to avoid non-invertible boundary matrices)
  DO site=1,2
     CALL copyTens(tempGamma, iGamma(site) % m,  (/iGamma(site) % LocalDim, minChi, minChi, minChi, minChi/))
     CALL allocate_peps_site(iGamma(site),  tempGamma)
  end DO

  !!  (3B) Cut dims of G,L to exclude any zero lambdas 
  !!       (i.e. to avoid non-invertible boundary matrices)
  DO bond=1,4
     CALL copyTens(tempLambda, iLambda(bond) % m,  (/minChi/))
     CALL allocate_lambda_site(iLambda(bond),      tempLambda)
  end DO

 END SUBROUTINE roundoff_ipeps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize complete iPEPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE initialize_rand_full_ipeps(iGamma, iLambda, localDim, chi)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  INTEGER,                         INTENT(IN)    :: localDim, chi

  !! Rand tensors
  COMPLEX(KIND=DP), ALLOCATABLE :: tmpLAM(:), tmpGAM(:,:,:,:,:)

  !! Num of sites
  INTEGER, PARAMETER :: N_sites = 2
  INTEGER, PARAMETER :: N_bonds = 4
  INTEGER            :: site, bond

  !! Norm fac
  COMPLEX(KIND=DP) :: fac

  !! Init random site tensors
  ALLOCATE(tmpLAM(chi), tmpGAM(localDim, chi, chi, chi, chi))
  CALL initialize_rand_tens(tmpGAM)
  CALL initialize_rand_tens(tmpLAM)

  !! Normalize lambda
  fac = SQRT(SUM(tmpLAM**2))
  tmpLAM = tmpLAM / fac
  tmpGAM = tmpGAM * SQRT(fac)**4

  !! Sort lambdas
  CALL sort_array(tmpLAM)

  !! Setup random iGamma 
  CALL allocate_empty_peps_block(iGamma, N_sites)
  DO site=1,N_sites
     CALL allocate_peps_site(iGamma(site), tmpGAM)
  END DO

  !! Setup random iLambda
  CALL allocate_empty_lambda_block(iLambda, N_bonds)
  DO bond=1,N_bonds
     CALL allocate_lambda_site(iLambda(bond), tmpLAM)
  END DO

 END SUBROUTINE initialize_rand_full_ipeps 


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BVEC: Initialize imps-imps boundary vector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Initialize bvec !!!
 SUBROUTINE initialize_bvec(bvec, Tmpo, DIR)

  COMPLEX(KIND=DP),   ALLOCATABLE, INTENT(INOUT)        :: bvec(:,:)
  TYPE(transfer_mpo),              INTENT(IN)           :: Tmpo(:)
  CHARACTER(LEN=*),                INTENT(IN)           :: DIR

  LOGICAL, PARAMETER :: init_eye_bvec = .TRUE.

  IF(init_eye_bvec) THEN

     !! Initialize eye bvec
     CALL initialize_eye_bvec(bvec,  Tmpo, DIR)
  ELSE

     !! Initialize random bvec
     CALL initialize_rand_bvec(bvec, Tmpo, DIR)
  end IF

 END SUBROUTINE initialize_bvec



 

 !!! Initialize bvec = Identity matrix !!!
 SUBROUTINE initialize_eye_bvec(bvec, Tmpo, DIR)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: bvec(:,:)
  TYPE(transfer_mpo),            INTENT(IN)    :: Tmpo(:)
  CHARACTER(LEN=*),              INTENT(IN)    :: DIR

  INTEGER :: dimB(2), minDimB

  !! Get bvec dims
  dimB = bvec_dims(Tmpo, DIR) 

  !! The smaller dim of dimB(1,2)
  minDimB = MIN(dimB(1), dimB(2)) 

  !! Construct Eye bvec
  CALL copyTens(bvec, matEye(minDimB), dimB)
  bvec = bvec/SUM(bvec**2)

 END SUBROUTINE initialize_eye_bvec





 !!! Initialize a random bvec for power method iteration !!!
 SUBROUTINE initialize_rand_bvec(bvec, Tmpo, DIR)

  COMPLEX(KIND=DP),   ALLOCATABLE, INTENT(INOUT) :: bvec(:,:)
  TYPE(transfer_mpo),              INTENT(IN)    :: Tmpo(:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! Dims & indices
  INTEGER :: i, j, site 
  INTEGER :: dimB(2), minDim

  !! Rand seed
  LOGICAL, PARAMETER :: use_const_seed = .TRUE. 
  REAL               :: cpu_times 
  INTEGER            :: seed 

  !! Generate rand nums using a const seed or a seed from cpu_time
  IF(use_const_seed) THEN
     seed = 86456
  ELSE
     CALL cpu_time(cpu_times)
     seed = NINT(cpu_times)
  end IF
  CALL srand(seed)

  !! Get dims of new bvec (= dims of Tmpo)
  dimB    = bvec_dims(Tmpo, DIR)
  minDim  = MIN(dimB(1), dimB(2)) 

  !! Initialize bvec = random vector
  ALLOCATE(bvec(dimB(1), dimB(2)))
  Do i=1,dimB(1)
    Do j=1,dimB(2)
       bvec(i,j) = CMPLX(rand(), 0d0)
     end Do
  end Do

 END SUBROUTINE initialize_rand_bvec





 !!! Find dims of bvec based on Tmpo !!!
 FUNCTION bvec_dims(Tmpo, DIR)

  TYPE(transfer_mpo), INTENT(IN) :: Tmpo(:)
  CHARACTER(LEN=*),   INTENT(IN) :: DIR
  INTEGER                        :: bvec_dims(2) 

  !! Size of Tmpo
  INTEGER :: N_sites

  N_sites = SIZE(Tmpo)

  !! Get dims of Tmpo
  SELECT CASE(DIR)
  CASE('R')
     !Right evec
     bvec_dims(1) = SIZE(Tmpo(N_sites) % D, 3)  
     bvec_dims(2) = SIZE(Tmpo(N_sites) % U, 3)    
  CASE('L')
     !Left evec
     bvec_dims(1) = SIZE(Tmpo(1) % D, 2)  
     bvec_dims(2) = SIZE(Tmpo(1) % U, 2)
  CASE DEFAULT
     CALL invalid_flag("bvec_dims -- invalid DIR ", DIR)
  end SELECT

 END FUNCTION bvec_dims

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INITIALIZE RANDOM TENSOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Initialize random 1D tensor !!!
 SUBROUTINE initialize_rand_tens_1D(tensA)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:) 

  !Dims & indices
  INTEGER :: i, dimA(1)

  !Rand seed
  LOGICAL, PARAMETER :: use_const_seed = .TRUE.
  REAL    :: cpu_times 
  INTEGER :: seed 

  !Generate rand nums using a const seed or a seed from cpu_time
  IF(use_const_seed) THEN
     seed = 86456
  ELSE
     CALL cpu_time(cpu_times)
     seed = NINT(cpu_times)
  end IF
  CALL srand(seed)

  !Get mpo dims
  dimA = shape(tensA)

  !Initialize random 3D tensor
  Do i=1,dimA(1)  
     tensA(i) = CMPLX(rand(), 0d0)
  end Do

 END SUBROUTINE initialize_rand_tens_1D







 !!! Initialize random 2D tensor !!!
 SUBROUTINE initialize_rand_tens_2D(tensA)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:) 

  !Dims & indices
  INTEGER :: i, j
  INTEGER :: dimA(2)

  !Rand seed
  LOGICAL, PARAMETER :: use_const_seed = .TRUE.
  REAL    :: cpu_times 
  INTEGER :: seed 

  !Generate rand nums using a const seed or a seed from cpu_time
  IF(use_const_seed) THEN
     seed = 86456
  ELSE
     CALL cpu_time(cpu_times)
     seed = NINT(cpu_times)
  end IF
  CALL srand(seed)

  !Get mpo dims
  dimA = shape(tensA)

  !Initialize random 3D tensor
  Do i=1,dimA(1)
    Do j=1,dimA(2)  
       tensA(i,j) = CMPLX(rand(), 0d0)
    end Do
  end Do

 END SUBROUTINE initialize_rand_tens_2D







 !!! Initialize random 3D tensor !!!
 SUBROUTINE initialize_rand_tens_3D(tensA)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:) 

  !Dims & indices
  INTEGER :: i, j, k
  INTEGER :: dimA(3)

  !Rand seed
  LOGICAL, PARAMETER :: use_const_seed = .TRUE.
  REAL    :: cpu_times 
  INTEGER :: seed 

  !Generate rand nums using a const seed or a seed from cpu_time
  IF(use_const_seed) THEN
     seed = 86456
  ELSE
     CALL cpu_time(cpu_times)
     seed = NINT(cpu_times)
  end IF
  CALL srand(seed)

  !Get mpo dims
  dimA = shape(tensA)

  !Initialize random 3D tensor
  Do i=1,dimA(1)
    Do j=1,dimA(2)
      Do k=1,dimA(3)    
         tensA(i,j,k) = CMPLX(rand(), 0d0)
      end Do
    end Do
  end Do

 END SUBROUTINE initialize_rand_tens_3D








 !!! Initialize random 4D tensor !!!
 SUBROUTINE initialize_rand_tens_4D(tensA)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:)

  !Dims & indices
  INTEGER :: i, j, k, l
  INTEGER :: dimA(4)

  !Rand seed
  LOGICAL, PARAMETER :: use_const_seed = .TRUE.
  REAL :: cpu_times 
  INTEGER :: seed 

  !Generate rand nums using a const seed or a seed from cpu_time
  IF(use_const_seed) THEN
     seed = 86456
  ELSE
     CALL cpu_time(cpu_times)
     seed = NINT(cpu_times)
  end IF
  CALL srand(seed)

  !Get mpo dims
  dimA = shape(tensA)

  !Initialize random 4D tensor
  Do i=1,dimA(1)
    Do j=1,dimA(2)
      Do k=1,dimA(3)
        Do l=1,dimA(4)       
           tensA(i,j,k,l) = CMPLX(rand(), 0d0)
        end Do
      end Do
    end Do
  end Do

 END SUBROUTINE initialize_rand_tens_4D






 !!! Initialize random 5D tensor !!!
 SUBROUTINE initialize_rand_tens_5D(tensA)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:,:)

  !Dims & indices
  INTEGER :: i, j, k, l, m
  INTEGER :: dimA(5)

  !Rand seed
  LOGICAL, PARAMETER :: use_const_seed = .TRUE.
  REAL    :: cpu_times 
  INTEGER :: seed 

  !Generate rand nums using a const seed or a seed from cpu_time
  IF(use_const_seed) THEN
     seed = 86456
  ELSE
     CALL cpu_time(cpu_times)
     seed = NINT(cpu_times)
  end IF
  CALL srand(seed)

  !Get mpo dims
  dimA = shape(tensA)

  !Initialize random 5D tensor
  Do i=1,dimA(1)
    Do j=1,dimA(2)
      Do k=1,dimA(3)
        Do l=1,dimA(4)   
          Do m=1,dimA(5)     
             tensA(i,j,k,l,m) = CMPLX(rand(), 0d0)
          end Do
        end Do
      end Do
    end Do
  end Do

 END SUBROUTINE initialize_rand_tens_5D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE iteration_helper
