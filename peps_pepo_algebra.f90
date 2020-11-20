MODULE peps_pepo_algebra

 USE utility
 USE definitions_mps_mpo
 USE definitions_peps_pepo
 USE mps_mpo_utility
 USE mps_mpo_algebra_finite
 USE mps_mpo_algebra_inf
 USE project_ipeps
 USE project_ipepo

 IMPLICIT NONE

 interface resize_peps_bonds
    module procedure resize_peps_GamLam
    module procedure resize_peps_PLAIN
 end interface resize_peps_bonds

 private resize_peps_GamLam, resize_peps_PLAIN

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (1) Multiply iPEPS by two-body MPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE mult_ipeps_2body_mpo(iGamma, iLambda, mpo, chi, eps, DIR)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT)        :: iGamma(:)           !iGamma of PEPS
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)        :: iLambda(:)          !iLambda of PEPS
  TYPE(block_mpo),                 INTENT(IN)           :: mpo(:)              !MPO Propagator
  INTEGER,                         INTENT(IN)           :: chi                 !SVD chi
  REAL(KIND=DP),                   INTENT(IN)           :: eps                 !SVD eps
  CHARACTER(LEN=*),                INTENT(IN)           :: DIR                 !Bond we're propagating
  
  !! iMPS projection 
  TYPE(block_mps),    ALLOCATABLE :: iGmps(:)
  TYPE(block_lambda), ALLOCATABLE :: iLmps(:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: Xa(:,:), Yb(:,:)
  
  !! (1) Project iPEPS onto iMPS
  !CALL project_ipeps_to_imps(iGmps, iLmps, iGamma, iLambda, DIR)
  CALL calc_factorized_ipeps(iGmps, iLmps, Xa, Yb, iGamma, iLambda, DIR)

  !! (2) Absorb central lambda 
  CALL LambdaMUL(iGmps(1) % m, iLmps(1) % m, iGmps(2) % m)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! (3) Multiply mpo_prop and iGmps 
  CALL simple_mult_mps_mpo(iGmps, mpo) 

  !! (4) Compute SVD of central iGmps bond
  CALL compute_svd_of_imps_bond(iGmps, iLmps(1), chi, eps, '12')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! (5) Restore iPEPS from iMPS projection
  !CALL unproject_ipeps_from_imps(iGamma, iLambda, iGmps, iLmps, DIR)
  CALL restore_unfactorized_ipeps(iGamma, iLambda, iGmps, iLmps, Xa, Yb, DIR)

 END SUBROUTINE mult_ipeps_2body_mpo





 !!! Mult iPEPS by one-site propagator !!!
 SUBROUTINE mult_ipeps_onsite_prop(iGamma, iLambda, prop, chi, eps)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGamma(:)           !iGamma of PEPS
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)          !iLambda of PEPS
  COMPLEX(KIND=DP),                INTENT(IN)    :: prop(:,:)           !Propagator
  INTEGER,                         INTENT(IN)    :: chi                 !SVD chi
  REAL(KIND=DP),                   INTENT(IN)    :: eps                 !SVD eps

  !! Multiply iPEPS by one-site propagator
  iGamma(1) % m = TENSMUL(iGamma(1) % m,  prop,  MULT='12')
  iGamma(2) % m = TENSMUL(iGamma(2) % m,  prop,  MULT='12')

 END SUBROUTINE mult_ipeps_onsite_prop




 !!! Mult iPEPS by dense propagator !!!
 SUBROUTINE mult_ipeps_dense_prop(iGamma, iLambda, prop_dense, chi, eps, DIR)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT)        :: iGamma(:)           !iGamma of PEPS
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)        :: iLambda(:)          !iLambda of PEPS
  COMPLEX(KIND=DP),                INTENT(IN)           :: prop_dense(:,:)     !Dense Propagator
  INTEGER,                         INTENT(IN)           :: chi                 !SVD chi
  REAL(KIND=DP),                   INTENT(IN)           :: eps                 !SVD eps
  CHARACTER(LEN=*),                INTENT(IN)           :: DIR                 !Bond we're propagating
  
  !! iMPS projection 
  TYPE(block_mps),    ALLOCATABLE :: iGmps(:)
  TYPE(block_lambda), ALLOCATABLE :: iLmps(:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: Xa(:,:), Yb(:,:)
  
  !! (1) Project iPEPS onto iMPS
  !CALL project_ipeps_to_imps(iGmps, iLmps, iGamma, iLambda, DIR)
  CALL calc_factorized_ipeps(iGmps, iLmps, Xa, Yb, iGamma, iLambda, DIR)

  !! (2) Absorb central lambda 
  CALL LambdaMUL(iGmps(1) % m, iLmps(1) % m, iGmps(2) % m)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! (3) Multiply dense propagator and iGmps, compute SVD
  CALL mult_imps_dense_prop(iGmps, iLmps(1), prop_dense, chi, eps, '12')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! (4) Restore iPEPS from iMPS projection
  !CALL unproject_ipeps_from_imps(iGamma, iLambda, iGmps, iLmps, DIR)
  CALL restore_unfactorized_ipeps(iGamma, iLambda, iGmps, iLmps, Xa, Yb, DIR)

 END SUBROUTINE mult_ipeps_dense_prop


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (2) Compute SVD of iPEPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute SVD of all iPEPS bonds !!!
 SUBROUTINE compute_svd_of_ipeps(iGamma, iLambda, chi, eps)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT)        :: iGamma(:)  !iGamma of PEPS
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)        :: iLambda(:) !iLambda of PEPS
  INTEGER,                         INTENT(IN)           :: chi        !SVD chi  
  REAL(KIND=DP),                   INTENT(IN)           :: eps        !SVD eps

  CALL compute_svd_of_ipeps_bond(iGamma, iLambda, chi, eps, DIR='S')
  CALL compute_svd_of_ipeps_bond(iGamma, iLambda, chi, eps, DIR='W')
  CALL compute_svd_of_ipeps_bond(iGamma, iLambda, chi, eps, DIR='N')
  CALL compute_svd_of_ipeps_bond(iGamma, iLambda, chi, eps, DIR='E')

 END SUBROUTINE compute_svd_of_ipeps




 !!! Compute SVD of a single iPEPS bond !!!
 SUBROUTINE compute_svd_of_ipeps_bond(iGamma, iLambda, chi, eps, DIR)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGamma(:)           !iGamma of PEPS
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)          !iLambda of PEPS
  INTEGER,                         INTENT(IN)    :: chi                 !SVD chi
  REAL(KIND=DP),                   INTENT(IN)    :: eps                 !SVD eps
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR                 !Bond we're propagating

  !! iMPS projection 
  TYPE(block_mps),    ALLOCATABLE :: iGmps(:)
  TYPE(block_lambda), ALLOCATABLE :: iLmps(:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: Xa(:,:), Yb(:,:)

  !! (1) Project iPEPS onto iMPS
  !CALL project_ipeps_to_imps(iGmps, iLmps, iGamma, iLambda, DIR)
  CALL calc_factorized_ipeps(iGmps, iLmps, Xa, Yb, iGamma, iLambda, DIR)

  !! (2) Absorb central lambda
  CALL LambdaMUL(iGmps(1) % m, iLmps(1) % m, iGmps(2) % m)

  !! (3) Compute SVD of central iGmps bond
  CALL compute_svd_of_imps_bond(iGmps, iLmps(1), chi, eps, '12')

  !! (4) Restore iPEPS from iMPS projection
  !CALL unproject_ipeps_from_imps(iGamma, iLambda, iGmps, iLmps, DIR)
  CALL restore_unfactorized_ipeps(iGamma, iLambda, iGmps, iLmps, Xa, Yb, DIR)

 END SUBROUTINE compute_svd_of_ipeps_bond

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (3) RESHAPE PEPS <-> PEPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Reshape PEPS to PEPO !!!
 SUBROUTINE reshape_PEPS_to_PEPO(pepo, peps)

  TYPE(block_pepo), ALLOCATABLE, INTENT(INOUT) :: pepo(:)
  TYPE(block_peps),              INTENT(IN)    :: peps(:)

  INTEGER :: site, N_sites
  INTEGER :: localDim
  
  !! Block size
  N_sites = SIZE(peps)

  !! PEPS local dim
  localDim = SQRTint(peps(1) % LocalDim)

  !! Allocate empty PEPO
  CALL allocate_empty_pepo_block(pepo, N_sites)

  !! Calc PEPO sites
  DO site=1,N_sites
     CALL allocate_pepo_site(pepo(site),  RESHAPE_6D(peps(site) % m, '12,3,4,5,6', (/ localDim, localDim /)))
  end DO

 END SUBROUTINE reshape_PEPS_to_PEPO





 !!! Reshape PEPO to PEPS !!!
 SUBROUTINE reshape_PEPO_to_PEPS(peps, pepo)

  TYPE(block_peps), ALLOCATABLE, INTENT(INOUT) :: peps(:)
  TYPE(block_pepo),              INTENT(IN)    :: pepo(:)

  INTEGER :: site, N_sites
  
  !! Block size
  N_sites = SIZE(pepo)

  !! Allocate empty PEPS
  CALL allocate_empty_peps_block(peps, N_sites)

  !! Calc PEPS sites
  DO site=1,N_sites
     CALL allocate_peps_site(peps(site),  RESHAPE_5D(pepo(site) % m, '12,3,4,5,6'))
  end DO

 END SUBROUTINE reshape_PEPO_to_PEPS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (4) Rescaling/normalizing PEPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Resize the bonds of a Gamma-Lambda iPEPS !!!
 SUBROUTINE resize_peps_GamLam(iGamma, iLambda, chi)

  TYPE(block_peps),   INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), INTENT(INOUT) :: iLambda(:)
  INTEGER,            INTENT(IN)    :: chi

  COMPLEX(KIND=DP), ALLOCATABLE :: tmpGam(:,:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tmpLam(:)

  INTEGER :: N_sites, N_bonds
  INTEGER :: site, bond

  N_sites = SIZE(iGamma)
  N_bonds = SIZE(iLambda)

  !! Resize Gammas
  DO site=1,N_sites
     CALL copyTens(tmpGam, iGamma(site) % m, (/iGamma(site) % LocalDim, chi, chi, chi, chi/))
     CALL allocate_peps_site(iGamma(site), tmpGam)
  end DO

  !! Resize Lambdas
  DO bond=1,N_bonds
     CALL copyTens(tmpLam, iLambda(bond) % m, (/chi/))
     CALL allocate_lambda_site(iLambda(bond), tmpLam)
  end DO

 END SUBROUTINE resize_peps_GamLam



 !!! Resize the bonds of a plain iPEPS !!!
 SUBROUTINE resize_peps_PLAIN(ipeps, chi)

  TYPE(block_peps),   INTENT(INOUT) :: ipeps(:)
  INTEGER,            INTENT(IN)    :: chi

  COMPLEX(KIND=DP), ALLOCATABLE :: tmpPeps(:,:,:,:,:)

  INTEGER :: N_sites, site

  N_sites = SIZE(ipeps)

  !! Resize iPEPS sites
  DO site=1,N_sites
     CALL copyTens(tmpPeps, ipeps(site) % m, (/ipeps(site) % LocalDim, chi, chi, chi, chi/))
     CALL allocate_peps_site(ipeps(site), tmpPeps)
  end DO

 END SUBROUTINE resize_peps_PLAIN







 !!! Multiply PEPS by a given value !!!
 SUBROUTINE mult_peps_by_value(peps_block, input_value)
 
  TYPE(block_peps), INTENT(INOUT) :: peps_block(:)
  COMPLEX(KIND=DP), INTENT(IN)    :: input_value

  !! Dims & indices
  INTEGER :: N_sites, site

  !! Size of mps_block
  N_sites = SIZE(peps_block)
  
  !! Multiply PEPS by value site-by-site
  Do site=1,N_sites
     peps_block(site) % m = input_value**(1.0/REAL(N_sites, KIND=DP)) * peps_block(site) % m
  end Do

 END SUBROUTINE mult_peps_by_value




 !!! Rescale PEPS block my max ABS element !!!
 SUBROUTINE rescale_peps_by_max_element(peps)
 
  TYPE(block_peps), INTENT(INOUT) :: peps(:)

  !! Dims & indices
  INTEGER :: N_sites, site

  !! Maximum abs element in mps_block
  REAL(KIND=DP), ALLOCATABLE :: maxEls(:)
  REAL(KIND=DP)              :: maxEl

  !! Size of mps_block
  N_sites = SIZE(peps)

  !! (1) Find max element on each peps site
  ALLOCATE(maxEls(N_sites))
  Do site=1,N_sites
     maxEls(site) = MAXVAL(ABS(peps(site) % m))
  end Do

  !! (2) Find max element in peps_block
  maxEl = MAXVAL(maxEls)
   
  !! (3) Rescale by [max element] -- i.e. make max element = 1
  CALL mult_peps_by_value(peps, 1.0d0/complx(maxEl))

 END SUBROUTINE rescale_peps_by_max_element




 SUBROUTINE equalize_ipeps_lambdas(iGamma, iLambda)

  TYPE(block_peps),   INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), INTENT(INOUT) :: iLambda(:)

  COMPLEX(KIND=DP)           :: fac
  REAL(KIND=DP), ALLOCATABLE :: maxlam(:)
  INTEGER                    :: site, N_sites

  N_sites = SIZE(iLambda)

  !! Get max lambda on each site
  ALLOCATE(maxlam(N_sites))
  DO site=1,N_sites
     maxlam(site) = MAXVAL(ABS(iLambda(site) % m)) 
  end DO

  !! Rescale lambdas to have the same max, normalize them so that Sum(iLambda**2) = 1
  DO site=1,N_sites
     iLambda(site) % m = iLambda(site) % m * PRODUCT(maxlam)**(0.25D0) / maxlam(site)
     !fac               = sqrt(SUM(iLambda(site) % m ** 2))
     !iLambda(site) % m = iLambda(site) % m / fac
     !iGamma(1) % m     = iGamma(1) % m * sqrt(fac)
     !iGamma(2) % m     = iGamma(2) % m * sqrt(fac)
  end DO

 END SUBROUTINE equalize_ipeps_lambdas






 !!! When PEPS is meant to represent, 
 !!! e.g. a vectorized Hermitian operator like rho, purify its Hermitian part 
 SUBROUTINE hermitize_peps(peps)

  TYPE(block_peps), INTENT(INOUT) :: peps(:)

  INTEGER :: N_sites, site

  N_sites = SIZE(peps)

  !! Hermitize each site of PEPS
  DO site=1,N_sites
     CALL allocate_peps_site(peps(site), 0.5D0*(peps(site) % m + CONJG(peps(site) % m)))
  end DO

 END SUBROUTINE hermitize_peps




 !!! Symmetrize iPEPS so that both sites have the same MAX element !!!
 SUBROUTINE symmetrize_ipeps_sites(ipeps)
 
  TYPE(block_peps), INTENT(INOUT) :: ipeps(:)

  !! Dims & indices
  INTEGER :: site

  !! Maximum abs element in mps_block
  REAL(KIND=DP) :: maxEls(2), fac

  !! Find max elements on each ipeps site
  Do site=1,2
     maxEls(site) = MAXVAL(ABS(ipeps(site) % m))
  end Do

  !! Rescale ipeps so that both sites have the same MAX element
  fac = sqrt(ABS(maxEls(2)/maxEls(1)))

  ipeps(1) % m = ipeps(1) % m * fac
  ipeps(2) % m = ipeps(2) % m / fac

 END SUBROUTINE symmetrize_ipeps_sites




 !!! Balance MPS bond using the proposal in Lubasch et al !!!
 SUBROUTINE balance_ipeps_all_bonds(iGamma, iLambda)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGamma(:)  
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:) 

  CALL balance_ipeps_bond(iGamma, iLambda, 'E')
  CALL balance_ipeps_bond(iGamma, iLambda, 'N')
  CALL balance_ipeps_bond(iGamma, iLambda, 'W')
  CALL balance_ipeps_bond(iGamma, iLambda, 'S')

 END SUBROUTINE balance_ipeps_all_bonds






 !!! Balance MPS bond using the proposal in Lubasch et al !!!
 SUBROUTINE balance_ipeps_bond(iGamma, iLambda, DIR)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGamma(:)  
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:) 
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR        !Bond to balance

  !! iMPS projection 
  TYPE(block_mps),    ALLOCATABLE :: iGmps(:)
  TYPE(block_lambda), ALLOCATABLE :: iLmps(:)

  !! Project iPEPS onto iMPS
  CALL project_ipeps_to_imps(iGmps, iLmps, iGamma, iLambda, DIR)

  !! Balance bond specified by DIR
  CALL balance_imps_bond(iGmps, iLmps(1)) 

  !! Project back onto iPEPS
  CALL unproject_ipeps_from_imps(iGamma, iLambda, iGmps, iLmps, DIR)

 END SUBROUTINE balance_ipeps_bond






 !!! Round off iPEPS tensors to reduce numerical noise !!!
 SUBROUTINE round_ipeps_sites(iGamma, iLambda, cutoff0)

  TYPE(block_peps),   INTENT(INOUT)        :: iGamma(:)
  TYPE(block_lambda), INTENT(INOUT)        :: iLambda(:)
  REAL(KIND=DP),      INTENT(IN), OPTIONAL :: cutoff0

  !! Indices
  INTEGER :: N_sites, site
  INTEGER :: N_bonds, bond

  !! Cutoff factor
  REAL(KIND=DP) :: cutoff

  !! Get num of iPEPS sites && bonds
  N_sites = SIZE(iGamma)
  N_bonds = SIZE(iLambda)

  !! Set cutoff factor
  cutoff = OptArg(cutoff0, 1.0d-08)

  !! Round/Filter iGamma
  DO site=1,N_sites
     CALL roundTens(iGamma(site) % m,  cutoff)
  end DO

  !! Round/Filter iLambda
  DO bond=1,N_bonds
     CALL roundTens(iLambda(bond) % m, cutoff)
  end DO
  
 END SUBROUTINE round_ipeps_sites




 !!! Compute trace of iPEPS assuming MF environment !!!
 SUBROUTINE calc_simplified_trace(C_norm, iGamma, iLambda, reduced_ipeps_site)

  COMPLEX(KIND=DP),   INTENT(OUT) :: C_norm
  TYPE(block_peps),   INTENT(IN)  :: iGamma(:)
  TYPE(block_lambda), INTENT(IN)  :: iLambda(:)

  INTERFACE 
    FUNCTION reduced_ipeps_site(ipeps, op1, op2)
      USE utility
      USE definitions_mps_mpo
      USE definitions_peps_pepo
      TYPE(block_peps), INTENT(IN)           :: ipeps                       !ipeps site 
      COMPLEX(KIND=DP), INTENT(IN), OPTIONAL :: op1(:,:), op2(:,:)          !operators applied to the site (we may apply two operators)
      COMPLEX(KIND=DP), ALLOCATABLE          :: reduced_ipeps_site(:,:,:,:) !result
    end FUNCTION reduced_ipeps_site
  end INTERFACE

  !! Nolambda iPEPS && its sites traced out
  TYPE(block_peps), ALLOCATABLE :: ipeps(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: ipepsA(:,:,:,:),  ipepsB(:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: ipepsAB(:,:,:,:), ipepsBA(:,:,:,:), ipepsC(:,:,:,:)

  !! Construct nolambda iPEPS (absorb all the lambdas)
  CALL absorb_ipeps_lambdas(ipeps, iGamma, iLambda)

  !! Trace out iPEPS sites: multiply 2x2 blocks together
  ipepsA = reduced_ipeps_site(ipeps(1))
  ipepsB = reduced_ipeps_site(ipeps(2))

  ipepsAB = TENSMUL(ipepsA,  ipepsB,  MULT='43', FUSE='(11,22)')
  ipepsBA = TENSMUL(ipepsB,  ipepsA,  MULT='43', FUSE='(11,22)')
  ipepsC  = TENSMUL(ipepsAB, ipepsBA, MULT='21', FUSE='(33,44)')

  !! Trace out over bond dims assuming MF environment 
  C_norm = TRACE(TRACE(ipepsC, '12'))

 END SUBROUTINE calc_simplified_trace

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 !!! Quasi-orthogonalize iPEPS: apply eye gate to all bonds until convergence !!!
 SUBROUTINE quasiorthogonalize_ipeps(iGamma, iLambda, chi, eps)

  USE datafile_utility, ONLY: setupConvDatafile

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGamma(:)   !iGamma of iPEPS
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)  !iLambda of iPEPS
  INTEGER,                         INTENT(IN)    :: chi         !SVD chi
  REAL(KIND=DP),                   INTENT(IN)    :: eps         !SVD eps

  !! Old iLambda for lambda convergence testing
  TYPE(block_lambda), ALLOCATABLE :: iLambdaOLD(:)

  !! Iteration params
  REAL(KIND=DP), PARAMETER :: eta     = 1.0D-04
  INTEGER,       PARAMETER :: N_steps = 20
  LOGICAL                  :: is_converged

  !! Indices
  INTEGER :: step
  INTEGER :: datafile

  !! Setup Lambda convergence file, initialize old Lambda
  CALL copy_lambda_block(iLambdaOLD, iLambda)
  CALL setupConvDatafile(datafile=datafile, descr="iPEPS_lambda_quasiorth", errstr="Err_lambda", valstr="Max_lambda")

  !! Initialize convergence tester
  is_converged = .FALSE.

  !! Iterate until lambdas converge
  DO step=1,N_steps
  
     !! Do SVD on all iPEPS bonds (single step)
     CALL compute_svd_of_ipeps(iGamma, iLambda, chi, eps)

     !! Test iLambda convergence
     CALL test_gen_lambda_convergence(is_converged, iLambda, iLambdaOLD, eta, datafile, step)
     IF(is_converged) EXIT
  end DO

  WRITE(datafile, *)
  WRITE(datafile, *)  
  CLOSE(datafile)

 END SUBROUTINE quasiorthogonalize_ipeps




 !!! Test if Lambdas have converged (general case, not only TEBD) !!!
 SUBROUTINE test_gen_lambda_convergence(is_converged, newLam, oldLam, eta, datafile, iter)

  LOGICAL,                         INTENT(INOUT) :: is_converged  !whether lambda has converged
  TYPE(block_lambda),              INTENT(IN)    :: newLam(:)     !new Lambdas
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: oldLam(:)     !old Lambdas
  REAL(KIND=DP),                   INTENT(IN)    :: eta           !prec
  INTEGER,                         INTENT(INOUT) :: datafile      !datafile for conv data
  INTEGER,                         INTENT(IN)    :: iter          !iter we are on

  INTEGER              :: N_sites, site  !! lambda sites
  LOGICAL, ALLOCATABLE :: lam_conv(:)    !! convergence tracking for individual lambdas

  !! Number of lambdas
  N_sites = SIZE(newLam) 

  !! Convergence tracking for diff lambdas
  ALLOCATE(lam_conv(N_sites)); lam_conv(:) = .FALSE.  

  !! Test Lambda convergence
  DO site=1,N_sites
     CALL test_convergence_vec(lam_conv(site), newLam(site) % m, oldLam(site) % m, eta, datafile, iter, PRINT_ON=.TRUE.)
  end DO
  is_converged = isTrueVec(lam_conv, 'AND')

  !! Set old --> old = new for next iter
  CALL copy_lambda_block(oldLam, newLam)

 END SUBROUTINE test_gen_lambda_convergence



END MODULE peps_pepo_algebra
