MODULE mps_mpo_algebra_inf

 USE utility
 USE definitions_mps_mpo
 USE mps_mpo_utility
 USE boundary_tensors
 USE mps_mpo_algebra_finite

 IMPLICIT NONE

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MULT iMPS by two-body MPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Apply mpo propagator to a given bond of iMPS !!!
 SUBROUTINE mult_imps_2body_mpo(iGamma, iLambda, mpo, chi, eps, SITEORD)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)  !iGamma of MPS
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:) !iLambda of MPS
  TYPE(block_mpo),                 INTENT(IN)    :: mpo(:)     !Two-body Mpo
  INTEGER,                         INTENT(IN)    :: chi        !SVD chi
  REAL(KIND=DP),                   INTENT(IN)    :: eps        !SVD eps
  CHARACTER(LEN=*),                INTENT(IN)    :: SITEORD    !Order of iMPS sites 
                                                               !(to specify which bond we're propagating)
  !! iMPS sites
  INTEGER :: sA, sB

  !! Determine iMPS sites 
  CALL get_imps_sites(sA, sB, SITEORD) 

  !! (1) Absorb iMPS lambdas
  CALL absorb_imps_lambdas(imps=iGamma, iLambda=iLambda, SITEORD=SITEORD, SYM='ALL')

  !! (2) Mult mpoProp * iMPS on bond = SITEORD
  CALL simple_mult_imps_2body_mpo(iGamma, mpo, SITEORD)

  !! (3) Compute SVD of bond = SITEORD -- restore the central lambda -- iLambda(sA)  
  CALL compute_svd_of_imps_bond(iGamma, iLambda(sA), chi, eps, SITEORD)

  !! (4) Restore external Lambdas -- iLambda(sB)
  CALL LambdaDIV(iGamma(sA) % m, iLambda(sB) % m,  '2')
  CALL LambdaDIV(iGamma(sB) % m, iLambda(sB) % m,  '3')
  CALL symmetrize_imps_sites(iGamma)

 END SUBROUTINE mult_imps_2body_mpo





 !!! Apply mpo propagator to a given bond of iMPO !!!
 SUBROUTINE mult_impo_twobody_gate(iGamma, iLambda, mpo, chi, eps, SITEORD)

  USE project_impo

  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)  !iGamma of iMPO
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:) !iLambda of iMPO
  TYPE(block_mpo),                 INTENT(IN)    :: mpo(:)     !Two-body Mpo
  INTEGER,                         INTENT(IN)    :: chi        !SVD chi
  REAL(KIND=DP),                   INTENT(IN)    :: eps        !SVD eps
  CHARACTER(LEN=*),                INTENT(IN)    :: SITEORD    !Order of iMPO sites 
                                                               !(to specify which bond we're propagating)
  !! Local iMPS copy
  TYPE(block_mps),  ALLOCATABLE :: imps(:)
  TYPE(block_lambda)            :: iLmps
  COMPLEX(KIND=DP), ALLOCATABLE :: Xa(:,:), Yb(:,:)

  !! (1) Project iMPO to iMPS
  !CALL project_impo_to_imps(imps, iGamma, iLambda, SITEORD)
  CALL calc_factorized_impo(imps, Xa, Yb, iGamma, iLambda, SITEORD)

  !! (2) Mult PROP * iMPS on bond = SITEORD
  CALL simple_mult_mps_mpo(imps, mpo) 

  !! (3) Compute SVD of bond = SITEORD -- restore the central lambda -- iLambda(sA)  
  CALL compute_svd_of_imps_bond(imps, iLmps, chi, eps, '12')

  !! (4) Symmetrize iMPS sites
  CALL symmetrize_imps_sites(imps)

  !! (5) Restore iMPO from iMPS projection
  !CALL unproject_impo_from_imps(iGamma, iLambda, imps, iLmps, SITEORD)
  CALL restore_unfactorized_impo(iGamma, iLambda, imps, iLmps, Xa, Yb, SITEORD)

 END SUBROUTINE mult_impo_twobody_gate





 !!! Apply mpo propagator to a given bond of iMPO !!!
 SUBROUTINE OLD_mult_impo_twobody_gate(iGamma, iLambda, mpo, chi, eps, SITEORD)

  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)  !iGamma of iMPO
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:) !iLambda of iMPO
  TYPE(block_mpo),                 INTENT(IN)    :: mpo(:)     !Two-body Mpo
  INTEGER,                         INTENT(IN)    :: chi        !SVD chi
  REAL(KIND=DP),                   INTENT(IN)    :: eps        !SVD eps
  CHARACTER(LEN=*),                INTENT(IN)    :: SITEORD    !Order of iMPO sites 
                                                               !(to specify which bond we're propagating)
  !! Local iMPO
  TYPE(block_mps), ALLOCATABLE :: iGmps(:)

  !! iMPO sites && dims
  INTEGER :: dims(2)
  INTEGER :: sA, sB

  !! Determine iMPO sites 
  CALL get_impo_sites(sA, sB, SITEORD)

  !! (1) Absorb iMPO lambdas
  CALL absorb_impo_lambdas(impo=iGamma, iLambda=iLambda, SITEORD=SITEORD, SYM='ALL')

  !! (2) Mult PROP * iMPO on bond = SITEORD
  CALL simple_mult_impo_twobody_gate(iGamma, mpo, SITEORD)

  !! (3) Reshape iMPO -> iMPS
  CALL reshape_MPO_to_MPS(iGmps, iGamma)

  !! (4) Compute SVD of bond = SITEORD -- restore the central lambda -- iLambda(sA)  
  CALL compute_svd_of_imps_bond(iGmps, iLambda(sA), chi, eps, SITEORD)

  !! (5) Restore external Lambdas -- iLambda(sB)
  CALL LambdaDIV(iGmps(sA) % m,  iLambda(sB) % m,  '2')
  CALL LambdaDIV(iGmps(sB) % m,  iLambda(sB) % m,  '3')
  CALL symmetrize_imps_sites(iGmps)

  !! (6) Reshape back iMPS -> iMPO
  CALL reshape_MPS_to_MPO(iGamma, iGmps, (/ iGamma(1) % Sdim, iGamma(1) % Ndim /))

 END SUBROUTINE OLD_mult_impo_twobody_gate






 !!! Apply mpo propagator to a given bond of iMPS !!!
 SUBROUTINE mult_purified_imps_2body_mpo(iGamma, iLambda, mpo, chi, eps, SITEORD)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)  !iGamma of MPS
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:) !iLambda of MPS
  TYPE(block_mpo),                 INTENT(IN)    :: mpo(:)     !Two-body Mpo
  INTEGER,                         INTENT(IN)    :: chi        !SVD chi
  REAL(KIND=DP),                   INTENT(IN)    :: eps        !SVD eps
  CHARACTER(LEN=*),                INTENT(IN)    :: SITEORD    !Order of iMPS sites 
                                                               !(to specify which bond we're propagating)
  !! iMPS sites
  INTEGER :: sA, sB

  !! Determine iMPS sites 
  CALL get_imps_sites(sA, sB, SITEORD) 

  !! (1) Absorb iMPS lambdas
  CALL absorb_imps_lambdas(imps=iGamma, iLambda=iLambda, SITEORD=SITEORD, SYM='ALL')

  !! (2) Mult mpoProp * iMPS on bond = SITEORD
  CALL BROKEN_simple_mult_purified_imps_2body_mpo(iGamma, mpo, SITEORD)

  !! (3) Compute SVD of bond = SITEORD -- restore the central lambda -- iLambda(sA)  
  CALL compute_svd_of_imps_bond(iGamma, iLambda(sA), chi, eps, SITEORD)

  !! (4) Restore external Lambdas -- iLambda(sB)
  CALL LambdaDIV(iGamma(sA) % m, iLambda(sB) % m,  '2')
  CALL LambdaDIV(iGamma(sB) % m, iLambda(sB) % m,  '3')
  CALL symmetrize_imps_sites(iGamma)

 END SUBROUTINE mult_purified_imps_2body_mpo




!!! Simple mult mpoProp * iMPS !!!
 SUBROUTINE BROKEN_simple_mult_purified_imps_2body_mpo(imps, mpo, SITEORD)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: imps(:)   !iMPS
  TYPE(block_mpo),              INTENT(IN)    :: mpo(:)    !Two-body Mpo
  CHARACTER(LEN=*),             INTENT(IN)    :: SITEORD 

  !! Full iMPO (restored from purified iMPS)   
  TYPE(block_mpo), ALLOCATABLE :: impoOUT(:)             

  !! iMPS sites && dims
  INTEGER :: localDim
  INTEGER :: sA, sB
  
  !! Determine iMPS sites
  CALL get_imps_sites(sA, sB, SITEORD)

  !! Reshape iMPS --> iMPO (dn-leg = real Hilbert space, up-leg = ancilla space)
  localDim = SQRTint(imps(1) % SNdim)
  CALL reshape_MPS_to_MPO(impoOUT, imps, (/ localDim, localDim /))

  !! Multiply MPO_PROP * iMPO
  CALL allocate_mpo_site(impoOUT(sA), TENSMUL(mpo(1) % m, impoOUT(sA) % m, MULT='21', FUSE='(33,44)'))
  CALL allocate_mpo_site(impoOUT(sB), TENSMUL(mpo(2) % m, impoOUT(sB) % m, MULT='21', FUSE='(33,44)'))

  !! Reshape iMPS --> iMPO
  CALL reshape_MPO_to_MPS(imps, impoOUT)

 END SUBROUTINE BROKEN_simple_mult_purified_imps_2body_mpo








 !!! Multiply iMPS by dense propagator, and compute SVD !!!
 SUBROUTINE mult_imps_dense_prop(imps, iLambda, prop_dense, chi0, eps0, SITEORD)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT)  :: imps(:)             !input/output imps
  TYPE(block_lambda),           INTENT(INOUT)  :: iLambda             !output only the central lambda returned by SVD
  COMPLEX(KIND=DP),             INTENT(IN)     :: prop_dense(:,:)     !Propagator
  INTEGER,                      INTENT(IN)     :: chi0
  REAL(KIND=DP),                INTENT(IN)     :: eps0
  CHARACTER(LEN=*),             INTENT(IN)     :: SITEORD

  !! iMPS sites
  INTEGER :: sA, sB

  !! iMPS sites combined into a single theta tensor
  COMPLEX(KIND=DP), ALLOCATABLE :: thetaTens(:,:,:,:)
  INTEGER                       :: dimT(4)

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:), VH(:,:), theta(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)
  INTEGER                       :: chi

  !! Determine iMPS sites
  CALL get_imps_sites(sA, sB, SITEORD)

  !! Construct thetaTens
  thetaTens = TENSMUL(imps(sA) % m, imps(sB) % m, '32')
  dimT = shape(thetaTens)

  !! Contract theta tensor with propagator
  thetaTens = RESHAPE_4D(TENSMUL(prop_dense, RESHAPE_2D(thetaTens, '12,34')), '12,34', (/dimT(1), dimT(2)/), (/dimT(3), dimT(4)/))

  !! Construct theta matrix
  theta = RESHAPE_2D(thetaTens, '13,24')

  !! Frobenius-normalize theta before SVD
  CALL normalize_matrix(theta)

  !! Compute SVD
  chi = chi0 
  CALL compute_lapack_svd(theta, U, VH, Sigma, chi, eps0)

  !! Copy back SVD results 
  CALL allocate_mps_site(imps(sA), RESHAPE_3D(U,  '12,3',  (/ imps(sA) % SNdim, imps(sA) % Wdim /)))
  CALL allocate_mps_site(imps(sB), RESHAPE_3D(VH, '2,13',  (/ imps(sB) % SNdim, imps(sB) % Edim /)))
  CALL allocate_lambda_site(iLambda, complx(Sigma))

 END SUBROUTINE mult_imps_dense_prop










 !!! Multiply iMPS by dense propagator, and compute SVD !!!
 SUBROUTINE ADVANCED_mult_imps_dense_prop(imps, iLambda, prop_dense, chi0, eps0, SITEORD)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT)  :: imps(:)             !input/output imps
  TYPE(block_lambda),           INTENT(INOUT)  :: iLambda             !output only the central lambda returned by SVD
  COMPLEX(KIND=DP),             INTENT(IN)     :: prop_dense(:,:,:,:) !Propagator
  INTEGER,                      INTENT(IN)     :: chi0
  REAL(KIND=DP),                INTENT(IN)     :: eps0
  CHARACTER(LEN=*),             INTENT(IN)     :: SITEORD

  !! iMPS sites
  INTEGER :: sA, sB

  !! Reshaped versions of propagator && MPS sites
  COMPLEX(KIND=DP), ALLOCATABLE :: prop_3D(:,:,:), mpsA(:,:), mpsB(:,:)
  INTEGER                       :: dimA(3), dimProp(4)

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:), VH(:,:), theta(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)
  INTEGER                       :: chi

  !! Determine iMPS sites
  CALL get_imps_sites(sA, sB, SITEORD)

  !! Get dims
  dimProp = shape(prop_dense)
  dimA    = shape(imps(sA) % m)

  !! Contract prop * mpsA
  prop_3D = RESHAPE_3D(prop_dense,  '12,3,4')
  mpsA    = RESHAPE_2D(imps(sA) % m,  '1,23')
  prop_3D = TENSMUL(prop_3D, mpsA, '21')

  !! Restore legs (2,3) of mpsA, then fuse legs (3,4) of prop_3D  swap (3,4) legs cause  they must correspond to (1,2) legs of mpsB
  prop_3D = RESHAPE_3D(RESHAPE_4D(prop_3D, '1,23,4', (/dimA(2), dimA(3)/)), '1,2,34')

  !! Contract prop * mpsB
  !! leg-3 of prop corresponds to leg-2 of mpsB, 
  !! and leg-4 of prop to leg-1 of mpsB --> so swap (1,2) before fusing so that both prop_3D and mpsB legs are combined in the same order
  mpsB    = RESHAPE_2D(TensTRANSPOSE(imps(sB) % m, '12'),  '12,3')
  prop_3D = TENSMUL(prop_3D, mpsB, '31')

  !! Reshape the final prop into theta
  theta = RESHAPE_2D(RESHAPE_4D(prop_3D, '12,3,4', (/dimProp(1), dimProp(2)/)), '13,24')

  !! Compute SVD
  chi = chi0 
  CALL compute_lapack_svd(theta, U, VH, Sigma, chi, eps0)

  !! Copy back SVD results 
  CALL allocate_mps_site(imps(sA), RESHAPE_3D(U,  '12,3',  (/ imps(sA) % SNdim, imps(sA) % Wdim /)))
  CALL allocate_mps_site(imps(sB), RESHAPE_3D(VH, '2,13',  (/ imps(sB) % SNdim, imps(sB) % Edim /)))
  CALL allocate_lambda_site(iLambda, complx(Sigma))

 END SUBROUTINE ADVANCED_mult_imps_dense_prop






 !!! Simple mult mpoProp * iMPS !!!
 SUBROUTINE simple_mult_imps_2body_mpo(imps, mpo, SITEORD)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: imps(:)   !iMPS
  TYPE(block_mpo),              INTENT(IN)    :: mpo(:)    !Two-body Mpo
  CHARACTER(LEN=*),             INTENT(IN)    :: SITEORD                 

  !! iMPS sites
  INTEGER :: sA, sB
  
  !! Determine iMPS sites
  CALL get_imps_sites(sA, sB, SITEORD)

  !! Multiply MPO_PROP * iMPS 
  CALL allocate_mps_site(imps(sA), TENSMUL(mpo(1) % m, imps(sA) % m, MULT='21', FUSE='(32,43)'))
  CALL allocate_mps_site(imps(sB), TENSMUL(mpo(2) % m, imps(sB) % m, MULT='21', FUSE='(32,43)'))

 END SUBROUTINE simple_mult_imps_2body_mpo





 !!! Simple mult gate_mpo * iMPO !!!
 SUBROUTINE simple_mult_impo_twobody_gate(impo, gate_mpo, SITEORD)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: impo(:)     !iMPO
  TYPE(block_mpo),              INTENT(IN)    :: gate_mpo(:) !Two-body Gate Mpo
  CHARACTER(LEN=*),             INTENT(IN)    :: SITEORD                 

  !! iMPO sites
  INTEGER :: sA, sB
  
  !! Determine iMPO sites
  CALL get_impo_sites(sA, sB, SITEORD)

  !! Multiply MPO_PROP * iMPO --> iMPO
  CALL allocate_mpo_site(impo(sA), TENSMUL(gate_mpo(1) % m, impo(sA) % m, MULT='21', FUSE='(33,44)'))
  CALL allocate_mpo_site(impo(sB), TENSMUL(gate_mpo(2) % m, impo(sB) % m, MULT='21', FUSE='(33,44)'))

 END SUBROUTINE simple_mult_impo_twobody_gate


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Kraus update for iMPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Apply Kraus propagator to iMPO !!!
 SUBROUTINE kraus_update_impo(impo, prop_kraus, chi, eps, SITEORD)

  TYPE(block_mpo),  ALLOCATABLE, INTENT(INOUT) :: impo(:)           !iMPO
  COMPLEX(KIND=DP),              INTENT(IN)    :: prop_kraus(:,:,:) !Kraus propagator
  INTEGER,                       INTENT(IN)    :: chi               !SVD chi
  REAL(KIND=DP),                 INTENT(IN)    :: eps               !SVD eps
  CHARACTER(LEN=*),              INTENT(IN)    :: SITEORD

  !! iMPO sites
  INTEGER :: sA, sB

  !! Determine iMPO sites
  CALL get_imps_sites(sA, sB, SITEORD)

  !! Kraus-propagate iMPO sites
  CALL kraus_update_impo_site(impo(sA), prop_kraus, chi, eps, '1')
  CALL kraus_update_impo_site(impo(sB), prop_kraus, chi, eps, '2')

 END SUBROUTINE kraus_update_impo




 !!! Apply Kraus propagator to iMPO site !!!
 SUBROUTINE kraus_update_impo_site(impo, prop_kraus, chi0, eps, SITE)

  TYPE(block_mpo),  INTENT(INOUT) :: impo              !iMPO site
  COMPLEX(KIND=DP), INTENT(IN)    :: prop_kraus(:,:,:) !Kraus propagator
  INTEGER,          INTENT(IN)    :: chi0              !SVD chi
  REAL(KIND=DP),    INTENT(IN)    :: eps               !SVD eps
  CHARACTER(LEN=*), INTENT(IN)    :: SITE              !Site flag

  COMPLEX(KIND=DP), ALLOCATABLE :: new_impo(:,:,:,:), theta(:,:)  
  COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:), VH(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)
  INTEGER                       :: chi

  CHARACTER(LEN=16) :: MULT
  INTEGER           :: dimM(4)

  !! Decide in which order to fuse Kraus legs
  SELECT CASE(SITE)
  CASE('1')
       MULT = '22(AB)'
  CASE('2')
       MULT = '22(BA)'
  CASE DEFAULT
       CALL invalid_flag("kraus_update_impo_site -- invalid site ", SITE)
  end SELECT

  !! (1) Contract iMPO site with Kraus propagator (transpose to put in the same config as Werner et al)
  new_impo = TensTRANSPOSE(impo % m, '12')
  !new_impo = TENSMUL(new_impo, prop_kraus, MULT=MULT)
  new_impo = TENSMUL(new_impo, prop_kraus, MULT='22') 
  dimM = shape(new_impo)

  !! (2) Reshape into theta matrix
  theta = RESHAPE_2D(new_impo, '1,234')

  !! (3) Compute SVD && truncate
  chi = chi0 
  CALL compute_lapack_svd(theta, U, VH, Sigma, chi, eps)

  !! (4) Copy back VH to iMPO (transpose back to restore from the config of Werner et al)
  CALL LambdaMUL(VH, complx(Sigma), MULT='1') 
  CALL allocate_mpo_site(impo, TensTRANSPOSE(RESHAPE_4D(VH, '1,234', (/ dimM(2), dimM(3), dimM(4) /)), '12'))

 END SUBROUTINE kraus_update_impo_site


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MULT && SVD iMPS-iMPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Multiply iMPS with binary iMPO !!!
 SUBROUTINE mult_imps_bimpo(iGamma, iLambda, impo, chi, eps)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  TYPE(block_mpo),                 INTENT(IN)    :: impo(:,:)
  INTEGER,                         INTENT(IN)    :: chi
  REAL(KIND=DP),                   INTENT(IN)    :: eps

  !! Mult && Canonicalize with impo(2,:) and with impo(1,:)
  CALL mult_imps_impo(iGamma, iLambda, impo(2,:), chi, eps)
  CALL mult_imps_impo(iGamma, iLambda, impo(1,:), chi, eps)

 END SUBROUTINE mult_imps_bimpo

  



 !!! Multiply iMPS with iMPO !!!
 SUBROUTINE mult_imps_impo(iGamma, iLambda, impo, chi, eps)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  TYPE(block_mpo),                 INTENT(IN)    :: impo(:)
  INTEGER,                         INTENT(IN)    :: chi
  REAL(KIND=DP),                   INTENT(IN)    :: eps

  !! Mult && Canonicalize
  CALL simple_mult_imps_impo(iGamma, iLambda, impo)
  CALL canonicalize_imps(iGamma, iLambda, chi, eps)

 END SUBROUTINE mult_imps_impo





 !!! Simple mult iMPS * iMPO (without SVD) !!!
 SUBROUTINE simple_mult_imps_impo(iGamma, iLambda, impo)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  TYPE(block_mpo),                 INTENT(IN)    :: impo(:)

  INTEGER :: N_sites, site

  N_sites = SIZE(iGamma)

  !! MPS*MPO mult at site=1,2
  DO site=1,N_sites
     CALL allocate_mps_site(iGamma(site), TENSMUL(impo(site) % m, iGamma(site) % m, MULT='21', FUSE='(32,43)'))
     CALL mult_lambda_site(iLambda(site), impo(site), '3')
  end DO
 
 END SUBROUTINE simple_mult_imps_impo




 !!! Multiply lambda site of lambda_block !!!
 SUBROUTINE mult_lambda_site(lambda_site, mpo_site, DIR)

  TYPE(block_lambda), INTENT(INOUT) :: lambda_site
  TYPE(block_mpo),    INTENT(IN)    :: mpo_site
  CHARACTER(LEN=*),   INTENT(IN)    :: DIR

  COMPLEX(KIND=DP), ALLOCATABLE :: lambdaMat(:,:), eyeTens(:,:)

  !! Setup eyeTens to represent mpo bond
  SELECT CASE(DIR)
  CASE('2') 
       eyeTens = matEye(mpo_site % Wdim)
  CASE('3')
       eyeTens = matEye(mpo_site % Edim)
  CASE DEFAULT
       CALL invalid_flag("mult_lambda_site: invalid DIR ", DIR)
  end SELECT

  !! Represent lambda as diagonal matrix
  lambdaMat = DiagonalAsMat(lambda_site % m)

  !! Tensor product = Eye x Lambda
  lambdaMat = TensKRON(eyeTens, lambdaMat)

  !! Write to lambda_site
  CALL allocate_lambda_site(lambda_site, DiagonalAsVec(lambdaMat))

 END SUBROUTINE mult_lambda_site





 !!! Compute SVD of iMPS (all bonds) !!!
 SUBROUTINE compress_imps(iGamma, iLambda, chi0, eps0)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  INTEGER,                         INTENT(IN)    :: chi0
  REAL(KIND=DP),                   INTENT(IN)    :: eps0

  !! Compute svd of imps bonds
  CALL compress_imps_bond(iGamma, iLambda, chi0, eps0, '12', 'absorb-lambda', 'restore-lambda')
  CALL compress_imps_bond(iGamma, iLambda, chi0, eps0, '21', 'absorb-lambda', 'restore-lambda')

 END SUBROUTINE compress_imps




 !!! Compress iMPS bond by SVD !!!
 SUBROUTINE compress_imps_bond(iGamma, iLambda, chi0, eps0, SITEORD, abs_lambda, res_lambda)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT)        :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)        :: iLambda(:)
  INTEGER,                         INTENT(IN)           :: chi0
  REAL(KIND=DP),                   INTENT(IN)           :: eps0
  CHARACTER(LEN=*),                INTENT(IN)           :: SITEORD
  CHARACTER(LEN=*),                INTENT(IN), OPTIONAL :: abs_lambda, res_lambda !Whether to absorb/restore external lambdas

  !! iMPS sites
  INTEGER :: sA, sB

  !! Determine iMPS sites 
  CALL get_imps_sites(sA, sB, SITEORD)

  !! Setup nolambda iMPS -- absorb all lambdas into gammas
  IF(PRESENT(abs_lambda)) CALL absorb_imps_lambdas(imps=iGamma, iLambda=iLambda, SITEORD=SITEORD, SYM='ALL')

  !! Calculate SVD of iMPS bond
  CALL compute_svd_of_imps_bond(iGamma, iLambda(sA), chi0, eps0, SITEORD)

  !! Restore external lambdas = iLambda(sB)
  IF(PRESENT(res_lambda)) THEN
     CALL LambdaDIV(iGamma(sA) % m, iLambda(sB) % m,  '2')
     CALL LambdaDIV(iGamma(sB) % m, iLambda(sB) % m,  '3')
  end IF

 END SUBROUTINE compress_imps_bond




 !!! Compute SVD of imps bond !!!
 SUBROUTINE compute_svd_of_imps_bond(imps, iLambda, chi0, eps0, SITEORD)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT)        :: imps(:)    !input/output imps
  TYPE(block_lambda),           INTENT(INOUT)        :: iLambda    !output only the central lambda returned by SVD
  INTEGER,                      INTENT(IN)           :: chi0
  REAL(KIND=DP),                INTENT(IN)           :: eps0
  CHARACTER(LEN=*),             INTENT(IN)           :: SITEORD

  !! iMPS sites
  INTEGER :: sA, sB

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:), VH(:,:), theta(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)
  INTEGER :: chi

  !! Determine iMPS sites
  CALL get_imps_sites(sA, sB, SITEORD)

  !! Compute SVD
  theta = TENSMUL(RESHAPE_2D(imps(sA)%m, '12,3'), RESHAPE_2D(imps(sB)%m, '2,13'))
  chi = chi0 
  CALL compute_lapack_svd(theta, U, VH, Sigma, chi, eps0)

  !! Copy back SVD results 
  CALL allocate_mps_site(imps(sA), RESHAPE_3D(U,  '12,3',  (/ imps(sA) % SNdim, imps(sA) % Wdim /)))
  CALL allocate_mps_site(imps(sB), RESHAPE_3D(VH, '2,13',  (/ imps(sB) % SNdim, imps(sB) % Edim /)))
  CALL allocate_lambda_site(iLambda, complx(Sigma))

 END SUBROUTINE compute_svd_of_imps_bond





 !!! Balance MPS bond using the proposal in Lubasch et al !!!
 SUBROUTINE balance_imps_bond(imps, iLambda)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: imps(:)    !input/output imps
  TYPE(block_lambda),              INTENT(INOUT) :: iLambda    !output only the central lambda returned by SVD

  !! Factorized MPS sites
  COMPLEX(KIND=DP), ALLOCATABLE :: mpsAA(:,:,:), mpsBB(:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: Xa(:,:), Yb(:,:)

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:), VH(:,:), theta(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)
  INTEGER                       :: chi

  !! Absorb Lambda
  CALL LambdaMUL(imps(1) % m, iLambda % m, imps(2) % m)

  !! Factorize MPS sites
  CALL LHS_factorize_Tens3D(mpsAA, Xa, imps(1) % m, '12,3')
  CALL RHS_factorize_Tens3D(mpsBB, Yb, imps(2) % m, '2,13')

  !! Compute SVD
  chi = -1 
  CALL compute_lapack_svd(TENSMUL(Xa, Yb), U, VH, Sigma, chi, eps=2.0D0)

  !! Copy back SVD results 
  CALL allocate_mps_site(imps(1),     TENSMUL(mpsAA, U,  '31'))
  CALL allocate_mps_site(imps(2),     TENSMUL(mpsBB, VH, '22'))
  CALL allocate_lambda_site(iLambda,  complx(Sigma))

 END SUBROUTINE balance_imps_bond

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SVD iMPO bond !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compress iMPO bond by SVD !!!
 SUBROUTINE compress_impo_bond(iGamma, iLambda, chi0, eps0, SITEORD, abs_lambda, res_lambda)

  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT)        :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)        :: iLambda(:)
  INTEGER,                         INTENT(IN)           :: chi0
  REAL(KIND=DP),                   INTENT(IN)           :: eps0
  CHARACTER(LEN=*),                INTENT(IN)           :: SITEORD
  CHARACTER(LEN=*),                INTENT(IN), OPTIONAL :: abs_lambda, res_lambda !Whether to absorb/restore external lambdas

  !! iMPO sites
  INTEGER :: sA, sB

  !! Determine iMPS sites 
  CALL get_impo_sites(sA, sB, SITEORD)

  !! Setup nolambda iMPS -- absorb all lambdas into gammas
  IF(PRESENT(abs_lambda)) CALL absorb_impo_lambdas(impo=iGamma, iLambda=iLambda, SITEORD=SITEORD, SYM='ALL')

  !! Calculate SVD of iMPO bond
  CALL compute_svd_of_impo_bond(iGamma, iLambda(sA), chi0, eps0, SITEORD)

  !! Restore external lambdas = iLambda(sB)
  IF(PRESENT(res_lambda)) THEN
     CALL LambdaDIV(iGamma(sA) % m, iLambda(sB) % m,  '3')
     CALL LambdaDIV(iGamma(sB) % m, iLambda(sB) % m,  '4')
  end IF

 END SUBROUTINE compress_impo_bond



 !!! Compute SVD of impo bond !!!
 SUBROUTINE compute_svd_of_impo_bond(impo, iLambda, chi0, eps0, SITEORD)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT)  :: impo(:)    !input/output impo
  TYPE(block_lambda),           INTENT(INOUT)  :: iLambda    !output only the central lambda returned by SVD
  INTEGER,                      INTENT(IN)     :: chi0
  REAL(KIND=DP),                INTENT(IN)     :: eps0
  CHARACTER(LEN=*),             INTENT(IN)     :: SITEORD

  !! Local iMPS && its sites
  TYPE(block_mps), ALLOCATABLE :: imps(:)
  INTEGER :: sA, sB

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:), VH(:,:), theta(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)
  INTEGER :: chi

  !! Reshape iMPO -> iMPS
  CALL reshape_MPO_to_MPS(imps, impo)

  !! Determine iMPS sites
  CALL get_imps_sites(sA, sB, SITEORD)

  !! Compute SVD
  theta = TENSMUL(RESHAPE_2D(imps(sA)%m, '12,3'), RESHAPE_2D(imps(sB)%m, '2,13'))
  chi = chi0 
  CALL compute_lapack_svd(theta, U, VH, Sigma, chi, eps0)

  !! Copy back SVD results 
  CALL allocate_mps_site(imps(sA), RESHAPE_3D(U,  '12,3',  (/ imps(sA) % SNdim, imps(sA) % Wdim /)))
  CALL allocate_mps_site(imps(sB), RESHAPE_3D(VH, '2,13',  (/ imps(sB) % SNdim, imps(sB) % Edim /)))
  CALL allocate_lambda_site(iLambda, complx(Sigma))

  !! Reshape iMPS -> iMPO
  CALL reshape_MPS_to_MPO(impo, imps, (/ impo(1) % Sdim, impo(1) % Ndim /))

 END SUBROUTINE compute_svd_of_impo_bond

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MULT iMPS-iMPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Contract iMPS-iMPS network to obtain expval = [result of contraction] !!!
 SUBROUTINE contract_imps_imps(expval, iG, iL, iGP, iLP, SITEORDin)

  COMPLEX(KIND=DP),   INTENT(OUT)          :: expval
  TYPE(block_mps),    INTENT(IN)           :: iG(:), iGP(:) 
  TYPE(block_lambda), INTENT(IN)           :: iL(:), iLP(:)
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: SITEORDin

  !! Boundary vecs
  COMPLEX(KIND=DP), ALLOCATABLE :: Rvec(:,:), Lvec(:,:) 
  CHARACTER(LEN=2)              :: SITEORD

  !! Treat factorized iMPS cases separately
  IF(is_one_imps_bond(iL) .OR. is_one_imps_bond(iLP)) THEN
     CALL contract_singleton_imps_imps(expval, iG, iL, iGP, iLP)
     RETURN
  end IF

  !! Determine site order
  SITEORD = OptArg(SITEORDin, '12')

  !! Get expval of contraction = eval of imps-imps transfer matrix
  CALL calc_boundary_vecs(expval, Rvec, Lvec, iG, iL, iGP, iLP, SITEORD)

 END SUBROUTINE contract_imps_imps




 !!! Normalize iMPS so that C_norm = 1 !!!
 SUBROUTINE normalize_imps(iGamma, iLambda, SITEORDin)

  TYPE(block_mps),    INTENT(INOUT)        :: iGamma(:)
  TYPE(block_lambda), INTENT(INOUT)        :: iLambda(:)
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: SITEORDin

  !! Local iMPS objects (HC Copy)
  TYPE(block_mps),    ALLOCATABLE :: iGammaP(:) 
  TYPE(block_lambda), ALLOCATABLE :: iLambdaP(:)
  COMPLEX(KIND=DP)                :: C_norm

  !! iMps sites
  INTEGER          :: sA, sB
  CHARACTER(LEN=2) :: SITEORD

  !! Determine iMps sites
  SITEORD = OptArg(SITEORDin, '12')
  CALL get_imps_sites(sA, sB, SITEORD)

  !! Rescale iMPS so that each site has the same max element 
  CALL symmetrize_imps_sites(iGamma)

  !! Normalize iMPS so that C_norm = 1 (take a sqrt to get the amount of normalization per iMPS)
  CALL copy_imps(iGammaP, iLambdaP, iGamma, iLambda, 'HC')
  CALL contract_imps_imps(C_norm, iGamma, iLambda, iGammaP, iLambdaP, SITEORD)
  CALL mult_mps_by_value(iGamma, C_norm**(-0.5d0)) 

 END SUBROUTINE normalize_imps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Handling iMPS with chi=1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Check if iMPS has bonds equal to one !!!
 FUNCTION is_one_imps_bond(iLambda)

  TYPE(block_lambda), INTENT(IN) :: iLambda(:)
  LOGICAL :: is_one_imps_bond

  is_one_imps_bond = (iLambda(1)%ChiDim .EQ. 1) .OR. (iLambda(2)%ChiDim .EQ. 1)

 END FUNCTION is_one_imps_bond



 !!! Contract two factorized iMPS
 !!! -- note that this is only used in iMPS routines for normalization or canonicalization purposes
 SUBROUTINE contract_singleton_imps_imps(expval, iG, iL, iGP, iLP)

  COMPLEX(KIND=DP),   INTENT(OUT) :: expval
  TYPE(block_mps),    INTENT(IN)  :: iG(:), iGP(:) 
  TYPE(block_lambda), INTENT(IN)  :: iL(:), iLP(:)

  !! Boundary vecs
  COMPLEX(KIND=DP), ALLOCATABLE :: Rvec(:), Lvec(:)
  COMPLEX(KIND=DP)              :: eval

  !! Nolambda iMPS
  TYPE(block_mps),  ALLOCATABLE :: imps(:), impsP(:) 
  INTEGER                       :: chi(2), chiP(2)

  !! Logical vars to check which bonds are one
  LOGICAL :: is_one_all_bonds, is_one_bond_12, is_one_bond_21

  !! MPS overlap
  COMPLEX(KIND=DP), ALLOCATABLE :: mps_overlap(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: overlap_1(:,:), overlap_2(:,:)

  !! (0) Get bond dims
  chi  = (/ iL(1)%ChiDim,  iL(2)%ChiDim  /)
  chiP = (/ iLP(1)%ChiDim, iLP(2)%ChiDim /)

  !! (1) Check which bonds are one
  is_one_all_bonds = (chi(1) .EQ. 1) .AND. (chi(2) .EQ. 1) .AND. (chiP(1) .EQ. 1) .AND. (chiP(2) .EQ. 1)
  is_one_bond_12   = (chi(1) .EQ. 1) .OR. (chiP(1) .EQ. 1)
  is_one_bond_21   = (chi(2) .EQ. 1) .OR. (chiP(2) .EQ. 1)

  !! (2) Otherwise check all cases with bonds equal to one
  IF(is_one_all_bonds) THEN 

      !! If all bonds dims = 1 --> take expval = 1, don't do any calculations 
      !! (NB. this is valid since it will only be used for normalization or canonicalization purposes 
      !!  as this routine is called only by iMPS procedures -- finite mps block contractions will never use this)
      expval = (1.0d0,0.0d0)
      RETURN

  ELSEIF(is_one_bond_12) THEN

      !! If only the '12' bond equals 1, construct nolambda imps such that singleton bond '12' is facing outwards (i.e. SITEORD = '21')
      CALL absorb_imps_lambdas(imps,  iG,  iL,  '21')
      CALL absorb_imps_lambdas(impsP, iGP, iLP, '21')

  ELSEIF(is_one_bond_21) THEN

      !! If only the '21' bond equals 1, construct nolambda imps such that singleton bond '21' is facing outwards (i.e. SITEORD = '12')
      CALL absorb_imps_lambdas(imps,  iG,  iL,  '12')
      CALL absorb_imps_lambdas(impsP, iGP, iLP, '12')
  ELSE
      !! In other cases, set expval to 1 to ensure we are not returning an empty variable
      expval = (1.0d0, 0.0d0)
      RETURN
  end IF

  !! (3) Compute overlap of two nolambda iMPS -- if either bond_12_is_one or bond_21_is_one, but not both
  overlap_1 = TENSMUL(imps(1) % m, impsP(1) % m, MULT='11', FUSE='(22,33)')
  overlap_2 = TENSMUL(imps(2) % m, impsP(2) % m, MULT='11', FUSE='(22,33)')
  mps_overlap = TENSMUL(overlap_1, overlap_2)

  !! (4) Find dominant eval of mps_overlap --> gives expval
  CALL dominant_eig(eval, Rvec, Lvec, mps_overlap)
  expval = eval/SUM(Rvec*Lvec)

 END SUBROUTINE contract_singleton_imps_imps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CANONICALIZE iMPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Canonicalize iMPS using the algorithm by Orus and Vidal !!!
 SUBROUTINE canonicalize_imps(iGamma, iLambda, chi, eps)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  INTEGER,                         INTENT(IN)    :: chi
  REAL(KIND=DP),                   INTENT(IN)    :: eps

  !! Boundary matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: X(:,:), Y(:,:)
 
  !! Do not canonicalize factorized iMPS
  IF(is_one_imps_bond(iLambda)) RETURN

  !! When using fixed chi, pretruncate iMPS to discard near-zero lambdas && rescale it so that each site has the same max element
  IF((eps .GT. 1.0D0) .OR. (eps .LT. 0.0D0)) THEN
      CALL compress_imps_bond(iGamma, iLambda, 2, 1.0d-12, '12', 'absorb-lambda', 'restore-lambda')
      CALL compress_imps_bond(iGamma, iLambda, 2, 1.0d-12, '21', 'absorb-lambda', 'restore-lambda')
  end IF

  WRITE(*,*) "CANONICALIZING iMPS"

  !! (1) Find the orthogonal boundary matrices of iMPS-iMPS network
  CALL calc_orthodecomposed_boundaries(X, Y, iGamma, iLambda, '12')

  !! (2) Transform iMPS into canonical form using bounds X,Y (to account for the environment of a repeated unit)
  CALL transform_GammaLambda_into_canonical_form(iGamma, iLambda, X, Y, chi, eps, '12') !--truncate here

  !! (3) Renormalize iMPS so that the dominant eval of TN = 1
  CALL normalize_imps(iGamma, iLambda)

 END SUBROUTINE canonicalize_imps





 !!! Transform iGamma, iLambda tensors to the iMPS canonical form !!!
 SUBROUTINE transform_GammaLambda_into_canonical_form(iGamma, iLambda, X, Y, chi0, eps, SITEORD)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  COMPLEX(KIND=DP),                INTENT(INOUT) :: X(:,:), Y(:,:)
  INTEGER,                         INTENT(IN)    :: chi0
  REAL(KIND=DP),                   INTENT(IN)    :: eps
  CHARACTER(LEN=*),                INTENT(IN)    :: SITEORD
 
  !! iMPS sites
  INTEGER :: sA, sB

  !! SVD matrices && params
  COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:), VH(:,:), theta(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)
  INTEGER :: chi

  !! Inverses of X,Y matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: X_inv(:,:), Y_inv(:,:)

  !! Determine iMps sites
  CALL get_imps_sites(sA, sB, SITEORD)

  !! (1) Compute the inverse of X, Y (do this before contracting X,Y with anything else!!!) -- DO NOT FILTER !!!
  X_inv = invert_matrix(X)
  Y_inv = invert_matrix(Y)

  !! (2) Construct theta = Y*X, compute SVD
  CALL LambdaMUL(Y, iLambda(sB) % m, X)
  theta = TENSMUL(Y, X) 
  chi = chi0
  CALL compute_lapack_svd(theta, U, VH, Sigma, chi, eps) 

  !! (3) Obtain canonicalized external iLambda, rescale so that SUM(iLambda**2) = 1
  CALL allocate_lambda_site(iLambda(sB), complx(Sigma)) 

  !! (4) Mult iMPS sites by VH * X^-1 and Y^-1 * U --> Obtain canonicalized iGamma
  CALL allocate_mps_site(iGamma(sA), TENSMUL(iGamma(sA) % m,  TENSMUL(VH, X_inv),  '2'))
  CALL allocate_mps_site(iGamma(sB), TENSMUL(iGamma(sB) % m,  TENSMUL(Y_inv, U),   '3'))
  
  !! (5) Compute SVD of the [A--B] bond
  CALL compress_imps_bond(iGamma, iLambda, chi0, eps, SITEORD, 'absorb-lambda', 'restore-lambda')

 END SUBROUTINE transform_GammaLambda_into_canonical_form

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE mps_mpo_algebra_inf
