MODULE observables_ipeps

 USE utility
 USE definitions_mps_mpo
 USE definitions_peps_pepo
 USE mps_mpo_utility

 USE project_ipeps
 USE peps_pepo_algebra
 USE TN_contractions

 IMPLICIT NONE

 interface construct_reduced_rho
    module procedure construct_reduced_rho_GamLam
    module procedure construct_reduced_rho_PLAIN
 end interface construct_reduced_rho

 interface calc_CTM_environment_ipeps
    module procedure calc_CTM_environment_ipeps_GamLam
    module procedure calc_CTM_environment_ipeps_PLAIN
 end interface calc_CTM_environment_ipeps

 private calc_CTM_environment_ipeps_GamLam, calc_CTM_environment_ipeps_PLAIN
 private construct_reduced_rho_GamLam, construct_reduced_rho_PLAIN

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Compute CTM environment for either plain or {Gamma, Lambda} iPEPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute environment of iPEPS unit cell !!!
 SUBROUTINE calc_CTM_environment_ipeps_PLAIN(Cmat, Tmat, ipeps, CTMRG, use_old_ctm)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)                !CTM Cmat
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)                !CTM Tmat
  TYPE(block_peps),                     INTENT(INOUT) :: ipeps(:)               !iPEPS
  TYPE(ctm_params),                     INTENT(IN)    :: CTMRG                  !CTM params
  LOGICAL,                              INTENT(IN)    :: use_old_ctm            !whether to use old CTM for initialization

  !! Bulk 2D network (TN2D)
  TYPE(block_mpo),    ALLOCATABLE :: TN2D(:)

  !! Transfer-MPO network tensors (Rvec, Lvec, Tmpo)
  TYPE(transfer_mpo), ALLOCATABLE :: Tmpo(:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: Rvec(:,:), Lvec(:,:)
  COMPLEX(KIND=DP)                :: C_norm

  !! If we're using previous CTM, normalize TN first
  IF(ALLOCATED(Cmat) .AND. ALLOCATED(Tmat) .AND. use_old_ctm) THEN

     CALL create_reduced_TN_from_ipeps(TN2D, ipeps)
     CALL contract_CTM_network(C_norm, Cmat, Tmat, TN2D)

     ipeps(1) % m = ipeps(1) % m / C_norm**(0.25D0)  
     ipeps(2) % m = ipeps(2) % m / C_norm**(0.25D0)  
  end IF

  !! Create reduced TN
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)  

  !! Compute CTMRG boundaries of TN
  CALL compute_CTMRG_boundaries(Cmat, Tmat, TN2D, CTMRG, use_old_ctm) 

 END SUBROUTINE calc_CTM_environment_ipeps_PLAIN





 !!! Compute environment of iPEPS unit cell !!!
 SUBROUTINE calc_CTM_environment_ipeps_GamLam(Cmat, Tmat, iGamma, iLambda, CTMRG, use_old_ctm)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)                !CTM Cmat
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)                !CTM Tmat
  TYPE(block_peps),                     INTENT(INOUT) :: iGamma(:)              !iPEPS gamma
  TYPE(block_lambda),                   INTENT(IN)    :: iLambda(:)             !iPEPS lambda
  TYPE(ctm_params),                     INTENT(IN)    :: CTMRG                  !CTM params
  LOGICAL,                              INTENT(IN)    :: use_old_ctm            !whether to use old CTM for initialization

  !! Bulk 2D network (TN2D)
  TYPE(block_peps),   ALLOCATABLE :: ipeps(:)
  TYPE(block_mpo),    ALLOCATABLE :: TN2D(:)

  !! Transfer-MPO network tensors (Rvec, Lvec, Tmpo)
  TYPE(transfer_mpo), ALLOCATABLE :: Tmpo(:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: Rvec(:,:), Lvec(:,:)
  COMPLEX(KIND=DP)                :: C_norm

  !! Create reduced TN
  CALL absorb_ipeps_lambdas(ipeps, iGamma, iLambda)
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)

  !! If we're using previous CTM, normalize TN first
  IF(ALLOCATED(Cmat) .AND. ALLOCATED(Tmat) .AND. use_old_ctm) THEN

     CALL contract_CTM_network(C_norm, Cmat, Tmat, TN2D)
     ipeps(1) % m = ipeps(1) % m / C_norm**(0.25D0)  
     ipeps(2) % m = ipeps(2) % m / C_norm**(0.25D0)   
  end IF

  !! Compute CTMRG boundaries of TN
  CALL compute_CTMRG_boundaries(Cmat, Tmat, TN2D, CTMRG, use_old_ctm) 

 END SUBROUTINE calc_CTM_environment_ipeps_GamLam





 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Construct reduced RHO for either plain or {Gamma, Lambda} iPEPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !!! By using environment of iPEPS unit cell construct reduced rho for calculating observables !!!
 SUBROUTINE construct_reduced_rho_PLAIN(rho_onesite, rho_twosite, Cmat, Tmat, ipeps, CTMRG, max_sep, opFlag)

  COMPLEX(KIND=DP),        ALLOCATABLE,  INTENT(INOUT)        :: rho_onesite(:,:)       !reduced rho matrix one-site
  COMPLEX(KIND=DP),        ALLOCATABLE,  INTENT(INOUT)        :: rho_twosite(:,:,:,:,:) !reduced rho matrix two-site
  TYPE(ctm_corner_type),                 INTENT(IN)           :: Cmat(:)                !CTM Cmat
  TYPE(ctm_transfer_type),               INTENT(IN)           :: Tmat(:)                !CTM Tmat
  TYPE(block_peps),                      INTENT(IN)           :: ipeps(:)               !iPEPS
  TYPE(ctm_params),                      INTENT(IN)           :: CTMRG                  !CTM params  
  INTEGER,                               INTENT(IN)           :: max_sep                !max separation for corr functions
  CHARACTER(LEN=*),                      INTENT(IN), OPTIONAL :: opFlag


  !! Bulk 2D network (TN2D)
  TYPE(block_mpo),    ALLOCATABLE :: TN2D(:)

  !! Transfer-MPO network tensors (Rvec, Lvec, Tmpo)
  TYPE(transfer_mpo), ALLOCATABLE :: Tmpo(:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: Rvec(:,:,:), Lvec(:,:,:)
  COMPLEX(KIND=DP)                :: C_norm

  !! Create reduced TN
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)

  !! Construct {Rvec, Lvec, Tmpo} from CTM tensors && 2D network
  CALL construct_Rvec_Tmpo_Lvec_network(Rvec, Tmpo, Lvec, Cmat, Tmat, TN2D, CTMRG)

  !! Calculate reduced rho (one-site and two-site) by contracting the above {Rvec, Lvec, Tmpo} network 
  CALL construct_rho_onesite(rho_onesite, ipeps, Rvec, Tmpo, Lvec, opFlag=OptArg(opFlag, 'op1'))
  CALL construct_rho_twosite(rho_twosite, ipeps, Rvec, Tmpo, Lvec, max_sep) 

 END SUBROUTINE construct_reduced_rho_PLAIN



 !!! By using environment of iPEPS unit cell construct reduced rho for calculating observables !!!
 SUBROUTINE construct_reduced_rho_GamLam(rho_onesite, rho_twosite, Cmat, Tmat, iGamma, iLambda, CTMRG, max_sep, opFlag)

  COMPLEX(KIND=DP),        ALLOCATABLE,  INTENT(INOUT)        :: rho_onesite(:,:)       !reduced rho matrix one-site
  COMPLEX(KIND=DP),        ALLOCATABLE,  INTENT(INOUT)        :: rho_twosite(:,:,:,:,:) !reduced rho matrix two-site
  TYPE(ctm_corner_type),                 INTENT(IN)           :: Cmat(:)                !CTM Cmat
  TYPE(ctm_transfer_type),               INTENT(IN)           :: Tmat(:)                !CTM Tmat
  TYPE(block_peps),                      INTENT(IN)           :: iGamma(:)              !iPEPS Gamma
  TYPE(block_lambda),                    INTENT(IN)           :: iLambda(:)             !iPEPS Lambda
  TYPE(ctm_params),                      INTENT(IN)           :: CTMRG                  !CTM params  
  INTEGER,                               INTENT(IN)           :: max_sep                !max separation for corr functions
  CHARACTER(LEN=*),                      INTENT(IN), OPTIONAL :: opFlag

  !! Bulk 2D network (TN2D)
  TYPE(block_peps),   ALLOCATABLE :: ipeps(:)
  TYPE(block_mpo),    ALLOCATABLE :: TN2D(:)

  !! Transfer-MPO network tensors (Rvec, Lvec, Tmpo)
  TYPE(transfer_mpo), ALLOCATABLE :: Tmpo(:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: Rvec(:,:,:), Lvec(:,:,:)
  COMPLEX(KIND=DP)                :: C_norm

  !! Create reduced TN
  CALL absorb_ipeps_lambdas(ipeps, iGamma, iLambda)
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)

  !! Construct {Rvec, Lvec, Tmpo} from CTM tensors && 2D network
  CALL construct_Rvec_Tmpo_Lvec_network(Rvec, Tmpo, Lvec, Cmat, Tmat, TN2D, CTMRG)

  !! Calculate reduced rho (one-site and two-site) by contracting the above {Rvec, Lvec, Tmpo} network 
  CALL construct_rho_onesite(rho_onesite, ipeps, Rvec, Tmpo, Lvec, opFlag=OptArg(opFlag, 'op1'))
  CALL construct_rho_twosite(rho_twosite, ipeps, Rvec, Tmpo, Lvec, max_sep) 

 END SUBROUTINE construct_reduced_rho_GamLam




 !!! Contruct reduced one-site rho !!!
 SUBROUTINE construct_rho_onesite(rho_onesite, ipeps, Rvec, Tmpo, Lvec, opFlag)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT)        :: rho_onesite(:,:)
  TYPE(block_peps),              INTENT(IN)           :: ipeps(:)
  COMPLEX(KIND=DP),              INTENT(IN)           :: Rvec(:,:,:), Lvec(:,:,:) 
  TYPE(transfer_mpo),            INTENT(IN)           :: Tmpo(:) 
  CHARACTER(LEN=*),              INTENT(IN), OPTIONAL :: opFlag

  !! Free operator sites of reduced rho
  TYPE(op_sites_type)           :: TN_rhoSites
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_tmp(:)
  INTEGER                       :: LocalDim, op_dim
  INTEGER                       :: i

  !! Get LocalDim, allocate rho_tmp
  LocalDim = ipeps(1) % LocalDim
  op_dim   = SQRTint(LocalDim) 

  CALL allocateTens(rho_tmp, (/ LocalDim /))
  
  !! Compute individual elements of rho_tmp one by one
  DO i=1,LocalDim
     CALL create_reduced_rho_sites_from_ipeps(TN_rhoSites, ipeps, i)
     CALL contract_TN_1P_obs(rho_tmp(i), Rvec, Tmpo, Lvec, TN_rhoSites, opFlag=OptArg(opFlag, 'op1'))
  end DO

  !! Reshape rho into matrix
  rho_onesite = RESHAPE_2D(rho_tmp, (/ op_dim, op_dim /))

 END SUBROUTINE construct_rho_onesite




 !!! Contruct reduced two-site rho !!!
 SUBROUTINE construct_rho_twosite(rho_twosite, ipeps, Rvec, Tmpo, Lvec, max_sep)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: rho_twosite(:,:,:,:,:)
  TYPE(block_peps),              INTENT(IN)    :: ipeps(:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: Rvec(:,:,:), Lvec(:,:,:) 
  TYPE(transfer_mpo),            INTENT(IN)    :: Tmpo(:)     
  INTEGER,                       INTENT(IN)    :: max_sep

  !! Free operator sites of reduced rho
  TYPE(op_sites_type)           :: TN_rhoSites
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_tmp(:,:,:)
  INTEGER                       :: LocalDim, op_dim
  INTEGER                       :: i, j, sep

  !! No calculation if max_sep = 0
  IF(max_sep .EQ. 0) RETURN

  !! Get LocalDim, allocate rho_tmp, excluding sep=0 (which we'll treat separately) -- set min_sep=1
  LocalDim = ipeps(1) % LocalDim
  op_dim   = SQRTint(LocalDim)  

  CALL allocateTens(rho_tmp, (/max_sep+1, LocalDim, LocalDim/))
  CALL allocateTens(rho_twosite, (/max_sep+1, op_dim, op_dim, op_dim, op_dim/))

  !! Compute individual elements of rho_twosite (for a range of separations)
  DO i=1,LocalDim
    DO j=1,LocalDim
       CALL create_reduced_rho_sites_from_ipeps(TN_rhoSites, ipeps, i, j)
       CALL contract_TN_2P_corr(rho_tmp(:,i,j), Rvec, Tmpo, Lvec, TN_rhoSites, min_sep=1, max_sep=max_sep)  
    end DO
  end DO

  !! Reshape rho into matrix
  DO sep=1,max_sep
     rho_twosite(sep+1,:,:,:,:) = RESHAPE_4D(rho_tmp(sep+1,:,:), '12,34', (/op_dim, op_dim/), (/op_dim, op_dim/))
  end DO

 END SUBROUTINE construct_rho_twosite

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPUTING OBSERVABLES DURING/AFTER TIME EVOLUTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute 1P observable !!!
 SUBROUTINE compute_obs_1P(expval, rho_onesite, op)

  COMPLEX(KIND=DP), INTENT(INOUT) :: expval
  COMPLEX(KIND=DP), INTENT(IN)    :: rho_onesite(:,:), op(:,:)

  !! Compute 1P observable
  expval = TRACE(TENSMUL(op, rho_onesite))

 END SUBROUTINE compute_obs_1P



 !!! Compute 2P observable !!!
 SUBROUTINE compute_obs_2P(expvals, rho_twosite, rho_onesite, ops, min_sep, max_sep)

  COMPLEX(KIND=DP), INTENT(INOUT) :: expvals(:)
  COMPLEX(KIND=DP), INTENT(IN)    :: rho_twosite(:,:,:,:,:), rho_onesite(:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: ops(:,:,:)
  INTEGER,          INTENT(IN)    :: min_sep, max_sep

  COMPLEX(KIND=DP), ALLOCATABLE   :: tmpRho(:,:)
  INTEGER                         :: sep

  !! Compute 2P observables (for different separations) -- sep=0 handled separately as onesite case
  DO sep = min_sep, max_sep

     IF(sep .EQ. 0) THEN 
        expvals(sep+1-min_sep) = TRACE(TENSMUL(ops(1,:,:), TENSMUL(ops(2,:,:),  rho_onesite)))
     ELSE
        tmpRho                 = TRACE(TENSMUL(rho_twosite(sep+1,:,:,:,:),  ops(1,:,:),  '12'),  '12')
        expvals(sep+1-min_sep) = TRACE(TENSMUL(ops(2,:,:), tmpRho))
     end IF

  end DO

 END SUBROUTINE compute_obs_2P


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPUTE LIOUVILLIAN ENERGY PER SITE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute Liouvillian energy per site !!!
 SUBROUTINE compute_energy_per_site(liouvillian_twosite, iGamma, iLambda, Cmat, Tmat, CTMRG, x, datafile, dump_parameters)

  USE datafile_utility

  COMPLEX(KIND=DP),        INTENT(IN)  :: liouvillian_twosite(:,:)
  TYPE(block_peps),        INTENT(IN)  :: iGamma(:)
  TYPE(block_lambda),      INTENT(IN)  :: iLambda(:)
  TYPE(ctm_corner_type),   INTENT(IN)  :: Cmat(:)    
  TYPE(ctm_transfer_type), INTENT(IN)  :: Tmat(:)  
  TYPE(ctm_params),        INTENT(IN)  :: CTMRG   
  INTEGER,                 INTENT(IN)  :: datafile     
  REAL(KIND=DP),           INTENT(IN)  :: x

  INTERFACE
     SUBROUTINE dump_parameters(HANDLE)
         INTEGER, INTENT(IN) :: HANDLE
     END SUBROUTINE dump_parameters
  end INTERFACE  

  !! Reduced rho && norm of rho TN (i.e. trace)
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_onesite(:,:), rho_twosite(:,:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_vec(:)
  COMPLEX(KIND=DP)              :: norm
  COMPLEX(KIND=DP)              :: energy

  !! Construct two-site rho with sep=1
  INTEGER,          PARAMETER   :: sep = 1

  !! Setup reduced rho
  CALL construct_reduced_rho(rho_onesite, rho_twosite, Cmat, Tmat, iGamma, iLambda, CTMRG, max_sep=sep)

  !! Construct rho vector from a two-site reduced rho with sep=1
  rho_vec = RESHAPE_1D(RESHAPE_2D(rho_twosite(sep+1,:,:,:,:), '12,34'))

  !! Get norm
  norm = SUM(CONJG(rho_vec) * rho_vec) 

  !! Contract <rho|L_{twosite}|rho>
  energy = SUM(CONJG(rho_vec) * TENSMUL(liouvillian_twosite, rho_vec))

  !! Print expval
  CALL write_summary_data_1P(datafile, x, energy/norm, norm)

 END SUBROUTINE compute_energy_per_site


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CREATE REDUCED NETWORK FROM RHO iPEPS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Create MPO representing reduced iPEPS network !!!
 SUBROUTINE create_reduced_TN_from_ipeps(TN2D, ipeps)

  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT) :: TN2D(:)
  TYPE(block_peps),                INTENT(IN)    :: ipeps(:)

  !! Indices && dims
  INTEGER :: site, N_sites

  !! TN size
  N_sites = SIZE(ipeps)

  !! Construct sites of reduced iPEPS network
  CALL allocate_empty_mpo_block(TN2D, 2)
  DO site=1,N_sites
     CALL allocate_mpo_site(TN2D(site), rho_tmat_site(ipeps(site)))
  end DO

 END SUBROUTINE create_reduced_TN_from_ipeps




 !!! Create operator sites of reduced iPEPS network !!!
 SUBROUTINE create_reduced_op_sites_from_ipeps(TN_OpSites, ipeps, op1, op2)

  TYPE(op_sites_type),             INTENT(INOUT)        :: TN_OpSites
  TYPE(block_peps),                INTENT(IN)           :: ipeps(:)
  COMPLEX(KIND=DP),                INTENT(IN)           :: op1(:,:)
  COMPLEX(KIND=DP),                INTENT(IN), OPTIONAL :: op2(:,:)
  
  !! Construct operator sites
  TN_OpSites % site1_Op1        = rho_tmat_site(ipeps(1), op1)

  IF(PRESENT(op2)) THEN
     TN_OpSites % site1_Op1_Op2 = rho_tmat_site(ipeps(1), op1, op2)
     TN_OpSites % site1_Op2     = rho_tmat_site(ipeps(1), op2)
     TN_OpSites % site2_Op2     = rho_tmat_site(ipeps(2), op2)
  end IF

 END SUBROUTINE create_reduced_op_sites_from_ipeps





 !!! [Reduced-rho iPEPS network] --> create untraced RHO_vec sites !!!
 SUBROUTINE create_reduced_rho_sites_from_ipeps(TN_OpSites, ipeps, i, j)

  TYPE(op_sites_type),             INTENT(INOUT)        :: TN_OpSites
  TYPE(block_peps),                INTENT(IN)           :: ipeps(:)
  INTEGER,                         INTENT(IN)           :: i
  INTEGER,                         INTENT(IN), OPTIONAL :: j

  !! Compute reduced TN sites-1,2 -- bare tensors, no operator applied, just rho traced out
  CALL copyTens(TN_OpSites % site1_Op1,         ipeps(1) % m(i,:,:,:,:)) 

  IF(PRESENT(j)) THEN
     CALL copyTens(TN_OpSites % site1_Op1_Op2,  ipeps(1) % m(j,:,:,:,:))
     CALL copyTens(TN_OpSites % site1_Op2,      ipeps(1) % m(j,:,:,:,:))
     CALL copyTens(TN_OpSites % site2_Op2,      ipeps(2) % m(j,:,:,:,:))
  end IF

 END SUBROUTINE create_reduced_rho_sites_from_ipeps





 !!! Construct a single site of reduced TN -- traced rho network !!!
 FUNCTION rho_tmat_site(ipeps, op, opP)

  TYPE(block_peps), INTENT(IN)           :: ipeps                  !ipeps site (nolambda)
  COMPLEX(KIND=DP), INTENT(IN), OPTIONAL :: op(:,:), opP(:,:)      !operators applied to the site (in case we apply two operators)
  COMPLEX(KIND=DP), ALLOCATABLE          :: rho_tmat_site(:,:,:,:) !result

  !! Local copy of iPEPS, and eye matrix
  COMPLEX(KIND=DP), ALLOCATABLE :: tempSite(:,:,:,:,:), eye(:,:)

  !! Create local copy of iPEPS site
  CALL copyTens(tempSite, ipeps % m)

  !! Create eye matrix
  IF(PRESENT(op) .OR. PRESENT(opP)) eye = matEye(SQRTint(ipeps % LocalDim))

  !! Apply op, opP, or both (if provided)
  IF(PRESENT(op))  tempSite = TENSMUL(tempSite, TensKRON(op,  eye), '12')
  IF(PRESENT(opP)) tempSite = TENSMUL(tempSite, TensKRON(opP, eye), '12')

  !! Take trace
  !! CALL copyTens(rho_tmat_site, rho_trace(tempSite))
  CALL copyTens(rho_tmat_site, rho_trace(tempSite))

 END FUNCTION rho_tmat_site





 !!! Calculate trace of a vectorized rho matrix !!!
 FUNCTION rho_trace(rho_site)

  COMPLEX(KIND=DP), INTENT(IN)  :: rho_site(:,:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_trace(:,:,:,:)

  INTEGER :: op_dim, i

  !! Determine operator dim
  op_dim = SQRTint(SIZE(rho_site, 1))

  !! Allocate rho_trace with identical bond dims to rho_site
  CALL allocateTens(rho_trace, shape(rho_site(1,:,:,:,:)))

  !! Take trace
  rho_trace = 0.0D0
  DO i=1,op_dim 
     rho_trace = rho_trace + rho_site(ICOM(i,i,op_dim),:,:,:,:)
  end DO

 END FUNCTION rho_trace

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Contract bare iPEPS-TN (i.e. no operators) embedded in CTM environment !!!
 SUBROUTINE contract_bare_TN(expval, Cmat, Tmat, ipeps, CTMRG, DIR)

  COMPLEX(KIND=DP),         INTENT(INOUT) :: expval    !expval = result of contraction
  TYPE(ctm_corner_type),    INTENT(IN)    :: Cmat(:)   !CTM Cmat
  TYPE(ctm_transfer_type),  INTENT(IN)    :: Tmat(:)   !CTM Tmat
  TYPE(block_peps),         INTENT(IN)    :: ipeps(:)  !iPEPS
  TYPE(ctm_params),         INTENT(IN)    :: CTMRG     !CTM params
  CHARACTER(LEN=*),         INTENT(IN)    :: DIR

  !! Bulk 2D network (TN2D)
  TYPE(block_mpo),    ALLOCATABLE :: TN2D(:)
  TYPE(block_mps),    ALLOCATABLE :: ENV2D(:)

  !! Transfer-MPO network tensors (Rvec, Lvec, Tmpo)
  TYPE(transfer_mpo), ALLOCATABLE :: Tmpo(:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: Rvec(:,:,:), Lvec(:,:,:)

  !! Indices && dims
  INTEGER :: s, sA, sB
  INTEGER :: xA, xB, yA, yB
  INTEGER, PARAMETER :: N_sites = 2

  !! Create reduced TN
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)

  !! Calculate single-bond environment of iPEPS
  CALL calc_single_bond_environment(ENV2D, Cmat, Tmat, TN2D, CTMRG, DIR=DIR, pbc=.FALSE.)

  !! Get iPEPS sites
  CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, DIR)

  SELECT CASE(DIR)
  CASE('S', 'N')
       CONTINUE
  CASE('W', 'E')
       DO s=1,N_sites
          CALL allocate_mpo_site(TN2D(s), TensROTATE(TN2D(s) % m, '+PI/2'))
       end DO
  CASE DEFAULT
       CALL invalid_flag("contract_bare_TN -- invalid DIR ", DIR)
  end SELECT

  !!! Construct transfer MPO from ENV2D && TN2D
  CALL allocate_empty_transfer_mpo(Tmpo, N_sites)
  DO s=1,N_sites
     CALL copyTens(Tmpo(s) % U,  ENV2D(s)   % m)
     CALL copyTens(Tmpo(s) % D,  ENV2D(s+2) % m)
  end DO

  CALL copyTens(Tmpo(1) % M,  TN2D(sA)    % m)
  CALL copyTens(Tmpo(2) % M,  TN2D(sB)    % m)

  !!! Construct boundary vectors from ENV2D
  CALL copyTens(Lvec,  ENV2D(5) % m)
  CALL copyTens(Rvec,  ENV2D(6) % m)

  !! Contract TN
  CALL mult_bvec_tmpo_bvec(expval, Rvec, Tmpo, Lvec)

 END SUBROUTINE contract_bare_TN






 !!! Compute environment of iPEPS unit cell, normalize by the TN contraction expval !!!
 SUBROUTINE normalize_ipeps_by_simple_trace(iGamma, iLambda, Cmat, Tmat, use_old_ctm)

  TYPE(block_peps),                     INTENT(INOUT) :: iGamma(:)              !iPEPS gamma
  TYPE(block_lambda),                   INTENT(IN)    :: iLambda(:)             !iPEPS lambda
  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)                !CTM Cmat
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)                !CTM Tmat
  LOGICAL,                              INTENT(IN)    :: use_old_ctm            !whether to use old CTM for initialization

  !! Bulk 2D network (TN2D)
  TYPE(block_peps),   ALLOCATABLE :: ipeps(:)
  TYPE(block_mpo),    ALLOCATABLE :: TN2D(:)

  !! CTM params
  TYPE(ctm_params) :: CTMRG  
  COMPLEX(KIND=DP) :: C_norm   

  !! Create a new instance of CTM
  CALL create_CTM_params(CTMRG, eps=2.0D0, chi=1, eta=1.0D-06, ctm_method='3')             

  !! Create reduced TN
  CALL absorb_ipeps_lambdas(ipeps, iGamma, iLambda)
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)

  !! Compute CTMRG boundaries of TN
  CALL compute_CTMRG_boundaries(Cmat, Tmat, TN2D, CTMRG, use_old_ctm) 
  CALL contract_CTM_network(C_norm, Cmat, Tmat, TN2D)

  iGamma(1) % m = iGamma(1) % m / C_norm**(0.25D0)
  iGamma(2) % m = iGamma(2) % m / C_norm**(0.25D0)

 END SUBROUTINE normalize_ipeps_by_simple_trace




 !!! Impose Hermiticity on RHO iPEPS !!!
 SUBROUTINE hermitize_rho_ipeps(iGamma, iLambda)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda),              INTENT(IN)    :: iLambda(:)

  TYPE(block_peps),   ALLOCATABLE :: ipeps(:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: rhoAA(:,:,:,:,:,:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: rhoBB(:,:,:,:,:,:)
  INTEGER                         :: dimHS

  INTEGER :: sA, sB
  INTEGER :: xA, xB, yA, yB

  !! Absorb iPEPS lambdas
  CALL absorb_ipeps_lambdas(ipeps, iGamma, iLambda)

  !! Get Hilbert dim
  dimHS = SQRTint(ipeps(1) % LocalDim)

  !! Reshape RHO iPEPS sites from vector-like into matrix-like format
  rhoAA = RESHAPE_6D(ipeps(1) % m, '12,3,4,5,6', (/dimHS, dimHS/))
  rhoBB = RESHAPE_6D(ipeps(2) % m, '12,3,4,5,6', (/dimHS, dimHS/))
  
  !! Impose Hermiticity
  rhoAA = 0.5D0 * (rhoAA + TensTRANSPOSE(CONJG(rhoAA), '12'))
  rhoBB = 0.5D0 * (rhoBB + TensTRANSPOSE(CONJG(rhoBB), '12')) 
  
  !! Reshape back, copy to the original RHO iPEPS
  CALL allocate_peps_site(iGamma(1),  RESHAPE_5D(rhoAA, '12,3,4,5,6')) 
  CALL allocate_peps_site(iGamma(2),  RESHAPE_5D(rhoBB, '12,3,4,5,6')) 

  !! Restore iPEPS lambdas
  CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, 'S')

  CALL LambdaDIV(iGamma(sA) % m, SQRT(iLambda(xA) % m), '5')
  CALL LambdaDIV(iGamma(sB) % m, SQRT(iLambda(xA) % m), '4')

  CALL LambdaDIV(iGamma(sA) % m, SQRT(iLambda(xB) % m), '4')
  CALL LambdaDIV(iGamma(sB) % m, SQRT(iLambda(xB) % m), '5')

  CALL LambdaDIV(iGamma(sA) % m, SQRT(iLambda(yA) % m), '3')
  CALL LambdaDIV(iGamma(sB) % m, SQRT(iLambda(yA) % m), '2')

  CALL LambdaDIV(iGamma(sA) % m, SQRT(iLambda(yB) % m), '2')
  CALL LambdaDIV(iGamma(sB) % m, SQRT(iLambda(yB) % m), '3')

 END SUBROUTINE hermitize_rho_ipeps




 !!! Calculate simplified observables without environment (i.e. assuming MF environment) !!!
 SUBROUTINE contract_simplified_obs1P(expval, iGamma, iLambda, op)

  COMPLEX(KIND=DP),   INTENT(INOUT) :: expval
  TYPE(block_peps),   INTENT(IN)    :: iGamma(:)
  TYPE(block_lambda), INTENT(IN)    :: iLambda(:)
  COMPLEX(KIND=DP),   INTENT(IN)    :: op(:,:)

  COMPLEX(KIND=DP) :: C_norm

  !! Contract 2D network --> get expval && norm
  CALL contract_simplified_network(expval, iGamma, iLambda, op)
  CALL contract_simplified_network(C_norm, iGamma, iLambda)

  !! Obtain normalized expval
  expval = expval/C_norm

 END SUBROUTINE contract_simplified_obs1P




 !!! Simplified contraction of 2D network without environment (i.e. assuming MF environment) !!!
 SUBROUTINE contract_simplified_network(expval, iGamma, iLambda, op)

  COMPLEX(KIND=DP),   INTENT(INOUT)        :: expval
  TYPE(block_peps),   INTENT(IN)           :: iGamma(:)
  TYPE(block_lambda), INTENT(IN)           :: iLambda(:)
  COMPLEX(KIND=DP),   INTENT(IN), OPTIONAL :: op(:,:)

  !! Plain iPEPS && single iPEPS site with operator
  TYPE(block_peps), ALLOCATABLE :: ipeps(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tmp(:,:,:,:,:,:)

  !! Reduced tensors of 2D network
  COMPLEX(KIND=DP), ALLOCATABLE :: AO(:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: A(:,:,:,:), B(:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: AB(:,:),    BA(:,:)

  !! Dims
  INTEGER :: localDim, opDim

  !! Absorb iPEPS lambdas
  CALL absorb_ipeps_lambdas(ipeps, iGamma, iLambda)

  !! Trace out iPEPS sites
  A = rho_tmat_site(ipeps(1))
  B = rho_tmat_site(ipeps(2))

  !! If operator is provided, apply it to iPEPS --> create reduced tens with op
  IF(PRESENT(op)) THEN

     !! Unvectorize rho site |rho> -> rho
     localDim = ipeps(1) % LocalDim
     opDim    = SQRTint(localDim) 
     tmp = RESHAPE_6D(ipeps(1) % m, '12,3,4,5,6', (/ opDim, opDim /))

     !! Apply op to rho matrix && take trace: AO = Tr(O*rho)
     AO  = TRACE(TENSMUL(tmp, op, MULT='12'), '12')
  ELSE
     CALL copyTens(AO, A)
  end IF

  !! Contract && trace 2x2 block of reduced TN sites (assuming MF environment)
  !! (including AO -- the site that can contain operator)
  AB     = TRACE(TENSMUL(AO, B, MULT='43', FUSE='(11,22)'), '34')
  BA     = TRACE(TENSMUL(B,  A, MULT='43', FUSE='(11,22)'), '34')
  expval = TRACE(TENSMUL(AB, BA))
  
 END SUBROUTINE contract_simplified_network

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE observables_ipeps
