MODULE imaginary_observables_ipeps

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

 interface resize_ipeps
    module procedure resize_ipeps_GamLam
    module procedure resize_ipeps_PLAIN
 end interface resize_ipeps

 interface compute_energy_per_site
    module procedure compute_energy_per_site_GamLam
    module procedure compute_energy_per_site_PLAIN
 end interface compute_energy_per_site

 private calc_CTM_environment_ipeps_GamLam,      calc_CTM_environment_ipeps_PLAIN
 private construct_reduced_rho_GamLam,           construct_reduced_rho_PLAIN
 private resize_ipeps_GamLam,                    resize_ipeps_PLAIN
 private compute_energy_per_site_GamLam,         compute_energy_per_site_PLAIN

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
  COMPLEX(KIND=DP)                :: C_norm

  !! If we're using previous CTM, normalize TN first
  IF(ALLOCATED(Cmat) .AND. ALLOCATED(Tmat)) THEN

     CALL create_reduced_TN_from_ipeps(TN2D, ipeps)
     CALL resize_transfer_mat(Tmat, TN2D)
     CALL contract_CTM_network(C_norm, Cmat, Tmat, TN2D)

     ipeps(1) % m = ipeps(1) % m / ABS(C_norm)**(0.125D0)
     ipeps(2) % m = ipeps(2) % m / ABS(C_norm)**(0.125D0)
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
  COMPLEX(KIND=DP)                :: C_norm

  !! If we're using previous CTM, normalize TN first
  IF(ALLOCATED(Cmat) .AND. ALLOCATED(Tmat) .AND. use_old_ctm) THEN

     CALL absorb_ipeps_lambdas(ipeps, iGamma, iLambda)
     CALL create_reduced_TN_from_ipeps(TN2D, ipeps)
     CALL contract_CTM_network(C_norm, Cmat, Tmat, TN2D)

     iGamma(1) % m = iGamma(1) % m / ABS(C_norm)**(0.125D0)
     iGamma(2) % m = iGamma(2) % m / ABS(C_norm)**(0.125D0)
  end IF

  !! Create reduced TN
  CALL absorb_ipeps_lambdas(ipeps, iGamma, iLambda)
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)

  !! Compute CTMRG boundaries of TN
  CALL compute_CTMRG_boundaries(Cmat, Tmat, TN2D, CTMRG, use_old_ctm) 

 END SUBROUTINE calc_CTM_environment_ipeps_GamLam

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Construct reduced RHO for either plain or {Gamma, Lambda} iPEPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! By using environment of iPEPS unit cell construct reduced rho for calculating observables !!!
 SUBROUTINE construct_reduced_rho_PLAIN(rho_onesite, rho_twosite, Cmat, Tmat, ipeps, CTMRG, max_sep)

  COMPLEX(KIND=DP),        ALLOCATABLE, INTENT(INOUT) :: rho_onesite(:,:)       !reduced rho matrix one-site
  COMPLEX(KIND=DP),        ALLOCATABLE, INTENT(INOUT) :: rho_twosite(:,:,:,:,:) !reduced rho matrix two-site
  TYPE(ctm_corner_type),                INTENT(IN)    :: Cmat(:)                !CTM Cmat
  TYPE(ctm_transfer_type),              INTENT(IN)    :: Tmat(:)                !CTM Tmat
  TYPE(block_peps),                     INTENT(IN)    :: ipeps(:)               !iPEPS
  TYPE(ctm_params),                     INTENT(IN)    :: CTMRG                  !CTM params
  INTEGER,                              INTENT(IN)    :: max_sep                !max separation for corr functions

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
  CALL construct_rho_onesite(rho_onesite, ipeps, Rvec, Tmpo, Lvec)
  CALL construct_rho_twosite(rho_twosite, ipeps, Rvec, Tmpo, Lvec, max_sep) 

 END SUBROUTINE construct_reduced_rho_PLAIN





 !!! By using environment of iPEPS unit cell construct reduced rho for calculating observables !!!
 SUBROUTINE construct_reduced_rho_GamLam(rho_onesite, rho_twosite, Cmat, Tmat, iGamma, iLambda, CTMRG, max_sep)

  COMPLEX(KIND=DP),        ALLOCATABLE, INTENT(INOUT) :: rho_onesite(:,:)       !reduced rho matrix one-site
  COMPLEX(KIND=DP),        ALLOCATABLE, INTENT(INOUT) :: rho_twosite(:,:,:,:,:) !reduced rho matrix two-site
  TYPE(ctm_corner_type),                INTENT(IN)    :: Cmat(:)                !CTM Cmat
  TYPE(ctm_transfer_type),              INTENT(IN)    :: Tmat(:)                !CTM Tmat
  TYPE(block_peps),                     INTENT(IN)    :: iGamma(:)              !iPEPS Gamma
  TYPE(block_lambda),                   INTENT(IN)    :: iLambda(:)             !iPEPS Lambda
  TYPE(ctm_params),                     INTENT(IN)    :: CTMRG                  !CTM params
  INTEGER,                              INTENT(IN)    :: max_sep                !max separation for corr functions

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
  CALL construct_rho_onesite(rho_onesite, ipeps, Rvec, Tmpo, Lvec)
  CALL construct_rho_twosite(rho_twosite, ipeps, Rvec, Tmpo, Lvec, max_sep) 

 END SUBROUTINE construct_reduced_rho_GamLam




 !!! Contruct reduced one-site rho !!!
 SUBROUTINE construct_rho_onesite(rho_onesite, ipeps, Rvec, Tmpo, Lvec)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: rho_onesite(:,:)
  TYPE(block_peps),              INTENT(IN)    :: ipeps(:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: Rvec(:,:,:), Lvec(:,:,:) 
  TYPE(transfer_mpo),            INTENT(IN)    :: Tmpo(:) 

  !! Free operator sites of reduced rho
  TYPE(op_sites_type) :: TN_psiSites
  INTEGER             :: i, j, LocalDim

  !! Get LocalDim, allocate rho_onesite
  LocalDim = ipeps(1) % LocalDim
  CALL allocateTens(rho_onesite, (/ LocalDim, LocalDim /))
  
  !! Compute individual elements of rho_onesite one by one
  DO i=1,LocalDim
    DO j=1,LocalDim
       CALL create_reduced_psi_sites_from_ipeps(TN_psiSites,  ipeps, i, j)
       CALL contract_TN_1P_obs(rho_onesite(i, j), Rvec, Tmpo, Lvec, TN_psiSites)       
    end DO
  end DO

 END SUBROUTINE construct_rho_onesite






 !!! Contruct reduced two-site rho !!!
 SUBROUTINE construct_rho_twosite(rho_twosite, ipeps, Rvec, Tmpo, Lvec, max_sep)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: rho_twosite(:,:,:,:,:)
  TYPE(block_peps),              INTENT(IN)    :: ipeps(:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: Rvec(:,:,:), Lvec(:,:,:) 
  TYPE(transfer_mpo),            INTENT(IN)    :: Tmpo(:)      
  INTEGER,                       INTENT(IN)    :: max_sep

  !! Free operator sites of reduced rho
  TYPE(op_sites_type) :: TN_psiSites
  INTEGER             :: i, j, k, l, LocalDim, sep

  !! Get LocalDim, allocate rho_twosite, excluding sep=0 (which we'll treat separately) -- set min_sep=1
  LocalDim = ipeps(1) % LocalDim
  CALL allocateTens(rho_twosite, (/max_sep+1, LocalDim, LocalDim, LocalDim, LocalDim/)) 

  !! Compute individual elements of rho_twosite (for a range of separations)
  DO i=1,LocalDim
    DO j=1,LocalDim
      DO k=1,LocalDim
        DO l=1,LocalDim
           CALL create_reduced_psi_sites_from_ipeps(TN_psiSites, ipeps, i, j, k, l)
           CALL contract_TN_2P_corr(rho_twosite(:, i, j, k, l), Rvec, Tmpo, Lvec, TN_psiSites, min_sep=1, max_sep=max_sep)  
        end DO
      end DO
    end DO
  end DO

 END SUBROUTINE construct_rho_twosite

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPUTING OBSERVABLES DURING/AFTER IMAGINARY TIME EVOLUTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

  COMPLEX(KIND=DP), ALLOCATABLE :: tmpRho(:,:)
  INTEGER :: sep

  !! Compute 2P observables (for different separations) -- sep=0 handled separately as onesite case
  DO sep = min_sep, max_sep

     IF(sep .EQ. 0) THEN 
        expvals(sep+1-min_sep) = TRACE(TENSMUL(ops(1,:,:), TENSMUL(ops(2,:,:), rho_onesite)))
     ELSE
        tmpRho                 = TRACE(TENSMUL(rho_twosite(sep+1,:,:,:,:),  ops(1,:,:),  '12'),  '12') 
        expvals(sep+1-min_sep) = TRACE(TENSMUL(ops(2,:,:), tmpRho))
     end IF

  end DO

 END SUBROUTINE compute_obs_2P


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPUTE ENERGY PER SITE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute energy per site !!!
 SUBROUTINE compute_energy_per_site_GamLam(hamiltonian_twosite, iGamma, iLambda, Cmat, Tmat, CTMRG, x, datafile, dump_parameters)

  USE datafile_utility

  COMPLEX(KIND=DP),        INTENT(IN)  :: hamiltonian_twosite(:,:)
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
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_vec(:,:)
  COMPLEX(KIND=DP)              :: norm
  COMPLEX(KIND=DP)              :: energy

  !! Construct two-site rho with sep=1
  INTEGER,          PARAMETER   :: sep = 1

  !! Setup reduced rho
  CALL construct_reduced_rho(rho_onesite, rho_twosite, Cmat, Tmat, iGamma, iLambda, CTMRG, max_sep=sep)

  !! Get norm of rho TN
  norm = TRACE(rho_onesite)

  !! Construct rho vector from a two-site reduced rho with sep=1
  rho_vec = RESHAPE_2D(rho_twosite(sep+1,:,:,:,:), '12,34')

  !! Contract <Psi|H_{twosite}|Psi>
  energy = TRACE(TENSMUL(hamiltonian_twosite, rho_vec))

  !! Print expval
  CALL write_summary_data_1P(datafile, x, energy/norm, norm)

 END SUBROUTINE compute_energy_per_site_GamLam




 !!! Compute energy per site !!!
 SUBROUTINE compute_energy_per_site_PLAIN(hamiltonian_twosite, ipeps, Cmat, Tmat, CTMRG, x, datafile, dump_parameters)

  USE datafile_utility

  COMPLEX(KIND=DP),        INTENT(IN)  :: hamiltonian_twosite(:,:)
  TYPE(block_peps),        INTENT(IN)  :: ipeps(:)
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
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_vec(:,:)
  COMPLEX(KIND=DP)              :: norm
  COMPLEX(KIND=DP)              :: energy

  !! Construct two-site rho with sep=1
  INTEGER,          PARAMETER   :: sep = 1

  !! Setup reduced rho
  CALL construct_reduced_rho(rho_onesite, rho_twosite, Cmat, Tmat, ipeps, CTMRG, max_sep=sep)

  !! Get norm of rho TN
  norm = TRACE(rho_onesite)

  !! Construct rho vector from a two-site reduced rho with sep=1
  rho_vec = RESHAPE_2D(rho_twosite(sep+1,:,:,:,:), '12,34')

  !! Contract <Psi|H_{twosite}|Psi>
  energy = TRACE(TENSMUL(hamiltonian_twosite, rho_vec))

  !! Print expval
  CALL write_summary_data_1P(datafile, x, energy/norm, norm)

 END SUBROUTINE compute_energy_per_site_PLAIN


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CREATE REDUCED NETWORK FROM PSI iPEPS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Create MPO representing reduced iPEPS network !!!
 SUBROUTINE create_reduced_TN_from_ipeps(TN2D, ipeps)

  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT)        :: TN2D(:)
  TYPE(block_peps),                INTENT(IN)           :: ipeps(:)

  !! Indices && dims
  INTEGER :: site, N_sites

  !! TN size
  N_sites = SIZE(ipeps)

  !! Construct sites of reduced iPEPS network
  CALL allocate_empty_mpo_block(TN2D, 2)
  DO site=1,N_sites
     CALL allocate_mpo_site(TN2D(site), psi_tmat_site(ipeps(site)))
  end DO

 END SUBROUTINE create_reduced_TN_from_ipeps




 !!! Create operator sites of reduced iPEPS network !!!
 SUBROUTINE create_reduced_op_sites_from_ipeps(TN_OpSites, ipeps, op1, op2)

  TYPE(op_sites_type),             INTENT(INOUT)        :: TN_OpSites
  TYPE(block_peps),                INTENT(IN)           :: ipeps(:)
  COMPLEX(KIND=DP),                INTENT(IN)           :: op1(:,:)
  COMPLEX(KIND=DP),                INTENT(IN), OPTIONAL :: op2(:,:)
  
  !! Construct operator sites
  TN_OpSites % site1_Op1        = psi_tmat_site(ipeps(1), op1)

  IF(PRESENT(op2)) THEN
     TN_OpSites % site1_Op1_Op2 = psi_tmat_site(ipeps(1), op1, op2)
     TN_OpSites % site1_Op2     = psi_tmat_site(ipeps(1), op2)
     TN_OpSites % site2_Op2     = psi_tmat_site(ipeps(2), op2)
  end IF

 END SUBROUTINE create_reduced_op_sites_from_ipeps





 !!! [Reduced-rho iPEPS network] --> create untraced |PSI><PSI| sites !!!
 SUBROUTINE create_reduced_psi_sites_from_ipeps(TN_OpSites, ipeps, i, j, k, l)

  TYPE(op_sites_type),             INTENT(INOUT)        :: TN_OpSites
  TYPE(block_peps),                INTENT(IN)           :: ipeps(:)
  INTEGER,                         INTENT(IN)           :: i, j
  INTEGER,                         INTENT(IN), OPTIONAL :: k, l

  !! Local HC copy
  TYPE(block_peps), ALLOCATABLE :: ipepsHC(:)

  !! Create a local HC copy of iPEPS
  CALL copy_peps_block(ipepsHC, ipeps, 'HC')

  !! Compute reduced TN sites-1,2 -- bare tensors, no operator applied, just rho traced out 
  !! (psi x psiHC --> consistent with psi_tmat_site())
  TN_OpSites % site1_Op1 = TensKRON(ipeps(1) % m(i,:,:,:,:), ipepsHC(1) % m(j,:,:,:,:))

  IF(PRESENT(k) .AND. PRESENT(l)) THEN

     CALL copyTens(TN_OpSites % site1_Op1_Op2,  TN_OpSites % site1_Op1)

     TN_OpSites % site1_Op2 = TensKRON(ipeps(1) % m(k,:,:,:,:), ipepsHC(1) % m(l,:,:,:,:)) 
     TN_OpSites % site2_Op2 = TensKRON(ipeps(2) % m(k,:,:,:,:), ipepsHC(2) % m(l,:,:,:,:)) 
  end IF

 END SUBROUTINE create_reduced_psi_sites_from_ipeps





 !!! [Reduced-rho iPEPO network] --> create untraced |PSI> [<Ancilla|Ancilla>] <PSI| sites !!!
 SUBROUTINE create_purified_psi_sites_from_ipepo(TN_OpSites, ipepo, i, j, k, l)

  TYPE(op_sites_type),             INTENT(INOUT)        :: TN_OpSites
  TYPE(block_pepo),                INTENT(IN)           :: ipepo(:)
  INTEGER,                         INTENT(IN)           :: i, j
  INTEGER,                         INTENT(IN), OPTIONAL :: k, l

  !! Local HC copy
  TYPE(block_pepo), ALLOCATABLE :: ipepoP(:)

  !! Create a local HC copy of iPEPO
  CALL copy_pepo_block(ipepoP, ipepo)

  !! Compute reduced TN sites-1,2 -- bare tensors, no operator applied, just rho traced out 
  !! (dn-leg = real Hilbert space, up-leg = ancilla space)
  TN_OpSites % site1_Op1        = TENSMUL(ipepo(1) % m(i,:,:,:,:,:),  CONJG(ipepoP(1) % m(j,:,:,:,:,:)), MULT='11', FUSE='(22,33,44,55)')

  IF(PRESENT(k) .AND. PRESENT(l)) THEN

     CALL copyTens(TN_OpSites % site1_Op1_Op2,  TN_OpSites % site1_Op1)

     TN_OpSites % site1_Op2     = TENSMUL(ipepo(1) % m(k,:,:,:,:,:),  CONJG(ipepoP(1) % m(l,:,:,:,:,:)), MULT='11', FUSE='(22,33,44,55)') 
     TN_OpSites % site2_Op2     = TENSMUL(ipepo(2) % m(k,:,:,:,:,:),  CONJG(ipepoP(2) % m(l,:,:,:,:,:)), MULT='11', FUSE='(22,33,44,55)') 
  end IF

 END SUBROUTINE create_purified_psi_sites_from_ipepo





 !!! Construct a single site of reduced TN -- <Psi|Psi> network !!!
 FUNCTION psi_tmat_site(ipeps, op, opP)

  TYPE(block_peps), INTENT(IN)           :: ipeps                  !ipeps site (nolambda)
  COMPLEX(KIND=DP), INTENT(IN), OPTIONAL :: op(:,:), opP(:,:)      !operators applied to the site (in case we apply two operators)
  COMPLEX(KIND=DP), ALLOCATABLE          :: psi_tmat_site(:,:,:,:) !result

  !! Intermediate tens produced by applying op to ipeps
  COMPLEX(KIND=DP), ALLOCATABLE :: tempSite(:,:,:,:,:) 

  IF(PRESENT(op) .AND. PRESENT(opP)) THEN
     !! Apply op and opP
     tempSite = TENSMUL(TENSMUL(ipeps % m,  op,  '12'),  opP,  '12')    

  ELSEIF(PRESENT(op)) THEN
     !! Apply op
     tempSite = TENSMUL(ipeps % m,  op,  '12') 
  ELSE
     CALL copyTens(tempSite, ipeps % m)
  end IF

  !! Take an inner product (psi x psi_HC)
  CALL copyTens(psi_tmat_site, TENSMUL(tempSite,  CONJG(ipeps % m),  MULT='11',  FUSE='(22,33,44,55)'))

 END FUNCTION psi_tmat_site

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SIMPLIFIED CONTRACTION ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute 1P observable (simplified version) !!!
 SUBROUTINE simplified_compute_obs_1P(expval, C_norm, ipeps, Cmat, Tmat, CTMRG, op)

  USE boundary_tensors

  COMPLEX(KIND=DP),        INTENT(INOUT) :: expval, C_norm 
  TYPE(block_peps),        INTENT(IN)    :: ipeps(:)
  TYPE(ctm_corner_type),   INTENT(IN)    :: Cmat(:)                
  TYPE(ctm_transfer_type), INTENT(IN)    :: Tmat(:)
  TYPE(ctm_params),        INTENT(IN)    :: CTMRG                 
  COMPLEX(KIND=DP),        INTENT(IN)    :: op(:,:)

  !! Bulk 2D network (TN2D)  && Transfer-MPO network tensors (Rvec, Lvec, Tmpo)
  TYPE(block_mpo),    ALLOCATABLE :: TN2D(:)
  TYPE(transfer_mpo), ALLOCATABLE :: Tmpo(:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: Rvec(:,:,:), Lvec(:,:,:)

  !! Local nolambda iPEPS && operator site
  TYPE(op_sites_type) :: TN_OpSites

  !! (1) Create reduced TN && Construct bounded MPS-biMPO-MPS network from CTM tensors
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)
  CALL construct_Rvec_Tmpo_Lvec_network(Rvec, Tmpo, Lvec, Cmat, Tmat, TN2D, CTMRG)

  !! (2) Compute 1P observable
  CALL create_reduced_op_sites_from_ipeps(TN_OpSites, ipeps, op)
  CALL contract_TN_1P_obs(expval, Rvec, Tmpo, Lvec, TN_opSites)

  !! (3) Compute TN norm (observable already normalized inside contract_TN_1P_obs)
  CALL mult_bvec_tmpo_bvec(C_norm, Rvec, Tmpo, Lvec)
  
 END SUBROUTINE simplified_compute_obs_1P





 !!! Contract bare TN !!!
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
       DO s=1,2
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


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 !!! Resize Gamma-Lambda iPEPS to a new bond dim, both the unit cell and its environment !!!
 SUBROUTINE resize_ipeps_GamLam(iGamma, iLambda, Cmat, Tmat, chi)

  TYPE(block_peps),        INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda),      INTENT(INOUT) :: iLambda(:)
  TYPE(ctm_corner_type),   INTENT(INOUT) :: Cmat(:)     !CTM Cmat
  TYPE(ctm_transfer_type), INTENT(INOUT) :: Tmat(:)     !CTM Tmat
  INTEGER,                 INTENT(IN)    :: chi         !Target chi

  !! Bulk 2D network (TN2D)
  TYPE(block_mpo),  ALLOCATABLE :: TN2D(:)
  TYPE(block_peps), ALLOCATABLE :: ipeps(:)

  !! Resize iPEPS
  CALL resize_peps_bonds(iGamma, iLambda, chi)

  !! Resize CTM
  CALL absorb_ipeps_lambdas(ipeps, iGamma, iLambda)
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)
  CALL resize_transfer_mat(Tmat, TN2D)

 END SUBROUTINE resize_ipeps_GamLam




 !!! Resize plain iPEPS to a new bond dim, both the unit cell and its environment !!!
 SUBROUTINE resize_ipeps_PLAIN(ipeps, Cmat, Tmat, chi, noise_str)

  TYPE(block_peps),        INTENT(INOUT)         :: ipeps(:)    !iPEPS
  TYPE(ctm_corner_type),   INTENT(INOUT)         :: Cmat(:)     !CTM Cmat
  TYPE(ctm_transfer_type), INTENT(INOUT)         :: Tmat(:)     !CTM Tmat
  INTEGER,                 INTENT(IN)            :: chi         !Target chi
  CHARACTER(LEN=*),        INTENT(IN),  OPTIONAL :: noise_str   !Incl noise str to add some noise

  !! Bulk 2D network (TN2D)
  TYPE(block_mpo),    ALLOCATABLE :: TN2D(:)

  !! Noise tensor
  COMPLEX(KIND=DP),   ALLOCATABLE :: eps_tens(:,:,:,:,:)

  !! Resize iPEPS
  CALL resize_peps_bonds(ipeps, chi)

  !! Resize CTM
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)
  CALL resize_transfer_mat(Tmat, TN2D)

  !! Optionally add noise
  IF (PRESENT(noise_str)) THEN
     ALLOCATE(eps_tens(ipeps(1) % LocalDim, chi, chi, chi, chi));  eps_tens(:,:,:,:,:) = 1.0D-02 
     ipeps(1) % m = ipeps(1) % m + eps_tens;  ipeps(1) % m = add_noise(ipeps(1) % m)
     ipeps(2) % m = ipeps(2) % m + eps_tens;  ipeps(2) % m = add_noise(ipeps(2) % m)
  end IF

 END SUBROUTINE resize_ipeps_PLAIN


END MODULE imaginary_observables_ipeps
