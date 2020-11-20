PROGRAM test_classical_ising_CTMRG

 USE utility
 USE datafile_utility, ONLY: setup_file_obs_1P, write_summary_data_1P, &
                           & setup_file_obs_2P, write_summary_data_2P

 USE definitions_mps_mpo
 USE mps_mpo_utility
 USE boundary_tensors
 USE mps_mpo_algebra_inf

 USE ctm_definitions
 USE corner_transfer_matrix

 USE TN_contractions
 USE classical_ising_2D

IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Physical param - loop over inverse temperature beta
  REAL(KIND=DP)            :: beta          
  REAL(KIND=DP), PARAMETER :: beta0      =  0.45d0    !Initial value of beta
  REAL(KIND=DP), PARAMETER :: delta_beta = -0.001d0   !Increment of beta in the loop
  INTEGER,       PARAMETER :: N_beta     =  25        !Number of beta values  
  INTEGER                  :: i                       !Loop index

  !New params && indices after readjusting the size of beta step
  REAL(KIND=DP), PARAMETER :: beta1       =  0.441d0
  REAL(KIND=DP), PARAMETER :: delta_beta1 = -0.0001d0
  INTEGER                  :: i1

  !! SVD && CTM params
  TYPE(ctm_params)  :: CTMRG
  REAL(KIND=DP)     :: eps 
  INTEGER           :: chi 
  REAL(KIND=DP)     :: eta
  CHARACTER(LEN=1)  :: ctm_method

  !! Separation vars for 2P corrs
  INTEGER :: sep
  INTEGER :: max_sep 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Spin-spin corrs && Magnetization expval (TN contraction && exact result)
  COMPLEX(KIND=DP), ALLOCATABLE :: expvals_2P(:), SiSj_exp(:)
  COMPLEX(KIND=DP)              :: expval, m_exact
  COMPLEX(KIND=DP)              :: C_norm, vec_norm

  !! 2D TN of classical Ising partition function && spin operator sites
  TYPE(block_mpo), ALLOCATABLE :: TN2D(:)
  TYPE(op_sites_type)          :: TN_Spin_Sites

  !! CTM objects
  TYPE(ctm_corner_type),   ALLOCATABLE :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE :: Tmat(:)

  !! RHS/LHS boundaries && Transfer MPO of infinite 2D network
  COMPLEX(KIND=DP),   ALLOCATABLE :: Rvec(:,:,:), Lvec(:,:,:)
  TYPE(transfer_mpo), ALLOCATABLE :: Tmpo(:)

  !! File indexing
  INTEGER :: datafile, corrs_datafile

  WRITE(*,*) "Reading CTMRG precision settings"
  READ(*,*) chi, eps, eta

  WRITE(*,*) "Reading CTMRG method selection"
  READ(*,*) ctm_method

  WRITE(*,*) "Reading max separation"
  READ(*,*) max_sep

  !! Create a new instance of CTM
  CALL create_CTM_params(CTMRG, eps, chi, eta, ctm_method)

  !! Determine the index at which we'll change the beta step size
  i1 = (beta1 - beta0)/delta_beta + 1

  !! Prepare files && array for writing data
  ALLOCATE(expvals_2P(max_sep + 1))
  CALL setup_file_obs_1P(datafile, "magnetization_summary", "beta", "m",  "exact sol m", "eval")
  CALL setup_file_obs_2P(corrs_datafile, descriptor="SiSj_corrs_summary",  param_str="beta", sepX_str="sep", &
                         & obs_str="SiSj", extra_str_1="C_norm", extra_str_2="exp SiSj")

  !! Loop over vals of beta
  DO i=1,N_beta 

     !! (0) Set temperature beta (Readjust the step size for small betas)
     IF(i .GT. i1) THEN
       beta = beta0 + (i1-1)*delta_beta + (i-i1)*delta_beta1
     ELSE
       beta = beta0 + (i-1)*delta_beta
     end IF

     !! (1) Create 2D network to represent classical Ising partition function
     CALL create_reduced_TN_from_classical_Z(TN2D, beta)
     CALL create_reduced_op_sites_from_classical_Z(TN_Spin_Sites, beta)

     !! (2) Calculate environment of infinite 2D tensor network
     CALL compute_CTMRG_boundaries(Cmat, Tmat, TN2D, CTMRG, use_old_ctm = (i .GT. 1))
     CALL construct_Rvec_Tmpo_Lvec_network(Rvec, Tmpo, Lvec, Cmat, Tmat, TN2D, CTMRG)

     !! (2.1) Calculate new norm of TN
     CALL mult_bvec_bvec(vec_norm,      Rvec,       Lvec)
     CALL mult_bvec_tmpo_bvec(C_norm,   Rvec, Tmpo, Lvec)

     !! (3)   Compute 1P magnetization
     CALL contract_TN_1P_obs(expval, Rvec, Tmpo, Lvec, TN_Spin_Sites)

     !! FIXME -- bruteforce calculation of expvals
     !CALL bruteforce_contraction(C_norm, Cmat, Tmat, TN2D)
     !CALL bruteforce_contraction(expval, Cmat, Tmat, TN2D, TN_Spin_Sites)
     !expval = expval/C_norm

     !! (3.1) Compute the exact solution -- 1P magnetization
     m_exact = (1 - sinh(2*beta)**(-4.0d0))**(1.0d0/8.0d0)

     !! (3.2) Write output to file -- 1P expvals
     CALL write_summary_data_1P(datafile, beta, expval, m_exact, C_norm)

     !! (4) Compute 2P spin-spin corrs
     CALL contract_TN_2P_corr(expvals_2P, Rvec, Tmpo, Lvec, TN_Spin_Sites, min_sep=0, max_sep=max_sep)

     !! (4.1) Write output to file -- 2P corrs
     CALL calc_exp_spin_corrs(SiSj_exp, expvals_2P, m_exact)
     DO sep=0,max_sep
        CALL write_summary_data_2P(corrs_datafile, beta, sep, expvals_2P(sep+1), C_norm, SiSj_exp(sep+1))
     end DO
     WRITE(corrs_datafile, *)
     WRITE(corrs_datafile, *)

  END DO

  !! Close files
  CLOSE(datafile)
  CLOSE(corrs_datafile)

CONTAINS

  !!! Get the exact solution for exponentially decaying correlations in classical Ising model !!!
  SUBROUTINE calc_exp_spin_corrs(sisj_exp, expvals_2P, m_exact)

   COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: SiSj_exp(:)
   COMPLEX(KIND=DP),              INTENT(IN)    :: expvals_2P(:)
   COMPLEX(KIND=DP),              INTENT(IN)    :: m_exact

   !! Parametric coefficients 
   COMPLEX(KIND=DP) :: A_2p, B_2p, xi

   !! Separations
   INTEGER :: max_sep, sep

   !! Get max sep
   max_sep = SIZE(expvals_2P) - 1

   !! Get parametric coefficients for exp decay function
   A_2p = ABS(expvals_2P(1) - m_exact**2) 
   B_2p = ABS(expvals_2P(2) - m_exact**2)
   xi = (log(A_2p) - log(B_2p))**(-1)

   !! Compute 2P corrs with exponential decay
   CALL allocateTens(SiSj_exp,  (/ max_sep+1 /))
   DO sep=0,max_sep
      SiSj_exp(sep+1) = A_2p * exp(-sep*1.0d0/xi) + m_exact**2
   end DO

  END SUBROUTINE calc_exp_spin_corrs





  !!! ROUTINES FOR TESTING/DEBUGGING NEW CALCULATION OF OBSERVABLES !!!

 !!! Construct boundary MPS using CTM tensors !!!
 SUBROUTINE construct_boundary_MPS(mps, Cmat, Tmat, boundFlag)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps(:)
  TYPE(ctm_corner_type),        INTENT(IN)    :: Cmat(:)
  TYPE(ctm_transfer_type),      INTENT(IN)    :: Tmat(:)
  CHARACTER(LEN=*),             INTENT(IN)    :: boundFlag

  !! TN length
  INTEGER, PARAMETER :: N_sites = 4

  !! Allocate empty boundary MPS
  CALL allocate_empty_mps_block(mps, N_sites)

  SELECT CASE(boundFlag)
  CASE('R')

       !! RHS (NORTH) boundary MPS
       CALL allocate_mps_site(mps(1),        increase_rank(              Cmat(4) % m,  '2'))
       CALL allocate_mps_site(mps(4),        increase_rank(TensTRANSPOSE(Cmat(3) % m), '3')) 

       CALL allocate_mps_site(mps(2),        Tmat(3) % T(2) % m) 
       CALL allocate_mps_site(mps(3),        Tmat(3) % T(1) % m)

  CASE('L')

       !! LHS (SOUTH) boundary MPS
       CALL allocate_mps_site(mps(1),        increase_rank(TensTRANSPOSE(Cmat(1) % m), '2'))
       CALL allocate_mps_site(mps(4),        increase_rank(              Cmat(2) % m,  '3'))

       CALL allocate_mps_site(mps(2),        TensTRANSPOSE(Tmat(1) % T(1) % m, '23'))  
       CALL allocate_mps_site(mps(3),        TensTRANSPOSE(Tmat(1) % T(2) % m, '23')) 

  CASE DEFAULT
       CALL invalid_flag("construct_boundary_MPS -- invalid boundFlag: ", boundFlag)
  end SELECT

 END SUBROUTINE construct_boundary_MPS




 !!! Construct transfer MPO using CTM tensors !!!
 SUBROUTINE construct_transfer_MPO(mpo, Tmat, TN2D)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: mpo(:,:)
  TYPE(ctm_transfer_type),      INTENT(IN)    :: Tmat(:)
  TYPE(block_mpo),              INTENT(IN)    :: TN2D(:)

  !! TN length
  INTEGER, PARAMETER :: N_sites = 4
  
  !! Allocate empty binary MPO block
  CALL allocate_empty_bimpo_block(mpo, N_sites)


  !! Setup MPO(2,:) -- the RHS (NORTH) part of transfer MPO
  CALL allocate_mpo_site(mpo(2,1),     TensROTATE(Tmat(4) % T(2) % m,  'CCW -PI/2',  '2'))
  CALL allocate_mpo_site(mpo(2,4),     TensROTATE(Tmat(2) % T(1) % m,  'CW +PI/2',   '2'))

  CALL allocate_mpo_site(mpo(2,2),     TN2D(2) % m) 
  CALL allocate_mpo_site(mpo(2,3),     TN2D(1) % m)


  !! Setup MPO(1,:) -- the LHS (SOUTH) part of transfer MPO
  CALL allocate_mpo_site(mpo(1,1),     TensROTATE(Tmat(4) % T(1) % m,  'CCW -PI/2',  '2'))
  CALL allocate_mpo_site(mpo(1,4),     TensROTATE(Tmat(2) % T(2) % m,  'CW +PI/2',   '2'))

  CALL allocate_mpo_site(mpo(1,2),     TN2D(1) % m) 
  CALL allocate_mpo_site(mpo(1,3),     TN2D(2) % m)

 END SUBROUTINE construct_transfer_MPO




 !!! Brute-force calculation of expectation values !!!
 SUBROUTINE bruteforce_contraction(expval, Cmat, Tmat, TN2D, TN_opSites)

  COMPLEX(KIND=DP),             INTENT(OUT)            :: expval
  TYPE(ctm_corner_type),        INTENT(IN)             :: Cmat(:)
  TYPE(ctm_transfer_type),      INTENT(IN)             :: Tmat(:)
  TYPE(block_mpo),              INTENT(IN)             :: TN2D(:)
  TYPE(op_sites_type),          INTENT(IN), OPTIONAL   :: TN_opSites       

  TYPE(block_mps),    ALLOCATABLE :: Rmps(:), Lmps(:)
  TYPE(block_mpo),    ALLOCATABLE :: mpo(:,:)

  !! Construct boundaries
  CALL construct_boundary_MPS(Rmps, Cmat, Tmat, 'R')
  CALL construct_boundary_MPS(Lmps, Cmat, Tmat, 'L')
  CALL construct_transfer_MPO(mpo,  Tmat, TN2D)

  !! Insert operator if provided
  IF(PRESENT(TN_opSites)) CALL copyTens(mpo(1,2) % m,  TN_OpSites % site1_Op1)

  !! Contract to find expval
  CALL simple_mult_mps_mpo(Rmps, mpo(2,:))
  CALL simple_mult_mps_mpo(Rmps, mpo(1,:))
  CALL contract_mps_mps(expval, Rmps, Lmps)

 END SUBROUTINE bruteforce_contraction

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM test_classical_ising_CTMRG
