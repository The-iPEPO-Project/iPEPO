MODULE CTM_semibound_routines

 USE utility
 USE datafile_utility, ONLY: setup_file_obs_1P, write_summary_data_1P

 USE definitions_mps_mpo
 USE mps_mpo_utility
 USE boundary_tensors

 USE ctm_definitions
 USE corner_transfer_matrix

 USE classical_ising_2D

 IMPLICIT NONE

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CTMRG routines for calculating semi-boundaries of TN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Main subroutine of CTMRG -- run CTMRG calculation to compute boundaries of an infinite 2D network !!!
 SUBROUTINE compute_CTMRG_semiboundaries(Cmat, Tmat, mps, TN2D, mpsP, CTMRG, use_old_ctm)

  USE datafile_utility, ONLY: setupConvDatafile

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2D(:)
  TYPE(block_mps),                      INTENT(IN)    :: mps(:), mpsP(:)
  TYPE(ctm_params),                     INTENT(IN)    :: CTMRG
  LOGICAL,                              INTENT(IN)    :: use_old_ctm

  !! CTM iteration && convergence params
  INTEGER            :: iter
  INTEGER, PARAMETER :: N_iters = 10000
  INTEGER            :: datafile
  LOGICAL            :: is_converged

  !! Singular vals of CTM
  COMPLEX(KIND=DP), ALLOCATABLE :: svals1(:), svals2(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: svals3(:), svals4(:)

  !! Initialize CTM (factorized network)
  CALL initialize_ctm_semibound(Cmat, Tmat, mps, TN2D, mpsP, chi=array(2,2,2,2), use_old_ctm=use_old_ctm)

  !! Setup file for tracking convergence
  CALL setupConvDatafile(datafile=datafile, descr="CTMRG_semibound", valstr="Max_Sval")

  !! Initialize convergence variable
  is_converged = .FALSE.

  mainloop: DO iter=1,N_iters

     WRITE(*,*)
     WRITE(*,*) "CTMRG, iter = ", iter
     WRITE(*,*)

     !! Single iteration of CTMRG
     CALL perform_CTMRG_semibound_one_iter(Cmat, Tmat, TN2D, CTMRG)

     !! Check convergence
     CALL check_CTM_convergence(is_converged, svals1, svals2, svals3, svals4, Cmat, CTMRG, datafile, iter)

     IF(is_converged) THEN  
        EXIT mainloop
     ELSEIF(iter .EQ. N_iters) THEN
        WRITE(*,*) "CTMRG HAS FAILED TO CONVERGE"
        STOP 
     end IF     
  end DO mainloop

  !! Close files
  CLOSE(datafile)
  WRITE(*,*) " CTMRG DONE "

 END SUBROUTINE compute_CTMRG_semiboundaries




 !!! Perform a single iteration of CTMRG semi-boundaries calculation !!!
 SUBROUTINE perform_CTMRG_semibound_one_iter(Cmat, Tmat, TN2Din, CTMRG)

  TYPE(ctm_corner_type),   INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),         INTENT(IN)    :: TN2Din(:)
  TYPE(ctm_params),        INTENT(IN)    :: CTMRG

  !! Local copy of 2D network (so we can swap sites after each CTM move)
  TYPE(block_mpo), ALLOCATABLE :: TN2D(:)

  !! Create a local copy of TN2D
  CALL copy_mpo_block(TN2D, TN2Din)

  !!! Evolve in W-dir (LHS bound) !!!
  CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR='W')
  CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR='W')

  !!! Evolve in E-dir (RHS bound) !!!
  CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR='E')
  CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR='E')

 END SUBROUTINE perform_CTMRG_semibound_one_iter

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INITIALIZE semi-boundary CTM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Initialize semi-boundary CTM !!!
 SUBROUTINE initialize_ctm_semibound(Cmat, Tmat, mps, TN2D, mpsP, chi, use_old_ctm)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT)  :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT)  :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)     :: TN2D(:)
  TYPE(block_mps),                      INTENT(IN)     :: mps(:), mpsP(:)
  INTEGER,                              INTENT(IN)     :: chi(:)
  LOGICAL,                              INTENT(IN)     :: use_old_ctm

  !! Initialize semi-boundary CTM, either by reusing an old CTM, or by initializing a new random CTM
  IF(use_old_ctm) THEN
     CALL resize_transfer_mat(Tmat, TN2D) 
  ELSE
     CALL initialize_random_ctm_semibound(Cmat, Tmat, mps, TN2D, mpsP, chi)
  end IF

 END SUBROUTINE initialize_ctm_semibound




 !!! Initialize random semi-boundary CTM !!!
 SUBROUTINE initialize_random_ctm_semibound(Cmat, Tmat, mps, TN2D, mpsP, chi)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2D(:)
  TYPE(block_mps),                      INTENT(IN)    :: mps(:), mpsP(:)
  INTEGER,                              INTENT(IN)    :: chi(:)

  INTEGER :: N_sites, N_sub
  INTEGER :: site, sub

  !Allocate new empty CTM
  CALL allocate_CTM_semibound(Cmat, Tmat, mps, TN2D, mpsP, chi)

  !Get allocated sizes of Cmat, Tmat blocks
  N_sites = SIZE(Tmat)
  N_sub   = SIZE(Tmat(1) % T)

  !Initialize random C-tensors
  DO site=1,N_sites
     CALL initialize_rand_tens(Cmat(site) % m)
  end DO

  !Initialize random T-tensors
  DO sub=1,N_sub
     CALL initialize_rand_tens(Tmat(2) % T(sub) % m)
     CALL initialize_rand_tens(Tmat(4) % T(sub) % m)
  end DO

 END SUBROUTINE initialize_random_ctm_semibound

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE semi-boundary CTM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Allocate semi-boundary CTM !!!
 SUBROUTINE allocate_CTM_semibound(Cmat, Tmat, mps, TN2D, mpsP, chi)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT)  :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT)  :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)     :: TN2D(:)
  TYPE(block_mps),                      INTENT(IN)     :: mps(:), mpsP(:)
  INTEGER,                              INTENT(IN)     :: chi(:)

  CALL allocate_ctm_semibound_corner_mat(Cmat,   mps,       mpsP, chi)
  CALL allocate_ctm_semibound_transfer_mat(Tmat, mps, TN2D, mpsP, chi)

 END SUBROUTINE allocate_CTM_semibound




 !!! Allocate new semi-boundary CTM corner matrix !!!
 SUBROUTINE allocate_ctm_semibound_corner_mat(Cmat, mps, mpsP, chi)

  TYPE(ctm_corner_type), ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(block_mps),                    INTENT(IN)    :: mps(:), mpsP(:)
  INTEGER,                            INTENT(IN)    :: chi(:)

  INTEGER, PARAMETER :: N_sites = 4
  INTEGER            :: site

  INTEGER, ALLOCATABLE :: dimC(:,:)

  !! Allocate empty corner matrix block
  CALL allocate_empty_ctm_corner_mat(Cmat, N_sites)

  !! Set CTM site dims
  CALL get_ctm_semibound_bondDims(dimC, mps, mpsP, chi, N_sites)

  !! Allocate individual corner nodes 
  DO site=1,N_sites
     CALL allocate_ctm_corner_site(Cmat(site), dimC(site,:))
  end DO

 END SUBROUTINE allocate_ctm_semibound_corner_mat




 !!! Allocate new semi-boundary CTM transfer matrix !!!
 SUBROUTINE allocate_ctm_semibound_transfer_mat(Tmat, mps, TN2D, mpsP, chi)

  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2D(:)
  TYPE(block_mps),                      INTENT(IN)    :: mps(:), mpsP(:)
  INTEGER,                              INTENT(IN)    :: chi(:)

  INTEGER, PARAMETER :: N_sites = 4
  INTEGER, PARAMETER :: N_sub   = 2
  INTEGER            :: site, sub

  INTEGER, ALLOCATABLE :: dimT(:,:,:)
  INTEGER, ALLOCATABLE :: localDims(:,:), bondDims(:,:)

  !! Check the number of subsites is correct
  CALL check_sizes_equal(N_sub, SIZE(TN2D), "allocate_ctm_transfer_mat: N_sub must equal size of TN2D")

  !! (1) Set CTM site dims
  CALL get_ctm_localDims(localDims, TN2D, N_sites)
  CALL get_ctm_semibound_bondDims(bondDims, mps, mpsP, chi, N_sites)

  ALLOCATE(dimT(N_sites, N_sub, 3))
  dimT(:,:,1)   = localDims
  dimT(:,:,2)   = bondDims 
  dimT(:,:,3)   = bondDims

  !! (2) Allocate empty transfer matrix block
  CALL allocate_empty_ctm_transfer_mat(Tmat, N_sites)

  !! (3) Allocate individual transfer sites -- LATERAL sites remain to be determined 
  CALL allocate_ctm_transfer_site(Tmat(2),  dimT(2, :, :))
  CALL allocate_ctm_transfer_site(Tmat(4),  dimT(4, :, :))

  !! (4) Allocate individual transfer sites -- UP/DOWN sites will contain mps && mpsP
  CALL allocate_ctm_transfer_site(Tmat(3),  mps,  flagUpDown='UP')
  CALL allocate_ctm_transfer_site(Tmat(1),  mpsP, flagUpDown='DN')

 END SUBROUTINE allocate_ctm_semibound_transfer_mat


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! GETTING CTM DIMENSIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Get CTM bond dims !!!
 SUBROUTINE get_ctm_semibound_bondDims(bondDims, mps, mpsP, chi, N_sites)

  INTEGER, ALLOCATABLE, INTENT(INOUT) :: bondDims(:,:)
  TYPE(block_mps),      INTENT(IN)    :: mps(:), mpsP(:)
  INTEGER,              INTENT(IN)    :: chi(:)
  INTEGER,              INTENT(IN)    :: N_sites

  INTEGER :: site

  !! Check chi(:) input has the correct size
  CALL check_sizes_equal(SIZE(chi), N_sites, "get_ctm_semibound_bondDims: size of chi must equal N_sites")

  !! Allocate bondDims
  ALLOCATE(bondDims(N_sites, 2))

  !! Read input chi(:) && MPS outer bond dims
  bondDims(1, :) = (/ mpsP(1) % Wdim,   chi(1)         /)
  bondDims(2, :) = (/ chi(2),           mpsP(2) % Edim /)
  bondDims(3, :) = (/ mps(2)  % Edim,   chi(3)         /)
  bondDims(4, :) = (/ chi(4),           mps(1)  % Wdim /)

 END SUBROUTINE get_ctm_semibound_bondDims

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE CTM_semibound_routines






PROGRAM test_mult_classical_ising_CTMRG

 USE utility
 USE datafile_utility, ONLY: setup_file_obs_1P, write_summary_data_1P

 USE definitions_mps_mpo
 USE mps_mpo_algebra_inf

 USE ctm_definitions
 USE corner_transfer_matrix
 USE TN_contractions

 USE CTM_semibound_routines

IMPLICIT NONE

 !! CTM objects
 TYPE(ctm_corner_type),   ALLOCATABLE :: Cmat(:)
 TYPE(ctm_transfer_type), ALLOCATABLE :: Tmat(:)

 !! RHS/LHS boundaries && Transfer MPO of infinite 2D network
 COMPLEX(KIND=DP),   ALLOCATABLE :: Rvec(:,:,:), Lvec(:,:,:)
 TYPE(transfer_mpo), ALLOCATABLE :: Tmpo(:), TmpoTransp(:)

 !! 2D TN of classical Ising partition function
 TYPE(block_mpo), ALLOCATABLE :: TN2D(:)

 !! 2D Ising Tmat && its transpose
 TYPE(block_mpo),    ALLOCATABLE :: Z_mpo(:,:), Z_mpoTransp(:,:)

 !! Right && Left Bounds
 TYPE(block_mps),    ALLOCATABLE :: Rmps(:),  Lmps(:)
 TYPE(block_mps),    ALLOCATABLE :: R_Gam(:), L_Gam(:)
 TYPE(block_lambda), ALLOCATABLE :: R_Lam(:), L_Lam(:)

 !! Copies of Right && Left Bounds for restoring at a later point
 TYPE(block_mps),    ALLOCATABLE :: R_Gam_Cpy(:), L_Gam_Cpy(:)
 TYPE(block_lambda), ALLOCATABLE :: R_Lam_Cpy(:), L_Lam_Cpy(:)

 !! Partition function Zval && temperature beta
 COMPLEX(KIND=DP)         :: Zval, ZvalTransp, ZvalBrute
 REAL(KIND=DP), PARAMETER :: beta = 0.4407d0

 !! File indexing
 INTEGER :: datafile

 !! SVD && CTM params
 TYPE(ctm_params)         :: CTMRG
 INTEGER,       PARAMETER :: chi = 1  
 REAL(KIND=DP), PARAMETER :: eps = 1.0d-15 
 REAL(KIND=DP), PARAMETER :: eta = 1.0d-15 


 !! Create a new instance of CTM
 CALL create_CTM_params(CTMRG, eps, chi, eta, ctm_method='1')


 !!!!!!!!!!!!!!!!!!!!!!!!!! Brute-force contraction of Z-network !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !! Construct Ising Tmat && Boundaries 
 CALL construct_ising_tmat(Z_mpo, beta)
 CALL construct_ising_boundVec(R_Gam, R_Lam, beta); CALL normalize_imps(R_Gam, R_Lam)
 CALL construct_ising_boundVec(L_Gam, L_Lam, beta); CALL normalize_imps(L_Gam, L_Lam)

 !! Create copies of boundary vecs (for transpose contraction)
 CALL copy_imps(R_Gam_Cpy, R_Lam_Cpy, R_Gam, R_Lam)
 CALL copy_imps(L_Gam_Cpy, L_Lam_Cpy, L_Gam, L_Lam)

 !! Brute-force calculation of Zval
 CALL calc_Zval_bruteforce(ZvalBrute, R_Gam, R_Lam, L_Gam, L_Lam, Z_mpo)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CTMRG contraction of Z-network !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !! (1) Create 2D network to represent classical Ising partition function
 CALL create_reduced_TN_from_classical_Z(TN2D, beta)
 CALL absorb_imps_lambdas(Rmps, R_Gam, R_Lam)
 CALL absorb_imps_lambdas(Lmps, L_Gam, L_Lam)

 !! (2) Compute CTMRG boundaries of TN
 CALL compute_CTMRG_semiboundaries(Cmat, Tmat, Rmps, TN2D, Lmps, CTMRG, use_old_ctm=.FALSE.)

 !! (3) Reconstruct TN as CTM-bounded MPS-MPO-MPS network
 CALL construct_Rvec_Tmpo_Lvec_network(Rvec, Tmpo, Lvec, Cmat, Tmat, TN2D, CTMRG)

 !! (4) Find the norm of TN
 CALL mult_bvec_tmpo_bvec(Zval, Rvec, Tmpo, Lvec)

 !! (5) Repeat the calculation with transposed transfer MPO
 !!CALL transpose_mpo_block(TmpoTransp, Tmpo)
 !!CALL mult_bvec_tmpo_bvec(ZvalTransp, Lvec, TmpoTransp, Rvec)

 !! (6) Print results
 WRITE(*,*) "Zval =       ", Zval
 !!WRITE(*,*) "ZvalTransp = ", ZvalTransp
 WRITE(*,*) "ZvalBrute =  ", ZvalBrute


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!! CONSTRUCT ISING BOUNDARY VECS -- FOR 2D LATTICE INFINITE ALONG X, FINITE ALONG Y !!!!!!!!!!!!!!!!!!!!!!!!!!!


 !Construct Ising transfer matrix
 SUBROUTINE construct_ising_tmat(Tmat, beta, N_sites_In, opsite_1, opsite_2)
 
  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT)        :: Tmat(:,:)
  REAL(KIND=DP),                INTENT(IN)           :: beta
  INTEGER,                      INTENT(IN), OPTIONAL :: N_sites_In
  INTEGER,                      INTENT(IN), OPTIONAL :: opsite_1, opsite_2

  TYPE(block_mpo), ALLOCATABLE :: TmatA(:), TmatB(:)
  INTEGER                      :: N_sites, site

  !Determine the number of sites
  N_sites = OptArg(N_sites_In, 2)

  !Allocate empty Tmat
  CALL allocate_empty_mpo_block(TmatA, N_sites)

  !Construct the sites of Tmat
  DO site=1,N_sites
     CALL allocate_mpo_site(TmatA(site), classical_ising_site(beta))
  end DO

  !Create a copy
  CALL copy_mpo_block(TmatB, TmatA)

  !Insert operators at given sites opsite_1, opsite_2 (if provided -- otherwise just a bare Tmat)
  IF(PRESENT(opsite_1)) CALL allocate_mpo_site(TmatB(opsite_1), classical_ising_site(beta, 'OP'))
  IF(PRESENT(opsite_2)) CALL allocate_mpo_site(TmatB(opsite_2), classical_ising_site(beta, 'OP'))

  !Combine TmatA, TmatB
  CALL combine_two_mpos_into_bimpo(Tmat, TmatA, TmatB)

 END SUBROUTINE construct_ising_tmat


 
 !Construct Ising boundary vec
 SUBROUTINE construct_ising_boundVec(Gamma_bound, Lambda_bound, beta, N_sites_In)
 
  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT)        :: Gamma_bound(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)        :: Lambda_bound(:)
  REAL(KIND=DP),                   INTENT(IN)           :: beta
  INTEGER,                         INTENT(IN), OPTIONAL :: N_sites_In

  INTEGER :: N_sites, site

  !Determine the number of sites
  N_sites = OptArg(N_sites_In, 2)

  !Allocate empty bound vec (gammas && lambdas)
  CALL allocate_empty_mps_block(Gamma_bound, N_sites)
  CALL allocate_empty_lambda_block(Lambda_bound, N_sites)

  !Construct the sites of boundary vec (gammas && lambdas)
  DO site=1,N_sites
     CALL allocate_mps_site(Gamma_bound(site), boundVec_site(beta))
     CALL allocate_lambda_site(Lambda_bound(site), DiagonalAsVec(matEye(2)))
  end DO

 END SUBROUTINE construct_ising_boundVec



 !Construct a single site of Ising boundary vec
 FUNCTION boundVec_site(beta)

  REAL(KIND=DP), INTENT(IN) :: beta
  COMPLEX(KIND=DP)          :: boundVec_site(2,2,2) !result

  INTEGER :: i,j,k !indices

  !Compute Tmat site element by element (i,j,k,l = spin-state indices indicating spin-states on adjacent sites)
  DO i=1,2
    DO j=1,2
      DO k=1,2
         boundVec_site(i,j,k) = boundVec_site_element(beta,i,j,k) 
      end DO
    end DO
  end DO

 END FUNCTION boundVec_site



 !!! Compute a single matrix element of a given Tmat site (indices correspond to spin states at adjacent sites) !!!
 FUNCTION boundVec_site_element(beta, i, j, k)

  REAL(KIND=DP), INTENT(IN) :: beta                  !temperature
  INTEGER,       INTENT(IN) :: i, j, k               !spin indices
  COMPLEX(KIND=DP)          :: boundVec_site_element !result = element at indices i,j,k,l

  COMPLEX(KIND=DP) :: Q_sqrt(2,2) !Q_sqrt matrix
  INTEGER          :: s           !local spin on our site

  !setup Q_sqrt matrix
  CALL setup_Q_sqrt(Q_sqrt, beta)

  !initialize
  boundVec_site_element = (0.0d0, 0.0d0)
     
  !Compute boundVec site element
  DO s=1,2
     boundVec_site_element = boundVec_site_element + Q_sqrt(i,s) * Q_sqrt(j,s) * Q_sqrt(k,s)
  end DO

 END FUNCTION boundVec_site_element



 !!! CONSTRUCT ISING BOUNDARY VECS -- FOR 2D LATTICE INFINITE ALONG X, FINITE ALONG Y !!!
 SUBROUTINE calc_Zval_bruteforce(Zval, R_Gamma, R_Lambda, L_Gamma, L_Lambda, Tmat)

  USE mps_mpo_utility

  COMPLEX(KIND=DP),   INTENT(OUT) :: Zval
  TYPE(block_mps),    INTENT(IN)  :: R_Gamma(:), L_Gamma(:)
  TYPE(block_lambda), INTENT(IN)  :: R_Lambda(:), L_Lambda(:)
  TYPE(block_mpo),    INTENT(IN)  :: Tmat(:,:)

  TYPE(block_mps), ALLOCATABLE :: R_mps(:), L_mps(:)

  !Local objects for deriving TransMat
  COMPLEX(KIND=DP), ALLOCATABLE :: TransMat(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: TransMat1(:,:), TransMat2(:,:) 
  COMPLEX(KIND=DP), ALLOCATABLE :: mpsTmp1(:,:,:), mpsTmp2(:,:,:)

  CALL absorb_imps_lambdas(R_mps, R_Gamma, R_Lambda, '12')
  CALL absorb_imps_lambdas(L_mps, L_Gamma, L_Lambda, '12')

  !Compute TransMat = Rank-2 version of Tmat
  mpsTmp1 = TENSMUL(Tmat(1,1) % m, TENSMUL(Tmat(2,1) % m, R_mps(1) % m, MULT='21', FUSE='(32,43)'), MULT='21', FUSE='(32,43)')
  !mpsTmp1 = TENSMUL(Tmat(2,1) % m, R_mps(1) % m, MULT='21', FUSE='(32,43)')
  TransMat1 = TENSMUL(mpsTmp1, L_mps(1) % m, MULT='11', FUSE='(22,33)')

  mpsTmp2 = TENSMUL(Tmat(1,2) % m, TENSMUL(Tmat(2,2) % m, R_mps(2) % m, MULT='21', FUSE='(32,43)'), MULT='21', FUSE='(32,43)')
  !mpsTmp2 = TENSMUL(Tmat(2,2) % m, R_mps(2) % m, MULT='21', FUSE='(32,43)')
  TransMat2 = TENSMUL(mpsTmp2, L_mps(2) % m, MULT='11', FUSE='(22,33)')

  TransMat = TENSMUL(TransMat1, TransMat2)

  !Get dominant eval = Zval
  CALL dominant_eig(eval=Zval, tensA=TransMat)
  !WRITE(*,*) "Zval brute force = ", Zval
   
 END SUBROUTINE calc_Zval_bruteforce

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM test_mult_classical_ising_CTMRG
