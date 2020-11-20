MODULE imaginary_full_update_ipeps

 USE utility
 USE definitions_mps_mpo
 USE definitions_peps_pepo
 USE mps_mpo_utility

 USE project_ipeps
 USE peps_pepo_algebra

 USE ctm_definitions
 USE corner_transfer_matrix
 USE TN_contractions
 USE imaginary_observables_ipeps, ONLY: create_reduced_TN_from_ipeps

 IMPLICIT NONE

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (1) iPEPS FULL UPDATE MAIN ROUTINE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Full-update iPEPS bond !!!
 SUBROUTINE full_update_ipeps_bond(ipeps, Cmat, Tmat, mpo_prop, chi, eps, CTMRG, DIR)

  USE imaginary_observables_ipeps, ONLY: calc_CTM_environment_ipeps

  TYPE(block_peps),        ALLOCATABLE,  INTENT(INOUT)  :: ipeps(:)        !iGamma of PEPS
  TYPE(ctm_corner_type),   ALLOCATABLE,  INTENT(INOUT)  :: Cmat(:)         !CTM Cmat
  TYPE(ctm_transfer_type), ALLOCATABLE,  INTENT(INOUT)  :: Tmat(:)         !CTM Tmat
  TYPE(block_mpo),                       INTENT(IN)     :: mpo_prop(:)     !MPO Propagator
  INTEGER,                               INTENT(IN)     :: chi             !SVD chi
  REAL(KIND=DP),                         INTENT(IN)     :: eps             !SVD eps
  TYPE(ctm_params),                      INTENT(IN)     :: CTMRG           !CTM params
  CHARACTER(LEN=*),                      INTENT(IN)     :: DIR             !Bond we're propagating
  
  !! Reduced iPEPS reduced tensors --> MPS sites && Residual tensors
  TYPE(block_mps),    ALLOCATABLE :: mps(:),   mpsOO(:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: Xa(:,:),  Yb(:,:)

  !! Norm matrix && reduced 2D network
  COMPLEX(KIND=DP),   ALLOCATABLE :: N_norm(:,:,:,:)
  TYPE(block_mpo),    ALLOCATABLE :: TN2D(:)
  COMPLEX(KIND=DP)                :: C_norm

  !! (0) Calc iPEPS environment
  CALL calc_CTM_environment_ipeps(Cmat, Tmat, ipeps, CTMRG, use_old_ctm = .TRUE.)

  !! (1) Calculate reduced (factorized) iPEPS
  CALL calc_factorized_ipeps(mps, Xa, Yb, ipeps, DIR)

  !! (2A) Calculate iPEPS norm matrix
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)
  CALL calc_norm_matrix(N_norm, Xa, Yb, TN2D, Cmat, Tmat, CTMRG, DIR)

  !! (2B) Normalize N_norm
  CALL simple_contract_norm_network(C_norm, N_norm, mps)
  N_norm = N_norm / C_norm
  N_norm = 0.5D0*(N_norm + TensTRANSPOSE(TensTRANSPOSE(CONJG(N_norm), '13'), '24'))

  !! (3) Apply propagator MPO to reduced PEPS sites 
  CALL simple_mult_mps_mpo(mps, mpo_prop) 

  !! (4) Make a copy of mps = mps * mpo_prop 
  !!  -- mps   = approximate sites we want to optimize
  !!  -- mpsOO = exact sites we keep fixed and use as ref pt for optimization
  !!  -- we want to minimize the distance d = || |mps> - |mpsOO> ||
  CALL copy_mps_block(mpsOO, mps)

  !! (5) Find initial guess for CG optimization = [result from simple update]  
  CALL initialize_ALS_ipeps(mps,   chi,    eps)

  !! (6) Run CG optimization
  CALL ALS_optimize_ipeps(mps, mpsOO, N_norm)

  !! Normalize optimized tensors
  mps(1) % m = mps(1) % m / matrix_frobenius_norm(TensTRANSPOSE(RESHAPE_2D(mps(1) % m, '1,23')))
  mps(2) % m = mps(2) % m / matrix_frobenius_norm(TensTRANSPOSE(RESHAPE_2D(mps(2) % m, '1,23')))

  !! (7) Recover full iPEPS from the reduced (factorized) one
  CALL restore_unfactorized_ipeps(ipeps, mps, Xa, Yb, DIR)

  !! Rescale iPEPS by their max elements
  ipeps(1) % m = ipeps(1) % m / MAXVAL(ABS(ipeps(1) % m))
  ipeps(2) % m = ipeps(2) % m / MAXVAL(ABS(ipeps(2) % m))

 END SUBROUTINE full_update_ipeps_bond

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (2) NORM MATRIX CALCULATIONS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Calculate norm matrix in PSI network !!!
 SUBROUTINE calc_norm_matrix(N_norm, Xa, Yb, TN2D, Cmat, Tmat, CTMRG, DIR)

  COMPLEX(KIND=DP),  ALLOCATABLE,   INTENT(INOUT)   :: N_norm(:,:,:,:)     !Norm matrix
  COMPLEX(KIND=DP),                 INTENT(IN)      :: Xa(:,:), Yb(:,:)    !X,Y residual nodes of factorized PEPS
  TYPE(block_mpo),                  INTENT(IN)      :: TN2D(:)             !2D network
  TYPE(ctm_corner_type),            INTENT(IN)      :: Cmat(:)             !CTM Cmat
  TYPE(ctm_transfer_type),          INTENT(IN)      :: Tmat(:)             !CTM Tmat
  TYPE(ctm_params),                 INTENT(IN)      :: CTMRG               !CTM params
  CHARACTER(LEN=*),                 INTENT(IN)      :: DIR                 !Bond we're updating

  !! Environment of 2D network && A,B-halfs of Norm matrix
  TYPE(block_mpo),  ALLOCATABLE  :: ENV2D(:)
  COMPLEX(KIND=DP), ALLOCATABLE  :: NhalfAA(:,:,:,:),  NhalfBB(:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE  :: NhalfAA_3D(:,:,:), NhalfBB_3D(:,:,:)

  !! (1) Construct environment of a single bond (given by DIR)
  CALL calc_PsiPsi_single_bond_environment(ENV2D, TN2D, Cmat, Tmat, CTMRG, DIR)

  !! (2) Calculate A-half && B-half of norm matrix 
  !!     (i.e. connected to site-A && site-B of iPEPS respectively)
  SELECT CASE(DIR)
  CASE('S', 'N')  !! S,N dirs (legsA = '1,234,5'; legsB = '1,4,235')

       NhalfAA = TENSMUL(ENV2D(3) % m,  TENSMUL(ENV2D(1) % m,  ENV2D(5) % m,  MULT='34', FUSE='(11,22)'),  MULT='43', FUSE='(11,22)')
       NhalfBB = TENSMUL(ENV2D(4) % m,  TENSMUL(ENV2D(2) % m,  ENV2D(6) % m,  MULT='43', FUSE='(11,22)'),  MULT='34', FUSE='(11,22)')

  CASE('E', 'W')  !! E,W dirs (legsA = '1,245,3'; legsB = '1,2,345')

       NhalfAA = TENSMUL(TENSMUL(ENV2D(5) % m,  ENV2D(1) % m,  MULT='43', FUSE='(11,22)'),  ENV2D(3) % m,  MULT='34', FUSE='(11,22)')
       NhalfBB = TENSMUL(TENSMUL(ENV2D(6) % m,  ENV2D(2) % m,  MULT='34', FUSE='(11,22)'),  ENV2D(4) % m,  MULT='43', FUSE='(11,22)')

  CASE DEFAULT
       CALL invalid_flag("calc_norm_matrix -- invalid DIR ", DIR)
  end SELECT

  !! (3) Contract A,B-halfs with Xa,Yb residual nodes of factorized iPEPS
  NhalfAA = TENSMUL(TENSMUL(NhalfAA,  Xa,  '11'),  CONJG(Xa),  '21')
  NhalfBB = TENSMUL(TENSMUL(NhalfBB,  Yb,  '12'),  CONJG(Yb),  '22') 

  !! (4) Contract A,B-halfs together --> gives a full norm matrix
  NhalfAA_3D = RESHAPE_3D(              NhalfAA,        '1,2,34')
  NhalfBB_3D = RESHAPE_3D(TensTRANSPOSE(NhalfBB, '34'), '1,2,34')
  N_norm     = TENSMUL(NhalfAA_3D,  NhalfBB_3D,  '33')

 END SUBROUTINE calc_norm_matrix





 !!! Construct single bond environment in PSI-PSI network 
 !!! -- calculate single bond environment && reshape it from MPS to MPO
 !!!    with two local dims (one connected to |PSI> and one to <PSI|)
 SUBROUTINE calc_PsiPsi_single_bond_environment(envMPO, TN2D, Cmat, Tmat, CTMRG, DIR)

  TYPE(block_mpo),    ALLOCATABLE,  INTENT(INOUT)   :: envMPO(:)           !Reshaped environment
  TYPE(block_mpo),                  INTENT(IN)      :: TN2D(:)             !2D network
  TYPE(ctm_corner_type),            INTENT(IN)      :: Cmat(:)             !CTM Cmat
  TYPE(ctm_transfer_type),          INTENT(IN)      :: Tmat(:)             !CTM Tmat
  TYPE(ctm_params),                 INTENT(IN)      :: CTMRG               !CTM params
  CHARACTER(LEN=*),                 INTENT(IN)      :: DIR                 !Bond we're updating

  !! Environment of 2D network (before reshaping)
  TYPE(block_mps), ALLOCATABLE  :: envMPS(:)

  !! Dims && indices
  INTEGER :: site, N_sites
  INTEGER :: LocalDim

  !! Construct environment of a single bond (given by DIR)
  CALL calc_single_bond_environment(envMPS, Cmat, Tmat, TN2D, CTMRG, DIR, pbc=.TRUE.)

  !! Size of envMPS
  N_sites  = SIZE(envMPS)

  !! Check if proposed reshaping gives consistent local dims (i.e. envMPS(site) % SNdim = LocalDim**2)
  DO site=1,N_sites
     LocalDim = INT(SQRT(1.0d0 * envMPS(site) % SNdim))
     CALL check_sizes_equal(LocalDim**2,  envMPS(site) % SNdim, "reshape_environment_to_PsiPsi -- inconsistent LocalDim ")
  end DO

  !! Reshape MPS-environment into MPO-environment (two local dims connected to |Psi> and <Psi|)
  CALL allocate_empty_mpo_block(envMPO, N_sites)
  DO site=1,N_sites
     CALL allocate_mpo_site(envMPO(site), RESHAPE_4D(envMPS(site) % m, '12,3,4', (/LocalDim, LocalDim/)))
  end DO

 END SUBROUTINE calc_PsiPsi_single_bond_environment






 !!! Calculate iPEPS R-matrix in R|v> = |S> optimization equation !!!
 SUBROUTINE calc_ipeps_Rmat(Rmat, V_site, N_norm, whichSITE)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: Rmat(:,:,:,:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: V_site(:,:,:) 
  COMPLEX(KIND=DP),              INTENT(IN)    :: N_norm(:,:,:,:)
  CHARACTER(LEN=*),              INTENT(IN)    :: whichSITE

  !! Intermediate tensors
  COMPLEX(KIND=DP), ALLOCATABLE :: VV(:,:,:), TmpN(:,:,:)

  SELECT CASE(whichSITE)
  CASE('A')

     !! Contract together V and V^HC, reshape Norm to be compatible with (V * V^HC)
     VV   = TensTRANSPOSE(RESHAPE_3D(TENSMUL(V_site,  CONJG(V_site), '11'), '1,3,24'),  '23') 
     TmpN =               RESHAPE_3D(N_norm,                                '1,3,24')

     !! Compute R = Norm * (V * V^HC)
     Rmat  = TENSMUL(TmpN, VV, '32')

  CASE('B')

     !! Contract together V and V^HC, reshape Norm to be compatible with (V * V^HC)
     VV   = TensTRANSPOSE(RESHAPE_3D(TENSMUL(V_site,  CONJG(V_site), '11'), '2,13,4'),  '23')
     TmpN =               RESHAPE_3D(N_norm,                                '2,13,4')

     !! Compute R = (V * V^HC) * Norm
     Rmat  = TENSMUL(VV, TmpN, '32')

  CASE DEFAULT
       CALL invalid_flag("calc_ipeps_Rmat -- invalid whichSITE ", whichSITE)
  end SELECT

 END SUBROUTINE calc_ipeps_Rmat






 !!! Calculate iPEPS S-vector for PSI-TN 
 !!! -- i.e. norm network with one site missing 
 !!! -- equivalent to iPEPS derivative wrt that site 
 SUBROUTINE calc_ipeps_Svec(Svec, AA, BB, VV, N_norm, whichSITE)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: Svec(:,:,:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: AA(:,:,:), BB(:,:,:) 
  COMPLEX(KIND=DP),              INTENT(IN)    :: VV(:,:,:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: N_norm(:,:,:,:)
  CHARACTER(LEN=*),              INTENT(IN)    :: whichSITE

  !! Intermediate tensors
  COMPLEX(KIND=DP), ALLOCATABLE :: SitesAB(:,:,:), Norm3D(:,:,:), TEMP(:,:,:,:)

  !! Contract A,B-sites of iPEPS with norm matrix --> gives TEMP intermediate tensor
  SitesAB = RESHAPE_3D(TENSMUL(AA, BB, '32'), '1,2,34')
  Norm3D  = TensTRANSPOSE(TensTRANSPOSE(RESHAPE_3D(N_norm, '12,3,4'),  '13'),  '12')
  TEMP    = TensTRANSPOSE(TENSMUL(SitesAB, Norm3D, '33'), '23')

  SELECT CASE(whichSITE)
  CASE('A')

     !! Contract VV with TEMP tensor --> gives Svec (equiv to derivative wrt missing site adjacent to VV)
     Svec = TENSMUL(RESHAPE_3D(TEMP, '1,3,24'), RESHAPE_2D(CONJG(VV), '2,13'), '32')

  CASE('B')

     !! Contract VV with TEMP tensor --> gives Svec (equiv to derivative wrt missing site adjacent to VV)
     Svec = TENSMUL(RESHAPE_3D(TEMP, '2,13,4'), RESHAPE_2D(CONJG(VV), '12,3'), '21')

  CASE DEFAULT
       CALL invalid_flag("calc_ipeps_Svec -- invalid whichSITE ", whichSITE)
  end SELECT

 END SUBROUTINE calc_ipeps_Svec




 



 !!! Fix gauge of norm matrix in iPEPS PSI-PSI network !!!
 SUBROUTINE gauge_fix_ipeps_norm(mps, Xa, Yb, N_norm)

  TYPE(block_mps),  ALLOCATABLE, INTENT(INOUT) :: mps(:)           !Reduced PEPS (takes form of MPS sites)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: Xa(:,:), Yb(:,:) !X,Y residual nodes of factorized PEPS
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: N_norm(:,:,:,:)  !Norm matrix

  !! Z matrix that orthodecomposes N_norm, and Matrices that factorize Z
  COMPLEX(KIND=DP), ALLOCATABLE :: Z(:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: Ul(:,:,:), L(:,:), Linv(:,:) 
  COMPLEX(KIND=DP), ALLOCATABLE :: Ur(:,:,:), R(:,:), Rinv(:,:)

  !! (1) Orthodecompose N_norm = Z * ZH
  CALL orthodecompose_PsiPsi_ipeps_norm(Z, N_norm)

  !! (2) Factorize Z = Ur * R, and Z = L * Ul
  CALL LHS_factorize_Tens3D(Ur, R, Z, '12,3') 
  CALL RHS_factorize_Tens3D(Ul, L, Z, '2,13') 

  Rinv = invert_matrix(R)
  Linv = invert_matrix(L)

  !! (3A) Transform Z -> Linv * Z * Rinv
  Z = TENSMUL(TENSMUL(Z, Linv, '22'),  Rinv, '31')

  !! (3B) Calculate new norm N = Z * ZH (from transformed Z)
  CALL copyTens(N_norm, TENSMUL(Z, CONJG(Z), '11'))

  !! (4) Transform mps sites: mps(1) -> L * mps(1), mps(2) -> mps(2) * R 
  CALL allocate_mps_site(mps(1),  TENSMUL(mps(1) % m,  L, '21'))
  CALL allocate_mps_site(mps(2),  TENSMUL(mps(2) % m,  R, '32'))

  !! (5) Transform residual sites Xa -> Linv * Xa, Yb -> Yb * Rinv
  CALL copyTens(Xa,  TENSMUL(Xa,  Linv, '22'))
  CALL copyTens(Yb,  TENSMUL(Yb,  Rinv, '11'))

  !! (6) Hermitize N_norm:
  N_norm = 0.5D0*(N_norm + TensTRANSPOSE(TensTRANSPOSE(CONJG(N_norm), '13'), '24'))

 END SUBROUTINE gauge_fix_ipeps_norm





 !!! Decompose norm matrix into orthonormal components: N_norm = Z Z^HC !!!
 SUBROUTINE orthodecompose_PsiPsi_ipeps_norm(Z, N_norm)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: Z(:,:,:)         !Z --> from factorization N=Z*ZH
  COMPLEX(KIND=DP),              INTENT(IN)    :: N_norm(:,:,:,:)  !Norm matrix

  !! Dims && indices
  INTEGER :: dimN(4), i

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: theta(:,:), U(:,:), VH(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)
  INTEGER                       :: chi

  !! Get dims of N_norm
  dimN = shape(N_norm)

  !! (1) Construct theta && impose hermiticity
  theta = RESHAPE_2D(N_norm, '12,34')
  CALL purify_hermitian_tens(theta)

  !! (2) Calculate SVD (We've VH = UH cause theta = hermitian)
  chi=-1
  CALL compute_lapack_svd(theta, U, VH, Sigma, chi, 2.0d0)

  !! (3A) Impose positivity: set all small -ve/+ve sigmas to zero
  DO i=1,SIZE(Sigma)
     IF(Sigma(i) .LT. 1.0D-08) Sigma(i) = 0.0D0
  end DO

  !! (3B) Absorb Sigmas into U, VH
  CALL LambdaMUL(U, complx(Sigma), VH)

  !! (4) Reshape U*SQRT(Sigma) into Z --> we've factorization N = Z * ZH
  Z = RESHAPE_3D(U, '1,23', dimN(1:2))

 END SUBROUTINE orthodecompose_PsiPsi_ipeps_norm


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (4) ALS optimization of ipeps !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE ALS_optimize_ipeps(mps, mpsOO, N_norm)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps(:)
  TYPE(block_mps),              INTENT(IN)    :: mpsOO(:)
  COMPLEX(KIND=DP),             INTENT(IN)    :: N_norm(:,:,:,:) 

  !! Distance, CG tolerance, sweeping params 
  COMPLEX(KIND=DP)         :: dist, dist_old, dist_00
  REAL(KIND=DP), PARAMETER :: tolSweep = 1.0D-10
  INTEGER,       PARAMETER :: N_iter   = 1000
  INTEGER,       PARAMETER :: min_iter = 10
  INTEGER                  :: iter

  !! Find initial distance
  CALL calc_ipeps_distance(dist, mps, mpsOO, N_norm)
  dist_old = dist
  dist_00  = dist

  IF(ABS(dist_00) .LT. tolSweep) dist_00 = 1.0D0
  
  !! Main iteration (alternating sweep)
  DO iter=1,N_iter

     !WRITE(*,*)
     !WRITE(*,*) "ALS iter = ", iter, " dist = ", dist

     !! Optimize sites-A,B of reduced iPEPS
     CALL do_ALS_step(mps, mpsOO, N_norm, 'A')
     CALL do_ALS_step(mps, mpsOO, N_norm, 'B')

     !! Calculate the distance (cost function)
     CALL calc_ipeps_distance(dist, mps, mpsOO, N_norm)

     !! Check convergence
     IF(iter .GE. min_iter) THEN
        IF(ABS(dist - dist_old)/ABS(dist_00) .LT. tolSweep) THEN
           WRITE(*,*) "ALS has converged, N = ", iter, " dist = ", dist
           EXIT 
        ELSEIF(iter .EQ. N_iter) THEN
           WRITE(*,*) "ALS has failed to converge, dist = ", dist, dist_old, dist_00
           STOP
        end IF
     end IF

     !! Update old distance
     dist_old = dist
  end DO

 END SUBROUTINE ALS_optimize_ipeps  





 !!! Do a single step of ALS sweep !!!
 SUBROUTINE do_ALS_step(mps, mpsOO, N_norm, whichSITE)

  USE lapack_wrappers, ONLY: lapack_gen_linsol

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps(:)
  TYPE(block_mps),              INTENT(IN)    :: mpsOO(:)
  COMPLEX(KIND=DP),             INTENT(IN)    :: N_norm(:,:,:,:) 
  CHARACTER(LEN=*),             INTENT(IN)    :: whichSITE

  COMPLEX(KIND=DP), ALLOCATABLE :: Rmat(:,:,:,:), Svec(:,:,:)
  INTEGER                       :: sA, sB

  COMPLEX(KIND=DP), ALLOCATABLE :: tmpx(:,:)
  INTEGER                       :: dimR(4), dimX(3)

  !! Get sites (site-A = the one we're optimizing; site-B = the one that is fixed)
  SELECT CASE(whichSITE)
  CASE('A')
       sA=1; sB=2
  CASE('B')
       sA=2; sB=1
  CASE DEFAULT
       CALL invalid_flag("do_ALS_step -- invalid whichSITE ", whichSITE)
  end SELECT

  !! (1) Calculate R-matrix && S-vector
  CALL calc_ipeps_Rmat(Rmat,                             mps(sB) % m, N_norm, whichSITE)
  CALL calc_ipeps_Svec(Svec, mpsOO(1) % m, mpsOO(2) % m, mps(sB) % m, N_norm, whichSITE)

  !! (2) Hermitize R-matrix 
  Rmat = 0.5D0 * (Rmat + TensTRANSPOSE(TensTRANSPOSE(CONJG(Rmat), '13'), '24'))

  !! (3) Prepare temp matrix to store the solution X
  dimX = shape(mps(sA) % m) 
  ALLOCATE(tmpx(dimX(2) * dimX(3), dimX(1)))

  !! (4) Solve linear system R*X = S, treating R, X, S as matrices
  CALL lapack_gen_linsol(tmpx,  TensTRANSPOSE(RESHAPE_2D(Rmat, '12,34')),  TensTRANSPOSE(RESHAPE_2D(Svec, '1,23')))
  
  !! (5) Write the solution X to reduced iPEPS
  CALL allocate_mps_site(mps(sA),  RESHAPE_3D(TensTRANSPOSE(tmpx), '1,23', (/dimX(2), dimX(3)/)))    

 END SUBROUTINE do_ALS_step




 !!! Calculate cost function -- distance between optimal and approximate iPEPS !!!
 SUBROUTINE calc_ipeps_distance(dist, mps, mpsOO, N_norm)

  COMPLEX(KIND=DP), INTENT(OUT) :: dist
  TYPE(block_mps),  INTENT(IN)  :: mps(:)
  TYPE(block_mps),  INTENT(IN)  :: mpsOO(:)
  COMPLEX(KIND=DP), INTENT(IN)  :: N_norm(:,:,:,:) 

  COMPLEX(KIND=DP)              :: dAA, dOO, dAO
  COMPLEX(KIND=DP), ALLOCATABLE :: vecAA(:,:,:), vecOO(:,:,:), vecAO(:,:,:)

  CALL calc_ipeps_Svec(vecAA,    mps(1) % m,    mps(2) % m,    mps(2) % m,  N_norm,  whichSITE='A')
  CALL calc_ipeps_Svec(vecOO,  mpsOO(1) % m,  mpsOO(2) % m,  mpsOO(2) % m,  N_norm,  whichSITE='A')
  CALL calc_ipeps_Svec(vecAO,  mpsOO(1) % m,  mpsOO(2) % m,    mps(2) % m,  N_norm,  whichSITE='A')

  dAA = SUM(RESHAPE_1D(CONJG(  mps(1) % m)) * RESHAPE_1D(vecAA))
  dOO = SUM(RESHAPE_1D(CONJG(mpsOO(1) % m)) * RESHAPE_1D(vecOO))
  dAO = SUM(RESHAPE_1D(CONJG(  mps(1) % m)) * RESHAPE_1D(vecAO))

  dist = dAA + dOO - dAO - CONJG(dAO)  

 END SUBROUTINE calc_ipeps_distance




 !!! Find initial guess for iPEPS CG iteration !!!
 SUBROUTINE initialize_ALS_ipeps(mps, chi, eps)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps(:)
  INTEGER,                      INTENT(IN)    :: chi
  REAL(KIND=DP),                INTENT(IN)    :: eps

  !! Lambda produced by SVD on reduced PEPS sites
  TYPE(block_lambda) :: lambda

  !! Do SVD on reduced iPEPS sites multiplied by propagator 
  !! --> gives SU result that we'll use as a starting pt for FU
  CALL compute_svd_of_imps_bond(mps, lambda, chi, eps, '12')

  !! Absorb lambda into reduced iPEPS sites 
  CALL LambdaMUL(mps(1) % m, lambda % m, mps(2) % m)

  !WRITE(*,*)
  !WRITE(*,*) "init ALS -- lambda "
  !CALL printVector(lambda % m)
  !WRITE(*,*)

 END SUBROUTINE initialize_ALS_ipeps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (5) Routines for testing iPEPS norm !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Contract iPEPS network using norm matrix !!!
 SUBROUTINE simple_contract_norm_network(C_norm, N_norm, mps)

  COMPLEX(KIND=DP),  INTENT(INOUT) :: C_norm
  COMPLEX(KIND=DP),  INTENT(IN)    :: N_norm(:,:,:,:)  !Norm matrix
  TYPE(block_mps),   INTENT(IN)    :: mps(:)           !Reduced PEPS (takes form of MPS sites)
 
  COMPLEX(KIND=DP), ALLOCATABLE :: overlap_11(:,:), overlap_22(:,:)

  overlap_11 = TENSMUL(mps(1) % m, CONJG(mps(1) % m), MULT='11', FUSE='(22,33)')
  overlap_22 = TENSMUL(mps(2) % m, CONJG(mps(2) % m), MULT='11', FUSE='(22,33)')

  C_norm = TRACE(TENSMUL(TENSMUL(overlap_11, overlap_22), TensTRANSPOSE(RESHAPE_2D(N_norm, '13,24'))))

 END SUBROUTINE simple_contract_norm_network





 !!! Test the norm of iPEPS network -- contract using CTM tensors, N_norm, Rmat, Svec !!!
 SUBROUTINE test_ipeps_norm(ipeps, Cmat, Tmat, CTMRG, DIR)

  USE imaginary_observables_ipeps, ONLY: calc_CTM_environment_ipeps, construct_reduced_rho, contract_bare_TN

  TYPE(block_peps),        ALLOCATABLE,  INTENT(INOUT)  :: ipeps(:)        !iGamma of PEPS
  TYPE(ctm_corner_type),   ALLOCATABLE,  INTENT(INOUT)  :: Cmat(:)         !CTM Cmat
  TYPE(ctm_transfer_type), ALLOCATABLE,  INTENT(INOUT)  :: Tmat(:)         !CTM Tmat
  TYPE(ctm_params),                      INTENT(IN)     :: CTMRG           !CTM params
  CHARACTER(LEN=*),                      INTENT(IN)     :: DIR             !Bond we're propagating

  !! Reduced iPEPS reduced tensors --> MPS sites && Residual tensors
  TYPE(block_mps),    ALLOCATABLE :: mps(:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: Xa(:,:),  Yb(:,:)

  !! Norm matrix && reduced 2D network
  COMPLEX(KIND=DP),   ALLOCATABLE :: N_norm(:,:,:,:)
  TYPE(block_mpo),    ALLOCATABLE :: TN2D(:)
  COMPLEX(KIND=DP)                :: C_norm

  !! Auxiliary vars
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_onesite(:,:),  rho_twosite(:,:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: Svec(:,:,:),       Rmat(:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: mps_tmp(:,:),      Rtmp(:,:)
  COMPLEX(KIND=DP)              :: dN, dS, dR, rho_norm


  !! Calculate reduced (factorized) iPEPS
  CALL calc_factorized_ipeps(mps, Xa, Yb, ipeps, DIR)

  !! Calculate iPEPS norm matrix
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)
  CALL calc_norm_matrix(N_norm, Xa, Yb, TN2D, Cmat, Tmat, CTMRG, DIR)

  !! Normalize N_norm
  CALL simple_contract_norm_network(C_norm, N_norm, mps)
  N_norm = N_norm/C_norm
  N_norm = 0.5D0*(N_norm + TensTRANSPOSE(TensTRANSPOSE(CONJG(N_norm), '13'), '24'))
  

  !! (1) Norm using trace
  CALL construct_reduced_rho(rho_onesite, rho_twosite, Cmat, Tmat, ipeps, CTMRG, max_sep=1)
  rho_norm = TRACE(rho_onesite)

  !! (2) Norm using CTM
  CALL contract_bare_TN(C_norm, Cmat, Tmat, ipeps, CTMRG, DIR)

  !! (3) Norm using N_norm
  CALL simple_contract_norm_network(dN,  N_norm,  mps)

  !! (4) Norm using Svec
  CALL calc_ipeps_Svec(Svec,    mps(1) % m,    mps(2) % m,    mps(2) % m,  N_norm,  whichSITE='A')

  dS = SUM(RESHAPE_1D(CONJG(mps(1) % m)) * RESHAPE_1D(Svec))

  !! (5) Norm using Rmat
  CALL calc_ipeps_Rmat(Rmat,    mps(2) % m,                                N_norm,  whichSITE='A')

  Rtmp    = TensTRANSPOSE(RESHAPE_2D(Rmat,       '12,34'))
  mps_tmp = TensTRANSPOSE(RESHAPE_2D(mps(1) % m, '1,23' ))
  dR      = TRACE(TENSMUL(TensTRANSPOSE(mps_tmp, 'HC'),  TENSMUL(Rtmp,  mps_tmp)))

  WRITE(*,*)
  WRITE(*,*) "TN norms -- C_norm ", C_norm
  WRITE(*,*) "TN norms --     dN ", dN
  WRITE(*,*) "TN norms --     dR ", dR
  WRITE(*,*) "TN norms --     dS ", dS 
  WRITE(*,*)

 END SUBROUTINE test_ipeps_norm

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE imaginary_full_update_ipeps
