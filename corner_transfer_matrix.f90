MODULE corner_transfer_matrix

 USE utility
 USE definitions_mps_mpo
 USE mps_mpo_utility
 USE ctm_definitions

 !!!!!!!!!!!!!!!!!!!!!!! CONTENTS: !!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! (1) CTMRG TOP-LEVEL ROUTINES (to manage the iteration process)
 !!
 !! (2) CTMRG SINGLE MOVE (one- && two- directional)
 !!
 !! (3) ABSORPTION, CONTRACTION, RENORMALIZATION
 !!
 !! (4) CALCULATE CTM ISOMETRY
 !!
 !! (5) CHECK CTM CONVERGENCE
 !!
 !! (6) RESCALE/NORMALIZE CTM TENSORS 
 !!
 !! (7) AUXILIARY ROUTINES: CHECK IF CTM HAS SINGLETON DIMS
 !!
 !! (8) CTM ALGEBRA -- LOW LEVEL ROUTINES 
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IMPLICIT NONE

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (1) CTMRG TOP-LEVEL ROUTINES (to manage the iteration process) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Main subroutine of CTMRG 
 !!! -- run CTMRG calculation to compute boundaries of an infinite 2D network
 SUBROUTINE compute_CTMRG_boundaries(Cmat, Tmat, TN2D, CTMRG, use_old_ctm)

  USE datafile_utility, ONLY: setupConvDatafile

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2D(:)
  TYPE(ctm_params),                     INTENT(IN)    :: CTMRG
  LOGICAL,                              INTENT(IN)    :: use_old_ctm

  !! CTM iteration && convergence params
  INTEGER            :: iter
  INTEGER, PARAMETER :: N_iters = 30000
  INTEGER            :: datafile
  LOGICAL            :: is_converged

  !! Singular vals of CTM
  COMPLEX(KIND=DP), ALLOCATABLE :: svals1(:), svals2(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: svals3(:), svals4(:)

  !! Initialize CTM (random network)
  CALL initialize_CTM(Cmat, Tmat, TN2D, chi=(/2,2,2,2/), init_rand_state='RAND', use_old_ctm=use_old_ctm)

  !! Setup file for tracking convergence
  CALL setupConvDatafile(datafile=datafile, descr="CTMRG", valstr="Max_Sval")

  !! Initialize convergence tester
  is_converged = .FALSE.

  mainloop: DO iter=1,N_iters

     !! Perform a single iteration of CTMRG
     CALL perform_CTMRG_iteration(Cmat, Tmat, TN2D, CTMRG)

     !! Check convergence
     CALL check_CTM_convergence(is_converged, svals1, svals2, svals3, svals4, Cmat, CTMRG, datafile, iter)

     IF(is_converged) THEN
        CALL print_ctm_final_message("CTMRG DONE -- ",       Cmat, iter, use_old_ctm, svals1, svals2, svals3, svals4)
        EXIT mainloop
     ELSEIF(iter .EQ. N_iters) THEN
        CALL print_ctm_final_message("CTMRG HAS FAILED -- ", Cmat, iter, use_old_ctm, svals1, svals2, svals3, svals4)
        STOP 
     end IF 
  
  end DO mainloop

  !! Close files
  CLOSE(datafile)

 END SUBROUTINE compute_CTMRG_boundaries




 !!! Perform a single iteration of two-directional CTMRG (consists of 4 moves) !!!
 SUBROUTINE perform_CTMRG_iteration(Cmat, Tmat, TN2Din, CTMRG)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2Din(:)
  TYPE(ctm_params),                     INTENT(IN)    :: CTMRG

  !! Local copy of 2D network (so we can swap sites after each CTM move)
  TYPE(block_mpo), ALLOCATABLE :: TN2D(:)

  !! Create a local copy of TN2D
  CALL copy_mpo_block(TN2D, TN2Din)

  !! Perform CTMRG iteration FIXME we've swapped W,E and S,N moves
  CALL perform_CTM_move_twoDir(Cmat, Tmat, TN2D, CTMRG, DIR = (/'W', 'E'/))
  CALL perform_CTM_move_twoDir(Cmat, Tmat, TN2D, CTMRG, DIR = (/'N', 'S'/))

  CALL perform_CTM_move_twoDir(Cmat, Tmat, TN2D, CTMRG, DIR = (/'W', 'E'/))
  CALL perform_CTM_move_twoDir(Cmat, Tmat, TN2D, CTMRG, DIR = (/'N', 'S'/))

  !! Normalize CTM to ensure numerical stability
  CALL rescale_CTM_tensors(Cmat, Tmat)

 END SUBROUTINE perform_CTMRG_iteration



 !!! Perform a single iteration of one-directional CTMRG (consists of 2 x 4 moves) !!!
 SUBROUTINE perform_CTMRG_iteration_oneDir(Cmat, Tmat, TN2Din, CTMRG)

  TYPE(ctm_corner_type),   INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),         INTENT(IN)    :: TN2Din(:)
  TYPE(ctm_params),        INTENT(IN)    :: CTMRG

  !! Local copy of 2D network (so we can swap sites after each CTM move)
  TYPE(block_mpo), ALLOCATABLE :: TN2D(:)

  !! Create a local copy of TN2D
  CALL copy_mpo_block(TN2D, TN2Din)

  !! Perform CTMRG iteration (NB. horizontal and vertical dirs must alternate, e.g. W-N-E-S, but NOT W-E-S-N)
  CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR='W')
  CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR='W')

  CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR='N')
  CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR='N')

  CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR='E')
  CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR='E')

  CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR='S')
  CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR='S')

  !! Normalize CTM to ensure numerical stability
  CALL rescale_CTM_tensors(Cmat, Tmat) 

 END SUBROUTINE perform_CTMRG_iteration_oneDir

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (2) CTMRG SINGLE MOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Perform a single one-directional move of CTM !!!
 SUBROUTINE perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR)

  TYPE(ctm_corner_type),   INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),         INTENT(INOUT) :: TN2D(:)
  TYPE(ctm_params),        INTENT(IN)    :: CTMRG
  CHARACTER(LEN=*),        INTENT(IN)    :: DIR

  !! Isometries for renormalization
  COMPLEX(KIND=DP),  ALLOCATABLE   :: Z(:,:), ZP(:,:)
  COMPLEX(KIND=DP),  ALLOCATABLE   :: W(:,:), WP(:,:)

  !! Rmat super corner matrices
  TYPE(ctm_corner_type), ALLOCATABLE :: Rmat(:)

  !! (1) Insert && Absorb extra block into CTM
  CALL absorb_extra_block_into_CTM(Cmat, Tmat, TN2D, DIR)

  !! (2A) Calc isometry with C-matrices (corners)
  CALL calc_CTM_isometry_oneDir(Z, ZP, Cmat, CTMRG, DIR)  

  !! (2B) Calc isometry with R-matrices (super corners)
  CALL copy_CTM_corner_mat(Rmat, Cmat)
  CALL self_contract_CTM_edge(Rmat, Tmat, DIR)
  CALL calc_CTM_isometry_oneDir(W, WP, Rmat, CTMRG, DIR) 

  !! (3) Renormalize CTM
  CALL renormalize_CTM(Cmat, Tmat, Z, ZP, W, WP, DIR)

  !! (4) The move swaps [lateral sites of Tmat] && [sites of 2D network]
  !!     -- i.e. sites on CTM blocks we've inserted
  CALL swap_CTM_block_sites(Tmat, TN2D, DIR) 

  !! (5) Check CTM dims
  CALL check_CTM_dims(Cmat, Tmat, TN2D, DIR)

 END SUBROUTINE perform_CTM_move_oneDir




 !!! Perform a single two-directional move of CTM !!!
 SUBROUTINE perform_CTM_move_twoDir(Cmat, Tmat, TN2D, CTMRG, DIR)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(INOUT) :: TN2D(:)
  TYPE(ctm_params),                     INTENT(IN)    :: CTMRG
  CHARACTER(LEN=*),                     INTENT(IN)    :: DIR(2)

  !! Orthogonal dir
  CHARACTER(LEN=1) :: orthDIR(2)

  !! Isometries for renormalization
  COMPLEX(KIND=DP),  ALLOCATABLE :: Z_11(:,:), ZP_11(:,:), Z_22(:,:), ZP_22(:,:)
  COMPLEX(KIND=DP),  ALLOCATABLE :: W_11(:,:), WP_11(:,:), W_22(:,:), WP_22(:,:)

  !! Updated CTM corner matrices, and CTM halfs
  TYPE(ctm_corner_type), ALLOCATABLE :: Rmat(:)
  COMPLEX(KIND=DP),      ALLOCATABLE :: HA_11(:,:), HB_11(:,:)
  COMPLEX(KIND=DP),      ALLOCATABLE :: HA_22(:,:), HB_22(:,:)

  INTEGER :: sA, sB, sC, sD
  INTEGER :: subA, subB

  !! Get the order of CTM sites on the edge specified by dirFlag
  CALL get_CTM_sites(sA=sA, sB=sB, sC=sC, sD=sD, subA=subA, subB=subB, DIR=DIR(1))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! (1) Insert && Absorb two extra blocks into CTM (two-directional absorption)
  CALL absorb_extra_block_into_CTM(Cmat, Tmat, TN2D, DIR(1))
  CALL absorb_extra_block_into_CTM(Cmat, Tmat, TN2D, DIR(2))

  !! (2) Construct R-matrices -- supercorners (Cmat --> Rmat)
  CALL copy_CTM_corner_mat(Rmat, Cmat)
  CALL self_contract_CTM_edge(Rmat, Tmat, DIR(1))
  CALL self_contract_CTM_edge(Rmat, Tmat, DIR(2))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! (3) Compute W-isometries for the middle cut
  CALL calc_CTM_isometry_twoDir(W_11, WP_11, Rmat, CTMRG, DIR(1))
  CALL calc_CTM_isometry_twoDir(W_22, WP_22, Rmat, CTMRG, DIR(2))

  !! (4A) Renormalize R-matrices before applying a shift
  orthDIR = orthogonalDIR(DIR)
  CALL renormalize_CTM_supercorners(Rmat, CTMRG, DIR=orthDIR(1))
  CALL renormalize_CTM_supercorners(Rmat, CTMRG, DIR=orthDIR(2))

  !! (4B) Apply a translationally-invariant shift to R-matrices 
  CALL self_contract_CTM_edge(Rmat, Tmat, DIR(1), is_inverse_order=.TRUE.)
  CALL self_contract_CTM_edge(Rmat, Tmat, DIR(2), is_inverse_order=.TRUE.)

  !! (4C) Compute isometries for the top/bottom cuts 
  CALL calc_CTM_isometry_twoDir(Z_11, ZP_11, Rmat, CTMRG, DIR(1))
  CALL calc_CTM_isometry_twoDir(Z_22, ZP_22, Rmat, CTMRG, DIR(2))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! (4) Renormalize CTM
  CALL renormalize_CTM(Cmat, Tmat, Z_11, ZP_11, W_11, WP_11, DIR(1))
  CALL renormalize_CTM(Cmat, Tmat, Z_22, ZP_22, W_22, WP_22, DIR(2))

  !! (5) The move swaps [lateral sites of Tmat] && [sites of 2D network]
  !!     -- i.e. sites on CTM blocks we've inserted
  CALL swap_CTM_block_sites(Tmat, TN2D, DIR(1))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! (6) Check CTM dims
  CALL check_CTM_dims(Cmat, Tmat, TN2D, DIR(1))
  CALL check_CTM_dims(Cmat, Tmat, TN2D, DIR(2))

 END SUBROUTINE perform_CTM_move_twoDir

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (3) ABSORPTION, CONTRACTION, RENORMALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Insert && Absorb an extra block of tensors into CTM !!!
 SUBROUTINE absorb_extra_block_into_CTM(Cmat, Tmat, TN2D, DIR)

  TYPE(ctm_corner_type),   INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),         INTENT(IN)    :: TN2D(:)
  CHARACTER(LEN=*),        INTENT(IN)    :: DIR

  INTEGER :: sA, sB, sD
  INTEGER :: subA, subB

  !! Get the order of CTM sites on the edge specified by dirFlag
  CALL get_CTM_sites(sA=sA, sB=sB, sD=sD, subA=subA, subB=subB, DIR=DIR)

  !! Multiply sites -- absorb extra block into CTM edge
  CALL mult_Cmat_Tmat(Cmat(sA),            Tmat(sD) % T(subA),    DIR='22')
  CALL mult_Tmat_TN  (Tmat(sA) % T(subA),  TN2D(subA),            DIR=DIR)
  CALL mult_Tmat_TN  (Tmat(sA) % T(subB),  TN2D(subB),            DIR=DIR)
  CALL mult_Cmat_Tmat(Cmat(sB),            Tmat(sB) % T(subB),    DIR='31')

 END SUBROUTINE absorb_extra_block_into_CTM




 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! Self contract an edge of CTM -- absorb T-tens into C-tens on that edge
 !!! to create an edge that consists of C-tens only. 
 !!!
 !!! NB. To avoid confusion: this operation is the same during both half-moves 'A,B' 
 !!!     since it involves only the edge itself ant not extra inserted blocks.
 !!!
 !!! NB. Consider renaming to cornerize_CTM_edge (or cornerize_[long/short]_CTM_edge)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE self_contract_CTM_edge(Cmat, Tmat, DIR, is_inverse_order)

  TYPE(ctm_corner_type),   INTENT(INOUT)        :: Cmat(:)
  TYPE(ctm_transfer_type), INTENT(INOUT)        :: Tmat(:)
  CHARACTER(LEN=*),        INTENT(IN)           :: DIR
  LOGICAL,                 INTENT(IN), OPTIONAL :: is_inverse_order

  INTEGER :: sA, sB, subA, subB
  LOGICAL :: is_inv

  !! Whether to absorb sites in inverse order (default to FALSE)
  is_inv = OptArg(is_inverse_order, .FALSE.)

  !! Get the order of CTM sites on the edge specified by dirFlag
  CALL get_CTM_sites(sA=sA, sB=sB, subA=subA, subB=subB, DIR=DIR)

  !! Swap subA, subB
  IF(is_inv) CALL swap(subA, subB)

  !! Multiply sites -- absorb transfer sites into corner sites
  CALL mult_Cmat_Tmat(Cmat(sA),  Tmat(sA) % T(subA),   DIR='31')
  CALL mult_Cmat_Tmat(Cmat(sB),  Tmat(sA) % T(subB),   DIR='22')

 END SUBROUTINE self_contract_CTM_edge





 !!! Renormalize an edge of CTM by inserting isometries !!!
 SUBROUTINE renormalize_CTM(Cmat, Tmat, Z, ZP, W, WP, DIR)

  TYPE(ctm_corner_type),   INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), INTENT(INOUT) :: Tmat(:)
  COMPLEX(KIND=DP),        INTENT(IN)    :: Z(:,:), ZP(:,:)
  COMPLEX(KIND=DP),        INTENT(IN)    :: W(:,:), WP(:,:)
  CHARACTER(LEN=*),        INTENT(IN)    :: DIR

  INTEGER :: sA, sB
  INTEGER :: subA, subB

  !! Do not renormalize if edge is singleton
  IF(is_singleton_CTM_edge(Cmat, DIR)) RETURN

  !! Get the order of CTM sites on edge specified by dirFlag (that we're renormalizing) 
  CALL get_CTM_sites(sA=sA, sB=sB, subA=subA, subB=subB, DIR=DIR)

  !! Insert ZP * Z
  CALL mult_Cmat_by_isometry(Cmat(sA),            ZP,   '12')   
  CALL mult_Tmat_by_isometry(Tmat(sA) % T(subA),  Z,    '31')  

  !! Insert WP * W
  CALL mult_Tmat_by_isometry(Tmat(sA) % T(subA),  WP,   '22')   
  CALL mult_Tmat_by_isometry(Tmat(sA) % T(subB),  W,    '31')  

  !! Insert ZP * Z 
  CALL mult_Tmat_by_isometry(Tmat(sA) % T(subB),  ZP,   '22')  
  CALL mult_Cmat_by_isometry(Cmat(sB),            Z,    '21')  

  !! Rescale CTM tensors on the edge that we've just updated && renormalized
  !! and then rescale all CTM tensors together to symmetrize their weights
  CALL rescale_CTM_tensors(Cmat, Tmat)

 END SUBROUTINE renormalize_CTM


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (4) CALCULATE CTM ISOMETRY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Calculate CTM isometry by contracting its tensors into two halfs
 !!! -- the edge being updated is specified by DIR
 SUBROUTINE calc_CTM_isometry_twoDir(Z, ZP, Cmat, CTMRG, DIR)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: Z(:,:), ZP(:,:)
  TYPE(ctm_corner_type),         INTENT(IN)    :: Cmat(:)
  TYPE(ctm_params),              INTENT(IN)    :: CTMRG
  CHARACTER(LEN=*),              INTENT(IN)    :: DIR

  !! Transforming operators of CTM halfs
  COMPLEX(KIND=DP), ALLOCATABLE :: theta(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: QA(:,:), QB(:,:)

  !! SVD matrices for [QB--QA] SVD --> we'll get their HC copies to calculate isometries
  COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:),  VH(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: V(:,:),  UH(:,:) 
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)
  INTEGER                       :: chi

  !! (1) Calculate transforming operators from CTM halfs (A,B)
  CALL calc_CTM_half_transform_op(QA, Cmat, DIR, ctm_half='A')
  CALL calc_CTM_half_transform_op(QB, Cmat, DIR, ctm_half='B')

  !! (2A) Construct theta = QB * QA (product of transforming op), normalize it 
  theta = TENSMUL(QB, QA); CALL normalize_matrix(theta)

  !! (2B) Calc SVD  --> find UH, V matrices 
  !!      (NB. do NOT filter small values at this stage, because it leads to substantial errors!)
  chi = CTMRG % chi
  CALL compute_lapack_svd(theta, U, VH, Sigma, chi, CTMRG % eps)
  UH = TensTRANSPOSE(U,  'HC')
  V  = TensTRANSPOSE(VH, 'HC')

  !! (3) Contract V -> V * Sigma**(-1/2), UH = Sigma**(-1/2) * UH
  CALL LambdaDIV(V,   complx(SQRT(Sigma)),  '2')
  CALL LambdaDIV(UH,  complx(SQRT(Sigma)),  '1')

  !! (4) Contract Z = QA * V, ZP = UH * QB --> gives isometries
  !!     i.e. overall we've:  
  !!      Z = QA * V * Sigma**(-1/2) 
  !!     ZP = Sigma**(-1/2) * UH * QB
  ZP = TENSMUL(UH, QB)
  Z  = TENSMUL(QA, V)

 END SUBROUTINE calc_CTM_isometry_twoDir





 !!! Calculate the transforming operator from a CTM half !!!
 SUBROUTINE calc_CTM_half_transform_op(Q, Cmat, DIR, ctm_half)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: Q(:,:)
  TYPE(ctm_corner_type),         INTENT(IN)    :: Cmat(:)
  CHARACTER(LEN=*),              INTENT(IN)    :: DIR
  CHARACTER(LEN=*),              INTENT(IN)    :: ctm_half

  !! CTM half -- theta matrix
  COMPLEX(KIND=DP), ALLOCATABLE :: theta(:,:)

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: VH(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)
  INTEGER                       :: chi

  !! (1) Contract one half of CTM (specified by A,B) into a single theta tensor
  !!  -- with a cut across [the edge we're renormalizing] && [the edge opposite to it]
  CALL calc_CTM_half(theta, Cmat, DIR, ctm_half)

  !! (2) Calc SVD of CTM half 
  !!     (default value of chi=-1 -- means 'no truncation' in fixed-chi mode)
  !!  -- discard small singular values (< 1.0d-08)
  chi=-1
  CALL compute_lapack_svd(theta, Q, VH, Sigma, chi, eps=2.0d0)  
  CALL filter_svd_matrices(Q, VH, Sigma, eps=1.0D-08)

  !! (3) Absorb Sigma into Q --> get transforming operator from a CTM half
  CALL LambdaMUL(Q, complx(Sigma), '2')   

  !! (4) If we're doing 'B' half, we must transpose Q, 
  !!     since theta of 'B' half is in transposed order  
  !!     to be consistent with the original algorithm in [Benedikt Bruognolo PhD Thesis] 
  IF(ctm_half .EQ. 'B') Q = TensTRANSPOSE(Q)

 END SUBROUTINE calc_CTM_half_transform_op





 !!! Calculate CTM halfs, each formed by two corner matrices !!!
 SUBROUTINE calc_CTM_half(theta, Cmat, DIR, ctm_half)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: theta(:,:)
  TYPE(ctm_corner_type),         INTENT(IN)    :: Cmat(:)
  CHARACTER(LEN=*),              INTENT(IN)    :: DIR
  CHARACTER(LEN=*),              INTENT(IN)    :: ctm_half

  !CTM sites
  INTEGER :: sA, sB, sC, sD

  !! Get the order of CTM sites as we renormalize an edge specified by DIR
  CALL get_CTM_sites(sA=sA, sB=sB, sC=sC, sD=sD, DIR=DIR)

  !! Contract a half of CTM into a single tensor 
  !! (A-half or B-half, specified by ctm_half flag -- correspond to U,D halfs in Bruognolo Thesis)
  SELECT CASE(ctm_half)
  CASE('A')
       theta =               TENSMUL(Cmat(sA) % m, Cmat(sD) % m)
  CASE('B')
       theta = TensTRANSPOSE(TENSMUL(Cmat(sC) % m, Cmat(sB) % m))
  CASE DEFAULT
       CALL invalid_flag("calc_CTM_half -- invalid ctm_half ", ctm_half)
  end SELECT

  !! Normalize theta matrix
  CALL normalize_matrix(theta)

 END SUBROUTINE calc_CTM_half



 !!! Renormalize CTM supercorners on edge specified by DIR -- compress the augmented chi !!!
 SUBROUTINE renormalize_CTM_supercorners(Cmat, CTMRG, DIR)

  TYPE(ctm_corner_type),   INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_params),        INTENT(IN)    :: CTMRG
  CHARACTER(LEN=*),        INTENT(IN)    :: DIR

  !! Isometries for renormalization
  COMPLEX(KIND=DP),  ALLOCATABLE  :: W(:,:), WP(:,:)

  !! CTM sites
  INTEGER :: sA,   sB
  INTEGER :: subA, subB

  !! Do not renormalize if edge is singleton
  IF(is_singleton_CTM_edge(Cmat, DIR)) RETURN

  !! Get the order of CTM sites on edge specified by dirFlag (that we're renormalizing) 
  CALL get_CTM_sites(sA=sA, sB=sB, subA=subA, subB=subB, DIR=DIR)

  !! Calc isometry
  CALL calc_CTM_isometry_twoDir(W, WP, Cmat, CTMRG, DIR)

  !! Insert WP * W
  CALL mult_Cmat_by_isometry(Cmat(sA),  WP, '12') 
  CALL mult_Cmat_by_isometry(Cmat(sB),  W,  '21') 

 END SUBROUTINE renormalize_CTM_supercorners




 !!! Calculate a product of Cmat matrices, with cut on the edge given by DIR !!!
 SUBROUTINE calc_CTM_corner_product(theta, Cmat, DIR)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: theta(:,:)
  TYPE(ctm_corner_type),         INTENT(IN)    :: Cmat(:)
  CHARACTER(LEN=*),              INTENT(IN)    :: DIR

  !! CTM sites
  INTEGER :: sA, sB, sC, sD

  !! Get the order of CTM sites on an edge specified by DIR (that we're renormalizing) 
  CALL get_CTM_sites(sA=sA, sB=sB, sC=sC, sD=sD, DIR=DIR)

  !! Compute product of Cmat in reverse cyclic order
  theta = TENSMUL(Cmat(sA) % m, TENSMUL(Cmat(sD) % m, TENSMUL(Cmat(sC) % m, Cmat(sB) % m)))

  !! Transpose if we're doing 'N' or 'E' edges
  IF((DIR .EQ. 'N') .OR. (DIR .EQ. 'E')) theta = TensTRANSPOSE(theta)

  !! Normalize theta matrix 
  CALL normalize_matrix(theta)

 END SUBROUTINE calc_CTM_corner_product




 !!! Calculate CTM isometry using superposition of two corners of the edge we've just updated !!!
 SUBROUTINE calc_CTM_isometry_oneDir(Z, ZP, Cmat, CTMRG, DIR)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: Z(:,:), ZP(:,:)
  TYPE(ctm_corner_type),         INTENT(IN)    :: Cmat(:)
  TYPE(ctm_params),              INTENT(IN)    :: CTMRG
  CHARACTER(LEN=*),              INTENT(IN)    :: DIR

  !! Corner matrices && SVD bond dim
  COMPLEX(KIND=DP), ALLOCATABLE :: M1(:,:), M2(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: theta(:,:)
  INTEGER :: chi

  !! CTM sites
  INTEGER :: sA, sB

  !! Get the order of CTM sites on an edge specified by DIR (that we're renormalizing) 
  CALL get_CTM_sites(sA=sA, sB=sB, DIR=DIR)

  !! (1) Get corner matrices
  CALL copyTens(M1,  Cmat(sA) % m)
  M2 = TensTRANSPOSE(Cmat(sB) % m)

  !! (2) Compute superposition of M1, M2 corner matrices
  theta = TENSMUL(M1, TensTRANSPOSE(M1, 'HC')) + CONJG(TENSMUL(M2, TensTRANSPOSE(M2, 'HC')))

  !! --- Normalize theta matrix
  CALL normalize_matrix(theta)

  !! (3) Compute SVD -- set [Z isometry] = U
  chi = CTMRG % chi 
  CALL compute_lapack_svd(theta, Z, chi, CTMRG % eps)

  !! (4) Output isometry
  ZP = TensTRANSPOSE(Z, 'HC')

 END SUBROUTINE calc_CTM_isometry_oneDir

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (5) CHECK CTM CONVERGENCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!! Check convergence of CTM singular value spectra !!!
  SUBROUTINE check_CTM_convergence(is_converged, svals1, svals2, svals3, svals4, Cmat, CTMRG, datafile, iter)

  LOGICAL,                       INTENT(INOUT) :: is_converged
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: svals1(:), svals2(:)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: svals3(:), svals4(:)
  TYPE(ctm_corner_type),         INTENT(IN)    :: Cmat(:)
  TYPE(ctm_params),              INTENT(IN)    :: CTMRG
  INTEGER,                       INTENT(INOUT) :: datafile
  INTEGER,                       INTENT(IN)    :: iter

  LOGICAL :: sval_converged(4)

  !! Initialize all to FALSE
  sval_converged(:) = .FALSE.

  !! Check convergence of svals
  CALL check_conv_single_CTM(sval_converged(1), svals1, Cmat, CTMRG, datafile, iter, DIR='S')
  CALL check_conv_single_CTM(sval_converged(2), svals2, Cmat, CTMRG, datafile, iter, DIR='E')
  CALL check_conv_single_CTM(sval_converged(3), svals3, Cmat, CTMRG, datafile, iter, DIR='N')
  CALL check_conv_single_CTM(sval_converged(4), svals4, Cmat, CTMRG, datafile, iter, DIR='W')

  !! Check if all of sval_converged(:) are true
  is_converged = isTrueVec(sval_converged, 'AND') 

 END SUBROUTINE check_CTM_convergence






 !!! Check convergence of a single CTM corner matrix (or a single product of corner matrices) !!!
 SUBROUTINE check_conv_single_CTM(is_converged, svals, Cmat, CTMRG, datafile, iter, DIR)

  USE iteration_helper, ONLY: test_convergence_vec

  LOGICAL,                       INTENT(INOUT) :: is_converged
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: svals(:)
  TYPE(ctm_corner_type),         INTENT(IN)    :: Cmat(:)
  TYPE(ctm_params),              INTENT(IN)    :: CTMRG
  INTEGER,                       INTENT(INOUT) :: datafile
  INTEGER,                       INTENT(IN)    :: iter
  CHARACTER(LEN=*),              INTENT(IN)    :: DIR

  !! Storage array of new singular vals (to be passed to next iteration)
  COMPLEX(KIND=DP), ALLOCATABLE :: new_svals(:)

  !! (1) Calculate spectrum of singular values for a given corner matrix or corner product
  CALL calc_CTM_spectrum_of_svals(new_svals, Cmat, DIR, CTMRG % ctm_method)
  
  !! (2) Test convergence of singular values of Cmat
  CALL test_convergence_vec(is_converged, new_svals, svals, CTMRG % eta, datafile, iter, PRINT_ON=.TRUE.) !PRINT_ON=.FALSE.) 

  !! (3) Update the recorded spectrum of singular values (for the next convergence check)
  CALL copyTens(svals, new_svals) 

 END SUBROUTINE check_conv_single_CTM







 !!! Calculate CTM spectrum of singular values for a given corner matrix or a product of corner matrices with a specified cut !!!
 SUBROUTINE calc_CTM_spectrum_of_svals(svals, Cmat, DIR, ctm_method)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: svals(:)
  TYPE(ctm_corner_type),         INTENT(IN)    :: Cmat(:)
  CHARACTER(LEN=*),              INTENT(IN)    :: DIR
  CHARACTER(LEN=*),              INTENT(IN)    :: ctm_method

  !! Principal matrix -- a corner matrix or a product of corner matrices 
  !! (either one is specified by DIR)
  COMPLEX(KIND=DP), ALLOCATABLE :: theta(:,:)

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:), VH(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)
  INTEGER                       :: chi

  !! The site of desired corner matrix, specified by DIR
  INTEGER :: sA, i

  !! (1) Get a specified corner matrix, 
  !!     or calculate a product of corner matrices with a specified cut,
  !!     or contract half of CTM TN into a single tensor, with a cut specified by DIR
  CALL get_CTM_sites(sA=sA, DIR=DIR)
  CALL copyTens(theta, Cmat(sA) % m)

  !CALL calc_CTM_corner_product(theta, Cmat, DIR)
  !CALL calc_CTM_half(theta, Cmat, DIR, ctm_half='A')

  !! (2) Frobenius-normalize theta before calculating the spectrum
  CALL normalize_matrix(theta)

  !! (3) Compute the spectrum of singular values for a given corner matrix
  chi=-1
  CALL compute_lapack_svd(theta, U, VH, Sigma, chi, eps=2.0d0)

  !! (4) Extract new singular values -- normalize && convert Sigma to complex format
  svals = complx(Sigma/sqrt(SUM(Sigma**2)))

 END SUBROUTINE calc_CTM_spectrum_of_svals

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (6) RESCALE/NORMALIZE CTM TENSORS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! RESCALE all CTM tensors to have the same max element !!!
 SUBROUTINE rescale_CTM_tensors(Cmat, Tmat)

  TYPE(ctm_corner_type),   INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), INTENT(INOUT) :: Tmat(:)

  !! CTM sites
  INTEGER :: site, sub 
  INTEGER :: N_sites, N_sub

  !! Get num of CTM sites && subsites
  N_sites = SIZE(Cmat)
  N_sub   = SIZE(Tmat(1) % T)

  !! Rescale C's
  DO site=1,N_sites
     CALL rescaleTensMaxel(Cmat(site) % m)
  end DO

  !! Rescale T's
  DO site=1,N_sites
    DO sub=1,N_sub
       CALL rescaleTensMaxel(Tmat(site) % T(sub) % m)
    end DO
  end DO

 END SUBROUTINE rescale_CTM_tensors

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (7) AUXILIARY ROUTINES: CHECK IF CTM HAS SINGLETON DIMS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Check if we've a singleton CTM !!!
 FUNCTION is_singleton_CTM(Cmat) result(is_singleton)

  TYPE(ctm_corner_type), INTENT(IN)  :: Cmat(:)
  LOGICAL                            :: is_singleton

  INTEGER :: s, N_sites
  LOGICAL :: is_singleton_corner(4)

  !Size of Cmat
  N_sites = SIZE(Cmat)

  !Initialize to FALSE
  is_singleton_corner(:) = .FALSE.

  !Check all CTM corner matrices
  DO s=1,N_sites
     is_singleton_corner(s) = is_singleton_matrix(Cmat(s) % m)
  end DO

  !Net result -- check if all corners are singleton
  is_singleton = isTrueVec(is_singleton_corner, 'AND')
  
 END FUNCTION is_singleton_CTM




 !!! Check if a matrix has at least one singleton leg !!!
 FUNCTION is_singleton_matrix(Mat) result(is_singleton)

  COMPLEX(KIND=DP), INTENT(IN)  :: Mat(:,:)
  LOGICAL                       :: is_singleton

  !Initialize to False
  is_singleton = .FALSE.

  !If matrix has at least one singleton leg, return True
  IF((SIZE(Mat,1) .EQ. 1) .OR. (SIZE(Mat,2) .EQ. 1)) is_singleton = .TRUE.
  
 END FUNCTION is_singleton_matrix




 !!! Check if a matrix has at least one singleton leg !!!
 FUNCTION is_singleton_CTM_edge(Cmat, DIR) result(is_singleton)

  TYPE(ctm_corner_type), INTENT(IN)  :: Cmat(:)
  CHARACTER(LEN=*),      INTENT(IN)  :: DIR
  LOGICAL                            :: is_singleton

  INTEGER :: sA, sB

  !Get sites of corner matrices on the edge = DIR
  CALL get_ctm_sites(sA=sA, sB=sB, DIR=DIR)

  !Default to FALSE
  is_singleton = .FALSE.

  !Check if corner matrices have at least one singleton leg
  IF((SIZE(Cmat(sA) % m, 1) .EQ. 1) .OR. &
                & (SIZE(Cmat(sB) % m, 2) .EQ. 1)) is_singleton = .TRUE. 
  
 END FUNCTION is_singleton_CTM_edge




 !!! Print CTM final message !!!
 SUBROUTINE print_ctm_final_message(msg, Cmat, iter, use_old_ctm, svals1, svals2, svals3, svals4)

  CHARACTER(LEN=*),       INTENT(IN) :: msg
  TYPE(ctm_corner_type),  INTENT(IN) :: Cmat(:)
  INTEGER,                INTENT(IN) :: iter
  LOGICAL,                INTENT(IN) :: use_old_ctm
  COMPLEX(KIND=DP),       INTENT(IN) :: svals1(:), svals2(:), svals3(:), svals4(:)
  
  WRITE(*,*)  
  WRITE(*,*) TRIM(ADJUSTL(msg))//" N_ITER = ", iter, " bond size = ", SIZE(Cmat(1) % m, 2), " use_old_ctm = ", use_old_ctm, &
              & "smallest sval", MIN(MINVAL(ABS(svals1)), MINVAL(ABS(svals2)), MINVAL(ABS(svals3)), MINVAL(ABS(svals4)))
  WRITE(*,*)  

 END SUBROUTINE print_ctm_final_message

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (8) CTM ALGEBRA -- LOW LEVEL ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Multiply Cmat * Tmat sites !!!
 SUBROUTINE mult_Cmat_Tmat(C_site, T_site, DIR)

  TYPE(ctm_corner_type), INTENT(INOUT) :: C_site
  TYPE(block_mps),       INTENT(IN)    :: T_site
  CHARACTER(LEN=*),      INTENT(IN)    :: DIR

  !! FUSE flag
  CHARACTER(LEN=16) :: FUSE

  !! Determine correct FUSE for a given DIR of contraction
  SELECT CASE(DIR)
  CASE('22')
       FUSE = '(11)'
  CASE('31')
       FUSE = '(12)'
  CASE DEFAULT
       CALL invalid_flag("mult_Cmat_Tmat -- invalid DIR ", DIR)
  end SELECT

  !! Contract C-tensor && T-tensor
  CALL allocate_ctm_corner_site(C_site, TENSMUL(T_site % m, C_site % m, MULT=DIR, FUSE=FUSE))
  
 END SUBROUTINE mult_Cmat_Tmat




 !!! Multiply Tmat * 2D network sites !!!
 SUBROUTINE mult_Tmat_TN(T_site, TN2D_site, DIR)

  TYPE(block_mps),  INTENT(INOUT) :: T_site
  TYPE(block_mpo),  INTENT(IN)    :: TN2D_site
  CHARACTER(LEN=*), INTENT(IN)    :: DIR

  !! MULT && FUSE flags
  CHARACTER(LEN=16) :: MULT, FUSE

  !! Determine correct MULT, FUSE for a given DIR of contraction
  SELECT CASE(DIR)
  CASE('S')
       MULT = '11'; FUSE = '(42,33)'
  CASE('N')
       MULT = '21'; FUSE = '(32,43)'
  CASE('W')
       MULT = '31'; FUSE = '(12,23)'
  CASE('E')
       MULT = '41'; FUSE = '(22,13)'
  CASE DEFAULT
       CALL invalid_flag("mult_Tmat_TN -- invalid DIR ", DIR)
  end SELECT

  !! Contract T-tensor && 2D network tensor
  CALL allocate_mps_site(T_site, TENSMUL(TN2D_site % m, T_site % m, MULT, FUSE))
  
 END SUBROUTINE mult_Tmat_TN




 !!! Multiply Cmat site by isometry !!!
 SUBROUTINE mult_Cmat_by_isometry(C_site, ZZ, DIR)

  TYPE(ctm_corner_type), INTENT(INOUT) :: C_site
  COMPLEX(KIND=DP),      INTENT(IN)    :: ZZ(:,:)
  CHARACTER(LEN=*),      INTENT(IN)    :: DIR

  !! Mult Cmat site by isometry -- contract legs '21' or '12' 
  !! (1st leg: C_site, 2nd leg: ZZ)
  CALL allocate_ctm_corner_site(C_site, TENSMUL(C_site % m, ZZ, DIR))
  
 END SUBROUTINE mult_Cmat_by_isometry




 !!! Multiply Tmat site by isometry !!!
 SUBROUTINE mult_Tmat_by_isometry(T_site, ZZ, DIR)

  TYPE(block_mps),  INTENT(INOUT) :: T_site
  COMPLEX(KIND=DP), INTENT(IN)    :: ZZ(:,:)
  CHARACTER(LEN=*), INTENT(IN)    :: DIR

  !! Mult Tmat site by isometry -- contract legs '22' or '31' 
  !! (1st leg: T_site, 2nd leg: ZZ)
  CALL allocate_mps_site(T_site, TENSMUL(T_site % m, ZZ, DIR))
  
 END SUBROUTINE mult_Tmat_by_isometry




 !!! Swap transfer sites on lateral edges && swap 2D network sites !!!
 SUBROUTINE swap_CTM_block_sites(Tmat, TN2D, DIR)

  TYPE(ctm_transfer_type), INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),         INTENT(INOUT) :: TN2D(:)
  CHARACTER(LEN=*),        INTENT(IN)    :: DIR

  INTEGER :: sB, sD
  INTEGER :: subA, subB

  !! Get the order of CTM sites on the edge specified by DIR
  CALL get_CTM_sites(sB=sB, sD=sD, DIR=DIR)

  !! Swap transfer sites on lateral edges (to the one being updated)
  CALL swap_mps_sites(Tmat(sD) % T(1), Tmat(sD) % T(2))
  CALL swap_mps_sites(Tmat(sB) % T(1), Tmat(sB) % T(2))

  !! Swap sites of 2D network
  CALL swap_mpo_sites(TN2D(1), TN2D(2))

 END SUBROUTINE swap_CTM_block_sites

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE corner_transfer_matrix
