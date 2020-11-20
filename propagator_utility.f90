MODULE propagator_utility

  USE utility
  USE definitions_mps_mpo
  USE simulation_parameters, ONLY: local_dim

  IMPLICIT NONE

  private
  public setup_ti_propagator_mpo,  setup_ti_propagator_dense
  public lindblad, commutator
  public sv_lindblad, sv_commutator, sv_anticommutator, sv_cross_commutator, sv_cross_anticommutator     
  public allocate_pair_mpo                        
  public mpo_add_operator, mpo_add_commutator, mpo_add_nll_term_2P,  mpo_add_nll_term_4P
  public canonicalize_propagator

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SETUP DENSE PROPAGATOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Construct TI dense propagator !!!
 SUBROUTINE setup_ti_propagator_dense(prop_dense, Hl, Hr, H_pair, dt, mps_vs_peps)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT)          :: prop_dense(:,:)
  COMPLEX(KIND=DP),              INTENT(IN)             :: Hl(:,:), Hr(:,:), H_pair(:,:)
  COMPLEX(KIND=DP),              INTENT(IN)             :: dt
  CHARACTER(LEN=*),              INTENT(IN),  OPTIONAL  :: mps_vs_peps

  COMPLEX(KIND=DP), ALLOCATABLE :: H_full(:,:), eye(:,:)
  REAL(KIND=DP)                 :: lfactor, rfactor
  INTEGER                       :: dimH

  !! Determine prefactors based on dimensionality of the problem (i.e. MPS vs PEPS)
  SELECT CASE(OptArg(mps_vs_peps, 'MPS'))
  CASE('MPS')
       lfactor = 0.50D0; rfactor = 0.50D0 
  CASE('PEPS')
       lfactor = 0.25D0; rfactor = 0.25D0 
  CASE DEFAULT
       CALL invalid_flag("setup_ti_propagator_dense -- invalid mps_vs_peps ", OptArg(mps_vs_peps, 'MPS'))
  end SELECT 

  !! Get dims, allocate H_full
  dimH = SIZE(Hl,1)

  !! Combine onsite && interaction terms
  eye = matEye(dimH)
  H_full = lfactor*TensKRON(Hl, eye) + rfactor*TensKRON(eye, Hr) + H_pair

  !! Compute exponential of H matrix
  prop_dense = matrix_exp(H_full, dt)  

 END SUBROUTINE setup_ti_propagator_dense

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SETUP MPO PROPAGATOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Construct TI MPO propagators !!!
 SUBROUTINE setup_ti_propagator_mpo(prop_mpo, super_H_site, H_pair, dt, mps_vs_peps, xi)

  TYPE(block_mpo),  ALLOCATABLE, INTENT(INOUT)          :: prop_mpo(:)
  COMPLEX(KIND=DP),              INTENT(IN)             :: super_H_site(:,:)
  TYPE(block_mpo),               INTENT(IN)             :: H_pair(:)
  COMPLEX(KIND=DP),              INTENT(IN)             :: dt
  CHARACTER(LEN=*),              INTENT(IN),  OPTIONAL  :: mps_vs_peps
  REAL(KIND=DP),                 INTENT(IN),  OPTIONAL  :: xi

  COMPLEX(KIND=DP), ALLOCATABLE :: prop_dense(:,:,:,:)
  TYPE(block_mpo),  ALLOCATABLE :: H_full(:)

  !! (1) combine onsite and two-body super H terms into H_full
  CALL combine_super_H_terms(H_full, super_H_site, super_H_site, H_pair, OptArg(mps_vs_peps, 'MPS'))

  !! (2) calculate dense propagator = exp(H_full*dt)
  CALL exponentiate_pair_mpo(prop_dense, H_full, dt)

  !! (3) canonicalize propagator: convert dense propagator into mpo propagator (with compressed bond dimension)
  CALL canonicalize_propagator(prop_mpo, prop_dense, xi)

 END SUBROUTINE setup_ti_propagator_mpo





 !!! Combine onsite and two-body terms of super H !!!
 SUBROUTINE combine_super_H_terms(H_full, Hl, Hr, H_pair, mps_vs_peps)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: H_full(:)
  COMPLEX(KIND=DP),             INTENT(IN)    :: Hl(:,:), Hr(:,:)
  TYPE(block_mpo),              INTENT(IN)    :: H_pair(:)
  CHARACTER(LEN=*),             INTENT(IN)    :: mps_vs_peps

  COMPLEX(KIND=DP), ALLOCATABLE :: eyeMat(:,:)
  REAL(KIND=DP)                 :: lfactor, rfactor
  INTEGER                       :: full_odim

  !! Determine prefactors based on dimensionality of the problem (i.e. MPS vs PEPS)
  SELECT CASE(mps_vs_peps)
  CASE('MPS')
       lfactor = 0.50D0; rfactor = 0.50D0 
  CASE('PEPS')
       lfactor = 0.25D0; rfactor = 0.25D0 
  CASE DEFAULT
       CALL invalid_flag("combine_super_H_terms -- invalid mps_vs_peps ", mps_vs_peps)
  end SELECT 

  !! Get operator bond dim
  full_odim = 2 + H_pair(1) % Edim

  !! Generate eye matrix (NB. Gell-Mann eye = same as the usual definition of eye)
  eyeMat = matEye(local_dim)
    
  !! (1) allocate H_full
  CALL allocate_pair_mpo(H_full, local_dim, full_odim)

  !! (2) copy onsite terms to H_full:
  H_full(1) % m(:,:,1,1) = lfactor*Hl
  H_full(2) % m(:,:,1,1) = eyeMat

  H_full(1) % m(:,:,1,2) = eyeMat
  H_full(2) % m(:,:,2,1) = rfactor*Hr

  !! (3) copy pair interaction terms:
  H_full(1) % m(:, :, 1:1, 3:full_odim) = H_pair(1) % m
  H_full(2) % m(:, :, 3:full_odim, 1:1) = H_pair(2) % m

 END SUBROUTINE combine_super_H_terms





 !!! Exponentiate pair MPO = exp(H_full*dt) !!!
 SUBROUTINE exponentiate_pair_mpo(prop, H_mpo, dt)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: prop(:,:,:,:)
  TYPE(block_mpo),               INTENT(IN)    :: H_mpo(:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: dt

  COMPLEX(KIND=DP), ALLOCATABLE :: H_tens(:,:,:,:), H_mat(:,:), Exp_Hdt(:,:)

  !(1) Combine H_mpo sites into a single H tensor
  H_tens = TENSMUL(H_mpo(1) % m, H_mpo(2) % m,  MULT='43', FUSE='(11,22)')

  !(2) Reshape H tensor into matrix
  H_mat  = reduce_rank(reduce_rank(H_tens))

  !(3) Compute exponential of H matrix
  Exp_Hdt = matrix_exp(H_mat, dt)
  
  !(4) Reshape matrix exponential (2D tensor) into dense propagator (4D tensor)
  prop = RESHAPE_4D(Exp_Hdt, '12,34',  (/ local_dim, local_dim /),  (/ local_dim, local_dim /))

 END SUBROUTINE exponentiate_pair_mpo

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!!!!! CROSS-TERM COMMUTATOR, COMMUTATOR, LINDBLAD (USING SIMPLE VECTORIZATION) !!!!!!!!!!!!!!!!!!!!!!

 !!! Calculate cross term in the commutator (sv -- simple vectorization) !!!
 FUNCTION sv_cross_commutator(op1, op2) result(commutator)

   COMPLEX(KIND=DP), INTENT(IN)  :: op1(:,:), op2(:,:)
   COMPLEX(KIND=DP), ALLOCATABLE :: commutator(:,:)

   COMPLEX(KIND=DP), ALLOCATABLE :: eye(:,:)

   !! Check op dims
   CALL check_sizes_equal(SIZE(op1,1), SIZE(op1,2), "sv_cross_commutator -- op1 must be a square matrix ")
   CALL check_sizes_equal(SIZE(op2,1), SIZE(op2,2), "sv_cross_commutator -- op2 must be a square matrix ")
   CALL check_sizes_equal(SIZE(op1,1), SIZE(op2,2), "sv_cross_commutator -- op1 and op2 must have equal dims ")

   !! Create eye operator
   eye = matEye(SIZE(op1,1))

   !! Calculate cross term
   commutator = &
   & TensKRON(TensKRON(op1, eye), TensKRON(op2, eye)) - TensKRON(TensKRON(eye, TRANSPOSE(op1)), TensKRON(eye, TRANSPOSE(op2)))

 END FUNCTION sv_cross_commutator



 !!! Calculate cross term in the commutator (sv -- simple vectorization) !!!
 FUNCTION sv_cross_anticommutator(op1, op2) result(commutator)

   COMPLEX(KIND=DP), INTENT(IN)  :: op1(:,:), op2(:,:)
   COMPLEX(KIND=DP), ALLOCATABLE :: commutator(:,:)

   COMPLEX(KIND=DP), ALLOCATABLE :: eye(:,:)

   !! Check op dims
   CALL check_sizes_equal(SIZE(op1,1), SIZE(op1,2), "sv_cross_anticommutator -- op1 must be a square matrix ")
   CALL check_sizes_equal(SIZE(op2,1), SIZE(op2,2), "sv_cross_anticommutator -- op2 must be a square matrix ")
   CALL check_sizes_equal(SIZE(op1,1), SIZE(op2,2), "sv_cross_anticommutator -- op1 and op2 must have equal dims ")

   !! Create eye operator
   eye = matEye(SIZE(op1,1))

   !! Calculate cross term
   commutator = &
   & TensKRON(TensKRON(op1, eye), TensKRON(op2, eye)) + TensKRON(TensKRON(eye, TRANSPOSE(op1)), TensKRON(eye, TRANSPOSE(op2)))

 END FUNCTION sv_cross_anticommutator




 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! NAME:   svcommutator (sv -- simple vectorization)
 ! VARIABLES:
 !           op - Some operator
 ! SYNOPSIS:
 ! Returns the super-operator corresponding to op - Transpose(op)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION sv_commutator(op) result(commutator)

   COMPLEX(KIND=DP), INTENT(IN)  :: op(:,:)
   COMPLEX(KIND=DP), ALLOCATABLE :: commutator(:,:)

   COMPLEX(KIND=DP), ALLOCATABLE :: eye(:,:)
   INTEGER                       :: dimO

   !! Check op dims
   CALL check_sizes_equal(SIZE(op,1), SIZE(op,2), "commutator -- op must be a square matrix ")

   !! Get dims, allocate commutator: op = [dimO, dimO] --> commutator = [dimO**2, dimO**2] 
   dimO = SIZE(op,1)
   ALLOCATE(commutator(dimO**2, dimO**2))

   !! Create eye operator
   eye = matEye(dimO)

   !! Calculate commutator
   commutator = TensKRON(op, eye) - TensKRON(eye, TRANSPOSE(op))

 END FUNCTION sv_commutator




 !!! Calculate anticommutator !!!
 FUNCTION sv_anticommutator(op) result(commutator)

   COMPLEX(KIND=DP), INTENT(IN)  :: op(:,:)
   COMPLEX(KIND=DP), ALLOCATABLE :: commutator(:,:)

   COMPLEX(KIND=DP), ALLOCATABLE :: eye(:,:)
   INTEGER                       :: dimO

   !! Check op dims
   CALL check_sizes_equal(SIZE(op,1), SIZE(op,2), "anticommutator -- op must be a square matrix ")

   !! Get dims, allocate commutator: op = [dimO, dimO] --> commutator = [dimO**2, dimO**2] 
   dimO = SIZE(op,1)
   ALLOCATE(commutator(dimO**2, dimO**2))

   !! Create eye operator
   eye = matEye(dimO)

   !! Calculate commutator
   commutator = TensKRON(op, eye) + TensKRON(eye, TRANSPOSE(op))

 END FUNCTION sv_anticommutator




 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! NAME:   svlindblad (sv -- simple vectorization)
 ! VARIABLES:
 !           X - jump operator
 ! SYNOPSIS:
 ! Produces the lindblad form for superoperator op 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION sv_lindblad(X) result(lindblad)

   COMPLEX(KIND=DP), INTENT(IN)  :: X(:,:)
   COMPLEX(KIND=DP), ALLOCATABLE :: lindblad(:,:)

   COMPLEX(KIND=DP), ALLOCATABLE :: Xp(:,:), eye(:,:)
   INTEGER                       :: dimX

   !! Check op dims
   CALL check_sizes_equal(SIZE(X,1), SIZE(X,2), "lindblad -- X must be a square matrix ")

   !! Get dims, allocate Lindbladian: X = [dimX, dimX] --> lindblad = [dimX**2, dimX**2] 
   dimX = SIZE(X,1)
   ALLOCATE(lindblad(dimX**2, dimX**2))

   !! Create eye && conjugate op
   eye = matEye(dimX)
   Xp  = TensTRANSPOSE(X, 'HC')

   !! Calculate Lindbladian
   lindblad = 2.0D0*TensKRON(X, TRANSPOSE(Xp)) - TensKRON(MATMUL(Xp, X), eye) - TensKRON(eye, TRANSPOSE(MATMUL(Xp, X)))
    
 END FUNCTION sv_lindblad


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMMUTATOR && LINDBLAD (USING GELL-MANN BASIS) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! NAME:   commutator
 ! VARIABLES:
 !           op - Some superoperator
 ! SYNOPSIS:
 ! Returns the super-operator corresponding to op - Transpose(op)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION commutator(op)

   COMPLEX(KIND=DP), INTENT(IN) :: op(:,:)
   COMPLEX(KIND=DP)             :: commutator(SIZE(op,1), SIZE(op,2))

   commutator = op - TRANSPOSE(op)

 END FUNCTION commutator




 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! NAME:   lindblad
 ! VARIABLES:
 !           X - Superoperator form of jump operator
 ! SYNOPSIS:
 ! Produces the lindblad form for superoperator op,
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 FUNCTION lindblad(X)

   COMPLEX(KIND=DP), INTENT(IN) :: X(:,:)
   COMPLEX(KIND=DP), DIMENSION(SIZE(X,1),SIZE(X,2)) :: lindblad, Xp
    
   Xp=TRANSPOSE(CONJG(X))

   lindblad = 2.0D0*MATMUL(X, TRANSPOSE(Xp)) - MATMUL(Xp,X) - TRANSPOSE(MATMUL(Xp,X))
    
 END FUNCTION lindblad

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PAIR MPOs (USING GELL-MANN BASIS) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Allocate pair mpo (two-body operator) !!!
 SUBROUTINE allocate_pair_mpo(pair_mpo, SNdim, EWdim)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: pair_mpo(:)
  INTEGER,                      INTENT(IN)    :: SNdim, EWdim

  CALL allocate_empty_mpo_block(pair_mpo, 2)

  CALL allocate_mpo_site(pair_mpo(1), SNdim, SNdim, 1,     EWdim)
  CALL allocate_mpo_site(pair_mpo(2), SNdim, SNdim, EWdim,     1)

 END SUBROUTINE allocate_pair_mpo



 !!! Add a commutator to pair_mpo !!!
 SUBROUTINE mpo_add_commutator(pair_mpo, ocnt, factor, op1, op2)

  TYPE(block_mpo),  INTENT(INOUT) :: pair_mpo(:)
  INTEGER,          INTENT(INOUT) :: ocnt
  COMPLEX(KIND=DP), INTENT(IN)    :: factor
  COMPLEX(KIND=DP), INTENT(IN)    :: op1(:,:), op2(:,:)

  INTEGER :: o 

  !Verify we have enough space
  CALL check_space_exists(ocnt+2, pair_mpo(1) % Edim, "Attempt to add too many elements to existing pair MPO")
    
  !Operators in the 1st term of commutator
  o=ocnt+1
  pair_mpo(1) % m(:,:,1,o) = sqrt(factor) * op1 
  pair_mpo(2) % m(:,:,o,1) = sqrt(factor) * op2 

  !Operators in the 2nd term of commutator. Prefactor of i means that we get requisite -1 in the product.
  o=ocnt+2
  pair_mpo(1) % m(:,:,1,o) = ii*sqrt(factor) * TRANSPOSE(op1) 
  pair_mpo(2) % m(:,:,o,1) = ii*sqrt(factor) * TRANSPOSE(op2)

  !Record where we got to.
  ocnt=o

 END SUBROUTINE mpo_add_commutator





 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! NAME:   mpo_add_operator
 ! VARIABLES:
 !           factor   - Prefactor, put sqrt of both on each
 !           op1, op2 - Operators for L&R sites
 ! SYNOPSIS:
 !           Add operator to pair_mpo
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE mpo_add_operator(pair_mpo, ocnt, factor, op1, op2)

  TYPE(block_mpo),  INTENT(INOUT) :: pair_mpo(:)
  INTEGER,          INTENT(INOUT) :: ocnt
  COMPLEX(KIND=DP), INTENT(IN)    :: factor
  COMPLEX(KIND=DP), INTENT(IN)    :: op1(:,:), op2(:,:)

  INTEGER :: o 

  !Verify we have enough space
  CALL check_space_exists(ocnt+1, pair_mpo(1) % Edim, "Attempt to add too many elements to existing pair MPO")
    
  !Add operators to MPO
  o=ocnt+1
  pair_mpo(1) % m(:,:,1,o) = sqrt(factor) * op1 
  pair_mpo(2) % m(:,:,o,1) = sqrt(factor) * op2 
    
  !Record where we got to.
  ocnt=o

 END SUBROUTINE mpo_add_operator




 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! NAME:   mpo_add_nll_term
 ! VARIABLES:
 !            factor    - Prefactor, put sqrt of both on each
 !            op1, op2 - Operators for site 1 and 2 resp., multiplied from the left
 !	      op3, op4 - Operators for site 1 (op3) and 2 (op4) multiplied from the right
 ! SYNOPSIS:
 !            Add 4-point nonlocal loss term
 ! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE mpo_add_nll_term_4P(pair_mpo, ocnt, factor, op1, op2, op3, op4)

  TYPE(block_mpo),  INTENT(INOUT) :: pair_mpo(:)
  INTEGER,          INTENT(INOUT) :: ocnt
  COMPLEX(KIND=DP), INTENT(IN)    :: factor
  COMPLEX(KIND=DP), INTENT(IN)    :: op1(:,:), op2(:,:), op3(:,:), op4(:,:)

  INTEGER :: o 

  !! Verify we have enough space
  CALL check_space_exists(ocnt+1, pair_mpo(1) % Edim, "Attempt to add too many elements to existing pair MPO")

  !! Add operators to MPO
  o=ocnt+1
  pair_mpo(1) % m(:,:,1,o) = sqrt(factor) * MATMUL(op1, TRANSPOSE(op3))
  pair_mpo(2) % m(:,:,o,1) = sqrt(factor) * MATMUL(op2, TRANSPOSE(op4))
    
  !! Record where we got to.
  ocnt=o

 END SUBROUTINE mpo_add_nll_term_4P





 !!! Add 2-point nonlocal loss term !!!
 !!! (i.e. we add operators op1, op2 acting on a pair of sites i,j to a pair_mpo) 
 SUBROUTINE mpo_add_nll_term_2P(pair_mpo, ocnt, fac, op1, op2)

   TYPE(block_mpo),   INTENT(INOUT) :: pair_mpo(:)
   INTEGER,           INTENT(INOUT) :: ocnt
   COMPLEX(KIND=DP),  INTENT(IN)    :: fac
   COMPLEX(KIND=DP),  INTENT(IN)    :: op1(:,:), op2(:,:)

   INTEGER :: o 
    
   !! Verify we have enough space
   CALL check_space_exists(ocnt+1, pair_mpo(1) % Edim, "Attempt to add too many elements to existing pair MPO")

   !! Add cross-term that consists of [op1 x op2]
   o=ocnt+1
   pair_mpo(1) % m(:,:,1,o) = sqrt(fac) * op1
   pair_mpo(2) % m(:,:,o,1) = sqrt(fac) * op2
    
   !! Record where we got to.
   ocnt = o

 END SUBROUTINE mpo_add_nll_term_2P

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CANONICALIZE MPO PROPAGATOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Canonicalize propagator, convert dense prop into mpo prop !!!
 SUBROUTINE canonicalize_propagator(prop_mpo, prop_dense, xi)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT)          :: prop_mpo(:)
  COMPLEX(KIND=DP),             INTENT(IN)             :: prop_dense(:,:,:,:)
  REAL(KIND=DP),                INTENT(IN),   OPTIONAL :: xi

  !SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: theta(:,:), U(:,:), VH(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: sigma(:) 

  !SVD params
  INTEGER, PARAMETER :: chi0 = 1    
  REAL(KIND=DP)      :: eps0

  !Indices & dims
  INTEGER            :: chi, i

  WRITE(*,'("*** Canonicalising propagator ***")')

  !Set SVD precision
  eps0 = OptArg(xi, 1.0D-15) !0.0D0)

  !Reshape into matrix (along EW axis)
  theta = RESHAPE_2D(prop_dense, '13,24')

  !Compute SVD (initialize chi first)
  chi=chi0
  CALL compute_lapack_svd(theta, U, VH, Sigma, chi, eps0)

  !Allocate prop_mpo
  CALL allocate_pair_mpo(prop_mpo, local_dim, chi)

  !Copy back SVD results
  CALL LambdaMUL(U, complx(sigma), VH)
  prop_mpo(1) % m = increase_rank(RESHAPE_3D(U,             '12,3',  (/ local_dim, local_dim /)), '3')
  prop_mpo(2) % m = increase_rank(RESHAPE_3D(TRANSPOSE(VH), '12,3',  (/ local_dim, local_dim /)), '4')

  !prop_mpo(1) % m = increase_rank(RESHAPE_3D(U,             '12,3',  (/ SIZE(prop_dense,1), SIZE(prop_dense,1) /)), '3')
  !prop_mpo(2) % m = increase_rank(RESHAPE_3D(TRANSPOSE(VH), '12,3',  (/ SIZE(prop_dense,1), SIZE(prop_dense,1) /)), '4')

  !Print sigma spectrum of the propagator
  DO i=1,SIZE(sigma)
     WRITE(*,*) "Canonicalize propagator - sigma: ", sigma(i), " at i = ", i
  end DO 

  WRITE(*,*) "Constructing new mpo propagator with bond dimension ", prop_mpo(1) % Edim

 END SUBROUTINE canonicalize_propagator

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE propagator_utility 
