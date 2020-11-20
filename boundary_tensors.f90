MODULE boundary_tensors

 USE utility
 USE definitions_mps_mpo
 USE mps_mpo_utility
 USE iteration_helper

 IMPLICIT NONE

 !! A module-level copy to be called by Arnoldi iterator (Fortran = pass by value = v.slow)
 TYPE(transfer_mpo), ALLOCATABLE :: Tmpo_arnoldi(:)

 interface mult_bvec_tmpo
    module procedure mult_bvec_tmpo_3D
    module procedure mult_bvec_tmpo_2D
 end interface mult_bvec_tmpo

 interface mult_bvec_bvec
    module procedure mult_bvec_bvec_3D
    module procedure mult_bvec_bvec_2D
 end interface mult_bvec_bvec

 interface mult_bvec_tmpo_bvec
    module procedure mult_bvec_tmpo_bvec_3D
    module procedure mult_bvec_tmpo_bvec_2D
 end interface mult_bvec_tmpo_bvec

 private mult_bvec_tmpo_3D,          mult_bvec_tmpo_2D
 private mult_bvec_bvec_3D,          mult_bvec_bvec_2D
 private mult_bvec_tmpo_bvec_3D,     mult_bvec_tmpo_bvec_2D
 private Tmpo_arnoldi

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Contraction routines for boundary vectors && transfer MPOs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Multiply bvec * Tmpo (3D version: Tmpo = MPS-MPO-MPS, bvec has 3 legs) !!!
 SUBROUTINE mult_bvec_tmpo_3D(bvec, Tmpo, MULT)

  COMPLEX(KIND=DP),   ALLOCATABLE, INTENT(INOUT) :: bvec(:,:,:)
  TYPE(transfer_mpo),              INTENT(IN)    :: Tmpo(:)
  CHARACTER(LEN=*),                INTENT(IN)    :: MULT

  !! Temp tensors
  COMPLEX(KIND=DP), ALLOCATABLE :: MID(:,:,:), UP(:,:), DN(:,:)

  !! Indices && dims
  INTEGER :: site, N_sites

  !! Get size of transfer MPO
  N_sites = SIZE(Tmpo)

  !! Contract bvec * Tmpo
  SELECT CASE(MULT)
  CASE('L')

       !! Contract from LHS
       DO site=1,N_sites

          MID = TENSMUL(Tmpo(site) % M,  bvec,  MULT='31',  FUSE='(12,23)')
          UP  = RESHAPE_2D(Tmpo(site) % U, '12,3')
          DN  = RESHAPE_2D(Tmpo(site) % D, '12,3')

          CALL copyTens(bvec, TENSMUL(TENSMUL(MID, UP, '31'), DN, '21'))
       end DO

  CASE('R')

       !! Contract from RHS  
       DO site=N_sites,1,-1

          MID = TENSMUL(Tmpo(site) % M,  bvec,  MULT='41',  FUSE='(12,23)')
          UP  = RESHAPE_2D(Tmpo(site) % U, '2,13')
          DN  = RESHAPE_2D(Tmpo(site) % D, '2,13')

          CALL copyTens(bvec, TENSMUL(TENSMUL(MID, UP, '32'), DN, '22'))
       end DO

  CASE DEFAULT
       CALL invalid_flag("mult_bvec_tmpo_3D -- invalid MULT ", MULT)
  end SELECT

 END SUBROUTINE mult_bvec_tmpo_3D






 !!! Multiply bvec * Tmpo (2D version: Tmpo = MPS-MPS, bvec has 2 legs) !!!
 SUBROUTINE mult_bvec_tmpo_2D(bvec, Tmpo, MULT)

  COMPLEX(KIND=DP),   ALLOCATABLE, INTENT(INOUT) :: bvec(:,:)
  TYPE(transfer_mpo),              INTENT(IN)    :: Tmpo(:)
  CHARACTER(LEN=*),                INTENT(IN)    :: MULT

  !! Temp tensors
  COMPLEX(KIND=DP), ALLOCATABLE :: UP(:,:), DN(:,:)

  !! Indices && dims
  INTEGER :: site, N_sites

  !! Get size of transfer MPO
  N_sites = SIZE(Tmpo)

  !! Contract bvec * Tmpo
  SELECT CASE(MULT)
  CASE('L')

       !! Contract from LHS
       DO site=1,N_sites

          UP  = TENSMUL(Tmpo(site) % U,  bvec, MULT='22',  FUSE='(11)')
          DN  = RESHAPE_2D(Tmpo(site) % D, '12,3')
          CALL copyTens(bvec, TENSMUL(UP, DN, '11'))
       end DO

  CASE('R')

       !! Contract from RHS  
       DO site=N_sites,1,-1

          UP  = TENSMUL(Tmpo(site) % U,  bvec, MULT='32',  FUSE='(11)')
          DN  = RESHAPE_2D(Tmpo(site) % D, '2,13')
          CALL copyTens(bvec, TENSMUL(DN, UP, '22'))
       end DO

  CASE DEFAULT
       CALL invalid_flag("mult_bvec_tmpo_2D -- invalid MULT ", MULT)
  end SELECT

 END SUBROUTINE mult_bvec_tmpo_2D





 !!! Contract bvec * bvec (3D version: bvec has 3 legs) !!!
 SUBROUTINE mult_bvec_bvec_3D(expval, bvecAA, bvecBB)

  COMPLEX(KIND=DP), INTENT(INOUT) :: expval
  COMPLEX(KIND=DP), INTENT(IN)    :: bvecAA(:,:,:), bvecBB(:,:,:)

  !! Contract bvecs
  expval = SUM(RESHAPE_1D(bvecAA) * RESHAPE_1D(bvecBB))

 END SUBROUTINE mult_bvec_bvec_3D



 !!! Contract bvec * bvec (2D version: bvec has 2 legs) !!!
 SUBROUTINE mult_bvec_bvec_2D(expval, bvecAA, bvecBB)

  COMPLEX(KIND=DP), INTENT(INOUT) :: expval
  COMPLEX(KIND=DP), INTENT(IN)    :: bvecAA(:,:), bvecBB(:,:)

  !! Contract bvecs
  expval = SUM(RESHAPE_1D(bvecAA) * RESHAPE_1D(bvecBB))

 END SUBROUTINE mult_bvec_bvec_2D





 !!! Contract bvec * Tmpo * bvec (3D version: Tmpo = MPS-MPO-MPS, bvec has 3 legs) !!!
 SUBROUTINE mult_bvec_tmpo_bvec_3D(expval, Rvec, Tmpo, Lvec)

  COMPLEX(KIND=DP),    INTENT(INOUT) :: expval
  COMPLEX(KIND=DP),    INTENT(IN)    :: Rvec(:,:,:), Lvec(:,:,:)
  TYPE(transfer_mpo),  INTENT(IN)    :: Tmpo(:)

  COMPLEX(KIND=DP), ALLOCATABLE :: RvecO(:,:,:)

  !! Calculate RvecO = Rvec * Tmpo
  CALL copyTens(RvecO, Rvec)
  CALL mult_bvec_tmpo(RvecO, Tmpo, MULT='R')

  !! Contract RvecO * Lvec
  expval = SUM(RESHAPE_1D(RvecO) * RESHAPE_1D(Lvec))

 END SUBROUTINE mult_bvec_tmpo_bvec_3D





 !!! Contract bvec * Tmpo * bvec (2D version: Tmpo = MPS-MPO-MPS, bvec has 2 legs) !!!
 SUBROUTINE mult_bvec_tmpo_bvec_2D(expval, Rvec, Tmpo, Lvec)

  COMPLEX(KIND=DP),    INTENT(INOUT) :: expval
  COMPLEX(KIND=DP),    INTENT(IN)    :: Rvec(:,:), Lvec(:,:)
  TYPE(transfer_mpo),  INTENT(IN)    :: Tmpo(:)

  COMPLEX(KIND=DP), ALLOCATABLE :: RvecO(:,:)

  !! Calculate RvecO = Rvec * Tmpo
  CALL copyTens(RvecO, Rvec)
  CALL mult_bvec_tmpo(RvecO, Tmpo, MULT='R')

  !! Contract RvecO * Lvec
  expval = SUM(RESHAPE_1D(RvecO) * RESHAPE_1D(Lvec))

 END SUBROUTINE mult_bvec_tmpo_bvec_2D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Calculate boundary tensors !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Calculate orthodecomposed boundary matrices !!!
 SUBROUTINE calc_orthodecomposed_boundaries(X, Y, iGamma, iLambda, SITEORD)

  COMPLEX(KIND=DP),   ALLOCATABLE, INTENT(OUT) :: X(:,:), Y(:,:)
  TYPE(block_mps),                 INTENT(IN)  :: iGamma(:)
  TYPE(block_lambda),              INTENT(IN)  :: iLambda(:)
  CHARACTER(LEN=*),                INTENT(IN)  :: SITEORD

  !! Local iMPS objects
  TYPE(block_mps),    ALLOCATABLE :: iGammaP(:)  !HC of iGamma
  TYPE(block_lambda), ALLOCATABLE :: iLambdaP(:) !HC of iLambda

  !! Boundary vectors
  COMPLEX(KIND=DP)              :: eval
  COMPLEX(KIND=DP), ALLOCATABLE :: Rvec(:,:), Lvec(:,:)

  !! Create HC copy of iMPS
  CALL copy_imps(iGammaP, iLambdaP, iGamma, iLambda, 'HC')
  
  !! Find boundary vecs
  CALL calc_boundary_vecs(eval, Rvec, Lvec, iGamma, iLambda, iGammaP, iLambdaP, SITEORD, 'CANONICAL') 

  !! Perform Hermitian purification of Rvec, Lvec
  CALL purify_hermitian_tens(Rvec); CALL purify_positive_tens(Rvec)
  CALL purify_hermitian_tens(Lvec); CALL purify_positive_tens(Lvec)

  !! RHS boundary, Rvec = X * X^HC
  CALL orthodecompose_matrix(TensTRANSPOSE(Rvec), VR=X)

  !! LHS boundary, Lvec^T = Y*HC * Y
  CALL orthodecompose_matrix(Lvec, VR=Y)
  Y = TensTRANSPOSE(Y, 'HC')

 END SUBROUTINE calc_orthodecomposed_boundaries




 !!! Construct Symmetric/Non-symmetric boundary vecs !!!
 SUBROUTINE calc_boundary_vecs(eval, Rvec, Lvec, iG, iL, iGP, iLP, SITEORD, calc_canonical_bounds)

  COMPLEX(KIND=DP),                INTENT(OUT)           :: eval 
  COMPLEX(KIND=DP),   ALLOCATABLE, INTENT(OUT), OPTIONAL :: Rvec(:,:), Lvec(:,:)
  TYPE(block_mps),                 INTENT(IN)            :: iG(:), iGP(:)
  TYPE(block_lambda),              INTENT(IN)            :: iL(:), iLP(:)
  CHARACTER(LEN=*),                INTENT(IN)            :: SITEORD
  CHARACTER(LEN=*),                INTENT(IN),  OPTIONAL :: calc_canonical_bounds !! Decides whether to compute 
                                                                                  !! canonical bounds or symmetric bounds
  !! Transfer matrix, bvec normalization constant
  TYPE(transfer_mpo), ALLOCATABLE :: Tmpo(:) 
  COMPLEX(KIND=DP)                :: C_norm

  !! Which method to use
  LOGICAL, PARAMETER :: use_arnoldi = .TRUE.

  IF(PRESENT(calc_canonical_bounds)) THEN

     !! Compute Rvec of RIGHT-normalized TN = [GA * LA * GB * LB]
     CALL construct_tmpo_from_imps_imps(Tmpo, iG, iL, iGP, iLP, SITEORD, 'R')
     IF(use_arnoldi) THEN
        CALL compute_bvec_arnoldi(eval, Rvec, Tmpo, 'R') 
     ELSE
        CALL compute_bvecs_power(eval=eval, Rvec=Rvec, Tmpo=Tmpo)
     end IF

     !! Compute Lvec of LEFT-normalized TN = [LB * GA * LA * GB]
     CALL construct_tmpo_from_imps_imps(Tmpo, iG, iL, iGP, iLP, SITEORD, 'L')
     IF(use_arnoldi) THEN
        CALL compute_bvec_arnoldi(eval, Lvec, Tmpo, 'L')
     ELSE
        CALL compute_bvecs_power(eval=eval, Lvec=Lvec, Tmpo=Tmpo)
     end IF

  ELSE
     !! Construct symmetrized transfer matrix [Sqrt(LB) * GA * LA * GB * Sqrt(LB)]
     CALL construct_tmpo_from_imps_imps(Tmpo, iG, iL, iGP, iLP, SITEORD) 

     !! Find boundary evecs && eval of Tmpo
     IF(use_arnoldi) THEN
        CALL compute_bvec_arnoldi(eval, Rvec, Tmpo, 'R')
        CALL compute_bvec_arnoldi(eval, Lvec, Tmpo, 'L')
     ELSE
        CALL compute_bvecs_power(eval, Rvec, Lvec, Tmpo)
     end IF
  end IF

 END SUBROUTINE calc_boundary_vecs




 !!! Construct transfer MPO from iMPS-iMPS inner product TN !!!
 SUBROUTINE construct_tmpo_from_imps_imps(Tmpo, iG, iL, iGP, iLP, SITEORD, DIR)

  TYPE(transfer_mpo), ALLOCATABLE, INTENT(INOUT)        :: Tmpo(:) 
  TYPE(block_mps),                 INTENT(IN)           :: iG(:), iGP(:)
  TYPE(block_lambda),              INTENT(IN)           :: iL(:), iLP(:)
  CHARACTER(LEN=*),                INTENT(IN)           :: SITEORD
  CHARACTER(LEN=*),                INTENT(IN), OPTIONAL :: DIR

  !! Local iMPS copy 
  TYPE(block_mps), ALLOCATABLE :: iMps(:), iMpsP(:)
  INTEGER                      :: sA, sB

  !! Determine iMps sites
  CALL get_imps_sites(sA, sB, SITEORD)

  !! Absorb iMPS lambdas
  IF(PRESENT(DIR)) THEN
     CALL absorb_imps_lambdas(iMps,  iG,  iL,  SITEORD, DIR)
     CALL absorb_imps_lambdas(iMpsP, iGP, iLP, SITEORD, DIR)
  ELSE
     CALL absorb_imps_lambdas(iMps,  iG,  iL,  SITEORD)
     CALL absorb_imps_lambdas(iMpsP, iGP, iLP, SITEORD)
  end IF

  !! Allocate empty Tmpo
  CALL allocate_empty_transfer_mpo(Tmpo, 2)

  !! Construct transfer mpo site-1 (= imps-imps siteA)
  CALL copyTens(Tmpo(1) % U,  iMps(sA)  % m)
  CALL copyTens(Tmpo(1) % D,  iMpsP(sA) % m)

  !! Construct transfer mpo site-2 (= imps-imps siteB)
  CALL copyTens(Tmpo(2) % U,  iMps(sB)  % m)
  CALL copyTens(Tmpo(2) % D,  iMpsP(sB) % m)

 END SUBROUTINE construct_tmpo_from_imps_imps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ARNOLDI METHOD routines for calculating the dominant boundary vecs (= evecs of Tmpo) !!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute boundary vecs using arnoldi method !!!
 SUBROUTINE compute_bvec_arnoldi(eval, bvec, Tmpo, DIR)

  USE arnoldi_wrappers, ONLY: mpo_arnoldi_largest_eig

  COMPLEX(KIND=DP),              INTENT(OUT) :: eval
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT) :: bvec(:,:)
  TYPE(transfer_mpo),            INTENT(IN)  :: Tmpo(:)
  CHARACTER(LEN=*),              INTENT(IN)  :: DIR

  !! Vectorized forms of bvecs
  COMPLEX(KIND=DP), ALLOCATABLE :: bvec_arnoldi(:)
  COMPLEX(KIND=DP)              :: EvNorm
  INTEGER                       :: dimB(2)

  !! Determine the shapes of bvecs (based on Tmpo shape)
  dimB = bvec_dims(Tmpo, DIR)

  !! (1) Initialize boundary vector, reshape into initial arnoldi evec
  CALL initialize_bvec(bvec, Tmpo, DIR)
  bvec_arnoldi = RESHAPE_1D(bvec)

  !! (2) Create a module-level copy solely for the purpose of Arnoldi iteration
  CALL copy_transfer_mpo(Tmpo_arnoldi, Tmpo) 

  !! (3) Find boundary vector using Arnoldi algorithm
  CALL mpo_arnoldi_largest_eig(app_tmpo_to_bvec, eval, bvec_arnoldi, DIR)

  !! (4) Convert the vector-form bvec into the matrix form bvec
  bvec = RESHAPE_2D(bvec_arnoldi, dimB)

 END SUBROUTINE compute_bvec_arnoldi




 !!! Apply-Op-To-Vec (Linear Operator) routine for Arnoldi !!!
 SUBROUTINE app_tmpo_to_bvec(bvecIn, bvecOut, DIR)

  USE iteration_helper, ONLY: bvec_dims

  COMPLEX(KIND=DP), INTENT(IN)  :: bvecIn(:)
  COMPLEX(KIND=DP), INTENT(OUT) :: bvecOut(:)
  CHARACTER(LEN=*), INTENT(IN)  :: DIR

  !! Local (unpacked) copy of bvec
  COMPLEX(KIND=DP), ALLOCATABLE :: bvec(:,:)
  INTEGER                       :: dimB(2)

  !! Get bvec dims
  dimB = bvec_dims(Tmpo_arnoldi, DIR)

  !! Unpack bvec_In: convert a vector-form bvec_In into a matrix-form bvec
  bvec = RESHAPE_2D(bvecIn, dimB)

  !! Mult tmpo * bvec
  CALL mult_bvec_tmpo(bvec, Tmpo_arnoldi, DIR)
 
  !! Reshape multiplied bvec into a vector again
  bvecOut = RESHAPE_1D(bvec)

 END SUBROUTINE app_tmpo_to_bvec

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! POWER METHOD routines for calculating dominant boundary vecs (= evecs of Tmpo) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute boundary vecs using power method !!!
 SUBROUTINE compute_bvecs_power(eval, Rvec, Lvec, Tmpo)

  USE datafile_utility, ONLY: setupConvDatafile

  COMPLEX(KIND=DP),              INTENT(OUT),   OPTIONAL :: eval
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: Rvec(:,:), Lvec(:,:)
  TYPE(transfer_mpo),            INTENT(IN)              :: Tmpo(:) 

  !! (1) Compute Rvec
  IF(PRESENT(Rvec)) THEN
     CALL initialize_bvec(Rvec, Tmpo, 'R')
     CALL normalize_bvec(Rvec)
     CALL compute_single_bvec_power(Rvec, Tmpo, 'R')
  end IF

  !! (2) Compute Lvec
  IF(PRESENT(Lvec)) THEN 
     CALL initialize_bvec(Lvec, Tmpo, 'L')
     CALL normalize_bvec(Lvec)
     CALL compute_single_bvec_power(Lvec, Tmpo, 'L')
  end IF

  !! (3) Compute eval
  IF(PRESENT(eval)) THEN
     IF(PRESENT(Rvec) .AND.        PRESENT(Lvec))  eval = bvecEval(Rvec, Tmpo,       Lvec,  'R')
     IF(PRESENT(Rvec) .AND. (.NOT. PRESENT(Lvec))) eval = bvecEval(Rvec, Tmpo, CONJG(Rvec), 'R')
     IF(PRESENT(Lvec) .AND. (.NOT. PRESENT(Rvec))) eval = bvecEval(Lvec, Tmpo, CONJG(Lvec), 'L')
  end IF

 END SUBROUTINE compute_bvecs_power





 !!! Compute a single boundary eigenvector of Tmpo !!!
 SUBROUTINE compute_single_bvec_power(bvec, Tmpo, DIR)

  USE datafile_utility, ONLY: setupConvDatafile

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: bvec(:,:)
  TYPE(transfer_mpo),            INTENT(IN)    :: Tmpo(:) 
  CHARACTER(LEN=*),              INTENT(IN)    :: DIR     

  !! Local vars -- old evec/eval, new eval, cost func
  COMPLEX(KIND=DP), ALLOCATABLE :: bvecOld(:,:)
  COMPLEX(KIND=DP)              :: eval, old_eval
  COMPLEX(KIND=DP)              :: costF

  !! Power method settings & params
  REAL(KIND=DP), PARAMETER :: eta = 1.0d-12 
  INTEGER,       PARAMETER :: max_iters = 5000
  INTEGER                  :: datafile, iter
  LOGICAL                  :: is_converged

  !! Setup convergence data file, init convergence tester
  CALL setupConvDatafile(datafile=datafile, descr="Bvec_Power_Method", valstr="costF", xstr="eval")
  is_converged = .FALSE.

  !! Power iteration
  powerloop:DO iter=1,max_iters

     !! Create copies of evecs
     CALL copyTens(bvecOld, bvec)

     !! Mult bvec by transfer matrix
     CALL mult_bvec_tmpo(bvec, Tmpo, DIR)

     !! Calc eval && cost function
     eval  = bvecEval(bvec, Tmpo, CONJG(bvec), DIR)
     costF = bvecCostFunc(bvec, bvecOld, C_fac=eval)

     !! Test convergence
     CALL test_cost_function(is_converged=is_converged, costF=costF, eta=eta, datafile=datafile, iter=iter, xval=eval)

     !! Normalize bvec (must do after convergence testing)
     CALL normalize_bvec(bvec)

     !! If converged: exit loop; If convergence failed after max_iterations: STOP
     IF(is_converged) THEN          
         EXIT powerloop
     ELSEIF(iter .EQ. max_iters) THEN
         WRITE(*,*) "Boundary Power Method: bvecs failed to converge after max_iters = ", iter         
         STOP
     end IF

     !! Pass eval from [iter-1] to [iter]
     old_eval = eval

  END DO powerloop

  !! Close files
  CLOSE(datafile) 

 END SUBROUTINE compute_single_bvec_power





 !!! Compute boundary vector cost function !!!
 FUNCTION bvecCostFunc(b1, b0_in, C_fac) result(costF)

  COMPLEX(KIND=DP), INTENT(IN)           :: b1(:,:), b0_in(:,:)
  COMPLEX(KIND=DP), INTENT(IN), OPTIONAL :: C_fac
  COMPLEX(KIND=DP)                       :: costF 

  !! Boudary vecs - local copies
  COMPLEX(KIND=DP), ALLOCATABLE :: b1_HC(:,:), b0_HC(:,:), b0(:,:)

  !! Bvec overlaps - local copies
  COMPLEX(KIND=DP) :: C_00, C_11, C_10, C_01

  !! Local copy of b0
  CALL copyTens(b0, b0_in)

  !! If provided, multiply b0 by a scalar constant C_fac
  IF(PRESENT(C_fac)) b0 = b0 * C_fac

  !! Create HC copies
  CALL copyTens(b0_HC, b0); b0_HC = CONJG(b0_HC)
  CALL copyTens(b1_HC, b1); b1_HC = CONJG(b1_HC)

  !! Calc bvec overlaps
  CALL mult_bvec_bvec(C_00, b0, b0_HC)
  CALL mult_bvec_bvec(C_11, b1, b1_HC)
  CALL mult_bvec_bvec(C_10, b1, b0_HC)
  CALL mult_bvec_bvec(C_01, b0, b1_HC)

  !! Calc cost function
  CostF = C_11 + C_00 - C_01 - C_10

 END FUNCTION bvecCostFunc





 !!! Compute eval of Tmpo, corresponding to bvec 
 !!! -- use aux_vec to take inner product && extract eval 
 FUNCTION bvecEval(bvec, Tmpo, aux_vec, DIR) result(eval)

  COMPLEX(KIND=DP),   INTENT(IN) :: bvec(:,:), aux_vec(:,:)
  TYPE(transfer_mpo), INTENT(IN) :: Tmpo(:)
  CHARACTER(LEN=*),   INTENT(IN) :: DIR 
  COMPLEX(KIND=DP)               :: eval 

  COMPLEX(KIND=DP), ALLOCATABLE :: bvecCpy(:,:)
  COMPLEX(KIND=DP) :: C_norm

  !! Calc eval = <aux_vec|T|bvec>
  CALL copyTens(bvecCpy, bvec)
  CALL mult_bvec_tmpo(bvecCpy, Tmpo, DIR)
  CALL mult_bvec_bvec(eval, aux_vec, bvecCpy)
    
  !! Normalize eval by C_norm = <aux_vec|bvec>
  CALL mult_bvec_bvec(C_norm, aux_vec, bvec)
  eval = eval/C_norm

 END FUNCTION bvecEval



 !!! Normalize a boundary vector !!!
 SUBROUTINE normalize_bvec(bvec)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: bvec(:,:)

  !! Local HC copy of bvec && norm const
  COMPLEX(KIND=DP), ALLOCATABLE :: bvecHC(:,:)
  COMPLEX(KIND=DP)              :: C_norm
  
  !! Create HC copy of bvec
  CALL copyTens(bvecHC, bvec)
  bvecHC = CONJG(bvecHC)

  !! Obtain norm constant 
  CALL mult_bvec_bvec(C_norm, bvec, bvecHC)

  !! Normalize bvec by C_norm
  bvec(:,:) = bvec(:,:)/sqrt(C_norm)

 END SUBROUTINE normalize_bvec

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE boundary_tensors
