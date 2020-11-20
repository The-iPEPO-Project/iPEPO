MODULE imaginary_observables_imps

 USE utility
 USE definitions_mps_mpo
 USE mps_mpo_utility
 USE boundary_tensors

 IMPLICIT NONE

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPUTING OBSERVABLES DURING/AFTER TIME EVOLUTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute boundaries of infinite MPS !!!
 SUBROUTINE compute_psi_imps_bounds(eval, Rvec, Lvec, iGamma, iLambda)

  COMPLEX(KIND=DP),              INTENT(INOUT) :: eval
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: Rvec(:,:), Lvec(:,:)
  TYPE(block_mps),               INTENT(IN)    :: iGamma(:)            !iGamma of MPS
  TYPE(block_lambda),            INTENT(IN)    :: iLambda(:)           !iLambda of MPS

  !! HC iMPS copy
  TYPE(block_mps),    ALLOCATABLE :: iGammaP(:) 
  TYPE(block_lambda), ALLOCATABLE :: iLambdaP(:) 

  !! Create HC copy of iMPS
  CALL copy_imps(iGammaP, iLambdaP, iGamma, iLambda, 'HC')

  !! Calc iMPS-iMPS bounds
  CALL calc_boundary_vecs(eval, Rvec, Lvec, iGamma, iLambda, iGammaP, iLambdaP, '12')

 END SUBROUTINE compute_psi_imps_bounds




 !!! Construct Tmat containing operator !!!
 SUBROUTINE construct_tmat_with_op(TmatO, iGamma, iLambda, op, opP, op_site)

  TYPE(transfer_mpo), ALLOCATABLE, INTENT(INOUT)        :: TmatO(:) 
  TYPE(block_mps),                 INTENT(IN)           :: iGamma(:)          !iGamma of MPS
  TYPE(block_lambda),              INTENT(IN)           :: iLambda(:)         !iLambda of MPS
  COMPLEX(KIND=DP),                INTENT(IN), OPTIONAL :: op(:,:), opP(:,:)
  INTEGER,                         INTENT(IN), OPTIONAL :: op_site
  
  !! HC iMPS and O*iMPS local copies
  TYPE(block_mps),    ALLOCATABLE :: iGammaP(:),  iGammaO(:) 
  TYPE(block_lambda), ALLOCATABLE :: iLambdaP(:), iLambdaO(:) 

  !! Create local copies
  CALL copy_imps(iGammaP, iLambdaP, iGamma, iLambda, 'HC')
  CALL copy_imps(iGammaO, iLambdaO, iGamma, iLambda)

  !! Apply Op at site=op_site
  IF(PRESENT(op) .AND. PRESENT(op_site) .AND. (.NOT. PRESENT(opP))) THEN

     CALL allocate_mps_site(iGammaO(op_site), TENSMUL(iGammaO(op_site) % m,  op,  '12'))

  ELSEIF(PRESENT(op) .AND. PRESENT(opP)) THEN

     CALL allocate_mps_site(iGammaO(1),  TENSMUL(iGammaO(1) % m,  op,  '12'))
     CALL allocate_mps_site(iGammaO(2),  TENSMUL(iGammaO(2) % m,  opP, '12'))
  end IF

  !! Construct TmatO (with operator)
  CALL construct_tmpo_from_imps_imps(TmatO, iGammaO, iLambdaO, iGammaP, iLambdaP, SITEORD='12')

 END SUBROUTINE construct_tmat_with_op




 !!! Compute 1P observable !!!
 SUBROUTINE compute_obs_1P(expval, op, iGamma, iLambda, Rvec, Lvec)

  COMPLEX(KIND=DP),   INTENT(INOUT) :: expval
  COMPLEX(KIND=DP),   INTENT(IN)    :: op(:,:)
  TYPE(block_mps),    INTENT(IN)    :: iGamma(:)            !iGamma of MPS
  TYPE(block_lambda), INTENT(IN)    :: iLambda(:)           !iLambda of MPS
  COMPLEX(KIND=DP),   INTENT(IN)    :: Rvec(:,:), Lvec(:,:) !Bounds of iMPS

  !! Local vars
  COMPLEX(KIND=DP),   ALLOCATABLE :: LvecO1(:,:)
  TYPE(transfer_mpo), ALLOCATABLE :: TmatO(:)
  COMPLEX(KIND=DP)                :: C_norm

  !! Create local copy of Lvec
  CALL copyTens(LvecO1, Lvec)

  !! (1) Construct TmatO (with Op)
  CALL construct_tmat_with_op(TmatO, iGamma, iLambda, op=op, op_site=1)

  !! (2) Contract with bounds --> get 1P EXPVAL
  CALL mult_bvec_tmpo(LvecO1, TmatO, 'L')
  CALL mult_bvec_bvec(expval, Rvec, LvecO1)

  !! (3) Normalize 1P EXPVAL
  CALL mult_bvec_bvec(C_norm, Rvec, Lvec)
  expval = expval/C_norm

 END SUBROUTINE compute_obs_1P




 !!! Compute 2P observable !!!
 SUBROUTINE compute_obs_2P(expvals, ops, iGamma, iLambda, Rvec, Lvec, min_sep, max_sep)

  COMPLEX(KIND=DP),   INTENT(INOUT) :: expvals(:)
  COMPLEX(KIND=DP),   INTENT(IN)    :: ops(:,:,:)
  TYPE(block_mps),    INTENT(IN)    :: iGamma(:)            !iGamma of MPS
  TYPE(block_lambda), INTENT(IN)    :: iLambda(:)           !iLambda of MPS
  COMPLEX(KIND=DP),   INTENT(IN)    :: Rvec(:,:), Lvec(:,:) !Bounds of iMPS
  INTEGER,            INTENT(IN)    :: min_sep, max_sep

  !! Local vars
  COMPLEX(KIND=DP),   ALLOCATABLE :: LvecO1(:,:),      LvecO1O2(:,:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: RvecO2_EVEN(:,:), RvecO2_ODD(:,:)
  TYPE(transfer_mpo), ALLOCATABLE :: Tmat(:)
  TYPE(transfer_mpo), ALLOCATABLE :: TmatO1(:),      TmatO1O2(:)
  TYPE(transfer_mpo), ALLOCATABLE :: TmatO2_EVEN(:), TmatO2_ODD(:)

  COMPLEX(KIND=DP) :: expval, C_norm !! Temp storage for expvals
  INTEGER          :: sep            !! Indices && dims

  !! Min sep greater than two is not supported
  CALL check_space_exists(min_sep, 2, "compute_obs_2P -- min_sep must be <= 2: ")

  !! Create local copies of Rvec, Lvec
  CALL copyTens(LvecO1,      Lvec);  CALL copyTens(LvecO1O2,   Lvec)
  CALL copyTens(RvecO2_EVEN, Rvec);  CALL copyTens(RvecO2_ODD, Rvec)

  !! (1) Construct bare Tmat and TmatO with Op for EVEN and ODD separations
  CALL construct_tmat_with_op(Tmat,         iGamma,  iLambda)
  CALL construct_tmat_with_op(TmatO1O2,     iGamma,  iLambda,  op=ops(1,:,:),   opP=ops(2,:,:))
  CALL construct_tmat_with_op(TmatO1,       iGamma,  iLambda,  op=ops(1,:,:),   op_site=1)
  CALL construct_tmat_with_op(TmatO2_EVEN,  iGamma,  iLambda,  op=ops(2,:,:),   op_site=1)
  CALL construct_tmat_with_op(TmatO2_ODD,   iGamma,  iLambda,  op=ops(2,:,:),   op_site=2)

  !! (2) Absorb TmatO1,TmatO2 into boundaries 
  CALL mult_bvec_tmpo(LvecO1,      TmatO1,      'L')
  CALL mult_bvec_tmpo(LvecO1O2,    TmatO1O2,    'L')
  CALL mult_bvec_tmpo(RvecO2_EVEN, TmatO2_EVEN, 'R')
  CALL mult_bvec_tmpo(RvecO2_ODD,  TmatO2_ODD,  'R')

  !! (3) Find evec norm
  CALL mult_bvec_bvec(C_norm, Rvec, Lvec)

  !! (4) Main loop: contract iMPS-O1-O2-iMPS for different SEP (for every even sep, we augment LvecO1 by a bare Tmat)
  DO sep=min_sep,max_sep

     IF(sep .EQ. 0) THEN

        CALL compute_obs_1P(expval, MATMUL(ops(1,:,:), ops(2,:,:)), iGamma, iLambda, Rvec, Lvec)

     ELSEIF(sep .EQ. 1) THEN

        CALL mult_bvec_bvec(expval, Rvec, LvecO1O2)

     ELSEIF(sep .EQ. 2) THEN

        CALL mult_bvec_bvec(expval, RvecO2_EVEN, LvecO1)

     ELSEIF(isEVEN(sep)) THEN

        CALL mult_bvec_tmpo(LvecO1,  Tmat,  'L')
        CALL mult_bvec_bvec(expval, RvecO2_EVEN, LvecO1)

     ELSEIF(isODD(sep)) THEN

        CALL mult_bvec_bvec(expval, RvecO2_ODD,  LvecO1)   
     ELSE
        WRITE(*,*) "compute_obs_2P -- invalid sep ", sep
        STOP    
     end IF   

     !! Write expval (not for sep=0 -- norm is already done in compute_obs_1P() routine)
     IF(sep .GT. 0) THEN
        expvals(sep-min_sep+1) = expval/C_norm
     ELSE
        expvals(sep-min_sep+1) = expval
     end IF
  end DO

 END SUBROUTINE compute_obs_2P

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE imaginary_observables_imps
