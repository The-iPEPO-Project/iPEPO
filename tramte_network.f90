MODULE tramte_network

 USE utility
 USE definitions_mps_mpo
 USE definitions_peps_pepo
 USE mps_mpo_utility
 USE mps_mpo_algebra_finite

 IMPLICIT NONE

 !! Operator sites
 TYPE op_sites_type
   COMPLEX(KIND=DP), ALLOCATABLE :: SA_Op1Op2(:,:,:,:) 
   COMPLEX(KIND=DP), ALLOCATABLE :: SA_Op1(:,:,:,:)
   COMPLEX(KIND=DP), ALLOCATABLE :: SA_Op2(:,:,:,:)
   COMPLEX(KIND=DP), ALLOCATABLE :: SB_Op2(:,:,:,:)
 END TYPE op_sites_type

CONTAINS 

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Contract TraMTE network to calculate 2P expvals !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Contract TN containing 2P correlator !!!
 SUBROUTINE contract_tramte_network(expvals, Rvec, Tmpo, Lvec, chi, eps, TN_opSites, min_sep, max_sep)

  COMPLEX(KIND=DP),    INTENT(INOUT) :: expvals(:)       !Results: expvals for each different sep
  TYPE(block_mps),     INTENT(IN)    :: Rvec(:)          !Rvec MPS
  TYPE(block_mps),     INTENT(IN)    :: Lvec(:)          !Lvec MPS
  TYPE(block_mpo),     INTENT(IN)    :: Tmpo(:,:)        !Transfer matrix MPO
  INTEGER,             INTENT(IN)    :: chi(2)           !SVD params
  REAL(KIND=DP),       INTENT(IN)    :: eps(2)           !SVD params
  TYPE(op_sites_type), INTENT(IN)    :: TN_opSites       !operator sites on TN
  INTEGER,             INTENT(IN)    :: min_sep, max_sep !sepX, sepY when computing corrs

  !! Transposed transfer MPO
  TYPE(block_mpo), ALLOCATABLE :: TmpoTransp(:,:)

  !! Right/Left boundary vectors
  TYPE(block_mps), ALLOCATABLE :: Lb(:),            LbO1(:)
  TYPE(block_mps), ALLOCATABLE :: RbO2_sep_EVEN(:), RbO2_sep_ODD(:)
  TYPE(block_mps), ALLOCATABLE :: RbO1O2_sep_00(:), RbO1O2_sep_11(:)
  
  !! Expvals && norm
  COMPLEX(KIND=DP) :: expval
  COMPLEX(KIND=DP) :: C_norm, vec_norm

  !! Indices
  INTEGER  :: sep           

  !! Find TN norm && Vector norm
  CALL contract_mps_mps(vec_norm,   Rvec,  Lvec)
  CALL contract_mps_mpo_mps(C_norm, mps1=Rvec, bimpo=Tmpo, mps2=Lvec, chi=chi, eps=eps)  !chi=(/2, 2/), eps=(/1.0D-12, 1.0D-12/))
  C_norm = C_norm / vec_norm

  !! Find transposed Tmpo
  CALL transpose_bimpo_block(TmpoTransp, Tmpo)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! (1A) Construct Left Bounds -- bare and with O1  
                      CALL calc_left_bound_OP(Lb,              Lvec, Tmpo, chi, eps)
  IF(max_sep .GT. 1)  CALL calc_left_bound_OP(LbO1,            Lvec, Tmpo, chi, eps, TN_opSites)

  !! (1B) Construct Right Bounds -- with O2 for EVEN/ODD sep (take sep=[2,3] for example)
  IF(max_sep .GT. 1)  CALL calc_right_bound_OP(RbO2_sep_EVEN,  Rvec, Tmpo, chi, eps, TN_opSites, sep=2)
  IF(max_sep .GT. 2)  CALL calc_right_bound_OP(RbO2_sep_ODD,   Rvec, Tmpo, chi, eps, TN_opSites, sep=3)

  !! (1C) Construct Right Bounds -- with O1,O2 on the same Tmpo for sep=[0,1]
                      CALL calc_right_bound_OP(RbO1O2_sep_00,  Rvec, Tmpo, chi, eps, TN_opSites, sep=0)
  IF(max_sep .GT. 0)  CALL calc_right_bound_OP(RbO1O2_sep_11,  Rvec, Tmpo, chi, eps, TN_opSites, sep=1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !! (2) Contract TN for different separations (for sepY > 2, augment LbO1 boundary as we go)
  mainloop:DO sep = min_sep, max_sep

     IF(sep .EQ. 0) THEN  

        !! (3) Contract TN for sep=0
        CALL contract_mps_mps(expval,    RbO1O2_sep_00,  Lb)

     ELSEIF(sep .EQ. 1) THEN

        !! (4) Contract TN for sep=1
        CALL contract_mps_mps(expval,    RbO1O2_sep_11,  Lb)

     ELSEIF(sep .EQ. 2) THEN

        !! (5) Contract TN for sep=EVEN: If [sep = 2] --> no need for padding, just contract the boundary vecs
        CALL contract_mps_mps(expval,    RbO2_sep_EVEN,  LbO1)

     ELSEIF(isEVEN(sep)) THEN

        !! (6) Contract TN for sep=EVEN: If [sep > 2] --> mult in extra padding Tmpo at each even sep --> augment LbO1 boundary
        CALL multiply_mps_bimpo(LbO1,    TmpoTransp,     chi,   eps)
        CALL contract_mps_mps(expval,    RbO2_sep_EVEN,  LbO1)

     ELSEIF(isODD(sep)) THEN

        !! (7) Contract TN for sep=ODD:  If [sep > 2]
        CALL contract_mps_mps(expval,    RbO2_sep_ODD,   LbO1)
     ELSE
        CALL invalid_value("contract_TN_2P_corr -- invalid sep ", sep)
     end IF

     !! Write expval
     expval                 = expval / vec_norm
     expvals(sep+1-min_sep) = expval / C_norm**(sep/2 + 1)

  END DO mainloop

 END SUBROUTINE contract_tramte_network

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Construct LHS/RHS boundaries of TraMTE network with operators absorbed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Construct LHS boundary with an operator-carrying Tmpo absorbed !!!
 SUBROUTINE calc_left_bound_OP(Lb, Lvec, Tmpo, chi, eps, TN_opSites, sep)

  TYPE(block_mps),  ALLOCATABLE, INTENT(INOUT)        :: Lb(:)
  TYPE(block_mps),               INTENT(IN)           :: Lvec(:)
  TYPE(block_mpo),               INTENT(IN)           :: Tmpo(:,:)
  INTEGER,                       INTENT(IN)           :: chi(2)               
  REAL(KIND=DP),                 INTENT(IN)           :: eps(2)  
  TYPE(op_sites_type),           INTENT(IN), OPTIONAL :: TN_opSites
  INTEGER,                       INTENT(IN), OPTIONAL :: sep

  !! The leftmost transfer matrix MPO
  TYPE(block_mpo), ALLOCATABLE :: TmpoLHS(:,:), TmpoLHSTransp(:,:)

  !! Create local copies of transfer MPO && Lvec
  CALL copy_bimpo_block(TmpoLHS, Tmpo)
  CALL copy_mps_block(Lb, Lvec)

  !! If operator is present: 
  IF(PRESENT(TN_opSites)) THEN

     !! -- insert operator into transfer mpo
     CALL insert_operator_site(TmpoLHS, TN_opSites, OP='1')

     !! -- transpose TmpoLHS for Lvec * TmpoLHS mult (must be after inserting op)
     CALL transpose_bimpo_block(TmpoLHSTransp, TmpoLHS)

     !! -- absorb [transfer mpo with operator] into Lvec
     CALL multiply_mps_bimpo(Lb,  TmpoLHSTransp,  chi,  eps)
  end IF

 END SUBROUTINE calc_left_bound_OP





 !!! Construct RHS boundary with an operator-carrying Tmpo absorbed !!!
 SUBROUTINE calc_right_bound_OP(Rb, Rvec, Tmpo, chi, eps, TN_opSites, sep)

  TYPE(block_mps),  ALLOCATABLE, INTENT(INOUT) :: Rb(:)
  TYPE(block_mps),               INTENT(IN)    :: Rvec(:)
  TYPE(block_mpo),               INTENT(IN)    :: Tmpo(:,:)
  INTEGER,                       INTENT(IN)    :: chi(2)               
  REAL(KIND=DP),                 INTENT(IN)    :: eps(2)   
  TYPE(op_sites_type),           INTENT(IN)    :: TN_opSites 
  INTEGER,                       INTENT(IN)    :: sep

  !! The rightmost transfer matrix MPO
  TYPE(block_mpo), ALLOCATABLE :: TmpoRHS(:,:)

  !! (1) Create local copies of transfer MPO && Rvec
  CALL copy_bimpo_block(TmpoRHS, Tmpo)
  CALL copy_mps_block(Rb, Rvec)

  !! (2) Insert operator(s) into TmpoRHS
  IF    (sep .EQ. 0)  THEN
                      CALL insert_operator_site(TmpoRHS, TN_opSites, OP='1-2-EVEN')
  ELSEIF(sep .EQ. 1)  THEN
                      CALL insert_operator_site(TmpoRHS, TN_opSites, OP='1-2-ODD')
  ELSEIF(isEVEN(sep)) THEN
                      CALL insert_operator_site(TmpoRHS, TN_opSites, OP='2-EVEN')
  ELSE
                      CALL insert_operator_site(TmpoRHS, TN_opSites, OP='2-ODD')
  end IF

  !! (3) Absorb [Tmpo with operator] into RHS boundary
  CALL multiply_mps_bimpo(Rb, TmpoRHS, chi, eps)

 END SUBROUTINE calc_right_bound_OP







 !!! Insert operator sites into transfer MPO !!!
 SUBROUTINE insert_operator_site(Tmpo, TN_opSites, OP)

  TYPE(block_mpo),     INTENT(INOUT) :: Tmpo(:,:)
  TYPE(op_sites_type), INTENT(IN)    :: TN_opSites
  CHARACTER(LEN=*),    INTENT(IN)    :: OP

  INTEGER :: N_sites

  N_sites = SIZE(Tmpo,2)

  SELECT CASE(OP)
  CASE('1-2-EVEN')
       CALL allocate_mpo_site(Tmpo(1, 2),        TN_OpSites % SA_Op1)
       CALL allocate_mpo_site(Tmpo(1, N_sites),  TN_OpSites % SA_Op2)
  CASE('1-2-ODD')
       CALL allocate_mpo_site(Tmpo(1, 2),        TN_OpSites % SA_Op1)
       CALL allocate_mpo_site(Tmpo(2, N_sites),  TN_OpSites % SB_Op2)
  CASE('1')
       CALL allocate_mpo_site(Tmpo(1, 2),        TN_OpSites % SA_Op1)
  CASE('2-EVEN')
       CALL allocate_mpo_site(Tmpo(1, N_sites),  TN_OpSites % SA_Op2)
  CASE('2-ODD')
       CALL allocate_mpo_site(Tmpo(2, N_sites),  TN_OpSites % SB_Op2)
  CASE DEFAULT
       CALL invalid_flag("insert_operator_site -- invalid OP ", OP)
  end SELECT

 END SUBROUTINE insert_operator_site


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!















 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Construct transfer matrix  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !!! Construct transfer matrix MPO of Trotter network !!!
 SUBROUTINE construct_transfer_matrix(Tmpo, Mps, TN2D, N_sites)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: Tmpo(:,:)
  TYPE(block_mps),              INTENT(IN)    :: Mps(:)
  TYPE(block_mpo),              INTENT(IN)    :: TN2D(:)
  INTEGER,                      INTENT(IN)    :: N_sites      

  !! Temp tensors
  COMPLEX(KIND=DP), ALLOCATABLE :: AA(:,:,:,:),   BB(:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: mpsA(:,:,:,:), mpsB(:,:,:,:)

  !! Indices
  INTEGER :: site
     
  !! Allocate empty boundary mps
  CALL allocate_empty_bimpo_block(Tmpo, N_sites)

  !! MPS sites
  mpsA = TensROTATE(Mps(1) % m,  'CCW -PI/2', '2')
  mpsB = TensROTATE(Mps(2) % m,  'CCW -PI/2', '2')

  !! MPO propagator sites
  AA   = TensROTATE(TN2D(1) % m, 'CCW -PI/2')
  BB   = TensROTATE(TN2D(2) % m, 'CCW -PI/2')


  !! (1) Write MPS sites to Tmpo
  CALL allocate_mpo_site(Tmpo(1,1),  mpsA(:,:,1:1,:))
  CALL allocate_mpo_site(Tmpo(2,1),  mpsB(:,:,1:1,:))

  !! (2) Bulk sites
  DO site = 2, N_sites-1

     IF(isEVEN(site)) THEN
        CALL allocate_mpo_site(Tmpo(1, site),  AA); CALL allocate_mpo_site(Tmpo(2, site),  BB)
     ELSE
        CALL allocate_mpo_site(Tmpo(1, site),  BB); CALL allocate_mpo_site(Tmpo(2, site),  AA)
     end IF
  end DO 

  !! (3) End sites
  CALL allocate_mpo_site(Tmpo(1, N_sites),  BB(:,:,:,1:1))
  CALL allocate_mpo_site(Tmpo(2, N_sites),  AA(:,:,:,1:1))

 END SUBROUTINE construct_transfer_matrix







 !!! Normalize <Lvec|Tmpo|Rvec> network !!!
 SUBROUTINE normalize_Rvec_Tmpo_Lvec_network(Rvec, Lvec, Tmpo, TN_opSites, C_norm)

  TYPE(block_mps),     ALLOCATABLE, INTENT(INOUT) :: Rvec(:), Lvec(:)
  TYPE(block_mpo),     ALLOCATABLE, INTENT(INOUT) :: Tmpo(:,:)
  TYPE(op_sites_type),              INTENT(INOUT) :: TN_opSites 
  COMPLEX(KIND=DP),    ALLOCATABLE, INTENT(INOUT) :: C_norm

  !! Norm constants && sites
  COMPLEX(KIND=DP) :: vec_norm
  INTEGER          :: N_sites

  !! Get size of Tmpo
  N_sites = SIZE(Tmpo, 2)

  !! Normalize Rvec, Lvec boundary vectors
  CALL contract_mps_mps(vec_norm, Rvec, Lvec)
  CALL mult_mps_by_value(Rvec, vec_norm**(-0.5D0))
  CALL mult_mps_by_value(Lvec, vec_norm**(-0.5D0))

  !! Restore unnormalized operator sites 
  !! (if they've been normalized before, i.e. C_norm already allocated)
  !! (C_norm is allocatable so we don't accidentally try to unnormalize before C_norm is computed for the 1st time)
  IF(.NOT. ALLOCATED(C_norm)) THEN
     ALLOCATE(C_norm)
  ELSE
     CALL mult_TN_opsites_by_value(TN_opSites, N_sites, C_norm)
  end IF

  !! Find new C_norm of <Lvec|Tmpo|Rvec>
  CALL contract_mps_mpo_mps(C_norm, mps1=Rvec, bimpo=Tmpo, mps2=Lvec, chi=(/2, 2/), eps=(/1.0D-12, 1.0D-12/))
  CALL mult_bimpo_by_value(Tmpo, C_norm**(-1))

  !! Normalize operator sites by new C_norm
  CALL mult_TN_opsites_by_value(TN_opSites, N_sites, C_norm**(-1))

 END SUBROUTINE normalize_Rvec_Tmpo_Lvec_network






 !!! Normalize TN operator sites !!!
 SUBROUTINE mult_TN_opsites_by_value(TN_opSites, N_sites, C_norm)

  TYPE(op_sites_type), INTENT(INOUT) :: TN_opSites
  INTEGER,             INTENT(IN)    :: N_sites 
  COMPLEX(KIND=DP),    INTENT(IN)    :: C_norm

  !! Normalize operator sites
  TN_OpSites % SA_Op1    = TN_OpSites % SA_Op1    * C_norm**(0.5D0/REAL(N_sites, KIND=DP))  
  TN_OpSites % SA_Op2    = TN_OpSites % SA_Op2    * C_norm**(0.5D0/REAL(N_sites, KIND=DP)) 
  TN_OpSites % SB_Op2    = TN_OpSites % SB_Op2    * C_norm**(0.5D0/REAL(N_sites, KIND=DP))
  TN_OpSites % SA_Op1Op2 = TN_OpSites % SA_Op1Op2 * C_norm**(0.5D0/REAL(N_sites, KIND=DP))

 END SUBROUTINE mult_TN_opsites_by_value

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Construct Trotter network !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !!! Create MPO representing Trotter network !!! (NOT conforming with CTM convention!)
 SUBROUTINE create_reduced_TN_from_trotter(TN2D, mpo)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT)  :: TN2D(:)
  TYPE(block_mpo),              INTENT(IN)     :: mpo(:)

  !! Allocate new 2D network mpo
  CALL allocate_empty_mpo_block(TN2D, 2)

  !! Construct sites of reduced Trotter network 
  CALL allocate_mpo_site(TN2D(1),  TENSMUL(mpo(2) % m, mpo(1) % m, MULT='21', FUSE='(33,44)'))
  CALL allocate_mpo_site(TN2D(2),  TENSMUL(mpo(1) % m, mpo(2) % m, MULT='21', FUSE='(33,44)'))

 END SUBROUTINE create_reduced_TN_from_trotter




 !!! Create operator sites of Trotter network !!!
 SUBROUTINE create_reduced_op_sites_trotter(TN_OpSites, TN2D, op1, op2, ROT)

  TYPE(op_sites_type),   INTENT(INOUT)         :: TN_OpSites
  TYPE(block_mpo),       INTENT(IN)            :: TN2D(:)
  COMPLEX(KIND=DP),      INTENT(IN)            :: op1(:,:), op2(:,:)
  CHARACTER(LEN=*),      INTENT(IN),  OPTIONAL :: ROT

  TN_OpSites % SA_Op1 = TENSMUL(TN2D(1) % m,  op1,         '21') 
  TN_OpSites % SA_Op2 = TENSMUL(TN2D(1) % m,  op2(1:1,:),  '12') 
  TN_OpSites % SB_Op2 = TENSMUL(TN2D(2) % m,  op2(1:1,:),  '12') 

  TN_OpSites % SA_Op1Op2 = TENSMUL(TENSMUL(TN2D(1) % m,  op2,  '21'),  op1,  '21')

  IF(PRESENT(ROT)) THEN
     TN_OpSites % SA_Op1    = TensROTATE(TN_OpSites % SA_Op1,    'CCW -PI/2')  
     TN_OpSites % SA_Op2    = TensROTATE(TN_OpSites % SA_Op2,    'CCW -PI/2')  
     TN_OpSites % SB_Op2    = TensROTATE(TN_OpSites % SB_Op2,    'CCW -PI/2')  
     TN_OpSites % SA_Op1Op2 = TensROTATE(TN_OpSites % SA_Op1Op2, 'CCW -PI/2')  
  end IF

 END SUBROUTINE create_reduced_op_sites_trotter

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Calculate factorized correlator !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Calculate factorized correlator !!!
 SUBROUTINE calc_factorized_corr(expval, iG, iL, op, opP)

  USE observables_imps

  COMPLEX(KIND=DP),   INTENT(INOUT)        :: expval
  TYPE(block_mps),    INTENT(IN)           :: iG(:)
  TYPE(block_lambda), INTENT(IN)           :: iL(:)
  COMPLEX(KIND=DP),   INTENT(IN)           :: op(:,:)
  COMPLEX(KIND=DP),   INTENT(IN), OPTIONAL :: opP(:,:)

  !! Constructing reduced RHO
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_onesite(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: Rvec(:), Lvec(:)
  COMPLEX(KIND=DP)              :: eval

  !! Temp storage for expval of 2nd operator (if present)
  COMPLEX(KIND=DP) :: temp

  !! Get boundary evecs && eval
  CALL compute_rho_imps_bounds(eval, Rvec, Lvec, iG, iL)

  !! Construct reduced rho for computing SS observables)
  CALL construct_rho_onesite(rho_onesite, iG, iL, Rvec, Lvec) 

  !! Calc expval
  CALL compute_obs_1P(expval, rho_onesite, op)
  expval = expval / eval

  !! Calc expval of the 2nd operator, if present
  IF(PRESENT(opP)) THEN
     CALL compute_obs_1P(temp, rho_onesite, opP)
     temp    = temp / eval
     expval  = expval * temp
  ELSE
     expval  = expval**2
  end IF

 END SUBROUTINE calc_factorized_corr

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE tramte_network
