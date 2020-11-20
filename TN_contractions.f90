MODULE TN_contractions

 USE utility
 USE definitions_mps_mpo
 USE definitions_peps_pepo
 USE mps_mpo_utility
 USE boundary_tensors

 USE ctm_definitions
 USE corner_transfer_matrix

 USE project_ipeps,        ONLY: absorb_ipeps_lambdas

 IMPLICIT NONE

 TYPE op_sites_type
   COMPLEX(KIND=DP), ALLOCATABLE :: site1_Op1(:,:,:,:)
   COMPLEX(KIND=DP), ALLOCATABLE :: site1_Op2(:,:,:,:), site2_Op2(:,:,:,:)
   COMPLEX(KIND=DP), ALLOCATABLE :: site1_Op1_Op2(:,:,:,:) 
 END TYPE op_sites_type

CONTAINS 

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Contract 2D network to calculate 1P and 2P expvals !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Contract TN containing 1P observable !!!
 SUBROUTINE contract_TN_1P_obs(expval, Rvec, Tmpo, Lvec, TN_opSites, opFlag)

  COMPLEX(KIND=DP),    INTENT(INOUT)        :: expval         !expval = result of contraction
  COMPLEX(KIND=DP),    INTENT(IN)           :: Rvec(:,:,:)    !boundary vecs -- right vec
  COMPLEX(KIND=DP),    INTENT(IN)           :: Lvec(:,:,:)    !boundary vecs -- left vec
  TYPE(transfer_mpo),  INTENT(IN)           :: Tmpo(:)        !transfer MPO
  TYPE(op_sites_type), INTENT(IN)           :: TN_opSites     !operator sites on TN
  CHARACTER(LEN=*),    INTENT(IN), OPTIONAL :: opFlag

  !! Local LHS boundary vector with operator absorbed inside
  COMPLEX(KIND=DP), ALLOCATABLE :: Lb(:,:,:)
  COMPLEX(KIND=DP)              :: C_norm

  !! Find TN norm
  CALL mult_bvec_tmpo_bvec(C_norm,  Rvec, Tmpo, Lvec)

  !! Construct Left Bound with or without operator
  CALL calc_left_bound_OP(Lb, Lvec, Tmpo, TN_opSites, opFlag=OptArg(opFlag, 'op1'))

  !! Contract TN
  CALL mult_bvec_bvec(expval, Rvec, Lb)
  expval = expval / C_norm

 END SUBROUTINE contract_TN_1P_obs







 !!! Contract TN containing 2P correlator !!!
 SUBROUTINE contract_TN_2P_corr(expvals, Rvec, Tmpo, Lvec, TN_opSites, min_sep, max_sep)

  COMPLEX(KIND=DP),    INTENT(INOUT) :: expvals(:)            !results: expvals for each different sep
  COMPLEX(KIND=DP),    INTENT(IN)    :: Rvec(:,:,:)           !boundary vecs -- right vec
  COMPLEX(KIND=DP),    INTENT(IN)    :: Lvec(:,:,:)           !boundary vecs -- left vec
  TYPE(transfer_mpo),  INTENT(IN)    :: Tmpo(:)               !transfer MPO
  TYPE(op_sites_type), INTENT(IN)    :: TN_opSites            !operator sites on TN
  INTEGER,             INTENT(IN)    :: min_sep, max_sep      !min && max sep when computing corrs

  !! Right/Left boundary vectors
  COMPLEX(KIND=DP), ALLOCATABLE :: Lb(:,:,:),            LbO1(:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: RbO2_sep_EVEN(:,:,:), RbO2_sep_ODD(:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: RbO1O2_sep_00(:,:,:), RbO1O2_sep_11(:,:,:)
  
  !! Expvals && norm
  COMPLEX(KIND=DP) :: expval, C_norm, vec_norm
  INTEGER          :: sep            

  !! Min sep greater than two is not supported
  CALL check_space_exists(min_sep, 2, "contract_TN_2P_corr -- min_sep must be <= 2: ")

  !! Find TN norm && Vector norm
  CALL mult_bvec_bvec(vec_norm, Rvec, Lvec)
  CALL mult_bvec_tmpo_bvec(C_norm, Rvec, Tmpo, Lvec)
  C_norm = C_norm / vec_norm


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! (1A) Construct Left Bounds -- bare and with O1
  CALL calc_left_bound_OP(Lb,              Lvec, Tmpo)
  CALL calc_left_bound_OP(LbO1,            Lvec, Tmpo, TN_opSites)

  !! (1B) Construct Right Bounds -- with O2 for EVEN/ODD sep (take sep=2,3 for example)
  CALL calc_right_bound_OP(RbO2_sep_EVEN,  Rvec, Tmpo, TN_opSites, sep=2)
  CALL calc_right_bound_OP(RbO2_sep_ODD,   Rvec, Tmpo, TN_opSites, sep=3)

  !! (1C) Construct Right Bounds -- with O1,O2 on the same Tmpo for sep=0,1
  CALL calc_right_bound_OP(RbO1O2_sep_00,  Rvec, Tmpo, TN_opSites, sep=0)
  CALL calc_right_bound_OP(RbO1O2_sep_11,  Rvec, Tmpo, TN_opSites, sep=1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !! (2) Contract TN for different separations (for sep > 2, augment LbO1 boundary as we go)
  separation_loop:DO sep = min_sep, max_sep

     IF(sep .EQ. 0) THEN  

        !! (3) Contract TN for sep=0
        CALL mult_bvec_bvec(expval,    RbO1O2_sep_00,  Lb)

     ELSEIF(sep .EQ. 1) THEN

        !! (4) Contract TN for sep=1
        CALL mult_bvec_bvec(expval,    RbO1O2_sep_11,  Lb)

     ELSEIF(sep .EQ. 2) THEN

        !! (5) Contract TN for sep=EVEN: If [sep = 2] --> no need for padding, just contract the boundary vecs
        CALL mult_bvec_bvec(expval,    RbO2_sep_EVEN,  LbO1)

     ELSEIF(isEVEN(sep)) THEN

        !! (6) Contract TN for sep=EVEN: If [sep > 2] --> mult in extra padding Tmpo at each even sep --> augment LbO1 boundary
        CALL mult_bvec_tmpo(LbO1, Tmpo, 'L')
        CALL mult_bvec_bvec(expval,    RbO2_sep_EVEN,  LbO1)

     ELSEIF(isODD(sep)) THEN

        !! (7) Contract TN for sep=ODD
        CALL mult_bvec_bvec(expval,    RbO2_sep_ODD,   LbO1)
     ELSE
        CALL invalid_value("contract_TN_2P_corr -- invalid sep ", sep)
     end IF

     !! Write expval
     expval         = expval / vec_norm
     expvals(sep+1) = expval / C_norm**(sep/2 + 1)

  END DO separation_loop

 END SUBROUTINE contract_TN_2P_corr

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!















 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Construct TN bounded by CTM, represented by Rvec -- Transfer MPO -- Lvec !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Construct LHS boundary with an operator-carrying Tmpo absorbed !!!
 SUBROUTINE calc_left_bound_OP(Lb, Lvec, Tmpo, TN_opSites, opFlag)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT)        :: Lb(:,:,:)
  COMPLEX(KIND=DP),              INTENT(IN)           :: Lvec(:,:,:)
  TYPE(transfer_mpo),            INTENT(IN)           :: Tmpo(:)
  TYPE(op_sites_type),           INTENT(IN), OPTIONAL :: TN_opSites
  CHARACTER(LEN=*),              INTENT(IN), OPTIONAL :: opFlag

  !! The leftmost transfer MPO
  TYPE(transfer_mpo), ALLOCATABLE :: TmpoLHS(:)

  !! Create local copies of transfer MPO && Lvec
  CALL copy_transfer_mpo(TmpoLHS, Tmpo)
  CALL copyTens(Lb, Lvec)

  !! If operator is present: 
  !! -- insert operator into transfer mpo 
  !! -- absorb [transfer mpo with operator] into Lvec
  IF(PRESENT(TN_opSites)) THEN
     CALL insert_operator_site(TmpoLHS, TN_opSites, opFlag=OptArg(opFlag, 'op1'))
     CALL mult_bvec_tmpo(Lb, TmpoLHS, 'L')
  end IF

 END SUBROUTINE calc_left_bound_OP






 !!! Construct RHS boundary with an operator-carrying Tmpo absorbed !!!
 SUBROUTINE calc_right_bound_OP(Rb, Rvec, Tmpo, TN_opSites, sep)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: Rb(:,:,:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: Rvec(:,:,:)
  TYPE(transfer_mpo),            INTENT(IN)    :: Tmpo(:)
  TYPE(op_sites_type),           INTENT(IN)    :: TN_opSites
  INTEGER,                       INTENT(IN)    :: sep

  !! The rightmost transfer MPO
  TYPE(transfer_mpo), ALLOCATABLE :: TmpoRHS(:)

  !! (1) Create local copies of transfer MPO && Rvec
  CALL copy_transfer_mpo(TmpoRHS, Tmpo)
  CALL copyTens(Rb, Rvec)

  !! (2) Insert operator(s) into TmpoRHS
  IF(sep .EQ. 0) THEN

     CALL insert_operator_site(TmpoRHS, TN_opSites, opFlag='op12')

  ELSEIF(sep .EQ. 1) THEN

     CALL insert_operator_site(TmpoRHS, TN_opSites, opFlag='op1')
     CALL insert_operator_site(TmpoRHS, TN_opSites, opFlag='op2-ODD')

  ELSEIF(isEVEN(sep)) THEN

     CALL insert_operator_site(TmpoRHS, TN_opSites, opFlag='op2-EVEN')

  ELSEIF(isODD(sep)) THEN

     CALL insert_operator_site(TmpoRHS, TN_opSites, opFlag='op2-ODD')
  ELSE
     CALL invalid_value("calc_right_bound_OP -- invalid sep ", sep)
  end IF

  !! (3) Absorb [Tmpo with operator] into RHS boundary
  CALL mult_bvec_tmpo(Rb, TmpoRHS, 'R')

 END SUBROUTINE calc_right_bound_OP






 !!! Construct boundary vecs && transfer MPO !!!
 SUBROUTINE construct_Rvec_Tmpo_Lvec_network(Rvec, Tmpo, Lvec, Cmat, Tmat, TN2DIn, CTMRG)

  COMPLEX(KIND=DP),   ALLOCATABLE, INTENT(INOUT) :: Rvec(:,:,:), Lvec(:,:,:)
  TYPE(transfer_mpo), ALLOCATABLE, INTENT(INOUT) :: Tmpo(:)
  TYPE(ctm_corner_type),           INTENT(IN)    :: Cmat(:)
  TYPE(ctm_transfer_type),         INTENT(IN)    :: Tmat(:)
  TYPE(block_mpo),                 INTENT(IN)    :: TN2DIn(:)
  TYPE(ctm_params),                INTENT(IN)    :: CTMRG  !CTM params          

  !! Environment of 2D network && Local copy of 2D network
  TYPE(block_mps), ALLOCATABLE  :: ENV2D(:)
  TYPE(block_mpo), ALLOCATABLE  :: TN2D(:)

  !! Indices && dims
  INTEGER :: s
  INTEGER, PARAMETER :: N_sites = 2

  !! Rescale 2D network
  CALL copy_mpo_block(TN2D, TN2DIn)

  !! Calculate single-bond environment of iPEPS
  CALL calc_single_bond_environment(ENV2D, Cmat, Tmat, TN2D, CTMRG, DIR='S', pbc=.FALSE.)

  !!! Construct transfer MPO from ENV2D && TN2D
  CALL allocate_empty_transfer_mpo(Tmpo, N_sites)
  DO s=1,N_sites
     CALL copyTens(Tmpo(s) % U,  ENV2D(s)   % m)
     CALL copyTens(Tmpo(s) % M,  TN2D(s)    % m)
     CALL copyTens(Tmpo(s) % D,  ENV2D(s+2) % m)
  end DO

  !!! Construct boundary vectors from ENV2D
  CALL copyTens(Lvec,  ENV2D(5) % m)
  CALL copyTens(Rvec,  ENV2D(6) % m)

  !! Normalize  <Lvec|Tmpo|Rvec> network
  CALL normalize_Rvec_Tmpo_Lvec_network(Rvec, Tmpo, Lvec)

 END SUBROUTINE construct_Rvec_Tmpo_Lvec_network





 !!! Normalize <Lvec|Tmpo|Rvec> network !!!
 SUBROUTINE normalize_Rvec_Tmpo_Lvec_network(Rvec, Tmpo, Lvec)

  COMPLEX(KIND=DP),   ALLOCATABLE, INTENT(INOUT) :: Rvec(:,:,:), Lvec(:,:,:)
  TYPE(transfer_mpo), ALLOCATABLE, INTENT(INOUT) :: Tmpo(:)

  !! Norm constants && sites
  COMPLEX(KIND=DP)   :: C_norm, vec_norm
  INTEGER            :: s, N_sites

  !! Get size of Tmpo
  N_sites = SIZE(Tmpo)

  !! Normalize Rvec, Lvec boundary vectors
  CALL mult_bvec_bvec(vec_norm,      Rvec,       Lvec)
  Rvec = Rvec / sqrt(vec_norm)
  Lvec = Lvec / sqrt(vec_norm)

  !! Find C_norm of <Lvec|Tmpo|Rvec>
  CALL mult_bvec_tmpo_bvec(C_norm,   Rvec, Tmpo, Lvec)

  !! Normalize UP && DN sites of Tmpo 
  !! -- since the central part MID will change when we insert an operator site
  !! -- and since UP && DN, like Rvec && Lvec correspond to CTM environment while MID corresponds to bulk TN 
  DO s=1,N_sites
     Tmpo(s) % U = Tmpo(s) % U / C_norm**(0.5D0/REAL(N_sites, KIND=DP)) 
     Tmpo(s) % D = Tmpo(s) % D / C_norm**(0.5D0/REAL(N_sites, KIND=DP))
  end DO

 END SUBROUTINE normalize_Rvec_Tmpo_Lvec_network


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!















 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Calculate environment of a single-bond in 2D network !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Calculate single-bond CTM environment of a general 2D network !!!
 SUBROUTINE calc_single_bond_environment(ENV2D, CmatIn, TmatIn, TN2DIn, CTMRG, DIR, pbc)

  TYPE(block_mps),         ALLOCATABLE, INTENT(INOUT) :: ENV2D(:)           !Environment
  TYPE(ctm_corner_type),                INTENT(IN)    :: CmatIn(:)          !Original CTM Cmat
  TYPE(ctm_transfer_type),              INTENT(IN)    :: TmatIn(:)          !Original CTM Tmat
  TYPE(block_mpo),                      INTENT(IN)    :: TN2DIn(:)          !2D network 
  TYPE(ctm_params),                     INTENT(IN)    :: CTMRG              !CTM params
  CHARACTER(LEN=*),                     INTENT(IN)    :: DIR                !bond whose environment we're calculating
  LOGICAL,                              INTENT(IN)    :: pbc                !whether to represent environment as periodic or standard MPS

  !! CTM && TN2D local copies 
  TYPE(ctm_corner_type),    ALLOCATABLE :: Cmat(:) 
  TYPE(ctm_transfer_type),  ALLOCATABLE :: Tmat(:)
  TYPE(block_mpo),          ALLOCATABLE :: TN2D(:) 

  !! CTM sites && CTMRG method
  INTEGER :: sA, sB
  INTEGER :: xA, xB, yA, yB
  CHARACTER(LEN=1) :: oppDIR

  !! (1) Create local copy of CTM && 2D network
  CALL copy_CTM(Cmat, Tmat, CmatIn, TmatIn)
  CALL copy_mpo_block(TN2D, TN2DIn)

  !! (2) Get iPEPS sites for a given bond && the opposite bond to the one being updated
  CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, DIR=DIR)
  oppDIR = oppositeDIR(DIR)


  !! (3) To obtain a single-link environment, absorb an extra block into 'N'-edge of CTM  
  !!     -- to maintain consistency, use the same CTM method as we did in CTM calculation
  !!     -- we must unswap CTM sites on lateral edges (to the one being updated) 
  !!        since they're swapped during every CTM move to reconfigure the network for the next move
  SELECT CASE(CTMRG % ctm_method)
  CASE('1')       
       CALL perform_CTM_move_oneDir(Cmat, Tmat, TN2D, CTMRG, DIR = oppDIR)
       CALL swap_CTM_block_sites(Tmat,          TN2D,        DIR = oppDIR)                
  CASE('2', '3')          
       CALL perform_CTM_move_twoDir(Cmat, Tmat, TN2D, CTMRG, DIR = (/DIR, oppDIR/)) 
       CALL swap_CTM_block_sites(Tmat,          TN2D,        DIR = oppDIR)
       CALL copy_CTM_edge(Cmat, Tmat, CmatIn, TmatIn,        DIR = DIR)                                     
  CASE DEFAULT
       CALL invalid_flag("calc_single_bond_environment -- invalid ctm_method ", CTMRG % ctm_method)
  end SELECT


  !! (4) Compute environment from CTM tensors (using iPEPS site indexing rather than the CTM one)
  CALL allocate_empty_mps_block(ENV2D, 6)

  SELECT CASE(DIR)
  CASE('S', 'N')

       CALL allocate_mps_site(ENV2D(5),  TENSMUL(TENSMUL(Tmat(4) % T(sA) % m,  Cmat(1) % m,  '22'),  Cmat(4) % m,  '31'))
       CALL allocate_mps_site(ENV2D(6),  TENSMUL(TENSMUL(Tmat(2) % T(sB) % m,  Cmat(3) % m,  '22'),  Cmat(2) % m,  '31')) 

       CALL allocate_mps_site(ENV2D(1),  Tmat(3) % T(2) % m)
       CALL allocate_mps_site(ENV2D(2),  Tmat(3) % T(1) % m)
       CALL allocate_mps_site(ENV2D(3),  Tmat(1) % T(1) % m)
       CALL allocate_mps_site(ENV2D(4),  Tmat(1) % T(2) % m)

  CASE('W', 'E')

       CALL allocate_mps_site(ENV2D(5),  TENSMUL(TENSMUL(Tmat(1) % T(sA) % m,  Cmat(2) % m,  '22'),  Cmat(1) % m,  '31'))
       CALL allocate_mps_site(ENV2D(6),  TENSMUL(TENSMUL(Tmat(3) % T(sB) % m,  Cmat(4) % m,  '22'),  Cmat(3) % m,  '31')) 

       CALL allocate_mps_site(ENV2D(1),  Tmat(4) % T(1) % m)
       CALL allocate_mps_site(ENV2D(2),  Tmat(4) % T(2) % m)
       CALL allocate_mps_site(ENV2D(3),  Tmat(2) % T(2) % m)
       CALL allocate_mps_site(ENV2D(4),  Tmat(2) % T(1) % m)
  CASE DEFAULT
       CALL invalid_flag("calc_single_bond_environment -- invalid DIR ", DIR)
  end SELECT


  !! (5) By default, environment is represented by PBC MPS
  !!     -- if we want a standard MPS instead , transpose bottom && right sites of environment
  IF(.NOT. pbc) THEN
     CALL allocate_mps_site(ENV2D(6),  TensTRANSPOSE(ENV2D(6) % m,  '23'))
     CALL allocate_mps_site(ENV2D(3),  TensTRANSPOSE(ENV2D(3) % m,  '23'))
     CALL allocate_mps_site(ENV2D(4),  TensTRANSPOSE(ENV2D(4) % m,  '23'))
  end IF
  
 END SUBROUTINE calc_single_bond_environment


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Construct transfer MPO with operators inserted !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Insert operator sites into transfer MPO
 SUBROUTINE insert_operator_site(Tmpo, TN_opSites, opFlag)

  TYPE(transfer_mpo),            INTENT(INOUT) :: Tmpo(:)
  TYPE(op_sites_type),           INTENT(IN)    :: TN_opSites
  CHARACTER(LEN=*),              INTENT(IN)    :: opFlag

  !! Allocate operator site on transfer MPO
  SELECT CASE(opFlag)
  CASE('op1')  
       CALL insert_operator_site_INTERNAL(Tmpo, TN_OpSites % site1_Op1,     opSite=1) 
  CASE('op2')  
       CALL insert_operator_site_INTERNAL(Tmpo, TN_OpSites % site1_Op1,     opSite=2) 
  CASE('op12')
       CALL insert_operator_site_INTERNAL(Tmpo, TN_OpSites % site1_Op1_Op2, opSite=1)
  CASE('op2-EVEN')
       CALL insert_operator_site_INTERNAL(Tmpo, TN_OpSites % site1_Op2,     opSite=1)
  CASE('op2-ODD')
       CALL insert_operator_site_INTERNAL(Tmpo, TN_OpSites % site2_Op2,     opSite=2)
  CASE DEFAULT
       CALL invalid_flag("insert_operator_site -- invalid opFlag", opFlag)
  end SELECT

 END SUBROUTINE insert_operator_site




 SUBROUTINE insert_operator_site_INTERNAL(Tmpo, opSiteTens, opSite)

  TYPE(transfer_mpo),            INTENT(INOUT) :: Tmpo(:)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(IN)    :: opSiteTens(:,:,:,:)
  INTEGER,                       INTENT(IN)    :: opSite

  !! Check if operator site tensor is allocated
  IF(call_exit(.NOT. ALLOCATED(opSiteTens), "insert_operator_site_INTERNAL: opSiteTens not allocated")) STOP

  !! Allocate operator site on transfer MPO
  CALL copyTens(Tmpo(opSite) % M, opSiteTens)

 END SUBROUTINE insert_operator_site_INTERNAL

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 







 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Contract unit cell of 2D Network surrounded by CTM environment !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Contract CTM network !!!
 SUBROUTINE contract_CTM_network(expval, CmatIn, TmatIn, TN2D)

  COMPLEX(KIND=DP),         INTENT(INOUT) :: expval      !expval = result of contraction
  TYPE(ctm_corner_type),    INTENT(IN)    :: CmatIn(:)   !CTM Cmat
  TYPE(ctm_transfer_type),  INTENT(IN)    :: TmatIn(:)   !CTM Tmat
  TYPE(block_mpo),          INTENT(IN)    :: TN2D(:)     !2D network

  TYPE(ctm_corner_type),   ALLOCATABLE :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE :: Tmat(:)
  COMPLEX(KIND=DP),        ALLOCATABLE :: theta(:,:)

  !! CTM sites
  INTEGER :: N_sites, N_sub

  !! Create CTM local copy
  CALL copy_CTM(Cmat, Tmat, CmatIn, TmatIn)
  N_sites = SIZE(Cmat)
  N_sub   = SIZE(Tmat(1) % T)

  !! Contract CTM network into four corners
  CALL absorb_extra_block_into_CTM(Cmat, Tmat, TN2D, 'S')
  CALL absorb_extra_block_into_CTM(Cmat, Tmat, TN2D, 'N')

  CALL self_contract_CTM_edge(Cmat, Tmat, 'S')
  CALL self_contract_CTM_edge(Cmat, Tmat, 'N')

  !! Calc CTM corner product (with cut placed anywhere -- doesn't matter cause we take a trace)
  CALL calc_CTM_corner_product(theta, Cmat, DIR='W')

  !! Compute expval (trace over product of corner matrices)
  expval = TRACE(theta)

 END SUBROUTINE contract_CTM_network


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE TN_contractions
