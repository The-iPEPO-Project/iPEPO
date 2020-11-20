MODULE project_ipepo

 USE utility
 USE definitions_mps_mpo
 USE definitions_peps_pepo
 USE mps_mpo_utility

 IMPLICIT NONE

 interface project_ipepo_to_impo
    module procedure project_ipepo_to_impo_GamLam
    module procedure project_ipepo_to_impo_PLAIN
 end interface project_ipepo_to_impo

 interface unproject_ipepo_from_impo
    module procedure unproject_ipepo_from_impo_GamLam
    module procedure unproject_ipepo_from_impo_PLAIN
 end interface unproject_ipepo_from_impo

 interface calc_factorized_ipepo
    module procedure calc_factorized_ipepo_GamLam
    module procedure calc_factorized_ipepo_PLAIN
 end interface calc_factorized_ipepo

 interface restore_unfactorized_ipepo
    module procedure restore_unfactorized_ipepo_GamLam
    module procedure restore_unfactorized_ipepo_PLAIN
 end interface restore_unfactorized_ipepo

 private !hides all items not listed in public statement 
 public project_ipepo_to_impo, unproject_ipepo_from_impo
 public calc_factorized_ipepo, restore_unfactorized_ipepo
 public absorb_ipepo_lambdas,  restore_ipepo_lambdas

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (1) PROJECT iPEPO TO iMPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Project ipepo onto impo: fuse Y-axis legs with X-axis legs in the XY plane !!!
 SUBROUTINE project_ipepo_to_impo_GamLam(iGmpo, iLmpo, iGpepo, iLpepo, DIR)

  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT) :: iGmpo(:)            
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLmpo(:)
  TYPE(block_pepo),                INTENT(IN)    :: iGpepo(:)            
  TYPE(block_lambda),              INTENT(IN)    :: iLpepo(:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! Local ipepo copy
  TYPE(block_pepo), ALLOCATABLE :: ipepo(:)

  !! iPEPO sites && reshaping flags
  CHARACTER(LEN=16) :: legsA, legsB
  INTEGER           :: sA, sB
  INTEGER           :: xA, xB, yA, yB

  !! Get iPEPO sites
  CALL get_ipepo_sites(sA, sB, xA, xB, yA, yB, DIR)

  !! Set reshaping flag
  SELECT CASE(DIR)
  CASE('S', 'N')
       legsA = '1,2,345,6'; legsB = '1,2,5,346'
  CASE('E', 'W') 
       legsA = '1,2,356,4'; legsB = '1,2,3,456'
  CASE DEFAULT
       CALL invalid_flag("project_ipepo_to_impo_GamLam -- invalid DIR ", DIR)
  end SELECT

  !! Absorb lambdas -- external lambdas lying on both Y-axis and X-axis (keep only one internal lambda)
  CALL absorb_ipepo_lambdas(ipepo, iGpepo, iLpepo, DIR, absX=.True., absY=.True.) 

  !! Reshape PEPO gammas into MPO gammas
  CALL allocate_empty_mpo_block(iGmpo, 2)
  CALL allocate_mpo_site(iGmpo(1),  RESHAPE_4D(ipepo(sA) % m, legsA))
  CALL allocate_mpo_site(iGmpo(2),  RESHAPE_4D(ipepo(sB) % m, legsB))

  !! Write PEPO central lambda to MPO lambda
  CALL allocate_empty_lambda_block(iLmpo, 1)
  CALL allocate_lambda_site(iLmpo(1), iLpepo(xA) % m)

 END SUBROUTINE project_ipepo_to_impo_GamLam






 !!! Project ipepo onto impo: fuse Y-axis legs with X-axis legs in the XY plane !!!
 SUBROUTINE project_ipepo_to_impo_PLAIN(impo, ipepo, DIR)

  TYPE(block_mpo),  ALLOCATABLE, INTENT(INOUT) :: impo(:)            
  TYPE(block_pepo),              INTENT(IN)    :: ipepo(:)            
  CHARACTER(LEN=*),              INTENT(IN)    :: DIR

  !! iPEPO sites && reshaping flags
  CHARACTER(LEN=16) :: legsA, legsB
  INTEGER           :: sA, sB
  INTEGER           :: xA, xB, yA, yB

  !! Get iPEPO sites
  CALL get_ipepo_sites(sA, sB, xA, xB, yA, yB, DIR)

  !! Set reshaping flag
  SELECT CASE(DIR)
  CASE('S', 'N')
       legsA = '1,2,345,6'; legsB = '1,2,5,346'
  CASE('E', 'W') 
       legsA = '1,2,356,4'; legsB = '1,2,3,456'
  CASE DEFAULT
       CALL invalid_flag("project_ipepo_to_impo_GamLam -- invalid DIR ", DIR)
  end SELECT

  !! Reshape PEPO into MPO
  CALL allocate_empty_mpo_block(impo, 2)
  CALL allocate_mpo_site(impo(1),  RESHAPE_4D(ipepo(sA) % m, legsA))
  CALL allocate_mpo_site(impo(2),  RESHAPE_4D(ipepo(sB) % m, legsB))

 END SUBROUTINE project_ipepo_to_impo_PLAIN

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (2) RESTORE iPEPO FROM iMPO PROJECTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Restore ipepo from impo projection: split external legs into Y-axis && X-axis legs in the XY plane !!!
 SUBROUTINE unproject_ipepo_from_impo_GamLam(iGpepo, iLpepo, iGmpo, iLmpo, DIR)

  TYPE(block_pepo),   ALLOCATABLE, INTENT(INOUT) :: iGpepo(:)            
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLpepo(:)
  TYPE(block_mpo),                 INTENT(IN)    :: iGmpo(:)            
  TYPE(block_lambda),              INTENT(IN)    :: iLmpo(:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! iPEPO sites && reshaping flags
  CHARACTER(LEN=16) :: legsA, legsB
  INTEGER           :: dimsA(3), dimsB(3)
  INTEGER           :: sA, sB
  INTEGER           :: xA, xB, yA, yB

  !! Get iPEPO sites 
  !! IMPORTANT: do not allocate new empty iPEPo blocks when restoring iPEPO from iMPO projection, 
  !! since old lambdas will be needed by restore_ipepo_lambdas()
  CALL get_ipepo_sites(sA, sB, xA, xB, yA, yB, DIR)

  !! Set reshaping flag && dims
  SELECT CASE(DIR)
  CASE('S', 'N')
       legsA = '1,2,345,6'; dimsA = (/ iGpepo(sA) % Sdim, iGpepo(sA) % Ndim, iGpepo(sA) % Wdim /)
       legsB = '1,2,5,346'; dimsB = (/ iGpepo(sB) % Sdim, iGpepo(sB) % Ndim, iGpepo(sB) % Edim /)                     
  CASE('E', 'W') 
       legsA = '1,2,356,4'; dimsA = (/ iGpepo(sA) % Sdim, iGpepo(sA) % Wdim, iGpepo(sA) % Edim /)
       legsB = '1,2,3,456'; dimsB = (/ iGpepo(sB) % Ndim, iGpepo(sB) % Wdim, iGpepo(sB) % Edim /)
  CASE DEFAULT
       CALL invalid_flag("unproject_ipepo_from_impo_GamLam -- invalid DIR ", DIR)
  end SELECT

  !! Reshape iMPO projection into iPEPO 
  CALL allocate_pepo_site(iGpepo(sA),    RESHAPE_6D(iGmpo(1) % m, legsA, dimsA))
  CALL allocate_pepo_site(iGpepo(sB),    RESHAPE_6D(iGmpo(2) % m, legsB, dimsB))
  CALL allocate_lambda_site(iLpepo(xA),  iLmpo(1) % m) 

  !! Restore lambdas -- i.e. external lambdas lying on both Y-axis and X-axis
  CALL restore_ipepo_lambdas(iGpepo, iLpepo, DIR, resX=.True., resY=.True.)

 END SUBROUTINE unproject_ipepo_from_impo_GamLam





 !!! Restore ipepo from impo projection: split external legs into Y-axis && X-axis legs in the XY plane !!!
 SUBROUTINE unproject_ipepo_from_impo_PLAIN(ipepo, impo, DIR)

  TYPE(block_pepo), ALLOCATABLE, INTENT(INOUT) :: ipepo(:)            
  TYPE(block_mpo),               INTENT(IN)    :: impo(:)            
  CHARACTER(LEN=*),              INTENT(IN)    :: DIR

  !! iPEPO sites && reshaping flags
  CHARACTER(LEN=16) :: legsA, legsB
  INTEGER           :: dimsA(3), dimsB(3)
  INTEGER           :: sA, sB
  INTEGER           :: xA, xB, yA, yB

  !! Get iPEPO sites 
  CALL get_ipepo_sites(sA, sB, xA, xB, yA, yB, DIR)

  !! Set reshaping flag && dims
  SELECT CASE(DIR)
  CASE('S', 'N')
       legsA = '1,2,345,6'; dimsA = (/ ipepo(sA) % Sdim, ipepo(sA) % Ndim, ipepo(sA) % Wdim /)
       legsB = '1,2,5,346'; dimsB = (/ ipepo(sB) % Sdim, ipepo(sB) % Ndim, ipepo(sB) % Edim /)                     
  CASE('E', 'W') 
       legsA = '1,2,356,4'; dimsA = (/ ipepo(sA) % Sdim, ipepo(sA) % Wdim, ipepo(sA) % Edim /)
       legsB = '1,2,3,456'; dimsB = (/ ipepo(sB) % Ndim, ipepo(sB) % Wdim, ipepo(sB) % Edim /)
  CASE DEFAULT
       CALL invalid_flag("unproject_ipepo_from_impo_PLAIN -- invalid DIR ", DIR)
  end SELECT

  !! Reshape iMPO projection into iPEPO
  CALL allocate_pepo_site(ipepo(sA),  RESHAPE_6D(impo(1) % m, legsA, dimsA))
  CALL allocate_pepo_site(ipepo(sB),  RESHAPE_6D(impo(2) % m, legsB, dimsB))

 END SUBROUTINE unproject_ipepo_from_impo_PLAIN

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (3) CALCULATE FACTORIZED iPEPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Project ipepo to impo, and factorize it as [impo --> imps * X-site] !!!
 SUBROUTINE calc_factorized_ipepo_GamLam(iGmps, iLmps, Xa, Yb, iGpepo, iLpepo, DIR)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGmps(:)            
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLmps(:)
  COMPLEX(KIND=DP),   ALLOCATABLE, INTENT(INOUT) :: Xa(:,:,:), Yb(:,:,:)
  TYPE(block_pepo),                INTENT(IN)    :: iGpepo(:)            
  TYPE(block_lambda),              INTENT(IN)    :: iLpepo(:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! Temporary iPEPO projection onto iMPO && temp factorized MPO tensors
  TYPE(block_mpo),  ALLOCATABLE :: iGamma_mpo(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: mps11(:,:,:), mps22(:,:,:)  

  !! (1) Project iPEPO onto iMPO
  CALL project_ipepo_to_impo(iGamma_mpo, iLmps, iGpepo, iLpepo, DIR)

  !! (2) Factorize each site of iMPO into impoA -> Xa * impsA, impoB -> impsB * Yb
  CALL allocate_empty_mps_block(iGmps, 2)

  CALL LHS_factorize_Tens4D(mps11,  Xa,  iGamma_mpo(1) % m,  '23,14')
  CALL allocate_mps_site(iGmps(1),  mps11)

  CALL RHS_factorize_Tens4D(mps22,  Yb,  iGamma_mpo(2) % m,  '13,24')
  CALL allocate_mps_site(iGmps(2),  mps22)

 END SUBROUTINE calc_factorized_ipepo_GamLam





 !!! Project ipepo to impo, and factorize it as [impo --> imps * X-site] !!!
 SUBROUTINE calc_factorized_ipepo_PLAIN(imps, Xa, Yb, ipepo, DIR)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: imps(:)            
  COMPLEX(KIND=DP),   ALLOCATABLE, INTENT(INOUT) :: Xa(:,:,:), Yb(:,:,:)
  TYPE(block_pepo),                INTENT(IN)    :: ipepo(:)            
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! Temporary iPEPO projection onto iMPO && temp factorized MPO tensors
  TYPE(block_mpo),  ALLOCATABLE :: impoTmp(:) 
  COMPLEX(KIND=DP), ALLOCATABLE :: mps11(:,:,:), mps22(:,:,:)  

  !! (1) Project iPEPO onto iMPO
  CALL project_ipepo_to_impo(impoTmp, ipepo, DIR)

  !! (2) Factorize each site of iMPO into impoA -> Xa * impsA, impoB -> impsB * Yb
  CALL allocate_empty_mps_block(imps, 2)

  CALL LHS_factorize_Tens4D(mps11,  Xa,  impoTmp(1) % m,  '23,14')
  CALL allocate_mps_site(imps(1),  mps11)

  CALL RHS_factorize_Tens4D(mps22,  Yb,  impoTmp(2) % m,  '13,24')
  CALL allocate_mps_site(imps(2),  mps22)

 END SUBROUTINE calc_factorized_ipepo_PLAIN

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (4) RESTORE UNFACTORIZED iPEPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Multiply factorized iMPS-SITE and X-SITE as [imps * X -> impo], 
 !!! and restore iPEPO from iMPO projection
 SUBROUTINE restore_unfactorized_ipepo_GamLam(iGpepo, iLpepo, iGmps, iLmps, Xa, Yb, DIR)

  TYPE(block_pepo),   ALLOCATABLE, INTENT(INOUT) :: iGpepo(:)            
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLpepo(:)
  TYPE(block_mps),                 INTENT(IN)    :: iGmps(:)            
  TYPE(block_lambda),              INTENT(IN)    :: iLmps(:)
  COMPLEX(KIND=DP),                INTENT(IN)    :: Xa(:,:,:), Yb(:,:,:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! Temporary iPEPO projection onto iMPO
  TYPE(block_mpo), ALLOCATABLE :: iGmpoTmp(:)  

  !! (1) Multiply factorized iPEPO sites (Xa * iGmps(1) --> iGmpo(1) gives SN dims wrong way round hence transpose)
  CALL allocate_empty_mpo_block(iGmpoTmp, 2)
  CALL allocate_mpo_site(iGmpoTmp(1), TensTRANSPOSE(TENSMUL(iGmps(1) % m, Xa, '23'), '12'))
  CALL allocate_mpo_site(iGmpoTmp(2),               TENSMUL(iGmps(2) % m, Yb, '32'))

  !! (2) Restore iPEPO from iMPO projection
  CALL unproject_ipepo_from_impo(iGpepo, iLpepo, iGmpoTmp, iLmps, DIR)

 END SUBROUTINE restore_unfactorized_ipepo_GamLam




 !!! Multiply factorized iMPS-SITE and X-SITE as [imps * X -> impo], 
 !!! and restore iPEPO from iMPO projection
 SUBROUTINE restore_unfactorized_ipepo_PLAIN(ipepo, imps, Xa, Yb, DIR)

  TYPE(block_pepo),   ALLOCATABLE, INTENT(INOUT) :: ipepo(:)            
  TYPE(block_mps),                 INTENT(IN)    :: imps(:)            
  COMPLEX(KIND=DP),                INTENT(IN)    :: Xa(:,:,:), Yb(:,:,:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! Temporary iPEPO projection onto iMPO
  TYPE(block_mpo), ALLOCATABLE :: impoTmp(:)  

  !! (1) Multiply factorized iPEPO sites (Xa * iGmps(1) --> iGmpo(1) gives SN dims wrong way round hence transpose)
  CALL allocate_empty_mpo_block(impoTmp, 2)
  CALL allocate_mpo_site(impoTmp(1), TensTRANSPOSE(TENSMUL(imps(1) % m, Xa, '23'), '12'))
  CALL allocate_mpo_site(impoTmp(2),               TENSMUL(imps(2) % m, Yb, '32'))

  !! (2) Restore iPEPO from iMPO projection
  CALL unproject_ipepo_from_impo(ipepo, impoTmp, DIR)

 END SUBROUTINE restore_unfactorized_ipepo_PLAIN

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (5) ABSORB/RESTORE iPEPO LAMBDAS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Absorb lambdas on Y-axis, external lambda on X-axis, or absorb all lambdas in symmetric fashion !!!
 SUBROUTINE absorb_ipepo_lambdas(ipepo, iGamma, iLambda, DIR, absX, absY)

  TYPE(block_pepo),   ALLOCATABLE, INTENT(INOUT)           :: ipepo(:)            !iPEPO with absorbed lambdas
  TYPE(block_pepo),                INTENT(IN),   OPTIONAL  :: iGamma(:)           !iPEPO gamma
  TYPE(block_lambda),              INTENT(IN)              :: iLambda(:)          !iPEPO lambda
  CHARACTER(LEN=*),                INTENT(IN),   OPTIONAL  :: DIR                 !central bond direction/orientation
  LOGICAL,                         INTENT(IN),   OPTIONAL  :: absX, absY          !whether to absorb lambdas lying on 
                                                                                  !tangential && normal axes
  !! iPEPO sites
  INTEGER :: sA, sB
  INTEGER :: xA, xB, yA, yB

  !! Copy iGamma to iPepo block (if iGamma not present, use iPEPO as input Gamma instead)
  IF(PRESENT(iGamma)) CALL copy_pepo_block(ipepo, iGamma) 

  !! Absorb lambdas
  IF(PRESENT(absX) .AND. PRESENT(absY) .AND. PRESENT(DIR)) THEN 

     !! Get iPEPO sites (LHS = S,W), (RHS = N,E)
     CALL get_ipepo_sites(sA, sB, xA, xB, yA, yB, DIR)

     !! Absorb external lambdas ONLY
     SELECT CASE(DIR)
     CASE('S', 'N')

         !! Absorb external lambdas on Y-axis (into siteA)
         IF(absX) CALL LambdaMUL(ipepo(sA) % m,  iLambda(xB) % m,  '5') !LHS
         IF(absY) CALL LambdaMUL(ipepo(sA) % m,  iLambda(yA) % m,  '4') !RHS 
         IF(absY) CALL LambdaMUL(ipepo(sA) % m,  iLambda(yB) % m,  '3') !LHS 

         !! Absorb external lambdas on Y-axis (into siteB)
         IF(absX) CALL LambdaMUL(ipepo(sB) % m,  iLambda(xB) % m,  '6') !RHS
         IF(absY) CALL LambdaMUL(ipepo(sB) % m,  iLambda(yB) % m,  '4') !RHS
         IF(absY) CALL LambdaMUL(ipepo(sB) % m,  iLambda(yA) % m,  '3') !LHS

     CASE('E', 'W')

         !! Absorb external lambdas on Y-axis (into siteA)
         IF(absX) CALL LambdaMUL(ipepo(sA) % m,  iLambda(xB) % m,  '3') !LHS
         IF(absY) CALL LambdaMUL(ipepo(sA) % m,  iLambda(yA) % m,  '6') !RHS
         IF(absY) CALL LambdaMUL(ipepo(sA) % m,  iLambda(yB) % m,  '5') !LHS

         !! Absorb external lambdas on Y-axis (into siteB)
         IF(absX) CALL LambdaMUL(ipepo(sB) % m,  iLambda(xB) % m,  '4') !RHS
         IF(absY) CALL LambdaMUL(ipepo(sB) % m,  iLambda(yB) % m,  '6') !RHS
         IF(absY) CALL LambdaMUL(ipepo(sB) % m,  iLambda(yA) % m,  '5') !LHS

     CASE DEFAULT
         CALL invalid_flag("absorb_ipepo_lambdas -- invalid DIR ", DIR)
     end SELECT

  ELSE

     !! Get iPEPO sites for the default configuration (centred on 'S' bond)
     CALL get_ipepo_sites(sA, sB, xA, xB, yA, yB, 'S')

     !! Absorb all lambdas in symmetric fashion -- on SN bonds 
     CALL LambdaMUL(ipepo(sA) % m,  iLambda(xA) % m,  ipepo(sB) % m,  '34')
     CALL LambdaMUL(ipepo(sB) % m,  iLambda(xB) % m,  ipepo(sA) % m,  '34')
 
     !! Absorb all lambdas in symmetric fashion -- on EW bonds
     CALL LambdaMUL(ipepo(sA) % m,  iLambda(yA) % m,  ipepo(sB) % m,  '56')
     CALL LambdaMUL(ipepo(sB) % m,  iLambda(yB) % m,  ipepo(sA) % m,  '56')
  end IF 

 END SUBROUTINE absorb_ipepo_lambdas






 !!! Restore external iPEPO lambdas: convert from nolambda form to Gamma-Lambda form !!!
 SUBROUTINE restore_ipepo_lambdas(iGamma, iLambda, DIR, resX, resY)

  TYPE(block_pepo),   INTENT(INOUT) :: iGamma(:)             !INPUT iPEPO: with all external Lambdas absorbed
  TYPE(block_lambda), INTENT(IN)    :: iLambda(:)            !INPUT iLambdas: input old external Lambdas so we can restore them
  CHARACTER(LEN=*),   INTENT(IN)    :: DIR                   !central bond direction/orientation
  LOGICAL,            INTENT(IN)    :: resX, resY            !external lambdas to be restored

  !! iPEPO sites
  INTEGER :: sA, sB
  INTEGER :: xA, xB, yA, yB

  !! Get iPEPO sites (LEFT = S,W --- RIGHT = N,E)
  CALL get_ipepo_sites(sA, sB, xA, xB, yA, yB, DIR)

  !! Restore external lambdas
  SELECT CASE(DIR)
  CASE('S', 'N')

      !! Restore external lambdas (from site-AA)
      IF(resX) CALL LambdaDIV(iGamma(sA) % m,  iLambda(xB) % m,  '5') !LHS
      IF(resY) CALL LambdaDIV(iGamma(sA) % m,  iLambda(yA) % m,  '4') !RHS
      IF(resY) CALL LambdaDIV(iGamma(sA) % m,  iLambda(yB) % m,  '3') !LHS

      !! Restore external lambdas (from site-BB)
      IF(resX) CALL LambdaDIV(iGamma(sB) % m,  iLambda(xB) % m,  '6') !RHS
      IF(resY) CALL LambdaDIV(iGamma(sB) % m,  iLambda(yB) % m,  '4') !RHS
      IF(resY) CALL LambdaDIV(iGamma(sB) % m,  iLambda(yA) % m,  '3') !LHS

  CASE('E', 'W')

      !! Restore external lambdas (from site-AA)
      IF(resX) CALL LambdaDIV(iGamma(sA) % m,  iLambda(xB) % m,  '3') !LHS
      IF(resY) CALL LambdaDIV(iGamma(sA) % m,  iLambda(yA) % m,  '6') !RHS
      IF(resY) CALL LambdaDIV(iGamma(sA) % m,  iLambda(yB) % m,  '5') !LHS

      !! Restore external lambdas (from site-BB)
      IF(resX) CALL LambdaDIV(iGamma(sB) % m,  iLambda(xB) % m,  '4') !RHS
      IF(resY) CALL LambdaDIV(iGamma(sB) % m,  iLambda(yB) % m,  '6') !RHS
      IF(resY) CALL LambdaDIV(iGamma(sB) % m,  iLambda(yA) % m,  '5') !LHS

  CASE('SYM-AA')

      !! Restore lambdas after symmetric absorption (on site-AA)
      CALL LambdaDIV(iGamma(sA) % m,  SQRT(iLambda(xA) % m),  '4') 
      CALL LambdaDIV(iGamma(sA) % m,  SQRT(iLambda(xB) % m),  '3') 
      CALL LambdaDIV(iGamma(sA) % m,  SQRT(iLambda(yA) % m),  '6') 
      CALL LambdaDIV(iGamma(sA) % m,  SQRT(iLambda(yB) % m),  '5')

  CASE('SYM-BB')

      !! Restore lambdas after symmetric absorption (on site-BB)
      CALL LambdaDIV(iGamma(sB) % m,  SQRT(iLambda(xB) % m),  '4') 
      CALL LambdaDIV(iGamma(sB) % m,  SQRT(iLambda(xA) % m),  '3') 
      CALL LambdaDIV(iGamma(sB) % m,  SQRT(iLambda(yB) % m),  '6')
      CALL LambdaDIV(iGamma(sB) % m,  SQRT(iLambda(yA) % m),  '5') 

  CASE DEFAULT
      CALL invalid_flag("restore_ipepo_lambdas -- invalid DIR ", DIR)
  end SELECT

 END SUBROUTINE restore_ipepo_lambdas

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE project_ipepo
