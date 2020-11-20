MODULE project_ipeps

 USE utility
 USE definitions_mps_mpo
 USE definitions_peps_pepo
 USE mps_mpo_utility

 IMPLICIT NONE

 interface project_ipeps_to_imps
    module procedure project_ipeps_to_imps_GamLam
    module procedure project_ipeps_to_imps_PLAIN
 end interface project_ipeps_to_imps

 interface unproject_ipeps_from_imps
    module procedure unproject_ipeps_from_imps_GamLam
    module procedure unproject_ipeps_from_imps_PLAIN
 end interface unproject_ipeps_from_imps

 interface calc_factorized_ipeps
    module procedure calc_factorized_ipeps_GamLam
    module procedure calc_factorized_ipeps_PLAIN
 end interface calc_factorized_ipeps

 interface restore_unfactorized_ipeps
    module procedure restore_unfactorized_ipeps_GamLam
    module procedure restore_unfactorized_ipeps_PLAIN
 end interface restore_unfactorized_ipeps

 private !hides all items not listed in public statement 
 public project_ipeps_to_imps, unproject_ipeps_from_imps
 public calc_factorized_ipeps, restore_unfactorized_ipeps
 public absorb_ipeps_lambdas,  restore_ipeps_lambdas

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (1) PROJECT iPEPS TO iMPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Project ipeps onto imps: fuse Y-axis legs with X-axis legs in the XY plane !!!
 SUBROUTINE project_ipeps_to_imps_GamLam(iGmps, iLmps, iGpeps, iLpeps, DIR)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGmps(:)            
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLmps(:)
  TYPE(block_peps),                INTENT(IN)    :: iGpeps(:)            
  TYPE(block_lambda),              INTENT(IN)    :: iLpeps(:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! Local ipeps copy
  TYPE(block_peps), ALLOCATABLE :: ipeps(:)

  !! iPEPS sites && reshaping flags
  CHARACTER(LEN=16) :: legsA, legsB
  INTEGER           :: sA, sB
  INTEGER           :: xA, xB, yA, yB

  !! Get iPEPS sites
  CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, DIR)

  !! Set reshaping flag
  SELECT CASE(DIR)
  CASE('S', 'N')
       legsA = '1,234,5'; legsB = '1,4,235'
  CASE('E', 'W') 
       legsA = '1,245,3'; legsB = '1,2,345'
  CASE DEFAULT
       CALL invalid_flag("project_ipeps_to_imps_GamLam -- invalid DIR ", DIR)
  end SELECT

  !! Absorb lambdas -- external lambdas lying on both Y-axis and X-axis (keep only one internal lambda)
  CALL absorb_ipeps_lambdas(ipeps, iGpeps, iLpeps, DIR, absX=.True., absY=.True.) 

  !! Reshape PEPS gammas into MPS gammas
  CALL allocate_empty_mps_block(iGmps, 2)
  CALL allocate_mps_site(iGmps(1),  RESHAPE_3D(ipeps(sA) % m, legsA))
  CALL allocate_mps_site(iGmps(2),  RESHAPE_3D(ipeps(sB) % m, legsB))

  !! Write PEPS central lambda to MPS lambda
  CALL allocate_empty_lambda_block(iLmps, 1)
  CALL allocate_lambda_site(iLmps(1), iLpeps(xA) % m)

 END SUBROUTINE project_ipeps_to_imps_GamLam




 !!! Project ipeps onto imps: fuse Y-axis legs with X-axis legs in the XY plane !!!
 SUBROUTINE NEW_project_ipeps_to_imps_GamLam(iGmps, iLmps, iGpeps, iLpeps, DIR)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGmps(:)            
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLmps(:)
  TYPE(block_peps),                INTENT(IN)    :: iGpeps(:)            
  TYPE(block_lambda),              INTENT(IN)    :: iLpeps(:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! Local ipeps copy
  TYPE(block_peps), ALLOCATABLE :: ipeps(:)

  !! iPEPS sites && reshaping flags
  CHARACTER(LEN=16) :: legsA, legsB
  INTEGER           :: sA, sB
  INTEGER           :: xA, xB, yA, yB

  !! Get iPEPS sites
  CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, DIR)

  !! Absorb lambdas -- external lambdas lying on both Y-axis and X-axis (keep only one internal lambda)
  CALL absorb_ipeps_lambdas(ipeps, iGpeps, iLpeps, DIR, absX=.True., absY=.True.) 

  !! If bond is vertical, rotate to make it horizontal
  IF((DIR .EQ. 'W') .OR. (DIR .EQ. 'E')) THEN
      CALL allocate_peps_site(ipeps(sA),  TensROTATE(ipeps(sA) % m, 'CW +PI/2'))
      CALL allocate_peps_site(ipeps(sB),  TensROTATE(ipeps(sB) % m, 'CW +PI/2'))
  end IF

  !! Reshape PEPS gammas into MPS gammas
  CALL allocate_empty_mps_block(iGmps, 2)
  CALL allocate_mps_site(iGmps(1),  RESHAPE_3D(ipeps(sA) % m, '1,234,5'))
  CALL allocate_mps_site(iGmps(2),  RESHAPE_3D(ipeps(sB) % m, '1,4,235'))

  !! Write PEPS central lambda to MPS lambda
  CALL allocate_empty_lambda_block(iLmps, 1)
  CALL allocate_lambda_site(iLmps(1), iLpeps(xA) % m)

 END SUBROUTINE NEW_project_ipeps_to_imps_GamLam






 !!! Project ipeps onto imps: fuse Y-axis legs with X-axis legs in the XY plane !!!
 SUBROUTINE project_ipeps_to_imps_PLAIN(imps, ipeps, DIR)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: imps(:)            
  TYPE(block_peps),                INTENT(IN)    :: ipeps(:)            
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! iPEPS sites && reshaping flags
  CHARACTER(LEN=16) :: legsA, legsB
  INTEGER           :: sA, sB
  INTEGER           :: xA, xB, yA, yB

  !! Get iPEPS sites
  CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, DIR)

  !! Set reshaping flag
  SELECT CASE(DIR)
  CASE('S', 'N')
       legsA = '1,234,5'; legsB = '1,4,235'
  CASE('E', 'W') 
       legsA = '1,245,3'; legsB = '1,2,345'
  CASE DEFAULT
       CALL invalid_flag("project_ipeps_to_imps_PLAIN -- invalid DIR ", DIR)
  end SELECT

  !! Reshape PEPS into MPS
  CALL allocate_empty_mps_block(imps, 2)
  CALL allocate_mps_site(imps(1),  RESHAPE_3D(ipeps(sA) % m, legsA))
  CALL allocate_mps_site(imps(2),  RESHAPE_3D(ipeps(sB) % m, legsB))

 END SUBROUTINE project_ipeps_to_imps_PLAIN

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (2) RESTORE iPEPS FROM iMPS PROJECTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Restore ipeps from imps projection: split external legs into Y-axis && X-axis legs in the XY plane !!!
 SUBROUTINE unproject_ipeps_from_imps_GamLam(iGpeps, iLpeps, iGmps, iLmps, DIR)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGpeps(:)            
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLpeps(:)
  TYPE(block_mps),                 INTENT(IN)    :: iGmps(:)            
  TYPE(block_lambda),              INTENT(IN)    :: iLmps(:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! iPEPS sites && reshaping flags
  CHARACTER(LEN=16) :: legsA, legsB
  INTEGER           :: dimsA(3), dimsB(3)
  INTEGER           :: sA, sB
  INTEGER           :: xA, xB, yA, yB

  !! Get iPEPS sites 
  !! IMPORTANT: do not allocate new empty iPEPS blocks when restoring iPEPS from iMPS projection, 
  !! since old lambdas will be needed by restore_ipeps_lambdas()
  CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, DIR)

  !! Set reshaping flag && dims
  SELECT CASE(DIR)
  CASE('S', 'N')
       legsA = '1,234,5'; dimsA = (/ iGpeps(sA) % Sdim, iGpeps(sA) % Ndim, iGpeps(sA) % Wdim /)
       legsB = '1,4,235'; dimsB = (/ iGpeps(sB) % Sdim, iGpeps(sB) % Ndim, iGpeps(sB) % Edim /)                     
  CASE('E', 'W') 
       legsA = '1,245,3'; dimsA = (/ iGpeps(sA) % Sdim, iGpeps(sA) % Wdim, iGpeps(sA) % Edim /)
       legsB = '1,2,345'; dimsB = (/ iGpeps(sB) % Ndim, iGpeps(sB) % Wdim, iGpeps(sB) % Edim /)
  CASE DEFAULT
       CALL invalid_flag("unproject_ipeps_from_imps_GamLam -- invalid DIR ", DIR)
  end SELECT

  !! Reshape iMPS projection into iPEPS 
  CALL allocate_peps_site(iGpeps(sA),    RESHAPE_5D(iGmps(1) % m, legsA, dimsA))
  CALL allocate_peps_site(iGpeps(sB),    RESHAPE_5D(iGmps(2) % m, legsB, dimsB))
  CALL allocate_lambda_site(iLpeps(xA),  iLmps(1) % m) 

  !! Restore lambdas -- i.e. external lambdas lying on both Y-axis and X-axis
  CALL restore_ipeps_lambdas(iGpeps, iLpeps, DIR, resX=.True., resY=.True.)

 END SUBROUTINE unproject_ipeps_from_imps_GamLam






 !!! Restore ipeps from imps projection: split external legs into Y-axis && X-axis legs in the XY plane !!!
 SUBROUTINE NEW_unproject_ipeps_from_imps_GamLam(iGpeps, iLpeps, iGmps, iLmps, DIR)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGpeps(:)            
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLpeps(:)
  TYPE(block_mps),                 INTENT(IN)    :: iGmps(:)            
  TYPE(block_lambda),              INTENT(IN)    :: iLmps(:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! iPEPS sites && reshaping flags
  CHARACTER(LEN=16) :: legsA, legsB
  INTEGER           :: dimsA(3), dimsB(3)
  INTEGER           :: sA, sB
  INTEGER           :: xA, xB, yA, yB

  !! Get iPEPS sites 
  !! IMPORTANT: do not allocate new empty iPEPS blocks when restoring iPEPS from iMPS projection, 
  !! since old lambdas will be needed by restore_ipeps_lambdas()
  CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, DIR)

  !! Set reshaping flag && dims
  SELECT CASE(DIR)
  CASE('S', 'N')
       dimsA = (/ iGpeps(sA) % Sdim, iGpeps(sA) % Ndim, iGpeps(sA) % Wdim /)
       dimsB = (/ iGpeps(sB) % Sdim, iGpeps(sB) % Ndim, iGpeps(sB) % Edim /)                     
  CASE('W', 'E') 
       dimsA = (/ iGpeps(sA) % Edim, iGpeps(sA) % Wdim, iGpeps(sA) % Sdim /)
       dimsB = (/ iGpeps(sB) % Edim, iGpeps(sB) % Wdim, iGpeps(sB) % Ndim /)
  CASE DEFAULT
       CALL invalid_flag("unproject_ipeps_from_imps_GamLam -- invalid DIR ", DIR)
  end SELECT

  !! Reshape iMPS projection into iPEPS 
  CALL allocate_peps_site(iGpeps(sA),    RESHAPE_5D(iGmps(1) % m, '1,234,5', dimsA))
  CALL allocate_peps_site(iGpeps(sB),    RESHAPE_5D(iGmps(2) % m, '1,4,235', dimsB))
  CALL allocate_lambda_site(iLpeps(xA),  iLmps(1) % m) 

  !! If bond is vertical, we had made it horizontal --> rotate back to make it vertical again
  IF((DIR .EQ. 'W') .OR. (DIR .EQ. 'E')) THEN
      CALL allocate_peps_site(iGpeps(sA), TensROTATE(iGpeps(sA) % m, 'CCW -PI/2'))
      CALL allocate_peps_site(iGpeps(sB), TensROTATE(iGpeps(sB) % m, 'CCW -PI/2'))
  end IF

  !! Restore lambdas -- i.e. external lambdas lying on both Y-axis and X-axis
  CALL restore_ipeps_lambdas(iGpeps, iLpeps, DIR, resX=.True., resY=.True.)

 END SUBROUTINE NEW_unproject_ipeps_from_imps_GamLam







 !!! Restore ipeps from imps projection: split external legs into Y-axis && X-axis legs in the XY plane !!!
 SUBROUTINE unproject_ipeps_from_imps_PLAIN(ipeps, imps, DIR)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: ipeps(:)            
  TYPE(block_mps),                 INTENT(IN)    :: imps(:)            
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! iPEPS sites && reshaping flags
  CHARACTER(LEN=16) :: legsA, legsB
  INTEGER           :: dimsA(3), dimsB(3)
  INTEGER           :: sA, sB
  INTEGER           :: xA, xB, yA, yB

  !! Get iPEPS sites 
  CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, DIR)

  !! Set reshaping flag && dims
  SELECT CASE(DIR)
  CASE('S', 'N')
       legsA = '1,234,5'; dimsA = (/ ipeps(sA) % Sdim, ipeps(sA) % Ndim, ipeps(sA) % Wdim /)
       legsB = '1,4,235'; dimsB = (/ ipeps(sB) % Sdim, ipeps(sB) % Ndim, ipeps(sB) % Edim /)                     
  CASE('E', 'W') 
       legsA = '1,245,3'; dimsA = (/ ipeps(sA) % Sdim, ipeps(sA) % Wdim, ipeps(sA) % Edim /)
       legsB = '1,2,345'; dimsB = (/ ipeps(sB) % Ndim, ipeps(sB) % Wdim, ipeps(sB) % Edim /)
  CASE DEFAULT
       CALL invalid_flag("unproject_ipeps_from_imps_PLAIN -- invalid DIR ", DIR)
  end SELECT

  !! Reshape iMPS projection into iPEPS
  CALL allocate_peps_site(ipeps(sA),  RESHAPE_5D(imps(1) % m, legsA, dimsA))
  CALL allocate_peps_site(ipeps(sB),  RESHAPE_5D(imps(2) % m, legsB, dimsB))

 END SUBROUTINE unproject_ipeps_from_imps_PLAIN

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (3) CALCULATE FACTORIZED iPEPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Project ipeps to imps, and factorize it as [imps --> imps * X-site] !!!
 SUBROUTINE calc_factorized_ipeps_GamLam(iGmps, iLmps, Xa, Yb, iGpeps, iLpeps, DIR)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGmps(:)            
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLmps(:)
  COMPLEX(KIND=DP),   ALLOCATABLE, INTENT(INOUT) :: Xa(:,:), Yb(:,:)
  TYPE(block_peps),                INTENT(IN)    :: iGpeps(:)            
  TYPE(block_lambda),              INTENT(IN)    :: iLpeps(:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! Temporary iPEPS projection onto iMPS && temp factorized MPS tensors
  TYPE(block_mps),  ALLOCATABLE :: iGmpsTmp(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: mps11(:,:,:), mps22(:,:,:)  

  !! (1) Project iPEPS onto iMPS
  CALL project_ipeps_to_imps(iGmpsTmp, iLmps, iGpeps, iLpeps, DIR)

  !! (2) Factorize each site of iMPS into impsA -> Xa * impsA, impsB -> impsB * Yb
  CALL allocate_empty_mps_block(iGmps, 2)

  CALL LHS_factorize_Tens3D(mps11,  Xa,  iGmpsTmp(1) % m,  '2,13')
  CALL allocate_mps_site(iGmps(1),  mps11)

  CALL RHS_factorize_Tens3D(mps22,  Yb,  iGmpsTmp(2) % m,  '12,3')
  CALL allocate_mps_site(iGmps(2),  mps22)

 END SUBROUTINE calc_factorized_ipeps_GamLam





 !!! Project ipeps to imps, and factorize it as [imps --> imps * X-site] !!!
 SUBROUTINE calc_factorized_ipeps_PLAIN(imps, Xa, Yb, ipeps, DIR)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: imps(:)            
  COMPLEX(KIND=DP),   ALLOCATABLE, INTENT(INOUT) :: Xa(:,:), Yb(:,:)
  TYPE(block_peps),                INTENT(IN)    :: ipeps(:)            
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! Temporary iPEPS projection onto iMPS && temp factorized MPS tensors
  TYPE(block_mps),  ALLOCATABLE :: impsTmp(:) 
  COMPLEX(KIND=DP), ALLOCATABLE :: mps11(:,:,:), mps22(:,:,:)  

  !! (1) Project iPEPS onto iMPS
  CALL project_ipeps_to_imps(impsTmp, ipeps, DIR)

  !! (2) Factorize each site of iMPS into impsA -> Xa * impsA, impsB -> impsB * Yb
  CALL allocate_empty_mps_block(imps, 2)

  CALL LHS_factorize_Tens3D(mps11,  Xa,  impsTmp(1) % m,  '2,13')
  CALL allocate_mps_site(imps(1),  mps11)

  CALL RHS_factorize_Tens3D(mps22,  Yb,  impsTmp(2) % m,  '12,3')
  CALL allocate_mps_site(imps(2),  mps22)

 END SUBROUTINE calc_factorized_ipeps_PLAIN

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (4) RESTORE UNFACTORIZED iPEPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Multiply factorized iMPS-SITE and X-SITE as [imps * X -> imps], 
 !!! and restore iPEPS from iMPS projection
 SUBROUTINE restore_unfactorized_ipeps_GamLam(iGpeps, iLpeps, iGmps, iLmps, Xa, Yb, DIR)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGpeps(:)            
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLpeps(:)
  TYPE(block_mps),                 INTENT(IN)    :: iGmps(:)            
  TYPE(block_lambda),              INTENT(IN)    :: iLmps(:)
  COMPLEX(KIND=DP),                INTENT(IN)    :: Xa(:,:), Yb(:,:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! Temporary iPEPS projection onto iMPS
  TYPE(block_mps), ALLOCATABLE :: iGmpsTmp(:)  

  !! (1) Multiply factorized iPEPS sites
  CALL allocate_empty_mps_block(iGmpsTmp, 2)
  CALL allocate_mps_site(iGmpsTmp(1), TENSMUL(iGmps(1) % m, Xa, '22'))
  CALL allocate_mps_site(iGmpsTmp(2), TENSMUL(iGmps(2) % m, Yb, '31'))

  !! (2) Restore iPEPS from iMPS projection
  CALL unproject_ipeps_from_imps(iGpeps, iLpeps, iGmpsTmp, iLmps, DIR)

 END SUBROUTINE restore_unfactorized_ipeps_GamLam




 !!! Multiply factorized iMPS-SITE and X-SITE as [imps * X -> imps], 
 !!! and restore iPEPS from iMPS projection
 SUBROUTINE restore_unfactorized_ipeps_PLAIN(ipeps, imps, Xa, Yb, DIR)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: ipeps(:)            
  TYPE(block_mps),                 INTENT(IN)    :: imps(:)            
  COMPLEX(KIND=DP),                INTENT(IN)    :: Xa(:,:), Yb(:,:)
  CHARACTER(LEN=*),                INTENT(IN)    :: DIR

  !! Temporary iPEPS projection onto iMPS
  TYPE(block_mps), ALLOCATABLE :: impsTmp(:)  

  !! (1) Multiply factorized iPEPS sites
  CALL allocate_empty_mps_block(impsTmp, 2)
  CALL allocate_mps_site(impsTmp(1), TENSMUL(imps(1) % m, Xa, '22'))
  CALL allocate_mps_site(impsTmp(2), TENSMUL(imps(2) % m, Yb, '31'))

  !! (2) Restore iPEPS from iMPS projection
  CALL unproject_ipeps_from_imps(ipeps, impsTmp, DIR)

 END SUBROUTINE restore_unfactorized_ipeps_PLAIN

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (5) ABSORB/RESTORE iPEPS LAMBDAS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Absorb lambdas on Y-axis, external lambda on X-axis, or absorb all lambdas in symmetric fashion !!!
 SUBROUTINE absorb_ipeps_lambdas(ipeps, iGamma, iLambda, DIR, absX, absY)

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT)           :: ipeps(:)            !iPEPS with absorbed lambdas
  TYPE(block_peps),                INTENT(IN),   OPTIONAL  :: iGamma(:)           !iPEPS gamma
  TYPE(block_lambda),              INTENT(IN)              :: iLambda(:)          !iPEPS lambda
  CHARACTER(LEN=*),                INTENT(IN),   OPTIONAL  :: DIR                 !central bond direction/orientation
  LOGICAL,                         INTENT(IN),   OPTIONAL  :: absX, absY          !whether to absorb lambdas lying on 
                                                                                  !tangential && normal axes
  !! iPEPS sites
  INTEGER :: sA, sB
  INTEGER :: xA, xB, yA, yB

  !! Copy iGamma to iPeps block (if iGamma not present, use iPEPS as input Gamma instead)
  IF(PRESENT(iGamma)) CALL copy_peps_block(ipeps, iGamma) 

  !! Absorb lambdas
  IF(PRESENT(absX) .AND. PRESENT(absY) .AND. PRESENT(DIR)) THEN 

     !! Get iPEPS sites (LHS = S,W), (RHS = N,E)
     CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, DIR)

     !! Absorb external lambdas ONLY
     SELECT CASE(DIR)
     CASE('S', 'N')

         !! Absorb external lambdas on Y-axis (into siteA)
         IF(absX) CALL LambdaMUL(ipeps(sA) % m,  iLambda(xB) % m,  '4') !LHS
         IF(absY) CALL LambdaMUL(ipeps(sA) % m,  iLambda(yA) % m,  '3') !RHS 
         IF(absY) CALL LambdaMUL(ipeps(sA) % m,  iLambda(yB) % m,  '2') !LHS 

         !! Absorb external lambdas on Y-axis (into siteB)
         IF(absX) CALL LambdaMUL(ipeps(sB) % m,  iLambda(xB) % m,  '5') !RHS
         IF(absY) CALL LambdaMUL(ipeps(sB) % m,  iLambda(yB) % m,  '3') !RHS
         IF(absY) CALL LambdaMUL(ipeps(sB) % m,  iLambda(yA) % m,  '2') !LHS

     CASE('E', 'W')

         !! Absorb external lambdas on Y-axis (into siteA)
         IF(absX) CALL LambdaMUL(ipeps(sA) % m,  iLambda(xB) % m,  '2') !LHS
         IF(absY) CALL LambdaMUL(ipeps(sA) % m,  iLambda(yA) % m,  '5') !RHS
         IF(absY) CALL LambdaMUL(ipeps(sA) % m,  iLambda(yB) % m,  '4') !LHS

         !! Absorb external lambdas on Y-axis (into siteB)
         IF(absX) CALL LambdaMUL(ipeps(sB) % m,  iLambda(xB) % m,  '3') !RHS
         IF(absY) CALL LambdaMUL(ipeps(sB) % m,  iLambda(yB) % m,  '5') !RHS
         IF(absY) CALL LambdaMUL(ipeps(sB) % m,  iLambda(yA) % m,  '4') !LHS

     CASE DEFAULT
         CALL invalid_flag("absorb_ipeps_lambdas -- invalid DIR ", DIR)
     end SELECT

  ELSE

     !! Get iPEPS sites for the default configuration (centred on 'S' bond)
     CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, 'S')

     !! Absorb all lambdas in symmetric fashion -- on SN bonds 
     CALL LambdaMUL(ipeps(sA) % m,  iLambda(xA) % m,  ipeps(sB) % m,  '23')
     CALL LambdaMUL(ipeps(sB) % m,  iLambda(xB) % m,  ipeps(sA) % m,  '23')
 
     !! Absorb all lambdas in symmetric fashion -- on EW bonds
     CALL LambdaMUL(ipeps(sA) % m,  iLambda(yA) % m,  ipeps(sB) % m,  '45')
     CALL LambdaMUL(ipeps(sB) % m,  iLambda(yB) % m,  ipeps(sA) % m,  '45')
  end IF 

 END SUBROUTINE absorb_ipeps_lambdas






 !!! Restore external iPEPS lambdas: convert from nolambda form to Gamma-Lambda form !!!
 SUBROUTINE restore_ipeps_lambdas(iGamma, iLambda, DIR, resX, resY)

  TYPE(block_peps),   INTENT(INOUT) :: iGamma(:)             !INPUT iPEPS: with all external Lambdas absorbed
  TYPE(block_lambda), INTENT(IN)    :: iLambda(:)            !INPUT iLambdas: input old external Lambdas so we can restore them
  CHARACTER(LEN=*),   INTENT(IN)    :: DIR                   !central bond direction/orientation
  LOGICAL,            INTENT(IN)    :: resX, resY            !external lambdas to be restored

  !! iPEPS sites
  INTEGER :: sA, sB
  INTEGER :: xA, xB, yA, yB

  !! Get iPEPS sites (LEFT = S,W --- RIGHT = N,E)
  CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, DIR)

  !! Restore external lambdas
  SELECT CASE(DIR)
  CASE('S', 'N')

      !! Restore external lambdas (from site-AA)
      IF(resX) CALL LambdaDIV(iGamma(sA) % m,  iLambda(xB) % m,  '4') !LHS
      IF(resY) CALL LambdaDIV(iGamma(sA) % m,  iLambda(yA) % m,  '3') !RHS
      IF(resY) CALL LambdaDIV(iGamma(sA) % m,  iLambda(yB) % m,  '2') !LHS

      !! Restore external lambdas (from site-BB)
      IF(resX) CALL LambdaDIV(iGamma(sB) % m,  iLambda(xB) % m,  '5') !RHS
      IF(resY) CALL LambdaDIV(iGamma(sB) % m,  iLambda(yB) % m,  '3') !RHS
      IF(resY) CALL LambdaDIV(iGamma(sB) % m,  iLambda(yA) % m,  '2') !LHS

  CASE('E', 'W')

      !! Restore external lambdas (from site-AA)
      IF(resX) CALL LambdaDIV(iGamma(sA) % m,  iLambda(xB) % m,  '2') !LHS
      IF(resY) CALL LambdaDIV(iGamma(sA) % m,  iLambda(yA) % m,  '5') !RHS
      IF(resY) CALL LambdaDIV(iGamma(sA) % m,  iLambda(yB) % m,  '4') !LHS

      !! Restore external lambdas (from site-BB)
      IF(resX) CALL LambdaDIV(iGamma(sB) % m,  iLambda(xB) % m,  '3') !RHS
      IF(resY) CALL LambdaDIV(iGamma(sB) % m,  iLambda(yB) % m,  '5') !RHS
      IF(resY) CALL LambdaDIV(iGamma(sB) % m,  iLambda(yA) % m,  '4') !LHS

  CASE('SYM-AA')

      !! Restore lambdas after symmetric absorption (on site-AA)
      CALL LambdaDIV(iGamma(sA) % m,  SQRT(iLambda(xA) % m),  '3') 
      CALL LambdaDIV(iGamma(sA) % m,  SQRT(iLambda(xB) % m),  '2') 
      CALL LambdaDIV(iGamma(sA) % m,  SQRT(iLambda(yA) % m),  '5') 
      CALL LambdaDIV(iGamma(sA) % m,  SQRT(iLambda(yB) % m),  '4')

  CASE('SYM-BB')

      !! Restore lambdas after symmetric absorption (on site-BB)
      CALL LambdaDIV(iGamma(sB) % m,  SQRT(iLambda(xB) % m),  '3') 
      CALL LambdaDIV(iGamma(sB) % m,  SQRT(iLambda(xA) % m),  '2') 
      CALL LambdaDIV(iGamma(sB) % m,  SQRT(iLambda(yB) % m),  '5')
      CALL LambdaDIV(iGamma(sB) % m,  SQRT(iLambda(yA) % m),  '4') 

  CASE DEFAULT
      CALL invalid_flag("restore_ipeps_lambdas -- invalid DIR ", DIR)
  end SELECT

 END SUBROUTINE restore_ipeps_lambdas

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE project_ipeps
