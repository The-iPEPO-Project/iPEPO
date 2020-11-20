MODULE project_impo

 USE utility
 USE definitions_mps_mpo
 USE mps_mpo_utility

 IMPLICIT NONE

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (1) PROJECT iMPO TO iMPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Project iMPO onto iMPS: fuse N-leg with either W or E legs
 SUBROUTINE project_impo_to_imps(imps, iGmpo, iLmpo, SITEORD)

  TYPE(block_mps),  ALLOCATABLE, INTENT(INOUT) :: imps(:)            
  TYPE(block_mpo),               INTENT(IN)    :: iGmpo(:) 
  TYPE(block_lambda),            INTENT(IN)    :: iLmpo(:)            
  CHARACTER(LEN=*),              INTENT(IN)    :: SITEORD

  !! Local iMPO copy
  TYPE(block_mpo), ALLOCATABLE :: impo(:)

  !! iMPO sites
  INTEGER :: sA, sB

  !! Get iMPO sites
  CALL get_impo_sites(sA, sB, SITEORD)

  !! Absorb iMPO lambdas
  CALL absorb_impo_lambdas(impo, iGmpo, iLmpo, SITEORD, SYM='ALL')
  !CALL copy_mpo_block(impo, iGmpo)

  !! Reshape MPO gammas into MPS gammas
  CALL allocate_empty_mps_block(imps, 2)
  CALL allocate_mps_site(imps(1),  RESHAPE_3D(impo(sA) % m, '1,23,4'))
  CALL allocate_mps_site(imps(2),  RESHAPE_3D(impo(sB) % m, '1,3,24'))

 END SUBROUTINE project_impo_to_imps




 !!! Restore iMPO from iMPS projection: split external legs into N-leg and either W or E leg !!!
 SUBROUTINE unproject_impo_from_imps(iGmpo, iLmpo, iGmps, iLmps, SITEORD)

  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT) :: iGmpo(:) 
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLmpo(:)           
  TYPE(block_mps),                 INTENT(IN)    :: iGmps(:)
  TYPE(block_lambda),              INTENT(IN)    :: iLmps           
  CHARACTER(LEN=*),                INTENT(IN)    :: SITEORD

  !! iMPO sites
  INTEGER :: sA, sB

  !! Get iMPO sites
  CALL get_impo_sites(sA, sB, SITEORD)

  !! Reshape iMPS projection into iMPO
  CALL allocate_mpo_site(iGmpo(sA),  RESHAPE_4D(iGmps(1) % m, '1,23,4', (/ iGmpo(sA) % Ndim, iGmpo(sA) % Wdim /)))
  CALL allocate_mpo_site(iGmpo(sB),  RESHAPE_4D(iGmps(2) % m, '1,3,24', (/ iGmpo(sB) % Ndim, iGmpo(sB) % Edim /)))

  !! Copy iLmps to iLmpo(sA)
  CALL allocate_lambda_site(iLmpo(sA), iLmps % m)

  !! Restore external iMPO lambdas
  CALL LambdaDIV(iGmpo(sA) % m, iLmpo(sB) % m,  '3')
  CALL LambdaDIV(iGmpo(sB) % m, iLmpo(sB) % m,  '4')

 END SUBROUTINE unproject_impo_from_imps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (2) CALCULATE FACTORIZED iMPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Project iMPO to iMPS, and factorize it as [imps --> imps * X-site] !!!
 SUBROUTINE calc_factorized_impo(imps, Xa, Yb, iGmpo, iLmpo, SITEORD)

  TYPE(block_mps),  ALLOCATABLE, INTENT(INOUT) :: imps(:)            
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: Xa(:,:), Yb(:,:)
  TYPE(block_mpo),               INTENT(IN)    :: iGmpo(:)            
  TYPE(block_lambda),            INTENT(IN)    :: iLmpo(:)
  CHARACTER(LEN=*),              INTENT(IN)    :: SITEORD

  !! MPO block with lambdas absorbed
  TYPE(block_mpo),  ALLOCATABLE :: impo(:)

  !! Temporary iMPO projection onto iMPS
  TYPE(block_mps),  ALLOCATABLE :: imps_tmp(:)  
  COMPLEX(KIND=DP), ALLOCATABLE :: mps11(:,:,:), mps22(:,:,:)  

  !! Project iMPO onto iMPS
  CALL project_impo_to_imps(imps_tmp, iGmpo, iLmpo, SITEORD)

  !! Factorize each site of iMPS into impsA -> Xa * impsA, impsB -> impsB * Yb
  CALL allocate_empty_mps_block(imps, 2)

  CALL LHS_factorize_Tens3D(mps11,  Xa,  imps_tmp(1) % m,  '2,13')
  CALL RHS_factorize_Tens3D(mps22,  Yb,  imps_tmp(2) % m,  '12,3')

  CALL allocate_mps_site(imps(1),  mps11)
  CALL allocate_mps_site(imps(2),  mps22)

 END SUBROUTINE calc_factorized_impo



 !!! Multiply factorized iMPS-SITE and X-SITE as [imps * X -> imps], and restore iMPO from iMPS projection
 SUBROUTINE restore_unfactorized_impo(iGmpo, iLmpo, iGmps, iLmps, Xa, Yb, SITEORD)

  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT) :: iGmpo(:)            
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLmpo(:)
  TYPE(block_mps),                 INTENT(IN)    :: iGmps(:) 
  TYPE(block_lambda),              INTENT(IN)    :: iLmps          
  COMPLEX(KIND=DP),                INTENT(IN)    :: Xa(:,:), Yb(:,:)
  CHARACTER(LEN=*),                INTENT(IN)    :: SITEORD

  !! Temporary iMPO projection onto iMPS
  TYPE(block_mps), ALLOCATABLE :: iGmpsTmp(:)  

  !! Multiply factorized iMPS sites
  CALL allocate_empty_mps_block(iGmpsTmp, 2)
  CALL allocate_mps_site(iGmpsTmp(1), TENSMUL(iGmps(1) % m, Xa, '22'))
  CALL allocate_mps_site(iGmpsTmp(2), TENSMUL(iGmps(2) % m, Yb, '31'))

  !! Restore iMPO from iMPS projection
  CALL unproject_impo_from_imps(iGmpo, iLmpo, iGmpsTmp, iLmps, SITEORD)

 END SUBROUTINE restore_unfactorized_impo

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE project_impo
