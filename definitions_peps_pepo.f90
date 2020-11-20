MODULE definitions_peps_pepo

 USE utility
 USE definitions_mps_mpo

 IMPLICIT NONE

 !! Define the PEPS type
 TYPE block_peps
   INTEGER :: LocalDim, Sdim, Ndim, Wdim, Edim
   COMPLEX(KIND=DP), ALLOCATABLE :: m(:,:,:,:,:) 
 END TYPE block_peps

 !! Define the PEPO type
 TYPE block_pepo
   INTEGER :: dnDim, upDim, Sdim, Ndim, Wdim, Edim
   COMPLEX(KIND=DP), ALLOCATABLE :: m(:,:,:,:,:,:) 
 END TYPE block_pepo

 interface allocate_peps_site
    module procedure allocate_peps_site_u_input_tensor
    module procedure allocate_peps_site_u_vec_of_dims
 end interface allocate_peps_site

 interface allocate_pepo_site
    module procedure allocate_pepo_site_u_input_tensor
    module procedure allocate_pepo_site_u_vec_of_dims
 end interface allocate_pepo_site

 private allocate_peps_site_u_input_tensor, allocate_peps_site_u_vec_of_dims

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE/COPY PEPS BLOCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!! allocate empty peps block without any sites in it !!!!
 SUBROUTINE allocate_empty_peps_block(peps_block, N_sites)

  TYPE(block_peps), ALLOCATABLE, INTENT(INOUT) :: peps_block(:)
  INTEGER, INTENT(IN) :: N_sites

  CALL deallocate_peps_block(peps_block)

  !allocate new empty block (new sites to be appended)
  ALLOCATE(peps_block(N_sites))

 END SUBROUTINE allocate_empty_peps_block



 !!!!!!!!!!! deallocate mps block !!!!!!!!!!!!!!!
 SUBROUTINE deallocate_peps_block(peps_block)

  TYPE(block_peps), ALLOCATABLE, INTENT(INOUT) :: peps_block(:)

  !dims & indices
  INTEGER :: N_sites, site

  N_sites = SIZE(peps_block)

  IF(ALLOCATED(peps_block)) THEN
    !if previously allocated, deallocate all sites first
    DO site=1,N_sites
         IF(ALLOCATED(peps_block(site)%m)) DEALLOCATE(peps_block(site)%m)
    end DO
    !then deallocate the block itself
    DEALLOCATE(peps_block)
  end IF

 END SUBROUTINE deallocate_peps_block



 !!!!!!!!!!!!!!!!!!!!! copy mps block !!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE copy_peps_block(peps_copy, peps_block, hc_flag)

   TYPE(block_peps), ALLOCATABLE, INTENT(OUT) :: peps_copy(:)
   TYPE(block_peps), INTENT(IN) :: peps_block(:)
   CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: hc_flag

   !dims & indices
   INTEGER :: N_sites, site

   N_sites = SIZE(peps_block)

   !Allocate all sites of the vector
   CALL allocate_empty_peps_block(peps_copy, N_sites)

   IF(PRESENT(hc_flag)) THEN
      !create hermitian conjugate of evec_original if requested
      Do site=1,N_sites
         CALL allocate_peps_site(peps_copy(site), CONJG(peps_block(site)%m))
      end Do
   ELSE
      !else just a simple copy of evec_original
      Do site=1,N_sites
         CALL allocate_peps_site(peps_copy(site), peps_block(site)%m)
      end Do
   end IF

 END SUBROUTINE copy_peps_block




 !!!!!!!!!!!!!!!!!!!!!!!!! copy iPEPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE copy_ipeps(iGamma_copy, iLambda_copy, iGamma, iLambda, hc_flag)

  TYPE(block_peps), ALLOCATABLE, INTENT(INOUT) :: iGamma_copy(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda_copy(:)
  TYPE(block_peps), INTENT(IN) :: iGamma(:)
  TYPE(block_lambda), INTENT(IN) :: iLambda(:)
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: hc_flag

  !Copy iGamma block (HC or simple copy)
  IF(PRESENT(hc_flag)) THEN
     CALL copy_peps_block(iGamma_copy, iGamma, hc_flag)
  ELSE
     CALL copy_peps_block(iGamma_copy, iGamma)
  end IF

  !Copy iLambda block
  CALL copy_lambda_block(iLambda_copy, iLambda)

 END SUBROUTINE copy_ipeps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE/COPY PEPO BLOCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!! allocate empty PEPO block without any sites in it !!!!
 SUBROUTINE allocate_empty_pepo_block(pepo, N_sites)

  TYPE(block_pepo), ALLOCATABLE, INTENT(INOUT) :: pepo(:)
  INTEGER,                       INTENT(IN)    :: N_sites

  CALL deallocate_pepo_block(pepo)

  !allocate new empty block (new sites to be appended)
  ALLOCATE(pepo(N_sites))

 END SUBROUTINE allocate_empty_pepo_block



 !!!!!!!!!!! deallocate PEPO block !!!!!!!!!!!!!!!
 SUBROUTINE deallocate_pepo_block(pepo)

  TYPE(block_pepo), ALLOCATABLE, INTENT(INOUT) :: pepo(:)

  !! Dims & indices
  INTEGER :: N_sites, site

  N_sites = SIZE(pepo)

  IF(ALLOCATED(pepo)) THEN

    !! If previously allocated, deallocate all sites first
    DO site=1,N_sites
       IF(ALLOCATED(pepo(site) % m)) DEALLOCATE(pepo(site) % m)
    end DO

    !! Then deallocate the block itself
    DEALLOCATE(pepo)
  end IF

 END SUBROUTINE deallocate_pepo_block



 !!!!!!!!!!!!!!!!!!!!! copy PEPO block !!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE copy_pepo_block(pepo_copy, pepo)

   TYPE(block_pepo), ALLOCATABLE, INTENT(OUT) :: pepo_copy(:)
   TYPE(block_pepo),              INTENT(IN)  :: pepo(:)

   !! Dims & indices
   INTEGER :: N_sites, site

   N_sites = SIZE(pepo)

   !! Allocate all sites 
   CALL allocate_empty_pepo_block(pepo_copy, N_sites)

   !! Copy sites
   Do site=1,N_sites
      CALL allocate_pepo_site(pepo_copy(site),  pepo(site) % m)
   end Do

 END SUBROUTINE copy_pepo_block




 !!!!!!!!!!!!!!!!!!!!!!!!! copy iPEPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE copy_ipepo(iGamma_copy, iLambda_copy, iGamma, iLambda)

  TYPE(block_pepo),   ALLOCATABLE, INTENT(INOUT) :: iGamma_copy(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda_copy(:)
  TYPE(block_pepo),                INTENT(IN)    :: iGamma(:)
  TYPE(block_lambda),              INTENT(IN)    :: iLambda(:)

  !! Copy iGamma, iLambda blocks
  CALL copy_pepo_block(iGamma_copy, iGamma)
  CALL copy_lambda_block(iLambda_copy, iLambda)

 END SUBROUTINE copy_ipepo

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE PEPS SITE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE allocate_peps_site_u_input_tensor(peps_site, input_tensor)

  TYPE(block_peps), INTENT(INOUT) :: peps_site
  COMPLEX(KIND=DP), INTENT(IN)    :: input_tensor(:,:,:,:,:)

  !set dimensions of peps_site
  peps_site%LocalDim = SIZE(input_tensor,1)
  peps_site%Sdim = SIZE(input_tensor,2)
  peps_site%Ndim = SIZE(input_tensor,3)
  peps_site%Wdim = SIZE(input_tensor,4)
  peps_site%Edim = SIZE(input_tensor,5)

  !allocate mps_site
  IF(ALLOCATED(peps_site%m)) DEALLOCATE(peps_site%m)
  ALLOCATE(peps_site%m(peps_site%LocalDim, peps_site%Sdim, peps_site%Ndim, peps_site%Wdim, peps_site%Edim))

  !Initialize mps_site
  peps_site%m = input_tensor

 END SUBROUTINE allocate_peps_site_u_input_tensor



 SUBROUTINE allocate_peps_site_u_vec_of_dims(peps_site, dim_peps)

  TYPE(block_peps), INTENT(INOUT) :: peps_site
  INTEGER,          INTENT(IN)    :: dim_peps(5)

  !set dimensions of mps_site
  peps_site%LocalDim = dim_peps(1)
  peps_site%Sdim = dim_peps(2)
  peps_site%Ndim = dim_peps(3)
  peps_site%Wdim = dim_peps(4)
  peps_site%Edim = dim_peps(5)

  !allocate mps_site
  IF(ALLOCATED(peps_site%m)) DEALLOCATE(peps_site%m)
  ALLOCATE(peps_site%m(peps_site%LocalDim, peps_site%Sdim, peps_site%Ndim, peps_site%Wdim, peps_site%Edim))

  !Initialize mps_site
  peps_site%m(:,:,:,:,:) = (0d0,0d0)

 END SUBROUTINE allocate_peps_site_u_vec_of_dims

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE PEPO SITE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE allocate_pepo_site_u_input_tensor(pepo, tensIN)

  TYPE(block_pepo), INTENT(INOUT) :: pepo
  COMPLEX(KIND=DP), INTENT(IN)    :: tensIN(:,:,:,:,:,:)

  !! Dims of input tensor
  INTEGER :: dimP(6)

  !! Get dims of tensIN
  dimP = shape(tensIN)

  !! Set dimensions of pepo site
  pepo % dnDim = dimP(1)
  pepo % upDim = dimP(2)
  pepo % Sdim  = dimP(3)
  pepo % Ndim  = dimP(4)
  pepo % Wdim  = dimP(5)
  pepo % Edim  = dimP(6)

  !! Allocate pepo site
  IF(ALLOCATED(pepo % m)) DEALLOCATE(pepo % m)
  ALLOCATE(pepo % m(pepo % dnDim, pepo % upDim, pepo % Sdim, pepo % Ndim, pepo % Wdim, pepo % Edim))

  !! Initialize pepo site
  pepo % m = tensIN

 END SUBROUTINE allocate_pepo_site_u_input_tensor




 SUBROUTINE allocate_pepo_site_u_vec_of_dims(pepo, dimP)

  TYPE(block_pepo), INTENT(INOUT) :: pepo
  INTEGER,          INTENT(IN)    :: dimP(6)

  !! Set dimensions of pepo site
  pepo % dnDim = dimP(1)
  pepo % upDim = dimP(2)
  pepo % Sdim  = dimP(3)
  pepo % Ndim  = dimP(4)
  pepo % Wdim  = dimP(5)
  pepo % Edim  = dimP(6)

  !! Allocate pepo site
  IF(ALLOCATED(pepo % m)) DEALLOCATE(pepo % m)
  ALLOCATE(pepo % m(pepo % dnDim, pepo % upDim, pepo % Sdim, pepo % Ndim, pepo % Wdim, pepo % Edim))

  !! Initialize pepo site
  pepo % m(:,:,:,:,:,:) = (0d0, 0d0)

 END SUBROUTINE allocate_pepo_site_u_vec_of_dims

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Auxiliary routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Determine iPEPS sites 
 !!! -- args = {gamma sites s(A,B), lambda sites x(A,B) && y(A,B), bond orientation DIR} 
 SUBROUTINE get_ipeps_sites(sA, sB, xA, xB, yA, yB, DIR)

  INTEGER,          INTENT(OUT) :: sA, sB
  INTEGER,          INTENT(OUT) :: xA, xB, yA, yB
  CHARACTER(LEN=*), INTENT(IN)  :: DIR

  !! Determine Gamma && Lambda sites for a given siteFlag && bondFlag
  SELECT CASE(DIR)
  CASE('W')

        sA = 1; sB = 2
        xA = 1; xB = 2; yA = 3; yB = 4

  CASE('S')

        sA = 1; sB = 2
        xA = 3; xB = 4; yA = 1; yB = 2

  CASE('E')

        sA = 2; sB = 1
        xA = 2; xB = 1; yA = 4; yB = 3

  CASE('N')

        sA = 2; sB = 1
        xA = 4; xB = 3; yA = 2; yB = 1

  CASE DEFAULT
        CALL invalid_flag("get_ipeps_sites -- invalid DIR ", DIR)
  end SELECT

 END SUBROUTINE get_ipeps_sites



 !!! Determine iPEPO sites 
 !!! -- same as get_ipeps_sites
 SUBROUTINE get_ipepo_sites(sA, sB, xA, xB, yA, yB, DIR)

  INTEGER,          INTENT(OUT) :: sA, sB
  INTEGER,          INTENT(OUT) :: xA, xB, yA, yB
  CHARACTER(LEN=*), INTENT(IN)  :: DIR

  CALL get_ipeps_sites(sA, sB, xA, xB, yA, yB, DIR)

 END SUBROUTINE get_ipepo_sites

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE definitions_peps_pepo
