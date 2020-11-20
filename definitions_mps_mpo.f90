MODULE definitions_mps_mpo

 USE utility

 IMPLICIT NONE

 !! Define the types for mps & mpo blocks 
 !! Dimensions: South, North (vertical), West, East (horizontal)
 TYPE block_mpo
   INTEGER :: Sdim, Ndim, Wdim, Edim
   COMPLEX(KIND=DP), ALLOCATABLE :: m(:,:,:,:)
 END TYPE block_mpo

 TYPE transfer_mpo
   COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:,:), M(:,:,:,:), D(:,:,:)
 end TYPE transfer_mpo

 TYPE block_mps
   INTEGER :: SNdim, Wdim, Edim
   COMPLEX(KIND=DP), ALLOCATABLE :: m(:,:,:) 
   LOGICAL :: is_mult = .FALSE.
 END TYPE block_mps

 TYPE block_lambda
    INTEGER :: ChiDim
    COMPLEX(KIND=DP), ALLOCATABLE :: m(:) 
 END TYPE block_lambda 

 interface allocate_mps_site
    module procedure allocate_mps_site_u_input_tensor
    module procedure allocate_mps_site_u_input_dims
    module procedure allocate_mps_site_u_vec_of_dims
 end interface allocate_mps_site

 interface allocate_lambda_site
    module procedure allocate_lambda_site_u_input_tensor
    module procedure allocate_lambda_site_u_input_dims
 end interface allocate_lambda_site

 interface allocate_mpo_site
    module procedure allocate_mpo_site_u_input_tensor
    module procedure allocate_mpo_site_u_input_dims
    module procedure allocate_mpo_site_u_vec_of_dims
 end interface allocate_mpo_site

 private allocate_mpo_site_u_input_tensor,    allocate_mpo_site_u_input_dims,    allocate_mpo_site_u_vec_of_dims
 private allocate_mps_site_u_input_tensor,    allocate_mps_site_u_input_dims,    allocate_mps_site_u_vec_of_dims
 private allocate_lambda_site_u_input_tensor, allocate_lambda_site_u_input_dims


CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE/COPY MPS BLOCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Allocate empty mps block without any sites in it !!!
 SUBROUTINE allocate_empty_mps_block(mps_block, N_sites)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:)
  INTEGER,                      INTENT(IN)    :: N_sites

  CALL deallocate_mps_block(mps_block)

  !allocate new empty block (new sites to be appended)
  ALLOCATE(mps_block(N_sites))

 END SUBROUTINE allocate_empty_mps_block



 !!! Deallocate mps block !!!
 SUBROUTINE deallocate_mps_block(mps_block)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:)

  !dims & indices
  INTEGER :: N_sites, site

  N_sites = SIZE(mps_block)

  IF(ALLOCATED(mps_block)) THEN
    !if previously allocated, deallocate all sites first
    DO site=1,N_sites
         IF(ALLOCATED(mps_block(site)%m)) DEALLOCATE(mps_block(site)%m)
    end DO
    !then deallocate the block itself
    DEALLOCATE(mps_block)
  end IF

 END SUBROUTINE deallocate_mps_block



 !!! Copy mps block !!!
 SUBROUTINE copy_mps_block(mps_copy, mps_block, hc_flag)

   TYPE(block_mps), ALLOCATABLE, INTENT(INOUT)        :: mps_copy(:)
   TYPE(block_mps),              INTENT(IN)           :: mps_block(:)
   CHARACTER(LEN=*),             INTENT(IN), OPTIONAL :: hc_flag

   !dims & indices
   INTEGER :: N_sites, site

   N_sites = SIZE(mps_block)

   !Allocate all sites of the vector
   CALL allocate_empty_mps_block(mps_copy, N_sites)

   IF(PRESENT(hc_flag)) THEN
      !create hermitian conjugate of evec_original if requested
      Do site=1,N_sites
         CALL allocate_mps_site(mps_copy(site), CONJG(mps_block(site)%m))
      end Do
   ELSE
      !else just a simple copy of evec_original
      Do site=1,N_sites
         CALL allocate_mps_site(mps_copy(site), mps_block(site)%m)
      end Do
   end IF

 END SUBROUTINE copy_mps_block

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!! ALLOCATE/COPY LAMBDA BLOCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Allocate empty lambda block !!!
 SUBROUTINE allocate_empty_lambda_block(lambda_block, N_sites)

  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: lambda_block(:)
  INTEGER, INTENT(IN) :: N_sites

  CALL deallocate_lambda_block(lambda_block)

  !allocate new empty block (new sites to be appended)
  ALLOCATE(lambda_block(N_sites))

 END SUBROUTINE allocate_empty_lambda_block



 !!! Deallocate lambda block !!!
 SUBROUTINE deallocate_lambda_block(lambda_block)

  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: lambda_block(:)

  !dims & indices
  INTEGER :: N_sites, site

  N_sites = SIZE(lambda_block)

  IF(ALLOCATED(lambda_block)) THEN
    !if previously allocated, deallocate all sites first
    DO site=1,N_sites
         IF(ALLOCATED(lambda_block(site)%m)) DEALLOCATE(lambda_block(site)%m)
    end DO
    !then deallocate the block itself
    DEALLOCATE(lambda_block)
  end IF

 END SUBROUTINE deallocate_lambda_block



 !!! Copy lambda block !!!
 SUBROUTINE copy_lambda_block(lambda_copy, lambda_block)

   TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: lambda_copy(:)
   TYPE(block_lambda),              INTENT(IN)    :: lambda_block(:)

   !dims & indices
   INTEGER :: N_sites, site

   N_sites = SIZE(lambda_block)

   !Allocate all sites of the vector
   CALL allocate_empty_lambda_block(lambda_copy, N_sites)

   !else just a simple copy of lambda_original
   Do site=1,N_sites
      CALL allocate_lambda_site(lambda_copy(site), lambda_block(site)%m)
   end Do

 END SUBROUTINE copy_lambda_block

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE/COPY MPO BLOCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !allocate empty mpo block without any sites in it
 SUBROUTINE allocate_empty_mpo_block(mpo_block, N_sites)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: mpo_block(:)
  INTEGER, INTENT(IN) :: N_sites

  CALL deallocate_mpo_block(mpo_block)

  !allocate new empty block (new sites to be appended)
  ALLOCATE(mpo_block(N_sites))

 END SUBROUTINE allocate_empty_mpo_block



 !deallocate mpo block
 SUBROUTINE deallocate_mpo_block(mpo_block)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: mpo_block(:)

  !dims & indices
  INTEGER :: N_sites, site

  N_sites = SIZE(mpo_block)

  IF(ALLOCATED(mpo_block)) THEN
    !if previously allocated, deallocate all sites first
    DO site=1,N_sites
       IF(ALLOCATED(mpo_block(site)%m)) DEALLOCATE(mpo_block(site)%m)
    end DO
    !then deallocate the block itself
    DEALLOCATE(mpo_block)
  end IF

 END SUBROUTINE deallocate_mpo_block



 !!! Copy mpo block !!!
 SUBROUTINE copy_mpo_block(mpo_copy, mpo_block, hc_flag)

   TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT)         :: mpo_copy(:)
   TYPE(block_mpo),              INTENT(IN)            :: mpo_block(:)
   CHARACTER(LEN=*),             INTENT(IN), OPTIONAL  :: hc_flag

   !! Dims & indices
   INTEGER :: N_sites, site

   N_sites = SIZE(mpo_block)

   !! Allocate all sites of the vector
   CALL allocate_empty_mpo_block(mpo_copy, N_sites)

   IF(PRESENT(hc_flag)) THEN
      !! Create hermitian conjugate of mpo if requested
      Do site=1,N_sites
         CALL allocate_mpo_site(mpo_copy(site),  TensTRANSPOSE(CONJG(mpo_block(site) % m), '12'))
      end Do
   ELSE
      !! Else just a simple copy of mpo
      Do site=1,N_sites
         CALL allocate_mpo_site(mpo_copy(site),  mpo_block(site) % m)
      end Do
   end IF

 END SUBROUTINE copy_mpo_block

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE/COPY BINARY MPO BLOCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !allocate empty mpo block without any sites in it
 SUBROUTINE allocate_empty_bimpo_block(mpo_block, N_sites)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: mpo_block(:,:)
  INTEGER,                      INTENT(IN)    :: N_sites

  CALL deallocate_bimpo_block(mpo_block)

  !allocate new empty block (new sites to be appended)
  ALLOCATE(mpo_block(2, N_sites))

 END SUBROUTINE allocate_empty_bimpo_block



 !deallocate mpo block
 SUBROUTINE deallocate_bimpo_block(mpo_block)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: mpo_block(:,:)

  !dims & indices
  INTEGER :: N_sites, site

  N_sites = SIZE(mpo_block,2)

  IF(ALLOCATED(mpo_block)) THEN
    !if previously allocated, deallocate all sites first
    DO site=1,N_sites
       IF(ALLOCATED(mpo_block(1,site)%m)) DEALLOCATE(mpo_block(1,site)%m)
       IF(ALLOCATED(mpo_block(2,site)%m)) DEALLOCATE(mpo_block(2,site)%m)
    end DO
    !then deallocate the block itself
    DEALLOCATE(mpo_block)
  end IF

 END SUBROUTINE deallocate_bimpo_block




 !combine two mpos into bin-mpo
 SUBROUTINE combine_two_mpos_into_bimpo(bimpo_block, mpo1, mpo2)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: bimpo_block(:,:)
  TYPE(block_mpo),              INTENT(IN)    :: mpo1(:), mpo2(:)

  INTEGER :: site, N_sites

  !assert sizes are equal
  IF(SIZE(mpo1) .NE. SIZE(mpo2)) THEN
      WRITE(*,*) "mpo1 and mpo2 must have equal size"
      STOP
  ELSE
      N_sites = SIZE(mpo1)
  end IF

  !create empty bin-mpo
  CALL allocate_empty_bimpo_block(bimpo_block, N_sites)

  !copy mpo1, mpo2 to bimpo(1,:), bimpo(2,:) respectively
  Do site=1,N_sites
     CALL allocate_mpo_site(bimpo_block(1,site), mpo1(site)%m)
     CALL allocate_mpo_site(bimpo_block(2,site), mpo2(site)%m)
  end Do

 END SUBROUTINE combine_two_mpos_into_bimpo




 !copy mpo block
 SUBROUTINE copy_bimpo_block(mpo_copy, mpo_block)

   TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: mpo_copy(:,:)
   TYPE(block_mpo),              INTENT(IN)    :: mpo_block(:,:)

   !dims & indices
   INTEGER :: N_sites, site, blk

   N_sites = SIZE(mpo_block,2)

   !Allocate all sites of the vector
   CALL allocate_empty_bimpo_block(mpo_copy, N_sites)

   !else just a simple copy of evec_original
   Do blk=1,2
      Do site=1,N_sites
         CALL allocate_mpo_site(mpo_copy(blk,site), mpo_block(blk,site)%m)
      end DO
   end Do

 END SUBROUTINE copy_bimpo_block

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE/COPY (INFINITE) iMPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !Allocate empty iMPS
 SUBROUTINE allocate_empty_imps(iGamma, iLambda, N_sites_in)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT)         :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)         :: iLambda(:)
  INTEGER,                         INTENT(IN),  OPTIONAL :: N_sites_in

  INTEGER :: N_sites

  !set input or default val
  N_sites = OptArg(N_sites_in, 2)
  
  !allocate empty iGamma && iLambda blocks
  CALL allocate_empty_mps_block(iGamma, N_sites)
  CALL allocate_empty_lambda_block(iLambda, N_sites)

 END SUBROUTINE allocate_empty_imps




 !Deallocate iMPS
 SUBROUTINE deallocate_imps(iGamma, iLambda)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)

  !Deallocate iGamma && iLambda blocks
  CALL deallocate_mps_block(iGamma)
  CALL deallocate_lambda_block(iLambda)

 END SUBROUTINE deallocate_imps





 !Copy iMPS
 SUBROUTINE copy_imps(iGamma_copy, iLambda_copy, iGamma, iLambda, hc_flag)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma_copy(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda_copy(:)
  TYPE(block_mps),                 INTENT(IN)    :: iGamma(:)
  TYPE(block_lambda),              INTENT(IN)    :: iLambda(:)
  CHARACTER(LEN=*),                INTENT(IN),  OPTIONAL  :: hc_flag

  !Copy iGamma block (HC or simple copy)
  IF(PRESENT(hc_flag)) THEN
     CALL copy_mps_block(iGamma_copy, iGamma, hc_flag)
  ELSE
     CALL copy_mps_block(iGamma_copy, iGamma)
  end IF

  !Copy iLambda block
  CALL copy_lambda_block(iLambda_copy, iLambda)

 END SUBROUTINE copy_imps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE/COPY (INFINITE) iMPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !Allocate empty iMPS
 SUBROUTINE allocate_empty_impo(iGamma, iLambda, N_sites_in)

  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT)          :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)          :: iLambda(:)
  INTEGER,                         INTENT(IN),   OPTIONAL :: N_sites_in

  INTEGER :: N_sites

  !set input or default val
  N_sites = OptArg(N_sites_in, 2)
  
  !allocate empty iGamma && iLambda blocks
  CALL allocate_empty_mpo_block(iGamma, N_sites)
  CALL allocate_empty_lambda_block(iLambda, N_sites)

 END SUBROUTINE allocate_empty_impo




 !Deallocate iMPS
 SUBROUTINE deallocate_impo(iGamma, iLambda)

  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)

  !Deallocate iGamma && iLambda blocks
  CALL deallocate_mpo_block(iGamma)
  CALL deallocate_lambda_block(iLambda)

 END SUBROUTINE deallocate_impo





 !!! Copy iMPO !!!
 SUBROUTINE copy_impo(iGamma_copy, iLambda_copy, iGamma, iLambda, hc_flag)

  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT)         :: iGamma_copy(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)         :: iLambda_copy(:)
  TYPE(block_mpo),                 INTENT(IN)            :: iGamma(:)
  TYPE(block_lambda),              INTENT(IN)            :: iLambda(:)
  CHARACTER(LEN=*),                INTENT(IN), OPTIONAL  :: hc_flag

  !Copy iGamma block (HC or simple copy)
  IF(PRESENT(hc_flag)) THEN
     CALL copy_mpo_block(iGamma_copy, iGamma, hc_flag)
  ELSE
     CALL copy_mpo_block(iGamma_copy, iGamma)
  end IF

  !Copy iLambda block
  CALL copy_lambda_block(iLambda_copy, iLambda)

 END SUBROUTINE copy_impo

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE/COPY TRANSFER MPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Allocate empty transfer MPO !!!
 SUBROUTINE allocate_empty_transfer_mpo(Tmpo, N_sites)

  TYPE(transfer_mpo), ALLOCATABLE, INTENT(INOUT) :: Tmpo(:)
  INTEGER,                         INTENT(IN)    :: N_sites

  CALL deallocate_transfer_mpo(Tmpo)

  !! Allocate new empty Tmpo
  ALLOCATE(Tmpo(N_sites))

 END SUBROUTINE allocate_empty_transfer_mpo




 !!! Deallocate transfer MPO !!!
 SUBROUTINE deallocate_transfer_mpo(Tmpo)

  TYPE(transfer_mpo), ALLOCATABLE, INTENT(INOUT) :: Tmpo(:)

  INTEGER :: N_sites, site

  N_sites = SIZE(Tmpo)

  IF(ALLOCATED(Tmpo)) THEN

     !if previously allocated, deallocate all sites first
     DO site=1,N_sites
        IF(ALLOCATED(Tmpo(site) % U)) DEALLOCATE(Tmpo(site) % U)
        IF(ALLOCATED(Tmpo(site) % D)) DEALLOCATE(Tmpo(site) % D)
        IF(ALLOCATED(Tmpo(site) % M)) DEALLOCATE(Tmpo(site) % M)
     end DO

     !then deallocate the block itself
     DEALLOCATE(Tmpo)
  end IF

 END SUBROUTINE deallocate_transfer_mpo




 !!! Copy transfer MPO !!!
 SUBROUTINE copy_transfer_mpo(Tmpo_copy, Tmpo)

  TYPE(transfer_mpo), ALLOCATABLE, INTENT(INOUT) :: Tmpo_copy(:)
  TYPE(transfer_mpo),              INTENT(IN)    :: Tmpo(:)

  !! Dims && indices
  INTEGER :: N_sites, site

  !! Get size of Tmpo
  N_sites = SIZE(Tmpo)

  !! Allocate new empty Tmpo
  CALL allocate_empty_transfer_mpo(Tmpo_copy, N_sites)

  !! Copy Tmpo sites
  DO site=1,N_sites
     CALL copyTens(Tmpo_copy(site) % U,  Tmpo(site) % U)
     CALL copyTens(Tmpo_copy(site) % D,  Tmpo(site) % D)
     CALL copyTens(Tmpo_copy(site) % M,  Tmpo(site) % M)     
  end DO

 END SUBROUTINE copy_transfer_mpo

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE MPS SITE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE allocate_mps_site_u_input_tensor(mps_site, input_tensor)

  TYPE(block_mps),  INTENT(INOUT) :: mps_site
  COMPLEX(KIND=DP), INTENT(IN)    :: input_tensor(:,:,:)

  !set dimensions of mps_site
  mps_site%SNdim = SIZE(input_tensor,1)
  mps_site%Wdim = SIZE(input_tensor,2)
  mps_site%Edim = SIZE(input_tensor,3)

  !allocate mps_site
  IF(ALLOCATED(mps_site%m)) DEALLOCATE(mps_site%m)
  ALLOCATE(mps_site%m(mps_site%SNdim, mps_site%Wdim, mps_site%Edim))

  !Initialize mps_site
  mps_site%m = input_tensor

 END SUBROUTINE allocate_mps_site_u_input_tensor



 SUBROUTINE allocate_mps_site_u_input_dims(mps_site, SNdim, Wdim, Edim)

  TYPE(block_mps), INTENT(INOUT) :: mps_site
  INTEGER, INTENT(IN) :: SNdim, Wdim, Edim
 
  !set dimensions of mps_site
  mps_site%SNdim = SNdim 
  mps_site%Wdim = Wdim
  mps_site%Edim = Edim

  !allocate mps_site
  IF(ALLOCATED(mps_site%m)) DEALLOCATE(mps_site%m)
  ALLOCATE(mps_site%m(mps_site%SNdim, mps_site%Wdim, mps_site%Edim))

  !Initialize mps_site
  mps_site%m(:,:,:)=(0d0,0d0)

 END SUBROUTINE allocate_mps_site_u_input_dims



 SUBROUTINE allocate_mps_site_u_vec_of_dims(mps_site, dim_mps)

  TYPE(block_mps), INTENT(INOUT) :: mps_site
  INTEGER, INTENT(IN) :: dim_mps(3)

  !set dimensions of mps_site
  mps_site%SNdim = dim_mps(1) 
  mps_site%Wdim = dim_mps(2) 
  mps_site%Edim = dim_mps(3) 

  !allocate mps_site
  IF(ALLOCATED(mps_site%m)) DEALLOCATE(mps_site%m)
  ALLOCATE(mps_site%m(mps_site%SNdim, mps_site%Wdim, mps_site%Edim))

  !Initialize mps_site
  mps_site%m(:,:,:)=(0d0,0d0)

 END SUBROUTINE allocate_mps_site_u_vec_of_dims

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE LAMBDA SITE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !allocate lambda site using an input tensor
 SUBROUTINE allocate_lambda_site_u_input_tensor(lambda_site, input_tensor)

  TYPE(block_lambda), INTENT(INOUT) :: lambda_site
  COMPLEX(KIND=DP), INTENT(IN) :: input_tensor(:)
 
  !set dimensions of lambda_site
  lambda_site%ChiDim = SIZE(input_tensor)

  !allocate mps_site
  IF(ALLOCATED(lambda_site%m)) DEALLOCATE(lambda_site%m)
  ALLOCATE(lambda_site%m(lambda_site%ChiDim))

  !Initialize lambda_site
  lambda_site%m = input_tensor

 END SUBROUTINE allocate_lambda_site_u_input_tensor



 !allocate lambda site using input dims, initialize it to zero
 SUBROUTINE allocate_lambda_site_u_input_dims(lambda_site, EWdim)

  TYPE(block_lambda), INTENT(INOUT) :: lambda_site
  INTEGER, INTENT(IN) :: EWdim
 
  !set dimensions of lambda_site
  lambda_site%ChiDim = EWdim 

  !allocate mps_site
  IF(ALLOCATED(lambda_site%m)) DEALLOCATE(lambda_site%m)
  ALLOCATE(lambda_site%m(lambda_site%ChiDim))

  !Initialize lambda_site
  lambda_site%m(:)=(0d0,0d0)

 END SUBROUTINE allocate_lambda_site_u_input_dims

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE MPO SITE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE allocate_mpo_site_u_input_tensor(mpo_site, input_tensor)

  TYPE(block_mpo), INTENT(INOUT) :: mpo_site
  COMPLEX(KIND=DP), INTENT(IN) :: input_tensor(:,:,:,:)

  !set dimensions of mpo_site
  mpo_site%Sdim = SIZE(input_tensor,1)
  mpo_site%Ndim = SIZE(input_tensor,2)
  mpo_site%Wdim = SIZE(input_tensor,3)
  mpo_site%Edim = SIZE(input_tensor,4)

  !allocate mpo_site
  IF(ALLOCATED(mpo_site%m)) DEALLOCATE(mpo_site%m)
  ALLOCATE(mpo_site%m(mpo_site%Sdim, mpo_site%Ndim, mpo_site%Wdim, mpo_site%Edim))

  !Initialize mpo_site
  mpo_site%m = input_tensor

 END SUBROUTINE allocate_mpo_site_u_input_tensor




 SUBROUTINE allocate_mpo_site_u_input_dims(mpo_site, Sdim, Ndim, Wdim, Edim)

  TYPE(block_mpo), INTENT(INOUT) :: mpo_site
  INTEGER, INTENT(IN) :: Sdim, Ndim, Wdim, Edim
 
  !set dimensions of mpo_site
  mpo_site%Sdim = Sdim 
  mpo_site%Ndim = Ndim 
  mpo_site%Wdim = Wdim
  mpo_site%Edim = Edim

  !allocate mpo_site
  IF(ALLOCATED(mpo_site%m)) DEALLOCATE(mpo_site%m)
  ALLOCATE(mpo_site%m(mpo_site%Sdim, mpo_site%Ndim, mpo_site%Wdim, mpo_site%Edim))

  !Initialize mpo_site
  mpo_site%m(:,:,:,:)=(0d0,0d0)

 END SUBROUTINE allocate_mpo_site_u_input_dims



 SUBROUTINE allocate_mpo_site_u_vec_of_dims(mpo_site, dim_mpo)

  TYPE(block_mpo), INTENT(INOUT) :: mpo_site
  INTEGER, INTENT(IN) :: dim_mpo(4)
 
  !set dimensions of mpo_site
  mpo_site%Sdim = dim_mpo(1)  
  mpo_site%Ndim = dim_mpo(2)  
  mpo_site%Wdim = dim_mpo(3) 
  mpo_site%Edim = dim_mpo(4) 

  !allocate mpo_site
  IF(ALLOCATED(mpo_site%m)) DEALLOCATE(mpo_site%m)
  ALLOCATE(mpo_site%m(mpo_site%Sdim, mpo_site%Ndim, mpo_site%Wdim, mpo_site%Edim))

  !Initialize mpo_site
  mpo_site%m(:,:,:,:)=(0d0,0d0)

 END SUBROUTINE allocate_mpo_site_u_vec_of_dims

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DETERMINE iMPS and compact MPO SITES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Determine the order of iMPS sites !!!
 SUBROUTINE get_imps_sites(siteA, siteB, SITEORD)

  INTEGER,          INTENT(OUT) :: siteA, siteB
  CHARACTER(LEN=*), INTENT(IN)  :: SITEORD
  
  SELECT CASE(SITEORD)
  CASE('12')
       siteA = 1; siteB = 2
  CASE('21')
       siteA = 2; siteB = 1
  CASE DEFAULT
       CALL invalid_flag("get_imps_sites -- invalid SITEORD ", SITEORD)
  end SELECT

 END SUBROUTINE get_imps_sites




 !!! Determine the order of iMPO sites (Same as get_imps_sites) !!!
 SUBROUTINE get_impo_sites(siteA, siteB, SITEORD)

  INTEGER,          INTENT(OUT) :: siteA, siteB
  CHARACTER(LEN=*), INTENT(IN)  :: SITEORD
  
  CALL get_imps_sites(siteA, siteB, SITEORD)

 END SUBROUTINE get_impo_sites




 !!! Inverse iMPS site Flag !!!
 FUNCTION invSITEORD(SITEORD) 

  CHARACTER(LEN=*), INTENT(IN) :: SITEORD
  CHARACTER(LEN=2)             :: invSITEORD

  SELECT CASE(SITEORD)
  CASE('12')
       invSITEORD = '21'
  CASE('21')
       invSITEORD = '12'
  CASE DEFAULT
       CALL invalid_flag("invSITEORD -- invalid SITEORD ", SITEORD)
  end SELECT

 END FUNCTION invSITEORD




 !!! Decode Compact MPO sites !!!
 FUNCTION getSiteMpo(site, N_sites) result(siteM)

   INTEGER, INTENT(IN) :: site, N_sites
   INTEGER :: siteM

   IF(site .EQ. 1) THEN
      siteM = 1
   ELSEIF(site .EQ. N_sites) THEN
      siteM = 4
   ELSEIF(mod(site,2) .EQ. 0) THEN
      siteM = 2
   ELSEIF(mod(site,2) .NE. 0) THEN
      siteM = 3
   end IF

 END FUNCTION getSiteMpo

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE definitions_mps_mpo
