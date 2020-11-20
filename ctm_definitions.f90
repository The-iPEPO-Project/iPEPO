MODULE ctm_definitions

 USE utility
 USE definitions_mps_mpo
 USE mps_mpo_utility

 !!!!!!!!!!!!!!!!!!!!!!! CONTENTS: !!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! (1) CREATE AN INSTANCE OF CTM PARAMS
 !!
 !! (2) ALLOCATE CTM 
 !!
 !!     (2A) Allocate empty CTM blocks 
 !!     (2B) Allocate a corner matrix site 
 !!     (2C) Allocate a transfer matrix site 
 !!
 !! (3) DEALLOCATE CTM
 !!
 !! (4) COPY CTM
 !!
 !! (5) CTM OUTPUT/RELOAD ROUTINES
 !!
 !! (6) INITIALIZE CTM TENSORS -- Random or Factorized Network 
 !!
 !! (7) GET/SET CTM DIMENSIONS
 !!
 !! (8) CHECK CTM DIMENSIONS 
 !!
 !! (9) GET/SET CTM SITES and SUBSITES
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IMPLICIT NONE

 !! CTM corner matrix type
 TYPE ctm_corner_type
   COMPLEX(KIND=DP), ALLOCATABLE :: m(:,:)
 END TYPE ctm_corner_type

 !! CTM transfer matrix type
 TYPE ctm_transfer_type
   TYPE(block_mps), ALLOCATABLE :: T(:)
 END TYPE ctm_transfer_type

 !! CTM instance, defined by a unique set of params
 TYPE ctm_params
   REAL(KIND=DP)    :: eps
   INTEGER          :: chi
   REAL(KIND=DP)    :: eta
   CHARACTER(LEN=1) :: ctm_method
 END TYPE ctm_params

 !! Allocate corner matrix site
 interface allocate_ctm_corner_site
    module procedure allocate_ctm_corner_site_u_input_dims
    module procedure allocate_ctm_corner_site_u_input_tens
 end interface allocate_ctm_corner_site

 !! Allocate transfer matrix site
 interface allocate_ctm_transfer_site
    module procedure allocate_ctm_transfer_site_u_input_dims
    module procedure allocate_ctm_transfer_site_u_input_mps
 end interface allocate_ctm_transfer_site

 private allocate_ctm_corner_site_u_input_dims,   allocate_ctm_corner_site_u_input_tens
 private allocate_ctm_transfer_site_u_input_dims, allocate_ctm_transfer_site_u_input_mps

CONTAINS 

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (1) CREATE AN INSTANCE OF CTM PARAMS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Create a new instance of CTM params !!!
 SUBROUTINE create_CTM_params(CTMRG, eps, chi, eta, ctm_method)

  TYPE(CTM_params),  INTENT(INOUT) :: CTMRG       !CTM params
  REAL(KIND=DP),     INTENT(IN)    :: eps         !SVD precision
  INTEGER,           INTENT(IN)    :: chi         !SVD chi
  REAL(KIND=DP),     INTENT(IN)    :: eta         !CTM convergence precision
  CHARACTER(LEN=*),  INTENT(IN)    :: ctm_method  !CTM method

  !! Set prec params
  CTMRG % eps = eps
  CTMRG % chi = chi
  CTMRG % eta = eta

  !! Check if ctm_method input is valid
  SELECT CASE(ctm_method)
  CASE('1', '2', '3')
       CONTINUE
  CASE DEFAULT
       CALL invalid_flag("create_CTM_params -- invalid ctm_method ", ctm_method)
  end SELECT

  !! Set CTM method
  CTMRG % ctm_method = ctm_method

 END SUBROUTINE create_CTM_params

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (2) ALLOCATE CTM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Allocate new CTM !!!
 SUBROUTINE allocate_CTM(Cmat, Tmat, TN2D, chi)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2D(:)
  INTEGER,                              INTENT(IN)    :: chi(:)

  CALL allocate_ctm_corner_mat(Cmat, chi)
  CALL allocate_ctm_transfer_mat(Tmat, TN2D, chi)

 END SUBROUTINE allocate_CTM



 !!! Allocate new CTM corner matrix !!!
 SUBROUTINE allocate_CTM_corner_mat(Cmat, chi)

  TYPE(ctm_corner_type), ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  INTEGER,                            INTENT(IN)    :: chi(:)

  INTEGER, PARAMETER :: N_sites = 4
  INTEGER            :: site

  INTEGER, ALLOCATABLE :: dimC(:,:)

  !! Allocate empty corner matrix block
  CALL allocate_empty_ctm_corner_mat(Cmat, N_sites)

  !! Set CTM site dims
  CALL get_ctm_bondDims(dimC, chi, N_sites)

  !! Allocate individual corner nodes 
  DO site=1,N_sites
     CALL allocate_ctm_corner_site(Cmat(site), dimC(site, :))
  end DO

 END SUBROUTINE allocate_CTM_corner_mat




 !!! Allocate new CTM transfer matrix !!!
 SUBROUTINE allocate_CTM_transfer_mat(Tmat, TN2D, chi)

  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2D(:)
  INTEGER,                              INTENT(IN)    :: chi(:)

  INTEGER, PARAMETER :: N_sites = 4
  INTEGER, PARAMETER :: N_sub   = 2
  INTEGER            :: site, sub

  INTEGER, ALLOCATABLE :: dimT(:,:,:)
  INTEGER, ALLOCATABLE :: localDims(:,:), bondDims(:,:)

  !! Check the number of subsites is correct
  CALL check_sizes_equal(N_sub, SIZE(TN2D), "allocate_ctm_transfer_mat: N_sub must equal size of TN2D")

  !! Set CTM site dims
  CALL get_ctm_localDims(localDims, TN2D, N_sites)
  CALL get_ctm_bondDims( bondDims,  chi,  N_sites)
  
  ALLOCATE(dimT(N_sites, N_sub, 3))
  dimT(:,:,1)   = localDims
  dimT(:,:,2)   = bondDims 
  dimT(:,:,3)   = bondDims 

  !! Allocate empty transfer matrix block
  CALL allocate_empty_ctm_transfer_mat(Tmat, N_sites)

  !! Allocate all transfer matrix sites (1,2,3,4) and subsites (1,2)
  DO site=1,N_sites
     CALL allocate_ctm_transfer_site(Tmat(site), dimT(site, :, :))
  end DO

 END SUBROUTINE allocate_CTM_transfer_mat


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (2A) Allocate empty CTM blocks !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Allocate empty corner matrix !!!
 SUBROUTINE allocate_empty_CTM_corner_mat(Cmat, N_sites)

  TYPE(ctm_corner_type), ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  INTEGER,                            INTENT(IN)    :: N_sites

  !! Deallocate Cmat (if previously allocated)
  CALL deallocate_ctm_corner_mat(Cmat)

  !! Allocate new empty corner block
  ALLOCATE(Cmat(N_sites))

 END SUBROUTINE allocate_empty_CTM_corner_mat



 !!! Allocate empty transfer matrix !!! 
 SUBROUTINE allocate_empty_CTM_transfer_mat(Tmat, N_sites)

  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  INTEGER,                              INTENT(IN)    :: N_sites

  !! Deallocate Tmat (if previously allocated)
  CALL deallocate_ctm_transfer_mat(Tmat)

  !! Allocate new empty transfer block
  ALLOCATE(Tmat(N_sites))

 END SUBROUTINE allocate_empty_CTM_transfer_mat

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (2B) Allocate a corner matrix site !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Allocate C-site using input dims !!!
 SUBROUTINE allocate_CTM_corner_site_u_input_dims(Cmat, dimC)

  TYPE(ctm_corner_type), INTENT(INOUT) :: Cmat
  INTEGER,               INTENT(IN)    :: dimC(2)

  CALL allocateTens(Cmat % m, dimC)

 END SUBROUTINE allocate_CTM_corner_site_u_input_dims




 !!! Allocate C-site using input tensor !!!
 SUBROUTINE allocate_CTM_corner_site_u_input_tens(Cmat, tensIn)

  TYPE(ctm_corner_type), INTENT(INOUT) :: Cmat
  COMPLEX(KIND=DP),      INTENT(IN)    :: tensIn(:,:)                 

  CALL copyTens(Cmat % m, tensIn)

 END SUBROUTINE allocate_CTM_corner_site_u_input_tens

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (2C) Allocate a transfer matrix site !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Allocate T-site using input dims !!!
 SUBROUTINE allocate_CTM_transfer_site_u_input_dims(Tmat, dimT)

  TYPE(ctm_transfer_type), INTENT(INOUT) :: Tmat
  INTEGER,                 INTENT(IN)    :: dimT(:,:)

  INTEGER :: sub, N_sub

  !! Get MPS block size
  N_sub = SIZE(dimT, 1)

  !! Allocate transfer site = MPS block
  CALL allocate_empty_mps_block(Tmat % T, N_sub)

  !! Allocate transfer matrix subsites (1,2)
  DO sub=1,N_sub
     CALL allocate_mps_site(Tmat % T(sub),  dimT(sub, :))
  end DO

 END SUBROUTINE allocate_CTM_transfer_site_u_input_dims




 !!! Allocate T-site using input tensor !!!
 SUBROUTINE allocate_CTM_transfer_site_u_input_mps(Tmat, mps, flagUpDown)

  TYPE(ctm_transfer_type), INTENT(INOUT) :: Tmat
  TYPE(block_mps),         INTENT(IN)    :: mps(:)
  CHARACTER(LEN=*),        INTENT(IN)    :: flagUpDown

  INTEGER :: sub, subX, N_sub

  !! Get MPS block size
  N_sub = SIZE(mps)

  !! Allocate transfer site = MPS block
  CALL allocate_empty_mps_block(Tmat % T, N_sub)

  SELECT CASE(flagUpDown)
  CASE('N', 'UP')

       !! Copy "UP"/"NORTH" MPS -- need to swap site order
       DO sub=1,N_sub
          subX = 1 + mod(sub, 2)
          CALL allocate_mps_site(Tmat % T(sub),  mps(subX) % m)
       end DO

  CASE('S', 'DN')

       !! Copy "DOWN"/"SOUTH" MPS -- site order already correct, 
       !! but need to transpose each MPS tensor (to swap their leg indices)
       DO sub=1,N_sub
          CALL allocate_mps_site(Tmat % T(sub),  TensTRANSPOSE(mps(sub) % m, '23'))
       end DO

  CASE DEFAULT
      CALL invalid_flag("allocate_ctm_transfer_site_u_input_mps -- invalid flagUpDown ", flagUpDown)
  end SELECT

 END SUBROUTINE allocate_CTM_transfer_site_u_input_mps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (3) DEALLOCATE CTM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Deallocate CTM !!!
 SUBROUTINE deallocate_CTM(Cmat, Tmat)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)

  !Deallocate C-tensors
  CALL deallocate_ctm_corner_mat(Cmat)

  !Deallocate T-tensors
  CALL deallocate_ctm_transfer_mat(Tmat)

 END SUBROUTINE deallocate_CTM



 !!! Deallocate CTM corner matrix !!!
 SUBROUTINE deallocate_CTM_corner_mat(Cmat)

  TYPE(ctm_corner_type), ALLOCATABLE, INTENT(INOUT) :: Cmat(:)

  INTEGER :: N_sites, site

  !Size of Cmat
  N_sites = SIZE(Cmat)

  IF(ALLOCATED(Cmat)) THEN

    !if previously allocated, deallocate individual sites of Cmat
    DO site=1,N_sites
       IF(ALLOCATED(Cmat(site) % m)) DEALLOCATE(Cmat(site) % m) 
    end DO

    !then deallocate Cmat itself
    DEALLOCATE(Cmat)
  end IF

 END SUBROUTINE deallocate_CTM_corner_mat




 !!! Deallocate CTM transfer matrix !!!
 SUBROUTINE deallocate_CTM_transfer_mat(Tmat)

  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)

  INTEGER :: N_sites, site

  !Size of Tmat
  N_sites = SIZE(Tmat)

  IF(ALLOCATED(Tmat)) THEN

     !Deallocate all MPO blocks contained in Tmat
     DO site=1,N_sites
        CALL deallocate_mps_block(Tmat(site) % T)
     end DO

     !Then deallocate Cmat itself
     DEALLOCATE(Tmat)
  end IF

 END SUBROUTINE deallocate_CTM_transfer_mat


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (4) COPY CTM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Copy CTM tensors !!!
 SUBROUTINE copy_CTM(Cmat_copy, Tmat_copy, Cmat, Tmat)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat_copy(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat_copy(:)
  TYPE(ctm_corner_type),                INTENT(IN)    :: Cmat(:)
  TYPE(ctm_transfer_type),              INTENT(IN)    :: Tmat(:)

  !Copy corner tensors
  CALL copy_ctm_corner_mat(Cmat_copy, Cmat)

  !Copy transfer tensors
  CALL copy_ctm_transfer_mat(Tmat_copy, Tmat)

 END SUBROUTINE copy_CTM



 !!! Copy CTM corner tensors !!!
 SUBROUTINE copy_CTM_corner_mat(Cmat_copy, Cmat)

  TYPE(ctm_corner_type), ALLOCATABLE, INTENT(INOUT) :: Cmat_copy(:)
  TYPE(ctm_corner_type),              INTENT(IN)    :: Cmat(:)

  INTEGER :: site, N_sites

  !Size of Tmat
  N_sites = SIZE(Cmat)

  !Allocate new empty corner matrix
  CALL allocate_empty_ctm_corner_mat(Cmat_copy, N_sites)

  !Copy all tensors contained in Cmat
  DO site=1,N_sites
     CALL copyTens(Cmat_copy(site) % m, Cmat(site) % m)
  end DO

 END SUBROUTINE copy_CTM_corner_mat





 !!! Copy CTM transfer tensors !!!
 SUBROUTINE copy_CTM_transfer_mat(Tmat_copy, Tmat)

  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat_copy(:)
  TYPE(ctm_transfer_type),              INTENT(IN)    :: Tmat(:)

  INTEGER :: site, N_sites

  !Size of Tmat
  N_sites = SIZE(Tmat)

  !Allocate new empty transfer matrix
  CALL allocate_empty_ctm_transfer_mat(Tmat_copy, N_sites)

  !Copy all MPO blocks contained in Tmat
  DO site=1,N_sites
     CALL copy_mps_block(Tmat_copy(site) % T, Tmat(site) % T)
  end DO

 END SUBROUTINE copy_CTM_transfer_mat


 




 !!! Copy CTM edge specified by DIR
 !!! -- either create a new set of {Cmat,Tmat} to store a copy of a single edge (if CTM copy not allocated)
 !!! -- write to an existing CTM (if CTM copy already allocated)
 SUBROUTINE copy_CTM_edge(Cmat_copy, Tmat_copy, Cmat, Tmat, DIR)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat_copy(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat_copy(:)
  TYPE(ctm_corner_type),                INTENT(IN)    :: Cmat(:)
  TYPE(ctm_transfer_type),              INTENT(IN)    :: Tmat(:)
  CHARACTER(LEN=*),                     INTENT(IN)    :: DIR

  !! CTM sites
  INTEGER :: sA, sB
  INTEGER :: N_sites

  !! Size of CTM
  N_sites = SIZE(Cmat)

  !! Get the order of CTM sites on the edge specified by DIR
  CALL get_ctm_sites(sA=sA, sB=sB, DIR=DIR)

  !! Allocate new CTM edge to store a copy (2 corner matrices, 1 transfer matrix) if Cmat, Tmat not allocated
  !! -- else, just write to an existing CTM
  IF(.NOT. ALLOCATED(Cmat_copy)) CALL allocate_empty_ctm_corner_mat(Cmat_copy,   N_sites)
  IF(.NOT. ALLOCATED(Tmat_copy)) CALL allocate_empty_ctm_transfer_mat(Tmat_copy, N_sites) 

  !! Copy of corner sites on the edge
  CALL copyTens(Cmat_copy(sA) % m,   Cmat(sA) % m)
  CALL copyTens(Cmat_copy(sB) % m,   Cmat(sB) % m)

  !! Copy of transfer site (two subsites) on the edge
  CALL copy_mps_block(Tmat_copy(sA) % T,   Tmat(sA) % T)

 END SUBROUTINE copy_CTM_edge


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (5) CTM OUTPUT/RELOAD ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! This routine sets up the file for printing CTM, then prints it using output_ctm routine !!!
 SUBROUTINE CTM_print_function(Cmat, Tmat, ctm_name)

  TYPE(ctm_corner_type),   INTENT(IN) :: Cmat(:)
  TYPE(ctm_transfer_type), INTENT(IN) :: Tmat(:)
  CHARACTER(LEN=*),        INTENT(IN) :: ctm_name

  CHARACTER :: dump_ctm_filename*256

  !! Determine file name for ipeps dump
  dump_ctm_filename=TRIM(ADJUSTL(ctm_name))//".dat"
  
  CALL output_ctm(dump_ctm_filename, Cmat, Tmat)

  !! Write to an index file where we dumped the state.
  WRITE(*,*) "Output to:", dump_ctm_filename
  WRITE(*,*)

 END SUBROUTINE CTM_print_function




 !!! Output CTM to file !!!
 SUBROUTINE output_CTM(filename, Cmat, Tmat)

  CHARACTER,               INTENT(IN) :: filename*256
  TYPE(ctm_corner_type),   INTENT(IN) :: Cmat(:)
  TYPE(ctm_transfer_type), INTENT(IN) :: Tmat(:)

  !! output file handle
  INTEGER, PARAMETER :: OUTPUT=90

  !! dims & indices
  INTEGER :: N_sites, site
  INTEGER :: N_sub, sub

  !! size of ipeps
  N_sites = SIZE(Cmat)
  N_sub   = SIZE(Tmat(1) % T)

  OPEN(OUTPUT, FILE=filename, form="unformatted", status="replace")

  !! Num of sites && subsites
  WRITE(OUTPUT) N_sites
  WRITE(OUTPUT) N_sub

  !! Size of each Cmat tensor
  DO site=1,N_sites
     WRITE(OUTPUT) shape(Cmat(site) % m)
  end DO

  !! Size of each Tmat tensor
  DO site=1,N_sites
    DO sub=1,N_sub
       WRITE(OUTPUT) shape(Tmat(site) % T(sub) % m)
    end DO
  end DO

  !! Output Cmat tensors
  DO site=1,N_sites
     WRITE(OUTPUT) Cmat(site) % m
  end DO

  !! Output Tmat tensors
  DO site=1,N_sites
    DO sub=1,N_sub
       WRITE(OUTPUT) Tmat(site) % T(sub) % m
    end DO
  end DO

  CLOSE(OUTPUT)
   
 END SUBROUTINE output_CTM




 !!! Reload CTM from file !!!
 SUBROUTINE reload_CTM(filename, Cmat, Tmat)

  CHARACTER,                            INTENT(IN)    :: filename*256
  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)

  !! Dims & indices
  INTEGER, ALLOCATABLE :: dimC(:,:), dimT(:,:,:)
  INTEGER :: N_sites, site
  INTEGER :: N_sub, sub
  INTEGER :: INPUT=90

  OPEN(INPUT, FILE=filename,form="unformatted",status="old")

  !! Number of CTM sites && subsites
  READ(INPUT) N_sites
  READ(INPUT) N_sub


  !!!!!!!!!!!!!!!!!!!!!! Read CTM shapes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Allocate shape arrays
  ALLOCATE(dimC(N_sites, 2), dimT(N_sites, N_sub, 3))

  !! Read shapes of Cmat from file
  DO site=1,N_sites
     READ(INPUT) dimC(site, :)
  end DO

  !! Read shapes of Tmat from file
  DO site=1,N_sites
    DO sub=1,N_sub
       READ(INPUT) dimT(site, sub, :)
    end DO
  end DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  !!!!!!!!!!!!!!!!!!!!! Allocate CTM tensors !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Allocate empty corner && transfer blocks
  CALL allocate_empty_ctm_corner_mat(Cmat,   N_sites)
  CALL allocate_empty_ctm_transfer_mat(Tmat, N_sites)

  !! Allocate CTM corner sites
  DO site=1,N_sites
     CALL allocate_ctm_corner_site(Cmat(site), dimC(site, :))
  end DO

  !! Allocate CTM transfer sites
  DO site=1,N_sites
     CALL allocate_ctm_transfer_site(Tmat(site), dimT(site, :, :))
  end DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !! Read Cmat tensors from file
  DO site=1,N_sites
     READ(INPUT) Cmat(site) % m
  end DO

  !Read Tmat tensors from file
  DO site=1,N_sites
    DO sub=1,N_sub
       READ(INPUT) Tmat(site) % T(sub) % m
    end DO
  end DO

  CLOSE(INPUT)
  WRITE(*,*) "Succesfully reloaded CTM from file: ", filename

 END SUBROUTINE reload_CTM

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!! (6) INITIALIZE CTM TENSORS -- Random or Factorized Network !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Initialize CTM !!!
 SUBROUTINE initialize_CTM(Cmat, Tmat, TN2D, chi, init_rand_state, use_old_CTM)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT)        :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT)        :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)           :: TN2D(:)
  INTEGER,                              INTENT(IN)           :: chi(:)
  CHARACTER(LEN=*),                     INTENT(IN), OPTIONAL :: init_rand_state
  LOGICAL,                              INTENT(IN)           :: use_old_CTM

  IF(use_old_CTM) THEN

     CALL resize_transfer_mat(Tmat, TN2D) 
    
  ELSEIF(PRESENT(init_rand_state)) THEN

     CALL initialize_symmetric_random_CTM(Cmat, Tmat, TN2D, chi)
     !CALL initialize_random_CTM(Cmat, Tmat, TN2D, chi)
     !CALL initialize_CTM_using_TN2D(Cmat, Tmat, TN2D)
  ELSE
     CALL initialize_factorized_CTM(Cmat, Tmat, TN2D)
  end IF

 END SUBROUTINE initialize_CTM





 !!! Recycle old CTM to be reused in a new initialization !!!
 SUBROUTINE resize_transfer_mat(Tmat, TN2D) 

  TYPE(ctm_transfer_type), INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),         INTENT(IN)    :: TN2D(:)

  INTEGER, ALLOCATABLE :: localDims(:,:)

  INTEGER :: N_sites, N_sub
  INTEGER :: site, sub

  N_sites = SIZE(Tmat)
  N_sub = SIZE(Tmat(1) % T)

  !!! (1) Adjust "local" dims of Transfer tensors !!!

  !! (A) Get new "local" dims of current 2D network
  ALLOCATE(localDims(N_sites, N_sub))
  CALL get_ctm_localDims(localDims, TN2D, N_sites)

  !! (B) Create new (resized) T-tensors
  DO site=1,N_sites
    DO sub=1,N_sub
       CALL resize_mps_site(Tmat(site) % T(sub), localDims(site, sub))
    end DO
  end DO

  !!! (2) Corner tensors remain unmodified !!!
  
 END SUBROUTINE resize_transfer_mat




 !!! Initialize factorized CTM !!!
 SUBROUTINE initialize_factorized_CTM(Cmat, Tmat, TN2D)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2D(:) 

  !Factorized state
  INTEGER :: chi(4)
  INTEGER :: N_sites, N_sub
  INTEGER :: site, sub

  !Allocate new CTM with singleton bond dim
  chi(:) = 1
  CALL allocate_CTM(Cmat, Tmat, TN2D, chi)

  !Get allocated sizes of Cmat, Tmat blocks
  N_sites = SIZE(Tmat)
  N_sub   = SIZE(Tmat(1) % T)

  !Initialize all factorized CTM elements to one 
  !(to get initial ctm that has some overlap with converged ctm -- cf. power meth initialization)
  DO site=1,N_sites
     Cmat(site) % m(:,:) = (1.0d0, 0.0d0)
  end DO

  DO site=1,N_sites
    DO sub=1,N_sub
       Tmat(site) % T(sub) % m(:,:,:) = (1.0d0, 0.0d0)
    end DO
  end DO

 END SUBROUTINE initialize_factorized_CTM




 !!! Initialize random CTM !!!
 SUBROUTINE initialize_random_CTM(Cmat, Tmat, TN2D, chi)

  USE iteration_helper, ONLY: initialize_rand_tens

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2D(:)
  INTEGER,                              INTENT(IN)    :: chi(:)

  INTEGER :: N_sites, N_sub
  INTEGER :: site, sub

  !Allocate new empty CTM
  CALL allocate_CTM(Cmat, Tmat, TN2D, chi)

  !Get allocated sizes of Cmat, Tmat blocks
  N_sites = SIZE(Tmat)
  N_sub   = SIZE(Tmat(1) % T)

  !Initialize random C-tensors
  DO site=1,N_sites
     CALL initialize_rand_tens(Cmat(site) % m)
  end DO

  !Initialize random T-tensors
  DO site=1,N_sites
    DO sub=1,N_sub
       CALL initialize_rand_tens(Tmat(site) % T(sub) % m)
    end DO
  end DO

 END SUBROUTINE initialize_random_CTM





 !!! Initialize random CTM !!!
 SUBROUTINE initialize_symmetric_random_CTM(Cmat, Tmat, TN2D, chi)

  USE iteration_helper, ONLY: initialize_rand_tens

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2D(:)
  INTEGER,                              INTENT(IN)    :: chi(:)

  INTEGER :: N_sites, N_sub
  INTEGER :: site, sub

  COMPLEX(KIND=DP), ALLOCATABLE :: CC_rand(:,:), TT_rand(:,:,:)

  !! Allocate new empty CTM
  CALL allocate_CTM(Cmat, Tmat, TN2D, chi)

  !! Get allocated sizes of Cmat, Tmat blocks
  N_sites = SIZE(Tmat)
  N_sub   = SIZE(Tmat(1) % T)

  !! Create random C-tensor && T-tensor
  CALL allocateTens(CC_rand, shape(Cmat(1) % m))
  CALL allocateTens(TT_rand, shape(Tmat(1) % T(1) % m))

  CALL initialize_rand_tens(CC_rand)
  CALL initialize_rand_tens(TT_rand)

  !! Impose rotational invariance
  CC_rand = 0.5D0 * (CC_rand + TensTRANSPOSE(CC_rand))
  TT_rand = 0.5D0 * (TT_rand + TensTRANSPOSE(TT_rand, '23'))

  !! Copy the same random C-tensor to all Cmat sites
  DO site=1,N_sites
     CALL copyTens(Cmat(site) % m, CC_rand)
  end DO

  !! Copy the same random T-tensor to all Tmat sites
  DO site=1,N_sites
    DO sub=1,N_sub
       CALL copyTens(Tmat(site) % T(sub) % m, TT_rand)
    end DO
  end DO

 END SUBROUTINE initialize_symmetric_random_CTM




 !!! Initialize CTM sites using 2D network !!! 
 SUBROUTINE initialize_CTM_using_TN2D(Cmat, Tmat, TN2D)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2D(:)

  INTEGER, PARAMETER :: N_sites = 4
  INTEGER, PARAMETER :: N_sub   = 2

  !! Initialize Cmat, Tmat 
  CALL initialize_Cmat_using_TN2D(Cmat, TN2D, N_sites)
  CALL initialize_Tmat_using_TN2D(Tmat, TN2D, N_sites, N_sub)
 
 END SUBROUTINE initialize_CTM_using_TN2D





 !!! Initialize Cmat site using 2D network !!!
 SUBROUTINE initialize_Cmat_using_TN2D(Cmat, TN2D, N_sites)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2D(:)
  INTEGER,                              INTENT(IN)    :: N_sites

  !! Allocate empty corner matrix block
  CALL allocate_empty_ctm_corner_mat(Cmat, N_sites)

  !! Initialize Cmat sites using 2D network 
  !! (if needed, transpose the resulting matrix after summing out the legs)
  CALL allocate_ctm_corner_site(Cmat(1),   TensTRANSPOSE(SUM(SUM(TN2D(1) % m, DIM=1), DIM=2)))  !! Sum over legs = {1, 3->2}
  CALL allocate_ctm_corner_site(Cmat(3),   TensTRANSPOSE(SUM(SUM(TN2D(1) % m, DIM=2), DIM=3)))  !! Sum over legs = {2, 4->3}
  CALL allocate_ctm_corner_site(Cmat(2),                 SUM(SUM(TN2D(2) % m, DIM=1), DIM=3))   !! Sum over legs = {1, 4->3}
  CALL allocate_ctm_corner_site(Cmat(4),                 SUM(SUM(TN2D(2) % m, DIM=2), DIM=2))   !! Sum over legs = {2, 3->2}

 END SUBROUTINE initialize_Cmat_using_TN2D






 !!! Initialize Tmat site using 2D network !!!
 SUBROUTINE initialize_Tmat_using_TN2D(Tmat, TN2D, N_sites, N_sub)

  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_mpo),                      INTENT(IN)    :: TN2D(:)
  INTEGER,                              INTENT(IN)    :: N_sites, N_sub

  !! Indices
  INTEGER :: site, sub 

  !! Check the number of subsites is correct
  CALL check_sizes_equal(N_sub, SIZE(TN2D), "allocate_ctm_transfer_mat -- N_sub must equal size of TN2D")

  !! Allocate empty transfer matrix with MPS subsites
  CALL allocate_empty_ctm_transfer_mat(Tmat, N_sites)
  DO site=1,N_sites
     CALL allocate_empty_mps_block(Tmat(site) % T, N_sub)
  end DO

  !! Initialize Tmat sites using 2D network 
  !! (rotate tens so that leg-1 is pointing inwards, while leg-2 is being summed out)
  DO sub=1,N_sub
     CALL allocate_mps_site(Tmat(2) % T(sub),  SUM(TensROTATE(TN2D(sub) % m,  'CW +PI/2' ),  DIM=2))
     CALL allocate_mps_site(Tmat(4) % T(sub),  SUM(TensROTATE(TN2D(sub) % m,  'CCW -PI/2'),  DIM=2))
     CALL allocate_mps_site(Tmat(1) % T(sub),  SUM(TensROTATE(TN2D(sub) % m,  'PI'       ),  DIM=2))
     CALL allocate_mps_site(Tmat(3) % T(sub),  SUM(           TN2D(sub) % m,                 DIM=2))
  end DO

 END SUBROUTINE initialize_Tmat_using_TN2D



 !!! Add noise to CTM tensors !!!
 SUBROUTINE add_noise_CTM(Cmat, Tmat)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)

  INTEGER :: N_sites, N_sub
  INTEGER :: site, sub

  N_sites = SIZE(Cmat)
  N_sub   = SIZE(Tmat(1) % T)

  !! Add noise to C's
  DO site=1,N_sites
     Cmat(site) % m = add_noise(Cmat(site) % m)
  end DO

  !! Add noise to T's
  DO site=1,N_sites
    DO sub=1,N_sub
       Tmat(site) % T(sub) % m = add_noise(Tmat(site) % T(sub) % m)
    end DO
  end DO

 END SUBROUTINE add_noise_CTM

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
















 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (7) GET/SET CTM DIMENSIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Get CTM "local" dims connecting it to 2D network !!!
 SUBROUTINE get_CTM_localDims(localDims, TN2D, N_sites)

  INTEGER, ALLOCATABLE, INTENT(OUT) :: localDims(:,:) 
  TYPE(block_mpo),      INTENT(IN)  :: TN2D(:) 
  INTEGER,              INTENT(IN)  :: N_sites

  INTEGER :: sub, N_sub

  !! Get num of subsites (each corresponding to 2D network sites they're connected to)
  N_sub = SIZE(TN2D)

  !! Allocate localDims
  ALLOCATE(localDims(N_sites, N_sub))

  !! Get "local" dims of CTM tensors 
  DO sub=1,N_sub
     localDims(1, sub) = TN2D(sub) % Sdim
     localDims(2, sub) = TN2D(sub) % Edim
     localDims(3, sub) = TN2D(sub) % Ndim
     localDims(4, sub) = TN2D(sub) % Wdim
  end DO

 END SUBROUTINE get_CTM_localDims




 !!! Get CTM bond dims !!!
 SUBROUTINE get_CTM_bondDims(bondDims, chi, N_sites)

  INTEGER, ALLOCATABLE, INTENT(INOUT) :: bondDims(:,:)
  INTEGER,              INTENT(IN)    :: chi(:)
  INTEGER,              INTENT(IN)    :: N_sites

  INTEGER :: site

  !! Check chi(:) input has the correct size
  CALL check_sizes_equal(SIZE(chi), N_sites, "get_ctm_bondDims: size of chi must equal N_sites")

  !! Allocate bondDims
  ALLOCATE(bondDims(N_sites, 2))

  !! Get bond dims of CTM tensors
  DO site=1,N_sites
     bondDims(site, :) = chi(site)
  end DO

 END SUBROUTINE get_CTM_bondDims

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (8) CHECK CTM DIMENSIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Check all dims of Cmat and Tmat in CTM !!!
 SUBROUTINE check_CTM_dims_ALL(Cmat, Tmat, TN2D)

  TYPE(ctm_corner_type),   INTENT(IN) :: Cmat(:)
  TYPE(ctm_transfer_type), INTENT(IN) :: Tmat(:)
  TYPE(block_mpo),         INTENT(IN) :: TN2D(:)

  CALL check_CTM_localDims(Tmat, TN2D)
  CALL check_CTM_bondDims(Cmat, Tmat, 'S')
  CALL check_CTM_bondDims(Cmat, Tmat, 'N')
  CALL check_CTM_bondDims(Cmat, Tmat, 'W')
  CALL check_CTM_bondDims(Cmat, Tmat, 'E')

 END SUBROUTINE check_CTM_dims_ALL




 !!! Check all dims of Cmat and Tmat on a given edge of CTM !!!
 SUBROUTINE check_CTM_dims(Cmat, Tmat, TN2D, DIR)

  TYPE(ctm_corner_type),   INTENT(IN) :: Cmat(:)
  TYPE(ctm_transfer_type), INTENT(IN) :: Tmat(:)
  TYPE(block_mpo),         INTENT(IN) :: TN2D(:)
  CHARACTER(LEN=*),        INTENT(IN) :: DIR

  CALL check_CTM_localDims(Tmat, TN2D)
  CALL check_CTM_bondDims(Cmat, Tmat, DIR)

 END SUBROUTINE check_CTM_dims


 

 !!! Check if CTM "local" dims match the dims of 2D TN tensors they're connected to !!!
 SUBROUTINE check_CTM_localDims(Tmat, TN2D)

  TYPE(ctm_transfer_type), INTENT(IN) :: Tmat(:)
  TYPE(block_mpo),         INTENT(IN) :: TN2D(:)

  INTEGER, ALLOCATABLE :: localDims(:,:)
  INTEGER              :: N_sites, N_sub
  INTEGER              :: site, sub

  !! Vars for error message
  CHARACTER(LEN=16)  :: char_site, char_sub
  CHARACTER(LEN=128) :: msg

  !! Get size num of sites && subsites in Tmat
  N_sites = SIZE(Tmat); N_sub = SIZE(Tmat(1) % T)

  CALL get_CTM_localDims(localDims, TN2D, N_sites)

  DO site=1,N_sites
    DO sub=1,N_sub  
        
       WRITE(char_site, '(I3)') site
       WRITE(char_sub,  '(I3)') sub

       CALL check_sizes_equal(SIZE(Tmat(site) % T(sub) % m, 1), &
            & localDims(site, sub), &
            & "check_ctm_localDims: mismatch at site "//TRIM(ADJUSTL(char_site))//", sub "//TRIM(ADJUSTL(char_sub)))
    end DO
  end DO
  
 END SUBROUTINE check_CTM_localDims






 !!! Check if CTM bond dims match for consecutive sites on an edge specified by dirFlag !!!
 SUBROUTINE check_CTM_bondDims(Cmat, Tmat, DIR)

  TYPE(ctm_corner_type),   INTENT(IN) :: Cmat(:)
  TYPE(ctm_transfer_type), INTENT(IN) :: Tmat(:)
  CHARACTER(LEN=*),        INTENT(IN) :: DIR
 
  !! Dims && indices
  INTEGER :: dimCa(2), dimCb(2), dimT(2,3)
  INTEGER :: sA, sB, subA, subB

  !! Get sites of CTM tensors for a given edge direction
  CALL get_ctm_sites(sA=sA, sB=sB, subA=subA, subB=subB, DIR=DIR) 

  !! Find dimensions of CTM sites on the edge specified by dirFlag
  dimCA        = shape(Cmat(sA) % m)
  dimCB        = shape(Cmat(sB) % m)
  dimT(subA,:) = shape(Tmat(sA) % T(subA) % m)
  dimT(subB,:) = shape(Tmat(sA) % T(subB) % m)
  
  !! Check if all consecutive sites have matching dims
  CALL check_sizes_equal(dimCA(1),      dimT(subA, 3), "check_ctm_bondDims: mismatch dimCA vs dimT,  dir "//TRIM(ADJUSTL(DIR)))
  CALL check_sizes_equal(dimT(subA, 2), dimT(subB, 3), "check_ctm_bondDims: mismatch dimT  vs dimT,  dir "//TRIM(ADJUSTL(DIR)))
  CALL check_sizes_equal(dimT(subB, 2), dimCB(2),      "check_ctm_bondDims: mismatch dimT  vs dimCB, dir "//TRIM(ADJUSTL(DIR)))

 END SUBROUTINE check_CTM_bondDims


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!















 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (9) GET/SET CTM SITES and SUBSITES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Determine the order of CTM sites && subsites !!!
 SUBROUTINE get_CTM_sites(sA, sB, sC, sD, subA, subB, DIR)

  INTEGER,          INTENT(OUT), OPTIONAL :: sA, sB, sC, sD
  INTEGER,          INTENT(OUT), OPTIONAL :: subA, subB
  CHARACTER(LEN=*), INTENT(IN)            :: DIR

  INTEGER :: sAA, sBB, sCC, sDD
  INTEGER :: subAA, subBB


  !! (1) CTM sites: get sAA from dirFlag, the rest follow by cyclic permutation !!
  SELECT CASE(DIR)
  CASE('S')
       sAA = 1; sBB = 2; sCC = 3; sDD = 4
  CASE('E')
       sAA = 2; sBB = 3; sCC = 4; sDD = 1
  CASE('N')
       sAA = 3; sBB = 4; sCC = 1; sDD = 2
  CASE('W')
       sAA = 4; sBB = 1; sCC = 2; sDD = 3
  CASE DEFAULT
       CALL invalid_flag("get_ctm_sites -- incorrect DIR ", DIR)
  end SELECT


  !! (2) CTM subsites: [subAA = 1, subBB = 2] if sAA = 1,3 !!
  !!                   [subAA = 2, subBB = 1] if sAA = 2,4 !!
  subAA = 2 - mod(sAA, 2)
  subBB = 1 + mod(sAA, 2)


  !! (3) Export to output !!
  IF(PRESENT(sA)) sA = sAA
  IF(PRESENT(sB)) sB = sBB
  IF(PRESENT(sC)) sC = sCC
  IF(PRESENT(sD)) sD = sDD

  IF(PRESENT(subA)) subA = subAA
  IF(PRESENT(subB)) subB = subBB

 END SUBROUTINE get_CTM_sites




 !!! Find opposite direction to the one specified by DIR !!!
 FUNCTION oppositeDIR(DIR)

  CHARACTER(LEN=*), INTENT(IN)  :: DIR
  CHARACTER(LEN=1)              :: oppositeDIR

  SELECT CASE(DIR)
  CASE('S')
       oppositeDIR = 'N'
  CASE('N')
       oppositeDIR = 'S'
  CASE('W')
       oppositeDIR = 'E'
  CASE('E')
       oppositeDIR = 'W'
  CASE DEFAULT
      CALL invalid_flag("oppositeDIR -- invalid DIR ", DIR)
  end SELECT

 END FUNCTION oppositeDIR



 !!! Find orthogonal directions to those specified by DIR !!!
 FUNCTION orthogonalDIR(DIR)

  CHARACTER(LEN=*), INTENT(IN)  :: DIR(2)
  CHARACTER(LEN=1)              :: orthogonalDIR(2)

  IF    ((DIR(1) .EQ. 'S') .AND. (DIR(2) .EQ. 'N')) THEN
                                                    orthogonalDIR = (/ 'W', 'E' /)
  ELSEIF((DIR(1) .EQ. 'N') .AND. (DIR(2) .EQ. 'S')) THEN
                                                    orthogonalDIR = (/ 'W', 'E' /)
  ELSEIF((DIR(1) .EQ. 'W') .AND. (DIR(2) .EQ. 'E')) THEN
                                                    orthogonalDIR = (/ 'S', 'N' /)
  ELSEIF((DIR(1) .EQ. 'E') .AND. (DIR(2) .EQ. 'W')) THEN
                                                    orthogonalDIR = (/ 'S', 'N' /)
  ELSE
      CALL invalid_flag("orthogonalDIR -- invalid DIR ", DIR(1), DIR(2))
  end IF

 END FUNCTION orthogonalDIR

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE ctm_definitions
