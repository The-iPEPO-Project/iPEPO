MODULE mps_mpo_algebra_finite

 USE utility
 USE definitions_mps_mpo
 USE mps_mpo_utility

 IMPLICIT NONE

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MULT && SVD MPS-MPO SITES/BONDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Multiply sites of mps & mpo !!!
 SUBROUTINE multiply_mps_and_mpo_site(mpsA, mpsB, mpoB, chi0, eps0)

  TYPE(block_mps), INTENT(INOUT) :: mpsA, mpsB
  TYPE(block_mpo), INTENT(IN)    :: mpoB 
  INTEGER,         INTENT(IN)    :: chi0
  REAL(KIND=DP),   INTENT(IN)    :: eps0

  !! Dims
  INTEGER :: chi

  !! SVD matrices & temp tensor to store contraction results
  COMPLEX(KIND=DP), ALLOCATABLE :: theta(:,:), U(:,:), Udag(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:)

  !! (1) Construct_theta
  theta = RESHAPE_2D(mpsA % m, '12,3')

  !! (2) Compute SVD
  chi = chi0
  CALL compute_lapack_svd(theta, U, chi, eps0)

  !! (3) Copy back svd results to MPS --> copy U to MPS & create Udag
  Udag = TensTRANSPOSE(U, 'HC')
  CALL allocate_mps_site(mpsA, RESHAPE_3D(U, '12,3', (/mpsA % SNdim, mpsA % Wdim/)))
  
  !! (4) Contract Udag -- theta -- {mpoB%m, mpsB%m}  ---> write to MPS site
  tensO = TENSMUL(TENSMUL(TENSMUL(mpoB % m, mpsB % m, MULT='21', FUSE='(32,43)'), theta, '22'), Udag, '22')
  CALL allocate_mps_site(mpsB, tensO)
  
 END SUBROUTINE multiply_mps_and_mpo_site




 !!! Compute SVD of finite MPS bond !!!
 SUBROUTINE calc_svd_of_mps_bond(mpsA, mpsB, chi0, eps0)

  TYPE(block_mps), INTENT(INOUT) :: mpsA, mpsB
  INTEGER,         INTENT(IN)    :: chi0
  REAL(KIND=DP),   INTENT(IN)    :: eps0

  !! Dims
  INTEGER :: chi

  !! SVD matrices & temp tensor to store contraction results
  COMPLEX(KIND=DP), ALLOCATABLE :: theta(:,:), U(:,:), Udag(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:)

  !! (1) Construct_theta
  theta = RESHAPE_2D(mpsA % m, '12,3')

  !! (2) Compute SVD
  chi = chi0
  CALL compute_lapack_svd(theta, U, chi, eps0) 

  !! (3) Copy back svd results to MPS --> copy U to MPS & create Udag
  Udag = TensTRANSPOSE(U, 'HC')
  CALL allocate_mps_site(mpsA, RESHAPE_3D(U, '12,3', (/ mpsA % SNdim, mpsA % Wdim /)))

  !! (4) Contract Udag -- theta -- mps(site+1)
  tensO = TENSMUL(TENSMUL(mpsB % m, theta, '22'), Udag, '22')
  CALL allocate_mps_site(mpsB, tensO)

 END SUBROUTINE calc_svd_of_mps_bond

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SWEEP && MULT THROUGH MPS-MPO BLOCK NETWORK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Perform LEFT SWEEP through mps-mpo network !!!
 SUBROUTINE left_sweep_mps_mpo_network(mps_block, mpo_block, is_multiplied, orth_centre, chi0, eps0)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:)
  TYPE(block_mpo),              INTENT(IN)    :: mpo_block(:)
  LOGICAL,                      INTENT(INOUT) :: is_multiplied(:)
  INTEGER,                      INTENT(IN)    :: orth_centre
  INTEGER,                      INTENT(IN)    :: chi0
  REAL(KIND=DP),                INTENT(IN)    :: eps0

  !indices & dims
  INTEGER :: site, N_sites
  INTEGER :: SiteA, SiteB, SiteMpoB

  N_sites = SIZE(mps_block)

  !If orth_centre = start of block --> exit subroutine
  IF(orth_centre .EQ. 1) RETURN

  !Initialize block mult at site=1
  IF(.NOT. is_multiplied(1)) CALL allocate_mps_site(mps_block(1), TENSMUL(mpo_block(1) % m, mps_block(1) % m, MULT='21', FUSE='(32,43)'))

  DO site=2,orth_centre

      IF(.NOT. is_multiplied(site)) THEN 

          !Multiply mps & mpo and perform SVD on each bond
          SiteA = (site-1); SiteB = site; SiteMpoB = site 
          CALL multiply_mps_and_mpo_site(mps_block(SiteA), mps_block(SiteB), mpo_block(SiteMpoB), chi0, eps0)
          is_multiplied(site) = .TRUE. 
      ELSE

          !If site has already been multiplied, we only need to do SVD & truncation
          CALL calc_svd_of_mps_bond(mps_block(site-1), mps_block(site), chi0, eps0)
      end IF 

  end DO

 END SUBROUTINE left_sweep_mps_mpo_network





 !!! Perform RIGHT SWEEP through mps-mpo network !!!
 SUBROUTINE right_sweep_mps_mpo_network(mps_block, mpo_block, is_multiplied, orth_centre, chi0, eps0)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:)
  TYPE(block_mpo),              INTENT(IN)    :: mpo_block(:)
  LOGICAL,         ALLOCATABLE, INTENT(INOUT) :: is_multiplied(:)
  INTEGER,                      INTENT(IN)    :: orth_centre
  INTEGER,                      INTENT(IN)    :: chi0
  REAL(KIND=DP),                INTENT(IN)    :: eps0

  TYPE(block_mpo), ALLOCATABLE :: mpo_inverse(:)

  !if orth_centre = start of block --> exit subroutine
  IF(orth_centre .EQ. 1) RETURN

  !Before R sweep --> invert mps_block & mpo_block
  CALL invert_mps_block(mps_block)
  CALL copy_mpo_block(mpo_inverse, mpo_block); CALL invert_mpo_block(mpo_inverse)
  CALL invert_array(is_multiplied)

  !Use L sweep routine, but with inverted (mps & mpo) network
  CALL left_sweep_mps_mpo_network(mps_block, mpo_inverse, is_multiplied, orth_centre, chi0, eps0)

  !After R sweep --> must revert mps_block & deallocate mpo_inverse
  CALL invert_mps_block(mps_block); CALL deallocate_mpo_block(mpo_inverse)
  CALL invert_array(is_multiplied)

 END SUBROUTINE right_sweep_mps_mpo_network

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SWEEP && SVD THROUGH MPS BLOCK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Canonicalize (FINITE) MPS BLOCK !!!
 SUBROUTINE canonicalize_mps_block(mps_block, orth_centre, chi, eps)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:)
  INTEGER,                      INTENT(IN)    :: orth_centre
  INTEGER,                      INTENT(IN)    :: chi(2)
  REAL(KIND=DP),                INTENT(IN)    :: eps(2)

  INTEGER :: N_sites

  !! Determine the length of mps
  N_sites = SIZE(mps_block)

  !! Left & Right sweeps
  CALL left_sweep_mps_block(mps_block,  orth_centre,               chi(1), eps(1))
  CALL right_sweep_mps_block(mps_block, N_sites - orth_centre + 1, chi(2), eps(2))

 END SUBROUTINE canonicalize_mps_block



 !!! Right SVD sweep through mps_block !!!
 SUBROUTINE right_sweep_mps_block(mps_block, orth_centre, chi, eps)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:)
  INTEGER,                      INTENT(IN)    :: orth_centre
  INTEGER,                      INTENT(IN)    :: chi
  REAL(KIND=DP),                INTENT(IN)    :: eps

  !if orth_centre = start of block --> exit subroutine
  IF(orth_centre .EQ. 1) RETURN

  !Before R sweep --> invert mps_block
  CALL invert_mps_block(mps_block)

  !Use L sweep routine, but with inverted mps_block
  CALL left_sweep_mps_block(mps_block, orth_centre, chi, eps)

  !After R sweep --> revert mps_block back
  CALL invert_mps_block(mps_block)

 END SUBROUTINE right_sweep_mps_block




 !!! Left SVD sweep through mps_block !!!
 SUBROUTINE left_sweep_mps_block(mps_block, orth_centre, chi, eps)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:)
  INTEGER,                      INTENT(IN)    :: orth_centre
  INTEGER,                      INTENT(IN)    :: chi
  REAL(KIND=DP),                INTENT(IN)    :: eps

  !indices & dims
  INTEGER :: site

  !if orth_centre = start of block --> exit subroutine
  IF(orth_centre .EQ. 1) RETURN

  !sweep through all sites of mps_block performing SVD & truncating
  DO site=2,orth_centre
     CALL calc_svd_of_mps_bond(mps_block(site-1), mps_block(site), chi, eps) 
  end DO

 END SUBROUTINE left_sweep_mps_block

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MULT && SVD MPS-MPO BLOCKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Mult mps_block * binary_mpo_block !!!
 SUBROUTINE multiply_mps_bimpo(mps_block, mpo_block, chi, eps)  

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT)         :: mps_block(:)   !MPS
  TYPE(block_mpo),              INTENT(IN)            :: mpo_block(:,:) !binary MPO
  INTEGER,                      INTENT(IN)            :: chi(2)         !svd params
  REAL(KIND=DP),                INTENT(IN)            :: eps(2)         !svd params        

  CALL multiply_mps_mpo(mps_block, mpo_block(2,:), chi, eps)
  CALL multiply_mps_mpo(mps_block, mpo_block(1,:), chi, eps)

 END SUBROUTINE multiply_mps_bimpo



 !!! Multiply MPS_BLOCK * MPO_BLOCK (full multiplication with SVD) !!!
 SUBROUTINE multiply_mps_mpo(mps_block, mpo_block, chi, eps, excludeSites)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT)        :: mps_block(:)    !MPS
  TYPE(block_mpo),              INTENT(IN)           :: mpo_block(:)    !MPO
  INTEGER,                      INTENT(IN)           :: chi(2)          !svd params
  REAL(KIND=DP),                INTENT(IN)           :: eps(2)          !svd params         
  INTEGER,                      INTENT(IN), OPTIONAL :: excludeSites(:) !optional exclude MPS sites from multiplication

  INTEGER              :: orth_centre, N_sites, i
  LOGICAL, ALLOCATABLE :: is_multiplied(:)

  !! Size of mps-mpo network && pos of orth_centre (mixed canonical: CEILING(0.5*N_sites))
  N_sites     = SIZE(mps_block) 
  orth_centre = N_sites

  !! Alloc & init [is_multiplied] array to track which mps-mpo sites have been multiplied
  ALLOCATE(is_multiplied(N_sites)); is_multiplied(1:N_sites) = .FALSE.  

  !! If any sites have been multiplied during an earlier process
  IF(PRESENT(excludeSites)) THEN
     DO i=1,SIZE(excludeSites)
        is_multiplied(excludeSites(i)) = .TRUE.
     end DO
  end IF

  !! (1) Sweep Left--to--Right through MPS-MPO network (the default choice)
  CALL left_sweep_mps_mpo_network(mps_block,  mpo_block, is_multiplied, orth_centre, chi(1), eps(1))
  !CALL right_sweep_mps_mpo_network(mps_block, mpo_block, is_multiplied, N_sites - orth_centre + 1, chi(1), eps(1))
  
  !! (2) Backward sweep Left--Right: 
  !!     canonicalize post-mult MPS --> [Oc on backward sweep] = opposite to [initial val of Oc]
  CALL canonicalize_mps_block(mps_block, N_sites - orth_centre + 1,  (/chi(2), chi(2)/),  (/eps(2), eps(2)/))

 END SUBROUTINE multiply_mps_mpo




 !!! Simple multiplication: MPO * MPS blocks (no SVD) !!!
 SUBROUTINE simple_mult_mps_mpo(mps_block, mpo_block)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:) 
  TYPE(block_mpo),              INTENT(IN)    :: mpo_block(:) 

  !! Indices
  INTEGER :: site, N_sites

  N_sites = SIZE(mpo_block)

  !! Multiply mpo_block * mps_block
  DO site=1,N_sites
     CALL allocate_mps_site(mps_block(site), TENSMUL(mpo_block(site) % m, mps_block(site) % m, MULT='21', FUSE='(32,43)'))
  end DO

 END SUBROUTINE simple_mult_mps_mpo




 !!! Combine two MPS vectors into a single MPO matrix !!!
 SUBROUTINE kronecker_mult_mps_mps(mpoCombo, mpsR, mpsL)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: mpoCombo(:)
  TYPE(block_mps),              INTENT(IN)    :: mpsR(:), mpsL(:)

  INTEGER :: site, N_sites

  N_sites = SIZE(mpsR)

  CALL allocate_empty_mpo_block(mpoCombo, N_sites)

  DO site=1,N_sites
     CALL allocate_mpo_site(mpoCombo(site), TensKRON(mpsR(site) % m, mpsL(site) % m))
  end DO

 END SUBROUTINE kronecker_mult_mps_mps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MULT && SVD MPO-MPO BLOCKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Simple multiplication: MPO * MPO blocks (no SVD) !!!
 SUBROUTINE simple_mult_MPO_MPO(mpoSELF, mpoOTHER, MULT)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: mpoSELF(:) 
  TYPE(block_mpo),              INTENT(IN)    :: mpoOTHER(:) 
  CHARACTER(LEN=*),             INTENT(IN)    :: MULT

  !! Indices
  INTEGER :: site, N_sites

  !! Get MPO size
  N_sites = SIZE(mpoSELF)

  !! Multiply mpoSELF * mpoOTHER
  SELECT CASE(MULT)
  CASE('R')
      DO site=1,N_sites
         CALL allocate_mpo_site(mpoSELF(site),  TENSMUL(mpoOTHER(site) % m,  mpoSELF(site) % m,   MULT='21',  FUSE='(33,44)'))
      end DO

  CASE('L')
      DO site=1,N_sites
         CALL allocate_mpo_site(mpoSELF(site),  TENSMUL(mpoSELF(site) % m,   mpoOTHER(site) % m,  MULT='21',  FUSE='(33,44)'))
      end DO

  CASE DEFAULT
      CALL invalid_flag("simple_mult_mpo_mpo -- invalid MULT ", MULT)
  end SELECT

 END SUBROUTINE simple_mult_MPO_MPO

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CONTRACT MPS--MPO--MPS, MPS--biMPO--MPS, and MPS--MPS networks !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! CONTRACT MPS-MPS NETWORK, FIND EXPVAL !!!
 SUBROUTINE contract_mps_mps(expval, mps1, mps2)

  COMPLEX (KIND=DP), INTENT(OUT) :: expval
  TYPE(block_mps),   INTENT(IN)  :: mps1(:), mps2(:)

  !! Overlap variable
  COMPLEX(KIND=DP), ALLOCATABLE :: mps_overlap(:,:)

  !! Indices
  INTEGER :: N_sites, site

  N_sites = SIZE(mps1)

  !! Initialize mps_overlap at site=1 of mps-mps network
  mps_overlap = TENSMUL(mps1(1) % m, mps2(1) % m, MULT='11', FUSE='(22,33)')

  !! Sweep through mps-mps network, zip the two MPS site-by-site & contract with mps_overlap
  !! the whole mps-mps network is contracted into a single tensor
  DO site=2,N_sites
     mps_overlap = TENSMUL(mps_overlap, TENSMUL(mps1(site) % m, mps2(site) % m, MULT='11', FUSE='(22,33)'))
  end DO

  !! Check if mps_overlap indeed has singleton dim
  IF((SIZE(mps_overlap,1) .NE. 1) .OR. (SIZE(mps_overlap,2) .NE. 1)) THEN
     WRITE(*,*) "Error in contract_two_mps: mps_overlap has dim > 1: ", shape(mps_overlap)
     STOP
  end IF

  !! Return output
  expval = mps_overlap(1,1)
  DEALLOCATE(mps_overlap)

 END SUBROUTINE contract_mps_mps




 !!! Contract mps-mpo-mps network, find expval !!!
 SUBROUTINE contract_mps_mpo_mps(expval, mps1, mpo, bimpo, mps2, chi, eps)

  COMPLEX(KIND=DP),  INTENT(OUT)          :: expval
  TYPE(block_mps),   INTENT(IN)           :: mps1(:)
  TYPE(block_mps),   INTENT(IN)           :: mps2(:)
  TYPE(block_mpo),   INTENT(IN), OPTIONAL :: mpo(:)  
  TYPE(block_mpo),   INTENT(IN), OPTIONAL :: bimpo(:,:)  
  INTEGER,           INTENT(IN)           :: chi(2) 
  REAL(KIND=DP),     INTENT(IN)           :: eps(2)

  !! Local copy of mps1
  TYPE(block_mps), ALLOCATABLE :: mpsC1(:)

  !! Create a local copy of mps1
  CALL copy_mps_block(mpsC1, mps1)

  !! Mult MPS-1 with MPO or biMPO
  IF(PRESENT(mpo)) THEN

       CALL multiply_mps_mpo(mpsC1,   mpo,   chi, eps)

  ELSEIF(PRESENT(bimpo)) THEN

       CALL multiply_mps_bimpo(mpsC1, bimpo, chi, eps)
  end IF

  !! Calc normalized expval
  CALL contract_mps_mps(expval, mpsC1, mps2)

 END SUBROUTINE contract_mps_mpo_mps




 !!! Normalize MPS !!!
 SUBROUTINE normalize_mps(mps_block, C_mps)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT)         :: mps_block(:)
  COMPLEX(KIND=DP),             INTENT(OUT), OPTIONAL :: C_mps 

  !HC of mps_block
  TYPE(block_mps), ALLOCATABLE :: mps_hc(:)

  !norm const
  COMPLEX(KIND=DP) :: C_norm

  !dims & indices
  INTEGER :: site
  
  !Obtain norm constant C_norm
  CALL copy_mps_block(mps_hc, mps_block, 'HC')
  CALL contract_mps_mps(C_norm, mps_block, mps_hc)

  !Normalize mps_block by C_norm
  CALL mult_mps_by_value(mps_block, C_norm**(-0.5))

  !Output C_norm if present
  IF(PRESENT(C_mps)) C_mps = C_norm

 END SUBROUTINE normalize_mps

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE mps_mpo_algebra_finite
