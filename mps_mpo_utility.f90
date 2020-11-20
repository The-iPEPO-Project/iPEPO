MODULE mps_mpo_utility

 USE utility
 USE definitions_mps_mpo

 IMPLICIT NONE

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ABSORB iMPS LAMBDAS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Absorb iMPS lambdas into gammas !!!
 SUBROUTINE absorb_imps_lambdas(imps, iGamma, iLambda, SITEORD, SYM)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT)        :: imps(:)
  TYPE(block_mps),              INTENT(IN), OPTIONAL :: iGamma(:)
  TYPE(block_lambda),           INTENT(IN)           :: iLambda(:)
  CHARACTER(LEN=*),             INTENT(IN), OPTIONAL :: SITEORD
  CHARACTER(LEN=*),             INTENT(IN), OPTIONAL :: SYM

  !! iMPS sites
  INTEGER :: sA, sB

  !! Copy iGamma to iMPS block (if iGamma not present, use iMPS as input gamma instead)
  IF(PRESENT(iGamma)) CALL copy_mps_block(imps, iGamma) 

  !! Determine siteA, siteB
  CALL get_imps_sites(sA, sB, OptArg(SITEORD, '12'))

  IF(PRESENT(SYM)) THEN

     IF(SYM .EQ. 'R') THEN

        !! Absorb iLambdas from East to West
        CALL LambdaMUL(imps(sA) % m, iLambda(sA) % m, '3')
        CALL LambdaMUL(imps(sB) % m, iLambda(sB) % m, '3')

     ELSEIF(SYM .EQ. 'L') THEN

        !! Absorb iLambdas from West to East
        CALL LambdaMUL(imps(sA) % m, iLambda(sB) % m, '2')
        CALL LambdaMUL(imps(sB) % m, iLambda(sA) % m, '2')

     ELSEIF(SYM .EQ. 'ALL') THEN

        !! Construct iL(B)*iG(A)*iL(A)*iG(B)*iL(B)
        CALL LambdaMUL(imps(sA) % m, iLambda(sA) % m, imps(sB) % m)
        CALL LambdaMUL(imps(sA) % m, iLambda(sB) % m, '2')
        CALL LambdaMUL(imps(sB) % m, iLambda(sB) % m, '3')

     ELSE
        CALL invalid_flag("absorb_imps_lambdas -- invalid SYM ", SYM)
     end IF

  ELSE
     !! Absorb iLambda(A) && iLambda(B) into imps(A,B) -- construct a symmetric combo
     CALL LambdaMUL(imps(sA) % m, iLambda(sA) % m, imps(sB) % m)
     CALL LambdaMUL(imps(sB) % m, iLambda(sB) % m, imps(sA) % m)
  end IF

 END SUBROUTINE absorb_imps_lambdas



 !!! Absorb iMPO lambdas into gammas !!!
 SUBROUTINE absorb_impo_lambdas(impo, iGamma, iLambda, SITEORD, SYM)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT)        :: impo(:)
  TYPE(block_mpo),              INTENT(IN), OPTIONAL :: iGamma(:)
  TYPE(block_lambda),           INTENT(IN)           :: iLambda(:)
  CHARACTER(LEN=*),             INTENT(IN), OPTIONAL :: SITEORD
  CHARACTER(LEN=*),             INTENT(IN), OPTIONAL :: SYM

  !! iMPO sites
  INTEGER :: sA, sB

  !! Copy iGamma to iMPO block (if iGamma not present, use iMPO as input gamma instead)
  IF(PRESENT(iGamma)) CALL copy_mpo_block(impo, iGamma) 

  !! Determine siteA, siteB
  CALL get_impo_sites(sA, sB, OptArg(SITEORD, '12'))

  IF(PRESENT(SYM)) THEN

     IF(SYM .EQ. 'R') THEN

        !! Absorb iLambdas from East to West
        CALL LambdaMUL(impo(sA) % m, iLambda(sA) % m, '4')
        CALL LambdaMUL(impo(sB) % m, iLambda(sB) % m, '4')

     ELSEIF(SYM .EQ. 'L') THEN

        !! Absorb iLambdas from West to East
        CALL LambdaMUL(impo(sA) % m, iLambda(sB) % m, '3')
        CALL LambdaMUL(impo(sB) % m, iLambda(sA) % m, '3')

     ELSEIF(SYM .EQ. 'ALL') THEN

        !! Construct iL(B)*iG(A)*iL(A)*iG(B)*iL(B)
        CALL LambdaMUL(impo(sA) % m, iLambda(sA) % m, impo(sB) % m)
        CALL LambdaMUL(impo(sA) % m, iLambda(sB) % m, '3')
        CALL LambdaMUL(impo(sB) % m, iLambda(sB) % m, '4')

     ELSE
        CALL invalid_flag("absorb_impo_lambdas -- invalid SYM ", SYM)
     end IF

  ELSE
     !! Absorb iLambda(A) && iLambda(B) into impo(A,B) -- construct a symmetric combo
     CALL LambdaMUL(impo(sA) % m, iLambda(sA) % m, impo(sB) % m)
     CALL LambdaMUL(impo(sB) % m, iLambda(sB) % m, impo(sA) % m)
  end IF

 END SUBROUTINE absorb_impo_lambdas


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TRANSPOSE OR INVERT MPS/MPO BLOCKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Take a transpose bi-MPO block (swap SN bonds, but keep EW bonds the same) !!!
 SUBROUTINE transpose_bimpo_block(mpo_transp, mpo_block) 

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: mpo_transp(:,:)
  TYPE(block_mpo),              INTENT(IN)    :: mpo_block(:,:)

  !! Temp mpo site for copying
  TYPE(block_mpo), ALLOCATABLE :: mpoTmp1(:), mpoTmp2(:)

  !! Transpose each component of binary MPO (swap block indices!)
  CALL transpose_mpo_block(mpoTmp1, mpo_block(2,:))
  CALL transpose_mpo_block(mpoTmp2, mpo_block(1,:))

  !! Combine two mpoTmp1, mpoTmp2 into binary-mpo = mpo_transp
  CALL combine_two_mpos_into_bimpo(mpo_transp, mpoTmp1, mpoTmp2)

 END SUBROUTINE transpose_bimpo_block






 !!! Take a transpose of mpo_block (swap SN bonds, but keep EW bonds the same) !!!
 SUBROUTINE transpose_mpo_block(mpo_transp, mpo_block, HC) 

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT)          :: mpo_transp(:)
  TYPE(block_mpo),              INTENT(IN)             :: mpo_block(:)
  CHARACTER(LEN=*),             INTENT(IN),   OPTIONAL :: HC

  !! Dims & indices
  INTEGER :: N_sites, site

  !! Allocate empty block
  N_sites = SIZE(mpo_block) 
  CALL allocate_empty_mpo_block(mpo_transp, N_sites)

  IF (.NOT. PRESENT(HC)) THEN

     !! Loop over mpo_block, transpose each site through '12' axes (i.e. swap South & North)
     DO site=1,N_sites
        CALL allocate_mpo_site(mpo_transp(site),        TensTRANSPOSE(mpo_block(site) % m, '12'))
     end DO
  ELSE

     !! Loop over mpo_block, transpose each site through '12' axes (i.e. swap South & North)
     DO site=1,N_sites
        CALL allocate_mpo_site(mpo_transp(site),  CONJG(TensTRANSPOSE(mpo_block(site) % m, '12')))
     end DO
  end IF

 END SUBROUTINE transpose_mpo_block





 !!! Invert MPO block 
 !!! (swap beginning and end + swap direction of EW bonds, but keep SN bonds the same) 
 SUBROUTINE invert_mpo_block(mpo_block, mpo_out)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT)           :: mpo_block(:)
  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: mpo_out(:)

  !! Temp MPO
  TYPE(block_mpo), ALLOCATABLE :: mpo_inverted(:)

  !! Indices
  INTEGER :: N_sites, site

  !! Allocate empty mpo_inverted
  N_sites = SIZE(mpo_block)
  CALL allocate_empty_mpo_block(mpo_inverted, N_sites)
   
  !! Swap mpo(site) and mpo(N_sites - site + 1), transpose each site through '34' axes
  DO site=1,N_sites
     CALL allocate_mpo_site(mpo_inverted(site), TensTRANSPOSE(mpo_block(N_sites-site+1) % m, '34'))
  end DO

  !! Write to output
  IF(PRESENT(mpo_out)) THEN
     CALL copy_mpo_block(mpo_out,   mpo_inverted)
  ELSE
     CALL copy_mpo_block(mpo_block, mpo_inverted)
  end IF

 END SUBROUTINE invert_mpo_block




 !!! Invert MPS block 
 !!! (swap beginning and end + swap direction of EW bonds, but keep SN bonds the same)
 SUBROUTINE invert_mps_block(mps_block, mps_out)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT)           :: mps_block(:)
  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: mps_out(:)

  !! Temp MPS
  TYPE(block_mps), ALLOCATABLE :: mps_inverted(:)

  !! Indices
  INTEGER :: N_sites, site

  !! Allocate empty mps_inverted 
  N_sites = SIZE(mps_block)
  CALL allocate_empty_mps_block(mps_inverted, N_sites)
   
  !! Swap mps(site) and mps(N_sites - site + 1), transpose each site through '23' axes
  DO site=1,N_sites
     CALL allocate_mps_site(mps_inverted(site), TensTRANSPOSE(mps_block(N_sites-site+1) % m, '23'))
  end DO

  !! Write to output 
  IF(PRESENT(mps_out)) THEN
     CALL copy_mps_block(mps_out,   mps_inverted)
  ELSE
     CALL copy_mps_block(mps_block, mps_inverted)
  end IF

 END SUBROUTINE invert_mps_block




 !!! Invert Lambda block !!!
 SUBROUTINE invert_lambda_block(lambda_block, lambda_out)

  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)           :: lambda_block(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: lambda_out(:)

  !! Temp inverted lambda block
  TYPE(block_lambda), ALLOCATABLE :: lambdaInverted(:) 

  !! Indices
  INTEGER :: N_sites, site

  !! Allocate empty lambda_inverted
  N_sites = SIZE(lambda_block)
  CALL allocate_empty_lambda_block(lambdaInverted, N_sites)
   
  !! Swap Lambda(site) and Lambda(N_sites - site + 1) -- no need to transpose cause Lambda = diagonal 
  DO site=1,N_sites
     CALL allocate_lambda_site(lambdaInverted(site),  lambda_block(N_sites-site+1) % m)
  end DO

  !! Write to output
  IF(PRESENT(lambda_out)) THEN
     CALL copy_lambda_block(lambda_out,   lambdaInverted)
  ELSE
     CALL copy_lambda_block(lambda_block, lambdaInverted)
  end IF

 END SUBROUTINE invert_lambda_block
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RESCALE MPS/MPO BLOCKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Multiply MPO by a given value !!!
 SUBROUTINE mult_mpo_by_value(mpo_block, input_value)
 
  TYPE(block_mpo),  INTENT(INOUT) :: mpo_block(:)
  COMPLEX(KIND=DP), INTENT(IN)    :: input_value

  !Dims & indices
  INTEGER :: N_sites, site

  !Size of MPO
  N_sites = SIZE(mpo_block)
  
  !Divide MPO by value going site-by-site
  Do site=1,N_sites
     mpo_block(site) % m = input_value**(1.0/REAL(N_sites, KIND=DP)) * mpo_block(site) % m
  end Do

 END SUBROUTINE mult_mpo_by_value



 !!! Multiply binary MPO by a given value !!!
 SUBROUTINE mult_bimpo_by_value(mpo_block, input_value)
 
  TYPE(block_mpo),  INTENT(INOUT) :: mpo_block(:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: input_value

  !! Dims & indices
  INTEGER :: N_sites, site, blk

  !! Size of binary MPO: 2 x N_sites
  N_sites = SIZE(mpo_block, 2)
  
  !! Divide MPO by value going site-by-site
  Do blk=1,2
    Do site=1,N_sites
       mpo_block(blk,site) % m = input_value**(1.0/REAL(2.0d0*N_sites, KIND=DP)) * mpo_block(blk,site) % m
    end Do
  end Do

 END SUBROUTINE mult_bimpo_by_value




 !!! Multiply mps by a given value !!!
 SUBROUTINE mult_mps_by_value(mps_block, input_value)
 
  TYPE(block_mps),  INTENT(INOUT) :: mps_block(:)
  COMPLEX(KIND=DP), INTENT(IN)    :: input_value

  !! Dims & indices
  INTEGER :: N_sites, site

  !! Size of mps_block
  N_sites = SIZE(mps_block)
  
  !! Divide evec by value going site-by-site
  Do site=1,N_sites
     mps_block(site) % m = input_value**(1.0/REAL(N_sites, KIND=DP)) * mps_block(site) % m
  end Do

 END SUBROUTINE mult_mps_by_value






 !!! Rescale MPO block my max ABS element !!!
 SUBROUTINE rescale_mpo_by_max_element(mpo_block)
 
  TYPE(block_mpo), INTENT(INOUT) :: mpo_block(:)

  !! Dims & indices
  INTEGER :: N_sites, site

  !! Maximum abs element in mpo_block
  REAL(KIND=DP), ALLOCATABLE :: maxEls(:)
  REAL(KIND=DP)              :: maxEl

  !! Size of mpo_block
  N_sites = SIZE(mpo_block)

  !! (1) Find max element on each mpo site
  ALLOCATE(maxEls(N_sites))
  Do site=1,N_sites
     maxEls(site) = MAXVAL(ABS(mpo_block(site) % m))
  end Do

  !! (2) Find max element in mpo_block
  maxEl = MAXVAL(maxEls)
   
  !! (3) Rescale by [max element] 
  CALL mult_mpo_by_value(mpo_block, 1.0d0/complx(maxEl))

 END SUBROUTINE rescale_mpo_by_max_element






 !!! Rescale MPS block my max ABS element !!!
 SUBROUTINE rescale_mps_by_max_element(mps_block)
 
  TYPE(block_mps), INTENT(INOUT) :: mps_block(:)

  !! Dims & indices
  INTEGER :: N_sites, site

  !! Maximum abs element in mps_block
  REAL(KIND=DP), ALLOCATABLE :: maxEls(:)
  REAL(KIND=DP)              :: maxEl

  !! Size of mps_block
  N_sites = SIZE(mps_block)

  !! (1) Find max element on each mps site
  ALLOCATE(maxEls(N_sites))
  Do site=1,N_sites
     maxEls(site) = MAXVAL(ABS(mps_block(site) % m))
  end Do

  !! (2) Find max element in mps_block
  maxEl = MAXVAL(maxEls)
   
  !! (3) Rescale by [max element] 
  CALL mult_mps_by_value(mps_block, 1.0d0/complx(maxEl))

 END SUBROUTINE rescale_mps_by_max_element




 !!! Rescale lambda so that Sum(iLambda**2) = 1 !!!
 SUBROUTINE rescale_lambda(iLambdaSite)

  TYPE(block_lambda), INTENT(INOUT) :: iLambdaSite
  COMPLEX(KIND=DP)                  :: fac

  fac = sqrt(SUM(iLambdaSite % m ** 2))

  !! Normalize lambda so that Sum(iLambda**2) = 1
  iLambdaSite % m = iLambdaSite % m / fac
  
 END SUBROUTINE rescale_lambda




 !!! Symmetrize iMPS so that both sites have the same MAX element !!!
 SUBROUTINE symmetrize_imps_sites(imps)
 
  TYPE(block_mps), INTENT(INOUT) :: imps(:)

  !! Dims & indices
  INTEGER :: site

  !! Maximum abs element in mps_block
  REAL(KIND=DP) :: maxEls(2), fac

  !! Find max elements on each imps site
  Do site=1,2
     maxEls(site) = MAXVAL(ABS(imps(site) % m))
  end Do

  !! Rescale imps so that both sites have the same MAX element
  fac = sqrt(ABS(maxEls(2) / maxEls(1)))
  imps(1) % m = imps(1) % m * fac
  imps(2) % m = imps(2) % m / fac

 END SUBROUTINE symmetrize_imps_sites





 !!! Equalize iMPS bond dims !!!
 SUBROUTINE equalize_imps_bonds(iGamma, iLambda)

  TYPE(block_mps),    ALLOCATABLE,  INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE,  INTENT(INOUT) :: iLambda(:)

  !! Dims && indices
  INTEGER :: chi(2), min_chi
  INTEGER :: site

  !! Get bond dims
  DO site=1,2
     chi(site) = iLambda(site) % ChiDim 
  end DO

  !! Check if already equal
  IF(chi(1) .EQ. chi(2)) RETURN

  !! Get min bond dim
  min_chi = MIN(chi(1), chi(2))

  !! Equalize bonds
  DO site=1,2
     CALL allocate_mps_site(iGamma(site),     iGamma(site)  % m(:, 1:min_chi, 1:min_chi)) 
     CALL allocate_lambda_site(iLambda(site), iLambda(site) % m(1:min_chi))
  end DO

 END SUBROUTINE equalize_imps_bonds




 !!! Equalize iMPO bond dims !!!
 SUBROUTINE equalize_impo_bonds(iGamma, iLambda)

  TYPE(block_mpo),    ALLOCATABLE,  INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE,  INTENT(INOUT) :: iLambda(:)

  !! Dims && indices
  INTEGER :: chi(2), min_chi
  INTEGER :: site

  !! Get bond dims
  DO site=1,2
     chi(site) = iLambda(site) % ChiDim 
  end DO

  !! Check if already equal
  IF(chi(1) .EQ. chi(2)) RETURN

  !! Get min bond dim
  min_chi = MIN(chi(1), chi(2))

  !! Equalize bonds
  DO site=1,2
     CALL allocate_mpo_site(iGamma(site),     iGamma(site)  % m(:, :, 1:min_chi, 1:min_chi)) 
     CALL allocate_lambda_site(iLambda(site), iLambda(site) % m(1:min_chi))
  end DO

 END SUBROUTINE equalize_impo_bonds

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MANIPULATING MPS/MPO SITES: FACTORIZE, RESHAPE, SWAP, etc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! LHS-SVD factorization theta = U * [UH * theta] !!!
 SUBROUTINE LHS_factorize_Tens3D(tensA_REDUCED, X, tensA, facFlag)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensA_REDUCED(:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: X(:,:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:)
  CHARACTER(LEN=*),              INTENT(IN)    :: facFlag

  !! Dims
  INTEGER :: dimA(3)

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:),  UH(:,:), theta(:,:)
  INTEGER                       :: chi,  sigma_dim

  !! Get dims of tensA
  dimA = shape(tensA)

  !! (1) Construct theta
  theta = RESHAPE_2D(tensA, facFlag)

  !! (2) Compute SVD
  chi = -1
  CALL compute_lapack_svd(theta, U, chi, 2.0d0)
  UH = TensTRANSPOSE(U, 'HC')

  !! (3) Copy back SVD results
  SELECT CASE(facFlag)
  CASE('2,13')

       CALL copyTens(X, U)
       tensA_REDUCED = RESHAPE_3D(TENSMUL(UH, theta), facFlag, (/ dimA(1), dimA(3)/))

  CASE('12,3')

       tensA_REDUCED = RESHAPE_3D(U,  facFlag, (/ dimA(1), dimA(2) /))
       X             = TENSMUL(UH, theta)

  CASE DEFAULT
       CALL invalid_flag("LHS_factorize_Tens3D -- invalid facFlag ", facFlag)
  end SELECT

 END SUBROUTINE LHS_factorize_Tens3D





 !!! RHS-SVD factorization theta = [theta * V] * VH !!!
 SUBROUTINE RHS_factorize_Tens3D(tensA_REDUCED, X, tensA, facFlag)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensA_REDUCED(:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: X(:,:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:)
  CHARACTER(LEN=*),              INTENT(IN)    :: facFlag

  !! Dims
  INTEGER :: dimA(3)

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: VH(:,:), V(:,:), U(:,:), theta(:,:)
  INTEGER                       :: chi, sigma_dim

  !! Get dims of tensA
  dimA = shape(tensA)

  !! (1) Construct theta
  theta = RESHAPE_2D(tensA, facFlag)

  !! (2) Compute SVD
  chi = -1
  CALL compute_lapack_svd(TensTRANSPOSE(theta), U, chi, 2.0d0)
  VH = TensTRANSPOSE(U)
  V  = TensTRANSPOSE(VH, 'HC')

  !! (3) Copy back SVD results
  SELECT CASE(facFlag)
  CASE('2,13')

       X             = TENSMUL(theta, V)
       tensA_REDUCED = RESHAPE_3D(VH, facFlag, (/ dimA(1), dimA(3) /))

  CASE('12,3')

       tensA_REDUCED = RESHAPE_3D(TENSMUL(theta, V),  facFlag, (/ dimA(1), dimA(2) /))
       CALL copyTens(X, VH)

  CASE DEFAULT
       CALL invalid_flag("RHS_factorize_Tens3D -- invalid facFlag ", facFlag)
  end SELECT

 END SUBROUTINE RHS_factorize_Tens3D





 !!! LHS-SVD factorization theta = U * [UH * theta] !!!
 SUBROUTINE LHS_factorize_Tens4D(tensA_REDUCED, X, tensA, facFlag)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensA_REDUCED(:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: X(:,:,:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:,:)
  CHARACTER(LEN=*),              INTENT(IN)    :: facFlag

  !! Dims
  INTEGER :: dimA(4)

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:),  UH(:,:), theta(:,:)
  INTEGER                       :: chi,  sigma_dim

  !! Get dims of tensA
  dimA = shape(tensA)

  !! (1) Construct theta
  theta = RESHAPE_2D(tensA, facFlag)

  !! (2) Compute SVD
  chi = -1
  CALL compute_lapack_svd(theta, U, chi, 2.0d0)
  UH = TensTRANSPOSE(U, 'HC')

  !! (3) Copy back SVD results
  SELECT CASE(facFlag)
  CASE('23,14')

       CALL copyTens(X, RESHAPE_3D(U, '12,3', (/ dimA(2), dimA(3) /)))
       tensA_REDUCED  = RESHAPE_3D(TENSMUL(UH, theta), '2,13', (/ dimA(1), dimA(4) /))

  CASE('13,24')

       tensA_REDUCED = RESHAPE_3D(U, '12,3', (/ dimA(1), dimA(3) /))
       X             = RESHAPE_3D(TENSMUL(UH, theta), '2,13', (/ dimA(2), dimA(4) /))

  CASE DEFAULT
       CALL invalid_flag("LHS_factorize_Tens4D -- invalid facFlag ", facFlag)
  end SELECT

 END SUBROUTINE LHS_factorize_Tens4D





 !!! RHS-SVD factorization theta = [theta * V] * VH !!!
 SUBROUTINE RHS_factorize_Tens4D(tensA_REDUCED, X, tensA, facFlag)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensA_REDUCED(:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: X(:,:,:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:,:)
  CHARACTER(LEN=*),              INTENT(IN)    :: facFlag

  !! Dims
  INTEGER :: dimA(4)

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: VH(:,:), V(:,:), U(:,:), theta(:,:)
  INTEGER                       :: chi, sigma_dim

  !! Get dims of tensA
  dimA = shape(tensA)

  !! (1) Construct theta
  theta = RESHAPE_2D(tensA, facFlag)

  !! (2) Compute SVD
  chi = -1
  CALL compute_lapack_svd(TensTRANSPOSE(theta), U, chi, 2.0d0)
  VH = TensTRANSPOSE(U)
  V  = TensTRANSPOSE(VH, 'HC')

  !! (3) Copy back SVD results
  SELECT CASE(facFlag)
  CASE('23,14')

       X             = RESHAPE_3D(TENSMUL(theta, V), '12,3', (/ dimA(2), dimA(3) /))
       tensA_REDUCED = RESHAPE_3D(VH, '2,13', (/ dimA(1), dimA(4) /))

  CASE('13,24')

       tensA_REDUCED  = RESHAPE_3D(TENSMUL(theta, V), '12,3', (/ dimA(1), dimA(3) /))
       CALL copyTens(X, RESHAPE_3D(VH, '2,13', (/ dimA(2), dimA(4) /)))

  CASE DEFAULT
       CALL invalid_flag("RHS_factorize_Tens4D -- invalid facFlag ", facFlag)
  end SELECT

 END SUBROUTINE RHS_factorize_Tens4D





 !!! Symmetric SVD factorization theta = [mps(1) -- mps(2)] !!!
 SUBROUTINE factorize_Tens4D(mps, theta, dimT, chi0, eps)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps(:)
  COMPLEX(KIND=DP),             INTENT(IN)    :: theta(:,:,:,:)
  INTEGER,                      INTENT(IN)    :: dimT(4)
  INTEGER,                      INTENT(IN)    :: chi0
  REAL(KIND=DP),                INTENT(IN)    :: eps

  !! SVD matrices
  COMPLEX(KIND=DP), ALLOCATABLE :: U(:,:), VH(:,:)
  REAL(KIND=DP),    ALLOCATABLE :: Sigma(:)
  INTEGER                       :: chi

  !! Allocate MPS
  CALL allocate_empty_mps_block(mps, 2)

  !! Factorize theta using SVD
  chi = chi0
  CALL compute_lapack_svd(RESHAPE_2D(theta, '13,24'), U, VH, Sigma, chi, eps)

  !! Copy SVD results to MPS sites
  CALL LambdaMUL(U, complx(Sigma), VH)
  CALL allocate_mps_site(mps(1),  RESHAPE_3D(U,  '12,3',  (/ dimT(1), dimT(3) /)))
  CALL allocate_mps_site(mps(2),  RESHAPE_3D(VH, '2,13',  (/ dimT(2), dimT(4) /)))

 END SUBROUTINE factorize_Tens4D









 !!! Transform TN sites (e.g. two-site MPS) into a vector !!!
 SUBROUTINE convert_TN_to_vec(vec, dimTN, mps)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: vec(:)
  INTEGER,                       INTENT(INOUT) :: dimTN(2,3)
  TYPE(block_mps),               INTENT(IN)    :: mps(:)

  !! Reshape MPS sites into vecs, concatenate them into a single vector
  CALL concatenate_vecs(vec, RESHAPE_1D(mps(1) % m), RESHAPE_1D(mps(2) % m))

  !! Record dims of MPS sites so we can reshape back
  dimTN(1,:) = shape(mps(1) % m)
  dimTN(2,:) = shape(mps(2) % m)

 END SUBROUTINE convert_TN_to_vec






 !!! Transform TN sites (e.g. two-site MPS) into a vector !!!
 SUBROUTINE convert_vec_to_TN(mps, vec, dimTN)

  TYPE(block_mps),   ALLOCATABLE, INTENT(INOUT) :: mps(:)
  COMPLEX(KIND=DP),               INTENT(IN)    :: vec(:)
  INTEGER,                        INTENT(IN)    :: dimTN(2,3)

  COMPLEX(KIND=DP), ALLOCATABLE :: vecA(:), vecB(:)

  !! Decatenate a two-component vec into vecA, vecB
  CALL decatenate_vecs(vecA, vecB, vec, product(dimTN(1,:)), product(dimTN(2,:)))

  !! Reshape vecA, vecB into MPS sites using dimTN provided
  CALL allocate_empty_mps_block(mps, 2)
  CALL allocate_mps_site(mps(1), RESHAPE_3D(vecA, dimTN(1,:)))
  CALL allocate_mps_site(mps(2), RESHAPE_3D(vecB, dimTN(2,:)))

 END SUBROUTINE convert_vec_to_TN





 !!! Resize MPS site using the input SNdim !!!
 SUBROUTINE resize_mps_site(mps, SNdim)

  TYPE(block_mps), INTENT(INOUT) :: mps
  INTEGER,         INTENT(IN)    :: SNdim

  COMPLEX(KIND=DP), ALLOCATABLE :: tmpSite(:,:,:)

  CALL copyTens(tmpSite, mps % m, (/SNdim, mps % Wdim, mps % Edim/))
  CALL allocate_mps_site(mps, tmpSite)

 END SUBROUTINE resize_mps_site





 !!! Swap two MPS sites !!!
 SUBROUTINE swap_mps_sites(mpsA, mpsB)
 
  TYPE(block_mps), INTENT(INOUT) :: mpsA, mpsB

  !! Temp storage for mpsA
  COMPLEX(KIND=DP), ALLOCATABLE :: tmpA(:,:,:)

  !! Copy mpsA to temp storage 
  CALL copyTens(tmpA, mpsA % m)

  !! Swap mpsA && mpsB
  CALL allocate_mps_site(mpsA, mpsB % m)
  CALL allocate_mps_site(mpsB, tmpA)

 END SUBROUTINE swap_mps_sites




 !!! Swap two MPO sites !!!
 SUBROUTINE swap_mpo_sites(mpoA, mpoB)
 
  TYPE(block_mpo), INTENT(INOUT) :: mpoA, mpoB

  !! Temp storage for mpoA
  COMPLEX(KIND=DP), ALLOCATABLE :: tmpA(:,:,:,:)

  !! Copy mpoA to temp storage 
  CALL copyTens(tmpA, mpoA % m)

  !! Swap mpoA && mpoB
  CALL allocate_mpo_site(mpoA, mpoB % m)
  CALL allocate_mpo_site(mpoB, tmpA)

 END SUBROUTINE swap_mpo_sites

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RESHAPE MPS <-> MPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Reshape MPS to MPO !!!
 SUBROUTINE reshape_MPS_to_MPO(mpo, mps, SNdims)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: mpo(:)
  TYPE(block_mps),              INTENT(IN)    :: mps(:)
  INTEGER,                      INTENT(IN)    :: SNdims(:)

  INTEGER :: site, N_sites
  INTEGER :: localDim
  INTEGER :: Sdim, Ndim
  
  !! Block size
  N_sites = SIZE(mps)

  !! Check dims
  CALL check_sizes_equal(product(SNdims), mps(1) % SNdim, "reshape_MPS_to_MPO -- mps && mpo dims must be consistent: ")

  !! Allocate empty MPO
  CALL allocate_empty_mpo_block(mpo, N_sites)

  !! Calc MPO sites
  DO site=1,N_sites
     CALL allocate_mpo_site(mpo(site),  RESHAPE_4D(mps(site) % m, '12,3,4', SNdims))
  end DO

 END SUBROUTINE reshape_MPS_to_MPO





 !!! Reshape MPO to MPS !!!
 SUBROUTINE reshape_MPO_to_MPS(mps, mpo)

  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps(:)
  TYPE(block_mpo),              INTENT(IN)    :: mpo(:)

  INTEGER :: site, N_sites
  
  !! Block size
  N_sites = SIZE(mpo)

  !! Allocate empty MPS
  CALL allocate_empty_mps_block(mps, N_sites)

  !! Calc MPS sites
  DO site=1,N_sites
     CALL allocate_mps_site(mps(site),  RESHAPE_3D(mpo(site) % m, '12,3,4'))
  end DO

 END SUBROUTINE reshape_MPO_to_MPS


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE mps_mpo_utility
