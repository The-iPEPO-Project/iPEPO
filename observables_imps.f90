MODULE observables_imps

 USE utility
 USE definitions_mps_mpo
 USE mps_mpo_utility

 IMPLICIT NONE

CONTAINS


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPUTING OBSERVABLES DURING/AFTER TIME EVOLUTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute 1P observable !!!
 SUBROUTINE compute_obs_1P(expval, rho_onesite, op)

  COMPLEX(KIND=DP), INTENT(INOUT) :: expval
  COMPLEX(KIND=DP), INTENT(IN)    :: rho_onesite(:), op(:,:)

  !! Compute 1P observable
  expval = SUM(op(1,:) * rho_onesite)

 END SUBROUTINE compute_obs_1P




 !!! Compute 2P observable !!!
 SUBROUTINE compute_obs_2P(expvals, rho_twosite, rho_onesite, ops, min_sep, max_sep)

  COMPLEX(KIND=DP), INTENT(INOUT) :: expvals(:)
  COMPLEX(KIND=DP), INTENT(IN)    :: rho_twosite(:,:,:), rho_onesite(:)
  COMPLEX(KIND=DP), INTENT(IN)    :: ops(:,:,:)
  INTEGER,          INTENT(IN)    :: min_sep, max_sep

  !! Indices
  INTEGER :: sep

  !! Compute 2P observables (for different separations) -- sep=0 handled separately as onesite case
  DO sep=min_sep,max_sep
     IF(sep .EQ. 0) THEN 
        expvals(sep+1-min_sep) = SUM(ops(1,1,:) * TENSMUL(ops(2,:,:), rho_onesite)) !NB. we mult in order: Tr(O1*O2*rho)
     ELSE
        expvals(sep+1-min_sep) = SUM(ops(1,1,:) * TENSMUL(rho_twosite(sep,:,:), ops(2,1,:)))
     end IF
  end DO

 END SUBROUTINE compute_obs_2P



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! COMPUTING OBSERVABLES FOR iMPS REDUCED TO iMPS -- A TEST CASE FOR TIME EVOLUTION !!!!!!!!!!!!!!!!!!!!!!!!

 !!! Compute boundaries of infinite MPS !!!
 SUBROUTINE compute_rho_imps_bounds(eval, Rvec, Lvec, iGamma, iLambda)

  COMPLEX(KIND=DP),              INTENT(INOUT) :: eval
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: Rvec(:), Lvec(:)
  TYPE(block_mps),               INTENT(IN)    :: iGamma(:)  !iGamma of MPS
  TYPE(block_lambda),            INTENT(IN)    :: iLambda(:) !iLambda of MPS

  !Local iMPS copy
  TYPE(block_mps), ALLOCATABLE :: imps(:) 

  !Create nolambda imps
  CALL absorb_imps_lambdas(imps, iGamma, iLambda, '12', 'R')

  !Get dominant eval && evecs of iG1-iL1-iG2-iL2 block
  CALL dominant_eig(eval, Rvec, Lvec, TENSMUL(imps(1) % m(1,:,:), imps(2) % m(1,:,:)))

 END SUBROUTINE compute_rho_imps_bounds





 !!! Construct reduced onesite && twosite rho, and dominant eval of TN unit cell !!!
 SUBROUTINE construct_rho_onesite(rho_onesite, iGamma, iLambda, Rvec, Lvec) 

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: rho_onesite(:) !onesite reduced rho
  TYPE(block_mps),               INTENT(IN)    :: iGamma(:)      !iGamma of MPS
  TYPE(block_lambda),            INTENT(IN)    :: iLambda(:)     !iLambda of MPS
  COMPLEX(KIND=DP),              INTENT(IN)    :: Rvec(:), Lvec(:)

  !Local iMPS copy
  TYPE(block_mps), ALLOCATABLE :: imps(:)

  !dims && indices
  INTEGER :: i, LocalDim

  !Construct nolambda iMPS
  CALL absorb_imps_lambdas(imps, iGamma, iLambda, '12', 'R')

  !Get LocalDim, Allocate rho_onesite
  LocalDim = iGamma(1)%SNdim
  CALL allocateTens(rho_onesite, (/LocalDim/))

  !Compute individual elements of rho_onesite one by one
  DO i=1,LocalDim
     rho_onesite(i) = SUM(Lvec * TENSMUL(TENSMUL(imps(1)%m(i,:,:), imps(2)%m(1,:,:)), Rvec))/SUM(Lvec*Rvec)
  end DO

 END SUBROUTINE construct_rho_onesite





 !!! Construct reduced onesite && twosite rho, and dominant eval of TN unit cell !!!
 SUBROUTINE construct_rho_twosite(rho_twosite, iGamma, iLambda, Rvec, Lvec, max_sep) 

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: rho_twosite(:,:,:) !twosite reduced rho
  TYPE(block_mps),               INTENT(IN)    :: iGamma(:)          !iGamma of MPS
  TYPE(block_lambda),            INTENT(IN)    :: iLambda(:)         !iLambda of MPS
  COMPLEX(KIND=DP),              INTENT(IN)    :: Rvec(:), Lvec(:)   !Boundary vecs
  INTEGER,                       INTENT(IN)    :: max_sep

  !! Local iMPS copy
  TYPE(block_mps), ALLOCATABLE :: imps(:)

  !! Tensors for calculating two-site reduced rho 
  COMPLEX(KIND=DP), ALLOCATABLE :: tensBulk(:,:) 
  COMPLEX(KIND=DP), ALLOCATABLE :: padding(:,:), lmat(:,:), rmat(:,:)

  !! Dims && indices
  INTEGER :: i, j, LocalDim, bdim
  INTEGER :: n, sep

  !! Absorb iMPS lambdas
  CALL absorb_imps_lambdas(imps, iGamma, iLambda, '12', 'R')

  !! Get dims, Allocate rho_twosite (use bdim = imps(2)%Edim since this is the relevant matrix dim)
  bdim = imps(2) % Edim
  LocalDim = imps(1) % SNdim
  CALL allocateTens(rho_twosite, (/max_sep, LocalDim, LocalDim/))

  !! Allocate all intermediate tensors 
  ALLOCATE(tensBulk(bdim,bdim), padding(bdim,bdim), lmat(bdim,bdim), rmat(bdim,bdim))

  !! Create two-site rho
  DO sep=1,max_sep

     !! Create padding between operator sites 
     padding = matEye(bdim)
     DO n=1,MAX((sep/2)-1, 0)
        padding = TENSMUL(padding, TENSMUL(imps(1)%m(1,:,:), imps(2)%m(1,:,:)))
     end DO

     !! Create two-site rho
     DO i=1,LocalDim
       DO j=1,LocalDim

          IF(sep .EQ. 1) THEN
             !! Special case: two-site block has both sites
             tensBulk = TENSMUL(imps(1)%m(i,:,:), imps(2)%m(j,:,:))
          ELSE
             !! Either two left or two right blocks, depending on parity. 
             !! If even separation we have a pattern ... [**] [i*] [**] [j*] [**] ...
             !! The above picture has separation 4). 
             !! If odd separation we have a pattern  ... [**] [i*] [**] [*j] [**] ... 
             !! The above picture has separation 5).
             !! Here * represents identity and the [] the two-site block.

             IF(isEVEN(sep)) THEN
                !! Even separation, both blocks have selected site left
                lmat= TENSMUL(imps(1)%m(i,:,:), imps(2)%m(1,:,:))
                rmat= TENSMUL(imps(1)%m(j,:,:), imps(2)%m(1,:,:))
             ELSE
                !! Odd separation, left has left, right has right.
                lmat= TENSMUL(imps(1)%m(i,:,:), imps(2)%m(1,:,:))
                rmat= TENSMUL(imps(1)%m(1,:,:), imps(2)%m(j,:,:))
             end IF

             !! Full Bulk Tensor
             tensBulk = TENSMUL(lmat, TENSMUL(padding, rmat))
          END IF
          
          !! Trace over the matrix product with projector.
          rho_twosite(sep,i,j) = SUM(Lvec * TENSMUL(tensBulk, Rvec))/SUM(Lvec*Rvec) 
       end DO
     end DO

  END DO

 END SUBROUTINE construct_rho_twosite

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALCULATE PURITY OF REDUCED RHO DURING TIME EVOLUTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !!! Calculate purity of reduced onesite rho !!!
 SUBROUTINE calc_purity_rho_onesite(purity, rho_onesite)

  COMPLEX(KIND=DP), INTENT(OUT) :: purity         !purity = Tr(rho**2)
  COMPLEX(KIND=DP), INTENT(IN)  :: rho_onesite(:) !onesite reduced rho

  !! Calculate purity
  CALL calc_purity_rho_vec(purity, rho_onesite)

 END SUBROUTINE calc_purity_rho_onesite




 !!! Calculate purity of reduced twosite rho !!!
 SUBROUTINE calc_purity_rho_twosite(purity, rho_twosite)

  COMPLEX(KIND=DP), INTENT(OUT) :: purity           !purity = Tr(rho**2)
  COMPLEX(KIND=DP), INTENT(IN)  :: rho_twosite(:,:) !twosite reduced rho

  COMPLEX(KIND=DP), ALLOCATABLE :: rho_vec(:)

  !! Combine two sites into a single vector
  rho_vec = RESHAPE_1D(rho_twosite)

  !! Calculate purity
  CALL calc_purity_rho_vec(purity, rho_vec)

 END SUBROUTINE calc_purity_rho_twosite




 !!! Compute purity of a given rho_vec !!!
 SUBROUTINE calc_purity_rho_vec(purity, rho_vec)

  COMPLEX(KIND=DP), INTENT(OUT) :: purity     !purity = Tr(rho**2)
  COMPLEX(KIND=DP), INTENT(IN)  :: rho_vec(:) !reduced rho

  INTEGER                       :: dimO(2), dimH
  REAL(KIND=DP)                 :: dimRho
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_mat(:,:)
  COMPLEX(KIND=DP)              :: norm

  !! Get Hilbert space dim
  dimRho = sqrt(1.0D0*SIZE(rho_vec))
  dimH = INT(dimRho)
  dimO(1) = dimH; dimO(2) = dimH

  !! Consistency check
  IF(PRODUCT(dimO) .NE. SIZE(rho_vec)) THEN
     WRITE(*,*) "calc_purity_rho_vec: wrong size of rho_vec ", dimH, dimRho 
     STOP
  end IF

  !! Reshape rho_vec into rho_mat
  rho_mat = RESHAPE_2D(rho_vec, dimO)

  !! Compute purity
  purity = TRACE(TENSMUL(rho_mat, rho_mat))

  !! Normalize purity
  norm = TRACE(rho_mat)
  purity = purity/norm

 END SUBROUTINE calc_purity_rho_vec

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE observables_imps
