MODULE matrix_functions

 USE basic_functions
 USE error_handling
 USE array_utility
 USE TENS_transpose
 USE TENS_mult
 USE TENS_mult_extension

 IMPLICIT NONE

CONTAINS

 !!! Diagonalize a matrix with LAPACK !!!
 SUBROUTINE diagonalize_matrix(evals, evecsR, evecsL, tensA_in)

  USE lapack_wrappers, ONLY: lapack_nsym_eigv_all, lapack_nsym_eig_all

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT)           :: evals(:)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT), OPTIONAL :: evecsR(:,:), evecsL(:,:)
  COMPLEX(KIND=DP),              INTENT(IN)            :: tensA_in(:,:)

  !! Local copy of the input tensA, output Evecs, and dim of tensA
  COMPLEX(KIND=DP), ALLOCATABLE :: tensA(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: R(:,:), L(:,:)
  INTEGER                       :: adim

  !! Verify square matrix
  CALL check_sizes_equal(SIZE(tensA_in,1), SIZE(tensA_in,2),  "diagonalize_matrix: tensA must be a square matrix ")

  !! Obtain dims
  adim = SIZE(tensA_in,1)

  !! Create a local copy of tensA_in
  CALL copyTens(tensA, tensA_in)

  !!!!!!!!!!!!!! Find evecs && evals of tensA !!!!!!!!!!!!!!!!!!!!!!!!
  IF(PRESENT(evecsL) .OR. PRESENT(evecsR)) THEN

     !! Allocate evals && evecs
     ALLOCATE(evals(adim), L(adim, adim), R(adim, adim))

     !! Filter out small values
     CALL filterTens(tensA)

     !! Eigendecomposition
     CALL lapack_nsym_eigv_all(tensA, evals, L, R) 

     !! According to lapack def: u^H A = lambda * u^H --> evecsL = contains {u} rather than {u^H}
     !! we want evecsL = contains {u^H}, so take HC of evecsL 
     L = TensTRANSPOSE(L, 'HC')
  ELSE

     !! Use lapack to compute evals of tensA
     ALLOCATE(evals(adim))
     CALL lapack_nsym_eig_all(tensA, evals) 
  end IF

  !! Copy to output
  IF(PRESENT(evecsL)) CALL copyTens(evecsL,  L)
  IF(PRESENT(evecsR)) CALL copyTens(evecsR,  R)
  
 END SUBROUTINE diagonalize_matrix






 !!! Find the dominant eval && evecs !!!
 SUBROUTINE dominant_eig(eval, Rvec, Lvec, tensA)

  COMPLEX(KIND=DP),              INTENT(OUT)           :: eval
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT), OPTIONAL :: Rvec(:), Lvec(:)
  COMPLEX(KIND=DP),              INTENT(IN)            :: tensA(:,:)

  !! Evals && Evecs after diagonalization
  COMPLEX(KIND=DP), ALLOCATABLE :: evecsR(:,:), evecsL(:,:), evals(:)
  INTEGER                       :: adim, locus(1), i
  COMPLEX(KIND=DP)              :: norm

  !! Rvec && Lvec
  IF(PRESENT(Rvec) .AND. PRESENT(Lvec)) THEN
     !! Diagonalize tensA - get evals && evecs
     CALL diagonalize_matrix(evals, evecsR, evecsL, tensA)
  ELSE
     !! Diagonalize tensA - get evals only
     CALL diagonalize_matrix(evals=evals, tensA_in=tensA)
  end IF

  !Find the largest eval
  locus = MAXLOC(ABS(evals))
  eval  = evals(locus(1))

  !Find the corresponding evecs
  IF(PRESENT(Rvec) .AND. PRESENT(Lvec)) THEN

     !Allocate evecs
     adim = SIZE(tensA,1)
     ALLOCATE(Rvec(adim), Lvec(adim))

     !Get dominant evecs
     Rvec = evecsR(:, locus(1))
     Lvec = evecsL(locus(1), :)

  end IF

 END SUBROUTINE dominant_eig




 !!! Orthodecompose a square matrix as M = V*V^H using eigendecomposition !!!
 SUBROUTINE orthodecompose_matrix(tensA, VR, VL)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT), OPTIONAL :: VR(:,:), VL(:,:) !output evec matrices
  COMPLEX(KIND=DP),              INTENT(IN)            :: tensA(:,:)       !input square matrix

  COMPLEX(KIND=DP), ALLOCATABLE :: evecsR(:,:), evecsL(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: evals(:)

  !Diagonalize tensA
  CALL diagonalize_matrix(evals, evecsR, evecsL, tensA)

  !Absorb evals into evecsR, evecsL (NB. We've tensA = evecsR * LambdaMat * evecsL = SUM(lambda * rvec * lvec))
  CALL LambdaMUL(evecsR, evals, evecsL)

  !Output VR if we want decomposition using Right Evecs, VL if we want decomposition using Left Evecs
  IF(PRESENT(VR)) CALL copyTens(VR, evecsR)
  IF(PRESENT(VL)) CALL copyTens(VL, evecsL)

 END SUBROUTINE orthodecompose_matrix






 !!! Find matrix inverse !!!
 FUNCTION invert_matrix(Xmat) result(Xinv) 

  USE lapack_wrappers, ONLY: lapack_gen_mat_inv

  COMPLEX(KIND=DP), INTENT(IN)  :: Xmat(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: Xinv(:,:)

  INTEGER :: dimX(2)

  dimX = shape(Xmat)

  !! Initialize Xinv
  ALLOCATE(Xinv(dimX(1), dimX(2)))
  Xinv = Xmat

  !! Compute Xinv using lapack
  CALL lapack_gen_mat_inv(Xinv)

 END FUNCTION invert_matrix





 !!! Find SQRT of matrix !!!
 FUNCTION SQRT_matrix(tensA) result(SQRTmat)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: SQRTmat(:,:)

  COMPLEX(KIND=DP), ALLOCATABLE :: D(:), V(:,:), Vinv(:,:)

  !! Diagonalize tensA 
  CALL diagonalize_matrix(evals=D, evecsR=V, tensA_in=tensA)

  !! Calc Vinv = V^-1
  Vinv = invert_matrix(V)

  !! Calc SQRT(tensA) = V * D^1/2 * V^-1
  CALL LambdaMUL(V, SQRT(D), Vinv)
  SQRTmat = TENSMUL(V, Vinv)

 END FUNCTION SQRT_matrix





 !!! Compute matrix exponential !!!
 FUNCTION matrix_exp_series(Mat, fac) result(MatExp)

  COMPLEX(KIND=DP), INTENT(IN)           :: Mat(:,:)
  COMPLEX(KIND=DP), INTENT(IN), OPTIONAL :: fac
  COMPLEX(KIND=DP), ALLOCATABLE          :: MatExp(:,:) 

  COMPLEX(KIND=DP)              :: fac0
  COMPLEX(KIND=DP), ALLOCATABLE :: eye(:,:), temp(:,:)
  INTEGER :: dimM(2), i

  fac0 = OptArg(fac, (1.0D0, 0.0D0))

  !! Generate eye matrix
  dimM = shape(Mat)
  eye  = matEye(dimM(1))

  !! Initialize matrix exponent with 1st order approx
  MatExp = eye + Mat*fac0

  !! Power loop
  temp = Mat*fac0
  DO i=2,2 !100
     temp   = TENSMUL(Mat*fac0, temp) / i
     MatExp = MatExp + temp
  end DO

 END FUNCTION matrix_exp_series




 


 !!! Compute matrix exponential !!!
 FUNCTION matrix_exp(Mat, fac) result(MatExp)

  COMPLEX(KIND=DP), INTENT(IN)           :: Mat(:,:)
  COMPLEX(KIND=DP), INTENT(IN), OPTIONAL :: fac
  COMPLEX(KIND=DP), ALLOCATABLE          :: MatExp(:,:)

  COMPLEX(KIND=DP)              :: tempVal
  COMPLEX(KIND=DP), ALLOCATABLE :: tempVec(:)
  COMPLEX(KIND=DP)              :: C_norm

  COMPLEX(KIND=DP), ALLOCATABLE :: evals(:), DiagMatExp(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: evecs(:,:), evecsInv(:,:), evecsL(:,:)
  COMPLEX(KIND=DP)              :: fac0

  !! Indices & dims
  INTEGER :: dimM, i

  fac0 = OptArg(fac, (1.0D0, 0.0D0))

  !! (1) Diagonalize input Mat
  CALL diagonalize_matrix(evals, evecs, evecsL, Mat)

  !! Filter out small values
  CALL filterTens(evecs)

  !! (2) Get inverted evec matrix
  evecsInv = invert_matrix(evecs)

  !! Filter out small values
  CALL filterTens(evals)
  CALL filterTens(evecs)
  CALL filterTens(evecsInv)

  !! (3) Exponentiate diagonal matrix
  dimM = SIZE(evals)
  ALLOCATE(DiagMatExp(dimM, dimM)); DiagMatExp(:,:) = (0.0D0, 0.0D0)
  DO i=1,dimM
     DiagMatExp(i,i) = exp(evals(i)*fac0)
  end DO

  !! Filter out small values
  CALL filterTens(DiagMatExp)

  !! (4) Combine Diagonal Exp Matrix with evec matrices --> Gives Matrix Exponential
  MatExp = TENSMUL(evecs, TENSMUL(DiagMatExp, evecsInv))

  !! Filter out small values
  CALL filterTens(MatExp)

  !DO i=1,dimM
     !C_norm = SQRT(SUM(evecs(:,i) * evecsL(i,:)))
     !evecs(:,i)  = evecs(:,i)  / C_norm
     !evecsL(i,:) = evecsL(i,:) / C_norm
     !evals(i) = evals(i) * C_norm**2
  !end DO

  !MatExp = TENSMUL(evecs, TENSMUL(DiagMatExp, evecsL))

 END FUNCTION matrix_exp





 !!! Normalize matrix using Frobenius norm !!!
 SUBROUTINE normalize_matrix(Mat)

  COMPLEX(KIND=DP), INTENT(INOUT) :: Mat(:,:)
  COMPLEX(KIND=DP)                :: F_norm

  F_norm = matrix_frobenius_norm(Mat)

  IF(ABS(F_norm) .GT. 0.0d0) Mat = Mat/F_norm

 END SUBROUTINE normalize_matrix





 !!! Calculate Frobenius norm of a matrix !!!
 FUNCTION matrix_frobenius_norm(Mat) result(F_norm)

  COMPLEX(KIND=DP), INTENT(IN) :: Mat(:,:)
  COMPLEX(KIND=DP)             :: F_norm

  F_norm = SQRT(TRACE(TENSMUL(Mat,  TensTRANSPOSE(Mat, 'HC'))))

 END FUNCTION matrix_frobenius_norm





 !!! Calculate L2 norm of matrix !!!
 FUNCTION matrix_L2_norm(Mat) result(L2_norm)

  COMPLEX(KIND=DP), INTENT(IN) :: Mat(:,:)
  COMPLEX(KIND=DP)             :: L2_norm

  COMPLEX(KIND=DP) :: eval

  CALL dominant_eig(eval=eval, tensA=TENSMUL(Mat,  TensTRANSPOSE(Mat, 'HC')))

  L2_norm = SQRT(eval)

 END FUNCTION matrix_L2_norm



 !! Calculate simple L2 norm
 FUNCTION L2_norm(vec) result(norm)

  COMPLEX(KIND=DP), INTENT(IN) :: vec(:)
  COMPLEX(KIND=DP)             :: norm

  norm = SQRT(SUM(CONJG(vec) * vec))           

 END FUNCTION L2_norm


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




 !!!!!!!!!!!!!!! Represent diagonal matrix as matrix or vector !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Represent diagonal matrix as a vector !!!
 FUNCTION DiagonalAsVec(diagMat) result(diagVec)

  COMPLEX(KIND=DP), INTENT(IN)  :: diagMat(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: diagVec(:)

  INTEGER :: dimA(2), i

  !Determine shape
  dimA = shape(diagMat)

  !Error if diagMat is not square
  CALL check_sizes_equal(dimA(1), dimA(2), "repDiagAsVec: diagMat must be square matrix ")

  !Allocate diagVec
  ALLOCATE(diagVec(dimA(1))); diagVec(:) = (0.0d0,0.0d0)

  !Write to diagVec
  DO i=1,dimA(1)
     diagVec(i) = diagMat(i,i)
  end DO

 END FUNCTION DiagonalAsVec





 !!! Represent diagonal matrix as a matrix !!!
 FUNCTION DiagonalAsMat(diagVec) result(diagMat)

  COMPLEX(KIND=DP), INTENT(IN)  :: diagVec(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: diagMat(:,:)

  INTEGER :: dimA, i

  !Determine shape
  dimA = SIZE(diagVec)

  !Allocate diagVec
  ALLOCATE(diagMat(dimA, dimA)); diagMat(:,:) = (0.0d0,0.0d0)

  !Write to diagVec
  DO i=1,dimA
     diagMat(i,i) = diagVec(i)
  end DO

 END FUNCTION DiagonalAsMat

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE matrix_functions
