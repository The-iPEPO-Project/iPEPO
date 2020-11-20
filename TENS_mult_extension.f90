MODULE TENS_mult_extension

 USE basic_functions
 USE error_handling
 USE array_utility
 USE TENS_helper

 implicit none

 interface LambdaMUL
    module procedure mult_2D_lambda
    module procedure mult_3D_lambda
    module procedure mult_4D_lambda
    module procedure mult_5D_lambda
    module procedure mult_6D_lambda
    module procedure mult_2D_lambda_2D
    module procedure mult_3D_lambda_3D
    module procedure mult_4D_lambda_4D
    module procedure mult_5D_lambda_5D
    module procedure mult_6D_lambda_6D
 end interface LambdaMUL

 interface LambdaDIV
    module procedure divide_2D_lambda
    module procedure divide_3D_lambda
    module procedure divide_4D_lambda
    module procedure divide_5D_lambda
    module procedure divide_6D_lambda
 end interface LambdaDIV

 interface TensKRON
    module procedure kronecker_4D_4D
    module procedure kronecker_3D_3D
    module procedure kronecker_2D_2D_result_2D
    module procedure kronecker_1D_1D
 end interface TensKRON

 interface TRACE
    module procedure trace_2D
    module procedure trace_3D
    module procedure trace_4D
    module procedure trace_6D
 end interface TRACE

 private !! hides all items not listed in public statement 
 public LambdaMUL, LambdaDIV, TensKRON, TRACE

CONTAINS


 !!!!!!!!!!!!!!!!!!!!!!!!!!!! LambdaMUL -- multiply tensA by lambda tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Contract 2D tensor and lambda tensor (a diagonal matrix represented as vector) !!!
 SUBROUTINE mult_2D_lambda(tensA, lambda, MULT) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)
  CHARACTER(LEN=*), INTENT(IN)    :: MULT

  REAL(KIND=DP), PARAMETER :: eps = 1.0d-08 
  COMPLEX(KIND=DP)         :: fac

  !! Dimensions & indices
  INTEGER :: dimA(2), dimL(1)
  INTEGER :: i
  
  !! Get dims
  dimA = shape(tensA); dimL = shape(lambda)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('1')

       CALL check_sizes_equal(dimL(1), dimA(1), "mult_2D_lambda -- must have dimL(1) = dimA(1) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(i,:) = fac * tensA(i,:)
       end DO
        
  CASE('2')

       CALL check_sizes_equal(dimL(1), dimA(2), "mult_2D_lambda -- must have dimL(1) = dimA(2) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,i) = tensA(:,i) * fac
       end DO

  CASE DEFAULT
       CALL invalid_flag("mult_2D_lambda -- invalid MULT", MULT)
  end SELECT
  
 END SUBROUTINE mult_2D_lambda




 !!! Contract 3D tensor and lambda tensor (a diagonal matrix represented as vector) !!!
 SUBROUTINE mult_3D_lambda(tensA, lambda, MULT) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)
  CHARACTER(LEN=*), INTENT(IN)    :: MULT

  REAL(KIND=DP), PARAMETER :: eps = 1.0d-08 
  COMPLEX(KIND=DP)         :: fac

  !! Dimensions & indices
  INTEGER :: dimA(3), dimL(1)
  INTEGER :: i
  
  !! Get dims
  dimA = shape(tensA); dimL = shape(lambda)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('2')

       CALL check_sizes_equal(dimL(1), dimA(2), "mult_3D_lambda -- must have dimL(1) = dimA(2) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,i,:) = fac * tensA(:,i,:)
       end DO
        
  CASE('3')

       CALL check_sizes_equal(dimL(1), dimA(3), "mult_3D_lambda -- must have dimL(1) = dimA(3) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,:,i) = tensA(:,:,i) * fac
       end DO

  CASE DEFAULT
       CALL invalid_flag("mult_3D_lambda -- invalid MULT", MULT)
  end SELECT
  
 END SUBROUTINE mult_3D_lambda





 !!! Contract 3D tensor and lambda tensor (a diagonal matrix represented as vector) !!!
 SUBROUTINE mult_4D_lambda(tensA, lambda, MULT) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)
  CHARACTER(LEN=*), INTENT(IN)    :: MULT

  REAL(KIND=DP), PARAMETER :: eps = 1.0d-08 
  COMPLEX(KIND=DP)         :: fac

  !! Dimensions & indices
  INTEGER :: dimA(4), dimL(1)
  INTEGER :: i
  
  !! Get dims
  dimA = shape(tensA); dimL = shape(lambda)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('3')

       CALL check_sizes_equal(dimL(1), dimA(3), "mult_4D_lambda -- must have dimL(1) = dimA(3) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,:,i,:) = fac * tensA(:,:,i,:)
       end DO
        
  CASE('4')

       CALL check_sizes_equal(dimL(1), dimA(4), "mult_4D_lambda -- must have dimL(1) = dimA(4) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,:,:,i) = tensA(:,:,:,i) * fac
       end DO

  CASE DEFAULT
       CALL invalid_flag("mult_4D_lambda -- invalid MULT", MULT)
  end SELECT
  
 END SUBROUTINE mult_4D_lambda






 !!! Contract 5D tensor and lambda tensor (a diagonal matrix represented as vector) !!!
 SUBROUTINE mult_5D_lambda(tensA, lambda, MULT) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)
  CHARACTER(LEN=*), INTENT(IN)    :: MULT

  REAL(KIND=DP), PARAMETER :: eps = 1.0d-08 
  COMPLEX(KIND=DP)         :: fac

  !! Dimensions & indices
  INTEGER :: dimA(5), dimL(1)
  INTEGER :: i
  
  !! Get dims
  dimA = shape(tensA); dimL = shape(lambda)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('2')

       CALL check_sizes_equal(dimL(1), dimA(2), "mult_5D_lambda -- must have dimL(1) = dimA(2) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,i,:,:,:) = fac * tensA(:,i,:,:,:)
       end DO
        
  CASE('3')

       CALL check_sizes_equal(dimL(1), dimA(3), "mult_5D_lambda -- must have dimL(1) = dimA(3) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,:,i,:,:) = tensA(:,:,i,:,:) * fac
       end DO

  CASE('4')

       CALL check_sizes_equal(dimL(1), dimA(4), "mult_5D_lambda -- must have dimL(1) = dimA(4) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,:,:,i,:) = fac * tensA(:,:,:,i,:)
       end DO
        
  CASE('5')

       CALL check_sizes_equal(dimL(1), dimA(5), "mult_5D_lambda -- must have dimL(1) = dimA(5) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,:,:,:,i) = tensA(:,:,:,:,i) * fac
       end DO

  CASE DEFAULT
       CALL invalid_flag("mult_5D_lambda -- invalid MULT", MULT)
  end SELECT
  
 END SUBROUTINE mult_5D_lambda




 !!! Contract 6D tensor and lambda tensor (a diagonal matrix represented as vector) !!!
 SUBROUTINE mult_6D_lambda(tensA, lambda, MULT) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:,:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)
  CHARACTER(LEN=*), INTENT(IN)    :: MULT

  REAL(KIND=DP), PARAMETER :: eps = 1.0d-08 
  COMPLEX(KIND=DP)         :: fac

  !! Dimensions & indices
  INTEGER :: dimA(6), dimL(1)
  INTEGER :: i
  
  !! Get dims
  dimA = shape(tensA); dimL = shape(lambda)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('3')

       CALL check_sizes_equal(dimL(1), dimA(3), "mult_6D_lambda -- must have dimL(1) = dimA(3) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,:,i,:,:,:) = fac * tensA(:,:,i,:,:,:)
       end DO
        
  CASE('4')

       CALL check_sizes_equal(dimL(1), dimA(4), "mult_6D_lambda -- must have dimL(1) = dimA(4) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,:,:,i,:,:) = tensA(:,:,:,i,:,:) * fac
       end DO

  CASE('5')

       CALL check_sizes_equal(dimL(1), dimA(5), "mult_6D_lambda -- must have dimL(1) = dimA(5) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,:,:,:,i,:) = fac * tensA(:,:,:,:,i,:)
       end DO
        
  CASE('6')

       CALL check_sizes_equal(dimL(1), dimA(6), "mult_6D_lambda -- must have dimL(1) = dimA(6) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .LT. eps) fac = 0.0d0
          tensA(:,:,:,:,:,i) = tensA(:,:,:,:,:,i) * fac
       end DO

  CASE DEFAULT
       CALL invalid_flag("mult_6D_lambda -- invalid MULT", MULT)
  end SELECT
  
 END SUBROUTINE mult_6D_lambda




 !!! Absorb lambda tensor into tensA, tensB !!!
 SUBROUTINE mult_2D_lambda_2D(tensA, lambda, tensB)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:), tensB(:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)

  !! Calculate [tensA * sqrt(lambda)] * [sqrt(lambda) * tensB]
  CALL mult_2D_lambda(tensA, SQRT(lambda), MULT='2') 
  CALL mult_2D_lambda(tensB, SQRT(lambda), MULT='1') 

 END SUBROUTINE mult_2D_lambda_2D




 !!! Absorb lambda tensor into tensA, tensB !!!
 SUBROUTINE mult_3D_lambda_3D(tensA, lambda, tensB)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:), tensB(:,:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)

  !! Calculate [tensA * sqrt(lambda)] * [sqrt(lambda) * tensB]
  CALL mult_3D_lambda(tensA, SQRT(lambda), MULT='3') 
  CALL mult_3D_lambda(tensB, SQRT(lambda), MULT='2') 

 END SUBROUTINE mult_3D_lambda_3D



 !!! Absorb lambda tensor into tensA, tensB !!!
 SUBROUTINE mult_4D_lambda_4D(tensA, lambda, tensB)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:), tensB(:,:,:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)

  !! Calculate [tensA * sqrt(lambda)] * [sqrt(lambda) * tensB]
  CALL mult_4D_lambda(tensA, SQRT(lambda), MULT='4') 
  CALL mult_4D_lambda(tensB, SQRT(lambda), MULT='3') 

 END SUBROUTINE mult_4D_lambda_4D



 !!! Absorb lambda tensor into tensA, tensB !!!
 SUBROUTINE mult_5D_lambda_5D(tensA, lambda, tensB, doubleMULT)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:,:), tensB(:,:,:,:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)
  CHARACTER(LEN=*), INTENT(IN)    :: doubleMULT

  SELECT CASE(doubleMULT)
  CASE('23')

      !! Calculate [tensA * sqrt(lambda)] * [sqrt(lambda) * tensB]
      CALL mult_5D_lambda(tensA, SQRT(lambda), MULT='3') 
      CALL mult_5D_lambda(tensB, SQRT(lambda), MULT='2') 

  CASE('45')

      !! Calculate [tensA * sqrt(lambda)] * [sqrt(lambda) * tensB]
      CALL mult_5D_lambda(tensA, SQRT(lambda), MULT='5') 
      CALL mult_5D_lambda(tensB, SQRT(lambda), MULT='4') 
 
  CASE DEFAULT
      CALL invalid_flag("mult_5D_lambda_5D -- invalid doubleMULT", doubleMULT)
  end SELECT

 END SUBROUTINE mult_5D_lambda_5D



 !!! Absorb lambda tensor into tensA, tensB !!!
 SUBROUTINE mult_6D_lambda_6D(tensA, lambda, tensB, doubleMULT)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:,:,:), tensB(:,:,:,:,:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)
  CHARACTER(LEN=*), INTENT(IN)    :: doubleMULT

  SELECT CASE(doubleMULT)
  CASE('34')

      !! Calculate [tensA * sqrt(lambda)] * [sqrt(lambda) * tensB]
      CALL mult_6D_lambda(tensA, SQRT(lambda), MULT='4') 
      CALL mult_6D_lambda(tensB, SQRT(lambda), MULT='3') 

  CASE('56')

      !! Calculate [tensA * sqrt(lambda)] * [sqrt(lambda) * tensB]
      CALL mult_6D_lambda(tensA, SQRT(lambda), MULT='6') 
      CALL mult_6D_lambda(tensB, SQRT(lambda), MULT='5') 
 
  CASE DEFAULT
      CALL invalid_flag("mult_6D_lambda_6D -- invalid doubleMULT", doubleMULT)
  end SELECT

 END SUBROUTINE mult_6D_lambda_6D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LambdaDIV -- divide tensA by lambda tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Divide 2D tensor by lambda tensor (a diagonal matrix represented as vector) !!!
 SUBROUTINE divide_2D_lambda(tensA, lambda, MULT) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)
  CHARACTER(LEN=*), INTENT(IN)    :: MULT

  REAL(KIND=DP), PARAMETER :: eps = 1.0d-08 
  COMPLEX(KIND=DP)         :: fac

  !! Dimensions & indices
  INTEGER :: dimA(2), dimL(1)
  INTEGER :: i
  
  !! Get dims
  dimA = shape(tensA); dimL = shape(lambda)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('1')

       CALL check_sizes_equal(dimL(1), dimA(1), "divide_2D_lambda -- must have dimL(1) = dimA(1) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(i,:) = tensA(i,:) / fac
          end IF
       end DO
        
  CASE('2')

       CALL check_sizes_equal(dimL(1), dimA(2), "divide_2D_lambda -- must have dimL(1) = dimA(2) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,i) = tensA(:,i) / fac
          end IF
       end DO

  CASE DEFAULT
       CALL invalid_flag("divide_2D_lambda -- invalid MULT ", MULT)
  end SELECT
  
 END SUBROUTINE divide_2D_lambda




 !!! Divide 3D tensor by lambda tensor (a diagonal matrix represented as vector) !!!
 SUBROUTINE divide_3D_lambda(tensA, lambda, MULT) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)
  CHARACTER(LEN=*), INTENT(IN)    :: MULT

  REAL(KIND=DP), PARAMETER :: eps = 1.0d-08 
  COMPLEX(KIND=DP)         :: fac

  !! Dimensions & indices
  INTEGER :: dimA(3), dimL(1)
  INTEGER :: i
  
  !! Get dims
  dimA = shape(tensA); dimL = shape(lambda)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('2')

       CALL check_sizes_equal(dimL(1), dimA(2), "divide_3D_lambda -- must have dimL(1) = dimA(2) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,i,:) = tensA(:,i,:) / fac
          end IF
       end DO
        
  CASE('3')

       CALL check_sizes_equal(dimL(1), dimA(3), "divide_3D_lambda -- must have dimL(1) = dimA(3) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,:,i) = tensA(:,:,i) / fac
          end IF
       end DO

  CASE DEFAULT
       CALL invalid_flag("divide_3D_lambda -- invalid MULT ", MULT)
  end SELECT
  
 END SUBROUTINE divide_3D_lambda





 !!! Divide 4D tensor by lambda tensor (a diagonal matrix represented as vector) !!!
 SUBROUTINE divide_4D_lambda(tensA, lambda, MULT) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)
  CHARACTER(LEN=*), INTENT(IN)    :: MULT

  REAL(KIND=DP), PARAMETER :: eps = 1.0d-08 
  COMPLEX(KIND=DP)         :: fac

  !! Dimensions & indices
  INTEGER :: dimA(4), dimL(1)
  INTEGER :: i
  
  !! Get dims
  dimA = shape(tensA); dimL = shape(lambda)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('3')

       CALL check_sizes_equal(dimL(1), dimA(3), "divide_4D_lambda -- must have dimL(1) = dimA(3) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,:,i,:) = tensA(:,:,i,:) / fac
          end IF
       end DO
        
  CASE('4')

       CALL check_sizes_equal(dimL(1), dimA(4), "divide_4D_lambda -- must have dimL(1) = dimA(4) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,:,:,i) = tensA(:,:,:,i) / fac
          end IF
       end DO

  CASE DEFAULT
       CALL invalid_flag("divide_4D_lambda -- invalid MULT ", MULT)
  end SELECT
  
 END SUBROUTINE divide_4D_lambda





 !!! Divide 5D tensor by lambda tensor (a diagonal matrix represented as vector) !!!
 SUBROUTINE divide_5D_lambda(tensA, lambda, MULT) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)
  CHARACTER(LEN=*), INTENT(IN)    :: MULT

  REAL(KIND=DP), PARAMETER :: eps = 1.0d-08 
  COMPLEX(KIND=DP)         :: fac

  !! Dimensions & indices
  INTEGER :: dimA(5), dimL(1)
  INTEGER :: i
  
  !! Get dims
  dimA = shape(tensA); dimL = shape(lambda)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('2')

       CALL check_sizes_equal(dimL(1), dimA(2), "divide_5D_lambda -- must have dimL(1) = dimA(2) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,i,:,:,:) = tensA(:,i,:,:,:) / fac
          end IF
       end DO
        
  CASE('3')

       CALL check_sizes_equal(dimL(1), dimA(3), "divide_5D_lambda -- must have dimL(1) = dimA(3) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,:,i,:,:) = tensA(:,:,i,:,:) / fac
          end IF
       end DO

  CASE('4')

       CALL check_sizes_equal(dimL(1), dimA(4), "divide_5D_lambda -- must have dimL(1) = dimA(4) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,:,:,i,:) = tensA(:,:,:,i,:) / fac
          end IF
       end DO
        
  CASE('5')

       CALL check_sizes_equal(dimL(1), dimA(5), "divide_5D_lambda -- must have dimL(1) = dimA(5) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,:,:,:,i) = tensA(:,:,:,:,i) / fac
          end IF
       end DO

  CASE DEFAULT
       CALL invalid_flag("divide_5D_lambda -- invalid MULT ", MULT)
  end SELECT
  
 END SUBROUTINE divide_5D_lambda





 !!! Divide 6D tensor by lambda tensor (a diagonal matrix represented as vector) !!!
 SUBROUTINE divide_6D_lambda(tensA, lambda, MULT) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:,:,:)
  COMPLEX(KIND=DP), INTENT(IN)    :: lambda(:)
  CHARACTER(LEN=*), INTENT(IN)    :: MULT

  REAL(KIND=DP), PARAMETER :: eps = 1.0d-08 
  COMPLEX(KIND=DP)         :: fac

  !! Dimensions & indices
  INTEGER :: dimA(6), dimL(1)
  INTEGER :: i
  
  !! Get dims
  dimA = shape(tensA); dimL = shape(lambda)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('3')

       CALL check_sizes_equal(dimL(1), dimA(3), "divide_6D_lambda -- must have dimL(1) = dimA(3) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,:,i,:,:,:) = tensA(:,:,i,:,:,:) / fac
          end IF
       end DO
        
  CASE('4')

       CALL check_sizes_equal(dimL(1), dimA(4), "divide_6D_lambda -- must have dimL(1) = dimA(4) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,:,:,i,:,:) = tensA(:,:,:,i,:,:) / fac
          end IF
       end DO

  CASE('5')

       CALL check_sizes_equal(dimL(1), dimA(5), "divide_6D_lambda -- must have dimL(1) = dimA(5) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,:,:,:,i,:) = tensA(:,:,:,:,i,:) / fac
          end IF
       end DO
        
  CASE('6')

       CALL check_sizes_equal(dimL(1), dimA(6), "divide_6D_lambda -- must have dimL(1) = dimA(6) ")

       DO i=1,dimL(1)
          fac = lambda(i)
          IF(abs(fac) .GT. eps) THEN
             tensA(:,:,:,:,:,i) = tensA(:,:,:,:,:,i) / fac
          end IF
       end DO

  CASE DEFAULT
       CALL invalid_flag("divide_6D_lambda -- invalid MULT ", MULT)
  end SELECT
  
 END SUBROUTINE divide_6D_lambda

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!! TensKRON -- calculate tensor Kronecker product !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Kronecker product 4D x 4D !!!
 FUNCTION kronecker_4D_4D(tensA, tensB) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:,:,:), tensB(:,:,:,:)
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:,:,:)

  INTEGER :: a1, a2, a3, a4
  INTEGER :: b1, b2, b3, b4
  INTEGER :: dimA(4), dimB(4)

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Outer product of tensors
  ALLOCATE(tensO(dimA(1)*dimB(1), dimA(2)*dimB(2), dimA(3)*dimB(3), dimA(4)*dimB(4)))

  DO a1=1,dimA(1) 
    DO b1=1,dimB(1) 

      DO a2=1,dimA(2)  
        DO b2=1,dimB(2)

          DO a3=1,dimA(3)
            DO b3=1,dimB(3)

              DO a4=1,dimA(4)
                DO b4=1,dimB(4)
                   
                   tensO(ICOM(a1,b1,dimB(1)), ICOM(a2,b2,dimB(2)), ICOM(a3,b3,dimB(3)), ICOM(a4,b4,dimB(4))) = &

                    & tensA(a1, a2, a3, a4) * tensB(b1, b2, b3, b4)
                end DO
              end DO
            end DO
          end DO
        end DO
      end DO
    end DO
  end DO

 END FUNCTION kronecker_4D_4D




 !!! Kronecker product 3D x 3D !!!
 FUNCTION kronecker_3D_3D(tensA, tensB) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:), tensB(:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:)

  INTEGER :: a1, a2, a3
  INTEGER :: b1, b2, b3
  INTEGER :: dimA(3), dimB(3)

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Outer product of tensors
  ALLOCATE(tensO(dimA(1), dimB(1), dimA(2)*dimB(2), dimA(3)*dimB(3)))

  DO a1=1,dimA(1) 
    DO b1=1,dimB(1) 

      DO a2=1,dimA(2)  
        DO b2=1,dimB(2)
          DO a3=1,dimA(3)
            DO b3=1,dimB(3)
               tensO(a1, b1, ICOM(a2,b2,dimB(2)), ICOM(a3,b3,dimB(3))) = tensA(a1, a2, a3) * tensB(b1, b2, b3)
             end DO
          end DO
        end DO
      end DO

    END DO
  END DO

 END FUNCTION kronecker_3D_3D



 !!! Kronecker product 2D x 2D --> 2D !!!
 FUNCTION kronecker_2D_2D_result_2D(tensA, tensB) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:), tensB(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:)

  INTEGER :: a1, a2, b1, b2
  INTEGER :: dimA(2), dimB(2)

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Outer product of tensors
  ALLOCATE(tensO(dimA(1)*dimB(1), dimA(2)*dimB(2)))

  DO a1=1,dimA(1) 
    DO b1=1,dimB(1) 
      DO a2=1,dimA(2)  
        DO b2=1,dimB(2)
           tensO(ICOM(a1,b1,dimB(1)), ICOM(a2,b2,dimB(2))) = tensA(a1, a2) * tensB(b1, b2)
        end DO
      end DO
    end DO
  end DO

 END FUNCTION kronecker_2D_2D_result_2D



 !!! Kronecker product 2D x 2D --> 4D !!!
 FUNCTION kronecker_2D_2D_result_4D(tensA, tensB) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:), tensB(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:)

  INTEGER :: a1, a2, b1, b2
  INTEGER :: dimA(2), dimB(2)

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Outer product of tensors
  ALLOCATE(tensO(dimA(1), dimB(1), dimA(2), dimB(2)))

  DO a1=1,dimA(1) 
    DO b1=1,dimB(1) 
      DO a2=1,dimA(2)  
        DO b2=1,dimB(2)
           tensO(a1, b1, a2, b2) = tensA(a1, a2) * tensB(b1, b2)
        end DO
      end DO
    end DO
  end DO

 END FUNCTION kronecker_2D_2D_result_4D




 !!! Kronecker product 1D x 1D --> 1D !!!
 FUNCTION kronecker_1D_1D(tensA, tensB) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:), tensB(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:)

  INTEGER :: a, b, dimA, dimB

  !! Get dims of tensA, tensB
  dimA = SIZE(tensA); dimB = SIZE(tensB)

  !! Outer product of tensors
  ALLOCATE(tensO(dimA*dimB))

  DO a=1,dimA 
    DO b=1,dimB 
       tensO(ICOM(a,b,dimB)) = tensA(a) * tensB(b)
    end DO
  end DO

 END FUNCTION kronecker_1D_1D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TRACE: trace out given legs of a tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Trace legs of 2D tensor !!!
 FUNCTION trace_2D(tensA) result(trace)

  COMPLEX(KIND=DP), INTENT(IN) :: tensA(:,:)
  COMPLEX(KIND=DP)             :: trace

  INTEGER :: i, dimA(2)

  dimA = shape(tensA)

  CALL check_sizes_equal(dimA(1), dimA(2), "trace_2D: traced dims must be equal")

  !! Get trace
  trace = 0.0d0
  DO i=1,dimA(1)
     trace = trace + tensA(i,i)
  end DO

 END FUNCTION trace_2D





 !!! Trace a pair of legs of 3D tensor !!!
 FUNCTION trace_3D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: legs
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:)

  INTEGER :: i, dimA(3)

  dimA = shape(tensA)

  SELECT CASE(legs)
  CASE('12')

       !! Get trace over axes-1,2
       CALL allocateTens(tensO, (/dimA(3)/))
       DO i=1,dimA(3)
          tensO(i) = TRACE(tensA(:,:,i))
       end DO

  CASE('23')

       !! Get trace over axes-2,3
       CALL allocateTens(tensO, (/dimA(1)/))
       DO i=1,dimA(1)
          tensO(i) = TRACE(tensA(i,:,:))
       end DO

  CASE('13')

       !! Get trace over axes-1,3
       CALL allocateTens(tensO, (/dimA(2)/))
       DO i=1,dimA(2)
          tensO(i) = TRACE(tensA(:,i,:))
       end DO

  CASE DEFAULT
       CALL invalid_flag("trace_3D -- invalid pair of legs to trace over ", legs)
  end SELECT 

 END FUNCTION trace_3D






 !!! Trace a pair of legs of 4D tensor !!!
 FUNCTION trace_4D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: legs
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:)

  INTEGER :: i, j, dimA(4)

  dimA = shape(tensA)

  SELECT CASE(legs)
  CASE('12')

       !! Get trace over axes-1,2
       CALL allocateTens(tensO,  (/dimA(3), dimA(4)/))
       DO i=1,dimA(3)
         DO j=1,dimA(4)
            tensO(i,j) = TRACE(tensA(:,:,i,j))
         end DO
       end DO

  CASE('34')

       !! Get trace over axes-3,4
       CALL allocateTens(tensO,  (/dimA(1), dimA(2)/))
       DO i=1,dimA(1)
         DO j=1,dimA(2)
            tensO(i,j) = TRACE(tensA(i,j,:,:))
         end DO
       end DO

  CASE('13')

       !! Get trace over axes-1,3
       CALL allocateTens(tensO,  (/dimA(2), dimA(4)/))
       DO i=1,dimA(2)
         DO j=1,dimA(4)
            tensO(i,j) = TRACE(tensA(:,i,:,j))
         end DO
       end DO

  CASE('24')

       !! Get trace over axes-2,4
       CALL allocateTens(tensO,  (/dimA(1), dimA(3)/))
       DO i=1,dimA(1)
         DO j=1,dimA(3)
            tensO(i,j) = TRACE(tensA(i,:,j,:))
         end DO
       end DO

  CASE DEFAULT
       CALL invalid_flag("trace_4D -- invalid pair of legs to trace over ", legs)
  end SELECT 

 END FUNCTION trace_4D







 !!! Trace a pair of legs of 6D tensor !!!
 FUNCTION trace_6D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:,:,:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: legs
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:)

  INTEGER :: i, j, dimA(6)

  dimA = shape(tensA)

  SELECT CASE(legs)
  CASE('12')

       !! Get trace over axes-1,2
       CALL allocateTens(tensO,  (/dimA(3), dimA(4), dimA(5), dimA(6)/))
       DO i=1,dimA(3)
         DO j=1,dimA(4)
            tensO(i,j,:,:) = TRACE(tensA(:,:,i,j,:,:), '12')
         end DO
       end DO

  CASE('34')

       !! Get trace over axes-3,4
       CALL allocateTens(tensO,  (/dimA(1), dimA(2), dimA(5), dimA(6)/))
       DO i=1,dimA(1)
         DO j=1,dimA(2)
            tensO(i,j,:,:) = TRACE(tensA(i,j,:,:,:,:), '12')
         end DO
       end DO

  CASE('56')

       !! Get trace over axes-5,6
       CALL allocateTens(tensO,  (/dimA(1), dimA(2), dimA(3), dimA(4)/))
       DO i=1,dimA(1)
         DO j=1,dimA(2)
            tensO(i,j,:,:) = TRACE(tensA(i,j,:,:,:,:), '34')
         end DO
       end DO

  CASE DEFAULT
       CALL invalid_flag("trace_6D -- invalid pair of legs to trace over ", legs)
  end SELECT 

 END FUNCTION trace_6D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE TENS_mult_extension
