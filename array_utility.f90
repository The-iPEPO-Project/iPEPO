MODULE array_utility

 USE basic_functions
 USE error_handling

 !!!!!!!!!!!!!!!!!!!!!!! CONTENTS: !!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! (1) ALLOCATE TENSOR (allocateTens)
 !!
 !! (2) COPY TENSOR (copyTens) 
 !!
 !! (3) GENERATE SOME ELEMENTARY ARRAYS
 !!
 !! (4) PRINTING ROUTINES
 !!
 !! (5) RESCALING ROUTINES
 !!
 !! (6) CONCATENATE/SPLIT/INVERT ARRAYS
 !!
 !! (7) PURIFICATION ROUTINES 
 !!
 !! (8) ROUNDING/FILTERING ROUTINES
 !!
 !! (9) ADD NOISE to tensor/value 
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IMPLICIT NONE

 interface allocateTens
    module procedure allocate_tens_6D
    module procedure allocate_tens_5D
    module procedure allocate_tens_4D
    module procedure allocate_tens_3D
    module procedure allocate_tens_2D
    module procedure allocate_tens_1D
 end interface allocateTens

 interface copyTens
    module procedure copy_tens_5D
    module procedure copy_tens_4D
    module procedure copy_tens_3D
    module procedure copy_tens_2D
    module procedure copy_tens_1D
 end interface copyTens

 interface array
    module procedure array_5D
    module procedure array_4D
    module procedure array_3D
    module procedure array_2D
    module procedure array_1D
 end interface array

 interface rescaleTensMaxel
    module procedure rescale_tens_by_max_element_3D
    module procedure rescale_tens_by_max_element_2D
 end interface rescaleTensMaxel

 interface maxElement
    module procedure max_element_5D
    module procedure max_element_4D
    module procedure max_element_3D
    module procedure max_element_2D
    module procedure max_element_1D
 end interface maxElement

 interface roundTens
    module procedure roundPEPS
    module procedure roundMpo
    module procedure roundMps
    module procedure roundMatrix
    module procedure roundArray
 end interface roundTens

 interface filterTens
    module procedure filterPeps
    module procedure filterMpo
    module procedure filterMps
    module procedure filterMatrix
    module procedure filterArray
 end interface filterTens

 interface add_noise
    module procedure add_noise_5D
    module procedure add_noise_4D
    module procedure add_noise_3D
    module procedure add_noise_2D
    module procedure add_noise_1D
    module procedure add_noise_Val
 end interface add_noise

 private
 public allocateTens, copyTens                  
 public array, matEye, vecOnes   
 public add_noise               
 public roundTens, roundValue, filterTens
 public rescaleTensMaxel, maxElement, concatenate_vecs, decatenate_vecs, invert_array
 public purify_hermitian_tens, purify_symmetric_tens, purify_positive_tens, removeZeroEvals
 public printMatrix, printVector

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (1) ALLOCATE TENSOR (allocateTens) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !Allocate 1D tensor
 SUBROUTINE allocate_tens_1D(tensA, dimA)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensA(:)
  INTEGER, INTENT(IN) :: dimA(1)

  !deallocate if previously allocated
  IF(ALLOCATED(tensA)) DEALLOCATE(tensA)

  !allocate tens
  ALLOCATE(tensA(dimA(1))); tensA(:) = (0.0d0,0.0d0)

 END SUBROUTINE allocate_tens_1D


 !Allocate 2D tensor
 SUBROUTINE allocate_tens_2D(tensA, dimA)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensA(:,:)
  INTEGER, INTENT(IN) :: dimA(2)

  !deallocate if previously allocated
  IF(ALLOCATED(tensA)) DEALLOCATE(tensA)

  !allocate tens
  ALLOCATE(tensA(dimA(1), dimA(2))); tensA(:,:) = (0.0d0,0.0d0)

 END SUBROUTINE allocate_tens_2D



 !Allocate 3D tensor
 SUBROUTINE allocate_tens_3D(tensA, dimA)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensA(:,:,:)
  INTEGER, INTENT(IN) :: dimA(3)

  !deallocate if previously allocated
  IF(ALLOCATED(tensA)) DEALLOCATE(tensA)

  !allocate tens
  ALLOCATE(tensA(dimA(1), dimA(2), dimA(3))); tensA(:,:,:) = (0.0d0,0.0d0)

 END SUBROUTINE allocate_tens_3D



 !Allocate 4D tensor
 SUBROUTINE allocate_tens_4D(tensA, dimA)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensA(:,:,:,:)
  INTEGER, INTENT(IN) :: dimA(4)

  !deallocate if previously allocated
  IF(ALLOCATED(tensA)) DEALLOCATE(tensA)

  !allocate tens
  ALLOCATE(tensA(dimA(1), dimA(2), dimA(3), dimA(4))); tensA(:,:,:,:) = (0.0d0,0.0d0)

 END SUBROUTINE allocate_tens_4D


 !Allocate 5D tensor
 SUBROUTINE allocate_tens_5D(tensA, dimA)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensA(:,:,:,:,:)
  INTEGER, INTENT(IN) :: dimA(5)

  !deallocate if previously allocated
  IF(ALLOCATED(tensA)) DEALLOCATE(tensA)

  !allocate tens
  ALLOCATE(tensA(dimA(1), dimA(2), dimA(3), dimA(4), dimA(5))); tensA(:,:,:,:,:) = (0.0d0,0.0d0)

 END SUBROUTINE allocate_tens_5D



 !Allocate 6D tensor
 SUBROUTINE allocate_tens_6D(tensA, dimA)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensA(:,:,:,:,:,:)
  INTEGER, INTENT(IN) :: dimA(6)

  !deallocate if previously allocated
  IF(ALLOCATED(tensA)) DEALLOCATE(tensA)

  !allocate tens
  ALLOCATE(tensA(dimA(1), dimA(2), dimA(3), dimA(4), dimA(5), dimA(6))); tensA(:,:,:,:,:,:) = (0.0d0,0.0d0)

 END SUBROUTINE allocate_tens_6D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (2) COPY TENSOR (copyTens) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !Copy 1D tensor
 SUBROUTINE copy_tens_1D(tensO, tensI, dimO)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensI(:)
  INTEGER,           INTENT(IN), OPTIONAL      :: dimO(1) !proposed dim of tensO

  INTEGER :: dimI(1), minDim(1)

  !get dims of tensI
  dimI = shape(tensI)

  !deallocate tensO if previously allocated
  IF(ALLOCATED(tensO)) DEALLOCATE(tensO)

  IF(.NOT. PRESENT(dimO)) THEN

     !allocate && copy tens
     CALL allocateTens(tensO, dimI)
     tensO = tensI
  ELSE

     !allocate tensO 
     CALL allocateTens(tensO, dimO)

     !insert tensI into tensO so that it fits in
     minDim = MIN(dimI, dimO)
     tensO(1:minDim(1)) = tensI(1:minDim(1))
  end IF

 END SUBROUTINE copy_tens_1D



 !Copy 2D tensor
 SUBROUTINE copy_tens_2D(tensO, tensI, dimO)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:)
  COMPLEX(KIND=DP), INTENT(IN) :: tensI(:,:)
  INTEGER, INTENT(IN), OPTIONAL :: dimO(2) !proposed dim of tensO

  INTEGER :: dimI(2), minDim(2)

  !get dims
  dimI = shape(tensI)

  !deallocate if previously allocated
  IF(ALLOCATED(tensO)) DEALLOCATE(tensO)

  IF(.NOT. PRESENT(dimO)) THEN

     !allocate && copy tens
     CALL allocateTens(tensO, dimI)
     tensO = tensI

  ELSE

     !allocate tensO 
     CALL allocateTens(tensO, dimO)

     !insert tensI into tensO so that it fits in
     minDim = MIN(dimI, dimO)
     tensO(1:minDim(1), 1:minDim(2)) = tensI(1:minDim(1), 1:minDim(2))

  end IF

 END SUBROUTINE copy_tens_2D



 !Copy 3D tensor
 SUBROUTINE copy_tens_3D(tensO, tensI, dimO)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:,:)
  COMPLEX(KIND=DP), INTENT(IN) :: tensI(:,:,:)
  INTEGER, INTENT(IN), OPTIONAL :: dimO(3) !proposed dim of tensO

  INTEGER :: dimI(3), minDim(3)

  !get dims
  dimI = shape(tensI)

  !deallocate if previously allocated
  IF(ALLOCATED(tensO)) DEALLOCATE(tensO)

  IF(.NOT. PRESENT(dimO)) THEN

     !allocate && copy tens
     CALL allocateTens(tensO, dimI)
     tensO = tensI
  ELSE

     !allocate tensO 
     CALL allocateTens(tensO, dimO)

     !insert tensI into tensO so that it fits in
     minDim = MIN(dimI, dimO)
     tensO(1:minDim(1), 1:minDim(2), 1:minDim(3)) = tensI(1:minDim(1), 1:minDim(2), 1:minDim(3))
  end IF

 END SUBROUTINE copy_tens_3D



 !Copy 4D tensor
 SUBROUTINE copy_tens_4D(tensO, tensI, dimO)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT)        :: tensO(:,:,:,:)
  COMPLEX(KIND=DP),              INTENT(IN)           :: tensI(:,:,:,:)
  INTEGER,                       INTENT(IN), OPTIONAL :: dimO(4)        !proposed dim of tensO

  INTEGER :: dimI(4), minDim(4)

  !get dims
  dimI = shape(tensI)

  !deallocate if previously allocated
  IF(ALLOCATED(tensO)) DEALLOCATE(tensO)

  IF(.NOT. PRESENT(dimO)) THEN

     !allocate && copy tens
     CALL allocateTens(tensO, dimI)
     tensO = tensI
  ELSE

     !allocate tensO 
     CALL allocateTens(tensO, dimO)

     !insert tensI into tensO so that it fits in
     minDim = MIN(dimI, dimO)
     tensO(1:minDim(1), 1:minDim(2), 1:minDim(3), 1:minDim(4)) = &

           & tensI(1:minDim(1), 1:minDim(2), 1:minDim(3), 1:minDim(4))
  end IF

 END SUBROUTINE copy_tens_4D




 !Copy 5D tensor
 SUBROUTINE copy_tens_5D(tensO, tensI, dimO)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:,:,:,:)
  COMPLEX(KIND=DP), INTENT(IN) :: tensI(:,:,:,:,:)
  INTEGER, INTENT(IN), OPTIONAL :: dimO(5) !proposed dim of tensO

  INTEGER :: dimI(5), minDim(5)

  !get dims
  dimI = shape(tensI)

  !deallocate if previously allocated
  IF(ALLOCATED(tensO)) DEALLOCATE(tensO)

  IF(.NOT. PRESENT(dimO)) THEN

     !allocate && copy tens
     CALL allocateTens(tensO, dimI)
     tensO = tensI
  ELSE

     !allocate tensO 
     CALL allocateTens(tensO, dimO)

     !insert tensI into tensO so that it fits in
     minDim = MIN(dimI, dimO)
     tensO(1:minDim(1), 1:minDim(2), 1:minDim(3), 1:minDim(4), 1:minDim(5)) = & 

                     & tensI(1:minDim(1), 1:minDim(2), 1:minDim(3), 1:minDim(4), 1:minDim(5))
  end IF

 END SUBROUTINE copy_tens_5D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (3) GENERATE SOME ELEMENTARY ARRAYS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Create integer array of dims / integer flags / etc !!!
 FUNCTION array_5D(a1, a2, a3, a4, a5) result(ARR)

  INTEGER, INTENT(IN) :: a1, a2, a3, a4, a5
  INTEGER :: ARR(5)

  ARR(1) = a1; ARR(2) = a2; ARR(3) = a3; ARR(4) = a4; ARR(5) = a5

 END FUNCTION array_5D



 FUNCTION array_4D(a1, a2, a3, a4) result(ARR)

  INTEGER, INTENT(IN) :: a1, a2, a3, a4
  INTEGER :: ARR(4)

  ARR(1) = a1; ARR(2) = a2; ARR(3) = a3; ARR(4) = a4

 END FUNCTION array_4D



 FUNCTION array_3D(a1, a2, a3) result(ARR)

  INTEGER, INTENT(IN) :: a1, a2, a3
  INTEGER :: ARR(3)

  ARR(1) = a1; ARR(2) = a2; ARR(3) = a3

 END FUNCTION array_3D



 FUNCTION array_2D(a1, a2) result(ARR)

  INTEGER, INTENT(IN) :: a1, a2
  INTEGER :: ARR(2)

  ARR(1) = a1; ARR(2) = a2

 END FUNCTION array_2D



 FUNCTION array_1D(a1) result(ARR)

  INTEGER, INTENT(IN) :: a1
  INTEGER :: ARR(1)

  ARR(1) = a1

 END FUNCTION array_1D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !!! Generate identity matrix eye !!!
 FUNCTION matEye(dimA) result(tensA)

  INTEGER, INTENT(IN) :: dimA
  COMPLEX(KIND=DP), ALLOCATABLE :: tensA(:,:)

  INTEGER :: i

  ALLOCATE(tensA(dimA, dimA)); tensA(:,:) = (0.0d0,0.0d0)

  DO i=1,dimA
     tensA(i,i) = (1.0d0,0.0d0)
  end DO

 END FUNCTION matEye



 !!! Generate a vector of ones !!!
 FUNCTION vecOnes(dimA) result(tensA)

  INTEGER, INTENT(IN) :: dimA
  COMPLEX(KIND=DP), ALLOCATABLE :: tensA(:)

  ALLOCATE(tensA(dimA)); tensA(:) = (1.0d0,0.0d0)

 END FUNCTION vecOnes

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (4) PRINTING ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Print vector !!!
 SUBROUTINE printVector(tensA_in, datafile)

  COMPLEX(KIND=DP),  INTENT(IN)           :: tensA_in(:)
  INTEGER,           INTENT(IN), OPTIONAL :: datafile
  
  COMPLEX(KIND=DP), ALLOCATABLE :: tensA(:)
  INTEGER                       :: dimA(1), i

  !! Create a local copy of tensA_in
  CALL copyTens(tensA, tensA_in)

  !! Get dims
  dimA = shape(tensA)

  IF(PRESENT(datafile)) THEN

     !! Print tensA
     DO i=1,dimA(1)
        !IF(ABS(tensA(i)) .GT. 5.0d-04) THEN
           WRITE(datafile,*) "tensA at ", i
           WRITE(datafile,*)  tensA(i)
        !end IF
     end DO
  ELSE

     !! Print tensA
     DO i=1,dimA(1)
        !IF(ABS(tensA(i)) .GT. 5.0d-04) THEN
           WRITE(*,*) "tensA at ", i
           WRITE(*,*)  tensA(i)
        !end IF
     end DO
  end IF

 END SUBROUTINE printVector






 !!! Print matrix !!!
 SUBROUTINE printMatrix(tensA_in, datafile)

  COMPLEX(KIND=DP),  INTENT(IN)           :: tensA_in(:,:)
  INTEGER,           INTENT(IN), OPTIONAL :: datafile
  
  COMPLEX(KIND=DP), ALLOCATABLE :: tensA(:,:)
  INTEGER                       :: dimA(2), i, j

  !! Create a local copy of tensA_in
  CALL copyTens(tensA, tensA_in)

  !! Get dims
  dimA = shape(tensA)

  !! Filter out small values
  !CALL filterMatrix(tensA)

  IF(PRESENT(datafile)) THEN

     !! Print tensA
     DO i=1,dimA(1)
       DO j=1,dimA(2)
          !IF(ABS(tensA(i,j)) .GT. 5.0d-04) THEN
             WRITE(datafile,*) "tensA at ", i, j
             WRITE(datafile,*)  tensA(i,j)
          !end IF
       end DO
     end DO
  ELSE

     !! Print tensA
     DO i=1,dimA(1)
       DO j=1,dimA(2)
          !IF(ABS(tensA(i,j)) .GT. 5.0d-04) THEN
             WRITE(*,*) "tensA at ", i, j
             WRITE(*,*)  tensA(i,j)
          !end IF
       end DO
     end DO
  end IF

 END SUBROUTINE printMatrix

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (5) RESCALING ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !!! Rescale 3D tensor by its max element !!!
 SUBROUTINE rescale_tens_by_max_element_3D(tensA)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:)

  REAL(KIND=DP) :: maxEl

  maxEl = MAXVAL(ABS(tensA))
  tensA = tensA * 1.0d0 * maxEl**(-1.0d0)

 END SUBROUTINE rescale_tens_by_max_element_3D



 !!! Rescale 2D tensor by its max element !!!
 SUBROUTINE rescale_tens_by_max_element_2D(tensA)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:)

  REAL(KIND=DP) :: maxEl

  maxEl = MAXVAL(ABS(tensA))
  tensA = tensA * 1.0d0 * maxEl**(-1.0d0)

 END SUBROUTINE rescale_tens_by_max_element_2D


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! max_element -- FIND MAX ELEMENT IN A COMPLEX ARRAY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !Finds maximum ABS element in a complex array tensA -- 5D version
 FUNCTION max_element_5D(tensA) result(max_element)

  COMPLEX(KIND=DP), INTENT(IN) :: tensA(:,:,:,:,:) !input array
  REAL(KIND=DP) :: max_element !result

  !Temp array tensO = ABS(tensA)
  REAL(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:,:)

  !Indices && dims
  INTEGER :: dimA(5), i, j, k, m, n

  !Get dims
  dimA = shape(tensA)
  
  !Transfer tensA to real array tensO = ABS(tensA)
  ALLOCATE(tensO(dimA(1), dimA(2), dimA(3), dimA(4), dimA(5)))
  Do i=1,dimA(1)
    Do j=1,dimA(2)
      Do k=1,dimA(3)
        Do m=1,dimA(4)
          Do n=1,dimA(5)
             tensO(i,j,k,m,n) = ABS(tensA(i,j,k,m,n))
          end Do
        end Do
      end Do
    end Do
  end Do

  !Find max element in tensO
  max_element = MAXVAL(tensO)

 END FUNCTION max_element_5D




 !Finds maximum ABS element in a complex array tensA -- 4D version
 FUNCTION max_element_4D(tensA) result(max_element)

  COMPLEX(KIND=DP), INTENT(IN) :: tensA(:,:,:,:) !input array
  REAL(KIND=DP) :: max_element !result

  !Temp array tensO = ABS(tensA)
  REAL(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:)

  !Indices && dims
  INTEGER :: dimA(4), i, j, k, l

  !Get dims
  dimA = shape(tensA)
  
  !Transfer tensA to real array tensO = ABS(tensA)
  ALLOCATE(tensO(dimA(1), dimA(2), dimA(3), dimA(4)))
  Do i=1,dimA(1)
    Do j=1,dimA(2)
      Do k=1,dimA(3)
        Do l=1,dimA(4)
           tensO(i,j,k,l) = ABS(tensA(i,j,k,l))
        end Do
      end Do
    end Do
  end Do

  !Find max element in tensO
  max_element = MAXVAL(tensO)

 END FUNCTION max_element_4D



 !Finds maximum ABS element in a complex array tensA -- 3D version
 FUNCTION max_element_3D(tensA) result(max_element)

  COMPLEX(KIND=DP), INTENT(IN) :: tensA(:,:,:) !input array
  REAL(KIND=DP) :: max_element !result

  !Temp array tensO = ABS(tensA)
  REAL(KIND=DP), ALLOCATABLE :: tensO(:,:,:)

  !Indices && dims
  INTEGER :: dimA(3), i, j, k

  !Get dims
  dimA = shape(tensA)
  
  !Transfer tensA to real array tensO = ABS(tensA)
  ALLOCATE(tensO(dimA(1), dimA(2), dimA(3)))
  Do i=1,dimA(1)
    Do j=1,dimA(2)
      Do k=1,dimA(3)
         tensO(i,j,k) = ABS(tensA(i,j,k))
      end Do
    end Do
  end Do

  !Find max element in tensO
  max_element = MAXVAL(tensO)

 END FUNCTION max_element_3D



 !Finds maximum ABS element in a complex array tensA -- 2D version
 FUNCTION max_element_2D(tensA) result(max_element)

  COMPLEX(KIND=DP), INTENT(IN) :: tensA(:,:) !input array
  REAL(KIND=DP) :: max_element !result

  !Temp array tensO = ABS(tensA)
  REAL(KIND=DP), ALLOCATABLE :: tensO(:,:)

  !Indices && dims
  INTEGER :: dimA(2), i, j

  !Get dims
  dimA = shape(tensA)
  
  !Transfer tensA to real array tensO = ABS(tensA)
  ALLOCATE(tensO(dimA(1), dimA(2)))
  Do i=1,dimA(1)
    Do j=1,dimA(2)
       tensO(i,j) = ABS(tensA(i,j))
    end Do
  end Do

  !Find max element in tensO
  max_element = MAXVAL(tensO)

 END FUNCTION max_element_2D



 !Finds maximum ABS element in a complex array tensA -- 1D version
 FUNCTION max_element_1D(tensA) result(max_element)

  COMPLEX(KIND=DP), INTENT(IN) :: tensA(:) !input array
  REAL(KIND=DP) :: max_element !result

  !Temp array tensO = ABS(tensA)
  REAL(KIND=DP), ALLOCATABLE :: tensO(:)

  !Indices && dims
  INTEGER :: dimA(1), i, j

  !Get dims
  dimA = shape(tensA)
  
  !Transfer tensA to real array tensO = ABS(tensA)
  ALLOCATE(tensO(dimA(1)))
  Do i=1,dimA(1)
     tensO(i) = ABS(tensA(i))
  end Do

  !Find max element in tensO
  max_element = MAXVAL(tensO)

 END FUNCTION max_element_1D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (6) CONCATENATE/SPLIT/INVERT ARRAYS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Concatenate vecA, vecB into vec !!!
 SUBROUTINE concatenate_vecs(vec, vecA, vecB)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: vec(:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: vecA(:)
  COMPLEX(KIND=DP),              INTENT(IN)    :: vecB(:)

  INTEGER :: Na, Nb

  !! Sizes of vecA, vecB
  Na = SIZE(vecA); Nb = SIZE(vecB)

  !! Construct vec = {vecA, vecB}
  ALLOCATE(vec(Na+Nb))
  vec(1    : Na)    = vecA
  vec(Na+1 : Na+Nb) = vecB

 END SUBROUTINE concatenate_vecs




 !!! Split vec into vecA, vecB of intended sizes Na, Nb !!!
 SUBROUTINE decatenate_vecs(vecA, vecB, vec, Na, Nb)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT)    :: vecA(:)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT)    :: vecB(:)
  COMPLEX(KIND=DP),              INTENT(IN)       :: vec(:)
  INTEGER,                       INTENT(IN)       :: Na, Nb

  !! Extract vecA, vecB from vec
  CALL copyTens(vecA,  vec(1    : Na))
  CALL copyTens(vecB,  vec(Na+1 : Na+Nb)) 

 END SUBROUTINE decatenate_vecs




 !!! Invert an array (swap beginning and end) !!!
 SUBROUTINE INVERT_array(Arr_in, Arr_out)

  LOGICAL, ALLOCATABLE, INTENT(INOUT)           :: Arr_in(:)
  LOGICAL, ALLOCATABLE, INTENT(INOUT), OPTIONAL :: Arr_out(:)

  !temp variable
  LOGICAL, ALLOCATABLE :: temp_arr(:)

  !indices
  INTEGER :: N_sites
  INTEGER :: site, swap_site, i

  !allocate an empty temp array
  N_sites = SIZE(Arr_in)
  ALLOCATE(temp_arr(N_sites))
   
  !Swap Array(site) and Array(N_sites - site + 1)
  DO site=1,N_sites
     temp_arr(site) = Arr_in(N_sites - site + 1)
  end DO

  !write to output
  IF(PRESENT(Arr_out)) THEN
     ALLOCATE(Arr_out(N_sites))
     Arr_out = temp_arr
  ELSE
     DEALLOCATE(Arr_in)
     ALLOCATE(Arr_in(N_sites))
     Arr_in = temp_arr
  end IF

  !clear memory
  DEALLOCATE(temp_arr)

 END SUBROUTINE INVERT_array


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (7) PURIFICATION ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! If tensA = meant to be Hermitian, we can purify it by removing the non-Hermitian part
 SUBROUTINE purify_hermitian_tens(tensA)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensA(:,:)

  !Local HC copy of tensA
  COMPLEX(KIND=DP), ALLOCATABLE :: tensA_HC(:,:)

  !Create local HC copy of tensA
  CALL copyTens(tensA_HC, tensA)
  tensA_HC = CONJG(TRANSPOSE(tensA_HC))

  !Purify tensA
  tensA = 0.5D0*(tensA + tensA_HC)

 END SUBROUTINE purify_hermitian_tens




 !!! If tensA = meant to be symmetric, we can purify it by removing the non-symmetric part
 SUBROUTINE purify_symmetric_tens(tensA)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensA(:,:)

  !! Local transposed copy of tensA
  COMPLEX(KIND=DP), ALLOCATABLE :: tensA_Transp(:,:)

  !! Create transposed copy of tensA
  CALL copyTens(tensA_Transp, tensA)
  tensA_Transp = TRANSPOSE(tensA_Transp)

  !! Purify tensA
  tensA = 0.5D0*(tensA + tensA_Transp)

 END SUBROUTINE purify_symmetric_tens




 !!! If tensA = meant to be positive, we can purify it by removing near-zero && slightly negative evals (see Orus review, p.45)
 SUBROUTINE purify_positive_tens(tensA, eps0)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT)        :: tensA(:,:)
  REAL(KIND=DP),                 INTENT(IN), OPTIONAL :: eps0

  !dims && indices
  INTEGER :: dimA(2)
  REAL(KIND=DP) :: eps

  !Set eps
  eps = OptArg(eps0, 2.0d-08)

  !Get dims of tensA
  dimA = shape(tensA)

  !Check if tensA = square matrix
  IF(dimA(1) .NE. dimA(2)) RETURN

  !Transform tensA by adding an infinitesimal constant proportional to the identity matrix
  tensA = tensA + eps * matEye(dimA(1))

 END SUBROUTINE purify_positive_tens



 !!! Auxiliary routine -- remove indices corresponding to zero evals !!!
 SUBROUTINE removeZeroEvals(evecsr, evecsl, evals)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: evecsr(:,:), evecsl(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: evals(:)

  COMPLEX(KIND=DP), ALLOCATABLE :: tempR(:,:), tempL(:,:), tempEv(:)
  
  INTEGER :: dimR(2), dimL(2), dimEv, chi
  INTEGER :: o, i
  REAL(KIND=DP), PARAMETER :: cutoff = 1.0d-14

  dimR = shape(evecsr); dimL = shape(evecsl); dimEv = SIZE(evals)
  
  !Count the num of nonzero evals
  chi = COUNT(abs(evals) .GT. cutoff)
  ALLOCATE(tempR(dimR(1), chi), tempEv(chi), tempL(chi, dimL(2)))

  !Init counter of nonzero lambdas
  o=0
 
  !Loop over all Lambdas, select only those that are nonzero
  DO i=1,dimEv
     IF(abs(evals(i)) .GT. cutoff) THEN 
        o=o+1
        tempR(:,o) = evecsr(:,i)
        tempEv(o) =  evals(i)
        tempL(o,:) = evecsl(i,:)
     end IF
  end DO

  !Allocate new evec/eval matrices
  DEALLOCATE(evecsr, evecsl, evals)
  ALLOCATE(evecsr(dimR(1), chi), evals(chi), evecsl(chi, dimL(2)))
  evecsr = tempR; evecsl = tempL; evals = tempEv

 END SUBROUTINE removeZeroEvals

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (8) ROUNDING/FILTERING ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Round PEPS site tensor !!!
 SUBROUTINE roundPEPS(tensA, roundFac)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:,:)
  REAL(KIND=DP),    INTENT(IN)    :: roundFac !rounding factor

  !dims && indices
  INTEGER :: dimA(5), i, j, k

  dimA = shape(tensA)

  !Round off
  DO i=1,dimA(1)
    DO j=1,dimA(2)
      DO k=1,dimA(3)
         CALL roundMatrix(tensA(i,j,k,:,:), roundFac)
      end DO
    end DO
  end DO

 END SUBROUTINE roundPEPS



 !!! Round MPO site tensor !!!
 SUBROUTINE roundMpo(tensA, roundFac)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:)
  REAL(KIND=DP),    INTENT(IN)    :: roundFac   !rounding factor

  !dims && indices
  INTEGER :: dimA(4), i, j

  dimA = shape(tensA)

  !Round off
  DO i=1,dimA(1)
    DO j=1,dimA(2)
       CALL roundMatrix(tensA(i,j,:,:), roundFac)
    end DO
  end DO

 END SUBROUTINE roundMpo



 !!! Round MPS site tensor !!!
 SUBROUTINE roundMps(tensA, roundFac)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:)
  REAL(KIND=DP),    INTENT(IN)    :: roundFac   !rounding factor

  !dims && indices
  INTEGER :: dimA(3), i

  dimA = shape(tensA)

  !Round off
  DO i=1,dimA(1)
     CALL roundMatrix(tensA(i,:,:), roundFac)
  end DO

 END SUBROUTINE roundMps



 !!! Round off values in a matrix !!!
 SUBROUTINE roundMatrix(tensA, roundFac) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:)
  REAL(KIND=DP),    INTENT(IN)    :: roundFac   !rounding factor

  !dims && indices
  INTEGER :: dimA(2), i, j

  dimA = shape(tensA)

  !Round off
  DO i=1,dimA(1)
    DO j=1,dimA(2)
       tensA(i,j) = roundValue(tensA(i,j), roundFac)
    end DO
  end DO
  
 END SUBROUTINE roundMatrix



 !!! Round off values in a vector array !!!
 SUBROUTINE roundArray(tensA, roundFac) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:)
  REAL(KIND=DP),    INTENT(IN)    :: roundFac   !rounding factor

  !dims && indices
  INTEGER :: dimA(1), i

  dimA = shape(tensA)

  !Round off
  DO i=1,dimA(1)
     tensA(i) = roundValue(tensA(i), roundFac)     
  end DO
  
 END SUBROUTINE roundArray



 !!! Round off a value !!!
 FUNCTION roundValue(valueIn, roundIn)

  COMPLEX(KIND=DP), INTENT(IN) :: valueIn      !input value
  REAL(KIND=DP),    INTENT(IN) :: roundIn      !rounding factor
  COMPLEX(KIND=DP)             :: roundValue   !result

  REAL(KIND=DP) :: roundFac
  REAL(KIND=DP) :: valueIn_Re, valueIn_Im
  REAL(KIND=DP) :: roundValue_Re, roundValue_Im

  !Adjust roundFac
  CALL adjust_round_factor(roundFac, roundIn)

  !get Re && Im parts of valueIn
  valueIn_Re = REAL(valueIn) 
  valueIn_Im = AIMAG(valueIn) 

  !round off Re && Im parts separately
  roundValue_Re = REAL(INT(valueIn_Re*RoundFac + 0.5d0), KIND=DP)/RoundFac
  roundValue_Im = REAL(INT(valueIn_Im*RoundFac + 0.5d0), KIND=DP)/RoundFac

  !combine Re && Im together
  roundValue = CMPLX(roundValue_Re, roundValue_Im, KIND=DP)

 END FUNCTION roundValue




 !!! Adjust rounding factor !!!
 SUBROUTINE adjust_round_factor(roundFac, roundIn)

  REAL(KIND=DP), INTENT(OUT) :: roundFac
  REAL(KIND=DP), INTENT(IN)  :: roundIn

  IF(roundIn .LT. 1.0d0) THEN
     roundFac = 1.0d0/roundIn
  ELSE
     roundFac = roundIn
  end IF

 END SUBROUTINE adjust_round_factor



 !filter peps
 SUBROUTINE filterPeps(tensA, cutoff_in)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:,:)
  REAL(KIND=DP), INTENT(IN), OPTIONAL :: cutoff_in

  INTEGER :: dimA(5), i, j

  dimA = shape(tensA)

  IF(PRESENT(cutoff_in)) THEN
     DO i=1,dimA(1)
       DO j=1,dimA(2)
          CALL filterTens(tensA(i,j,:,:,:), cutoff_in)
       end DO
     end DO
  ELSE
     DO i=1,dimA(1)
       DO j=1,dimA(2)
          CALL filterTens(tensA(i,j,:,:,:))
       end DO
     end DO
  end IF

 END SUBROUTINE filterPeps


 !filter mpo
 SUBROUTINE filterMpo(tensA, cutoff_in)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:,:)
  REAL(KIND=DP), INTENT(IN), OPTIONAL :: cutoff_in

  INTEGER :: dimA(4), i, j

  dimA = shape(tensA)

  IF(PRESENT(cutoff_in)) THEN
     DO i=1,dimA(1)
       DO j=1,dimA(2)
          CALL filterTens(tensA(i,j,:,:), cutoff_in)
       end DO
     end DO
  ELSE
     DO i=1,dimA(1)
       DO j=1,dimA(2)
          CALL filterTens(tensA(i,j,:,:))
       end DO
     end DO
  end IF

 END SUBROUTINE filterMpo




 !filter mps
 SUBROUTINE filterMps(tensA, cutoff_in)

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:,:)
  REAL(KIND=DP), INTENT(IN), OPTIONAL :: cutoff_in

  INTEGER :: dimA(3), i

  dimA = shape(tensA)

  IF(PRESENT(cutoff_in)) THEN
     DO i=1,dimA(1)
        CALL filterTens(tensA(i,:,:), cutoff_in)
     end DO
  ELSE
     DO i=1,dimA(1)
        CALL filterTens(tensA(i,:,:))
     end DO
  end IF

 END SUBROUTINE filterMps





 !filter out small values in a matrix
 SUBROUTINE filterMatrix(tensA, cutoff_in) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:,:)
  REAL(KIND=DP), INTENT(IN), OPTIONAL :: cutoff_in

  REAL(KIND=DP) :: cutoff 
  INTEGER :: dimA(2), i, j

  IF(PRESENT(cutoff_in)) THEN
     cutoff = cutoff_in
  ELSE
     cutoff = 1.0d-08
  end IF

  dimA = shape(tensA)

  DO i=1,dimA(1)
    DO j=1,dimA(2)
       tensA(i,j) = filterValue(tensA(i,j), cutoff)       
    end DO
  end DO
  
 END SUBROUTINE filterMatrix





 !filter out small values in a vector array
 SUBROUTINE filterArray(tensA, cutoff_in) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: tensA(:)
  REAL(KIND=DP), INTENT(IN), OPTIONAL :: cutoff_in

  REAL(KIND=DP) :: cutoff 
  INTEGER :: dimA(1), i

  IF(PRESENT(cutoff_in)) THEN
     cutoff = cutoff_in
  ELSE
     cutoff = 1.0d-08
  end IF

  dimA = shape(tensA)

  DO i=1,dimA(1)
     tensA(i) = filterValue(tensA(i), cutoff)      
  end DO
  
 END SUBROUTINE filterArray





 !filter out small values
 FUNCTION filterValue(valueIn, cutoff_in)

  COMPLEX(KIND=DP), INTENT(IN) :: valueIn !input value
  REAL(KIND=DP), INTENT(IN), OPTIONAL :: cutoff_in !filtering cutoff
  COMPLEX(KIND=DP) :: filterValue !result

  REAL(KIND=DP) :: valueIn_Re, valueIn_Im
  REAL(KIND=DP) :: valueOUT_Re, valueOUT_Im
  REAL(KIND=DP) :: cutoff 

  !set cutoff
  IF(PRESENT(cutoff_in)) THEN
     cutoff = cutoff_in
  ELSE
     cutoff = 1.0d-08
  end IF

  !get Re && Im parts of valueIn
  valueIn_Re = REAL(valueIn) 
  valueIn_Im = AIMAG(valueIn) 

  !filter out Re && Im parts separately
  IF(ABS(valueIn_Re) .LT. cutoff) THEN
     valueOUT_Re = 0.0d0
  ELSE
     valueOUT_Re = valueIn_Re
  end IF 

  IF(ABS(valueIn_Im) .LT. cutoff) THEN
     valueOUT_Im = 0.0d0
  ELSE
     valueOUT_Im = valueIn_Im
  end IF

  !combine Re && Im together
  filterValue = CMPLX(valueOUT_Re, valueOUT_Im, KIND=DP)

 END FUNCTION filterValue

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (9) ADD NOISE to tensor/value !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Add noise to 5D tensor !!!
 FUNCTION add_noise_5D(tensA) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:,:)

  INTEGER :: dimA(5), i, j, k, l, m

  !! Get dims of tensA
  dimA = shape(tensA)

  !! Copy tensA to tensO
  CALL copyTens(tensO, tensA)

  !! Add noise to tensO
  DO i=1,dimA(1)
     DO j=1,dimA(2)
        DO k=1,dimA(3)
           DO l=1,dimA(4)
              DO m=1,dimA(5)
                 tensO(i,j,k,l,m) = add_noise(tensO(i,j,k,l,m))
              end DO
           end DO
        end DO
     end DO
  end DO

 END FUNCTION add_noise_5D




 !!! Add noise to 4D tensor !!!
 FUNCTION add_noise_4D(tensA) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:)

  INTEGER :: dimA(4), i, j, k, l

  !! Get dims of tensA
  dimA = shape(tensA)

  !! Copy tensA to tensO
  CALL copyTens(tensO, tensA)

  !! Add noise to tensO
  DO i=1,dimA(1)
     DO j=1,dimA(2)
        DO k=1,dimA(3)
           DO l=1,dimA(4)
              tensO(i,j,k,l) = add_noise(tensO(i,j,k,l))
           end DO
        end DO
     end DO
  end DO

 END FUNCTION add_noise_4D




 !!! Add noise to 3D tensor !!!
 FUNCTION add_noise_3D(tensA) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:)

  INTEGER :: dimA(3), i, j, k

  !! Get dims of tensA
  dimA = shape(tensA)

  !! Copy tensA to tensO
  CALL copyTens(tensO, tensA)

  !! Add noise to tensO
  DO i=1,dimA(1)
     DO j=1,dimA(2)
        DO k=1,dimA(3)
           tensO(i,j,k) = add_noise(tensO(i,j,k))
        end DO
     end DO
  end DO

 END FUNCTION add_noise_3D




 !!! Add noise to 2D tensor !!!
 FUNCTION add_noise_2D(tensA) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:)

  INTEGER :: dimA(2), i, j

  !! Get dims of tensA
  dimA = shape(tensA)

  !! Copy tensA to tensO
  CALL copyTens(tensO, tensA)

  !! Add noise to tensO
  DO i=1,dimA(1)
     DO j=1,dimA(2)
        tensO(i,j) = add_noise(tensO(i,j))
     end DO
  end DO

 END FUNCTION add_noise_2D




 !!! Add noise to 1D tensor !!!
 FUNCTION add_noise_1D(tensA) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:)

  INTEGER :: dimA(1), i

  !! Get dims of tensA
  dimA = shape(tensA)

  !! Copy tensA to tensO
  CALL copyTens(tensO, tensA)

  !! Add noise to tensO
  DO i=1,dimA(1)
     tensO(i) = add_noise(tensO(i))
  end DO

 END FUNCTION add_noise_1D



 !!! Add noise to value !!!
 FUNCTION add_noise_Val(val) result(valOut)

  COMPLEX(KIND=DP), INTENT(IN) :: val
  COMPLEX(KIND=DP)             :: valOut

  !! Real/Imag parts
  REAL(KIND=DP) :: val_Re,    val_Im
  REAL(KIND=DP) :: valOUT_Re, valOUT_Im

  !! Rand seed
  REAL    :: cpu_times 
  INTEGER :: seed 

  !! Noise level
  REAL(KIND=DP), PARAMETER :: eps = 1.0d-02

  !! Start random num generator
  CALL cpu_time(cpu_times)
  seed = NINT(cpu_times)
  CALL srand(seed)

  !! Get Re/Im part of val
  val_Re = REAL(val) 
  val_Im = AIMAG(val)

  !! Add noise
  valOut_Re = val_Re + eps * val_Re * rand()
  valOut_Im = val_Im + eps * val_Im * rand()

  !! Combine into valOut
  valOut = CMPLX(valOUT_Re, valOUT_Im, KIND=DP)

 END FUNCTION add_noise_Val

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE array_utility
