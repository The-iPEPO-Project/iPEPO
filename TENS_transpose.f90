MODULE TENS_transpose

 USE basic_functions
 USE error_handling
 USE array_utility
 USE TENS_helper

 implicit none

 interface TensTRANSPOSE
    module procedure transpose_6D
    module procedure transpose_4D
    module procedure transpose_3D
    module procedure transpose_2D
 end interface TensTRANSPOSE

 interface TensROTATE
    module procedure rotate_5D
    module procedure rotate_4D
    module procedure rotate_3D
 end interface TensROTATE

 private !! hides all items not listed in public statement 
 public TensTRANSPOSE, TensROTATE

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TRANSPOSE a general tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Transpose 4D tensor !!!
 FUNCTION transpose_6D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:,:,:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: legs
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:,:,:)

  INTEGER :: dimA(6), i, j, k, l
  
  dimA = shape(tensA)

  SELECT CASE(legs)
  CASE('12')

     ALLOCATE(tensO(dimA(2), dimA(1), dimA(3), dimA(4), dimA(5), dimA(6)))

     DO i=1,dimA(3)
       DO j=1,dimA(4)
         DO k=1,dimA(5)
           DO l=1,dimA(6)
              tensO(:,:,i,j,k,l) = TRANSPOSE(tensA(:,:,i,j,k,l))
           end DO
         end DO
       end DO
     end DO

  CASE DEFAULT
     CALL invalid_flag("transpose_6D -- invalid legs ", legs)
  end SELECT
 
 END FUNCTION transpose_6D





 !!! Transpose 4D tensor !!!
 FUNCTION transpose_4D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: legs
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:)

  INTEGER :: dimA(4), i, j
  
  dimA = shape(tensA)

  SELECT CASE(legs)
  CASE('12')

     ALLOCATE(tensO(dimA(2), dimA(1), dimA(3), dimA(4)))

     DO i=1,dimA(3)
       DO j=1,dimA(4)
          tensO(:,:,i,j) = TRANSPOSE(tensA(:,:,i,j))
       end DO
     end DO

  CASE('13')

     ALLOCATE(tensO(dimA(3), dimA(2), dimA(1), dimA(4)))

     DO i=1,dimA(2)
       DO j=1,dimA(4)
          tensO(:,i,:,j) = TRANSPOSE(tensA(:,i,:,j))
       end DO
     end DO

  CASE('14')

     ALLOCATE(tensO(dimA(4), dimA(2), dimA(3), dimA(1)))

     DO i=1,dimA(2)
       DO j=1,dimA(3)
          tensO(:,i,j,:) = TRANSPOSE(tensA(:,i,j,:))
       end DO
     end DO

  CASE('23')

     ALLOCATE(tensO(dimA(1), dimA(3), dimA(2), dimA(4)))

     DO i=1,dimA(1)
       DO j=1,dimA(4)
          tensO(i,:,:,j) = TRANSPOSE(tensA(i,:,:,j))
       end DO
     end DO

  CASE('24')

     ALLOCATE(tensO(dimA(1), dimA(4), dimA(3), dimA(2)))

     DO i=1,dimA(1)
       DO j=1,dimA(3)
          tensO(i,:,j,:) = TRANSPOSE(tensA(i,:,j,:))
       end DO
     end DO

  CASE('34')

     ALLOCATE(tensO(dimA(1), dimA(2), dimA(4), dimA(3)))

     DO i=1,dimA(1)
       DO j=1,dimA(2)
          tensO(i,j,:,:) = TRANSPOSE(tensA(i,j,:,:))
       end DO
     end DO

  CASE DEFAULT
     CALL invalid_flag("transpose_4D -- invalid legs ", legs)
  end SELECT
 
 END FUNCTION transpose_4D





 !!! Transpose 3D tensor !!!
 FUNCTION transpose_3D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: legs
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:)

  INTEGER :: dimA(3), i, j
  
  dimA = shape(tensA)

  SELECT CASE (legs)
  CASE('23')

     ALLOCATE(tensO(dimA(1), dimA(3), dimA(2)))

     DO i=1,dimA(1)
        tensO(i,:,:) = TRANSPOSE(tensA(i,:,:))
     end DO

  CASE('12') 

     ALLOCATE(tensO(dimA(2), dimA(1), dimA(3)))

     DO i=1,dimA(3)
        tensO(:,:,i) = TRANSPOSE(tensA(:,:,i))
     end DO

  CASE('13')

     ALLOCATE(tensO(dimA(3), dimA(2), dimA(1)))

     DO i=1,dimA(2)
        tensO(:,i,:) = TRANSPOSE(tensA(:,i,:))
     end DO

  CASE DEFAULT
     CALL invalid_flag("transpose_3D -- invalid legs ", legs)
  end SELECT
 
 END FUNCTION transpose_3D



 !!! Transpose 2D tensor (HC or Simple transpose) !!!
 FUNCTION transpose_2D(tensA, HCflag) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)           :: tensA(:,:)
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HCflag
  COMPLEX(KIND=DP), ALLOCATABLE          :: tensO(:,:)

  INTEGER :: dimA(2)
  
  dimA = shape(tensA)

  ALLOCATE(tensO(dimA(2), dimA(1)))

  IF(PRESENT(HCflag)) THEN
     tensO(:,:) = CONJG(TRANSPOSE(tensA))
  ELSE
     tensO(:,:) = TRANSPOSE(tensA) 
  end IF

 END FUNCTION transpose_2D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ROTATE tensor by +PI/2, -PI/2, PI angles !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Rotate 5D tensor: rotateFlag = [not present = CW], CCW, PI !!!
 FUNCTION rotate_5D(tensA, rotFlagIn) result(tensO)
  
  COMPLEX(KIND=DP), INTENT(IN)           :: tensA(:,:,:,:,:)
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: rotFlagIn
  COMPLEX(KIND=DP), ALLOCATABLE          :: tensO(:,:,:,:,:)

  !! Dims && indices
  INTEGER           :: dimA(5), i, j, k
  CHARACTER(LEN=16) :: rotFlag

  !! Set default value of rotateFlag
  rotFlag = OptArg(rotFlagIn, '+PI/2')

  !! Get dims
  dimA=shape(tensA)

  !! Rotate tensor
  SELECT CASE(rotFlag)
  CASE('+PI/2', 'CW', 'CW +PI/2')

       !! Rotate by +PI/2 (clockwise) 
       ALLOCATE(tensO(dimA(1), dimA(5), dimA(4), dimA(2), dimA(3)))

       DO i=1,dimA(4)
         DO j=1,dimA(5)
            tensO(:,j,i,:,:) = tensA(:,:,:,i,j)
         end DO
       end DO

  CASE('-PI/2', 'CCW', 'CCW -PI/2')

       !! Rotate by -PI/2 (counter-clockwise) 
       ALLOCATE(tensO(dimA(1), dimA(4), dimA(5), dimA(3), dimA(2)))

       DO i=1,dimA(2)
         DO j=1,dimA(3)
            tensO(:,:,:,j,i) = tensA(:,i,j,:,:)
         end DO
       end DO

  CASE('PI')

       !! Rotate by PI 
       ALLOCATE(tensO(dimA(1), dimA(3), dimA(2), dimA(5), dimA(4)))

       DO k=1,dimA(1)
         DO i=1,dimA(2)
           DO j=1,dimA(3)
              tensO(k,j,i,:,:) = TRANSPOSE(tensA(k,i,j,:,:))
           end DO
         end DO
       end DO

  CASE DEFAULT
       CALL invalid_flag("rotate_5D -- invalid rotFlag ", rotFlag)
  end SELECT

 END FUNCTION rotate_5D





 !!! Rotate 4D tensor: rotateFlag = [not present = CW], CCW, PI !!!
 FUNCTION rotate_4D(tensA, rotFlagIn) result(tensO)
  
  COMPLEX(KIND=DP), INTENT(IN)           :: tensA(:,:,:,:)
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: rotFlagIn
  COMPLEX(KIND=DP), ALLOCATABLE          :: tensO(:,:,:,:)

  !! Dims && indices
  INTEGER           :: dimA(4), i, j
  CHARACTER(LEN=16) :: rotFlag

  !! Set default value of rotateFlag
  rotFlag = OptArg(rotFlagIn, '+PI/2')

  !! Get dims
  dimA=shape(tensA)

  !! Rotate tensor
  SELECT CASE(rotFlag)
  CASE('+PI/2', 'CW', 'CW +PI/2')

       !! Rotate by +PI/2 (clockwise) 
       ALLOCATE(tensO(dimA(4), dimA(3), dimA(1), dimA(2)))

       DO i=1,dimA(3)
         DO j=1,dimA(4)
            tensO(j,i,:,:) = tensA(:,:,i,j)
         end DO
       end DO

  CASE('-PI/2', 'CCW', 'CCW -PI/2')

       !! Rotate by -PI/2 (counter-clockwise) 
       ALLOCATE(tensO(dimA(3), dimA(4), dimA(2), dimA(1)))

       DO i=1,dimA(1)
         DO j=1,dimA(2)
            tensO(:,:,j,i) = tensA(i,j,:,:)
         end DO
       end DO

  CASE('PI')

       !! Rotate by PI 
       ALLOCATE(tensO(dimA(2), dimA(1), dimA(4), dimA(3)))

       DO i=1,dimA(1)
         DO j=1,dimA(2)
            tensO(j,i,:,:) = TRANSPOSE(tensA(i,j,:,:))
         end DO
       end DO

  CASE DEFAULT
       CALL invalid_flag("rotate_4D -- invalid rotFlag ", rotFlag)
  end SELECT

 END FUNCTION rotate_4D




 !!! Rotate 3D tensor by PI/2 clockwise 
 !!! -- add an extra leg at position specified by [extra_leg] to convert 3D into 4D
 !!! -- then rotate by an angle specified by rotFlag
 FUNCTION rotate_3D(tensA, rotFlag, extra_leg) result(tensO)
  
  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:,:)
  CHARACTER(LEN=*),  INTENT(IN)  :: rotFlag, extra_leg
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:,:,:)

  !! Dims && indices
  INTEGER :: dimA(3), i

  !! Get dims
  dimA=shape(tensA)

  !! Rotate tensor
  SELECT CASE(rotFlag)
  CASE('+PI/2', 'CW', 'CW +PI/2')

       SELECT CASE(extra_leg)
       CASE('2')  

            !! Rotate by +PI/2 (clockwise) 
            ALLOCATE(tensO(dimA(3), dimA(2), dimA(1), 1))
            DO i=1,dimA(1)
               tensO(:,:,i,1) = TRANSPOSE(tensA(i,:,:))
            end DO

       CASE('1') 

            !! Rotate by PI/2 (clockwise) 
            ALLOCATE(tensO(dimA(3), dimA(2), 1, dimA(1)))
            DO i=1,dimA(1)
               tensO(:,:,1,i) = TRANSPOSE(tensA(i,:,:))
            end DO

       CASE DEFAULT
            CALL invalid_flag("rotate_3D -- invalid extra_leg ", extra_leg)
       end SELECT


  CASE('-PI/2', 'CCW', 'CCW -PI/2')

       SELECT CASE(extra_leg)
       CASE('2')  

            !! Rotate by -PI/2 (counterclockwise) 
            ALLOCATE(tensO(dimA(2), dimA(3), 1, dimA(1)))
            DO i=1,dimA(1)
               tensO(:,:,1,i) = tensA(i,:,:)
            end DO

       CASE('1') 

            !! Rotate by -PI/2 (counterclockwise) 
            ALLOCATE(tensO(dimA(2), dimA(3), dimA(1), 1))
            DO i=1,dimA(1)
               tensO(:,:,i,1) = tensA(i,:,:)
            end DO

       CASE DEFAULT
            CALL invalid_flag("rotate_3D -- invalid extra_leg ", extra_leg)
       end SELECT

  CASE('PI')

       SELECT CASE(extra_leg)
       CASE('1')  

            !! Rotate by PI
            ALLOCATE(tensO(1, dimA(1), dimA(3), dimA(2)))
            DO i=1,dimA(1)
               tensO(1,i,:,:) = TRANSPOSE(tensA(i,:,:))
            end DO

       CASE('2') 

            !! Rotate by PI
            ALLOCATE(tensO(dimA(1), 1, dimA(3), dimA(2)))
            DO i=1,dimA(1)
               tensO(i,1,:,:) = TRANSPOSE(tensA(i,:,:))
            end DO

       CASE DEFAULT
            CALL invalid_flag("rotate_3D -- invalid extra_leg ", extra_leg)
       end SELECT

  CASE DEFAULT
       CALL invalid_flag("rotate_3D -- invalid rotFlag ", rotFlag)
  end SELECT

 END FUNCTION rotate_3D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE TENS_transpose
