MODULE TENS_helper

 USE basic_functions
 USE error_handling
 USE array_utility

 implicit none

 interface fuse_legs
    module procedure fuse_legs_6D_to_4D
    module procedure fuse_legs_5D_to_3D
    module procedure fuse_legs_5D_to_2D
    module procedure fuse_legs_4D_to_3D
    module procedure fuse_legs_4D_to_2D
    module procedure fuse_legs_4D_to_2D_oneside
    module procedure fuse_legs_3D_to_2D
 end interface fuse_legs

 interface split_legs
    module procedure split_legs_4D_to_6D
    module procedure split_legs_3D_to_5D
    module procedure split_legs_2D_to_5D
    module procedure split_legs_3D_to_4D
    module procedure split_legs_2D_to_4D
    module procedure split_legs_2D_to_4D_oneside
    module procedure split_legs_2D_to_3D
 end interface split_legs

 interface increase_rank
    module procedure increase_rank_1D_to_2D
    module procedure increase_rank_2D_to_3D
    module procedure increase_rank_3D_to_4D
 end interface increase_rank

 interface reduce_rank
    module procedure reduce_rank_4D
    module procedure reduce_rank_3D
    module procedure reduce_rank_2D
 end interface reduce_rank

 private !! hides all items not listed in public statement 
 public ICOM
 public fuse_legs, split_legs
 public increase_rank, reduce_rank

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ICOM -- calculate combo index for fused legs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Get combo index ICOM !!!
 RECURSIVE FUNCTION ICOM(o1, o2, dimO2, o3, dimO3, o4, dimO4) result(I_COM)

  INTEGER, INTENT(IN)           :: o1, o2, dimO2
  INTEGER, INTENT(IN), OPTIONAL :: o3, dimO3
  INTEGER, INTENT(IN), OPTIONAL :: o4, dimO4
  INTEGER                       :: I_COM 

  IF(PRESENT(o4) .AND. PRESENT(dimO4)) THEN

     !! Quadruple index
     I_COM = o4 + (ICOM(o1, o2, dimO2, o3, dimO3) - 1) * dimO4

  ELSEIF(PRESENT(o3) .AND. PRESENT(dimO3)) THEN

     !! Triple index
     I_COM = o3 + (ICOM(o1, o2, dimO2) - 1) * dimO3
  ELSE

     !! Double index
     I_COM = o2 + (o1-1) * dimO2
  end IF

 END FUNCTION ICOM




 !!! Get combo index ICOM -- new (untested) version !!!
 FUNCTION ICOM_NEW(o1, o2, dimO2, o3, dimO3, o4, dimO4) result(ICOM)

  INTEGER, INTENT(IN)           :: o1, o2, dimO2
  INTEGER, INTENT(IN), OPTIONAL :: o3, dimO3
  INTEGER, INTENT(IN), OPTIONAL :: o4, dimO4
  INTEGER                       :: ICOM 

  IF(PRESENT(o4) .AND. PRESENT(dimO4)) THEN

     !! Quadruple index
     ICOM = o4  +  (o3-1) * dimO4  & 
                                   & +  (o2-1) * dimO3 * dimO4  &
                                   & +  (o1-1) * dimO2 * dimO3 * dimO4   

  ELSEIF(PRESENT(o3) .AND. PRESENT(dimO3)) THEN

     !! Triple index
     ICOM = o3  +  (o2-1) * dimO3  & 
                                   & +  (o1-1) * dimO2 * dimO3   
  ELSE

     !! Double index
     ICOM = o2  +  (o1-1) * dimO2
  end IF

 END FUNCTION ICOM_NEW

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FUSE tensor legs inside a loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Fuse legs 6D --> 4D !!!
 SUBROUTINE fuse_legs_6D_to_4D(tensO, tensA, legs, aTOT, a1, a2, a3)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:,:,:)       
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:,:,:,:)   
  CHARACTER(LEN=*),              INTENT(IN)    :: legs      
  INTEGER,                       INTENT(IN)    :: aTOT, a1, a2, a3   

  SELECT CASE(legs)
  CASE('1,2,345,6') 
       tensO(:, :, aTOT, :) = tensA(:, :, a1, a2, a3, :)

  CASE('1,2,5,346')
       tensO(:, :, :, aTOT) = tensA(:, :, a1, a2, :, a3)

  CASE('1,2,356,4')
       tensO(:, :, aTOT, :) = tensA(:, :, a1, :, a2, a3)

  CASE('1,2,3,456')
       tensO(:, :, :, aTOT) = tensA(:, :, :, a1, a2, a3)

  CASE DEFAULT
       CALL invalid_flag("fuse_legs_6D_to_4D -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE fuse_legs_6D_to_4D



 !!! Fuse legs 5D --> 2D !!!
 SUBROUTINE fuse_legs_5D_to_2D(tensO, tensA, legs, aTOT, a1, a2, a3, a4)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:)         
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:,:,:)       
  CHARACTER(LEN=*),              INTENT(IN)    :: legs            
  INTEGER,                       INTENT(IN)    :: aTOT, a1, a2, a3, a4  

  SELECT CASE(legs)
  CASE('1234,5') 
       tensO(aTOT, :) = tensA(a1, a2, a3, a4, :)

  CASE('4,1235') 
       tensO(:, aTOT) = tensA(a1, a2, a3, :, a4) 

  CASE('1245,3')
       tensO(aTOT, :) = tensA(a1, a2, :, a3, a4)

  CASE('2,1345') 
       tensO(:, aTOT) = tensA(a1, :, a2, a3, a4)

  CASE('1,2345') 
       tensO(:, aTOT) = tensA(:, a1, a2, a3, a4)

  CASE DEFAULT
       CALL invalid_flag("fuse_legs_5D_to_2D -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE fuse_legs_5D_to_2D



 !!! Fuse legs 5D --> 3D !!!
 SUBROUTINE fuse_legs_5D_to_3D(tensO, tensA, legs, aTOT, a1, a2, a3)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:,:)       
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:,:,:)   
  CHARACTER(LEN=*),              INTENT(IN)    :: legs      
  INTEGER,                       INTENT(IN)    :: aTOT, a1, a2, a3   

  SELECT CASE(legs)
  CASE('1,234,5') 
       tensO(:, aTOT, :) = tensA(:, a1, a2, a3, :)

  CASE('1,4,235')
       tensO(:, :, aTOT) = tensA(:, a1, a2, :, a3)

  CASE('1,245,3')
       tensO(:, aTOT, :) = tensA(:, a1, :, a2, a3)

  CASE('1,2,345')
       tensO(:, :, aTOT) = tensA(:, :, a1, a2, a3)

  CASE('123,4,5') 
       tensO(aTOT, :, :) = tensA(a1, a2, a3, :, :)

  CASE('145,2,3')
       tensO(aTOT, :, :) = tensA(a1, : , :, a2, a3)

  CASE DEFAULT
       CALL invalid_flag("fuse_legs_5D_to_3D -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE fuse_legs_5D_to_3D






 !!! Fuse legs 4D --> 3D !!!
 SUBROUTINE fuse_legs_4D_to_3D(tensO, tensA, legs, aTOT, a1, a2)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:,:)       
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:,:)   
  CHARACTER(LEN=*),              INTENT(IN)    :: legs      
  INTEGER,                       INTENT(IN)    :: aTOT, a1, a2  

  SELECT CASE(legs)
  CASE('12,3,4') 
       tensO(aTOT, :, :) = tensA(a1, a2, :, :)

  CASE('1,2,34')
       tensO(:, :, aTOT) = tensA(:, :, a1, a2)

  CASE('2,3,14')
       tensO(:, :, aTOT) = tensA(a1, :, :, a2)

  CASE('2,13,4')
       tensO(:, aTOT, :) = tensA(a1, :, a2, :)

  CASE('1,3,24') 
       tensO(:, :, aTOT) = tensA(:, a1, :, a2)

  CASE('1,23,4')
       tensO(:, aTOT, :) = tensA(:, a1, a2, :)

  CASE DEFAULT
       CALL invalid_flag("fuse_legs_4D_to_3D -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE fuse_legs_4D_to_3D





 !!! Fuse legs 4D --> 2D !!!
 SUBROUTINE fuse_legs_4D_to_2D(tensO, tensA, legs, aTOT, a1, a2, bTOT, b1, b2)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:)       
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:,:)   
  CHARACTER(LEN=*),              INTENT(IN)    :: legs      
  INTEGER,                       INTENT(IN)    :: aTOT, a1, a2 
  INTEGER,                       INTENT(IN)    :: bTOT, b1, b2 

  !! [a1,a2] --> aTOT, [b1,b2] --> bTOT
  SELECT CASE(legs)
  CASE('13,24') 
       tensO(aTOT, bTOT) = tensA(a1, b1, a2, b2)

  CASE('12,34')
       tensO(aTOT, bTOT) = tensA(a1, a2, b1, b2)

  CASE DEFAULT
       CALL invalid_flag("fuse_legs_4D_to_2D -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE fuse_legs_4D_to_2D




 !!! Fuse legs 4D --> 2D !!!
 SUBROUTINE fuse_legs_4D_to_2D_oneside(tensO, tensA, legs, aTOT, a1, a2, a3)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:)       
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:,:)   
  CHARACTER(LEN=*),              INTENT(IN)    :: legs      
  INTEGER,                       INTENT(IN)    :: aTOT, a1, a2, a3 

  !! [a1,a2,a3] --> aTOT
  SELECT CASE(legs)
  CASE('1,234') 
       tensO(:, aTOT) = tensA(:, a1, a2, a3)

  CASE DEFAULT
       CALL invalid_flag("fuse_legs_4D_to_2D_oneside -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE fuse_legs_4D_to_2D_oneside





 !!! Fuse legs 3D --> 2D !!!
 SUBROUTINE fuse_legs_3D_to_2D(tensO, tensA, legs, aTOT, a1, a2)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:)       
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:)   
  CHARACTER(LEN=*),              INTENT(IN)    :: legs      
  INTEGER,                       INTENT(IN)    :: aTOT, a1, a2  

  SELECT CASE(legs)
  CASE('12,3') 
       tensO(aTOT, :) = tensA(a1, a2, :)

  CASE('2,13')
       tensO(:, aTOT) = tensA(a1, :, a2)

  CASE('1,23')
       tensO(:, aTOT) = tensA(:, a1, a2)

  CASE DEFAULT
       CALL invalid_flag("fuse_legs_3D_to_2D -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE fuse_legs_3D_to_2D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SPLIT tensor legs inside a loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Split legs 4D --> 6D !!!
 SUBROUTINE split_legs_4D_to_6D(tensO, tensA, legs, oTOT, o1, o2, o3)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:,:,:,:,:)   
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:,:)      
  CHARACTER(LEN=*),              INTENT(IN)    :: legs        
  INTEGER,                       INTENT(IN)    :: oTOT, o1, o2, o3  

  SELECT CASE(legs)
  CASE('1,2,345,6')  
       tensO(:,:,o1,o2,o3,:) = tensA(:, :, oTOT, :)

  CASE('1,2,5,346') 
       tensO(:,:,o1,o2,:,o3) = tensA(:, :, :, oTOT)

  CASE('1,2,356,4') 
       tensO(:,:,o1,:,o2,o3) = tensA(:, :, oTOT, :)

  CASE('1,2,3,456') 
       tensO(:,:,:,o1,o2,o3) = tensA(:, :, :, oTOT)

  CASE DEFAULT
       CALL invalid_flag("split_legs_4D_to_6D -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE split_legs_4D_to_6D



 !!! Split legs 2D --> 5D !!!
 SUBROUTINE split_legs_2D_to_5D(tensO, tensA, legs, oTOT, o1, o2, o3, o4)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:,:,:,:)      
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:)           
  CHARACTER(LEN=*),              INTENT(IN)    :: legs           
  INTEGER,                       INTENT(IN)    :: oTOT, o1, o2, o3, o4 

  SELECT CASE(legs)
  CASE('1234,5')
       tensO(o1,o2,o3,o4,:) = tensA(oTOT, :)

  CASE('4,1235') 
       tensO(o1,o2,o3,:,o4) = tensA(:, oTOT)

  CASE('1245,3') 
       tensO(o1,o2,:,o3,o4) = tensA(oTOT, :)

  CASE('2,1345') 
       tensO(o1,:,o2,o3,o4) = tensA(:, oTOT)

  CASE('1,2345') 
       tensO(:,o1,o2,o3,o4) = tensA(:, oTOT)

  CASE DEFAULT
       CALL invalid_flag("split_legs_2D_to_5D -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE split_legs_2D_to_5D




 !!! Split legs 3D --> 5D !!!
 SUBROUTINE split_legs_3D_to_5D(tensO, tensA, legs, oTOT, o1, o2, o3)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:,:,:,:)   
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:)      
  CHARACTER(LEN=*),              INTENT(IN)    :: legs        
  INTEGER,                       INTENT(IN)    :: oTOT, o1, o2, o3  

  SELECT CASE(legs)
  CASE('1,234,5')  
       tensO(:,o1,o2,o3,:) = tensA(:, oTOT, :)

  CASE('1,4,235') 
       tensO(:,o1,o2,:,o3) = tensA(:, :, oTOT)

  CASE('1,245,3') 
       tensO(:,o1,:,o2,o3) = tensA(:, oTOT, :)

  CASE('1,2,345') 
       tensO(:,:,o1,o2,o3) = tensA(:, :, oTOT)

  CASE('123,4,5') 
       tensO(o1,o2,o3,:,:) = tensA(oTOT, :, :)

  CASE('145,2,3') 
       tensO(o1,:,:,o2,o3) = tensA(oTOT, :, :)

  CASE DEFAULT
       CALL invalid_flag("split_legs_3D_to_5D -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE split_legs_3D_to_5D




 !!! Split legs 3D --> 4D !!!
 SUBROUTINE split_legs_3D_to_4D(tensO, tensA, legs, oTOT, o1, o2)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:,:,:)   
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:,:)      
  CHARACTER(LEN=*),              INTENT(IN)    :: legs        
  INTEGER,                       INTENT(IN)    :: oTOT, o1, o2 

  SELECT CASE(legs)
  CASE('12,3,4') 
       tensO(o1,o2,:,:) = tensA(oTOT, :, :)

  CASE('1,2,34')
       tensO(:,:,o1,o2) = tensA(:, :, oTOT)

  CASE('2,3,14') 
       tensO(o1,:,:,o2) = tensA(:, :, oTOT)

  CASE('2,13,4') 
       tensO(o1,:,o2,:) = tensA(:, oTOT, :)

  CASE('1,3,24') 
       tensO(:,o1,:,o2) = tensA(:, :, oTOT)

  CASE('1,23,4') 
       tensO(:,o1,o2,:) = tensA(:, oTOT, :)

  CASE DEFAULT
       CALL invalid_flag("split_legs_3D_to_4D -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE split_legs_3D_to_4D




 !!! Split legs 2D --> 4D !!!
 SUBROUTINE split_legs_2D_to_4D(tensO, tensA, legs, aTOT, a1, a2, bTOT, b1, b2)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:,:,:)   
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:)      
  CHARACTER(LEN=*),              INTENT(IN)    :: legs        
  INTEGER,                       INTENT(IN)    :: aTOT, a1, a2 
  INTEGER,                       INTENT(IN)    :: bTOT, b1, b2 

  !! aTOT --> [a1,a2], bTOT --> [b1,b2]
  SELECT CASE(legs)
  CASE('13,24') 
        tensO(a1, b1, a2, b2) = tensA(aTOT, bTOT)

  CASE('12,34')
        tensO(a1, a2, b1, b2) = tensA(aTOT, bTOT)

  CASE DEFAULT
        CALL invalid_flag("split_legs_2D_to_4D -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE split_legs_2D_to_4D



 !!! Split legs 2D --> 4D !!!
 SUBROUTINE split_legs_2D_to_4D_oneside(tensO, tensA, legs, aTOT, a1, a2, a3)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:,:,:)   
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:)  
  CHARACTER(LEN=*),              INTENT(IN)    :: legs      
  INTEGER,                       INTENT(IN)    :: aTOT, a1, a2, a3 

  !! aTOT --> [a1,a2,a3]
  SELECT CASE(legs)
  CASE('1,234') 
       tensO(:, a1, a2, a3) = tensA(:, aTOT) 

  CASE DEFAULT
       CALL invalid_flag("split_legs_2D_to_4D_oneside -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE split_legs_2D_to_4D_oneside





 !!! Split legs 2D --> 3D !!!
 SUBROUTINE split_legs_2D_to_3D(tensO, tensA, legs, oTOT, o1, o2)

  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: tensO(:,:,:)   
  COMPLEX(KIND=DP),              INTENT(IN)    :: tensA(:,:)      
  CHARACTER(LEN=*),              INTENT(IN)    :: legs        
  INTEGER,                       INTENT(IN)    :: oTOT, o1, o2 

  SELECT CASE(legs)
  CASE('12,3') 
       tensO(o1, o2, :) = tensA(oTOT, :)

  CASE('2,13')
       tensO(o1, :, o2) = tensA(:, oTOT)

  CASE('1,23')
       tensO(:, o1, o2) = tensA(:, oTOT)

  CASE DEFAULT
       CALL invalid_flag("split_legs_2D_to_3D -- invalid legs ", legs)
  end SELECT

 END SUBROUTINE split_legs_2D_to_3D


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!!!!!! Increase tensor rank: insert an extra singleton leg at a specified position !!!!!!!!!!!!!!!!!

 !!! Convert 2D tensor into 3D tensor by inserting an extra singleton leg !!!
 FUNCTION increase_rank_1D_to_2D(tens_1d, which_leg) result(tens_2d)

  COMPLEX(KIND=DP), INTENT(IN) :: tens_1d(:)
  CHARACTER(LEN=*), INTENT(IN) :: which_leg

  COMPLEX(KIND=DP), ALLOCATABLE :: tens_2d(:,:)
  INTEGER                       :: dim2D(1)

  dim2D = shape(tens_1d)

  IF(which_leg .EQ. '1') THEN 

     ALLOCATE(tens_2d(1, dim2D(1)))
     tens_2d(1,:) = tens_1d

  ELSEIF(which_leg .EQ. '2') THEN

     ALLOCATE(tens_2d(dim2D(1), 1))
     tens_2d(:,1) = tens_1d

  ELSE
     CALL invalid_flag("Error: insert_leg has invalid value", which_leg)
  end IF

 END FUNCTION increase_rank_1D_to_2D




 !!! Convert 2D tensor into 3D tensor by inserting an extra singleton leg !!!
 FUNCTION increase_rank_2D_to_3D(tens_2d, which_leg) result(tens_3d)

  COMPLEX(KIND=DP), INTENT(IN) :: tens_2d(:,:)
  CHARACTER(LEN=*), INTENT(IN) :: which_leg

  COMPLEX(KIND=DP), ALLOCATABLE :: tens_3d(:,:,:)
  INTEGER                       :: dim2D(2)

  dim2D = shape(tens_2d)

  IF(which_leg .EQ. '1') THEN 

     ALLOCATE(tens_3d(1, dim2D(1), dim2D(2)))
     tens_3d(1,:,:) = tens_2d

  ELSEIF(which_leg .EQ. '2') THEN

     ALLOCATE(tens_3d(dim2D(1), 1, dim2D(2)))
     tens_3d(:,1,:) = tens_2d

  ELSEIF(which_leg .EQ. '3') THEN 

     ALLOCATE(tens_3d(dim2D(1), dim2D(2), 1))
     tens_3d(:,:,1) = tens_2d

  ELSE
     CALL invalid_flag("Error: insert_leg has invalid value", which_leg)
  end IF

 END FUNCTION increase_rank_2D_to_3D




 !!! Convert 3D tensor into 4D tensor by inserting an extra singleton leg !!!
 FUNCTION increase_rank_3D_to_4D(tens_3d, which_leg) result(tens_4d)

  COMPLEX(KIND=DP), INTENT(IN) :: tens_3d(:,:,:)
  CHARACTER,        INTENT(IN) :: which_leg*1

  COMPLEX(KIND=DP), ALLOCATABLE :: tens_4d(:,:,:,:)
  INTEGER :: dim3D(3)

  dim3D = shape(tens_3d)

  IF(which_leg .EQ. '1') THEN 

     ALLOCATE(tens_4d(1, dim3D(1), dim3D(2), dim3D(3)))
     tens_4d(1,:,:,:) = tens_3d

  ELSEIF(which_leg .EQ. '2') THEN

     ALLOCATE(tens_4d(dim3D(1), 1, dim3D(2), dim3D(3)))
     tens_4d(:,1,:,:) = tens_3d

  ELSEIF(which_leg .EQ. '3') THEN 

     ALLOCATE(tens_4d(dim3D(1), dim3D(2), 1, dim3D(3)))
     tens_4d(:,:,1,:) = tens_3d

  ELSEIF(which_leg .EQ. '4') THEN 

     ALLOCATE(tens_4d(dim3D(1), dim3D(2), dim3D(3), 1))
     tens_4d(:,:,:,1) = tens_3d

  ELSE
     CALL invalid_flag("Error: insert_leg has invalid value", which_leg)
  end IF

 END FUNCTION increase_rank_3D_to_4D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!! Reduce tensor rank: remove the first encountered singleton leg !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Removes the first encountered singleton leg from a 4D tensor !!!
 FUNCTION reduce_rank_4D(tens_4D) result(tens_3D)

  COMPLEX(KIND=DP), INTENT(IN)  :: tens_4D(:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tens_3D(:,:,:)

  IF(SIZE(tens_4D,1) .EQ. 1) THEN

     ALLOCATE(tens_3D(SIZE(tens_4D,2), SIZE(tens_4D,3), SIZE(tens_4D,4)))
     tens_3D(:,:,:) = tens_4D(1,:,:,:)

  ELSEIF(SIZE(tens_4D,2) .EQ. 1) THEN

     ALLOCATE(tens_3D(SIZE(tens_4D,1), SIZE(tens_4D,3), SIZE(tens_4D,4)))
     tens_3D(:,:,:) = tens_4D(:,1,:,:)

  ELSEIF(SIZE(tens_4D,3) .EQ. 1) THEN

     ALLOCATE(tens_3D(SIZE(tens_4D,1), SIZE(tens_4D,2), SIZE(tens_4D,4)))
     tens_3D(:,:,:) = tens_4D(:,:,1,:)

  ELSEIF(SIZE(tens_4D,4) .EQ. 1) THEN

     ALLOCATE(tens_3D(SIZE(tens_4D,1), SIZE(tens_4D,2), SIZE(tens_4D,3)))
     tens_3D(:,:,:) = tens_4D(:,:,:,1)
 
  ELSE
     WRITE(*,*) "Error: no singleton legs to remove"
     STOP
  end IF 

 END FUNCTION reduce_rank_4D



 !!! Removes the first encountered singleton leg from a 3D tensor !!!
 FUNCTION reduce_rank_3D(tens_3D) result(tens_2D)

  COMPLEX(KIND=DP), INTENT(IN)  :: tens_3D(:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tens_2D(:,:)

  IF(SIZE(tens_3D,1) .EQ. 1) THEN

     ALLOCATE(tens_2D(SIZE(tens_3D,2), SIZE(tens_3D,3)))
     tens_2D(:,:) = tens_3D(1,:,:)

  ELSEIF(SIZE(tens_3D,2) .EQ. 1) THEN

     ALLOCATE(tens_2D(SIZE(tens_3D,1), SIZE(tens_3D,3)))
     tens_2D(:,:) = tens_3D(:,1,:)

  ELSEIF(SIZE(tens_3D,3) .EQ. 1) THEN

     ALLOCATE(tens_2D(SIZE(tens_3D,1), SIZE(tens_3D,2)))
     tens_2D(:,:) = tens_3D(:,:,1)
 
  ELSE
     WRITE(*,*) "Error: no singleton legs to remove"
     STOP
  end IF 

 END FUNCTION reduce_rank_3D



 !!! Removes the first encountered singleton leg from a 2D tensor !!!
 FUNCTION reduce_rank_2D(tens_2D) result(tens_1D)

  COMPLEX(KIND=DP), INTENT(IN)  :: tens_2D(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tens_1D(:)

  INTEGER :: dimA(2)

  !get shape of tens_2D
  dimA = shape(tens_2D)

  IF(dimA(1) .EQ. 1) THEN

     ALLOCATE(tens_1D(dimA(2)))
     tens_1D(:) = tens_2D(1,:)

  ELSEIF(dimA(2) .EQ. 1) THEN

     ALLOCATE(tens_1D(dimA(1)))
     tens_1D(:) = tens_2D(:,1)

  ELSE
     WRITE(*,*) "Error: no singleton legs to remove"
     STOP
  end IF 

 END FUNCTION reduce_rank_2D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE TENS_helper
