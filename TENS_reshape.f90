MODULE TENS_reshape

 USE basic_functions
 USE error_handling
 USE array_utility
 USE TENS_helper

 implicit none

 interface RESHAPE_6D
    module procedure reshape_5D_to_6D
    module procedure reshape_4D_to_6D
 end interface RESHAPE_6D

 interface RESHAPE_5D
    module procedure reshape_6D_to_5D
    module procedure reshape_3D_to_5D
    module procedure reshape_2D_to_5D
 end interface RESHAPE_5D

 interface RESHAPE_4D
    module procedure reshape_6D_to_4D
    module procedure reshape_3D_to_4D
    module procedure reshape_2D_to_4D_double_ind
    module procedure reshape_2D_to_4D_triple_ind
    module procedure reshape_1D_to_4D
 end interface RESHAPE_4D

 interface RESHAPE_3D
    module procedure reshape_5D_to_3D
    module procedure reshape_4D_to_3D
    module procedure reshape_2D_to_3D
    module procedure reshape_1D_to_3D
 end interface RESHAPE_3D

 interface RESHAPE_2D
    module procedure reshape_5D_to_2D
    module procedure reshape_4D_to_2D
    module procedure reshape_3D_to_2D
    module procedure reshape_1D_to_2D
 end interface RESHAPE_2D

 interface RESHAPE_1D
    module procedure reshape_5D_to_1D
    module procedure reshape_4D_to_1D
    module procedure reshape_3D_to_1D
    module procedure reshape_2D_to_1D
 end interface RESHAPE_1D

 private !! hides all items not listed in public statement 
 public RESHAPE_6D, RESHAPE_5D, RESHAPE_4D, RESHAPE_3D, RESHAPE_2D, RESHAPE_1D

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RESHAPE_6D -- reshape into 6D tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Reshape 5D tensor into 6D tensor !!!
 FUNCTION reshape_5D_to_6D(tensA, legs, dimLEG) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:,:,:,:) 
  CHARACTER(LEN=*),  INTENT(IN)  :: legs
  INTEGER,           INTENT(IN)  :: dimLEG(2)
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:,:,:,:,:)   

  !! Dims & indices
  INTEGER :: dimA(5), dimO(6)
  INTEGER :: i, j

  !! Get dims
  dimA = shape(tensA)

  !! (0) CONSISTENCY CHECKS
  SELECT CASE(legs)
  CASE('12,3,4,5,6')  
       CALL check_sizes_equal(dimA(1), product(dimLEG), "reshape_5D_to_6D -- tensA && tensO must have compatible dims ")      
  CASE DEFAULT
       CALL invalid_flag("reshape_5D_to_6D -- invalid legs ", legs)
  end SELECT

  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('12,3,4,5,6')  
       dimO = (/ dimLEG(1), dimLEG(2), dimA(2), dimA(3), dimA(4), dimA(5) /) 
  CASE DEFAULT
       CALL invalid_flag("reshape_5D_to_6D -- invalid legs ", legs)
  end SELECT

  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  Do i=1,dimLEG(1)
    Do j=1,dimLEG(2)  
       tensO(i,j,:,:,:,:) = tensA(ICOM(i,j,dimLEG(2)),:,:,:,:)   
    end Do
  end Do

 END FUNCTION reshape_5D_to_6D




 !!! Reshape 4D tensor into 6D tensor !!!
 FUNCTION reshape_4D_to_6D(tensA, legs, dimLEG) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN) :: tensA(:,:,:,:)
  CHARACTER(LEN=*),  INTENT(IN) :: legs
  INTEGER,           INTENT(IN) :: dimLEG(3)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:,:,:)

  !! Dims & indices 
  INTEGER :: dimA(4), dimO(6)
  INTEGER :: o1, o2, o3

  !! Get dims
  dimA = shape(tensA)

  !! (0) CONSISTENCY CHECKS
  SELECT CASE(legs)
  CASE('1,2,345,6') 
        CALL check_sizes_equal(dimA(3), product(dimLEG), "reshape_4D_to_6D -- tensA && tensO must have compatible dims ")
  CASE('1,2,5,346')  
        CALL check_sizes_equal(dimA(4), product(dimLEG), "reshape_4D_to_6D -- tensA && tensO must have compatible dims ")
  CASE('1,2,356,4')  
        CALL check_sizes_equal(dimA(3), product(dimLEG), "reshape_4D_to_6D -- tensA && tensO must have compatible dims ")
  CASE('1,2,3,456')  
        CALL check_sizes_equal(dimA(4), product(dimLEG), "reshape_4D_to_6D -- tensA && tensO must have compatible dims ")
  CASE DEFAULT
        CALL invalid_flag("reshape_4D_to_6D -- invalid legs ", legs)
  end SELECT


  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('1,2,345,6') 
        dimO = (/ dimA(1),   dimA(2),   dimLEG(1), dimLEG(2), dimLEG(3), dimA(4) /)
  CASE('1,2,5,346') 
        dimO = (/ dimA(1),   dimA(2),   dimLEG(1), dimLEG(2), dimA(3),   dimLEG(3) /)
  CASE('1,2,356,4') 
        dimO = (/ dimA(1),   dimA(2),   dimLEG(1), dimA(4),   dimLEG(2), dimLEG(3) /)
  CASE('1,2,3,456')  
        dimO = (/ dimA(1),   dimA(2),   dimA(3),   dimLEG(1), dimLEG(2), dimLEG(3) /)
  CASE DEFAULT
        CALL invalid_flag("reshape_4D_to_6D -- invalid legs ", legs)
  end SELECT


  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  Do o1=1,dimLEG(1)
    Do o2=1,dimLEG(2)
      Do o3=1,dimLEG(3)
         CALL split_legs(tensO, tensA, legs, ICOM(o1,o2,dimLEG(2),o3,dimLEG(3)), o1, o2, o3)
      end Do
    end Do 
  end Do

 END FUNCTION reshape_4D_to_6D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RESHAPE_5D -- reshape into 5D tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Reshape 6D tensor into 5D tensor !!!
 FUNCTION reshape_6D_to_5D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:,:,:,:,:) 
  CHARACTER(LEN=*),  INTENT(IN)  :: legs
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:,:,:,:)   

  !! Dims & indices
  INTEGER :: dimA(6), dimO(5), dimCOM(2)
  INTEGER :: i, j

  !! Get dims
  dimA = shape(tensA)

  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('12,3,4,5,6')  
     dimCOM = (/ dimA(1), dimA(2) /)
       dimO = (/ product(dimCOM), dimA(3), dimA(4), dimA(5), dimA(6) /) 
  CASE DEFAULT
       CALL invalid_flag("reshape_6D_to_5D -- invalid legs ", legs)
  end SELECT

  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  Do i=1,dimCOM(1)
    Do j=1,dimCOM(2)  
       tensO(ICOM(i,j,dimCOM(2)),:,:,:,:) = tensA(i,j,:,:,:,:)   
    end Do
  end Do

 END FUNCTION reshape_6D_to_5D




 !!! Reshape 3D tensor into 5D tensor !!!
 FUNCTION reshape_3D_to_5D(tensA, legs, dimLEG) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN) :: tensA(:,:,:)
  CHARACTER(LEN=*),  INTENT(IN) :: legs
  INTEGER,           INTENT(IN) :: dimLEG(3)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:,:)

  !! Dims & indices 
  INTEGER :: dimA(3), dimO(5)
  INTEGER :: o1, o2, o3

  !! Get dims
  dimA = shape(tensA)

  !! (0) CONSISTENCY CHECKS
  SELECT CASE(legs)
  CASE('1,234,5') 
        CALL check_sizes_equal(dimA(2), product(dimLEG), "reshape_3D_to_5D -- tensA && tensO must have compatible dims ")
  CASE('1,4,235')  
        CALL check_sizes_equal(dimA(3), product(dimLEG), "reshape_3D_to_5D -- tensA && tensO must have compatible dims ")
  CASE('1,245,3')  
        CALL check_sizes_equal(dimA(2), product(dimLEG), "reshape_3D_to_5D -- tensA && tensO must have compatible dims ")
  CASE('1,2,345')  
        CALL check_sizes_equal(dimA(3), product(dimLEG), "reshape_3D_to_5D -- tensA && tensO must have compatible dims ")
  CASE('123,4,5') 
        CALL check_sizes_equal(dimA(1), product(dimLEG), "reshape_3D_to_5D -- tensA && tensO must have compatible dims ")
  CASE('145,2,3') 
        CALL check_sizes_equal(dimA(1), product(dimLEG), "reshape_3D_to_5D -- tensA && tensO must have compatible dims ")
  CASE DEFAULT
        CALL invalid_flag("reshape_3D_to_5D -- invalid legs ", legs)
  end SELECT


  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('1,234,5') 
        dimO = (/ dimA(1),   dimLEG(1), dimLEG(2), dimLEG(3), dimA(3) /)
  CASE('1,4,235') 
        dimO = (/ dimA(1),   dimLEG(1), dimLEG(2), dimA(2),   dimLEG(3) /)
  CASE('1,245,3') 
        dimO = (/ dimA(1),   dimLEG(1), dimA(3),   dimLEG(2), dimLEG(3) /)
  CASE('1,2,345')  
        dimO = (/ dimA(1),   dimA(2),   dimLEG(1), dimLEG(2), dimLEG(3) /)
  CASE('123,4,5')  
        dimO = (/ dimLEG(1), dimLEG(2), dimLEG(3), dimA(2),   dimA(3) /)
  CASE('145,2,3') 
        dimO = (/ dimLEG(1), dimA(2),   dimA(3),   dimLEG(2), dimLEG(3) /)
  end SELECT


  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  Do o1=1,dimLEG(1)
    Do o2=1,dimLEG(2)
      Do o3=1,dimLEG(3)
         CALL split_legs(tensO, tensA, legs, ICOM(o1,o2,dimLEG(2),o3,dimLEG(3)), o1, o2, o3)
      end Do
    end Do 
  end Do

 END FUNCTION reshape_3D_to_5D









 !!! Reshape 2D tensor into 5D tensor !!!
 FUNCTION reshape_2D_to_5D(tensA, legs, dimLEG) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:)
  CHARACTER(LEN=*),  INTENT(IN)  :: legs
  INTEGER,           INTENT(IN)  :: dimLEG(4)
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:,:,:,:)

  !! Dims & indices
  INTEGER :: dimA(2), dimO(5)
  INTEGER :: o1, o2, o3, o4

  !! Get dims
  dimA = shape(tensA)

  !! (0) CONSISTENCY CHECKS
  SELECT CASE(legs)
  CASE('1234,5') 
       CALL check_sizes_equal(dimA(1), product(dimLEG), "reshape_2D_to_5D -- tensA && tensO must have compatible dims ")
  CASE('4,1235') 
       CALL check_sizes_equal(dimA(2), product(dimLEG), "reshape_2D_to_5D -- tensA && tensO must have compatible dims ")
  CASE('1245,3') 
       CALL check_sizes_equal(dimA(1), product(dimLEG), "reshape_2D_to_5D -- tensA && tensO must have compatible dims ")
  CASE('2,1345') 
       CALL check_sizes_equal(dimA(2), product(dimLEG), "reshape_2D_to_5D -- tensA && tensO must have compatible dims ")
  CASE('1,2345') 
       CALL check_sizes_equal(dimA(2), product(dimLEG), "reshape_2D_to_5D -- tensA && tensO must have compatible dims ")
  CASE DEFAULT
       CALL invalid_flag("reshape_2D_to_5D -- invalid legs ", legs)
  end SELECT


  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('1234,5')
       dimO = (/ dimLEG(1), dimLEG(2), dimLEG(3), dimLEG(4), dimA(2)   /)
  CASE('4,1235') 
       dimO = (/ dimLEG(1), dimLEG(2), dimLEG(3), dimA(1),   dimLEG(4) /)
  CASE('1245,3')
       dimO = (/ dimLEG(1), dimLEG(2), dimA(2),   dimLEG(3), dimLEG(4) /)
  CASE('2,1345') 
       dimO = (/ dimLEG(1), dimA(1),   dimLEG(2), dimLEG(3), dimLEG(4) /)
  CASE('1,2345') 
       dimO = (/ dimA(1),   dimLEG(1), dimLEG(2), dimLEG(3), dimLEG(4) /) 
  CASE DEFAULT
       CALL invalid_flag("reshape_2D_to_5D -- invalid legs ", legs)
  end SELECT


  !! (2) RESHAPE LOOP
  CALL allocateTens(tensO, dimO)

  Do o1=1,dimLEG(1)
    Do o2=1,dimLEG(2)
      Do o3=1,dimLEG(3)
        Do o4=1,dimLEG(4)
           CALL split_legs(tensO, tensA, legs, ICOM(o1,o2,dimLEG(2),o3,dimLEG(3),o4,dimLEG(4)), o1, o2, o3, o4)
        end Do
      end Do
    end Do 
  end Do

 END FUNCTION reshape_2D_to_5D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RESHAPE_4D -- reshape into 4D tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !!! Reshape 4D tensor into 6D tensor !!!
 FUNCTION reshape_6D_to_4D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:,:,:,:,:)
  CHARACTER(LEN=*),  INTENT(IN)  :: legs
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:,:,:)

  !! Dims & indices 
  INTEGER :: dimA(6), dimO(4), dimCOM(3)
  INTEGER :: a1, a2, a3

  !! Get dims
  dimA = shape(tensA)

  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('1,2,345,6')  !! Combine [axis-3,4,5 of tensA] into [axis-3 of tensO]

     dimCOM = (/ dimA(3), dimA(4), dimA(5) /)
       dimO = (/ dimA(1), dimA(2), product(dimCOM), dimA(6) /)

  CASE('1,2,5,346')  !! Combine [axis-3,4,6 of tensA] into [axis-4 of tensO]
     
     dimCOM = (/ dimA(3), dimA(4), dimA(6) /)
       dimO = (/ dimA(1), dimA(2), dimA(5), product(dimCOM) /)

  CASE('1,2,356,4')  !! Combine [axis-3,5,6 of tensA] into [axis-3 of tensO]

     dimCOM = (/ dimA(3), dimA(5), dimA(6) /)
       dimO = (/ dimA(1), dimA(2), product(dimCOM), dimA(4) /)

  CASE('1,2,3,456')  !! Combine [axis-4,5,6 of tensA] into [axis-4 of tensO]
     
     dimCOM = (/ dimA(4), dimA(5), dimA(6) /)
       dimO = (/ dimA(1), dimA(2), dimA(3), product(dimCOM) /)

  CASE DEFAULT
     CALL invalid_flag("reshape_5D_to_3D -- invalid legs ", legs)
  end SELECT


  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  Do a1=1,dimCOM(1)
    Do a2=1,dimCOM(2)
      Do a3=1,dimCOM(3)
         CALL fuse_legs(tensO, tensA, legs, ICOM(a1,a2,dimCOM(2),a3,dimCOM(3)), a1, a2, a3)
      end Do
    end Do 
  end Do

 END FUNCTION reshape_6D_to_4D



 !!! Reshape 3D tensor into 4D tensor !!!
 FUNCTION reshape_3D_to_4D(tensA, legs, dimLEG) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:)    
  CHARACTER(LEN=*), INTENT(IN)  :: legs
  INTEGER,          INTENT(IN)  :: dimLEG(2)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:)  

  !! Dims & indices
  INTEGER :: dimA(3), dimO(4)
  INTEGER :: i, j

  !! Get dims
  dimA = shape(tensA)

  !! (0) CONSISTENCY CHECKS
  SELECT CASE(legs)
  CASE('12,3,4')  
       CALL check_sizes_equal(dimA(1), product(dimLEG), "reshape_4D_to_3D -- tensA && tensO must have compatible dims ")
  CASE('1,2,34')
       CALL check_sizes_equal(dimA(3), product(dimLEG), "reshape_4D_to_3D -- tensA && tensO must have compatible dims ")      
  CASE('2,3,14')
       CALL check_sizes_equal(dimA(3), product(dimLEG), "reshape_4D_to_3D -- tensA && tensO must have compatible dims ")       
  CASE('2,13,4')
       CALL check_sizes_equal(dimA(2), product(dimLEG), "reshape_4D_to_3D -- tensA && tensO must have compatible dims ")       
  CASE('1,3,24')
       CALL check_sizes_equal(dimA(3), product(dimLEG), "reshape_4D_to_3D -- tensA && tensO must have compatible dims ")       
  CASE('1,23,4')
       CALL check_sizes_equal(dimA(2), product(dimLEG), "reshape_4D_to_3D -- tensA && tensO must have compatible dims ")       
  CASE DEFAULT
       CALL invalid_flag("reshape_3D_to_4D -- invalid legs ", legs)
  end SELECT


  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('12,3,4')  
       dimO = (/ dimLEG(1), dimLEG(2), dimA(2),   dimA(3) /) 
  CASE('1,2,34')
       dimO = (/ dimA(1),   dimA(2),   dimLEG(1), dimLEG(2) /) 
  CASE('2,3,14')
       dimO = (/ dimLEG(1), dimA(1),   dimA(2),   dimLEG(2) /) 
  CASE('2,13,4')
       dimO = (/ dimLEG(1), dimA(1),   dimLEG(2), dimA(3) /) 
  CASE('1,3,24')
       dimO = (/ dimA(1),   dimLEG(1), dimA(2),   dimLEG(2) /) 
  CASE('1,23,4')
       dimO = (/ dimA(1),   dimLEG(1), dimLEG(2), dimA(3) /) 
  CASE DEFAULT
      CALL invalid_flag("reshape_3D_to_4D -- invalid legs ", legs)
  end SELECT


  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  Do i=1,dimLEG(1)
    Do j=1,dimLEG(2)
       CALL split_legs(tensO, tensA, legs, ICOM(i,j,dimLEG(2)), i, j)
    end Do
  end Do

 END FUNCTION reshape_3D_to_4D







 !!! Reshape 2D tensor into 4D tensor
 FUNCTION reshape_2D_to_4D_double_ind(tensA, legs, dimLEG_AA, dimLEG_BB) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:) 
  CHARACTER(LEN=*),  INTENT(IN)  :: legs
  INTEGER,           INTENT(IN)  :: dimLEG_AA(2), dimLEG_BB(2) 
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:,:,:)

  !! Dims & indices
  INTEGER :: dimA(2), dimO(4)
  INTEGER :: ia, ja, ib, jb

  !! Get dims
  dimA = shape(tensA)

  !! (0) CONSISTENCY CHECKS
  SELECT CASE(legs)
  CASE('13,24', '23,14', '12,34')
       CALL check_sizes_equal(dimA(1), product(dimLEG_AA), "reshape_2D_to_4D_double_ind -- tensA && tensO must have compatible dims ")
       CALL check_sizes_equal(dimA(2), product(dimLEG_BB), "reshape_2D_to_4D_double_ind -- tensA && tensO must have compatible dims ")
  CASE DEFAULT
       CALL invalid_flag("reshape_2D_to_4D_double_ind -- invalid legs ", legs) 
  end SELECT

  
  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('13,24')
       dimO = (/ dimLEG_AA(1), dimLEG_BB(1), dimLEG_AA(2), dimLEG_BB(2) /)
  CASE('23,14')
       dimO = (/ dimLEG_BB(1), dimLEG_AA(1), dimLEG_AA(2), dimLEG_BB(2) /)
  CASE('12,34')
       dimO = (/ dimLEG_AA(1), dimLEG_AA(2), dimLEG_BB(1), dimLEG_BB(2) /)     
  CASE DEFAULT
       CALL invalid_flag("reshape_2D_to_4D_double_ind -- invalid legs ", legs) 
  end SELECT


  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  DO ia=1,dimLEG_AA(1)
    DO ja=1,dimLEG_AA(2)
      DO ib=1,dimLEG_BB(1)
        DO jb=1,dimLEG_BB(2)
           CALL split_legs(tensO, tensA, legs, ICOM(ia,ja,dimLEG_AA(2)), ia, ja, ICOM(ib,jb,dimLEG_BB(2)), ib, jb)
        end DO
      end DO
    end DO
  end DO

 END FUNCTION reshape_2D_to_4D_double_ind






 !!! Reshape 2D tensor into 4D tensor
 FUNCTION reshape_2D_to_4D_triple_ind(tensA, legs, dimLEG) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:) 
  CHARACTER(LEN=*),  INTENT(IN)  :: legs
  INTEGER,           INTENT(IN)  :: dimLEG(3) 
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:,:,:)

  !! Dims & indices
  INTEGER :: dimA(2), dimO(4)
  INTEGER :: i,j,k

  !! Get dims
  dimA = shape(tensA)

  !! (0) CONSISTENCY CHECKS
  SELECT CASE(legs)
  CASE('1,234')
       CALL check_sizes_equal(dimA(2), product(dimLEG), "reshape_2D_to_4D_triple_ind -- tensA && tensO must have compatible dims ")
  CASE DEFAULT
       CALL invalid_flag("reshape_2D_to_4D_triple_ind -- invalid legs ", legs) 
  end SELECT

  
  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('1,234')
       dimO = (/ dimA(1), dimLEG(1), dimLEG(2), dimLEG(3) /)    
  CASE DEFAULT
       CALL invalid_flag("reshape_2D_to_4D_triple_ind -- invalid legs ", legs) 
  end SELECT


  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  DO i=1,dimLEG(1)
    DO j=1,dimLEG(2)
      DO k=1,dimLEG(3)
         CALL split_legs(tensO, tensA, legs, ICOM(i, j, dimLEG(2), k, dimLEG(3)), i, j, k)
      end DO
    end DO
  end DO

 END FUNCTION reshape_2D_to_4D_triple_ind





 !!! Reshape vector into 4D tensor !!!
 FUNCTION reshape_1D_to_4D(tensA, dimO) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:)
  INTEGER,           INTENT(IN)  :: dimO(4)
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:,:,:) 

  !! Indices && dims (of given dims)
  INTEGER :: o1, o2, o3, o4

  !! Check if dims match
  CALL check_sizes_equal(product(dimO), SIZE(tensA), "reshape_1D_to_4D: sizes of tensO and tensA must match ")

  !! Reshape 1D into 4D
  ALLOCATE(tensO(dimO(1), dimO(2), dimO(3), dimO(4)))

  DO o1=1,dimO(1)
    DO o2=1,dimO(2)
      DO o3=1,dimO(3)
        DO o4=1,dimO(4)
           tensO(o1, o2, o3, o4) = tensA(ICOM(o1, o2, dimO(2), o3, dimO(3), o4, dimO(4)))
        end DO
      end DO
    end DO
  end DO

 END FUNCTION reshape_1D_to_4D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RESHAPE_3D -- reshape into 3D tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Reshape 5D tensor into 3D tensor !!!
 FUNCTION reshape_5D_to_3D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)   :: tensA(:,:,:,:,:)
  CHARACTER(LEN=*),  INTENT(IN)   :: legs
  COMPLEX(KIND=DP),  ALLOCATABLE  :: tensO(:,:,:)

  !! Dims & indices 
  INTEGER :: dimA(5), dimO(3), dimCOM(3)
  INTEGER :: a1, a2, a3

  !! Get dims
  dimA = shape(tensA)

  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('1,234,5')  !! Combine [axis-2,3,4 of tensA] into [axis-2 of tensO]

     dimCOM = (/ dimA(2), dimA(3), dimA(4) /)
       dimO = (/ dimA(1), product(dimCOM), dimA(5) /)

  CASE('1,4,235')  !! Combine [axis-2,3,5 of tensA] into [axis-3 of tensO]
     
     dimCOM = (/ dimA(2), dimA(3), dimA(5) /)
       dimO = (/ dimA(1), dimA(4), product(dimCOM) /)

  CASE('1,245,3')  !! Combine [axis-2,4,5 of tensA] into [axis-2 of tensO]

     dimCOM = (/ dimA(2), dimA(4), dimA(5) /)
       dimO = (/ dimA(1), product(dimCOM), dimA(3) /)

  CASE('1,2,345')  !! Combine [axis-3,4,5 of tensA] into [axis-3 of tensO]
     
     dimCOM = (/ dimA(3), dimA(4), dimA(5) /)
       dimO = (/ dimA(1), dimA(2), product(dimCOM) /)

  CASE('123,4,5')  !! Combine [axis-1,2,3 of tensA] into [axis-1 of tensO]

     dimCOM = (/ dimA(1), dimA(2), dimA(3) /)
       dimO = (/ product(dimCOM), dimA(4), dimA(5) /)

  CASE('145,2,3')  !! Combine [axis-1,4,5 of tensA] into [axis-1 of tensO]
     
     dimCOM = (/ dimA(1), dimA(4), dimA(5) /)
       dimO = (/ product(dimCOM), dimA(2), dimA(3) /)

  CASE DEFAULT
     CALL invalid_flag("reshape_5D_to_3D -- invalid legs ", legs)
  end SELECT


  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  Do a1=1,dimCOM(1)
    Do a2=1,dimCOM(2)
      Do a3=1,dimCOM(3)
         CALL fuse_legs(tensO, tensA, legs, ICOM(a1,a2,dimCOM(2),a3,dimCOM(3)), a1, a2, a3)
      end Do
    end Do 
  end Do

 END FUNCTION reshape_5D_to_3D





 !!! Reshape 4D tensor into 3D tensor !!!
 FUNCTION reshape_4D_to_3D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:,:,:) 
  CHARACTER(LEN=*),  INTENT(IN)  :: legs
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:,:)   

  !! Dims & indices
  INTEGER :: dimA(4), dimO(3), dimCOM(2)
  INTEGER :: i, j

  !! Get dims
  dimA = shape(tensA)

  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('12,3,4')  

    dimCOM = (/ dimA(1), dimA(2) /)
      dimO = (/ product(dimCOM), dimA(3), dimA(4) /) 

  CASE('1,2,34')

    dimCOM = (/ dimA(3), dimA(4) /)
      dimO = (/ dimA(1), dimA(2), product(dimCOM) /) 

  CASE('2,3,14')

    dimCOM = (/ dimA(1), dimA(4) /)
      dimO = (/ dimA(2), dimA(3), product(dimCOM) /) 

  CASE('2,13,4')

    dimCOM = (/ dimA(1), dimA(3) /)
      dimO = (/ dimA(2), product(dimCOM), dimA(4) /) 

  CASE('1,3,24')

    dimCOM = (/ dimA(2), dimA(4) /)
      dimO = (/ dimA(1), dimA(3), product(dimCOM) /) 

  CASE('1,23,4')

    dimCOM = (/ dimA(2), dimA(3) /)
      dimO = (/ dimA(1), product(dimCOM), dimA(4) /) 

  CASE DEFAULT
    CALL invalid_flag("reshape_4D_to_3D -- invalid legs ", legs)
  end SELECT


  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  Do i=1,dimCOM(1)
    Do j=1,dimCOM(2)     
       CALL fuse_legs(tensO, tensA, legs, ICOM(i,j,dimCOM(2)), i, j)
    end Do
  end Do

 END FUNCTION reshape_4D_to_3D






 !!! Reshape 2D tensor into 3D tensor !!!
 FUNCTION reshape_2D_to_3D(tensA, legs, dimLEG) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN) :: tensA(:,:)
  CHARACTER(LEN=*),  INTENT(IN) :: legs
  INTEGER,           INTENT(IN) :: dimLEG(2)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:)

  !! Dims & indices
  INTEGER :: dimA(2), dimO(3)
  INTEGER :: i, j

  !! Get dims
  dimA = shape(tensA)

  !! (0) CONSISTENCY CHECKS
  SELECT CASE(legs)
  CASE('12,3')
       CALL check_sizes_equal(dimA(1), product(dimLEG), "reshape_2D_to_3D -- tensA && tensO must have compatible dims ")
  CASE('2,13')
       CALL check_sizes_equal(dimA(2), product(dimLEG), "reshape_2D_to_3D -- tensA && tensO must have compatible dims ")
  CASE('1,23')
       CALL check_sizes_equal(dimA(2), product(dimLEG), "reshape_2D_to_3D -- tensA && tensO must have compatible dims ")
  CASE DEFAULT
       CALL invalid_flag("reshape_2D_to_3D -- invalid legs ", legs) 
  end SELECT


  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('12,3')
       dimO = (/ dimLEG(1), dimLEG(2), dimA(2)   /)
  CASE('2,13')
       dimO = (/ dimLEG(1), dimA(1),   dimLEG(2) /)       
  CASE('1,23')
       dimO = (/ dimA(1),   dimLEG(1), dimLEG(2) /)       
  CASE DEFAULT
       CALL invalid_flag("reshape_2D_to_3D -- invalid legs ", legs) 
  end SELECT


  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  DO i=1,dimLEG(1)
    DO j=1,dimLEG(2)
       CALL split_legs(tensO, tensA, legs, ICOM(i,j,dimLEG(2)), i, j)
    end DO
  end DO

 END FUNCTION reshape_2D_to_3D






 !!! Reshape vector into 3D tensor !!!
 FUNCTION reshape_1D_to_3D(tensA, dimO) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:)
  INTEGER,           INTENT(IN)  :: dimO(3)
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:,:) 

  !! Indices && dims (of given dims)
  INTEGER :: o1, o2, o3

  !! Check if dims match
  CALL check_sizes_equal(product(dimO), SIZE(tensA), "reshape_1D_to_3D: sizes of tensO and tensA must match ")

  !! Reshape 1D into 3D
  ALLOCATE(tensO(dimO(1), dimO(2), dimO(3)))

  DO o1=1,dimO(1)
    DO o2=1,dimO(2)
      DO o3=1,dimO(3)
         tensO(o1, o2, o3) = tensA(ICOM(o1, o2, dimO(2), o3, dimO(3)))
      end DO
    end DO
  end DO

 END FUNCTION reshape_1D_to_3D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RESHAPE_2D -- reshape into 2D tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Reshape 5D tensor into 2D tensor !!!
 FUNCTION reshape_5D_to_2D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN) :: tensA(:,:,:,:,:)
  CHARACTER(LEN=*),  INTENT(IN) :: legs
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:)

  !! Dims & indices 
  INTEGER :: dimA(5), dimCOM(4), dimO(2)
  INTEGER :: a1, a2, a3, a4

  !! Get dims
  dimA = shape(tensA)

  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('1234,5')  !Combine [axis-1,2,3,4 of tensA] into [axis-1 of tensO]

       dimCOM = (/ dimA(1), dimA(2), dimA(3), dimA(4) /)
         dimO = (/ product(dimCOM), dimA(5) /)     

  CASE('4,1235')  !Combine [axis-1,2,3,5 of tensA] into [axis-2 of tensO]

       dimCOM = (/ dimA(1), dimA(2), dimA(3), dimA(5) /)
         dimO = (/ dimA(4), product(dimCOM) /)

  CASE('1245,3')  !Combine [axis-1,2,4,5 of tensA] into [axis-1 of tensO]

       dimCOM = (/ dimA(1), dimA(2), dimA(4), dimA(5) /)
         dimO = (/ product(dimCOM), dimA(3) /)

  CASE('2,1345')  !Combine [axis-1,3,4,5 of tensA] into [axis-2 of tensO]

       dimCOM = (/ dimA(1), dimA(3), dimA(4), dimA(5) /)
         dimO = (/ dimA(2), product(dimCOM) /)

  CASE('1,2345')  !Combine [axis-2,3,4,5 of tensA] into [axis-2 of tensO]

       dimCOM = (/ dimA(2), dimA(3), dimA(4), dimA(5) /)
         dimO = (/ dimA(1), product(dimCOM) /)

  CASE DEFAULT
       CALL invalid_flag("reshape_5D_to_2D -- invalid legs ", legs)
  end SELECT


  !! (2) RESHAPE LOOP
  CALL allocateTens(tensO, dimO)

  Do a1=1,dimCOM(1)
    Do a2=1,dimCOM(2)
      Do a3=1,dimCOM(3)
        Do a4=1,dimCOM(4)
           CALL fuse_legs(tensO, tensA, legs, ICOM(a1,a2,dimCOM(2),a3,dimCOM(3),a4,dimCOM(4)), a1, a2, a3, a4)
        end Do
      end Do
    end Do 
  end Do

 END FUNCTION reshape_5D_to_2D 






 !!! Reshape 4D tensor into 2D tensor !!!
 FUNCTION reshape_4D_to_2D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:,:,:) 
  CHARACTER(LEN=*),  INTENT(IN)  :: legs
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:) 

  !! Dims & indices
  INTEGER :: dimA(4), dimO(2)
  INTEGER :: dimCOM_AA(2), dimCOM_BB(2)
  INTEGER :: ia, ja, ib, jb

  INTEGER :: dimCOM(3)
  INTEGER :: i,j,k

  !! Get dims
  dimA = shape(tensA)

  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('13,24')
       dimCOM_AA = (/ dimA(1), dimA(3) /)
       dimCOM_BB = (/ dimA(2), dimA(4) /)
            dimO = (/ product(dimCOM_AA), product(dimCOM_BB) /)
  CASE('23,14')
       dimCOM_AA = (/ dimA(2), dimA(3) /)
       dimCOM_BB = (/ dimA(1), dimA(4) /)
            dimO = (/ product(dimCOM_AA), product(dimCOM_BB) /)
  CASE('12,34')
       dimCOM_AA = (/ dimA(1), dimA(2) /)
       dimCOM_BB = (/ dimA(3), dimA(4) /)
            dimO = (/ product(dimCOM_AA), product(dimCOM_BB) /)   
  CASE('1,234')
          dimCOM = (/ dimA(2), dimA(3), dimA(4) /)
            dimO = (/ dimA(1), product(dimCOM) /)    
  CASE DEFAULT
       CALL invalid_flag("reshape_4D_to_2D -- invalid legs ", legs) 
  end SELECT


  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  SELECT CASE(legs)
  CASE('13,24', '12,34', '23,14')

      DO ia=1,dimCOM_AA(1)
        DO ja=1,dimCOM_AA(2)
          DO ib=1,dimCOM_BB(1)
            DO jb=1,dimCOM_BB(2)
               CALL fuse_legs(tensO, tensA, legs, ICOM(ia,ja,dimCOM_AA(2)), ia, ja, ICOM(ib,jb,dimCOM_BB(2)), ib, jb)
            end DO
          end DO
        end DO
      end DO

  CASE('1,234')

      DO i=1,dimCOM(1)
        DO j=1,dimCOM(2)
          DO k=1,dimCOM(3)
             CALL fuse_legs(tensO, tensA, legs, ICOM(i, j, dimCOM(2), k, dimCOM(3)), i, j, k)
          end DO
        end DO
      end DO

  CASE DEFAULT
       CALL invalid_flag("reshape_4D_to_2D -- invalid legs ", legs) 
  end SELECT

 END FUNCTION reshape_4D_to_2D





 !!! Reshape 3D tensor into 2D tensor !!!
 FUNCTION reshape_3D_to_2D(tensA, legs) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:,:) 
  CHARACTER(LEN=*),  INTENT(IN)  :: legs
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:)   

  !! Dims & indices
  INTEGER :: dimA(3), dimCOM(2), dimO(2)
  INTEGER :: i, j

  !! Get dims
  dimA = shape(tensA)

  !! (1) SET DIMS
  SELECT CASE(legs)
  CASE('12,3')
       dimCOM = (/ dimA(1),          dimA(2) /)
         dimO = (/ product(dimCOM),  dimA(3) /)
  CASE('2,13')
       dimCOM = (/ dimA(1),  dimA(3) /)
         dimO = (/ dimA(2),  product(dimCOM) /)       
  CASE('1,23')
       dimCOM = (/ dimA(2),  dimA(3) /)
         dimO = (/ dimA(1),  product(dimCOM) /)       
  CASE DEFAULT
       CALL invalid_flag("reshape_3D_to_2D -- invalid legs ", legs) 
  end SELECT

  !! (2) RESHAPING LOOP
  CALL allocateTens(tensO, dimO)

  DO i=1,dimCOM(1)
    DO j=1,dimCOM(2)
       CALL fuse_legs(tensO, tensA, legs, ICOM(i,j,dimCOM(2)), i, j)
    end DO
  end DO
 
 END FUNCTION reshape_3D_to_2D





 !!! Reshape vector into 2D tensor !!!
 FUNCTION reshape_1D_to_2D(tensA, dimO) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:)
  INTEGER,           INTENT(IN)  :: dimO(2)
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:) 

  !! Indices && dims
  INTEGER :: o1, o2

  !Check if dims match
  CALL check_sizes_equal(product(dimO), SIZE(tensA), "reshape_1D_to_2D: sizes of tensO and tensA must match ")

  !! Reshape 1D into 2D
  ALLOCATE(tensO(dimO(1), dimO(2)))

  DO o1=1,dimO(1)
    DO o2=1,dimO(2)
       tensO(o1, o2) = tensA(ICOM(o1, o2, dimO(2)))
    end DO
  end DO

 END FUNCTION reshape_1D_to_2D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RESHAPE_1D -- reshape into 1D tensor !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Reshape 4D tensor into vector !!!
 FUNCTION reshape_5D_to_1D(tensA) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:,:,:,:)
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:) 

  !! Dims && indices
  INTEGER :: dimA(5)

  !! Get dims
  dimA = shape(tensA)

  !! Reshape tensor
  ALLOCATE(tensO(product(dimA)))

  tensO = RESHAPE(tensA, (/ PRODUCT(dimA) /))

 END FUNCTION reshape_5D_to_1D




 !!! Reshape 4D tensor into vector !!!
 FUNCTION reshape_4D_to_1D(tensA) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:,:,:)
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:) 

  !! Dims && indices
  INTEGER :: dimA(4), a1, a2, a3, a4

  !! Get dims
  dimA = shape(tensA)

  !! Reshape tensor
  ALLOCATE(tensO(product(dimA)))

  DO a1=1,dimA(1)
    DO a2=1,dimA(2)
      DO a3=1,dimA(3)
        DO a4=1,dimA(4)
           tensO(ICOM(a1, a2, dimA(2), a3, dimA(3), a4, dimA(4))) = tensA(a1, a2, a3, a4)
        end DO
      end DO
    end DO
  end DO

 END FUNCTION reshape_4D_to_1D




 !!! Reshape 3D tensor into vector !!!
 FUNCTION reshape_3D_to_1D(tensA) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:,:)
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:) 

  !! Dims && indices
  INTEGER :: dimA(3), a1, a2, a3

  !! Get dims
  dimA = shape(tensA)

  !! Reshape tensor
  ALLOCATE(tensO(product(dimA)))

  DO a1=1,dimA(1)
    DO a2=1,dimA(2)
      DO a3=1,dimA(3)
         tensO(ICOM(a1, a2, dimA(2), a3, dimA(3))) = tensA(a1, a2, a3)
      end DO
    end DO
  end DO

 END FUNCTION reshape_3D_to_1D





 !!! Reshape 2D tensor into vector !!!
 FUNCTION reshape_2D_to_1D(tensA) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:)
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:) 

  !! Dims && indices
  INTEGER :: dimA(2), a1, a2

  !! Get dims
  dimA = shape(tensA)

  !! Reshape 1D into 2D
  ALLOCATE(tensO(product(dimA)))

  DO a1=1,dimA(1)
    DO a2=1,dimA(2)
       tensO(ICOM(a1, a2, dimA(2))) = tensA(a1, a2)
    end DO
  end DO

 END FUNCTION reshape_2D_to_1D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE TENS_reshape
