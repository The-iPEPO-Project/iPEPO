MODULE TENS_mult

 USE basic_functions
 USE error_handling
 USE array_utility
 USE TENS_helper
 USE TENS_reshape
 USE TENS_transpose

 implicit none

 interface TENSMUL

    module procedure contract_4D_4D_result_4D
    module procedure contract_4D_3D_result_3D
    module procedure contract_4D_2D_result_4D

    module procedure contract_4D_3D_result_4D
    module procedure contract_3D_3D_result_4D
    module procedure contract_3D_3D_result_2D
    module procedure contract_3D_2D_result_3D
    module procedure contract_3D_2D_result_2D

    module procedure contract_2D_2D_result_2D_transp
    module procedure contract_2D_2D_result_2D
    module procedure contract_2D_1D_result_1D
    module procedure contract_1D_2D_result_1D

    module procedure contract_5D_2D_result_5D
    module procedure contract_5D_4D_result_5D
    module procedure contract_5D_5D_result_4D

    module procedure contract_6D_2D_result_6D
    module procedure contract_6D_5D_result_5D
    module procedure contract_6D_6D_result_6D

 end interface TENSMUL

 private !! hides all items not listed in public statement 
 public TENSMUL

CONTAINS

 !!! Contract 4D and 4D tensors !!!
 FUNCTION contract_4D_4D_result_4D(tensA, tensBIN, MULT, FUSE) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)   :: tensA(:,:,:,:), tensBIN(:,:,:,:)
  CHARACTER(LEN=*),  INTENT(IN)   :: MULT, FUSE
  COMPLEX(KIND=DP),  ALLOCATABLE  :: tensO(:,:,:,:)

  !! Local copies
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensB(:,:,:,:)
  
  !! Dims & indices
  INTEGER :: dimA(4), dimB(4)
  INTEGER :: ia, ib, ja, jb

  !! Verify flags
  SELECT CASE(MULT)
  CASE('21', '12', '11', '22')

       SELECT CASE(FUSE)
       CASE('(33,44)', '(34,43)')
            CONTINUE
       CASE DEFAULT 
            CALL invalid_flag("contract_4D_4D_result_4D -- invalid FUSE ", FUSE, MULT)
       end SELECT

  CASE('43', '34', '33', '44')

       SELECT CASE(FUSE)
       CASE('(11,22)', '(12,21)')
            CONTINUE
       CASE DEFAULT 
            CALL invalid_flag("contract_4D_4D_result_4D -- invalid FUSE ", FUSE, MULT)
       end SELECT

  CASE DEFAULT
            CALL invalid_flag("contract_4D_4D_result_4D -- invalid MULT ", MULT)
  end SELECT


  !! Transpose tensors if specified by MULT, FUSE flags
  CALL copyTens(tensB, tensBIN)
  IF(MULT .EQ. '11')       tensB = TensTRANSPOSE(tensB, '12') 
  IF(MULT .EQ. '22')       tensB = TensTRANSPOSE(tensB, '12')
  IF(MULT .EQ. '33')       tensB = TensTRANSPOSE(tensB, '34') 
  IF(MULT .EQ. '44')       tensB = TensTRANSPOSE(tensB, '34')
  IF(FUSE .EQ. '(12,21)')  tensB = TensTRANSPOSE(tensB, '12')
  IF(FUSE .EQ. '(34,43)')  tensB = TensTRANSPOSE(tensB, '34')

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('21','22')

     ALLOCATE(tensO(dimA(1), dimB(2), dimA(3)*dimB(3), dimA(4)*dimB(4)))

     DO ia=1,dimA(3)
       DO ib=1,dimB(3)     
         DO ja=1,dimA(4)
           DO jb=1,dimB(4)
              tensO(:, :, ICOM(ia,ib,dimB(3)), ICOM(ja,jb,dimB(4))) = TENSMUL(tensA(:,:,ia,ja), tensB(:,:,ib,jb))
           end DO
         end DO
       end DO
     end DO

  CASE('12','11')

     ALLOCATE(tensO(dimB(1), dimA(2), dimA(3)*dimB(3), dimA(4)*dimB(4)))

     DO ia=1,dimA(3)
       DO ib=1,dimB(3)     
         DO ja=1,dimA(4)
           DO jb=1,dimB(4)
              tensO(:, :, ICOM(ia,ib,dimB(3)), ICOM(ja,jb,dimB(4))) = TENSMUL(tensB(:,:,ib,jb), tensA(:,:,ia,ja))
           end DO
         end DO
       end DO
     end DO

  CASE('43','44')

     ALLOCATE(tensO(dimA(1)*dimB(1), dimA(2)*dimB(2), dimA(3), dimB(4)))

     DO ia=1,dimA(1)
       DO ib=1,dimB(1)     
         DO ja=1,dimA(2)
           DO jb=1,dimB(2)
              tensO(ICOM(ia,ib,dimB(1)), ICOM(ja,jb,dimB(2)), :, :) = TENSMUL(tensA(ia,ja,:,:), tensB(ib,jb,:,:))
           end DO
         end DO
       end DO
     end DO

  CASE('34','33')

     ALLOCATE(tensO(dimA(1)*dimB(1), dimA(2)*dimB(2), dimB(3), dimA(4)))

     DO ia=1,dimA(1)
       DO ib=1,dimB(1)     
         DO ja=1,dimA(2)
           DO jb=1,dimB(2)
              tensO(ICOM(ia,ib,dimB(1)), ICOM(ja,jb,dimB(2)), :, :) = TENSMUL(tensB(ib,jb,:,:), tensA(ia,ja,:,:))
           end DO
         end DO
       end DO
     end DO

  CASE DEFAULT
     CALL invalid_flag("contract_4D_4D_result_4D -- invalid MULT ", MULT) 
  end SELECT

 END FUNCTION contract_4D_4D_result_4D





 !!! Contract 4D and 3D tensors !!!
 FUNCTION contract_4D_3D_result_3D(tensA, tensBIN, MULT, FUSE) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)    :: tensA(:,:,:,:), tensBIN(:,:,:)
  CHARACTER(LEN=*),  INTENT(IN)    :: MULT, FUSE
  COMPLEX(KIND=DP),  ALLOCATABLE   :: tensO(:,:,:)

  !! Local copies
  COMPLEX(KIND=DP), ALLOCATABLE :: tensB(:,:,:)
  
  !! Dims & indices
  INTEGER :: dimA(4), dimB(3)
  INTEGER :: ia, ib, ja, jb

  !! Verify flags
  SELECT CASE(MULT)
  CASE('21','11')

       SELECT CASE(FUSE)
       CASE('(32,43)', '(42,33)')
            CONTINUE
       CASE DEFAULT 
            CALL invalid_flag("contract_4D_3D_result_3D -- invalid FUSE ", FUSE, MULT)
       end SELECT

  CASE('31','41')

       SELECT CASE(FUSE)
       CASE('(12,23)', '(22,13)')
            CONTINUE
       CASE DEFAULT 
            CALL invalid_flag("contract_4D_3D_result_3D -- invalid FUSE ", FUSE, MULT)
       end SELECT

  CASE DEFAULT
            CALL invalid_flag("contract_4D_3D_result_3D -- invalid MULT ", MULT)
  end SELECT


  !! Transpose tensors if specified by MULT, FUSE flags
  CALL copyTens(tensB, tensBIN)
  IF(FUSE .EQ. '(22,13)')  tensB = TensTRANSPOSE(tensB, '23')
  IF(FUSE .EQ. '(42,33)')  tensB = TensTRANSPOSE(tensB, '23')

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('21')

     ALLOCATE(tensO(dimA(1), dimA(3)*dimB(2), dimA(4)*dimB(3)))

     DO ia=1,dimA(3)
       DO ib=1,dimB(2)     
         DO ja=1,dimA(4)
           DO jb=1,dimB(3)
              tensO(:, ICOM(ia,ib,dimB(2)), ICOM(ja,jb,dimB(3))) = TENSMUL(tensA(:,:,ia,ja), tensB(:,ib,jb))
           end DO
         end DO
       end DO
     end DO

  CASE('11')

     ALLOCATE(tensO(dimA(2), dimA(3)*dimB(2), dimA(4)*dimB(3)))

     DO ia=1,dimA(3)
       DO ib=1,dimB(2)     
         DO ja=1,dimA(4)
           DO jb=1,dimB(3)
              tensO(:, ICOM(ia,ib,dimB(2)), ICOM(ja,jb,dimB(3))) = TENSMUL(tensB(:,ib,jb), tensA(:,:,ia,ja))
           end DO
         end DO
       end DO
     end DO

  CASE('41')

     ALLOCATE(tensO(dimA(3), dimA(1)*dimB(2), dimA(2)*dimB(3)))

     DO ia=1,dimA(1)
       DO ib=1,dimB(2)     
         DO ja=1,dimA(2)
           DO jb=1,dimB(3)
              tensO(:, ICOM(ia,ib,dimB(2)), ICOM(ja,jb,dimB(3))) = TENSMUL(tensA(ia,ja,:,:), tensB(:,ib,jb))
           end DO
         end DO
       end DO
     end DO

  CASE('31')

     ALLOCATE(tensO(dimA(4), dimA(1)*dimB(2), dimA(2)*dimB(3)))

     DO ia=1,dimA(1)
       DO ib=1,dimB(2)     
         DO ja=1,dimA(2)
           DO jb=1,dimB(3)
              tensO(:, ICOM(ia,ib,dimB(2)), ICOM(ja,jb,dimB(3))) = TENSMUL(tensB(:,ib,jb), tensA(ia,ja,:,:))
           end DO
         end DO
       end DO
     end DO

  CASE DEFAULT
     CALL invalid_flag("contract_4D_3D_result_3D -- invalid MULT ", MULT) 
  end SELECT

  !! In [contract_4D_3D] we set leg order according to that of tensB
  !! --> so if tensB was transposed, we should now transpose back tensO to restore the original leg order
  IF(FUSE .EQ. '(22,13)')  tensO = TensTRANSPOSE(tensO, '23')
  IF(FUSE .EQ. '(42,33)')  tensO = TensTRANSPOSE(tensO, '23')

 END FUNCTION contract_4D_3D_result_3D







 !!! Contract 4D and 3D tensors !!!
 FUNCTION contract_4D_3D_result_4D(tensA, tensB, MULT) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)  :: tensA(:,:,:,:), tensB(:,:,:)
  CHARACTER(LEN=*),  INTENT(IN)  :: MULT
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensO(:,:,:,:)
  
  !! Dims & indices
  INTEGER :: dimA(4), dimB(3)
  INTEGER :: ia, ib, ja, jb 

  !! Verify flags
  SELECT CASE(MULT)
  CASE('22', '22(AB)', '22(BA)')
       CONTINUE
  CASE DEFAULT
       CALL invalid_flag("contract_4D_3D_result_4D -- invalid MULT ", MULT)
  end SELECT

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('22', '22(AB)')

      ALLOCATE(tensO(dimA(1)*dimB(1), dimB(3), dimA(3), dimA(4)))

      DO ia=1,dimA(1)
        DO ib=1,dimB(1)     
           tensO(ICOM(ia,ib,dimB(1)),:,:,:) = TENSMUL(tensA(ia,:,:,:), tensB(ib,:,:), MULT='11')
        end DO
      end DO

  CASE('22(BA)')

      ALLOCATE(tensO(dimA(1)*dimB(1), dimB(3), dimA(3), dimA(4)))

      DO ib=1,dimB(1)
        DO ia=1,dimA(1)     
           tensO(ICOM(ib,ia,dimA(1)),:,:,:) = TENSMUL(tensA(ia,:,:,:), tensB(ib,:,:), MULT='11')
        end DO
      end DO

  CASE DEFAULT
     CALL invalid_flag("contract_4D_3D_result_4D -- invalid MULT ", MULT) 
  end SELECT

 END FUNCTION contract_4D_3D_result_4D




 !!! Contract 3D and 3D tensors into 4D tensor !!!
 FUNCTION contract_3D_3D_result_4D(tensA, tensBIN, MULT) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)   :: tensA(:,:,:), tensBIN(:,:,:)
  CHARACTER(LEN=*),  INTENT(IN)   :: MULT
  COMPLEX(KIND=DP),  ALLOCATABLE  :: tensO(:,:,:,:)

  !! Local copies
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensB(:,:,:)

  !! Dims & indices
  INTEGER :: dimA(3), dimB(3)
  INTEGER :: ia, ib, ja, jb

  !! Verify flags
  SELECT CASE(MULT)
  CASE('23',  '22')
       CONTINUE
  CASE('32',  '33')
       CONTINUE
  CASE('11')
       CONTINUE
  CASE DEFAULT
       CALL invalid_flag("contract_3D_3D_result_4D -- invalid MULT ", MULT)
  end SELECT


  !! Transpose tensors if specified by MULT, FUSE flags
  CALL copyTens(tensB, tensBIN)
  IF(MULT .EQ. '33')     tensB = TensTRANSPOSE(tensB, '23') 
  IF(MULT .EQ. '22')     tensB = TensTRANSPOSE(tensB, '23')

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('32', '33')

       ALLOCATE(tensO(dimA(1), dimB(1), dimA(2), dimB(3)))
       DO ia=1,dimA(1)
         DO ib=1,dimB(1)
            tensO(ia, ib, :, :) = TENSMUL(tensA(ia,:,:), tensB(ib,:,:))
         end DO
       end DO

  CASE('23', '22')

       ALLOCATE(tensO(dimA(1), dimB(1), dimB(2), dimA(3)))
       DO ia=1,dimA(1)
         DO ib=1,dimB(1)
            tensO(ia, ib, :, :) = TENSMUL(tensB(ib,:,:), tensA(ia,:,:))
         end DO
       end DO

  CASE('11')

       ALLOCATE(tensO(dimA(2), dimA(3), dimB(2), dimB(3)))
       DO ia=1,dimA(2)
         DO ib=1,dimB(2)
           DO ja=1,dimA(3)     
             DO jb=1,dimB(3)
                tensO(ia, ja, ib, jb) = SUM(tensA(:, ia, ja) * tensB(:, ib, jb))
             end DO
           end DO
         end DO
       end DO

  CASE DEFAULT
     CALL invalid_flag("contract_3D_3D_result_4D -- invalid MULT ", MULT) 
  end SELECT

 END FUNCTION contract_3D_3D_result_4D






 !!! Contract 3D and 3D tensors into 2D tensor !!!
 FUNCTION contract_3D_3D_result_2D(tensA, tensBIN, MULT, FUSE) result(tensO)

  COMPLEX(KIND=DP),  INTENT(IN)   :: tensA(:,:,:), tensBIN(:,:,:)
  CHARACTER(LEN=*),  INTENT(IN)   :: MULT, FUSE
  COMPLEX(KIND=DP),  ALLOCATABLE  :: tensO(:,:)

  !! Local copies
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensB(:,:,:)

  !! Dims & indices
  INTEGER :: dimA(3), dimB(3)
  INTEGER :: ia, ib, ja, jb

  !! Verify flags, create local copies
  SELECT CASE(MULT)
  CASE('11')

       SELECT CASE(FUSE)
       CASE('(22,33)', '(23,32)')
            CONTINUE
       CASE DEFAULT 
            CALL invalid_flag("contract_3D_3D_result_2D -- invalid FUSE ", FUSE, MULT)
       end SELECT

  CASE DEFAULT
       CALL invalid_flag("contract_3D_3D_result_2D -- invalid MULT ", MULT)
  end SELECT


  !! Transpose tensors if specified by MULT, FUSE flags
  CALL copyTens(tensB, tensBIN)
  IF(FUSE .EQ. '(23,32)')  tensB = TensTRANSPOSE(tensB, '23') 

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  ALLOCATE(tensO(dimA(2)*dimB(2), dimA(3)*dimB(3)))

  DO ia=1,dimA(2)
    DO ib=1,dimB(2)
      DO ja=1,dimA(3)     
         DO jb=1,dimB(3)
            tensO(ICOM(ia,ib,dimB(2)), ICOM(ja,jb,dimB(3))) = SUM(tensA(:, ia, ja) * tensB(:, ib, jb))
         end DO
      end DO
    end DO
  end DO

 END FUNCTION contract_3D_3D_result_2D






 !!! Contract 3D and 2D tensors into 3D tensor !!!
 FUNCTION contract_3D_2D_result_3D(tensA, tensBIN, MULT) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:), tensBIN(:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: MULT
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:)

  !! Local copies
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensB(:,:)

  !! Dimensions & indices
  INTEGER :: dimA(3), dimB(2)
  INTEGER :: i, j

  !! Verify flags, create local copies
  SELECT CASE(MULT)
  CASE('12',  '11')
       CONTINUE  
  CASE('22',  '21')
       CONTINUE 
  CASE('31',  '32')
       CONTINUE 
  CASE DEFAULT
       CALL invalid_flag("contract_3D_2D_result_3D -- invalid MULT ", MULT)
  end SELECT


  !! Transpose tensors if specified by MULT, FUSE flags
  CALL copyTens(tensB, tensBIN)
  IF(MULT .EQ. '11')  tensB = TensTRANSPOSE(tensB) 
  IF(MULT .EQ. '21')  tensB = TensTRANSPOSE(tensB) 
  IF(MULT .EQ. '32')  tensB = TensTRANSPOSE(tensB) 

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  SELECT CASE(MULT)
  CASE('12','11')

     ALLOCATE(tensO(dimB(1), dimA(2), dimA(3)))
     DO i=1,dimA(2)
        DO j=1,dimA(3)
           tensO(:,i,j) = TENSMUL(tensB(:,:), tensA(:,i,j))
        end DO
     end DO
  
  CASE('22','21')

     ALLOCATE(tensO(dimA(1), dimB(1), dimA(3)))
     DO i=1,dimA(1)
        tensO(i,:,:) = TENSMUL(tensB(:,:), tensA(i,:,:))
     end DO

  CASE('31','32')

     ALLOCATE(tensO(dimA(1), dimA(2), dimB(2)))
     DO i=1,dimA(1)
        tensO(i,:,:) = TENSMUL(tensA(i,:,:), tensB(:,:))
     end DO
 
  CASE DEFAULT
     CALL invalid_flag("contract_3D_2D_result_3D -- invalid MULT ", MULT)
  end SELECT

 END FUNCTION contract_3D_2D_result_3D







 !!! Contract 2D and 3D tensors into 2D tensor !!!
 FUNCTION contract_3D_2D_result_2D(tensA, tensBIN, MULT, FUSE) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:), tensBIN(:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: MULT, FUSE
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:)

  !! Local copies
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensB(:,:)

  !! Dimensions & indices
  INTEGER :: dimA(3), dimB(2)
  INTEGER :: ia, ib

  !! Verify flags, create local copies
  SELECT CASE(MULT)
  CASE('22', '21', '31', '32')

       SELECT CASE(FUSE)
       CASE('(11)', '(12)')
            CONTINUE
       CASE DEFAULT 
            CALL invalid_flag("contract_3D_2D_result_2D -- invalid FUSE ", FUSE, MULT)
       end SELECT
 
  CASE DEFAULT
       CALL invalid_flag("contract_3D_2D_result_2D -- invalid MULT ", MULT)
  end SELECT


  !! Transpose tensors if specified by MULT, FUSE flags
  CALL copyTens(tensB, tensBIN)
  IF(MULT .EQ. '21')  tensB = TensTRANSPOSE(tensB) 
  IF(MULT .EQ. '32')  tensB = TensTRANSPOSE(tensB) 

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  SELECT CASE(MULT)
  CASE('22', '21')

       ALLOCATE(tensO(dimA(1)*dimB(1), dimA(3)))
       DO ia=1,dimA(1)
         DO ib=1,dimB(1)
            tensO(ICOM(ia,ib,dimB(1)), :) = TENSMUL(tensB(ib, :), tensA(ia, :, :))
         end DO
       end DO

  CASE('31', '32')

       ALLOCATE(tensO(dimA(2), dimA(1)*dimB(2)))
       DO ia=1,dimA(1)
         DO ib=1,dimB(2)
            tensO(:, ICOM(ia,ib,dimB(2))) = TENSMUL(tensA(ia, :, :), tensB(:, ib))
         end DO
       end DO

  CASE DEFAULT
      CALL invalid_flag("contract_3D_2D_result_2D -- invalid MULT ", MULT)
  end SELECT
  
 END FUNCTION contract_3D_2D_result_2D





 !!! Contract 2D && 2D tensors !!!
 FUNCTION contract_2D_2D_result_2D_transp(tensA, tensBIN, MULT) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:), tensBIN(:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: MULT
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:)

  !! Local copies
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensB(:,:)

  !! Verify flags, create local copies
  SELECT CASE(MULT)
  CASE('12', '21', '22', '11')
       CONTINUE
  CASE DEFAULT
       CALL invalid_flag("contract_2D_2D_result_2D_transp -- invalid MULT ", MULT)
  end SELECT

  !! Transpose tensors if specified by MULT, FUSE flags
  CALL copyTens(tensB, tensBIN)
  IF(MULT .EQ. '22')  tensB = TensTRANSPOSE(tensB) 
  IF(MULT .EQ. '11')  tensB = TensTRANSPOSE(tensB) 

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('21','22')
       tensO = TENSMUL(tensA, tensB)
  CASE('11','12')
       tensO = TENSMUL(tensB, tensA)
  CASE DEFAULT
       CALL invalid_flag("contract_2D_2D_result_2D_transp -- invalid MULT ", MULT)
  end SELECT

 END FUNCTION contract_2D_2D_result_2D_transp






 !!! Contract 2D && 2D tensors !!!
 FUNCTION contract_2D_2D_result_2D(tensA, tensB) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:), tensB(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:)

  !! Dimensions & indices
  INTEGER :: dimA(2), dimB(2)
  
  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Check matrix dims are consistent
  CALL check_sizes_equal(dimA(2), dimB(1), "contract_2D_2D_result_2D -- tensA && tensB must have compatible dims ")

  !! Allocate & initialize tensO
  ALLOCATE(tensO(dimA(1), dimB(2)))

  !! Contract tensors
  tensO = MATMUL(tensA, tensB)

 END FUNCTION contract_2D_2D_result_2D







 !!! Contract 2D && 1D tensors !!!
 FUNCTION contract_2D_1D_result_1D(tensA, tensB) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:), tensB(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:)

  !! Dimensions & indices
  INTEGER :: dimA(2), dimB(1)
  
  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Check matrix dims are consistent
  CALL check_sizes_equal(dimA(2), dimB(1), "contract_2D_1D_result_1D: tensA && tensB must have compatible dims ")

  !! Allocate tensO
  ALLOCATE(tensO(dimA(1)))

  !! Contract tensors
  tensO = MATMUL(tensA, tensB)

 END FUNCTION contract_2D_1D_result_1D







 !!! Contract 1D && 2D tensors !!!
 !! NB. argument names in contract_1D_2D MUST BE DISTINCT from those in contract_2D_1D
 !! for generic function interface to work properly -- hence tensAA, tensBB.
 FUNCTION contract_1D_2D_result_1D(tensAA, tensBB) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensAA(:), tensBB(:,:) !DO NOT RENAME to tensA, tensB (see above)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:)

  !! Dimensions & indices
  INTEGER :: dimA(1), dimB(2)
  
  !! Get dims of tensA, tensB
  dimA = shape(tensAA); dimB = shape(tensBB)

  !! Check matrix dims are consistent
  CALL check_sizes_equal(dimA(1), dimB(1), "contract_1D_2D_result_1D: tensA && tensB must have compatible dims ")

  !! Allocate tensO
  ALLOCATE(tensO(dimB(2)))

  !! Contract tensors
  tensO = MATMUL(tensAA, tensBB)

 END FUNCTION contract_1D_2D_result_1D









 !!! Contract 4D && 2D tensors !!!
 FUNCTION contract_4D_2D_result_4D(tensA, tensBIN, MULT) result(tensO) 

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:,:), tensBIN(:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: MULT
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:)

  !! Local copies
  COMPLEX(KIND=DP), ALLOCATABLE :: tensB(:,:)

  !! Dims & indices
  INTEGER :: dimA(4), dimB(2)
  INTEGER :: i, j

  !! Verify flags, create local copies
  SELECT CASE(MULT)
  CASE('12',  '11')
       CONTINUE 
  CASE('21',  '22')
       CONTINUE  
  CASE('32',  '31')
       CONTINUE  
  CASE('41',  '42')
       CONTINUE   
  CASE DEFAULT
       CALL invalid_flag("contract_4D_2D_result_4D -- invalid MULT ", MULT)
  end SELECT

  !! Transpose tensors if specified by MULT, FUSE flags
  CALL copyTens(tensB, tensBIN)
  IF(MULT .EQ. '11')  tensB = TensTRANSPOSE(tensB) 
  IF(MULT .EQ. '22')  tensB = TensTRANSPOSE(tensB) 
  IF(MULT .EQ. '31')  tensB = TensTRANSPOSE(tensB) 
  IF(MULT .EQ. '42')  tensB = TensTRANSPOSE(tensB)

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('12','11')

     ALLOCATE(tensO(dimB(1), dimA(2), dimA(3), dimA(4)))
     DO i=1,dimA(3)
       DO j=1,dimA(4)
          tensO(:,:,i,j) = TENSMUL(tensB(:,:), tensA(:,:,i,j))
       end DO
     end DO

  CASE('21','22')

     ALLOCATE(tensO(dimA(1), dimB(2), dimA(3), dimA(4)))
     DO i=1,dimA(3)
       DO j=1,dimA(4)
          tensO(:,:,i,j) = TENSMUL(tensA(:,:,i,j), tensB(:,:))
       end DO
     end DO

  CASE('32','31')

     ALLOCATE(tensO(dimA(1), dimA(2), dimB(1), dimA(4)))
     DO i=1,dimA(1)
       DO j=1,dimA(2)
          tensO(i,j,:,:) = TENSMUL(tensB(:,:), tensA(i,j,:,:))
       end DO
     end DO

  CASE('41','42')

     ALLOCATE(tensO(dimA(1), dimA(2), dimA(3), dimB(2)))
     DO i=1,dimA(1)
       DO j=1,dimA(2)
          tensO(i,j,:,:) = TENSMUL(tensA(i,j,:,:), tensB(:,:))
       end DO
     end DO

  CASE DEFAULT
     CALL invalid_flag("contract_4D_2D_result_4D -- invalid MULT ", MULT)
  end SELECT

 END FUNCTION contract_4D_2D_result_4D







 !!! Contract 5D && 2D tensors !!!
 FUNCTION contract_5D_2D_result_5D(tensA, tensBIN, MULT) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)   :: tensA(:,:,:,:,:), tensBIN(:,:)
  CHARACTER(LEN=*), INTENT(IN)   :: MULT
  COMPLEX(KIND=DP), ALLOCATABLE  :: tensO(:,:,:,:,:)
  
  !! Local copies
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensB(:,:)
  COMPLEX(KIND=DP),  ALLOCATABLE :: tempTens(:,:,:)

  !! Dims & indices
  INTEGER :: dimA(5), dimB(2)
  INTEGER :: i, j

  !! Verify flags, create local copies
  SELECT CASE(MULT)
  CASE('12',  '11')
       CONTINUE
  CASE('22',  '21')
       CONTINUE
  CASE('31',  '32')
       CONTINUE
  CASE('42',  '41')
       CONTINUE
  CASE('51',  '52')
       CONTINUE
  CASE DEFAULT
       CALL invalid_flag("contract_5D_2D_result_5D -- invalid MULT ", MULT)
  end SELECT

  !! Transpose tensors if specified by MULT, FUSE flags
  CALL copyTens(tensB, tensBIN)
  IF(MULT .EQ. '11')  tensB = TensTRANSPOSE(tensB) 
  IF(MULT .EQ. '21')  tensB = TensTRANSPOSE(tensB) 
  IF(MULT .EQ. '32')  tensB = TensTRANSPOSE(tensB) 
  IF(MULT .EQ. '41')  tensB = TensTRANSPOSE(tensB) 
  IF(MULT .EQ. '52')  tensB = TensTRANSPOSE(tensB)

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('12','11')

     ALLOCATE(tensO(dimB(1), dimA(2), dimA(3), dimA(4), dimA(5)))
     DO i=1,dimA(2)
       DO j=1,dimA(3)
          tempTens = TENSMUL(tensA(:,i,j,:,:), tensB, '12')
          tensO(:,i,j,:,:) = tempTens 
       end DO
     end DO

  CASE('22','21')

     ALLOCATE(tensO(dimA(1), dimB(1), dimA(3), dimA(4), dimA(5)))
     DO i=1,dimA(4)
       DO j=1,dimA(5)
          tempTens = TENSMUL(tensA(:,:,:,i,j), tensB, '22') 
          tensO(:,:,:,i,j) = tempTens 
       end DO
     end DO  

  CASE('31','32')

     ALLOCATE(tensO(dimA(1), dimA(2), dimB(2), dimA(4), dimA(5)))
     DO i=1,dimA(4)
       DO j=1,dimA(5)
          tempTens = TENSMUL(tensA(:,:,:,i,j), tensB, '31')
          tensO(:,:,:,i,j) = tempTens 
       end DO
     end DO  

  CASE('42','41')

     ALLOCATE(tensO(dimA(1), dimA(2), dimA(3), dimB(1), dimA(5)))
     DO i=1,dimA(2)
       DO j=1,dimA(3)
          tempTens = TENSMUL(tensA(:,i,j,:,:), tensB, '22') 
          tensO(:,i,j,:,:) = tempTens 
       end DO
     end DO  

  CASE('51','52')

     ALLOCATE(tensO(dimA(1), dimA(2), dimA(3), dimA(4), dimB(2)))
     DO i=1,dimA(2)
       DO j=1,dimA(3)
          tempTens = TENSMUL(tensA(:,i,j,:,:), tensB, '31') 
          tensO(:,i,j,:,:) = tempTens 
       end DO
     end DO  

  CASE DEFAULT
     CALL invalid_flag("contract_5D_2D_result_5D -- invalid MULT ", MULT)
  end SELECT 

 END FUNCTION contract_5D_2D_result_5D







 !!! Contract 5D && 4D tensors !!!
 FUNCTION contract_5D_4D_result_5D(tensA, tensBIN, MULT, FUSE) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:,:,:), tensBIN(:,:,:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: MULT, FUSE
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:,:)
  
  !! Local copies
  COMPLEX(KIND=DP), ALLOCATABLE :: tensB(:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tempTens(:,:,:)

  !! Dims & indices
  INTEGER :: dimA(5), dimB(4)
  INTEGER :: i, j

  !! Verify flags, create local copies
  SELECT CASE(MULT)
  CASE('12', '11')

       SELECT CASE(FUSE)
       CASE('(43,54)', '(23,34)')
            CONTINUE
       CASE DEFAULT 
            CALL invalid_flag("contract_5D_4D_result_5D -- invalid FUSE ", FUSE, MULT)
       end SELECT
 
  CASE DEFAULT
       CALL invalid_flag("contract_5D_4D_result_5D -- invalid MULT ", MULT)
  end SELECT


  !! Transpose tensors if specified by MULT, FUSE flags
  CALL copyTens(tensB, tensBIN)
  IF(MULT .EQ. '11')  tensB = TensTRANSPOSE(tensB, '12') 

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Contract tensors
  SELECT CASE(FUSE)
  CASE('(43,54)')
 
     ALLOCATE(tensO(dimB(1), dimA(2), dimA(3), dimA(4)*dimB(3), dimA(5)*dimB(4)))
     DO i=1,dimA(2)
       DO j=1,dimA(3)
          tempTens = TENSMUL(tensB, tensA(:,i,j,:,:), MULT='21', FUSE='(32,43)')
          tensO(:,i,j,:,:) = tempTens 
       end DO
     end DO

  CASE('(23,34)') 

     ALLOCATE(tensO(dimB(1), dimA(2)*dimB(3), dimA(3)*dimB(4), dimA(4), dimA(5)))
     DO i=1,dimA(4)
       DO j=1,dimA(5)
          tempTens = TENSMUL(tensB, tensA(:,:,:,i,j), MULT='21', FUSE='(32,43)')
          tensO(:,:,:,i,j) = tempTens 
       end DO
     end DO  

  CASE DEFAULT
     CALL invalid_flag("contract_5D_4D_result_5D -- invalid FUSE ", FUSE)
  end SELECT 

 END FUNCTION contract_5D_4D_result_5D







 !!! Contract 5D && 5D tensors !!!
 FUNCTION contract_5D_5D_result_4D(tensA, tensB, MULT, FUSE) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:,:,:), tensB(:,:,:,:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: MULT, FUSE
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:)

  !! tensA and tensB reshaped into 3D tensors 
  COMPLEX(KIND=DP), ALLOCATABLE :: tensA_3D(:,:,:), tensB_3D(:,:,:)

  !! Dims & indices
  INTEGER :: dimA(5), ia, ja, ka, la
  INTEGER :: dimB(5), ib, jb, kb, lb

  !! Verify flags
  SELECT CASE(MULT)
  CASE('54')

       SELECT CASE(FUSE)
       CASE('(33,22)')
            CONTINUE
       CASE DEFAULT 
            CALL invalid_flag("contract_5D_5D_result_4D -- invalid FUSE ", FUSE, MULT)
       end SELECT

  CASE('32')

       SELECT CASE(FUSE)
       CASE('(44,55)')
            CONTINUE
       CASE DEFAULT 
            CALL invalid_flag("contract_5D_5D_result_4D -- invalid FUSE ", FUSE, MULT)
       end SELECT

  CASE('11')

       SELECT CASE(FUSE)
       CASE('(22,33,44,55)')
            CONTINUE
       CASE DEFAULT 
            CALL invalid_flag("contract_5D_5D_result_4D -- invalid FUSE ", FUSE, MULT)
       end SELECT

  CASE DEFAULT
            CALL invalid_flag("contract_5D_5D_result_4D -- invalid MULT ", MULT)
  end SELECT


  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)


  !! Contract tensors
  SELECT CASE(MULT)
  CASE('54') 

     !! Reshape 5D into 3D
     tensA_3D = RESHAPE_3D(tensA, '1,234,5')  
     tensB_3D = RESHAPE_3D(tensB, '1,4,235') 

     !! Contract
     tensO = TENSMUL(tensA_3D, tensB_3D, '32')

  CASE('32')

     !! Reshape 5D into 3D
     tensA_3D = RESHAPE_3D(tensA, '1,245,3')  
     tensB_3D = RESHAPE_3D(tensB, '1,2,345')

     !! Contract
     tensO = TENSMUL(tensA_3D, tensB_3D, '32')

  CASE('11')

     ALLOCATE(tensO(dimA(2)*dimB(2), dimA(3)*dimB(3), dimA(4)*dimB(4), dimA(5)*dimB(5)))

     DO ia=1,dimA(2)
       DO ib=1,dimB(2)     
         DO ja=1,dimA(3)
           DO jb=1,dimB(3)
             DO ka=1,dimA(4)
               DO kb=1,dimB(4)     
                 DO la=1,dimA(5)
                   DO lb=1,dimB(5)

                      tensO(ICOM(ia,ib,dimB(2)), ICOM(ja,jb,dimB(3)), ICOM(ka,kb,dimB(4)), ICOM(la,lb,dimB(5))) &
                      & = SUM(tensA(:,ia,ja,ka,la) * tensB(:,ib,jb,kb,lb))                           
                   end DO
                 end DO
               end DO
             end DO
           end DO
         end DO
       end DO
     end DO

  CASE DEFAULT
     CALL invalid_flag("contract_5D_5D_result_4D -- invalid MULT ", MULT)
  end SELECT

 END FUNCTION contract_5D_5D_result_4D








 !!! Contract 6D && 2D tensors !!!
 FUNCTION contract_6D_2D_result_6D(tensA, tensBIN, MULT) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)   :: tensA(:,:,:,:,:,:), tensBIN(:,:)
  CHARACTER(LEN=*), INTENT(IN)   :: MULT
  COMPLEX(KIND=DP), ALLOCATABLE  :: tensO(:,:,:,:,:,:)
  
  !! Local copies
  COMPLEX(KIND=DP),  ALLOCATABLE :: tensB(:,:)
  COMPLEX(KIND=DP),  ALLOCATABLE :: tempTens(:,:,:,:)

  !! Dims & indices
  INTEGER :: dimA(6), dimB(2)
  INTEGER :: i, j

  !! Verify flags, create local copies
  SELECT CASE(MULT)
  CASE('12',  '11')
       CONTINUE
  CASE('22',  '21')
       CONTINUE
  CASE DEFAULT
       CALL invalid_flag("contract_6D_2D_result_6D -- invalid MULT ", MULT)
  end SELECT

  !! Transpose tensors if specified by MULT, FUSE flags
  CALL copyTens(tensB, tensBIN)
  IF(MULT .EQ. '11')  tensB = TensTRANSPOSE(tensB) 
  IF(MULT .EQ. '22')  tensB = TensTRANSPOSE(tensB) 

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('12','11')

     ALLOCATE(tensO(dimB(1), dimA(2), dimA(3), dimA(4), dimA(5), dimA(6)))
     DO i=1,dimA(3)
       DO j=1,dimA(4)
          tempTens = TENSMUL(tensA(:,:,i,j,:,:), tensB, '12')
          tensO(:,:,i,j,:,:) = tempTens 
       end DO
     end DO

  CASE('22','21')

     ALLOCATE(tensO(dimA(1), dimB(2), dimA(3), dimA(4), dimA(5), dimA(6)))
     DO i=1,dimA(3)
       DO j=1,dimA(4)
          tempTens = TENSMUL(tensA(:,:,i,j,:,:), tensB, '21')
          tensO(:,:,i,j,:,:) = tempTens 
       end DO
     end DO  

  CASE DEFAULT
     CALL invalid_flag("contract_6D_2D_result_6D -- invalid MULT ", MULT)
  end SELECT 

 END FUNCTION contract_6D_2D_result_6D






 !!! Contract 6D && 5D tensors !!!
 FUNCTION contract_6D_5D_result_5D(tensA, tensB, MULT, FUSE) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:,:,:,:), tensB(:,:,:,:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: MULT, FUSE
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:,:)

  !! Dims & indices
  INTEGER :: dimA(6), ia, ja, ka, la
  INTEGER :: dimB(5), ib, jb, kb, lb

  !! Verify flags
  SELECT CASE(MULT)
  CASE('21', '11')
       SELECT CASE(FUSE)
       CASE('(32,43,54,65)')
            CONTINUE
       CASE DEFAULT 
            CALL invalid_flag("contract_6D_5D_result_5D -- invalid FUSE ", FUSE, MULT)
       end SELECT
  CASE DEFAULT
       CALL invalid_flag("contract_6D_5D_result_5D -- invalid MULT ", MULT)
  end SELECT

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('21')

     ALLOCATE(tensO(dimA(1), dimA(3)*dimB(2), dimA(4)*dimB(3), dimA(5)*dimB(4), dimA(6)*dimB(5)))

     DO ia=1,dimA(3)
       DO ib=1,dimB(2)     
         DO ja=1,dimA(4)
           DO jb=1,dimB(3)
             DO ka=1,dimA(5)
               DO kb=1,dimB(4)     
                 DO la=1,dimA(6)
                   DO lb=1,dimB(5)

                      tensO(:, ICOM(ia,ib,dimB(2)), ICOM(ja,jb,dimB(3)), ICOM(ka,kb,dimB(4)), ICOM(la,lb,dimB(5))) &
                      & = TENSMUL(tensA(:,:,ia,ja,ka,la),  tensB(:,ib,jb,kb,lb))                           
                   end DO
                 end DO
               end DO
             end DO
           end DO
         end DO
       end DO
     end DO

  CASE('11')

     ALLOCATE(tensO(dimA(2), dimA(3)*dimB(2), dimA(4)*dimB(3), dimA(5)*dimB(4), dimA(6)*dimB(5)))

     DO ia=1,dimA(3)
       DO ib=1,dimB(2)     
         DO ja=1,dimA(4)
           DO jb=1,dimB(3)
             DO ka=1,dimA(5)
               DO kb=1,dimB(4)     
                 DO la=1,dimA(6)
                   DO lb=1,dimB(5)

                      tensO(:, ICOM(ia,ib,dimB(2)), ICOM(ja,jb,dimB(3)), ICOM(ka,kb,dimB(4)), ICOM(la,lb,dimB(5))) &
                      & = TENSMUL(tensB(:,ib,jb,kb,lb),  tensA(:,:,ia,ja,ka,la))                           
                   end DO
                 end DO
               end DO
             end DO
           end DO
         end DO
       end DO
     end DO

  CASE DEFAULT
     CALL invalid_flag("contract_6D_5D_result_5D -- invalid MULT ", MULT)
  end SELECT

 END FUNCTION contract_6D_5D_result_5D




 !!! Contract 6D && 6D tensors !!!
 FUNCTION contract_6D_6D_result_6D(tensA, tensB, MULT, FUSE) result(tensO)

  COMPLEX(KIND=DP), INTENT(IN)  :: tensA(:,:,:,:,:,:), tensB(:,:,:,:,:,:)
  CHARACTER(LEN=*), INTENT(IN)  :: MULT, FUSE
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:,:,:,:,:,:)

  !! Dims & indices
  INTEGER :: dimA(6), ia, ja, ka, la
  INTEGER :: dimB(6), ib, jb, kb, lb

  !! Verify flags
  SELECT CASE(MULT)
  CASE('21')
       SELECT CASE(FUSE)
       CASE('(33,44,55,66)')
            CONTINUE
       CASE DEFAULT 
            CALL invalid_flag("contract_6D_6D_result_6D -- invalid FUSE ", FUSE, MULT)
       end SELECT
  CASE DEFAULT
       CALL invalid_flag("contract_6D_6D_result_6D -- invalid MULT ", MULT)
  end SELECT

  !! Get dims of tensA, tensB
  dimA = shape(tensA); dimB = shape(tensB)

  !! Contract tensors
  SELECT CASE(MULT)
  CASE('21')

     ALLOCATE(tensO(dimA(1), dimB(2), dimA(3)*dimB(3), dimA(4)*dimB(4), dimA(5)*dimB(5), dimA(6)*dimB(6)))

     DO ia=1,dimA(3)
       DO ib=1,dimB(3)     
         DO ja=1,dimA(4)
           DO jb=1,dimB(4)
             DO ka=1,dimA(5)
               DO kb=1,dimB(5)     
                 DO la=1,dimA(6)
                   DO lb=1,dimB(6)

                      tensO(:, :, ICOM(ia,ib,dimB(3)), ICOM(ja,jb,dimB(4)), ICOM(ka,kb,dimB(5)), ICOM(la,lb,dimB(6))) &
                      & = TENSMUL(tensA(:,:,ia,ja,ka,la),  tensB(:,:,ib,jb,kb,lb)) 
                          
                   end DO
                 end DO
               end DO
             end DO
           end DO
         end DO
       end DO
     end DO

  CASE DEFAULT
     CALL invalid_flag("contract_6D_6D_result_6D -- invalid MULT ", MULT)
  end SELECT

 END FUNCTION contract_6D_6D_result_6D

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE TENS_mult
