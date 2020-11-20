MODULE error_handling 

 USE basic_functions

 IMPLICIT NONE

 private
 public check_sizes_equal, check_space_exists !! argument/flag checks
 public invalid_flag, invalid_value
 public call_exit

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHECKING FUNCTION ARGUMENTS && FLAGS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Stop && throw error if an invalid flag is passed !!!
 SUBROUTINE invalid_flag(msg, incorrectFlag, incorrectFlag_1)

  CHARACTER(LEN=*), INTENT(IN)           :: msg
  CHARACTER(LEN=*), INTENT(IN)           :: incorrectFlag
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: incorrectFlag_1

  !! Stop && throw error if an incorrect flag is passed
  IF(PRESENT(incorrectFlag_1)) THEN

     WRITE(*,*) msg, incorrectFlag, " ", incorrectFlag_1 
     STOP
  ELSE
     WRITE(*,*) msg, incorrectFlag 
     STOP
  end IF

 END SUBROUTINE invalid_flag



 !!! Stop && throw error if an invalid value is passed !!!
 SUBROUTINE invalid_value(msg, incorrectVal, incorrectVal_1)

  CHARACTER(LEN=*), INTENT(IN)           :: msg
  INTEGER,          INTENT(IN)           :: incorrectVal
  INTEGER,          INTENT(IN), OPTIONAL :: incorrectVal_1

  !! Stop && throw error if an incorrect flag is passed
  IF(PRESENT(incorrectVal_1)) THEN
     WRITE(*,*) msg, incorrectVal, " ", incorrectVal_1 
     STOP
  ELSE
     WRITE(*,*) msg, incorrectVal
     STOP
  end IF

 END SUBROUTINE invalid_value



 !!! Check sizes are consistent !!!
 SUBROUTINE check_sizes_equal(sizeA, sizeB, msg)

  INTEGER,          INTENT(IN) :: sizeA, sizeB
  CHARACTER(LEN=*), INTENT(IN) :: msg

  !! Consistency check
  IF(sizeA .NE. sizeB) THEN
     WRITE(*,*) msg, sizeA, sizeB
     STOP
  end IF

 END SUBROUTINE check_sizes_equal



 !!! Check space exists !!!
 SUBROUTINE check_space_exists(desired_size, max_size, msg)

  INTEGER,          INTENT(IN) :: desired_size, max_size
  CHARACTER(LEN=*), INTENT(IN) :: msg

  !! Consistency check
  IF(desired_size .GT. max_size) THEN
     WRITE(*,*) msg, desired_size, max_size
     STOP
  end IF

 END SUBROUTINE check_space_exists


 !!! Function that triggers loop exit if exit_cond is true !!!
 !!! && prints an appropriate exit message                 !!! 
 FUNCTION call_exit(exit_cond, msg)

  LOGICAL,          INTENT(IN) :: exit_cond
  CHARACTER(LEN=*), INTENT(IN) :: msg
  LOGICAL                      :: call_exit

  call_exit = exit_cond

  IF(call_exit) WRITE(*,*) msg
  
 END FUNCTION call_exit


END MODULE error_handling
