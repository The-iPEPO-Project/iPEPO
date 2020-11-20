MODULE basic_functions

 IMPLICIT NONE

 INTEGER,          PARAMETER :: DP=KIND(1.0d0)
 COMPLEX(KIND=DP), PARAMETER :: ii=(0.0D0, 1.0D0)
 COMPLEX(KIND=DP), PARAMETER :: one=(1.0D0, 0.0D0)

 INTEGER :: DUMP_INDEX

 interface OptArg
    module procedure OptArg_Char
    module procedure OptArg_Log
    module procedure OptArg_Int
    module procedure OptArg_Re
    module procedure OptArg_CMPLX
 end interface OptArg

 interface setOptVar
    module procedure setOptVar_Char
    module procedure setOptVar_Log
    module procedure setOptVar_Int
    module procedure setOptVar_Re
    module procedure setOptVar_CMPLX
 end interface setOptVar

 interface complx
    module procedure cmplxArr
    module procedure cmplxVal
 end interface complx

 interface logN
    module procedure logN_cmplx
    module procedure logN_real
 end interface logN

 interface swap
    module procedure swapInt
    module procedure swapRe
    module procedure swapCmplx
 end interface swap

 private
 public DP, ii, DUMP_INDEX
 public OptArg, setOptVar, complx
 public isEVEN, isODD, isTrueVec, isTrueMat !basic logical functions
 public logN, SQRTint
 public swap, sort_array

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! READING OPTIONAL ARGUMENTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !set optional args (CHARACTER)
 FUNCTION OptArg_Char(argIn, argDefault) result(arg)

  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: argIn
  CHARACTER(LEN=*), INTENT(IN) :: argDefault

  CHARACTER :: arg*128 

  IF(PRESENT(argIn)) THEN
     arg = argIn !set to input val
  ELSE
     arg = argDefault !set to default val
  end IF
 
 END FUNCTION OptArg_Char



 !set optional args (LOGICAL)
 FUNCTION OptArg_Log(argIn, argDefault) result(arg)

  LOGICAL, INTENT(IN), OPTIONAL :: argIn
  LOGICAL, INTENT(IN) :: argDefault

  LOGICAL :: arg 

  IF(PRESENT(argIn)) THEN
     arg = argIn !set to input val
  ELSE
     arg = argDefault !set to default val
  end IF
 
 END FUNCTION OptArg_Log




 !set optional args (INTEGER)
 FUNCTION OptArg_Int(argIn, argDefault) result(arg)

  INTEGER, INTENT(IN), OPTIONAL :: argIn
  INTEGER, INTENT(IN) :: argDefault

  INTEGER :: arg 

  IF(PRESENT(argIn)) THEN
     arg = argIn !set to input val
  ELSE
     arg = argDefault !set to default val
  end IF
 
 END FUNCTION OptArg_Int




 !set optional args (REAL)
 FUNCTION OptArg_Re(argIn, argDefault) result(arg)

  REAL(KIND=DP), INTENT(IN), OPTIONAL :: argIn
  REAL(KIND=DP), INTENT(IN) :: argDefault

  REAL(KIND=DP) :: arg

  IF(PRESENT(argIn)) THEN
     arg = argIn !set to input val
  ELSE
     arg = argDefault !set to default val
  end IF
 
 END FUNCTION OptArg_Re




 !set optional args (COMPLEX)
 FUNCTION OptArg_CMPLX(argIn, argDefault) result(arg)

  COMPLEX(KIND=DP), INTENT(IN), OPTIONAL :: argIn
  COMPLEX(KIND=DP), INTENT(IN) :: argDefault

  COMPLEX(KIND=DP) :: arg

  IF(PRESENT(argIn)) THEN
     arg = argIn !set to input val
  ELSE
     arg = argDefault !set to default val
  end IF
 
 END FUNCTION OptArg_CMPLX

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! setOptVar: SETTING OPTIONAL ARGUMENTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !Set opt variable, make sure it exists
 SUBROUTINE setOptVar_Char(var, varIn, msg)

  CHARACTER(LEN=*), INTENT(INOUT), OPTIONAL :: var
  CHARACTER(LEN=*), INTENT(IN)              :: varIn
  CHARACTER(LEN=*), INTENT(IN)              :: msg

  IF(PRESENT(var)) THEN
     var = varIn !set to input val
  ELSE
     WRITE(*,*) msg
     STOP
  end IF

 END SUBROUTINE setOptVar_Char



 !Set opt variable, make sure it exists
 SUBROUTINE setOptVar_Log(var, varIn, msg)

  LOGICAL,          INTENT(INOUT), OPTIONAL :: var
  LOGICAL,          INTENT(IN)              :: varIn
  CHARACTER(LEN=*), INTENT(IN)              :: msg

  IF(PRESENT(var)) THEN
     var = varIn !set to input val
  ELSE
     WRITE(*,*) msg
     STOP
  end IF

 END SUBROUTINE setOptVar_Log



 !Set opt variable, make sure it exists
 SUBROUTINE setOptVar_Int(var, varIn, msg) 

  INTEGER,          INTENT(INOUT), OPTIONAL :: var
  INTEGER,          INTENT(IN)              :: varIn
  CHARACTER(LEN=*), INTENT(IN)              :: msg

  IF(PRESENT(var)) THEN
     var = varIn !set to input val
  ELSE
     WRITE(*,*) msg
     STOP
  end IF

 END SUBROUTINE setOptVar_Int



 !Set opt variable, make sure it exists
 SUBROUTINE setOptVar_Re(var, varIn, msg)

  REAL(KIND=DP),    INTENT(INOUT), OPTIONAL :: var
  REAL(KIND=DP),    INTENT(IN)              :: varIn
  CHARACTER(LEN=*), INTENT(IN)              :: msg

  IF(PRESENT(var)) THEN
     var = varIn !set to input val
  ELSE
     WRITE(*,*) msg
     STOP
  end IF

 END SUBROUTINE setOptVar_Re



 !Set opt variable, make sure it exists
 SUBROUTINE setOptVar_CMPLX(var, varIn, msg)

  COMPLEX(KIND=DP), INTENT(INOUT), OPTIONAL :: var
  COMPLEX(KIND=DP), INTENT(IN)              :: varIn
  CHARACTER(LEN=*), INTENT(IN)              :: msg

  IF(PRESENT(var)) THEN
     var = varIn !set to input val
  ELSE
     WRITE(*,*) msg
     STOP
  end IF

 END SUBROUTINE setOptVar_CMPLX

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BASIC LOGICAL FUNCTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !Check if num is odd
 FUNCTION isODD(num_in)

  INTEGER, INTENT(IN) :: num_in !input number (index, site, etc)
  LOGICAL :: isODD !result

  isODD = .NOT. isEVEN(num_in)

 END FUNCTION isODD



 !Check if num is even
 FUNCTION isEVEN(num_in)
 
  INTEGER, INTENT(IN) :: num_in !input number (index, site, etc)
  LOGICAL :: isEVEN !result

  IF(mod(num_in,2) .EQ. 0) THEN
     isEVEN = .TRUE.
  ELSE
     isEVEN = .FALSE.
  end IF
       
 END FUNCTION isEVEN




 !Check if the whole vector is true
 FUNCTION isTrueVec(vec, log_op)

  LOGICAL,          INTENT(IN) :: vec(:)     !input logical vec
  CHARACTER(LEN=*), INTENT(IN) :: log_op     !logical operator to be used
  LOGICAL                      :: isTrueVec  !result

  !vector size
  INTEGER :: nvec

  !count elements that satisfy the condition
  IF(log_op .EQ. 'AND') THEN

     !check if all elements satisfy the condition
     nvec = SIZE(vec)
     isTrueVec = (COUNT(vec) .EQ. nvec)

  ELSEIF(log_op .EQ. 'OR') THEN

     !check if at least one element satisfies the condition
     isTrueVec = (COUNT(vec) .GT. 0)
  ELSE
     WRITE(*,*) "isTrueVec -- invalid log_op ", log_op
     STOP
  end IF

 END FUNCTION isTrueVec




 !Check if the whole matrix is true
 FUNCTION isTrueMat(mat, log_op)

  LOGICAL,          INTENT(IN) :: mat(:,:)   !input logical vec
  CHARACTER(LEN=*), INTENT(IN) :: log_op     !logical operator to be used
  LOGICAL                      :: isTrueMat  !result

  !matrix size (total num of elements)
  INTEGER :: nmat

  !count elements that satisfy the condition
  IF(log_op .EQ. 'AND') THEN

     !check if all elements satisfy the condition
     nmat = SIZE(mat)
     isTrueMat = (COUNT(mat) .EQ. nmat)

  ELSEIF(log_op .EQ. 'OR') THEN

     !check if at least one element satisfies the condition
     isTrueMat = (COUNT(mat) .GT. 0)
  ELSE
     WRITE(*,*) "isTrueMat -- invalid log_op ", log_op
     STOP
  end IF

 END FUNCTION isTrueMat

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!! Complexify values and arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !Complexify a real array
 FUNCTION cmplxArr(tensA) result(tensO)

  REAL(KIND=DP),    INTENT(IN)  :: tensA(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: tensO(:)

  INTEGER :: i, dimA

  dimA = SIZE(tensA)

  ALLOCATE(tensO(dimA))
  DO i=1,dimA
     tensO(i) = CMPLX(tensA(i), 0.0d0, KIND=DP)
  end DO

 END FUNCTION cmplxArr


 !Complexify a real value
 FUNCTION cmplxVal(valIn) 

  REAL(KIND=DP),   INTENT(IN) :: valIn
  COMPLEX(KIND=DP)            :: cmplxVal

  cmplxVal = CMPLX(valIn, 0.0d0, KIND=DP)

 END FUNCTION cmplxVal

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Some simple mathematical functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !Find logN(x) for any log-base N
 FUNCTION logN_cmplx(x, n) result(logN)

  COMPLEX(KIND=DP), INTENT(IN) :: x
  REAL(KIND=DP),    INTENT(IN) :: n
  COMPLEX(KIND=DP)             :: logN

  logN = log(x)/log(n)

 END FUNCTION logN_cmplx


 FUNCTION logN_real(x, n) result(logN)

  REAL(KIND=DP), INTENT(IN) :: x
  REAL(KIND=DP), INTENT(IN) :: n
  REAL(KIND=DP)             :: logN

  logN = log(x)/log(n)

 END FUNCTION logN_real


 
 FUNCTION SQRTint(x) 

  INTEGER, INTENT(IN) :: x
  INTEGER             :: SQRTint

  SQRTint = NINT(SQRT(1.0D0 * x))

  IF(SQRTint * SQRTint .NE. x) THEN
     WRITE(*,*) "SQRTint -- SQRTint * SQRTint and x must match ", SQRTint * SQRTint, x
     STOP
  end IF

 END FUNCTION SQRTint


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Routines for swapping/reordering variables !!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Swap two integers !!!
 SUBROUTINE swapInt(varA, varB) 

  INTEGER, INTENT(INOUT) :: varA, varB
  INTEGER :: tempA

  !! temporary storage for varA
  tempA = varA
  
  !! swap varA, varB
  varA = varB
  varB = tempA

 END SUBROUTINE swapInt



 !!! Swap two complex vars !!!
 SUBROUTINE swapRe(varA, varB) 

  REAL(KIND=DP), INTENT(INOUT) :: varA, varB
  REAL(KIND=DP) :: tempA

  !! temporary storage for varA
  tempA = varA
  
  !! swap varA, varB
  varA = varB
  varB = tempA

 END SUBROUTINE swapRe




 !!! Swap two complex vars !!!
 SUBROUTINE swapCmplx(varA, varB) 

  COMPLEX(KIND=DP), INTENT(INOUT) :: varA, varB
  COMPLEX(KIND=DP) :: tempA

  !! temporary storage for varA
  tempA = varA
  
  !! swap varA, varB
  varA = varB
  varB = tempA

 END SUBROUTINE swapCmplx



 !!! Sort array using a simple bubble sort !!!
 SUBROUTINE sort_array(A)

  COMPLEX(KIND=DP), INTENT(INOUT) :: A(:)

  INTEGER :: i, imax(1), dimA

  dimA = SIZE(A)

  DO i=1,dimA
     imax = MAXLOC(ABS(A(i:dimA)))
     CALL swap(A(i), A(i-1+imax(1))) 
  end DO

 END SUBROUTINE sort_array

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE basic_functions
