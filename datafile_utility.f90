MODULE datafile_utility

 USE basic_functions
 USE error_handling

 IMPLICIT NONE

 !! The following is for keeping track of which files we have opened.
 INTEGER, PARAMETER :: FIRST_HANDLE=300

 !! Initialise last to be one fewer, as equal values indicate a single file exists.
 INTEGER :: LAST_HANDLE=FIRST_HANDLE-1

 !! File handles for printing observables
 INTEGER, ALLOCATABLE :: TRACE_OUT_1P(:), TRACE_OUT_2P(:,:)
 INTEGER, ALLOCATABLE ::  SUMM_OUT_1P(:), SUMM_OUT_2P(:)

 !! Operators to take expectations of during trace
 COMPLEX(KIND=DP), ALLOCATABLE :: TRACE_OP_1P(:,:,:), TRACE_OP_2P(:,:,:,:)
 COMPLEX(KIND=DP), ALLOCATABLE ::  SUMM_OP_1P(:,:,:),  SUMM_OP_2P(:,:,:,:)

 INTEGER :: num_1p_trace=0,   num_2p_trace=0
 INTEGER :: num_1p_summary=0, num_2p_summary=0

CONTAINS

 !!!!! Setup convergence datafile !!!!!!
 SUBROUTINE setupConvDatafile(datafile, descr, parstr, errstr, valstr, xstr, xlog1, xlog2, chi, N_sites)

  INTEGER,          INTENT(INOUT)        :: datafile
  CHARACTER(LEN=*), INTENT(IN)           :: descr, valstr
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: parstr, errstr, xstr, xlog1, xlog2 !parstr = param str, xstr = extra val str, xlog = extra logical str
  INTEGER,          INTENT(IN), OPTIONAL :: chi, N_sites

  CHARACTER :: char_chi*16, char_N*16, char_err*16

  !! Setup data file (check if chi && N_sites variables are to be included in the filename)
  IF(PRESENT(chi) .AND. PRESENT(N_sites)) THEN 
      WRITE(char_chi, '(I3)') chi 
      WRITE(char_N,   '(I6)') N_sites
      CALL prepare_datafile(datafile, TRIM(ADJUSTL(descr))//"_N="//TRIM(ADJUSTL(char_N))//"_chi="//TRIM(ADJUSTL(char_chi)), descr)
  ELSE
      CALL prepare_datafile(datafile, TRIM(ADJUSTL(descr)), descr)
  end IF

  !! Set error str
  char_err = OptArg(errstr, "Err_"//TRIM(ADJUSTL(valstr)))

  !!! Write heading to file !!!

  !! Spacing from the top
  WRITE(datafile, '("#")')

  !! Check which strings are to be included in the file header (need to check this one by one...)
  IF(PRESENT(parstr)) THEN

      if(PRESENT(xstr) .AND. PRESENT(xlog1) .AND. PRESENT(xlog2)) then
                                 WRITE(datafile, '("# iter",3X,A12,10X,A12,5X,A12,24X,A12,4X,A12,4X,A12)') parstr, char_err, valstr, xstr, xlog1, xlog2
      elseif(PRESENT(xstr) .AND. PRESENT(xlog1)) then
                                 WRITE(datafile, '("# iter",3X,A12,10X,A12,5X,A12,24X,A12,4X,A12)') parstr, char_err, valstr, xstr, xlog1
      elseif(PRESENT(xlog1) .AND. PRESENT(xlog2)) then
                                 WRITE(datafile, '("# iter",3X,A12,10X,A12,5X,A12,24X,A12,4X,A12)') parstr, char_err, valstr, xlog1, xlog2
      elseif(PRESENT(xstr)) then
                                 WRITE(datafile, '("# iter",3X,A12,10X,A12,5X,A12,24X,A12)') parstr, char_err, valstr, xstr
      elseif(PRESENT(xlog1)) then
                                 WRITE(datafile, '("# iter",3X,A12,10X,A12,5X,A12,24X,A12)') parstr, char_err, valstr, xlog1
      else
                                 WRITE(datafile, '("# iter",3X,A12,10X,A12,5X,A12)') parstr, char_err, valstr
      endif
  ELSE
      if(PRESENT(xstr) .AND. PRESENT(xlog1) .AND. PRESENT(xlog2)) then
                                 WRITE(datafile, '("# iter",8X,A12,8X,A12,24X,A12,4X,A12,4X,A12)') char_err, valstr, xstr, xlog1, xlog2
      elseif(PRESENT(xstr) .AND. PRESENT(xlog1)) then
                                 WRITE(datafile, '("# iter",8X,A12,8X,A12,24X,A12,4X,A12)') char_err, valstr, xstr, xlog1
      elseif(PRESENT(xlog1) .AND. PRESENT(xlog2)) then
                                 WRITE(datafile, '("# iter",8X,A12,8X,A12,24X,A12,4X,A12)') char_err, valstr, xlog1, xlog2
      elseif(PRESENT(xstr)) then
                                 WRITE(datafile, '("# iter",8X,A12,8X,A12,24X,A12)') char_err, valstr, xstr
      elseif(PRESENT(xlog1)) then
                                 WRITE(datafile, '("# iter",8X,A12,8X,A12,24X,A12)') char_err, valstr, xlog1
      else
                                 WRITE(datafile, '("# iter",8X,A12,8X,A12)') char_err, valstr
      endif
  END IF

  !! Spacing from the bottom
  WRITE(datafile, '("#")')

 END SUBROUTINE setupConvDatafile





 !!! Setup file for writing 1P observables !!!
 SUBROUTINE setup_file_obs_1P(datafile, descriptor, param_str, obs_str, extra_str_1, extra_str_2)

  INTEGER,          INTENT(INOUT)        :: datafile
  CHARACTER(LEN=*), INTENT(IN)           :: descriptor, param_str, obs_str
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extra_str_1, extra_str_2

  !! Setup datafile
  CALL prepare_datafile(datafile, TRIM(ADJUSTL(descriptor)), descriptor)

  !! Write heading to file
  IF(PRESENT(extra_str_1) .AND. PRESENT(extra_str_2)) THEN

     WRITE(datafile, '("#",A12,8X,A12,24X,A12,8X,A12)') param_str, obs_str, extra_str_1, extra_str_2
     WRITE(datafile, '("#")')

  ELSEIF(PRESENT(extra_str_1)) THEN

     WRITE(datafile, '("#",A12,8X,A12,24X,A12)') param_str, obs_str, extra_str_1
     WRITE(datafile, '("#")')
  ELSE

     WRITE(datafile, '("#",A12,10X,A12)') param_str, obs_str
     WRITE(datafile, '("#")')
  end IF
       
 END SUBROUTINE setup_file_obs_1P








 !!! Setup file for writing 2P observables !!!
 SUBROUTINE setup_file_obs_2P(datafile, descriptor, param_str, sepX_str, sepY_str, obs_str, extra_str_1, extra_str_2)

  INTEGER,          INTENT(INOUT)        :: datafile
  CHARACTER(LEN=*), INTENT(IN)           :: descriptor, param_str, sepX_str, obs_str
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: sepY_str, extra_str_1, extra_str_2

  !! Setup datafile
  CALL prepare_datafile(datafile, TRIM(ADJUSTL(descriptor)), descriptor)

  !! Write heading to file
  IF(PRESENT(sepY_str) .AND. PRESENT(extra_str_1) .AND. PRESENT(extra_str_2)) THEN

     WRITE(datafile, '("#",A12,8X,A12,8X,A12,8X,A12,24X,A12,8X,A12)') param_str, sepX_str, sepY_str, obs_str, extra_str_1, extra_str_2
     WRITE(datafile, '("#")')

  ELSEIF(PRESENT(sepY_str) .AND. PRESENT(extra_str_1)) THEN

     WRITE(datafile, '("#",A12,8X,A12,8X,A12,8X,A12,24X,A12)') param_str, sepX_str, sepY_str, obs_str, extra_str_1
     WRITE(datafile, '("#")')

  ELSEIF(PRESENT(extra_str_1) .AND. PRESENT(extra_str_2)) THEN

     WRITE(datafile, '("#",A12,8X,A12,8X,A12,24X,A12,8X,A12)') param_str, sepX_str, obs_str, extra_str_1, extra_str_2
     WRITE(datafile, '("#")')

  ELSEIF(PRESENT(sepY_str)) THEN

     WRITE(datafile, '("#",A12,8X,A12,8X,A12,8X,A12)') param_str, sepX_str, sepY_str, obs_str
     WRITE(datafile, '("#")')

  ELSEIF(PRESENT(extra_str_1)) THEN

     WRITE(datafile, '("#",A12,8X,A12,8X,A12,24X,A12)') param_str, sepX_str, obs_str, extra_str_1
     WRITE(datafile, '("#")')
  ELSE

     WRITE(datafile, '("#",A12,8X,A12,10X,A12)') param_str, sepX_str, obs_str
     WRITE(datafile, '("#")')
  end IF
       
 END SUBROUTINE setup_file_obs_2P







 !!!!!!!! Write convergence data !!!!!!!!!!!
 SUBROUTINE write_conv_data(datafile, iter, ipar, errval, val, xval, xlog1, xlog2)

  INTEGER,          INTENT(IN)           :: datafile
  INTEGER,          INTENT(IN)           :: iter
  REAL(KIND=DP),    INTENT(IN), OPTIONAL :: ipar
  REAL(KIND=DP),    INTENT(IN)           :: errval
  COMPLEX(KIND=DP), INTENT(IN)           :: val
  COMPLEX(KIND=DP), INTENT(IN), OPTIONAL :: xval          !xval = extra val
  LOGICAL,          INTENT(IN), OPTIONAL :: xlog1, xlog2  !xlog = extra logical var

  !! Write output to file
  IF(PRESENT(ipar)) THEN

      if(PRESENT(xval) .AND. PRESENT(xlog1) .AND. PRESENT(xlog2)) then
                                  WRITE(datafile, '(I7,4X,D12.5,8X,D12.5,2X,2(D12.5,X),5X,2(D12.5,X),5X,L3,5X,L3)') iter, ipar, errval, val, xval, xlog1, xlog2
      elseif(PRESENT(xval) .AND. PRESENT(xlog1)) then
                                  WRITE(datafile, '(I7,4X,D12.5,8X,D12.5,2X,2(D12.5,X),5X,2(D12.5,X),5X,L3)') iter, ipar, errval, val, xval, xlog1
      elseif(PRESENT(xlog1) .AND. PRESENT(xlog2)) then
                                  WRITE(datafile, '(I7,4X,D12.5,8X,D12.5,2X,2(D12.5,X),5X,L3,5X,L3)') iter, ipar, errval, val, xval, xlog1, xlog2
      elseif(PRESENT(xval)) then
                                  WRITE(datafile, '(I7,4X,D12.5,8X,D12.5,2X,2(D12.5,X),5X,2(D12.5,X))') iter, ipar, errval, val, xval
      elseif(PRESENT(xlog1)) then
                                  WRITE(datafile, '(I7,4X,D12.5,8X,D12.5,2X,2(D12.5,X),5X,L3)') iter, ipar, errval, val, xlog1
      else
                                  WRITE(datafile, '(I7,4X,D12.5,8X,D12.5,2X,2(D12.5,X))') iter, ipar, errval, val
      endif

  ELSE

      if(PRESENT(xval) .AND. PRESENT(xlog1) .AND. PRESENT(xlog2)) then
                                  WRITE(datafile, '(I7,8X,D12.5,2X,2(D12.5,X),5X,2(D12.5,X),5X,L3,5X,L3)') iter, errval, val, xval, xlog1, xlog2
      elseif(PRESENT(xval) .AND. PRESENT(xlog1)) then
                                  WRITE(datafile, '(I7,8X,D12.5,2X,2(D12.5,X),5X,2(D12.5,X),5X,L3)') iter, errval, val, xval, xlog1
      elseif(PRESENT(xlog1) .AND. PRESENT(xlog2)) then
                                  WRITE(datafile, '(I7,8X,D12.5,2X,2(D12.5,X),5X,L3,5X,L3)') iter, errval, val, xval, xlog1, xlog2
      elseif(PRESENT(xval)) then
                                  WRITE(datafile, '(I7,8X,D12.5,2X,2(D12.5,X),5X,2(D12.5,X))') iter, errval, val, xval
      elseif(PRESENT(xlog1)) then
                                  WRITE(datafile, '(I7,8X,D12.5,2X,2(D12.5,X),5X,L3)') iter, errval, val, xlog1
      else
                                  WRITE(datafile, '(I7,8X,D12.5,2X,2(D12.5,X))') iter, errval, val
      endif
  END IF

 END SUBROUTINE write_conv_data



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   add_1p_trace
  ! VARIABLES:
  !           prefix     - Filename prefix
  !           op         - Operator to store
  !           descriptor - File descriptor
  ! SYNOPSIS:
  ! Add a single site expectation to the list of operators.  This routine
  ! ensures that filenames stay correlated with the operators measured.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE add_1p_conv_trace(op, descriptor, param_str, obs_str, extra_str1)

    COMPLEX(KIND=DP), INTENT(IN)           :: op(:,:)
    CHARACTER(LEN=*), INTENT(IN)           :: descriptor, param_str, obs_str
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extra_str1

    CHARACTER :: full_descriptor*48

    !! (1) Increment the file counter, and check we do not overrun
    num_1p_trace = num_1p_trace + 1
    
    IF(num_1p_trace .GT. SIZE(TRACE_OUT_1P)) THEN
       WRITE(*,*) "Not enough one point files allocated"
       STOP
    end IF

    !! Create a full descriptor for the new datafile/observable
    full_descriptor = "exp_1P_"//TRIM(ADJUSTL(descriptor))

    !! (2) Setup new file, record the file handle it gives, write file heading to file
    IF(PRESENT(extra_str1)) THEN
       CALL setupConvDatafile(datafile=TRACE_OUT_1P(num_1p_trace), descr=full_descriptor, parstr=param_str, valstr=obs_str, xstr=extra_str1)
    ELSE
       CALL setupConvDatafile(datafile=TRACE_OUT_1P(num_1p_trace), descr=full_descriptor, parstr=param_str, valstr=obs_str)
    end IF

    !! (3) Record the operator
    TRACE_OP_1P(num_1p_trace, :, :) = op
       
  END SUBROUTINE add_1p_conv_trace



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   add_2p_trace
  ! VARIABLES:
  !           prefix     - Filename prefix
  !           op1, op2   - Operators to store
  !           descriptor - FIle descriptor
  ! SYNOPSIS:
  ! Add a two site expectation to the list of operators.  This routine
  ! ensures that filenames stay correlated with the operators measured.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE add_2p_conv_trace(op2, op1, min_sep, max_sep, descriptor, param_str, obs_str, extra_str1)

    COMPLEX(KIND=DP), INTENT(IN)           :: op2(:,:), op1(:,:) !NB input in reverse order: we'll calc TRACE(op2*op1*RHO) but store as {op1,op2}
    INTEGER,          INTENT(IN)           :: min_sep, max_sep
    CHARACTER(LEN=*), INTENT(IN)           :: descriptor, param_str, obs_str
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extra_str1

    INTEGER :: sep
    CHARACTER :: full_descriptor*48, sep_str*4

    !! (1) Increment the  counter, and check we do not overrun
    num_2p_trace = num_2p_trace + 1
    
    IF(num_2p_trace .GT. SIZE(TRACE_OUT_2P,1)) THEN
       WRITE(*,*) "Not enough two point files allocated"
       STOP
    end IF

    !! (2) Create trace file for each sep
    DO sep = min_sep, max_sep

       !! Create a full descriptor for each separation
       WRITE(sep_str,'(I4)') sep
       full_descriptor = "corr_2P_"//TRIM(ADJUSTL(descriptor))//"_sep="//TRIM(ADJUSTL(sep_str))

       !! Setup new file, record the file handle it gives, write file heading to file
       IF(PRESENT(extra_str1)) THEN
          CALL setupConvDatafile(datafile=TRACE_OUT_2P(num_2p_trace, sep), descr=full_descriptor, parstr=param_str, valstr=obs_str, xstr=extra_str1)
       ELSE
          CALL setupConvDatafile(datafile=TRACE_OUT_2P(num_2p_trace, sep), descr=full_descriptor, parstr=param_str, valstr=obs_str)
       end IF
    end DO

    !! (3) Record the operators
    TRACE_OP_2P(num_2p_trace, 1, :, :) = op1
    TRACE_OP_2P(num_2p_trace, 2, :, :) = op2
    
  END SUBROUTINE add_2p_conv_trace




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   add_1p_summary
  ! VARIABLES:
  !           suffix     - Filename suffix
  !           op         - Operator to store
  !           descriptor - FIle descriptor
  !           xstr       - Name of what is varying
  !           prefix     - (Optional) : Prefix for filenames
  ! SYNOPSIS:
  ! Add a single site expectation to the list of operators.  This routine
  ! ensures that filenames stay correlated with the operators measured.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE add_1p_summary(op, descriptor, param_str, obs_str, extra_str1)

    COMPLEX(KIND=DP), INTENT(IN)           :: op(:,:)
    CHARACTER(LEN=*), INTENT(IN)           :: descriptor, param_str, obs_str
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extra_str1

    !! Increment the file counter, and check we do not overrun
    num_1p_summary = num_1p_summary + 1
    
    !! Open new file, and record the file handle it gives. Note augmentation of strings
    IF(num_1p_summary .GT. SIZE(SUMM_OUT_1P)) THEN
       WRITE(*,*) "Not enough one point files allocated"
       STOP
    end IF

    !! Open new file, and record the file handle it gives. Note the augmentation of strings
    CALL prepare_datafile(SUMM_OUT_1P(num_1p_summary), "summary_exp_1P_"//TRIM(ADJUSTL(descriptor)), TRIM(ADJUSTL(descriptor)))

    !! Write heading to file
    IF(PRESENT(extra_str1)) THEN

       WRITE(SUMM_OUT_1P(num_1p_summary), '("#",A12,8X,A12,24X,A12)') param_str, obs_str, extra_str1
       WRITE(SUMM_OUT_1P(num_1p_summary), '("#")')
    ELSE

       WRITE(SUMM_OUT_1P(num_1p_summary), '("#",A12,10X,A12)') param_str, obs_str
       WRITE(SUMM_OUT_1P(num_1p_summary), '("#")')
    end IF

    !! Record the operator
    SUMM_OP_1P(num_1p_summary, :, :) = op
       
  END SUBROUTINE add_1p_summary




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   add_2p_summary
  ! VARIABLES:
  !           suffix     - Filename suffix
  !           op1, op2   - Operators to store
  !           descriptor - FIle descriptor
  !           xstr       - Name of what is varying
  !           prefix     - (Optional) : Prefix for filenames
  ! SYNOPSIS:
  ! Add a two site expectation to the list of operators.  This routine
  ! ensures that filenames stay correlated with the operators measured.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE add_2p_summary(op1, op2, descriptor, param_str, sep_str, obs_str, extra_str1) !! FIXME -- swap op1, op2 order!!!

    COMPLEX(KIND=DP), INTENT(IN)           :: op1(:,:), op2(:,:)
    CHARACTER(LEN=*), INTENT(IN)           :: descriptor, param_str, sep_str, obs_str
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: extra_str1

    !! Increment the  counter, and check we do not overrun
    num_2p_summary = num_2p_summary + 1
    
    IF(num_2p_summary .GT. SIZE(SUMM_OUT_2P)) THEN
       WRITE(*,*) "Not enough two point files allocated"
       STOP
    end IF

    !! Open new file, and record the file handle it gives. Note the augmentation of strings
    CALL prepare_datafile(SUMM_OUT_2P(num_2p_summary), "summary_corr_2P_"//TRIM(ADJUSTL(descriptor)), TRIM(ADJUSTL(descriptor)))

    !! Write heading to each file
    IF(PRESENT(extra_str1)) THEN

       WRITE(SUMM_OUT_2P(num_2p_summary), '("#",A12,8X,A12,24X,A12,8X,A12)') param_str, sep_str, obs_str, extra_str1
       WRITE(SUMM_OUT_2P(num_2p_summary), '("#")')
    ELSE

       WRITE(SUMM_OUT_2P(num_2p_summary), '("#",A12,8X,A12,24X,A12)') param_str, sep_str, obs_str
       WRITE(SUMM_OUT_2P(num_2p_summary), '("#")')
    end IF

    !! Record the operators
    SUMM_OP_2P(num_2p_summary, 1, :, :) = op1
    SUMM_OP_2P(num_2p_summary, 2, :, :) = op2
    
  END SUBROUTINE add_2p_summary



 !!! Write data to 1P file !!!
 SUBROUTINE write_summary_data_1P(datafile, param, expval, xval1, xval2)

  INTEGER,          INTENT(IN)           :: datafile
  REAL(KIND=DP),    INTENT(IN)           :: param
  COMPLEX(KIND=DP), INTENT(IN)           :: expval
  COMPLEX(KIND=DP), INTENT(IN), OPTIONAL :: xval1, xval2
  
  !! Write output to file
  IF(PRESENT(xval1) .AND. PRESENT(xval2)) THEN

     WRITE(datafile, '(D12.5,5X,2(D12.5,X),5X,2(D12.5,X),5X,2(D12.5,X))') param, expval, xval1, xval2

  ELSEIF(PRESENT(xval1)) THEN

     WRITE(datafile, '(D12.5,5X,2(D12.5,X),5X,2(D12.5,X))') param, expval, xval1
  ELSE

     WRITE(datafile, '(D12.5,5X,2(D12.5,X))') param, expval
  end IF

 END SUBROUTINE write_summary_data_1P



 !!! Write data to 2P file !!!
 SUBROUTINE write_summary_data_2P(datafile, param, sep, expval, xval1, xval2)

  INTEGER,          INTENT(IN)           :: datafile
  REAL(KIND=DP),    INTENT(IN)           :: param
  INTEGER,          INTENT(IN)           :: sep
  COMPLEX(KIND=DP), INTENT(IN)           :: expval
  COMPLEX(KIND=DP), INTENT(IN), OPTIONAL :: xval1, xval2

  !! Write output to file
  IF(PRESENT(xval1) .AND. PRESENT(xval2)) THEN

     WRITE(datafile, '(D12.5,5X,I5,5X,2(D12.5,X),5X,2(D12.5,X),5X,2(D12.5,X))') param, sep, expval, xval1, xval2

  ELSEIF(PRESENT(xval1)) THEN

     WRITE(datafile, '(D12.5,5X,I5,5X,2(D12.5,X),5X,2(D12.5,X))') param, sep, expval, xval1
  ELSE

     WRITE(datafile, '(D12.5,5X,I5,5X,2(D12.5,X))') param, sep, expval
  end IF

 END SUBROUTINE write_summary_data_2P



 !!! Set up a summary file in SUMMARY_FILES(:) array for recording dynamical expectation values at different times & separations
 !!! The definition of t_steps (as a function of t_steps_in) assumes we're using TraMTE for calculating dynamical observables
 SUBROUTINE setup_files_dynamics(SUMMARY_FILES, descriptor, tsteps_raw, max_sep, dump_params)

  INTERFACE
     SUBROUTINE dump_params(HANDLE)
         INTEGER, INTENT(IN) :: HANDLE
     END SUBROUTINE dump_params
  end INTERFACE

  INTEGER, ALLOCATABLE, INTENT(INOUT) :: SUMMARY_FILES(:)
  CHARACTER(LEN=*),     INTENT(IN)    :: descriptor
  INTEGER,              INTENT(IN)    :: tsteps_raw, max_sep

  INTEGER :: sep, tsteps
  CHARACTER :: char_sep*16, char_tsteps*16, summary_filename*32, full_descriptor*32

  !! Allocate expval files
  ALLOCATE(SUMMARY_FILES(max_sep + 1))

  !! Number of t-steps
  tsteps = (tsteps_raw - 1)/2
  WRITE(char_tsteps,'(I6)') tsteps 

  !! Open a new file, and record its file handle (a number corresponding to the file) in the array SUMMARY_FILES(sep+1)
  DO sep=0,max_sep   
   
      !! Setup a new file for each different separation
      WRITE(char_sep,'(I5)') sep

      IF(tsteps .LT. 10) THEN 

           summary_filename="summary_"//TRIM(ADJUSTL(descriptor))//"_sep="//TRIM(ADJUSTL(char_sep))//"_t=000"//TRIM(ADJUSTL(char_tsteps))

      ELSEIF((tsteps .LT. 100) .AND. (tsteps .GE. 10)) THEN 

           summary_filename="summary_"//TRIM(ADJUSTL(descriptor))//"_sep="//TRIM(ADJUSTL(char_sep))//"_t=00"//TRIM(ADJUSTL(char_tsteps))

      ELSEIF((tsteps .LT. 1000) .AND. (tsteps .GE. 100)) THEN 

           summary_filename="summary_"//TRIM(ADJUSTL(descriptor))//"_sep="//TRIM(ADJUSTL(char_sep))//"_t=0"//TRIM(ADJUSTL(char_tsteps))
      ELSE
           summary_filename="summary_"//TRIM(ADJUSTL(descriptor))//"_sep="//TRIM(ADJUSTL(char_sep))//"_t="//TRIM(ADJUSTL(char_tsteps))
      end IF

      !! full descriptor for each file
      full_descriptor=TRIM(ADJUSTL(descriptor))//", |i-j| = "//TRIM(ADJUSTL(char_sep))

      !! prepare datafile
      CALL prepare_datafile(SUMMARY_FILES(sep+1), summary_filename, full_descriptor)
  end DO

  !! print params to all files
  CALL print_params(dump_params)

  !! print headings to all expval files
  DO sep=0,max_sep
      WRITE(SUMMARY_FILES(sep+1), '("# t_steps",5X,A12,20X,A12)') "expval", "trace"
      WRITE(SUMMARY_FILES(sep+1), '("#")')
  end DO
 
 END SUBROUTINE setup_files_dynamics



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BASIC ROUTINES FOR MANAGING DATAFILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Prepare a data file: sets up && opens a file for writing data
 !!! HANDLE = File handle to open, prefix = Prefix for filename, descriptor = Descriptor to write in file
 SUBROUTINE prepare_datafile(HANDLE, prefix, descriptor)

  INTEGER,          INTENT(OUT) :: HANDLE
  CHARACTER(LEN=*), INTENT(IN)  :: prefix, descriptor

  CHARACTER :: filename*64
  LOGICAL :: file_exists
  
  !! Increment file handle
  LAST_HANDLE=LAST_HANDLE+1

  !! Set the file handle to return.
  HANDLE=LAST_HANDLE

  !! Open the file
  filename=TRIM(ADJUSTL(prefix))//".txt"
  INQUIRE(FILE=filename, EXIST=file_exists) 

  IF(file_exists) THEN
      !! If file exists, append data to the end of the file
      OPEN(HANDLE, FILE=filename, status='old', action='write', form='formatted', position="append")
  ELSE
      !! Otherwise create a new file, write descriptor to the new file
      OPEN(HANDLE, FILE=filename)
      WRITE(HANDLE,'("#",A64)') ADJUSTL(descriptor)
  end IF

 END SUBROUTINE prepare_datafile



 !!! Write params to datafiles !!!
 SUBROUTINE print_params(dump_params)

  INTERFACE
     SUBROUTINE dump_params(HANDLE)
         INTEGER, INTENT(IN) :: HANDLE
     END SUBROUTINE dump_params
  end INTERFACE

  INTEGER :: file_handle

  CALL dump_params(6)
  DO file_handle=FIRST_HANDLE, LAST_HANDLE
     CALL dump_params(file_handle)
  end DO
 
 END SUBROUTINE print_params



 !!! Close datafiles !!!
 SUBROUTINE close_datafiles()

  INTEGER :: file_handle

  !! Closing all files
  DO file_handle=FIRST_HANDLE, LAST_HANDLE
      CLOSE(file_handle)
  end DO

 END SUBROUTINE close_datafiles

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE datafile_utility
