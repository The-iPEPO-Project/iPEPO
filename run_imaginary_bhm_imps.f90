PROGRAM run_imaginary_bhm_imps

  USE utility
  USE datafile_utility
  USE definitions_mps_mpo

  USE simulation_parameters
  USE psi_sun_basis
  USE bhm_imaginary_propagator

  USE imaginary_observables_imps
  USE imaginary_time_evolution_imps

  IMPLICIT NONE

  !Psi state
  TYPE(block_mps),    ALLOCATABLE :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE :: iLambda(:)

  !Mpo propagator
  TYPE(block_mpo), ALLOCATABLE :: mpo_prop(:)

  !Simulation params
  TYPE(tebd_params) :: tebd
  REAL(KIND=DP)     :: eps, eta
  INTEGER           :: chi
  
  !Physical params
  REAL(KIND=DP)    :: g, U, J
  COMPLEX(KIND=DP) :: nval, g2val

  !Initialization, Reloading && Output Printing vars
  COMPLEX(KIND=DP), ALLOCATABLE :: psi_vec(:)
  REAL(KIND=DP), ALLOCATABLE :: s0(:)
  CHARACTER :: read_state_filename*256, read_state_str*1
  LOGICAL   :: read_state

  !Range of correlations (summary files)
  INTEGER :: min_sep, max_sep  

  !TEBD params - Main iteration
  COMPLEX(KIND=DP)   :: dt
  REAL(KIND=DP)      :: real_dt
  INTEGER            :: N_steps, test_interval
  INTEGER, PARAMETER :: test_min_sep = 1, test_max_sep = 4

  !Parameter loop vars
  REAL(KIND=DP) :: x, x0, dx
  INTEGER       :: nx, ix
  CHARACTER     :: xname*16

  !cpu times
  REAL(KIND=DP) :: cpu_times(4)

  WRITE(*,*) "Physical parameters: g, U, J "
  READ(*,*) g, U, J

  WRITE(*,*) "Hilbert space size"
  READ(*,*) basis_size

  WRITE(*,*) "N_steps, timestep, test_interval"
  READ(*,*) N_steps, real_dt, test_interval
  dt = CMPLX(real_dt, 0.0D0)

  WRITE(*,*) "range of correlations"
  READ(*,*) min_sep, max_sep

  WRITE(*,*) "TEBD and SVD precision settings"
  READ(*,*) eta, chi, eps

  WRITE(*,*) "Physical param to vary, Initial val, Increment, Number of vals"
  READ(*,*) xname, x0, dx, nx

  WRITE(*,*) "Reload existing state?"
  READ(*,*) read_state_str

  IF(SCAN(read_state_str,"yYtT") .GT. 0) THEN
     read_state=.TRUE.
     WRITE(*,*) "Reading state from file.  Filename?"
     READ(*,*) read_state_filename
  ELSE 
     read_state=.FALSE.
     ALLOCATE(s0(basis_size))
     WRITE(*,*) "Initial (product) state"
     READ(*,*) s0
  end IF

  !Set TEBD params
  CALL create_tebd_params(tebd, & 
                          & eps, chi, eta, & 
                          & dt, N_steps, test_interval, test_min_sep, test_max_sep)
  
  !Set up operators -- MUST be called before initializing psi state
  CALL setup_basic_operators()

  !PSI state: either reload, or initialize to a product state
  IF(read_state) THEN
     CALL initialize_psi_iMPS(iGamma, iLambda, filename=read_state_filename)
  ELSE
     !Prepare initial psi state
     CALL create_initial_state_vec(psi_vec, s0) 
     CALL initialize_psi_iMPS(iGamma, iLambda, psi_vec)
  end IF

  !Setup observables (operator lists and files)
  CALL setup_observables()

  !!!!!! Loop over a specified parameter !!!!!!!
  WRITE(*,*) "Starting parameter loop"
 
  DO ix=1,nx

     !Set the value of a variable, according to xname
     IF(nx .GT. 1) THEN
        x = x0 + (ix-1)*dx
     ELSE
        x = x0
     end IF

     !Decide which param to use
     CALL read_xname(x, xname, J, "J", U, "U")

     !Create propagator
     CALL setup_propagator(mpo_prop, g, U, J, dt)

     !Trotter time evolution to obtain NESS
     CALL time_evolve(iGamma, iLambda, mpo_prop, tebd)

     !Calculate summary observables and output SS (steady state)
     CALL calc_tebd_summary(iGamma, iLambda, x, min_sep, max_sep, g2_index=2, dump_parameters=dump_parameters)

  END DO

  !Close output files
  CALL close_datafiles()

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! NAME:   dump_parameters
 ! VARIABLES:
 !           HANDLE - File handle
 ! SYNOPSIS:
 ! Dump parameters to a file, or to the screen.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE dump_parameters(HANDLE)

    INTEGER, INTENT(IN) :: HANDLE
    
    WRITE(HANDLE, '("# local_dim: ",I5)'),  local_dim
    WRITE(HANDLE, '("# timestep : ", 2(D12.5,X))'),  dt
    WRITE(HANDLE, '("# g, U, J   : ", 3(E12.5,4X))'),  U, J

 END SUBROUTINE dump_parameters



 !!!!!!!!!!!!!!!!!!! Setup observables (operator lists and files) !!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE setup_observables()

  !Setup convergence files to record observables during time evolution
  ALLOCATE(TRACE_OUT_1P(4), TRACE_OP_1P(4, local_dim, local_dim))

  CALL add_1p_conv_trace(dn_op,     "a",  "time", "a",  "Norm")
  CALL add_1p_conv_trace(number_op, "n",  "time", "n",  "Norm")
  CALL add_1p_conv_trace(g2bare_op, "G2", "time", "G2", "Norm")
  CALL add_1p_conv_trace(eye,       "Id", "time", "Id", "Norm")

  ALLOCATE(TRACE_OUT_2P(3, test_min_sep : test_max_sep), TRACE_OP_2P(3, 2, local_dim, local_dim))

  CALL add_2p_conv_trace(up_op,     dn_op,     test_min_sep, test_max_sep, "ada", "time", "ada", "Norm")
  CALL add_2p_conv_trace(number_op, number_op, test_min_sep, test_max_sep, "ninj", "time", "ninj", "Norm")
  CALL add_2p_conv_trace(eye,       eye,       test_min_sep, test_max_sep, "IdId", "time", "IdId", "Norm")

  !get dump index -- DUMP_NUM = first file created after the convergence files
  CALL prepare_datafile(DUMP_NUM, "dump_ind", "Index of states for reloading")

  !Additionally setup summary files to record SS observables
  ALLOCATE(SUMM_OUT_1P(3), SUMM_OP_1P(3, local_dim, local_dim))

  CALL add_1p_summary(dn_op,     "a",  xname, "a",  "Norm")
  CALL add_1p_summary(number_op, "n",  xname, "n",  "Norm")
  CALL add_1p_summary(g2bare_op, "G2", xname, "G2", "Norm")

  ALLOCATE(SUMM_OUT_2P(2), SUMM_OP_2P(2, 2, local_dim, local_dim))

  CALL add_2p_summary(up_op,     dn_op,     "ada", xname, "ada", "Norm")
  CALL add_2p_summary(number_op, number_op, "nn",  xname, "nn",  "Norm")

  !Extra file for post-processed g2
  CALL setup_file_obs_2P(datafile=G2_SUMMARY, descriptor="g2_post_processed", param_str=xname, sepX_str="|i-j|", obs_str="g2", extra_str_1="nval")

 END SUBROUTINE setup_observables

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM run_imaginary_bhm_imps
