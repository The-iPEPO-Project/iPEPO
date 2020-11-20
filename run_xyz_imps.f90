PROGRAM run_xyz_imps

  USE utility
  USE datafile_utility
  USE definitions_mps_mpo

  USE simulation_parameters
  USE xyz_propagator_pauli_basis

  USE observables_imps
  USE time_evolution_imps

  IMPLICIT NONE

  !Rho state
  TYPE(block_mps),    ALLOCATABLE :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE :: iLambda(:)

  !Mpo propagator
  TYPE(block_mpo), ALLOCATABLE :: mpo_prop(:)

  !Simulation params
  TYPE(tebd_params) :: tebd
  REAL(KIND=DP)     :: eps, eta
  INTEGER           :: chi
  
  !Physical params
  REAL(KIND=DP) :: kappa, Jx, Jy, Jz

  !Initialization, Reloading && Output Printing vars
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_vec(:)
  REAL(KIND=DP) :: s0(3)
  CHARACTER :: read_state_filename*256, read_state_str*1
  LOGICAL   :: read_state

  !Range of correlations (summary files)
  INTEGER :: min_sep, max_sep  

  !TEBD params - Main iteration
  COMPLEX(KIND=DP)   :: dt
  REAL(KIND=DP)      :: real_dt
  INTEGER            :: N_steps, test_interval
  INTEGER, PARAMETER :: test_min_sep = 1, test_max_sep = 2

  !Parameter loop vars
  REAL(KIND=DP) :: x, x0, dx
  INTEGER       :: nx, ix
  CHARACTER     :: xname*16

  !cpu times
  REAL(KIND=DP) :: cpu_times(4)

  WRITE(*,*) "kappa, Jx, Jy, Jz"
  READ(*,*) kappa, Jx, Jy, Jz

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
     WRITE(*,*) "Initial (product) state, sx,sy,sz"
     READ(*,*) s0
  end IF

  CALL cpu_time(cpu_times(1))

  !Set TEBD params
  CALL create_tebd_params(tebd, & 
                          & eps, chi, eta, & 
                          & dt, N_steps, test_interval, test_min_sep, test_max_sep)
  
  !Set up operators -- MUST be called before initializing rho state
  CALL setup_basic_operators()

  !RHO state: either reload, or initialize to a product state
  IF(read_state) THEN
     CALL initialize_rho_iMPS(iGamma, iLambda, filename=read_state_filename)
  ELSE
     CALL create_initial_state_vec(rho_vec, s0) 
     CALL initialize_rho_iMPS(iGamma, iLambda, rho_vec)
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
     CALL read_xname(x, xname, Jx, "Jx", Jy, "Jy", Jz, "Jz")

     !Create propagator
     CALL setup_propagator(mpo_prop, Jx, Jy, Jz, kappa, dt, 'MPS')

     !Trotter time evolution to obtain NESS
     CALL time_evolve(iGamma, iLambda, mpo_prop, tebd)

     !Calculate summary observables and output SS (steady state)
     CALL calc_tebd_summary(iGamma, iLambda, x, min_sep, max_sep, dump_parameters=dump_parameters)

  END DO

  CALL cpu_time(cpu_times(2))

  !Simulation time
  WRITE(*,*) "TEBD time = ", cpu_times(2) - cpu_times(1)

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

    USE simulation_parameters

    INTEGER, INTENT(IN) :: HANDLE
    
    WRITE(HANDLE,*)
    WRITE(HANDLE,*)
    WRITE(HANDLE, '("# local_dim: ", I5)') local_dim
    WRITE(HANDLE, '("# timestep : ", 2(D12.5,X))'),  dt
    WRITE(HANDLE, '("# kappa= ", D12.5, " J=", 3(D12.5,X))'), kappa, Jx, Jy, Jz

 END SUBROUTINE dump_parameters



 !!!!!!!!!!!!!!!!!!! Setup observables (operator lists and files) !!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE setup_observables()

  !Setup convergence files to record observables during time evolution
  ALLOCATE(TRACE_OUT_1P(4), TRACE_OP_1P(4, local_dim, local_dim))

  CALL add_1p_conv_trace(sigma_x, "Sx", "time", "Sx", "Trace")
  CALL add_1p_conv_trace(sigma_y, "Sy", "time", "Sy", "Trace")
  CALL add_1p_conv_trace(sigma_z, "Sz", "time", "Sz", "Trace")
  CALL add_1p_conv_trace(eye,     "Id", "time", "Id", "Trace")

  ALLOCATE(TRACE_OUT_2P(5, test_min_sep : test_max_sep), TRACE_OP_2P(5, 2, local_dim, local_dim))

  CALL add_2p_conv_trace(sigma_x, sigma_x, test_min_sep, test_max_sep, "SxSx", "time", "SxSx", "Trace")
  CALL add_2p_conv_trace(sigma_y, sigma_y, test_min_sep, test_max_sep, "SySy", "time", "SySy", "Trace")
  CALL add_2p_conv_trace(sigma_x, sigma_y, test_min_sep, test_max_sep, "SxSy", "time", "SxSy", "Trace")
  CALL add_2p_conv_trace(sigma_y, sigma_x, test_min_sep, test_max_sep, "SySx", "time", "SySx", "Trace")
  CALL add_2p_conv_trace(eye,     eye,     test_min_sep, test_max_sep, "IdId", "time", "IdId", "Trace")

  !get dump index -- DUMP_NUM = first file created after the convergence files
  CALL prepare_datafile(DUMP_NUM, "dump_ind", "Index of states for reloading")

  !Additionally setup summary files to record SS observables
  ALLOCATE(SUMM_OUT_1P(3), SUMM_OP_1P(3, local_dim, local_dim))

  CALL add_1p_summary(sigma_x, "Sx", xname, "Sx", "Trace")
  CALL add_1p_summary(sigma_y, "Sy", xname, "Sy", "Trace")
  CALL add_1p_summary(sigma_z, "Sz", xname, "Sz", "Trace")

  ALLOCATE(SUMM_OUT_2P(4), SUMM_OP_2P(4, 2, local_dim, local_dim))

  CALL add_2p_summary(sigma_x, sigma_x, "SxSx", xname, "SxSx", "Trace")
  CALL add_2p_summary(sigma_y, sigma_y, "SySy", xname, "SySy", "Trace")
  CALL add_2p_summary(sigma_x, sigma_y, "SxSy", xname, "SxSy", "Trace")
  CALL add_2p_summary(sigma_y, sigma_x, "SySx", xname, "SySx", "Trace")

 END SUBROUTINE setup_observables

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM run_xyz_imps
