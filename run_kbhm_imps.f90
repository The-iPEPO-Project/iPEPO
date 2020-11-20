PROGRAM run_kbhm_imps

  USE utility
  USE datafile_utility
  USE definitions_mps_mpo

  USE simulation_parameters
  USE sun_basis
  USE kbhm_propagator

  USE observables_imps
  USE time_evolution_imps

  IMPLICIT NONE

  !! Rho state
  TYPE(block_mps),    ALLOCATABLE :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE :: iLambda(:)

  !! Mpo propagator
  TYPE(block_mpo),    ALLOCATABLE :: mpo_prop(:)

  !! Simulation params
  TYPE(tebd_params) :: tebd
  REAL(KIND=DP)     :: eps, eta, xi
  INTEGER           :: chi
  
  !! Physical params
  REAL(KIND=DP)    :: delta, U, f, kappa, modJ, argJ
  COMPLEX(KIND=DP) :: J
  COMPLEX(KIND=DP) :: nval, g2val

  !! Initialization, Reloading && Output Printing vars
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_vec(:), rho_mat(:,:)
  CHARACTER :: read_state_filename*256, read_state_str*1
  LOGICAL   :: read_state

  !! Range of correlations (summary files)
  INTEGER :: min_sep, max_sep  

  !! TEBD params - Main iteration
  INTEGER, PARAMETER :: N_tebd = 5
  INTEGER            :: ll, ll0, HANDLE 

  COMPLEX(KIND=DP)   :: dt(N_tebd)
  REAL(KIND=DP)      :: real_dt(N_tebd)
  INTEGER            :: N_steps(N_tebd), test_interval(N_tebd)

  INTEGER, PARAMETER :: test_min_sep = 1, test_max_sep = 4
  INTEGER, PARAMETER :: exit_level(N_tebd) = (/1, 2, 2, 2, 2/)

  !! Parameter loop vars
  REAL(KIND=DP) :: x, x0, dx
  INTEGER       :: nx, ix
  CHARACTER     :: xname*16

  !! cpu times
  REAL(KIND=DP) :: cpu_times(4)

  WRITE(*,*) "delta, U, f, kappa, |J|, arg(J) = k "
  READ(*,*) delta, U, f, kappa, modJ, argJ

  WRITE(*,*) "Hilbert space size"
  READ(*,*) basis_size

  WRITE(*,*) "TEBD -- N_steps"
  READ(*,*) N_steps

  WRITE(*,*) "TEBD -- timestep"
  READ(*,*) real_dt
  dt = CMPLX(real_dt, 0.0D0)
  WRITE(*,*) "Complex dt = ", dt

  WRITE(*,*) "TEBD -- test_interval"
  READ(*,*) test_interval

  WRITE(*,*) "range of correlations"
  READ(*,*) min_sep, max_sep

  WRITE(*,*) "TEBD and SVD precision settings"
  READ(*,*) eta, chi, eps

  WRITE(*,*) "Propagator precision"
  READ(*,*) xi

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
  end IF
  
  !! Set up operators -- MUST be called before initializing rho state
  CALL setup_basic_operators()

  !! RHO state: either reload, or initialize to a product state
  IF(read_state) THEN
     CALL initialize_rho_iMPS(iGamma, iLambda, filename=read_state_filename)
  ELSE
     !! Prepare initial rho matrix: rep product state as supervector (start with an empty state)
     ALLOCATE(rho_mat(basis_size, basis_size)); rho_mat=0.0D0; rho_mat(1,1) = 1.0D0

     CALL create_initial_state_vec(rho_vec, rho_mat, represent) 
     CALL initialize_rho_iMPS(iGamma, iLambda, rho_vec=rho_vec) !, chi=1)
  end IF

  !! Setup observables (operator lists and files)
  CALL setup_observables()

  !!!!!! Loop over a specified parameter !!!!!!!
  WRITE(*,*) "Starting parameter loop"
 
  DO ix=1,nx

     !! Set the value of a variable, according to xname
     IF(nx .GT. 1) THEN
        x = x0 + (ix-1)*dx
     ELSE
        x = x0
     end IF

     !! Decide which param to use
     CALL read_xname(x, xname, f, "f", modJ, "J", argJ, "k")

     !! Combine modulus and phase to get j
     J=modJ*exp(ii*argJ)

     !! REAL-TIME TEBD propagation using decreasing values of timestep dt 
     DO ll=1,1 !N_tebd

        !! Create a new instance of TEBD
        CALL create_TEBD_params(TEBD, & 
                                & eps, chi, eta, & 
                                & dt(ll), N_steps(ll), test_interval(ll), test_min_sep, test_max_sep, exit_level(ll))
        !! Create propagator
        CALL setup_propagator(mpo_prop, delta, f, U, J, kappa, 0.5D0*TEBD % dt, xi=xi)

        !! Trotter time evolution to obtain NESS
        CALL time_evolve(iGamma, iLambda, mpo_prop, TEBD)

        !! Calculate summary observables and output SS (steady state)
        CALL calc_tebd_summary(iGamma, iLambda, x, min_sep, max_sep, g2_index=2, dump_parameters=dump_parameters)

        !! Add linebreak to datafiles 
        DO HANDLE=FIRST_HANDLE, DUMP_NUM
           WRITE(HANDLE, *)
           WRITE(HANDLE, *)
        end DO
     END DO
  END DO

  !! Close output files
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
    WRITE(HANDLE, '("# hopping  : ", 2(D12.5,X))'),  J
    WRITE(HANDLE, '("# delta,U  : ", 2(D12.5,4X))'),  delta,U
    WRITE(HANDLE, '("# f,kappa  : ", 2(D12.5,4X))'),  f,kappa

 END SUBROUTINE dump_parameters



 !!!!!!!!!!!!!!!!!!! Setup observables (operator lists and files) !!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE setup_observables()

   !Setup convergence files to record observables during time evolution
   ALLOCATE(TRACE_OUT_1P(4), TRACE_OP_1P(4, local_dim, local_dim))

   CALL add_1p_conv_trace(dn_op,     "a",  "time", "a",  "Trace")
   CALL add_1p_conv_trace(number_op, "n",  "time", "n",  "Trace")
   CALL add_1p_conv_trace(g2bare_op, "G2", "time", "G2", "Trace")
   CALL add_1p_conv_trace(eye,       "Id", "time", "Id", "Trace")

   ALLOCATE(TRACE_OUT_2P(3, test_min_sep : test_max_sep), TRACE_OP_2P(3, 2, local_dim, local_dim))

   CALL add_2p_conv_trace(up_op,     dn_op,     test_min_sep, test_max_sep, "ada", "time", "ada", "Trace")
   CALL add_2p_conv_trace(number_op, number_op, test_min_sep, test_max_sep, "ninj", "time", "ninj", "Trace")
   CALL add_2p_conv_trace(eye,       eye,       test_min_sep, test_max_sep, "IdId", "time", "IdId", "Trace")

   !get dump index -- DUMP_NUM = first file created after the convergence files
   CALL prepare_datafile(DUMP_NUM, "dump_ind", "Index of states for reloading")

   !Additionally setup summary files to record SS observables
   ALLOCATE(SUMM_OUT_1P(3), SUMM_OP_1P(3, local_dim, local_dim))

   CALL add_1p_summary(dn_op,     "a",  xname, "a",  "Trace")
   CALL add_1p_summary(number_op, "n",  xname, "n",  "Trace")
   CALL add_1p_summary(g2bare_op, "G2", xname, "G2", "Trace")

   ALLOCATE(SUMM_OUT_2P(2), SUMM_OP_2P(2, 2, local_dim, local_dim))

   CALL add_2p_summary(up_op,     dn_op,     "ada", xname, "ada", "Trace")
   CALL add_2p_summary(number_op, number_op, "nn",  xname, "nn",  "Trace")

   !Extra file for post-processed g2
   CALL setup_file_obs_2P(datafile=G2_SUMMARY, descriptor="g2_post_processed", param_str=xname, sepX_str="|i-j|", obs_str="g2", extra_str_1="nval")

 END SUBROUTINE setup_observables

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM run_kbhm_imps
