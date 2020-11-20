PROGRAM run_xyz_ipeps

  USE utility
  USE datafile_utility
  USE definitions_mps_mpo
  USE definitions_peps_pepo

  USE simulation_parameters
  USE xyz_propagator

  USE observables_ipeps
  USE time_evolution_ipeps

  IMPLICIT NONE

  !! Rho state
  TYPE(block_peps),   ALLOCATABLE :: ipeps(:)
  TYPE(block_peps),   ALLOCATABLE :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE :: iLambda(:)

  !! CTM environment
  TYPE(ctm_corner_type),   ALLOCATABLE :: Cmat(:), CC(:)
  TYPE(ctm_transfer_type), ALLOCATABLE :: Tmat(:), TT(:)
  
  !! Propagator
  !TYPE(block_mpo),  ALLOCATABLE :: prop(:)
  COMPLEX(KIND=DP),  ALLOCATABLE :: prop(:,:)

  !! Simulation params
  TYPE(tebd_params) :: TEBD
  TYPE(ctm_params)  :: CTMRG
  REAL(KIND=DP)     :: eps, epsENV
  REAL(KIND=DP)     :: eta, etaENV
  INTEGER           :: chi, chiENV
  
  !! Physical params
  REAL(KIND=DP) :: kappa, Jx, Jy, Jz

  !! Initialization, Reloading && Output Printing vars
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_vec(:)
  REAL(KIND=DP) :: s0(3)
  CHARACTER     :: read_state_filename*256, read_state_str*1
  CHARACTER     :: read_ctm_filename*256  
  LOGICAL       :: read_state

  !! Range of correlations (summary files)
  INTEGER :: min_sep, max_sep  

  !! TEBD params - Main iteration
  INTEGER, PARAMETER :: N_tebd = 5
  INTEGER            :: j, j0, HANDLE 

  COMPLEX(KIND=DP)   :: dt(N_tebd)
  REAL(KIND=DP)      :: real_dt(N_tebd)
  INTEGER            :: N_steps(N_tebd), test_interval(N_tebd)

  INTEGER, PARAMETER :: test_min_sep = 1, test_max_sep = 2
  INTEGER, PARAMETER :: exit_level(N_tebd) = (/1, 2, 2, 2, 2/) !(/0, 0, 0, 0, 0/)

  !! Parameter loop vars
  REAL(KIND=DP) :: x, x0, dx
  INTEGER       :: nx, ix
  CHARACTER     :: xname*16

  !! Cpu times
  REAL(KIND=DP) :: cpu_times(4)

  !! String vars for SS output filename
  CHARACTER :: SS_state_filename*256, char_chi*16, char_x*16

  WRITE(*,*) "kappa, Jx, Jy, Jz"
  READ(*,*) kappa, Jx, Jy, Jz

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

  WRITE(*,*) "TEBD precision settings"
  READ(*,*) eta, chi, eps

  WRITE(*,*) "Environment precision settings"
  READ(*,*) etaENV, chiENV, epsENV

  WRITE(*,*) "Physical param to vary, Initial val, Increment, Number of vals"
  READ(*,*) xname, x0, dx, nx

  WRITE(*,*) "Reload existing state?"
  READ(*,*) read_state_str

  IF(SCAN(read_state_str,"yYtT") .GT. 0) THEN

     read_state=.TRUE.
     WRITE(*,*) "Reading state from file. Filename?"
     READ(*,*) read_state_filename
     READ(*,*) read_ctm_filename
  ELSE
     read_state=.FALSE.
     WRITE(*,*) "Initial (product) state, sx,sy,sz"
     READ(*,*) s0
  end IF

  CALL cpu_time(cpu_times(1))

  !! Create a new instance of CTM
  CALL create_CTM_params(CTMRG, epsENV, chiENV, etaENV, ctm_method='3')
  
  !! Set up operators -- MUST be called before initializing rho state
  CALL setup_basic_operators('RHO BASIS')

  !! RHO state: either reload, or initialize to a product state
  IF(read_state) THEN
     CALL initialize_rho_iPEPS(iGamma, iLambda, local_dim, filename=read_state_filename)
     CALL reload_rho_ipeps_CTM(Cmat, Tmat, iGamma, iLambda, read_ctm_filename) 
  ELSE
     CALL initialize_rho_iPEPS(iGamma, iLambda, local_dim, chi=chi)
     !! CALL create_initial_state_vec(rho_vec, s0)
     !! CALL initialize_rho_iPEPS(iGamma, iLambda, local_dim, rho_vec=rho_vec)
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
     CALL read_xname(x, xname, Jy, "Jy", kappa, "kappa")

     !! REAL-TIME TEBD propagation using decreasing values of timestep dt 
     DO j=1,3 !N_tebd

        !! Create a new instance of TEBD
        CALL create_TEBD_params(TEBD, & 
                                      & eps, chi, eta, & 
                                      & dt(j), N_steps(j), test_interval(j), test_min_sep, test_max_sep, exit_level(j))
        !! Create propagator
        CALL setup_propagator(prop, Jx, Jy, Jz, kappa, TEBD % dt, 'PEPS')

        !! Trotter time evolution to obtain NESS
        CALL time_evolve(iGamma, iLambda, prop, TEBD)

        !! Calculate final environment
        CALL calc_CTM_environment_ipeps(Cmat, Tmat, iGamma, iLambda, CTMRG, use_old_ctm = read_state .OR. (ix .GT. 1) .OR. (j .GT. 1)) 

        !! Calculate summary observables and output SS (steady state) 
        CALL calc_TEBD_summary(iGamma, iLambda, Cmat, Tmat, CTMRG, x, min_sep, max_sep, dump_parameters=dump_parameters)

        !! Add linebreak to datafiles 
        DO HANDLE=FIRST_HANDLE, DUMP_NUM
           WRITE(HANDLE, *)
           WRITE(HANDLE, *)
        end DO
     END DO
  END DO

  CALL cpu_time(cpu_times(2))

  !! Simulation time
  WRITE(*,*) "TEBD time = ", cpu_times(2) - cpu_times(1)

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

    USE simulation_parameters

    INTEGER, INTENT(IN) :: HANDLE
    
    WRITE(HANDLE,*)
    WRITE(HANDLE,*)
    WRITE(HANDLE, '("# local_dim: ", I5)') local_dim
    WRITE(HANDLE, '("# timestep : ", 6(D12.5,X))'),  dt
    WRITE(HANDLE, '("# kappa= ", D12.5, ", gx= ", D12.5, " J=", 3(D12.5,X))'), kappa, Jx, Jy, Jz
    WRITE(HANDLE, '("# eps: ", D12.5, "  chi : ", I5, "  eta : ", D12.5)')          TEBD % eps,   TEBD % chi,   TEBD % eta
    WRITE(HANDLE, '("# epsENV: ", D12.5, "  chiENV : ", I5, "  etaENV : ", D12.5)') CTMRG % eps,  CTMRG % chi,  CTMRG % eta
    WRITE(HANDLE, '("# Lambda MIN/MAX: ", D12.5)') lambda_min_max(iLambda)
    WRITE(HANDLE,*)
    WRITE(HANDLE,*)

 END SUBROUTINE dump_parameters

 
 !!!!!!!!!!!!!!!!!!! Setup observables (operator lists and files) !!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE setup_observables()

  !! Setup convergence files to record observables during time evolution
  ALLOCATE(TRACE_OUT_1P(4), TRACE_OP_1P(4, hs_dim, hs_dim))

  CALL add_1p_conv_trace(sigma_x, "Sx", "time", "Sx", "Trace")
  CALL add_1p_conv_trace(sigma_y, "Sy", "time", "Sy", "Trace")
  CALL add_1p_conv_trace(sigma_z, "Sz", "time", "Sz", "Trace")
  CALL add_1p_conv_trace(eye,     "Id", "time", "Id", "Trace")

  ALLOCATE(TRACE_OUT_2P(5, test_min_sep : test_max_sep), TRACE_OP_2P(5, 2, hs_dim, hs_dim))

  CALL add_2p_conv_trace(sigma_x, sigma_x, test_min_sep, test_max_sep, "SxSx", "time", "SxSx", "Trace")
  CALL add_2p_conv_trace(sigma_y, sigma_y, test_min_sep, test_max_sep, "SySy", "time", "SySy", "Trace")
  CALL add_2p_conv_trace(sigma_x, sigma_y, test_min_sep, test_max_sep, "SxSy", "time", "SxSy", "Trace")
  CALL add_2p_conv_trace(sigma_y, sigma_x, test_min_sep, test_max_sep, "SySx", "time", "SySx", "Trace")
  CALL add_2p_conv_trace(eye,     eye,     test_min_sep, test_max_sep, "IdId", "time", "IdId", "Trace")

  !! Get dump index -- DUMP_NUM = first file created after the convergence files
  CALL prepare_datafile(DUMP_NUM, "dump_ind", "Index of states for reloading")

  !! Setup summary files to record SS observables
  ALLOCATE(SUMM_OUT_1P(4), SUMM_OP_1P(4, hs_dim, hs_dim))

  CALL add_1p_summary(sigma_x, "Sx", xname, "Sx", "Trace")
  CALL add_1p_summary(sigma_y, "Sy", xname, "Sy", "Trace")
  CALL add_1p_summary(sigma_z, "Sz", xname, "Sz", "Trace")
  CALL add_1p_summary(eye,     "Id", xname, "Id", "Trace")

  ALLOCATE(SUMM_OUT_2P(4), SUMM_OP_2P(4, 2, hs_dim, hs_dim))

  CALL add_2p_summary(sigma_x, sigma_x, "SxSx", xname, "SxSx", "Trace")
  CALL add_2p_summary(sigma_y, sigma_y, "SySy", xname, "SySy", "Trace")
  CALL add_2p_summary(sigma_x, sigma_y, "SxSy", xname, "SxSy", "Trace")
  CALL add_2p_summary(sigma_y, sigma_x, "SySx", xname, "SySx", "Trace")

 END SUBROUTINE setup_observables

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM run_xyz_ipeps
