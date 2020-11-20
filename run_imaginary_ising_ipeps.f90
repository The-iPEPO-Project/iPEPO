PROGRAM run_imaginary_ising_ipeps

  USE utility
  USE datafile_utility
  USE definitions_mps_mpo
  USE definitions_peps_pepo

  USE simulation_parameters
  USE ising_imaginary_propagator

  USE imaginary_observables_ipeps
  USE imaginary_time_evolution_ipeps

  USE imaginary_full_update_ipeps

  IMPLICIT NONE

  !! Psi state
  TYPE(block_peps),   ALLOCATABLE :: ipeps(:)
  TYPE(block_peps),   ALLOCATABLE :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE :: iLambda(:)

  !! CTM environment
  TYPE(ctm_corner_type),   ALLOCATABLE :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE :: Tmat(:)
  
  !! Mpo propagator
  TYPE(block_mpo), ALLOCATABLE :: mpo_prop(:)

  !! Simulation params
  TYPE(tebd_params) :: tebd
  TYPE(ctm_params)  :: CTMRG
  REAL(KIND=DP)     :: eps, epsENV
  REAL(KIND=DP)     :: eta, etaENV
  INTEGER           :: chi, chiENV
  
  !! Physical params
  REAL(KIND=DP) :: g, delta

  !! Initialization, Reloading && Output Printing vars
  COMPLEX(KIND=DP), ALLOCATABLE :: psi_vec(:)
  REAL(KIND=DP) :: s0(2)
  CHARACTER     :: read_state_filename*256, read_state_str*1
  CHARACTER     :: read_ctm_filename*256  
  LOGICAL       :: read_state

  !! Range of correlations (summary files)
  INTEGER :: min_sep, max_sep  

  !! TEBD params - Main iteration
  INTEGER, PARAMETER :: N_tebd = 5
  INTEGER            :: j, HANDLE 

  COMPLEX(KIND=DP)   :: dt(N_tebd)
  REAL(KIND=DP)      :: real_dt(N_tebd)
  INTEGER            :: N_steps(N_tebd), test_interval(N_tebd)

  INTEGER, PARAMETER :: test_min_sep = 1, test_max_sep = 2
  INTEGER, PARAMETER :: exit_level(N_tebd) = (/2, 2, 2, 2, 2/)

  !! Parameter loop vars
  REAL(KIND=DP) :: x, x0, dx
  INTEGER       :: nx, ix
  CHARACTER     :: xname*16

  !! Cpu times
  REAL(KIND=DP) :: cpu_times(4)

  !! String vars for SS output filename
  CHARACTER :: SS_state_filename*256, char_chi*16, char_x*16

  WRITE(*,*) "g, delta"
  READ(*,*) g, delta

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

  !! Set CTM params
  CALL create_CTM_params(CTMRG, epsENV, chiENV, etaENV, ctm_method='3')

  !! Set up operators -- MUST be called before initializing Psi state
  CALL setup_basic_operators()

  !! PSI state: either reload, or initialize to a product state
  IF(read_state) THEN
     CALL initialize_psi_iPEPS(iGamma, iLambda, local_dim, filename=read_state_filename)
     CALL reload_psi_ipeps_CTM(Cmat, Tmat, iGamma, iLambda, read_ctm_filename)
     CALL resize_ipeps(iGamma, iLambda, Cmat, Tmat, chi)
  ELSE
     ALLOCATE(psi_vec(2)) 
     CALL initialize_rand_tens(psi_vec)
     CALL initialize_psi_iPEPS(iGamma, iLambda, local_dim, psi_vec, chi=chi)
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
     CALL read_xname(x, xname, g, "g", delta, "delta")

     WRITE(*,*)
     WRITE(*,*) "Starting x = ", x

     !! (1) Simple Update TEBD propagation using decreasing values of timestep dt 
     DO j=1,N_tebd 

        !! Create a new instance of TEBD
        CALL create_TEBD_params(TEBD, & 
                                & eps, chi, eta, & 
                                & dt(j), N_steps(j), test_interval(j), test_min_sep, test_max_sep, exit_level(j))

        !! Create propagator
        CALL setup_propagator(mpo_prop, g, delta, TEBD % dt, 'PEPS')

        !! Trotter time evolution to obtain GS
        CALL time_evolve(iGamma, iLambda, mpo_prop, TEBD)

        !! Add linebreak to datafiles 
        DO HANDLE=FIRST_HANDLE, DUMP_NUM
           WRITE(HANDLE, *)
           WRITE(HANDLE, *)
        end DO
     END DO

     !! (2) Calculate environment of iPEPS obtained from SU
     CALL calc_CTM_environment_ipeps(Cmat, Tmat, iGamma, iLambda, CTMRG, use_old_ctm = read_state .OR. (ix .GT. 1)) 

     !! (3) Calculate summary observables and output SS (steady state) 
     CALL calc_TEBD_summary(iGamma, iLambda, Cmat, Tmat, CTMRG, x, min_sep, max_sep, dump_parameters=dump_parameters)

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
    WRITE(HANDLE, '("# local_dim: ", I5, "bond_dim : ", I5)') local_dim, chi
    WRITE(HANDLE, '("# timestep : ", 2(D12.5,X))'),  dt
    WRITE(HANDLE, '("# g= ", D12.5, " delta=", D12.5)') g, delta
    WRITE(HANDLE,*)
    WRITE(HANDLE,*)

 END SUBROUTINE dump_parameters



 !!!!!!!!!!!!!!!!!!! Setup observables (operator lists and files) !!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE setup_observables()

  !! Setup convergence files to record observables during time evolution 
  ALLOCATE(TRACE_OUT_1P(4), TRACE_OP_1P(4, local_dim, local_dim))

  CALL add_1p_conv_trace(sigma_x, "Sx", "time", "Sx", "Norm")
  CALL add_1p_conv_trace(sigma_y, "Sy", "time", "Sy", "Norm")
  CALL add_1p_conv_trace(sigma_z, "Sz", "time", "Sz", "Norm")
  CALL add_1p_conv_trace(eye,     "Id", "time", "Id", "Norm")

  ALLOCATE(TRACE_OUT_2P(4, test_min_sep : test_max_sep), TRACE_OP_2P(4, 2, local_dim, local_dim))

  CALL add_2p_conv_trace(sigma_z, sigma_z, test_min_sep, test_max_sep, "SzSz", "time", "SzSz", "Norm")
  CALL add_2p_conv_trace(sigma_x, sigma_x, test_min_sep, test_max_sep, "SxSx", "time", "SxSx", "Norm")
  CALL add_2p_conv_trace(sigma_y, sigma_y, test_min_sep, test_max_sep, "SySy", "time", "SySy", "Norm")
  CALL add_2p_conv_trace(eye,     eye,     test_min_sep, test_max_sep, "IdId", "time", "IdId", "Norm")

  !! Get dump index -- DUMP_NUM = first file created after the convergence files
  CALL prepare_datafile(DUMP_NUM, "dump_ind", "Index of states for reloading")

  !! Setup summary files to record SS observables
  ALLOCATE(SUMM_OUT_1P(4), SUMM_OP_1P(4, local_dim, local_dim))

  CALL add_1p_summary(sigma_x, "Sx", xname, "Sx", "Norm")
  CALL add_1p_summary(sigma_y, "Sy", xname, "Sy", "Norm")
  CALL add_1p_summary(sigma_z, "Sz", xname, "Sz", "Norm")
  CALL add_1p_summary(eye,     "Id", xname, "Id", "Norm")

  ALLOCATE(SUMM_OUT_2P(3), SUMM_OP_2P(3, 2, local_dim, local_dim))

  CALL add_2p_summary(sigma_z, sigma_z, "SzSz", xname, "SzSz", "Norm")
  CALL add_2p_summary(sigma_x, sigma_x, "SxSx", xname, "SxSx", "Norm")
  CALL add_2p_summary(sigma_y, sigma_y, "SySy", xname, "SySy", "Norm")

 END SUBROUTINE setup_observables

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM run_imaginary_ising_ipeps
