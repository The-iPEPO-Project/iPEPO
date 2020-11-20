PROGRAM run_kbhm_tramte_evecs

 USE utility
 USE datafile_utility
 USE definitions_mps_mpo
 USE mps_mpo_utility
 USE tramte_network
 USE power_method

 USE time_evolution_imps, ONLY: initialize_rho_imps
 USE mps_peps_INOUT_files, ONLY: reload_mps
 USE simulation_parameters
 USE kbhm_propagator

 IMPLICIT NONE

 !! MPS RHO NESS
 TYPE(block_mps),     ALLOCATABLE :: mps(:)
 TYPE(block_mps),     ALLOCATABLE :: iG(:)
 TYPE(block_lambda),  ALLOCATABLE :: iL(:)

 !! Propagator && 2D Trotter network
 TYPE(block_mpo),     ALLOCATABLE :: prop_mpo(:)
 TYPE(block_mpo),     ALLOCATABLE :: TN2D(:)

 !! Tmat && evecs
 TYPE(block_mpo),     ALLOCATABLE :: Tmat(:,:)
 TYPE(block_mps),     ALLOCATABLE :: Rvec(:), Lvec(:)

 !! Physical params, timestep
 REAL(KIND=DP)    :: delta, U, f, kappa, modJ, argJ
 COMPLEX(KIND=DP) :: J
 COMPLEX(KIND=DP) :: nval, g2val
 COMPLEX(KIND=DP) :: dt
 REAL(KIND=DP)    :: real_dt

 !! Simulation params
 TYPE(power_params) :: POWER
 REAL(KIND=DP)      :: eps(2)
 INTEGER            :: chi(2)
 REAL(KIND=DP)      :: eta
 REAL(KIND=DP)      :: xi

 !! Start & end points of time evolution, and time interval between evec calculations
 !! (NB calc_interval must be an even number to ensure convergence!)
 INTEGER :: N_start, N_end, deltaN
 INTEGER :: N_steps

 !! Reloading evecs
 CHARACTER :: read_evecs_str*1
 LOGICAL   :: read_evecs
 CHARACTER :: filename_Rvec*256, filename_Lvec*256, evecFlag*16
 LOGICAL   :: use_last_evec

 !! Reloading MPS NESS
 CHARACTER(LEN=256) :: filename_RHO

 !! Cpu times
 REAL(KIND=DP) :: cputimes(2)

 !! Starting time
 CALL cpu_time(cputimes(1))

 !! Physical params
 WRITE(*,*) "delta, U, f, kappa, |J|, arg(J) = k "
 READ(*,*) delta, U, f, kappa, modJ, argJ

 !! Combine modulus and phase to get j
 J=modJ*exp(ii*argJ)

 !! Hilbert space
 WRITE(*,*) "Hilbert space size"
 READ(*,*) basis_size

 !! Timestep
 WRITE(*,*) "Timestep"
 READ(*,*) real_dt
 dt = CMPLX(real_dt, 0.0D0)
 WRITE(*,*) "Complex dt = ", dt

 !! Start, end, interval of time evolution
 WRITE(*,*) "N_start, N_end, deltaN"
 READ(*,*) N_start, N_end, deltaN

 !! TraMTE precision params
 WRITE(*,*) "TraMTE precision -- eps"
 READ(*,*) eps
 WRITE(*,*) "TraMTE precision -- chi"
 READ(*,*) chi
 WRITE(*,*) "TraMTE precision -- eta"
 READ(*,*) eta

 !! Propagator precision
 WRITE(*,*) "Propagator precision -- xi"
 READ(*,*) xi

 !! Reading NESS RHO from file
 WRITE(*,*) "Reading state from file. Filename?"
 READ(*,*) filename_RHO

 !! Evec we're calculating
 WRITE(*,*) "Evecs to calculate?"
 READ(*,*) evecFlag

 !! Reload evecs from file?
 WRITE(*,*) "Reload existing state?"
 READ(*,*) read_evecs_str

 IF(read_evecs_str .EQ. "Y") THEN

     !! Reload evecs from file
     read_evecs=.TRUE.
     WRITE(*,*) "Reload Rvec & Lvec from file"
     READ(*,*) filename_Rvec, filename_Lvec 

     !! If reloading evecs:
     !! -- reload evecs of size = N_start
     !! -- increment them to  N_start --> N_start + deltaN
     N_start = N_start + deltaN
 ELSE
     read_evecs=.FALSE.
 end IF

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !! Set power method params
 CALL create_power_params(POWER, eps, chi, eta)

 !! Set up operators -- MUST be called before constructing propagator
 CALL setup_basic_operators()

 !! Reload RHO state (will become N-bound of CTM)
 CALL initialize_rho_imps(iG,  iL, filename=filename_RHO)
 CALL absorb_imps_lambdas(mps, iG, iL)

 !! Setup propagator and the corresponding 2D Trotter network
 CALL setup_propagator(prop_mpo, delta, f, U, J, kappa, 0.5D0*dt, xi)
 CALL create_reduced_TN_from_trotter(TN2D, prop_mpo)


 !!!!!!!!!!!!!!!!! Evolve in time computing boundary evecs every [deltaN] steps !!!!!!!!!!!!!!!!!!!!!!!
 DO N_steps = N_start, N_end, deltaN 

    !! Construct transfer matrix MPO
    CALL construct_transfer_matrix(Tmat, MPS, TN2D, N_sites=N_steps)

    !! Decide whether to initialize using the last converged evec
    IF(read_evecs .AND. (N_steps .EQ. N_start)) THEN 

        !! Reload boundary evecs
        SELECT CASE(evecFlag)
        CASE('Rvec') 
             CALL reload_mps(filename_Rvec, Tmat(2,:), Rvec)
        CASE('Lvec')
             CALL reload_mps(filename_Lvec, Tmat(2,:), Lvec)
        CASE DEFAULT
             CALL reload_mps(filename_Rvec, Tmat(2,:), Rvec)
             CALL reload_mps(filename_Lvec, Tmat(2,:), Lvec)
        end SELECT

        !! Since we're reloading last evec --> set use_last_evec = .TRUE.
        use_last_evec = .TRUE.
    ELSE
        !! Increment current evec if the below is true
        use_last_evec = (N_steps .GT. N_start) .AND. ((N_steps - deltaN) .GT. 2)
    end IF

    !! Compute evecs of Tmat using MPS power method
    CALL calc_mpsevecs_power_meth(Rvec, Lvec, Tmat, POWER, use_last_evec, N_sites=N_steps, dN=deltaN, evecFlag=evecFlag) 
 end DO

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !! Finishing time
 CALL cpu_time(cputimes(2))

 WRITE(*,*) "total time = ", cputimes(2) - cputimes(1)

CONTAINS

END PROGRAM run_kbhm_tramte_evecs
