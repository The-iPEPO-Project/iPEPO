PROGRAM run_ising_tramte_obs

 USE utility
 USE datafile_utility
 USE definitions_mps_mpo
 USE mps_mpo_utility
 USE tramte_network

 USE time_evolution_imps,  ONLY: initialize_rho_imps
 USE mps_peps_INOUT_files, ONLY: reload_mps
 USE simulation_parameters
 USE ising_propagator_pauli_basis

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
 REAL(KIND=DP)    :: g, delta, kappa
 COMPLEX(KIND=DP) :: dt
 REAL(KIND=DP)    :: real_dt

 !! Precision params
 REAL(KIND=DP)  :: eps(2)
 INTEGER        :: chi(2)
 REAL(KIND=DP)  :: xi

 !! Current point of time evolution
 INTEGER :: N_steps


 !! Reloading evecs
 CHARACTER :: read_evecs_str*1
 LOGICAL   :: read_evecs
 CHARACTER :: filename_Rvec*256, filename_Lvec*256, evecFlag*16

 !! Reloading MPS NESS
 CHARACTER(LEN=256) :: filename_RHO

 !! Observables
 COMPLEX(KIND=DP), ALLOCATABLE :: expvals(:)
 COMPLEX(KIND=DP)              :: expSS
 TYPE(op_sites_type)           :: TN_OpSites
 CHARACTER(LEN=16)             :: corr_str, N_str, descr
 INTEGER                       :: min_sep, max_sep
 INTEGER                       :: i

 !! Cpu times
 REAL(KIND=DP) :: cputimes(2)

 !! Starting time
 CALL cpu_time(cputimes(1))

 !! Physical params
 WRITE(*,*) "kappa, g, delta"
 READ(*,*) kappa, g, delta

 !! Timestep
 WRITE(*,*) "Timestep"
 READ(*,*) real_dt
 dt = CMPLX(real_dt, 0.0D0)
 WRITE(*,*) "Complex dt = ", dt

 !! Calculate observables after N_steps
 WRITE(*,*) "Number of timesteps"
 READ(*,*) N_steps

 !! Sep for correlators
 WRITE(*,*) "min_separation, max_separation?"
 READ(*,*) min_sep, max_sep

 !! TraMTE precision params
 WRITE(*,*) "TraMTE precision -- eps"
 READ(*,*) eps
 WRITE(*,*) "TraMTE precision -- chi"
 READ(*,*) chi

 !! Propagator precision
 WRITE(*,*) "Propagator precision -- xi"
 READ(*,*) xi

 !! Reading NESS && conj NESS from file
 WRITE(*,*) "Reading state from file. Filename?"
 READ(*,*) filename_RHO

 !! Reload evecs from file
 WRITE(*,*) "Reload Rvec & Lvec from file"
 READ(*,*) filename_Rvec, filename_Lvec 

 !! Correlator to calculate
 WRITE(*,*) "Correlator to calculate"
 READ(*,*) corr_str

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !! Set up operators -- MUST be called before constructing propagator
 CALL setup_basic_operators()

 !! Reload RHO state (will become N-bound of CTM)
 CALL initialize_rho_imps(iG,  iL, filename=filename_RHO)
 CALL absorb_imps_lambdas(mps, iG, iL)

 !! Setup propagator and the corresponding 2D Trotter network
 CALL setup_propagator(prop_mpo, g, delta, kappa, 0.0D0, dt, xi)
 CALL create_reduced_TN_from_trotter(TN2D, prop_mpo)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! setup datafiles and observables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ALLOCATE(SUMM_OUT_1P(1), SUMM_OUT_2P(1), SUMM_OP_2P(1, 2, local_dim, local_dim),  expvals(max_sep+1-min_sep))

  WRITE(N_str,'(I6)') N_steps
  descr = TRIM(ADJUSTL(corr_str))//"_N="//TRIM(ADJUSTL(N_str))

  SELECT CASE(corr_str)
        CASE('Sxt_Sx0')
            CALL add_2p_summary(sigma_x,  sigma_x,             descr,  'timestep', 'sep', 'Sx(t)Sx(0)', 'Sx(0)')
        CASE('Sx0_Sxt')
            CALL add_2p_summary(sigma_x,  TRANSPOSE(sigma_x),  descr,  'timestep', 'sep', 'Sx(0)Sx(t)', 'Sx(0)')
        CASE('Syt_Sy0')
            CALL add_2p_summary(sigma_y,  sigma_y,             descr,  'timestep', 'sep', 'Sy(t)Sy(0)', 'Sy(0)')
        CASE('Sy0_Syt')
            CALL add_2p_summary(sigma_y,  TRANSPOSE(sigma_y),  descr,  'timestep', 'sep', 'Sy(0)Sy(t)', 'Sy(0)')
        CASE('Szt_Sz0')
            CALL add_2p_summary(sigma_z,  sigma_z,             descr,  'timestep', 'sep', 'Sz(t)Sz(0)', 'Sz(0)')
        CASE('Sz0_Szt')
            CALL add_2p_summary(sigma_z,  TRANSPOSE(sigma_z),  descr,  'timestep', 'sep', 'Sz(0)Sz(t)', 'Sz(0)')
        CASE('Idt_Id0')
            CALL add_2p_summary(eye,      eye,                 descr,  'timestep', 'sep', 'Id(t)Id(0)', 'Id(0)')
        CASE('Id0_Idt')
            CALL add_2p_summary(eye,      TRANSPOSE(eye),      descr,  'timestep', 'sep', 'Id(0)Id(t)', 'Id(0)')
        CASE('upt_dn0')
            CALL add_2p_summary(up,       down,                descr,  'timestep', 'sep', 'up(t)dn(0)', 'dn(0)')
        CASE('dnt_up0')
            CALL add_2p_summary(down,     up,                  descr,  'timestep', 'sep', 'dn(t)up(0)', 'dn(0)')
        CASE DEFAULT        
            CALL invalid_flag("run_ising_tramte_obs -- invalid corr_str ", corr_str)
        end SELECT

 !! Create operator sites to be inserted into 2D Trotter network
 !! (NB. we apply operators in order SUMM(1,1,:,:) * SUMM(1,2,:,:) * rho --> i.e. set op1 = SUMM(1,2,:,:), op2 = SUMM(1,1,:,:))
 CALL create_reduced_op_sites_trotter(TN_OpSites, TN2D, SUMM_OP_2P(1,2,:,:), SUMM_OP_2P(1,1,:,:), ROT='PI/2-ROTATED')

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !! Construct transfer matrix MPO && TN operator sites
 CALL construct_transfer_matrix(Tmat, MPS, TN2D, N_steps)

 !! Reload boundary evecs
 CALL reload_mps(filename_Rvec, Tmat(2,:), Rvec)
 CALL reload_mps(filename_Lvec, Tmat(2,:), Lvec)

 !! Contract TraMTE network --> get expvals at different separations
 CALL contract_tramte_network(expvals, Rvec, Tmat, Lvec, chi, eps, TN_OpSites, min_sep, max_sep)

 !! Calc factorized correlator
 CALL calc_factorized_corr(expSS, iG, iL, SUMM_OP_2P(1,2,:,:))

 !! Write expvals to file
 DO i=1,SIZE(expvals)
    CALL write_summary_data_2P(SUMM_OUT_2P(1), param=ABS(dt)*(N_steps/2), sep=i-1+min_sep, expval=expvals(i) - expSS, xval1=expSS)
 end DO

 !! Close output files
 CALL close_datafiles()

CONTAINS

END PROGRAM run_ising_tramte_obs
