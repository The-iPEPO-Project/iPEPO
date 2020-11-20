MODULE time_evolution_ipeps

 USE utility
 USE definitions_mps_mpo
 USE definitions_peps_pepo
 USE mps_peps_INOUT_files
 USE mps_mpo_algebra_inf

 USE peps_pepo_algebra
 USE project_ipeps

 USE simulation_parameters
 USE observables_ipeps
 USE tebd_callback_routines_ipeps

 USE TN_contractions

 IMPLICIT NONE

 interface time_evolve
    module procedure time_evolve_simple
    module procedure time_evolve_simple_MPO
    module procedure time_evolve_simple_with_obs
 end interface time_evolve

 interface calc_tebd_summary
    module procedure calc_tebd_summary_PLAIN
    module procedure calc_tebd_summary_GamLam
 end interface calc_tebd_summary

 interface output_converged_ipeps
    module procedure output_converged_ipeps_PLAIN
    module procedure output_converged_ipeps_GamLam
 end interface output_converged_ipeps

 private time_evolve_simple,  time_evolve_simple_MPO,  time_evolve_simple_with_obs

 private calc_tebd_summary_PLAIN,       calc_tebd_summary_GamLam
 private output_converged_ipeps_PLAIN,  output_converged_ipeps_GamLam

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEBD algorithm for evolving iPEPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! TEBD simple update -- lambda convergence testing, no environment && observables !!!
 SUBROUTINE time_evolve_simple(iGamma, iLambda, prop_dense, TEBD) 

  USE datafile_utility, ONLY: setupConvDatafile, write_conv_data

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGamma(:)            !iGamma of iPEPS
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)           !iLambda of iPEPS
  COMPLEX(KIND=DP),                INTENT(IN)    :: prop_dense(:,:)      !Propagator
  TYPE(tebd_params),               INTENT(IN)    :: TEBD                 !TEBD params

  !! TEBD iteration vars
  INTEGER          :: step, i
  COMPLEX(KIND=DP) :: time, time_offset
  INTEGER          :: datafile, datafile_Sop
  LOGICAL          :: is_converged

  !! Old iLambda for lambda convergence testing
  TYPE(block_peps),   ALLOCATABLE :: iGammaOLD(:)
  TYPE(block_lambda), ALLOCATABLE :: iLambdaOLD(:)
  REAL(KIND=DP)                   :: lamErr(3)
  COMPLEX(KIND=DP)                :: Sop

  !! Norm const
  COMPLEX(KIND=DP) :: C_norm

  !! Initialize time && convergence tester
  time_offset = 0.0D0
  time = 0.0D0
  is_converged = .FALSE.

  !! Setup Lambda convergence file, initialize old Lambda
  CALL copy_lambda_block(iLambdaOLD, iLambda)
  CALL setupConvDatafile(datafile=datafile, descr="iPEPS_TEBD_lambda", errstr="Err_lambda", valstr="Max_lambda", xstr="dt", xlog1="use_old_ipeps")

  !! Setup file for operator entanglement entropy
  CALL setupConvDatafile(datafile=datafile_Sop, descr="Sop_entanglement", valstr="Sop")

  WRITE(*,*) " Starting iPEPS TEBD "

  !!! TIME EVOLUTION -- THE MAIN LOOP (using 1st order Trotter Decomposition) !!!
  mainloop:DO step = 1, (TEBD % N_steps)

    !! WRITE(*,*)
    !! WRITE(*,*) "TIMESTEP = ", step
    !! WRITE(*,*)

    !! (1) Evolve iPEPS bonds
    CALL mult_ipeps_dense_prop(iGamma, iLambda, prop_dense, TEBD % chi, TEBD % eps,  'S')
    CALL mult_ipeps_dense_prop(iGamma, iLambda, prop_dense, TEBD % chi, TEBD % eps,  'W')
    CALL mult_ipeps_dense_prop(iGamma, iLambda, prop_dense, TEBD % chi, TEBD % eps,  'N')
    CALL mult_ipeps_dense_prop(iGamma, iLambda, prop_dense, TEBD % chi, TEBD % eps,  'E')

    !! (2) Rescale iPEPS sites
    iGamma(1) % m = iGamma(1) % m / MAXVAL(ABS(iGamma(1) % m))
    iGamma(2) % m = iGamma(2) % m / MAXVAL(ABS(iGamma(2) % m))

    !! (3) Increment time by a single step after applying propagator to all bonds
    time = step * (TEBD % dt) + time_offset

    !! (4) Test convergence of TEBD iLambda spectrum
    IF((TEBD % test_interval .EQ. 1) .OR. (MOD(step, TEBD % test_interval) .EQ. 0)) THEN 
        CALL test_lambda_convergence(is_converged, iLambda, iLambdaOLD, TEBD, datafile, step, ABS(time))
    end IF

    !! (5) If converged --> exit TEBD loop
    IF(is_converged .AND. (TEBD % exit_level .GE. 1)) THEN
       EXIT mainloop
    ELSEIF((step .EQ. TEBD % N_steps) .AND. (TEBD % exit_level .GE. 2)) THEN
       WRITE(*,*) "iPEPS TEBD HAS FAILED TO CONVERGE "
       STOP
    end IF

  end DO mainloop

  CLOSE(datafile_Sop)
  WRITE(datafile, *)
  WRITE(datafile, *)
  CLOSE(datafile)
  WRITE(*,*)
  WRITE(*,*) "FINISHED TEBD"
  WRITE(*,*)

 END SUBROUTINE time_evolve_simple







 !!! TEBD simple update -- lambda convergence testing, no environment && observables !!!
 SUBROUTINE time_evolve_simple_MPO(iGamma, iLambda, prop, TEBD) 

  USE datafile_utility, ONLY: setupConvDatafile, write_conv_data

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGamma(:)   !iGamma of iPEPS
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)  !iLambda of iPEPS
  TYPE(block_mpo),                 INTENT(IN)    :: prop(:)     !Propagator
  TYPE(tebd_params),               INTENT(IN)    :: TEBD        !TEBD params

  !! TEBD iteration vars
  INTEGER          :: step, i
  COMPLEX(KIND=DP) :: time, time_offset
  INTEGER          :: datafile, datafile_Sop
  LOGICAL          :: is_converged

  !! Old iLambda for lambda convergence testing
  TYPE(block_lambda), ALLOCATABLE :: iLambdaOLD(:)
  COMPLEX(KIND=DP)                :: Sop

  !! Initialize time && convergence tester
  time_offset = 0.0D0
  time = 0.0D0
  is_converged = .FALSE.

  !! Setup Lambda convergence file, initialize old Lambda
  CALL copy_lambda_block(iLambdaOLD, iLambda)
  CALL setupConvDatafile(datafile=datafile, descr="iPEPS_TEBD_lambda", errstr="Err_lambda", valstr="Max_lambda", xstr="dt", xlog1="use_old_ipeps")

  !! Setup file for operator entanglement entropy
  CALL setupConvDatafile(datafile=datafile_Sop, descr="Sop_entanglement", valstr="Sop")

  WRITE(*,*) " Starting iPEPS TEBD "

  !!! TIME EVOLUTION -- THE MAIN LOOP (using 1st order Trotter Decomposition) !!!
  mainloop:DO step = 1, (TEBD % N_steps)

    !! WRITE(*,*)
    !! WRITE(*,*) "TIMESTEP = ", step
    !! WRITE(*,*)

    !! (1,2) Evolve iPEPS bonds
    CALL mult_ipeps_2body_mpo(iGamma, iLambda, prop, TEBD % chi, TEBD % eps,  'S')
    CALL mult_ipeps_2body_mpo(iGamma, iLambda, prop, TEBD % chi, TEBD % eps,  'W')
    CALL mult_ipeps_2body_mpo(iGamma, iLambda, prop, TEBD % chi, TEBD % eps,  'N')
    CALL mult_ipeps_2body_mpo(iGamma, iLambda, prop, TEBD % chi, TEBD % eps,  'E')

    !! (2) Rescale iPEPS sites
    iGamma(1) % m = iGamma(1) % m / MAXVAL(ABS(iGamma(1) % m))
    iGamma(2) % m = iGamma(2) % m / MAXVAL(ABS(iGamma(2) % m))

    !! (3) Increment time by a single step after applying propagator to all bonds
    time = step * (TEBD % dt) + time_offset

    !! (4) Test convergence of TEBD iLambda spectrum
    IF((TEBD % test_interval .EQ. 1) .OR. (MOD(step, TEBD % test_interval) .EQ. 0)) THEN 
        CALL test_lambda_convergence(is_converged, iLambda, iLambdaOLD, TEBD, datafile, step, ABS(time))
    end IF

    !! (5) If converged --> exit TEBD loop
    IF(is_converged .AND. (TEBD % exit_level .GE. 1)) THEN
       EXIT mainloop
    ELSEIF((step .EQ. TEBD % N_steps) .AND. (TEBD % exit_level .GE. 2)) THEN
       WRITE(*,*) "iPEPS TEBD HAS FAILED TO CONVERGE "
       STOP
    end IF

  end DO mainloop

  CLOSE(datafile_Sop)
  WRITE(datafile, *)
  WRITE(datafile, *)
  CLOSE(datafile)
  WRITE(*,*)
  WRITE(*,*) "FINISHED TEBD"
  WRITE(*,*)

 END SUBROUTINE time_evolve_simple_MPO




 !!! TEBD simple update -- but with environment && observables used in convergence testing !!!
 SUBROUTINE time_evolve_simple_with_obs(iGamma, iLambda, Cmat, Tmat, prop_dense, TEBD, CTMRG, use_old_ctm)

  TYPE(block_peps),        ALLOCATABLE, INTENT(INOUT) :: iGamma(:)        !iGamma of iPEPS
  TYPE(block_lambda),      ALLOCATABLE, INTENT(INOUT) :: iLambda(:)       !iLambda of iPEPS
  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)          !CTM Cmat
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)          !CTM Tmat
  COMPLEX(KIND=DP),                     INTENT(IN)    :: prop_dense(:,:)  !Propagator
  TYPE(tebd_params),                    INTENT(IN)    :: TEBD             !TEBD params
  TYPE(ctm_params),                     INTENT(IN)    :: CTMRG            !CTM params
  LOGICAL,                              INTENT(IN)    :: use_old_ctm      !Whether to use old CTM for initialization

  !! TEBD iteration vars
  INTEGER          :: step
  COMPLEX(KIND=DP) :: time, time_offset
  INTEGER          :: datafile
  LOGICAL          :: is_converged
  LOGICAL          :: use_old_env

  !! Reduced rho && observables
  COMPLEX(KIND=DP),   ALLOCATABLE :: rho_onesite(:,:),  rho_twosite(:,:,:,:,:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: new_obs_1P(:),     new_obs_2P(:,:)
  COMPLEX(KIND=DP),   ALLOCATABLE :: old_obs_1P(:),     old_obs_2P(:,:)

  !! Norm const
  COMPLEX(KIND=DP) :: C_norm

  !! Initialize time && convergence tester
  time_offset = 0.0D0
  time = 0.0D0
  is_converged = .FALSE.

  WRITE(*,*) " Starting iPEPS TEBD "

  !!! TIME EVOLUTION -- THE MAIN LOOP (using 1st order Trotter Decomposition) !!!
  mainloop:DO step = 1, (TEBD % N_steps)

    !WRITE(*,*)
    !WRITE(*,*) "TIMESTEP = ", step
    !WRITE(*,*)

    !! (1) Evolve iPEPS bonds (--- apply propagators in a symmetrized fashion to ensure stability -- is it really the right order???)
    CALL mult_ipeps_dense_prop(iGamma, iLambda, prop_dense, TEBD % chi, TEBD % eps,  'S')
    CALL mult_ipeps_dense_prop(iGamma, iLambda, prop_dense, TEBD % chi, TEBD % eps,  'W')
    CALL mult_ipeps_dense_prop(iGamma, iLambda, prop_dense, TEBD % chi, TEBD % eps,  'N')
    CALL mult_ipeps_dense_prop(iGamma, iLambda, prop_dense, TEBD % chi, TEBD % eps,  'E')

    !! (2) Rescale iPEPS sites
    iGamma(1) % m = iGamma(1) % m / matrix_frobenius_norm(TensTRANSPOSE(RESHAPE_2D(iGamma(1) % m, '1,2345')))
    iGamma(2) % m = iGamma(2) % m / matrix_frobenius_norm(TensTRANSPOSE(RESHAPE_2D(iGamma(2) % m, '1,2345')))

    !! (3) Increment time by a single step after applying propagator to all bonds 
    time = step * (TEBD % dt) + time_offset

    !! (4) Compute environment && test convergence 
    IF((TEBD % test_interval .EQ. 1) .OR. (MOD(step, TEBD % test_interval) .EQ. 0)) THEN 

        !! Whether to reuse old environment
        use_old_env = use_old_ctm .OR. (step .GT. TEBD % test_interval)

        !! Compute iPEPS environment
        CALL calc_CTM_environment_ipeps(Cmat, Tmat, iGamma, iLambda, CTMRG, use_old_ctm=use_old_env)

        !! Construct reduced RHO for computing observables
        CALL construct_reduced_rho(rho_onesite, rho_twosite, Cmat, Tmat, iGamma, iLambda, CTMRG, TEBD % max_sep)

        !! Test convergence (using both 1P observables && 2P correlations)
        CALL test_tebd_convergence(is_converged, new_obs_2P, old_obs_2P, new_obs_1P, old_obs_1P, &
                                                 & rho_onesite, rho_twosite, TEBD, step, ABS(time), use_old_env)
    end IF

    !! (5) If converged --> exit TEBD loop
    IF(is_converged .AND. (TEBD % exit_level .GE. 1)) THEN
       EXIT mainloop
    ELSEIF((step .EQ. TEBD % N_steps) .AND. (TEBD % exit_level .GE. 2)) THEN
       WRITE(*,*) "iPEPS TEBD HAS FAILED TO CONVERGE "
       STOP
    end IF

  end DO mainloop

  WRITE(*,*)
  WRITE(*,*) "FINISHED TEBD"
  WRITE(*,*)
 
 END SUBROUTINE time_evolve_simple_with_obs

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Auxiliary routines for TEBD time evolution  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Calculate summary observables and output SS (steady state) !!!
 SUBROUTINE calc_tebd_summary_GamLam(iGamma, iLambda, Cmat, Tmat, CTMRG, x, min_sep, max_sep, g2_index, dump_parameters)

  USE datafile_utility

  TYPE(block_peps),        INTENT(IN)           :: iGamma(:)
  TYPE(block_lambda),      INTENT(IN)           :: iLambda(:)
  TYPE(ctm_corner_type),   INTENT(IN)           :: Cmat(:)    
  TYPE(ctm_transfer_type), INTENT(IN)           :: Tmat(:)  
  TYPE(ctm_params),        INTENT(IN)           :: CTMRG        
  REAL(KIND=DP),           INTENT(IN)           :: x
  INTEGER,                 INTENT(IN)           :: min_sep, max_sep
  INTEGER,                 INTENT(IN), OPTIONAL :: g2_index

  INTERFACE
     SUBROUTINE dump_parameters(HANDLE)
         INTEGER, INTENT(IN) :: HANDLE
     END SUBROUTINE dump_parameters
  end INTERFACE  

  !Reduced rho && norm of rho TN (i.e. trace)
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_onesite_s1(:,:), rho_onesite_s2(:,:), rho_twosite(:,:,:,:,:)
  COMPLEX(KIND=DP)              :: norm

  !Observables (temp storage)
  COMPLEX(KIND=DP)              :: expval_1P_avg
  COMPLEX(KIND=DP)              :: expval_1P_s1  !1P obs
  COMPLEX(KIND=DP)              :: expval_1P_s2  !1P obs
  COMPLEX(KIND=DP), ALLOCATABLE :: expvals_2P(:) !2P obs
  COMPLEX(KIND=DP)              :: g2val, nval   !store G2 in g2val, store nval to normalize G2 (if requested)

  !string vars for SS output filename
  CHARACTER :: SS_state_filename*256, CTM_filename*256, char_chi*16, char_chiENV*16, char_x*16

  !indices && dims
  INTEGER :: unit, i, sep

  !! Construct reduced rho for computing SS observables
  CALL construct_reduced_rho(rho_onesite_s1, rho_twosite, Cmat, Tmat, iGamma, iLambda, CTMRG, max_sep=0)
  CALL construct_reduced_rho(rho_onesite_s2, rho_twosite, Cmat, Tmat, iGamma, iLambda, CTMRG, max_sep=0)

  !! Get norm of rho TN (i.e. trace)
  norm = TRACE(rho_onesite_s1)

  !! Single site expectations (divide by norm)
  DO i=1,num_1p_summary

     CALL compute_obs_1P(expval_1P_s1, rho_onesite_s1, SUMM_OP_1p(i,:,:))
     CALL compute_obs_1P(expval_1P_s2, rho_onesite_s2, SUMM_OP_1p(i,:,:))

     expval_1P_avg = CMPLX(SIGN(0.5D0 * (ABS(expval_1P_s1) + ABS(expval_1P_s2)), REAL(expval_1P_s1)))

     CALL write_summary_data_1P(SUMM_OUT_1p(i), x, expval_1P_avg/norm, norm)
  end DO

  !! Two-site expectations (divide by norm) 
  !DO i=1,num_2p_summary

     !! Compute observables
     !CALL allocateTens(expvals_2P, (/max_sep+1-min_sep/)) !(allocate storage for summary (SS) observables)
     !CALL compute_obs_2P(expvals_2P, rho_twosite, rho_onesite_s1, SUMM_OP_2p(i,:,:,:), min_sep, max_sep)

     !! Normalize observables
     !DO sep=min_sep,max_sep

        !! Normalize 
        !expvals_2P(sep+1-min_sep) = expvals_2P(sep+1-min_sep)/(norm**(sep/2+1))
        !CALL write_summary_data_2P(SUMM_OUT_2p(i), x, sep, expvals_2P(sep+1-min_sep), norm)
     !end DO
  !end DO

  !! Add blank lines to all two-point summary files.
  !DO i=1,num_2p_summary
     !WRITE(SUMM_OUT_2p(i),*)
  !end DO

  !! Dump the converged SS and its CTM Environment
  WRITE(char_x,   '(F6.3)') x
  WRITE(char_chi,   '(I5)') iGamma(1) % Wdim 
  WRITE(char_chiENV,'(I5)') SIZE(Cmat(1) % m, 1)

  SS_state_filename="steady_state_x="//TRIM(ADJUSTL(char_x))//"_chi="//TRIM(ADJUSTL(char_chi))
  CALL peps_print_function(iGamma, iLambda, SS_state_filename)

  CTM_filename="CTM_x="//TRIM(ADJUSTL(char_x))//"_chi="//TRIM(ADJUSTL(char_chi))//"_chiENV="//TRIM(ADJUSTL(char_chiENV))
  CALL CTM_print_function(Cmat, Tmat, CTM_filename)

  !! Write params to all trace files (as a footer), so we know when parameters changed
  DO UNIT=FIRST_HANDLE, DUMP_NUM
     CALL dump_parameters(UNIT)
  end DO

  !! Print iPEPS bonds
  CALL print_ipeps_bonds(iGamma, iLambda)

 END SUBROUTINE calc_tebd_summary_GamLam







 !!! Calculate summary observables and output SS (steady state) !!!
 SUBROUTINE calc_tebd_summary_PLAIN(ipeps, Cmat, Tmat, CTMRG, x, min_sep, max_sep, g2_index, dump_parameters)

  USE datafile_utility

  TYPE(block_peps),        INTENT(IN)           :: ipeps(:)
  TYPE(ctm_corner_type),   INTENT(IN)           :: Cmat(:)    
  TYPE(ctm_transfer_type), INTENT(IN)           :: Tmat(:)  
  TYPE(ctm_params),        INTENT(IN)           :: CTMRG        
  REAL(KIND=DP),           INTENT(IN)           :: x
  INTEGER,                 INTENT(IN)           :: min_sep, max_sep
  INTEGER,                 INTENT(IN), OPTIONAL :: g2_index

  INTERFACE
     SUBROUTINE dump_parameters(HANDLE)
         INTEGER, INTENT(IN) :: HANDLE
     END SUBROUTINE dump_parameters
  end INTERFACE  

  !Reduced rho && norm of rho TN (i.e. trace)
  COMPLEX(KIND=DP), ALLOCATABLE :: rho_onesite(:,:), rho_twosite(:,:,:,:,:)
  COMPLEX(KIND=DP)              :: norm

  !Observables (temp storage)
  COMPLEX(KIND=DP)              :: expval_1P     !1P obs
  COMPLEX(KIND=DP), ALLOCATABLE :: expvals_2P(:) !2P obs
  COMPLEX(KIND=DP)              :: g2val, nval   !store G2 in g2val, store nval to normalize G2 (if requested)

  !string vars for SS output filename
  CHARACTER :: SS_state_filename*256, CTM_filename*256, char_chi*16, char_chiENV*16, char_x*16

  !indices && dims
  INTEGER :: unit, i, sep

  !! Construct reduced rho for computing SS observables
  CALL construct_reduced_rho(rho_onesite, rho_twosite, Cmat, Tmat, ipeps, CTMRG, max_sep=0)

  !! Get norm of rho TN (i.e. trace)
  norm = TRACE(rho_onesite)

  !! Single site expectations (divide by norm)
  DO i=1,num_1p_summary

     CALL compute_obs_1P(expval_1P, rho_onesite, SUMM_OP_1p(i,:,:))
     expval_1P = expval_1P/norm
     CALL write_summary_data_1P(SUMM_OUT_1p(i), x, expval_1P, norm)
     
     !Record expectation of number for post-processing below
     IF(PRESENT(g2_index)) THEN
         IF(i .EQ. g2_index) nval=expval_1P 
     end IF
  end DO

  !! Two-site expectations (divide by norm) 
  !DO i=1,num_2p_summary

     !Compute observables
     !CALL allocateTens(expvals_2P, setShape(max_sep+1-min_sep)) !(allocate storage for summary (SS) observables)
     !CALL compute_obs_2P(expvals_2P, rho_twosite, rho_onesite, SUMM_OP_2p(i,:,:,:), min_sep, max_sep)

     !Normalize observables
     !DO sep=min_sep,max_sep

        !Normalize 
        !expvals_2P(sep+1-min_sep) = expvals_2P(sep+1-min_sep) / (norm**(sep/2+1))
        !CALL write_summary_data_2P(SUMM_OUT_2p(i), x, sep, expvals_2P(sep+1-min_sep), norm)
        
        !Extra summary file, G2 post processing. This should give <a^dagger_i a^dagger_j a_j a_i>
        !IF(PRESENT(g2_index)) THEN
           !IF(i .EQ. g2_index) THEN
              !g2val=expvals_2P(sep+1-min_sep); IF(sep .EQ. 0) g2val = g2val - nval
              !CALL write_summary_data_2P(G2_SUMMARY, x, sep, g2val, nval)
           !end IF
        !end IF
     !end DO
  !end DO

  !! Add blank lines to all two-point summary files.
  !DO i=1,num_2p_summary
     !WRITE(SUMM_OUT_2p(i),*)
  !end DO

  !! Dump the converged SS and its CTM Environment
  WRITE(char_x,   '(F6.3)') x
  WRITE(char_chi,   '(I5)') ipeps(1) % Wdim 
  WRITE(char_chiENV,'(I5)') SIZE(Cmat(1) % m, 1)

  SS_state_filename="FU_steady_state_x="//TRIM(ADJUSTL(char_x))//"_chi="//TRIM(ADJUSTL(char_chi))
  CALL peps_print_function(ipeps, SS_state_filename)

  CTM_filename="CTM_x="//TRIM(ADJUSTL(char_x))//"_chi="//TRIM(ADJUSTL(char_chi))//"_chiENV="//TRIM(ADJUSTL(char_chiENV))
  CALL CTM_print_function(Cmat, Tmat, CTM_filename)

  !! Write params to all trace files (as a footer), so we know when parameters changed
  DO UNIT=FIRST_HANDLE, DUMP_NUM
     CALL dump_parameters(UNIT)
  end DO

 END SUBROUTINE calc_tebd_summary_PLAIN







 !!! Output converged iPEPS -- Gamma-Lambda form !!!
 SUBROUTINE output_converged_ipeps_GamLam(iGamma, iLambda, Cmat, Tmat, x)

  TYPE(block_peps),        INTENT(IN)           :: iGamma(:)
  TYPE(block_lambda),      INTENT(IN)           :: iLambda(:)
  TYPE(ctm_corner_type),   INTENT(IN)           :: Cmat(:)    
  TYPE(ctm_transfer_type), INTENT(IN)           :: Tmat(:)      
  REAL(KIND=DP),           INTENT(IN)           :: x

  !! String vars for SS output filename
  CHARACTER :: SS_state_filename*256, CTM_filename*256, char_chi*16, char_chiENV*16, char_x*16

  !! Set string vars
  WRITE(char_x,   '(F6.3)') x
  WRITE(char_chi,   '(I5)') iGamma(1) % Wdim 
  WRITE(char_chiENV,'(I5)') SIZE(Cmat(1) % m, 1)

  !! Output converged iPEPS
  SS_state_filename="steady_state_x="//TRIM(ADJUSTL(char_x))//"_chi="//TRIM(ADJUSTL(char_chi))
  CALL peps_print_function(iGamma, iLambda, SS_state_filename)

  !! Output CTM environment
  CTM_filename="CTM_x="//TRIM(ADJUSTL(char_x))//"_chi="//TRIM(ADJUSTL(char_chi))//"_chiENV="//TRIM(ADJUSTL(char_chiENV))
  CALL CTM_print_function(Cmat, Tmat, CTM_filename)

 END SUBROUTINE output_converged_ipeps_GamLam





 !!! Output converged iPEPS -- plain (no-lambda) form !!!
 SUBROUTINE output_converged_ipeps_PLAIN(ipeps, Cmat, Tmat, x)

  TYPE(block_peps),        INTENT(IN)           :: ipeps(:)
  TYPE(ctm_corner_type),   INTENT(IN)           :: Cmat(:)    
  TYPE(ctm_transfer_type), INTENT(IN)           :: Tmat(:)      
  REAL(KIND=DP),           INTENT(IN)           :: x

  !! String vars for SS output filename
  CHARACTER :: SS_state_filename*256, CTM_filename*256, char_chi*16, char_chiENV*16, char_x*16

  !! Set string vars
  WRITE(char_x,   '(F6.3)') x
  WRITE(char_chi,   '(I5)') ipeps(1) % Wdim 
  WRITE(char_chiENV,'(I5)') SIZE(Cmat(1) % m, 1)

  !! Output converged iPEPS
  SS_state_filename="steady_state_x="//TRIM(ADJUSTL(char_x))//"_chi="//TRIM(ADJUSTL(char_chi))
  CALL peps_print_function(ipeps, SS_state_filename)

  !! Output CTM environment
  CTM_filename="CTM_x="//TRIM(ADJUSTL(char_x))//"_chi="//TRIM(ADJUSTL(char_chi))//"_chiENV="//TRIM(ADJUSTL(char_chiENV))
  CALL CTM_print_function(Cmat, Tmat, CTM_filename)

 END SUBROUTINE output_converged_ipeps_PLAIN





 !!! Print iPEPS bond dims and Lambdas !!!
 SUBROUTINE print_ipeps_bonds(iGamma, iLambda)

  TYPE(block_peps),   INTENT(IN) :: iGamma(:)
  TYPE(block_lambda), INTENT(IN) :: iLambda(:)

  INTEGER :: i, site

  WRITE(*,*) "Gamma(1,2) dims ",  shape(iGamma(1)%m),  shape(iGamma(2)%m)
  WRITE(*,*) "Lambda(1,2) dims ", shape(iLambda(1)%m), shape(iLambda(2)%m)
  WRITE(*,*) "Lambda(3,4) dims ", shape(iLambda(3)%m), shape(iLambda(4)%m)

  WRITE(*,*) 
  DO i=1,iLambda(1)%ChiDim
     WRITE(*,*) "Lambda(1) -- ", iLambda(1)%m(i), " at i = ", i
  end DO
  WRITE(*,*)
  DO i=1,iLambda(2)%ChiDim
     WRITE(*,*) "Lambda(2) -- ", iLambda(2)%m(i), " at i = ", i
  end DO
  WRITE(*,*)
  DO i=1,iLambda(3)%ChiDim
     WRITE(*,*) "Lambda(3) -- ", iLambda(3)%m(i), " at i = ", i
  end DO
  WRITE(*,*)
  DO i=1,iLambda(4)%ChiDim
     WRITE(*,*) "Lambda(4) -- ", iLambda(4)%m(i), " at i = ", i
  end DO
  WRITE(*,*)

 END SUBROUTINE print_ipeps_bonds

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














 !!!!!!!!!!!!!!!!!!!! Initialization routines, specifically for density matrix TEBD (state_parameters must be available) !!!!!!!!!!!!!!!

 !!! Initialize rho PEPS state !!!
 SUBROUTINE initialize_rho_iPEPS(iGamma, iLambda, local_dim, chi, rho_vec, filename)

  USE iteration_helper
  !USE simulation_parameters

  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT)           :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT)           :: iLambda(:)
  INTEGER,                         INTENT(IN)              :: local_dim
  INTEGER,                         INTENT(IN),  OPTIONAL   :: chi
  COMPLEX(KIND=DP),                INTENT(IN),  OPTIONAL   :: rho_vec(:)
  CHARACTER(LEN=*),                INTENT(IN),  OPTIONAL   :: filename

  !! indices && dims
  INTEGER :: site

  !! size of Gamma && Lambda blocks
  INTEGER, PARAMETER :: N_sites = 2
  INTEGER, PARAMETER :: N_bonds = 4

  !! Local dim defaults to -1 (see simulation_parameters.f90) so we can detect if we did not set it yet.
  IF(local_dim .EQ. -1) THEN
     WRITE(*,*) "Allocate state called before basis state setup"
     WRITE(*,*) "Must call routine setting up basis to determine local_dim"
     STOP
  end IF

  IF(PRESENT(filename)) THEN

     !! Reload PEPS rho state
     CALL reload_peps(filename, local_dim, iGamma, iLambda)

  ELSEIF(PRESENT(rho_vec) .AND. (.NOT. PRESENT(chi))) THEN

     !! Initialize factorized iPEPS with pre-determined rho_vec
     CALL initialize_factorized_peps(iGamma, local_dim, N_sites, rho_vec) 
     CALL initialize_factorized_lambda(iLambda, N_bonds) 
  ELSE

     !! Initialize random iPEPS
     CALL initialize_rand_full_ipeps(iGamma, iLambda, local_dim, OptArg(chi, 1))
  end IF

 END SUBROUTINE initialize_rho_iPEPS




 !!! Reload CTM environment of rho iPEPS !!!
 SUBROUTINE reload_rho_ipeps_CTM(Cmat, Tmat, iGamma, iLambda, filename)

  TYPE(ctm_corner_type),   ALLOCATABLE, INTENT(INOUT) :: Cmat(:)
  TYPE(ctm_transfer_type), ALLOCATABLE, INTENT(INOUT) :: Tmat(:)
  TYPE(block_peps),                     INTENT(IN)    :: iGamma(:)
  TYPE(block_lambda),                   INTENT(IN)    :: iLambda(:)
  CHARACTER,                            INTENT(IN)    :: filename*256

  !! Reduced TN
  TYPE(block_peps), ALLOCATABLE :: ipeps(:)
  TYPE(block_mpo),  ALLOCATABLE :: TN2D(:) 

  !! Reload CTM tensors from file
  CALL reload_ctm(filename, Cmat, Tmat)

  !! Check CTM dims (in particular, check if CTM local dims match iPEPS bond dims)
  CALL absorb_ipeps_lambdas(ipeps, iGamma, iLambda)
  CALL create_reduced_TN_from_ipeps(TN2D, ipeps)

  CALL check_ctm_dims(Cmat, Tmat, TN2D, DIR='S')
  CALL check_ctm_dims(Cmat, Tmat, TN2D, DIR='N')
  CALL check_ctm_dims(Cmat, Tmat, TN2D, DIR='W')
  CALL check_ctm_dims(Cmat, Tmat, TN2D, DIR='E')

 END SUBROUTINE reload_rho_ipeps_CTM



 !!! Calculate operator entanglement entropy of 2 x 2 PEPS block with PBC !!!
 SUBROUTINE calc_op_entropy(Sop, iLambda)

  COMPLEX(KIND=DP),   INTENT(INOUT) :: Sop
  TYPE(block_lambda), INTENT(IN)    :: iLambda(:) 

  Sop = 4.0D0*Sop_bound(iLambda(2)) + 4.0D0*Sop_bound(iLambda(4))

  !Sop = 0.5D0*(4.0D0*Sop_bound(iLambda(2)) + 4.0D0*Sop_bound(iLambda(4)) + 4.0D0*Sop_bound(iLambda(1)) + 4.0D0*Sop_bound(iLambda(3)))

 END SUBROUTINE calc_op_entropy




 !!! Calculate a single boundary contribution to operator entanglement entropy of 2 x 2 PEPS block with PBC !!!
 FUNCTION Sop_bound(iLambda) result(Sop)

  TYPE(block_lambda), INTENT(IN) :: iLambda 
  COMPLEX(KIND=DP)               :: Sop

  COMPLEX(KIND=DP) :: temp
  INTEGER :: i, chi

  !! Get bond dim of iLambda, allocate temp storages
  chi = SIZE(iLambda % m) 

  !! Compute Sop contribution from a single boundary
  Sop = 0.0D0
  DO i=1,chi
     temp = (iLambda % m(i))**2
     Sop  = Sop - temp * logN(temp, 2.0D0)
  end DO

 END FUNCTION Sop_bound

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE time_evolution_ipeps
