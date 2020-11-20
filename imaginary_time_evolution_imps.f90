MODULE imaginary_time_evolution_imps

 USE utility
 USE definitions_mps_mpo
 USE mps_peps_INOUT_files
 USE mps_mpo_algebra_inf

 USE simulation_parameters
 USE imaginary_observables_imps
 USE imaginary_tebd_callback_routines_imps
 
 IMPLICIT NONE

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEBD algorithm for evolving iMPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE time_evolve(iGamma, iLambda, mpo_prop, TEBD)

  USE datafile_utility, ONLY: setupConvDatafile

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)   !iGamma of MPS
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)  !iLambda of MPS
  TYPE(block_mpo),                 INTENT(IN)    :: mpo_prop(:) !Propagator
  TYPE(tebd_params),               INTENT(IN)    :: TEBD        !TEBD params

  !! TEBD iteration vars
  INTEGER          :: step
  COMPLEX(KIND=DP) :: time, time_offset
  INTEGER          :: datafile
  LOGICAL          :: is_converged

  !! Old iLambda for lambda convergence testing
  TYPE(block_lambda), ALLOCATABLE :: iLambdaOLD(:)

  !! Reduced rho && observables 
  COMPLEX(KIND=DP), ALLOCATABLE :: new_obs_1P(:),   old_obs_1P(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: new_obs_2P(:,:), old_obs_2P(:,:)

  !! Boundaries && eval of iMPS network
  COMPLEX(KIND=DP), ALLOCATABLE :: Rvec(:,:), Lvec(:,:)
  COMPLEX(KIND=DP)              :: eval
  
  !! Timings
  REAL(KIND=DP) :: cpu_times(4)

  !! Initialize time && convergence tester
  is_converged = .FALSE.
  time_offset = 0.0D0
  time = 0.0D0

  !! Setup Lambda convergence file, initialize old Lambda
  CALL copy_lambda_block(iLambdaOLD, iLambda)
  CALL setupConvDatafile(datafile=datafile, descr="iMPS_TEBD_lambda", errstr="Err_lambda", valstr="Max_lambda", xstr="dt")

  WRITE(*,*) " Starting iMPS Imaginary TEBD "

  !!! TIME EVOLUTION -- THE MAIN LOOP (using 2nd order Trotter Decomposition) !!!
  mainloop:DO step = 1, (TEBD % N_steps)

    CALL cpu_time(cpu_times(1))

    !! Evolve iMPS bonds 
    CALL mult_imps_2body_mpo(iGamma, iLambda, mpo_prop, TEBD % chi, TEBD % eps, '12')
    CALL mult_imps_2body_mpo(iGamma, iLambda, mpo_prop, TEBD % chi, TEBD % eps, '21')
    CALL mult_imps_2body_mpo(iGamma, iLambda, mpo_prop, TEBD % chi, TEBD % eps, '21')
    CALL mult_imps_2body_mpo(iGamma, iLambda, mpo_prop, TEBD % chi, TEBD % eps, '12')

    CALL cpu_time(cpu_times(2))

    !! Increment time by a single step after applying propagator to all bonds
    time = step * (TEBD % dt) + time_offset

    !! Test convergence (using single-site observables)
    IF((TEBD % test_interval .EQ. 1) .OR. (MOD(step, TEBD % test_interval) .EQ. 0)) THEN

        WRITE(*,*) "Finished TIMESTEP = ", step, " Bond dims =  ", iLambda(1) % ChiDim, iLambda(2) % ChiDim

        !! Calc iMPS--iMPS_HC bounds 
        CALL compute_psi_imps_bounds(eval, Rvec, Lvec, iGamma, iLambda)

        CALL cpu_time(cpu_times(3))

        !! Test convergence && print observables
        CALL test_lambda_convergence(is_converged, iLambda, iLambdaOLD, TEBD, datafile, step, ABS(time))
        !CALL test_tebd_convergence(is_converged, new_obs_2P, old_obs_2P, new_obs_1P, old_obs_1P, &
                                   !& iGamma, iLambda, Rvec, Lvec, eval, TEBD, step, ABS(time))

        !! Normalize iGamma by eval = <Psi|Psi> norm constant
        CALL mult_mps_by_value(iGamma, eval**(-0.5d0))

        CALL cpu_time(cpu_times(4))

        !! If converged --> exit TEBD loop
        IF(is_converged) THEN  
              CALL equalize_imps_bonds(iGamma, iLambda) !! Equalize iMPS bonds (we do it here not to interfere with convergence)
              EXIT mainloop
        ELSEIF(step .EQ. TEBD % N_steps) THEN
              WRITE(*,*) "TEBD: SS HAS FAILED TO CONVERGE"
              STOP
        end IF
    end IF

  end DO mainloop

  CLOSE(datafile)
 
 END SUBROUTINE time_evolve

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Auxiliary routines for TEBD time evolution  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Calculate summary observables and output SS (steady state) !!!
 SUBROUTINE calc_tebd_summary(iGamma, iLambda, x, min_sep, max_sep, g2_index, dump_parameters)

  USE datafile_utility

  TYPE(block_mps),    INTENT(IN) :: iGamma(:)
  TYPE(block_lambda), INTENT(IN) :: iLambda(:)
  REAL(KIND=DP),      INTENT(IN) :: x
  INTEGER,            INTENT(IN) :: min_sep, max_sep
  INTEGER,            INTENT(IN), OPTIONAL :: g2_index

  INTERFACE
     SUBROUTINE dump_params(HANDLE)
         INTEGER, INTENT(IN) :: HANDLE
     END SUBROUTINE dump_params
  end INTERFACE  

  !! Reduced rho && boundaries for computing observables in infinite TN
  COMPLEX(KIND=DP), ALLOCATABLE :: Rvec(:,:), Lvec(:,:)
  COMPLEX(KIND=DP)              :: eval

  !! Observables (temp storage)
  COMPLEX(KIND=DP)              :: expval_1P     !1P obs
  COMPLEX(KIND=DP), ALLOCATABLE :: expvals_2P(:) !2P obs
  COMPLEX(KIND=DP)              :: g2val, nval   !store G2 in g2val, store nval to normalize G2 (if requested)

  !! String vars for SS output filename
  CHARACTER :: SS_state_filename*256, char_chi*16, char_x*16

  !! Indices && dims
  INTEGER :: unit, i, sep

  !! Construct reduced rho for computing SS observables
  CALL compute_psi_imps_bounds(eval, Rvec, Lvec, iGamma, iLambda)

  !! Single site expectations (renormalize by eval)
  DO i=1,num_1p_summary

     CALL compute_obs_1P(expval_1P, SUMM_OP_1p(i,:,:), iGamma, iLambda, Rvec, Lvec)
     expval_1P = expval_1P/eval
     CALL write_summary_data_1P(SUMM_OUT_1p(i), x, expval_1P, eval)
     
     !! Record expectation of number for post-processing below
     IF(PRESENT(g2_index)) THEN
         IF(i .EQ. g2_index) nval=expval_1P 
     end IF
  end DO


  !! Two-site expectations (renormalize by eval) 
  DO i=1,num_2p_summary

     !! Compute observables
     CALL allocateTens(expvals_2P,  (/max_sep+1-min_sep/)) !(allocate storage for summary (SS) observables)
     CALL compute_obs_2P(expvals_2P, SUMM_OP_2p(i,:,:,:), iGamma, iLambda, Rvec, Lvec, min_sep, max_sep)

     !! Normalize observables
     DO sep=min_sep,max_sep

        !! Normalize 
        expvals_2P(sep+1-min_sep) = expvals_2P(sep+1-min_sep)/(eval**(sep/2+1))
        CALL write_summary_data_2P(SUMM_OUT_2p(i), x, sep, expvals_2P(sep+1-min_sep), eval)
        
        !! Extra summary file, G2 post processing. This should give <a^dagger_i a^dagger_j a_j a_i>
        IF(PRESENT(g2_index)) THEN
           IF(i .EQ. g2_index) THEN
              g2val=expvals_2P(sep+1-min_sep); IF(sep .EQ. 0) g2val = g2val - nval
              CALL write_summary_data_2P(G2_SUMMARY, x, sep, g2val, nval)
           end IF
        end IF
     end DO
  end DO

  !! Add blank lines to all two-point summary files.
  DO i=1,num_2p_summary
     WRITE(SUMM_OUT_2p(i),*)
  end DO
  IF(PRESENT(g2_index)) WRITE(G2_SUMMARY,*)

  !! Dump the converged SS
  WRITE(char_x,'(F5.2)') x
  WRITE(char_chi,'(I5)') iLambda(1)%ChiDim

  SS_state_filename="steady_state_x="//TRIM(ADJUSTL(char_x))//"_chi="//TRIM(ADJUSTL(char_chi))
  CALL mps_print_function(iGamma, iLambda, SS_state_filename)

  !! Write params to all trace files (as a footer), so we know when parameters changed
  DO UNIT=FIRST_HANDLE, DUMP_NUM
     CALL dump_parameters(UNIT)
  end DO

  !! Print iMPS bonds
  CALL print_imps_bonds(iGamma, iLambda)

 END SUBROUTINE calc_tebd_summary




 !!! Print iMPS bond dims and Lambdas !!!
 SUBROUTINE print_imps_bonds(iGamma, iLambda)

  TYPE(block_mps),    INTENT(IN) :: iGamma(:)
  TYPE(block_lambda), INTENT(IN) :: iLambda(:)

  INTEGER :: i, site

  WRITE(*,*) "Gamma(1,2) dims ",  shape(iGamma(1)%m),  shape(iGamma(2)%m)
  WRITE(*,*) "Lambda(1,2) dims ", shape(iLambda(1)%m), shape(iLambda(2)%m)

  WRITE(*,*) 
  DO i=1,iLambda(1)%ChiDim
     WRITE(*,*) "Lambda(1) -- ", iLambda(1)%m(i), " at i = ", i
  end DO
  WRITE(*,*)
  DO i=1,iLambda(2)%ChiDim
     WRITE(*,*) "Lambda(2) -- ", iLambda(2)%m(i), " at i = ", i
  end DO
  WRITE(*,*)

 END SUBROUTINE print_imps_bonds

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





 !!!!!!!!!!!!!!!!!!!! Initialization routines, specifically for |Psi> imaginary TEBD (simulation_parameters must be available) !!!!!!!!!!!!!!!!!

 !!! Initialize Psi MPS state !!!
 SUBROUTINE initialize_psi_iMPS(iGamma, iLambda, psi_vec, filename)

  USE simulation_parameters
  USE iteration_helper

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  COMPLEX(KIND=DP),   INTENT(IN),  OPTIONAL      :: psi_vec(:)
  CHARACTER(LEN=*),   INTENT(IN),  OPTIONAL      :: filename

  TYPE(block_mpo), ALLOCATABLE :: dummy_mpo(:)
  INTEGER :: site

  !! size of Gamma && Lambda blocks
  INTEGER, PARAMETER :: N_sites = 2
  INTEGER, PARAMETER :: N_bonds = 2

  !! Local dim defaults to -1 (see state_parameters.f90) so we can detect if we did not set it yet.
  IF(local_dim .EQ. -1) THEN
     WRITE(*,*) "Allocate state called before basis state setup"
     WRITE(*,*) "Must call routine setting up basis to determine local_dim"
     STOP
  end IF 

  IF(PRESENT(filename)) THEN

     !! Construct a dummy mpo block for allocating mps sites with correct local dim
     CALL allocate_empty_mpo_block(dummy_mpo, N_sites)

     DO site=1,N_sites
        CALL allocate_mpo_site(dummy_mpo(site), local_dim, local_dim, 1, 1)
     end DO

     !! Reload MPS Psi state
     CALL reload_mps(filename, dummy_mpo, iGamma, iLambda)

  ELSEIF(PRESENT(psi_vec)) THEN

     !! init to product state (using input psi_vec)
     CALL initialize_factorized_mps(iGamma, local_dim, N_sites, psi_vec)
     CALL initialize_factorized_lambda(iLambda, N_bonds)
  ELSE

     !! init to product state
     CALL initialize_factorized_mps(iGamma, local_dim, N_sites)
     CALL initialize_factorized_lambda(iLambda, N_bonds)
  end IF

 END SUBROUTINE initialize_psi_iMPS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Converting |Psi> to pure state density matrix |Psi><Psi| !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Convert Psi iMPS to Pure Rho iMPS !!!
 SUBROUTINE convert_psi_to_pure_rho(rhoG, rhoL, psiG, psiL)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: rhoG(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: rhoL(:)
  TYPE(block_mps),    ALLOCATABLE, INTENT(IN)    :: psiG(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(IN)    :: psiL(:)

  !! Nolambda Psi MPS, rho MPO
  TYPE(block_mps),    ALLOCATABLE :: psiMps(:), psiMpsHC(:)
  TYPE(block_mpo),    ALLOCATABLE :: rhoMpo(:)

  WRITE(*,*) "Converting Psi to Rho... "

  !! (1) Create nolambda Psi iMPS and its HC copy 
  CALL absorb_imps_lambdas(psiMps, psiG, psiL)
  CALL copy_mps_block(psiMpsHC, psiMps, 'HC')

  !! (2) Tensor product of Psi iMPS and its HC --> gives rho iMPO
  CALL kronecker_mult_mps_mps(rhoMpo, psiMps, psiMpsHC)

  !! (3) Vectorize iMPO --> gives rho iMPS
  CALL reshape_mpo_to_mps(rhoG, rhoMpo)

  !! (4) Convert nolambda iMPS to GL (Vidal form) iMPS
  CALL convert_imps_to_GL_form(rhoG, rhoL, canonical_flag='CANONICAL')

  WRITE(*,*) "Finished Psi-Rho conversion "

 END SUBROUTINE convert_psi_to_pure_rho




 !!! Reshape MPO to MPS !!! FIXME OBSOLETE -- use new routine in mps_mpo_utility
 !SUBROUTINE reshape_mpo_to_mps(mps_block, mpo_block)

  !TYPE(block_mps), ALLOCATABLE, INTENT(INOUT) :: mps_block(:)
  !TYPE(block_mpo),              INTENT(IN)    :: mpo_block(:)

  !INTEGER :: site, N_sites

  !N_sites = SIZE(mpo_block)

  !CALL allocate_empty_mps_block(mps_block, N_sites)

  !! Reshape each MPO site into MPS site
  !DO site=1,N_sites
     !CALL allocate_mps_site(mps_block(site), RESHAPE_3D(mpo_block(site) % m, '12,3,4'))
  !end DO

 !END SUBROUTINE reshape_mpo_to_mps






 !!! Convert nolambda iMPS to GL (Vidal form) iMPS !!!
 SUBROUTINE convert_imps_to_GL_form(iGamma, iLambda, imps, canonical_flag)

  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  TYPE(block_mps),    INTENT(IN),  OPTIONAL      :: imps(:)
  CHARACTER(LEN=*),   INTENT(IN),  OPTIONAL      :: canonical_flag

  INTEGER :: site

  !! Allocate empty iLambda block
  CALL allocate_empty_lambda_block(iLambda, 2)

  !! Local copy of nolambda iMPS
  IF(PRESENT(imps)) CALL copy_mps_block(iGamma, imps)

  !! Convert nolambda iMPS to GL iMPS 
  IF(PRESENT(canonical_flag)) THEN

     !! Create initial iLambda = Id matrix
     DO site=1,2
        CALL allocate_lambda_site(iLambda(site), vecOnes(iGamma(site)%Edim))
     end DO

     !! Canonicalize iMPS
     CALL canonicalize_imps(iGamma, iLambda, chi=-1, eps=2.0d0)
  ELSE

     !! (1) SVD(12) --> Extract iLambda(2)
     CALL compute_svd_of_imps_bond(iGamma, iLambda(1), chi0=-1, eps0=2.0d0, SITEORD='12')

     !! (2) Absorb iLambda(1)
     CALL LambdaMUL(iGamma(1) % m, iLambda(1) % m, iGamma(2) % m)

     !! (3) SVD(21) --> Extract iLambda(2)
     CALL compute_svd_of_imps_bond(iGamma, iLambda(2), chi0=-1, eps0=2.0d0, SITEORD='21')

     !! (4) Restore iLambda(1)
     CALL LambdaDIV(iGamma(1) % m, SQRT(iLambda(1) % m), '3')
     CALL LambdaDIV(iGamma(2) % m, SQRT(iLambda(1) % m), '2')
  end IF

 END SUBROUTINE convert_imps_to_GL_form

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE imaginary_time_evolution_imps
