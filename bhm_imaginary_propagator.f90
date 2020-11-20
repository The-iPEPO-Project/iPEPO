MODULE bhm_imaginary_propagator

  USE utility
  USE definitions_mps_mpo
  USE psi_sun_basis
  USE propagator_utility

  IMPLICIT NONE

  !!!!! Hamiltonian Propagator for BHM !!!!!

CONTAINS

  !!! Setup BHM propagator !!!
  SUBROUTINE setup_propagator(propagator_mpo, g, U, J, dt, mps_vs_peps, xi)
    
    USE simulation_parameters, ONLY: local_dim !Use system parameters   

    TYPE(block_mpo),  ALLOCATABLE, INTENT(OUT)           :: propagator_mpo(:)
    REAL(KIND=DP),                 INTENT(IN)            :: g, U, J
    COMPLEX(KIND=DP),              INTENT(IN)            :: dt
    CHARACTER(LEN=*),              INTENT(IN),  OPTIONAL :: mps_vs_peps
    REAL(KIND=DP),                 INTENT(IN),  OPTIONAL :: xi

    !! Hamiltonian sites
    COMPLEX(KIND=DP) :: super_H_site(local_dim, local_dim)
    COMPLEX(KIND=DP) :: H_pair(local_dim**2, local_dim**2)

    !! Dense propagator
    COMPLEX(KIND=DP), ALLOCATABLE :: prop_2D(:,:), prop_4D(:,:,:,:)

    !! (1) Setup onsite term. 
    !!     NB, all terms in our Hamiltonian would be negative if not for exp(-H*tau) 
    super_H_site = 0.0D0
    super_H_site = g * number_op - U * g2bare_op

    !! (2) Setup two-body term.
    H_pair = 0.0D0
    H_pair = J * TensKRON(up_op, dn_op) + J * TensKRON(dn_op, up_op)

    !! (3) Set up dense propagator given all sites are identical.
    CALL setup_ti_propagator_dense(prop_2D, super_H_site, super_H_site, H_pair, dt, mps_vs_peps)

    !! (4) Reshape into 4D
    prop_4D = RESHAPE_4D(prop_2D, '12,34', (/ local_dim, local_dim /), (/ local_dim, local_dim /))

    !! (5) Construct MPO propagator
    CALL canonicalize_propagator(propagator_mpo, prop_4D, xi)

  END SUBROUTINE setup_propagator




  !!! Setup BHM propagator !!!
  SUBROUTINE setup_propagator_dense(prop, g, U, J, dt, mps_vs_peps)
    
    USE simulation_parameters, ONLY: local_dim !Use system parameters   

    COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT) :: prop(:,:)
    REAL(KIND=DP),                 INTENT(IN)  :: g, U, J
    COMPLEX(KIND=DP),              INTENT(IN)  :: dt
    CHARACTER(LEN=*),              INTENT(IN)  :: mps_vs_peps

    COMPLEX(KIND=DP) :: super_H_site(local_dim, local_dim)
    COMPLEX(KIND=DP) :: H_pair(local_dim**2, local_dim**2)

    !! (1) Setup onsite term. 
    !!     NB, all terms in our Hamiltonian would be negative if not for exp(-H*tau) 
    super_H_site = 0.0D0
    super_H_site = g * number_op - U * g2bare_op

    !! (2) Setup two-body term.
    H_pair = 0.0D0
    H_pair = J * TensKRON(up_op, dn_op) + J * TensKRON(dn_op, up_op)

    !! (3) Set up propagators given all sites are identical.
    CALL setup_ti_propagator_dense(prop, super_H_site, super_H_site, H_pair, dt, mps_vs_peps)

  END SUBROUTINE setup_propagator_dense





  !!! Compute two-site Hamiltonian !!!
  SUBROUTINE setup_hamiltonian_twosite(hamiltonian_twosite, g, U, J, mps_vs_peps)

    USE simulation_parameters, ONLY: local_dim !! Use system parameters  

    COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT)          :: hamiltonian_twosite(:,:)
    REAL(KIND=DP),                 INTENT(IN)           :: g, U, J
    CHARACTER(LEN=*),              INTENT(IN), OPTIONAL :: mps_vs_peps

    COMPLEX(KIND=DP) :: super_H_site(local_dim, local_dim)
    COMPLEX(KIND=DP) :: H_pair(local_dim**2, local_dim**2)
    REAL(KIND=DP)    :: lfactor, rfactor

    !! (1) Setup onsite term. 
    super_H_site = 0.0D0
    super_H_site = -g * number_op + U * g2bare_op

    !! (2) Setup two-body term.
    H_pair = 0.0D0
    H_pair = -J * TensKRON(up_op, dn_op) - J * TensKRON(dn_op, up_op)

    !! Determine prefactors based on dimensionality of the problem (i.e. MPS vs PEPS)
    SELECT CASE(OptArg(mps_vs_peps, 'MPS'))
    CASE('MPS')
       lfactor = 0.50D0; rfactor = 0.50D0 
    CASE('PEPS')
       lfactor = 0.25D0; rfactor = 0.25D0 
    CASE DEFAULT
       CALL invalid_flag("setup_ti_propagator_dense -- invalid mps_vs_peps ", OptArg(mps_vs_peps, 'MPS'))
    end SELECT 

    !! Combine onsite && interaction terms
    ALLOCATE(hamiltonian_twosite(local_dim**2, local_dim**2))
    hamiltonian_twosite = lfactor * TensKRON(super_H_site, eye) + rfactor * TensKRON(eye, super_H_site) + H_pair

  END SUBROUTINE setup_hamiltonian_twosite

END MODULE bhm_imaginary_propagator
