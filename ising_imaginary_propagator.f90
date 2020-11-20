MODULE ising_imaginary_propagator

  USE utility
  USE definitions_mps_mpo
  USE psi_pauli_basis
  USE propagator_utility

  IMPLICIT NONE

CONTAINS

  !!! setup Ising propagator !!!
  SUBROUTINE setup_propagator(propagator_mpo, g, delta, dt, mps_vs_peps)
    
    USE simulation_parameters, ONLY: local_dim !Use system parameters  

    TYPE(block_mpo),  ALLOCATABLE, INTENT(OUT) :: propagator_mpo(:)
    REAL(KIND=DP),                 INTENT(IN)  :: g, delta
    COMPLEX(KIND=DP),              INTENT(IN)  :: dt
    CHARACTER(LEN=*),              INTENT(IN)  :: mps_vs_peps

    COMPLEX(KIND=DP)             :: super_H_site(local_dim, local_dim)
    TYPE(block_mpo), ALLOCATABLE :: H_pair(:)
    INTEGER                      :: ocnt

    !! (1) Setup onsite term. 
    !! NB, all terms in our Hamiltonian would be negative if not for exp(-H*tau) 
    super_H_site = 0.0D0
    super_H_site = g * sigma_x

    !! (2) Setup two-body term.
    CALL allocate_pair_mpo(H_pair, local_dim, 1)

    ocnt = 0
    CALL mpo_add_operator(H_pair, ocnt, (1.0D0, 0.0D0), sigma_z, sigma_z)

    !! (3) Set up propagators given all sites are identical.
    CALL setup_ti_propagator_mpo(propagator_mpo, super_H_site, H_pair, dt, mps_vs_peps)

  END SUBROUTINE setup_propagator




  !!! Compute two-site Hamiltonian !!!
  SUBROUTINE setup_hamiltonian_twosite(hamiltonian_twosite, g, delta, mps_vs_peps)

    USE simulation_parameters, ONLY: local_dim !! Use system parameters  

    COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT)          :: hamiltonian_twosite(:,:)
    REAL(KIND=DP),                 INTENT(IN)           :: g, delta
    CHARACTER(LEN=*),              INTENT(IN), OPTIONAL :: mps_vs_peps

    COMPLEX(KIND=DP) :: super_H_site(local_dim, local_dim)
    COMPLEX(KIND=DP) :: H_pair(local_dim**2, local_dim**2)
    REAL(KIND=DP)    :: lfactor, rfactor

    !! (1) Setup onsite term. 
    super_H_site = 0.0D0
    super_H_site = -g * sigma_x 

    !! (2) Setup two-body term.
    H_pair = 0.0D0
    H_pair = -TensKRON(sigma_z, sigma_z) 

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


END MODULE ising_imaginary_propagator
