MODULE ising_propagator_pauli_basis

  USE utility
  USE definitions_mps_mpo
  USE propagator_utility
  USE pauli_basis

  IMPLICIT NONE

CONTAINS

  !!! Setup MPO propagator with operators transformed to Pauli basis !!!
  SUBROUTINE setup_propagator(propagator_mpo, g, delta, kappa1, kappa2, dt, xi, mps_vs_peps)
    
    USE simulation_parameters, ONLY: local_dim !! Use system parameters  

    TYPE(block_mpo),  ALLOCATABLE, INTENT(OUT)          :: propagator_mpo(:)
    REAL(KIND=DP),                 INTENT(IN)           :: g, delta, kappa1, kappa2
    COMPLEX(KIND=DP),              INTENT(IN)           :: dt
    REAL(KIND=DP),                 INTENT(IN), OPTIONAL :: xi
    CHARACTER(LEN=*),              INTENT(IN), OPTIONAL :: mps_vs_peps

    COMPLEX(KIND=DP)             :: super_H_site(local_dim, local_dim)
    TYPE(block_mpo), ALLOCATABLE :: H_pair(:)
    INTEGER                      :: ocnt

    !! (1) Setup onsite term. 
    !! NB, all terms in our Hamiltonian would be negative, hence the +ii factor from (-ii)*(-1)
    !! NB, commutator and lindblad are declared in propagator_utility.
    super_H_site = 0.0D0
    super_H_site = + ii * g * commutator(sigma_z) + kappa1 * lindblad(down) + kappa2 * lindblad(up)

    !! (2) Setup two-body term.
    CALL allocate_pair_mpo(H_pair, local_dim, 4)

    ocnt = 0
    CALL mpo_add_commutator(H_pair, ocnt, ii*0.5D0*(1+delta), sigma_x, sigma_x)
    CALL mpo_add_commutator(H_pair, ocnt, ii*0.5D0*(1-delta), sigma_y, sigma_y)

    !! (3) Set up propagators given all sites are identical.
    CALL setup_ti_propagator_mpo(propagator_mpo, super_H_site, H_pair, dt, xi=xi)

  END SUBROUTINE setup_propagator

END MODULE ising_propagator_pauli_basis
