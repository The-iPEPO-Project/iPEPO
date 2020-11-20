MODULE xyz_propagator_pauli_basis

  USE utility
  USE definitions_mps_mpo
  USE propagator_utility
  USE pauli_basis

  IMPLICIT NONE

CONTAINS

  !!! Setup MPO propagator with operators transformed to Pauli basis  !!! 
  SUBROUTINE setup_propagator(propagator_mpo, Jx, Jy, Jz, kappa, dt, mps_vs_peps)
    
    USE simulation_parameters, ONLY: local_dim !! Use system parameters  

    TYPE(block_mpo),  ALLOCATABLE, INTENT(OUT)          :: propagator_mpo(:)
    REAL(KIND=DP),                 INTENT(IN)           :: Jx, Jy, Jz, kappa
    COMPLEX(KIND=DP),              INTENT(IN)           :: dt
    CHARACTER(LEN=*),              INTENT(IN), OPTIONAL :: mps_vs_peps

    COMPLEX(KIND=DP)             :: super_H_site(local_dim, local_dim)
    TYPE(block_mpo), ALLOCATABLE :: H_pair(:)
    INTEGER                      :: ocnt

    !! (1) Setup onsite term. 
    !! NB, all terms in our Hamiltonian would be negative, hence the +ii factor from (-ii)*(-1)
    !! NB, commutator and lindblad are declared in propagator_utility.
    super_H_site = 0.0D0
    super_H_site = 0.5D0 * kappa * lindblad(down)

    !! (2) Setup two-body term.
    CALL allocate_pair_mpo(H_pair, local_dim, 6)

    ocnt = 0
    CALL mpo_add_commutator(H_pair, ocnt, ii*Jx, sigma_x, sigma_x)
    CALL mpo_add_commutator(H_pair, ocnt, ii*Jy, sigma_y, sigma_y)
    CALL mpo_add_commutator(H_pair, ocnt, ii*Jz, sigma_z, sigma_z)

    !! (3) Set up propagators given all sites are identical.
    CALL setup_ti_propagator_mpo(propagator_mpo, super_H_site, H_pair, dt, mps_vs_peps)

  END SUBROUTINE setup_propagator


END MODULE xyz_propagator_pauli_basis
