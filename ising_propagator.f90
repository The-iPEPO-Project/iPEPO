MODULE ising_propagator

  USE utility
  USE definitions_mps_mpo
  USE propagator_utility
  USE psi_pauli_basis

  IMPLICIT NONE

  interface setup_propagator
    module procedure setup_propagator_mpo_v1
    module procedure setup_propagator_dense_v1
  end interface setup_propagator

  private setup_propagator_mpo_v1, setup_propagator_dense_v1

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!! Setup Ising propagator (model-1 with g, delta) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!! Setup MPO propagator !!!
  SUBROUTINE setup_propagator_mpo_v1(propagator_mpo, g, kappa1, kappa2, dt, mps_vs_peps, xi)
    
    USE simulation_parameters, ONLY: local_dim !! Use system parameters  

    TYPE(block_mpo),  ALLOCATABLE, INTENT(INOUT)        :: propagator_mpo(:)
    REAL(KIND=DP),                 INTENT(IN)           :: g, kappa1, kappa2
    COMPLEX(KIND=DP),              INTENT(IN)           :: dt
    CHARACTER(LEN=*),              INTENT(IN)           :: mps_vs_peps
    REAL(KIND=DP),                 INTENT(IN), OPTIONAL :: xi

    COMPLEX(KIND=DP), ALLOCATABLE :: prop_2D(:,:), prop_4D(:,:,:,:)
    COMPLEX(KIND=DP)              :: super_H_site(local_dim, local_dim)
    COMPLEX(KIND=DP)              :: H_pair(local_dim**2, local_dim**2)

    !! (1) Setup onsite term. 
    !!     NB, all terms in our Hamiltonian would be negative, hence the -ii factor from (-ii)*(1)
    !!     NB, commutator and lindblad are declared in propagator_utility.
    super_H_site = 0.0D0 
    super_H_site = ii * g * sv_commutator(sigma_z) + 0.5D0 * kappa1 * sv_lindblad(down) 

    !! (2) Setup two-body term.
    H_pair = 0.0D0
    H_pair = ii * sv_cross_commutator(sigma_x, sigma_x)

    !! (3) Set up dense propagator given all sites are identical.
    CALL setup_ti_propagator_dense(prop_2D, super_H_site, super_H_site, H_pair, dt, mps_vs_peps)

    !! (4) Reshape into 4D
    prop_4D = RESHAPE_4D(prop_2D, '12,34', (/ local_dim, local_dim /), (/ local_dim, local_dim /))

    !! (5) Construct MPO propagator
    CALL canonicalize_propagator(propagator_mpo, prop_4D, xi)
    
  END SUBROUTINE setup_propagator_mpo_v1



  !!! Setup dense propagator !!!
  SUBROUTINE setup_propagator_dense_v1(propagator_dense, g, kappa1, kappa2, dt, mps_vs_peps)
    
    USE simulation_parameters, ONLY: local_dim !! Use system parameters  

    COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT) :: propagator_dense(:,:)
    REAL(KIND=DP),                 INTENT(IN)  :: g, kappa1, kappa2
    COMPLEX(KIND=DP),              INTENT(IN)  :: dt
    CHARACTER(LEN=*),              INTENT(IN)  :: mps_vs_peps

    COMPLEX(KIND=DP) :: super_H_site(local_dim, local_dim)
    COMPLEX(KIND=DP) :: H_pair(local_dim**2, local_dim**2)

    !! (1) Setup onsite term. 
    !!     NB, all terms in our Hamiltonian would be negative, hence the ii factor from (-ii)*(-1)
    !!     NB, commutator and lindblad are declared in propagator_utility.
    super_H_site =  0.0D0
    super_H_site = ii * g * sv_commutator(sigma_z) + 0.5D0 * kappa1 * sv_lindblad(down) 

    !! (2) Setup two-body term.
    H_pair =  0.0D0
    H_pair = ii * sv_cross_commutator(sigma_x, sigma_x)

    !! (3) Set up propagators given all sites are identical.
    CALL setup_ti_propagator_dense(propagator_dense, super_H_site, super_H_site, H_pair, dt, mps_vs_peps)

  END SUBROUTINE setup_propagator_dense_v1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE ising_propagator
