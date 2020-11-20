MODULE xyz_propagator

  USE utility
  USE definitions_mps_mpo
  USE propagator_utility
  USE psi_pauli_basis

  IMPLICIT NONE

  interface setup_propagator
    module procedure setup_propagator_mpo
    module procedure setup_propagator_dense
  end interface setup_propagator

  private setup_propagator_mpo, setup_propagator_dense

CONTAINS

  !!! Compute two-site Liouvillian !!!
  SUBROUTINE setup_liouvillian_twosite(liouvillian_twosite, Jx, Jy, Jz, kappa, mps_vs_peps)

    USE simulation_parameters, ONLY: local_dim !! Use system parameters  

    COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT)          :: liouvillian_twosite(:,:)
    REAL(KIND=DP),                 INTENT(IN)           :: Jx, Jy, Jz, kappa
    CHARACTER(LEN=*),              INTENT(IN), OPTIONAL :: mps_vs_peps

    COMPLEX(KIND=DP) :: super_H_site(local_dim, local_dim)
    COMPLEX(KIND=DP) :: H_pair(local_dim**2, local_dim**2)
    COMPLEX(KIND=DP) :: super_eye(local_dim, local_dim)
    REAL(KIND=DP)    :: lfactor, rfactor

    !! (1) Setup onsite term. 
    !!     NB, all terms in our Hamiltonian would be negative, hence the -ii factor from (-ii)*(1)
    !!     NB, commutator and lindblad are declared in propagator_utility.
    super_H_site = 0.0D0 
    super_H_site = 0.5D0 * kappa * sv_lindblad(down)  

    !! (2) Setup two-body term.
    H_pair = 0.0D0
    H_pair =   - ii * Jx * sv_cross_commutator(sigma_x, sigma_x) &
             & - ii * Jy * sv_cross_commutator(sigma_y, sigma_y) & 
             & - ii * Jz * sv_cross_commutator(sigma_z, sigma_z)

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
    ALLOCATE(liouvillian_twosite(local_dim**2, local_dim**2))
    super_eye = matEye(SIZE(super_H_site,1))
    liouvillian_twosite = lfactor * TensKRON(super_H_site, super_eye) + rfactor * TensKRON(super_eye, super_H_site) + H_pair

  END SUBROUTINE setup_liouvillian_twosite



  !!! Setup dense propagator !!!
  SUBROUTINE setup_propagator_dense(propagator_dense, Jx, Jy, Jz, kappa, dt, mps_vs_peps)
    
    USE simulation_parameters, ONLY: local_dim !! Use system parameters  

    COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT)          :: propagator_dense(:,:)
    REAL(KIND=DP),                 INTENT(IN)           :: Jx, Jy, Jz, kappa
    COMPLEX(KIND=DP),              INTENT(IN)           :: dt
    CHARACTER(LEN=*),              INTENT(IN), OPTIONAL :: mps_vs_peps

    COMPLEX(KIND=DP) :: super_H_site(local_dim, local_dim)
    COMPLEX(KIND=DP) :: H_pair(local_dim**2, local_dim**2)

    !! (1) Setup onsite term. 
    !!     NB, all terms in our Hamiltonian would be negative, hence the -ii factor from (-ii)*(1)
    !!     NB, commutator and lindblad are declared in propagator_utility.
    super_H_site = 0.0D0 
    super_H_site = 0.5D0 * kappa * sv_lindblad(down)

    !! (2) Setup two-body term.
    H_pair = 0.0D0
    H_pair =   - ii * Jx * sv_cross_commutator(sigma_x, sigma_x) &
             & - ii * Jy * sv_cross_commutator(sigma_y, sigma_y) & 
             & - ii * Jz * sv_cross_commutator(sigma_z, sigma_z) 

    !! (3) Set up propagators given all sites are identical.
    CALL setup_ti_propagator_dense(propagator_dense, super_H_site, super_H_site, H_pair, dt, mps_vs_peps)
    
  END SUBROUTINE setup_propagator_dense



  !!! Setup MPO propagator !!!
  SUBROUTINE setup_propagator_mpo(propagator_mpo, Jx, Jy, Jz, kappa, dt, mps_vs_peps, xi)

    USE simulation_parameters, ONLY: local_dim !! Use system parameters  

    TYPE(block_mpo),  ALLOCATABLE, INTENT(INOUT)        :: propagator_mpo(:)
    REAL(KIND=DP),                 INTENT(IN)           :: Jx, Jy, Jz, kappa
    COMPLEX(KIND=DP),              INTENT(IN)           :: dt
    CHARACTER(LEN=*),              INTENT(IN), OPTIONAL :: mps_vs_peps
    REAL(KIND=DP),                 INTENT(IN), OPTIONAL :: xi

    COMPLEX(KIND=DP), ALLOCATABLE :: prop_2D(:,:), prop_4D(:,:,:,:)
    COMPLEX(KIND=DP)              :: super_H_site(local_dim, local_dim)
    COMPLEX(KIND=DP)              :: H_pair(local_dim**2, local_dim**2)

    !! (1) Setup onsite term. 
    !!     NB, all terms in our Hamiltonian would be negative, hence the -ii factor from (-ii)*(1)
    !!     NB, commutator and lindblad are declared in propagator_utility.
    super_H_site = 0.0D0 
    super_H_site = 0.5D0 * kappa * sv_lindblad(down)

    !! (2) Setup two-body term.
    H_pair = 0.0D0
    H_pair =   - ii * Jx * sv_cross_commutator(sigma_x, sigma_x) &
             & - ii * Jy * sv_cross_commutator(sigma_y, sigma_y) & 
             & - ii * Jz * sv_cross_commutator(sigma_z, sigma_z) 

    !! (3) Set up dense propagator given all sites are identical.
    CALL setup_ti_propagator_dense(prop_2D, super_H_site, super_H_site, H_pair, dt, mps_vs_peps)

    !! (4) Reshape into 4D
    prop_4D = RESHAPE_4D(prop_2D, '12,34', (/ local_dim, local_dim /), (/ local_dim, local_dim /))

    !! (5) Construct MPO propagator
    CALL canonicalize_propagator(propagator_mpo, prop_4D, xi)
    
  END SUBROUTINE setup_propagator_mpo

END MODULE xyz_propagator
