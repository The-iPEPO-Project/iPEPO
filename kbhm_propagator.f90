MODULE kbhm_propagator

  USE utility
  USE definitions_mps_mpo
  USE sun_basis
  USE propagator_utility

  IMPLICIT NONE

  !! Bose Hubbard model with coherent pumping with a given phase pattern. 
  !! This is implemented by a phase factor in the hopping between adjacent sites.

CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_propagator
  ! VARIABLES:
  !           delta - Detuning of pump from cavity, i.e omega_p - omega_c,
  !                   so that in rotating frame, get - Delta in Heff
  !           f     - Strength of drive
  !           U     - Boson nonlinearity
  !           J     - Complex hopping amplitude
  !           kappa - Photon loss
  !           dt
  ! SYNOPSIS:
  ! Propagator for coherent driven BHM, with a phase twist between sites.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_propagator(propagator_mpo, delta, f, U, J, kappa, dt, xi)
    
    USE simulation_parameters, ONLY: local_dim !Use system parameters  

    TYPE(block_mpo),  ALLOCATABLE, INTENT(OUT)          :: propagator_mpo(:)
    REAL(KIND=DP),                 INTENT(IN)           :: delta, f, U, kappa
    COMPLEX(KIND=DP),              INTENT(IN)           :: J, dt
    REAL(KIND=DP),                 INTENT(IN), OPTIONAL :: xi

    COMPLEX(KIND=DP)             :: super_H_site(local_dim, local_dim)
    TYPE(block_mpo), ALLOCATABLE :: H_pair(:)
    INTEGER :: ocnt

    !! (1) Setup onsite term. 
    !! NB, all terms in our Hamiltonian would be negative, hence the +ii factor from (-ii)*(-1)
    !! NB, commutator and lindblad are declared in propagator_utility.
    super_H_site = 0.0D0
    super_H_site = -ii*commutator(delta*number_op + 0.5*U*g2bare_op + 0.5*f*(up_op + dn_op)) + 0.5*kappa*lindblad(dn_op)

    !! (2) Setup two-body term.
    CALL allocate_pair_mpo(H_pair, local_dim, 4)

    ocnt = 0
    CALL mpo_add_commutator(H_pair, ocnt, ii      *J,  up_op, dn_op)
    CALL mpo_add_commutator(H_pair, ocnt, ii*CONJG(J), dn_op, up_op)

    !! (3) Set up propagators given all sites are identical.
    CALL setup_ti_propagator_mpo(propagator_mpo, super_H_site, H_pair, dt, xi=xi)

  END SUBROUTINE setup_propagator

END MODULE kbhm_propagator
