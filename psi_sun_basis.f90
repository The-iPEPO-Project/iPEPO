MODULE psi_sun_basis

  USE utility

  IMPLICIT NONE

  LOGICAL :: setup_complete=.FALSE.

  !! Default value;  can be adjusted before call to setup matrix.
  !! Actual basis size used to check if something changed the basis
  !! size since initialization.
  INTEGER          :: basis_size=3
  INTEGER, PRIVATE :: actual_basis_size=0

  !! Matrices used in Hamiltonian etc.
  COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: &
       & up_op, dn_op, number_op, g2bare_op, eye, upup_op, dndn_op
  
  !! This module should not be directly included except via the propagator.
  !! This code obsolotes all other basis codes (su3, pauli....). Note
  !! that basis_size is adjustable
  
CONTAINS

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_basic_operators
  ! SYNOPSIS:
  ! Initialise the basic operators
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_basic_operators(rhoFlag)

    USE simulation_parameters, ONLY: local_dim, hs_dim

    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: rhoFlag

    INTEGER :: i, j
    REAL (KIND=DP) :: norm

    !! Set local_dim
    local_dim = hs_dim

    IF(PRESENT(rhoFlag)) local_dim = local_dim**2

    !! Allocate module level storage, clearing if necessary.
    IF (ALLOCATED(eye)) THEN
       DEALLOCATE(up_op, dn_op, number_op, g2bare_op, eye)
    end IF
    ALLOCATE(  up_op(hs_dim, hs_dim), &
         &     dn_op(hs_dim, hs_dim), &
         & number_op(hs_dim, hs_dim), &
         & g2bare_op(hs_dim, hs_dim), &
         &       eye(hs_dim, hs_dim), &
	 &   upup_op(hs_dim, hs_dim), &
	 &   dndn_op(hs_dim, hs_dim))
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Now allocate representations of operators.
    
    !! Identity matrix
    eye=0.0D0
    DO i=1, hs_dim
       eye(i,i) = 1.0D0
    end DO
    
    !! Initialize basic operators.
    !! Note that we have choosen to have increasing occupation with
    !! increasing index in the matrix rep space, i.e. the bottom right
    !! component is the highest occupation.  This only matters when
    !! we come to identify meaning of states. This way round is natural
    !! for SU(N) but for some reason SU(2) normally other way...
    up_op=0.0D0
    DO i=2,hs_dim
       up_op(i,i-1) = sqrt(1.0D0*(i-1))
    end DO
    

    !! Derive other operators from this.  Hermitian conjugate gives
    !! other operator (Transpose alone moves it to other side).
    dn_op=CONJG(TRANSPOSE(up_op))
    number_op=MATMUL(up_op,dn_op)

    !! dndn_op is two photon loss operator, upup_op two photon gain
    dndn_op=MATMUL(dn_op,dn_op)
    upup_op=MATMUL(up_op,up_op)

    !! For G2 and interactions:
    g2bare_op = MATMUL(up_op, MATMUL(number_op, dn_op))

    !! Record that we have finished
    basis_size        = hs_dim
    actual_basis_size = hs_dim
    setup_complete    = .TRUE.

  END SUBROUTINE setup_basic_operators
  
  
END MODULE psi_sun_basis
