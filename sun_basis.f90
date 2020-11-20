MODULE sun_basis

  USE utility
  USE gell_mann_basis

  IMPLICIT NONE

  LOGICAL :: setup_complete=.FALSE.

  ! Default value;  can be adjusted before call to setup matrix.
  ! Actual basis size used to check if something changed the basis
  ! size since initialization.
  INTEGER :: basis_size=3
  INTEGER, PRIVATE :: actual_basis_size=0

  ! Matrices used in Hamiltonian etc.
  COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: &
       & up_op, dn_op, number_op, g2bare_op, eye, upup_op, dndn_op, dn_rho_up
  
  ! This module should not be directly included except via the propagator.
  ! This code obsolotes all other basis codes (su3, pauli....). Note
  ! that basis_size is adjustable
  
CONTAINS

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_basic_operators
  ! SYNOPSIS:
  ! Initialise the basic operators
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_basic_operators()

    USE simulation_parameters, ONLY:  local_dim

    COMPLEX (KIND=DP), DIMENSION(basis_size,basis_size) :: uptemp, dntemp

    INTEGER :: i, j, index
    REAL (KIND=DP) :: norm

    ! Set up Gell Mann Basis States
    CALL setup_gell_mann_basis(basis_size)

    ! Allocate module level storage, clearing if necessary.
    IF (ALLOCATED(eye)) THEN
       DEALLOCATE(up_op,dn_op,number_op,g2bare_op,eye)
    end IF
    ALLOCATE(  up_op(local_dim, local_dim), &
         &     dn_op(local_dim, local_dim), &
         & number_op(local_dim, local_dim), &
         & g2bare_op(local_dim, local_dim), &
         &       eye(local_dim, local_dim), &
	 &   upup_op(local_dim, local_dim), &
	 &   dndn_op(local_dim, local_dim), &
         &   dn_rho_up(local_dim, local_dim))
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now allocate reprsentations of operators.
    
    ! Identity matrix:  Note that the representation of identity
    ! is identity.  i.e. if one calculates:
    !     0.5*Trace(b(i)*one*b(j))
    ! one gets 0.5*2*delta_{ij}; it is the latter we are using here.
    eye=0.0D0
    DO i=1, local_dim
       eye(i,i) = 1.0D0
    end DO
    
    ! Initialize basic operators, by taking trace of bare operator.
    ! Note that we have choosen to have increasing occupation with
    ! increasing index in the matrix rep space, i.e. the bottom right
    ! component is the highest occupation.  This only matters when
    ! we come to identify meaning of states. This way round is natural
    ! for SU(N) but for some reason SU(2) normally other way...

    uptemp=0.0D0;
    DO i=2, basis_size
       uptemp(i, i-1)=sqrt(1.0D0*(i-1))
    end DO
    

    DO i=1, local_dim
       DO j=1, local_dim
          up_op(i,j) = 0.5D0*&
               & MTRACE(MATMUL(basis(i,:,:), MATMUL(uptemp, basis(j,:,:))))
       end DO
    end DO

    ! Setup dn_rho_up: a compound of dn_op to Left-mult rho and up_op to Right-mult rho
    dntemp = CONJG(TRANSPOSE(uptemp))
    DO i=1, local_dim
       DO j=1, local_dim
          dn_rho_up(i,j) = 0.5D0*&
               & MTRACE(MATMUL(MATMUL(basis(i,:,:), MATMUL(dntemp, basis(j,:,:))), uptemp))
       end DO
    end DO

    ! Derive other operators from this.  Hermitian conjugate gives
    ! other operator (Transpose alone moves it to other side).
    dn_op=CONJG(TRANSPOSE(up_op))
    number_op=MATMUL(up_op,dn_op)

    ! dndn_op is two photon loss operator, upup_op two photon gain
    dndn_op=MATMUL(dn_op,dn_op)
    upup_op=MATMUL(up_op,up_op)

    ! For G2 and interactions:
    g2bare_op = MATMUL(up_op, MATMUL(number_op, dn_op))

    ! Record that we have finished
    actual_basis_size=basis_size
    setup_complete=.TRUE.

  end SUBROUTINE setup_basic_operators
  
  
END MODULE sun_basis
