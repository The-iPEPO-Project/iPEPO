MODULE pauli_basis

  USE utility

  IMPLICIT NONE

  ! Basis matrices (Pauli matrices)
  COMPLEX (KIND=DP), DIMENSION(4,2,2) :: basis

  ! Matrices used in Hamiltonian
  COMPLEX (KIND=DP), DIMENSION(4,4) :: eye, sigma_x, sigma_y, sigma_z, up, down

  !This module should not be directly included except via the propagator,
  !i.e. this is code to be reused for all 2x2 problems,
 
CONTAINS
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_basic_operators
  ! SYNOPSIS:
  ! Initialise the basic operators
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_basic_operators()

    USE simulation_parameters, ONLY: local_dim, trace_basis

    COMPLEX (KIND=DP), PARAMETER :: &
         & ii=(0.0D0, 1.0D0), nul=(0.0D0, 0.0D0), one=(1.0D0,0.0D0)

    ! Size of basis and size of local_dim
    local_dim=4
    trace_basis=2.0D0
    
    basis = 0.0D0
    ! ID
    basis(1,1,1) = one
    basis(1,2,2) = one
    ! SX
    basis(2,1,2) = one
    basis(2,2,1) = one
    ! SY
    basis(3,1,2) = -ii
    basis(3,2,1) = ii
    ! SZ
    basis(4,1,1) =  one
    basis(4,2,2) = -one
    
    !initialize using rshape to 4x4 case
    !spin operators in 4-by-4 basis. These are for application from the
    !left. Application from the right, i.e. post-multiplication is
    !obtained by transpose (NOT h.c.)

    eye = RESHAPE(&
         & (/  one,  nul,  nul,  nul, &
         &     nul,  one,  nul,  nul, &
         &     nul,  nul,  one,  nul, &
         &     nul,  nul,  nul,  one/),&
         &   (/ 4, 4/))
    sigma_x=RESHAPE(&
         & (/  nul,  one,  nul,  nul, &
         &     one,  nul,  nul,  nul, &
         &     nul,  nul,  nul, +ii , &
         &     nul,  nul, -ii ,  nul/),&
         &   (/ 4, 4/))
    sigma_y=RESHAPE(&
         & (/  nul,  nul,  one,  nul, &
         &     nul,  nul,  nul, -ii , &
         &     one,  nul,  nul,  nul , &
         &     nul, +ii ,  nul,  nul/),&
         &   (/ 4, 4/))
    sigma_z=RESHAPE(&
         & (/  nul,  nul,  nul,  one, &
         &     nul,  nul, +ii ,  nul, &
         &     nul, -ii ,  nul,  nul , &
         &     one,  nul,  nul,  nul/),&
         &   (/ 4, 4/))

    up   = 0.5d0*(sigma_x + ii*sigma_y)
    down = 0.5d0*(sigma_x - ii*sigma_y)
   
  END SUBROUTINE setup_basic_operators

END MODULE pauli_basis
