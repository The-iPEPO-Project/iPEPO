MODULE psi_pauli_basis

 USE utility

 IMPLICIT NONE

  !spin operators 
  COMPLEX(KIND=DP), DIMENSION(2,2) :: eye, sigma_x, sigma_y, sigma_z, up, down 
  
CONTAINS

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_basic_operators
  ! SYNOPSIS:
  ! Initialise the basic operators
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_basic_operators(rhoFlag)

  USE simulation_parameters, ONLY: local_dim, hs_dim

  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: rhoFlag

  COMPLEX(KIND=DP), PARAMETER :: ii=(0.0D0, 1.0D0), nul=(0.0D0, 0.0D0), one=(1.0D0,0.0D0)

  !Set local dim
  local_dim = 2
  hs_dim    = local_dim

  IF(PRESENT(rhoFlag)) THEN
     hs_dim    = local_dim
     local_dim = local_dim**2
  end IF 

  !Initialize
  eye=0.0D0; sigma_x=0.0D0; sigma_y=0.0D0; sigma_z=0.0D0
    
  !ID
  eye(1,1) = one
  eye(2,2) = one

  !SX
  sigma_x(1,2) = one
  sigma_x(2,1) = one

  !SY
  sigma_y(1,2) = -ii
  sigma_y(2,1) =  ii

  !SZ
  sigma_z(1,1) =  one
  sigma_z(2,2) = -one

  up   = 0.5D0*(sigma_x + ii*sigma_y)
  down = 0.5D0*(sigma_x - ii*sigma_y)

 END SUBROUTINE setup_basic_operators





 !!! Setup two-site vectorized basis !!!
 SUBROUTINE setup_basis_twosite(basis_vecs)

   USE simulation_parameters, ONLY: local_dim, hs_dim  !! Use system parameters  

   COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: basis_vecs(:,:)

   COMPLEX(KIND=DP), ALLOCATABLE :: basis(:,:,:), tmp(:)
   INTEGER                       :: i,j

   !! Put pauli basis matrices in one var
   ALLOCATE(basis(local_dim, hs_dim, hs_dim))
   basis(1,:,:) = eye
   basis(2,:,:) = sigma_x
   basis(3,:,:) = sigma_y
   basis(4,:,:) = sigma_z

   !! Construct array of basis vecs -- to access each vector use: basis_vecs(k,:)
   ALLOCATE(basis_vecs(local_dim**2, local_dim**2))
   DO i=1,local_dim
     DO j=1,local_dim
        basis_vecs(ICOM(i,j,local_dim), :) = RESHAPE_1D(TensKRON(basis(i,:,:), basis(j,:,:)))
     end DO
   end DO
    
  END SUBROUTINE setup_basis_twosite





  !!! Setup one-site vectorized basis !!!
  SUBROUTINE setup_basis_onesite(basis_vecs)

   USE simulation_parameters, ONLY: local_dim, hs_dim  !! Use system parameters  

   COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: basis_vecs(:,:)

   COMPLEX(KIND=DP), ALLOCATABLE :: basis(:,:,:), tmp(:)
   INTEGER                       :: i

   !! Put pauli basis matrices in one var
   ALLOCATE(basis(local_dim, hs_dim, hs_dim))
   basis(1,:,:) = eye
   basis(2,:,:) = sigma_x
   basis(3,:,:) = sigma_y
   basis(4,:,:) = sigma_z

   !! Construct array of basis vecs -- to access each vector use: basis_vecs(k,:)
   ALLOCATE(basis_vecs(local_dim, local_dim))
   DO i=1,local_dim
      basis_vecs(i,:) = SQRT(0.5D0)*RESHAPE_1D(basis(i,:,:))
   end DO
    
  END SUBROUTINE setup_basis_onesite


END MODULE psi_pauli_basis
