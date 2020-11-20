MODULE gell_mann_basis

  IMPLICIT NONE
  INTEGER,PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  ! Bare gell mann basis, without operators, re-used for double basis,
  ! and for sun etc.

  ! Basis matrices (for explicit representation
  COMPLEX (KIND=DP), ALLOCATABLE :: basis(:,:,:)
  INTEGER ::  hs_size
    
CONTAINS
  
  ! Helper function for taking trace of given matrix
  FUNCTION mtrace(matrix)
    COMPLEX(KIND=DP) :: matrix(:,:), mtrace
    INTEGER :: i, d
    d=SIZE(matrix,1)

    mtrace=SUM( (/ (matrix(i,i), i=1, d) /))
  end FUNCTION mtrace

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   represent
  ! VARIABLES:
  !           matrix
  ! SYNOPSIS:
  ! Find representation of a given density matrix
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION represent(matrix)

    COMPLEX(KIND=DP), INTENT(IN) :: matrix(hs_size,hs_size)
    COMPLEX(KIND=DP) :: represent(hs_size*hs_size)

    INTEGER :: i

    IF (.NOT. ALLOCATED(basis)) THEN
       WRITE(*,*) "Must setup basis before calling represent"
       STOP
    end IF

    DO i=1, hs_size*hs_size
       represent(i) = 0.5D0*mTRACE(MATMUL(matrix, basis(i,:,:)))
    end DO

  end FUNCTION represent


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   setup_gell_mann_basis
  ! VARIABLES:
  !           p_hs_size - Hilbert space size
  ! SYNOPSIS:
  ! Set up Gell Mann basis states
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_gell_mann_basis(p_hs_size)

    USE simulation_parameters, ONLY: trace_basis, local_dim 

    COMPLEX (KIND=DP), PARAMETER :: ii=(0.0D0, 1.0D0)
    INTEGER, INTENT(IN) ::p_hs_size
    REAL (KIND=DP) :: norm

    INTEGER :: i, j, index


    ! Set module hilbert space size, and local dim
    hs_size=p_hs_size
    local_dim=hs_size*hs_size
    trace_basis=sqrt(2.0D0*hs_size)


    ! Clear and allocate basis
    IF (ALLOCATED(basis)) DEALLOCATE(basis)
    ALLOCATE(basis(local_dim, hs_size, hs_size))
    basis=0.0D0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set up basis states

    ! Identity
    DO i=1, hs_size
       basis(1, i, i) = sqrt(2.0D0/hs_size)
    end DO

    ! Other diagonal matrices
    DO i=1, hs_size-1
       ! i is number of equal elements of diagonal, which may go up to
       ! one less than the matrix size (if all matrix elements equal,
       ! this would be identity, which is included above.)

       ! index is which basis to write, norm is normalisation
       index=i+1
       norm=sqrt(2.0/(i*(i+1)))

       ! Set elements such that orthogonal to ID, and to each other.
       DO j=1, i
          basis(index,j,j) = norm
       end DO
       basis(index,i+1,i+1) = - i*norm
       
    end DO

    ! Set current index, up to hs_size elements set, now need off
    ! diagonal
    index=hs_size

    DO i=1, hs_size
       DO j=i+1, hs_size
          index=index+1
          !Symmetric case:
          basis(index,i,j)=1.0
          basis(index,j,i)=1.0

          index=index+1
          ! Antisymmetric case
          basis(index,i,j)=-ii
          basis(index,j,i)=+ii

       end DO
    end DO
    
  END SUBROUTINE setup_gell_mann_basis


END MODULE gell_mann_basis
