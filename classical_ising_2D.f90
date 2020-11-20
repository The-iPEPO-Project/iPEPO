MODULE classical_ising_2D

 USE utility
 USE definitions_mps_mpo
 USE definitions_peps_pepo
 USE mps_mpo_utility

 USE ctm_definitions
 USE corner_transfer_matrix

 USE project_ipeps,        ONLY: absorb_ipeps_lambdas
 USE TN_contractions

 IMPLICIT NONE

CONTAINS 

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Constructing tensors of 2D classical Ising partition function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! Create MPO representing reduced classical Z network !!!
 SUBROUTINE create_reduced_TN_from_classical_Z(TN2D, beta)

  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT) :: TN2D(:)
  REAL(KIND=DP),                INTENT(IN)    :: beta

  INTEGER :: site

  !! Allocate empty Tmat
  CALL allocate_empty_mpo_block(TN2D, 2)

  !! Construct the sites of Tmat
  DO site=1,2
     CALL allocate_mpo_site(TN2D(site), classical_ising_site(beta))
  end DO

 END SUBROUTINE create_reduced_TN_from_classical_Z





 !!! Create operator sites of reduced classical Z network !!!
 SUBROUTINE create_reduced_op_sites_from_classical_Z(TN_OpSites, beta)

  TYPE(op_sites_type), INTENT(INOUT) :: TN_OpSites
  REAL(KIND=DP),       INTENT(IN)    :: beta

  !! Construct operator sites -- single-site expvals
  TN_OpSites % site1_Op1     = classical_ising_site(beta, 'OP')

  !! Construct operator sites -- two-site correlators
  TN_OpSites % site1_Op1_Op2 = classical_ising_site(beta, 'OP', 'OP')
  TN_OpSites % site1_Op2     = classical_ising_site(beta, 'OP')
  TN_OpSites % site2_Op2     = classical_ising_site(beta, 'OP')

 END SUBROUTINE create_reduced_op_sites_from_classical_Z




 !!! Construct a single site of Ising transfer matrix !!!
 FUNCTION classical_ising_site(beta, spinOpFlag, spinOpFlag1) result(tmat_site)

  USE array_utility

  REAL(KIND=DP),    INTENT(IN)           :: beta
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: spinOpFlag
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: spinOpFlag1

  COMPLEX(KIND=DP) :: tmat_site(2,2,2,2) !! result

  INTEGER :: i,j,k,l !! indices

  !! Compute Tmat site element by element (i,j,k,l = spin-state indices indicating spin-states on adjacent sites)
  DO i=1,2
    DO j=1,2
      DO k=1,2
        DO l=1,2

           IF(PRESENT(spinOpFlag) .AND. PRESENT(spinOpFlag1)) THEN
              !! Compute Tmat site with double spin operator
              tmat_site(i,j,k,l) = classical_ising_site_element(beta,i,j,k,l,spinOpFlag,spinOpFlag1) 
           ELSEIF(PRESENT(spinOpFlag)) THEN
              !! Compute Tmat site with spin operator
              tmat_site(i,j,k,l) = classical_ising_site_element(beta,i,j,k,l,spinOpFlag) 
           ELSE
              !! Compute bare Tmat site
              tmat_site(i,j,k,l) = classical_ising_site_element(beta,i,j,k,l) 
           end IF

        end DO
      end DO
    end DO
  end DO

 END FUNCTION classical_ising_site




 !!! Compute a single matrix element of a given Tmat site (indices correspond to spin states at adjacent sites) !!!
 FUNCTION classical_ising_site_element(beta, i, j, k, l, spinOpFlag, spinOpFlag1) result(tmat_site_element)

  REAL(KIND=DP),    INTENT(IN)           :: beta       !temperature
  INTEGER,          INTENT(IN)           :: i, j, k, l !spin indices
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: spinOpFlag
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: spinOpFlag1

  COMPLEX(KIND=DP) :: tmat_site_element !result = element at indices i,j,k,l

  COMPLEX(KIND=DP) :: Q_sqrt(2,2) !Q_sqrt matrix
  INTEGER          :: s           !local spin on our site

  !! setup Q_sqrt matrix
  CALL setup_Q_sqrt(Q_sqrt, beta)

  !! initialize
  tmat_site_element = (0.0d0, 0.0d0)

  !!! Construct Tmat site with or without operator !!!
  IF(PRESENT(spinOpFlag) .AND. PRESENT(spinOpFlag1)) THEN

     !! Compute Tmat site elements with spin operator -- Sum over all spin values of spinS
     DO s=1,2
        tmat_site_element = tmat_site_element + spinS(s) * spinS(s) * Q_sqrt(i,s) * Q_sqrt(j,s) * Q_sqrt(k,s) * Q_sqrt(l,s)
     end DO

  ELSEIF(PRESENT(spinOpFlag)) THEN

     !! Compute Tmat site elements with spin operator -- Sum over all spin values of spinS
     DO s=1,2
        tmat_site_element = tmat_site_element + spinS(s) * Q_sqrt(i,s) * Q_sqrt(j,s) * Q_sqrt(k,s) * Q_sqrt(l,s)
     end DO

  ELSE

     !! Compute bare Tmat site elements
     DO s=1,2
        tmat_site_element = tmat_site_element + Q_sqrt(i,s) * Q_sqrt(j,s) * Q_sqrt(k,s) * Q_sqrt(l,s)
     end DO
  end IF

 END FUNCTION classical_ising_site_element




 !!! Compute sqrt(Q) bond function !!!
 SUBROUTINE setup_Q_sqrt(Q_sqrt, beta)

  COMPLEX(KIND=DP), INTENT(OUT) :: Q_sqrt(2,2) !output sqrt of Q matrix
  REAL(KIND=DP),    INTENT(IN)  :: beta        !temperature

  COMPLEX(KIND=DP) :: ev1, ev2

  !! sqrt of evals of Q matrix
  ev1 = 0.5d0*sqrt(exp(beta) + exp(-beta))
  ev2 = 0.5d0*sqrt(exp(beta) - exp(-beta))

  !! write to Q_sqrt
  Q_sqrt(:,:) = (0.0d0,0.0d0)
  Q_sqrt(1,1) = ev1 + ev2; Q_sqrt(1,2) = ev1 - ev2
  Q_sqrt(2,1) = ev1 - ev2; Q_sqrt(2,2) = ev1 + ev2

 END SUBROUTINE setup_Q_sqrt



 !!! Find spin given spin-state index s !!!
 FUNCTION spinS(s)

  INTEGER, INTENT(IN) :: s     !input spin-state index
  INTEGER             :: spinS !result = spin-state

  IF(s .EQ. 1) THEN
     spinS = 1
  ELSEIF(s .EQ. 2) THEN
     spinS = -1
  ELSE
     WRITE(*,*) "spinS: invalid value of s: ", s
     STOP
  end IF

 END FUNCTION spinS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE classical_ising_2D
