PROGRAM test_linsol

 USE utility
 USE lapack_wrappers, ONLY: lapack_gen_linsol
 USE definitions_mps_mpo
 USE mps_mpo_utility

 IMPLICIT NONE

  COMPLEX(KIND=DP), ALLOCATABLE :: Amat(:,:), bvec(:), bvecPP(:), Bmat(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: Xvec(:,:)
  INTEGER,          PARAMETER   :: N = 4

  COMPLEX(KIND=DP), ALLOCATABLE :: Rmat(:,:,:,:), Svec(:,:,:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! SPECIFY LINEAR SYSTEM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Setup a quadratic problem specified by Amat, bvec
  ALLOCATE(Amat(N,N), bvec(N), bvecPP(N), Bmat(N,2))
  Amat(1, 1:4) = (/  3,  1,   0,  0 /)
  Amat(2, 1:4) = (/  1,  4,   1,  3 /)
  Amat(3, 1:4) = (/  0,  1,  10,  0 /)
  Amat(4, 1:4) = (/  0,  3,   0,  3 /)

  bvec(1:4)    = (/  1,  1,   1,  1 /)
  bvecPP(1:4)  = (/  1,  3,   2, -1 /)
  Bmat(:,1)    = bvec
  Bmat(:,2)    = bvecPP


  

  !ALLOCATE(Amat(N,N), bvec(N))
  !Amat(1, 1:8) = (/  3,    0,    1,    0,    0,    0,    0,    0 /)
  !Amat(2, 1:8) = (/  0,    6,    0,    2,    0,    0,    0,    0 /)
  !Amat(3, 1:8) = (/  1,    0,    4,    0,    1,    0,    3,    0 /)
  !Amat(4, 1:8) = (/  0,    2,    0,    8,    0,    2,    0,    6 /)
  !Amat(5, 1:8) = (/  0,    0,    1,    0,   10,    0,    0,    0 /)
  !Amat(6, 1:8) = (/  0,    0,    0,    2,    0,   20,    0,    0 /)
  !Amat(7, 1:8) = (/  0,    0,    3,    0,    0,    0,    3,    0 /)
  !Amat(8, 1:8) = (/  0,    0,    0,    6,    0,    0,    0,    6 /)
  !bvec(1:8)    = (/  1,    1,    1,    1,    1,    1,    1,    1 /)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Simple lapack linsol !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !ALLOCATE(Xvec(N,2))
  !CALL lapack_gen_linsol(Xvec,  Amat,  increase_rank(bvecPP, '2'))
  !CALL lapack_gen_linsol(Xvec,  Amat,  Bmat)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Linsol using Rmat, Svec !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Construct Rmat
  Rmat = RESHAPE_4D(TensTRANSPOSE(Amat), '12,34', (/2, 2/), (/2, 2/))

  !! Construct Svec
  ALLOCATE(Svec(2,2,2))
  Svec(1,:,:) = RESHAPE_2D(bvec,    (/2, 2/)) 
  Svec(2,:,:) = RESHAPE_2D(bvecPP,  (/2, 2/)) 

  !! Solve for X
  ALLOCATE(Xvec(N,2))
  CALL lapack_gen_linsol(Xvec,  TensTRANSPOSE(RESHAPE_2D(Rmat, '12,34')),  TensTRANSPOSE(RESHAPE_2D(Svec, '1,23')))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  WRITE(*,*) "LINSOL result: "
  !CALL printMatrix(Xvec)
  CALL printVector(Xvec(:,1))
  WRITE(*,*)
  CALL printVector(Xvec(:,2))

CONTAINS

END PROGRAM test_linsol
