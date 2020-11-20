MODULE arnoldi_wrappers

  IMPLICIT NONE

  INTEGER, PRIVATE, PARAMETER :: DP=KIND(1.0D0)

  ! This contains routines associated with ARPACK arnoldi evaluatlion.
  ! Routines for calculating application of MPO to object are in
  ! the utility funciton.  There is an equivalent version of this 
  ! using NAG

  ! 0=off, 1=minimum, 2=full
  INTEGER :: DBG_ARNOLDI=0 

CONTAINS

 !!! Find the dominant eval && evec using Arnoldi !!!
 SUBROUTINE mpo_arnoldi_largest_eig(app_op_to_vec, eval, evec, which_evec)

  USE arnoldi_utility, ONLY: check_dim_match, macheps, tol_factor

  !! Apply-Op-To-Vec routine
  INTERFACE
    SUBROUTINE app_op_to_vec(vi, vo, vecFlag)
         COMPLEX(KIND(1.0D0)), INTENT(IN)  :: vi(:)
         COMPLEX(KIND(1.0D0)), INTENT(OUT) :: vo(:)
         CHARACTER(LEN=*),     INTENT(IN)  :: vecFlag
    end SUBROUTINE app_op_to_vec
  end INTERFACE

  !! Dominant Eval && Evec to be computed
  COMPLEX(KIND=DP),              INTENT(INOUT) :: eval
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: evec(:)
  CHARACTER(LEN=*),              INTENT(IN)    :: which_evec


  !! Number of eigenvalues we are trying to find.
  INTEGER :: nev

  !! Output evec for testing the residuals
  COMPLEX(KIND=DP), ALLOCATABLE :: evec_out(:)
  
  !! evec norm const && location of the min phase element
  COMPLEX(KIND=DP) :: C_norm
  INTEGER          :: minPhiLoc(1)

  !! Dimension of problem (of matrix), max number singular vals
  INTEGER :: ncv, lworkl, nconv, ev, ev_sel

  !! If reseting nev, record whether it unchanged.  If we had to shrink nev this time, 
  !! we do not let it grow again, on this step, as otherwise we get into tight loops.
  LOGICAL   :: nev_unchanged, cease_growing_nev
  INTEGER   :: irevcm, info, iparam(11), ipntr(14)
  INTEGER   :: Ndim
  CHARACTER :: bmat*1, which*2

  COMPLEX(KIND=DP), ALLOCATABLE :: resid(:), v(:,:), workd(:), ax(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: evals(:), workl(:), workev(:)
  LOGICAL,          ALLOCATABLE :: selection(:)
  REAL(KIND=DP),    ALLOCATABLE :: rwork(:)
  REAL(KIND=DP)                 :: tol, shift, norm

  REAL(KIND=DP) :: timers(3), otime, ctime

  EXTERNAL :: znaupd

  !! Find size of evec (Our matrix = [Ndim x Ndim], vec = [Ndim])
  nev=1     

  !! Size of evecs
  Ndim=SIZE(evec)

  !! Number of Arnoldi vectors
  ncv=MIN(3*nev+1,Ndim)

  !! Work array size
  lworkl=4*ncv*(ncv+5)
  !! Tolerance in solver
  tol=1.0D-14 !1.0D-10
  !! Shift, unused as simple problem
  shift=0.0D0

  !! In case we had to re-run after changing nev, need to free this up
  IF(ALLOCATED(resid)) DEALLOCATE(Resid, v, ax, evals, workd, workl, rwork, selection, workev)

  !! Work spaces etc. ax, used in multiplication.
  ALLOCATE(resid(Ndim), v(Ndim,ncv), ax(Ndim), evals(nev+1), evec_out(Ndim))
  ALLOCATE(workd(3*Ndim), workl(lworkl), rwork(ncv))
  ALLOCATE(selection(ncv), workev(2*ncv))

  !! Initial vec
  IF(ALLOCATED(evec)) resid = evec

  !Reverse communication flags etc
  irevcm=0
  info=1 !0        ! Reset residuals (CHANGE THIS BACK IF YOU WANNA RANDOM INITIAL VECS)

  iparam(1)=1      ! Use exact shifts
  iparam(3)=1000*nev ! Max iterations
  iparam(7)=1      ! Mode 1, i.e. simple, not generalized.
  bmat='I'         ! Also required for mode1.
  which='LM'       ! Want largest magnitude vals.
  selection=.TRUE.  ! Unreferenced as we ask for all.

  CALL cpu_time(otime); timers=0
  IF(DBG_ARNOLDI .GE. 1) WRITE(*,'(" ARNOLDI(ARPACK): nev=",I3," Ndim=",I5)') nev, Ndim

  !! Run Arnoldi Loop with ZNAUPD && App-Op-To-Vec until convergence
  arnoldi_loop: DO

       !! Call ZNAUPD to build Arnoldi basis
       CALL znaupd(irevcm, &
               & bmat, Ndim, which, nev, tol, resid, ncv, v, Ndim,  &
               & iparam, ipntr, workd, workl, lworkl, rwork, info)

       CALL largest_eig_timer(1)

       !! Check if no errors in ZNAUPD
       IF(info .NE. 0) THEN
          WRITE(*,*) "Error code in znaupd. Code:", info, " irevcm ", irevcm
          EXIT arnoldi_loop
       end IF

       !! Apply op to vec, or exit, or throw error
       IF(irevcm .EQ. 99) THEN

          !! Once irevcm = 99, it means convergence has been achieved -- Exit Arnoldi Loop
          EXIT arnoldi_loop

       ELSEIF(irevcm .EQ. -1 .OR. irevcm .EQ. 1) THEN

          !! Apply Linear operator to vector
          CALL app_op_to_vec(workd(ipntr(1):ipntr(1)+Ndim-1), workd(ipntr(2):ipntr(2)+Ndim-1), which_evec)
          CALL largest_eig_timer(2)

       ELSE 
          WRITE(*,*) "Unexpected request in Reverse communication"
       end IF

  end DO arnoldi_loop

  !WRITE(*,*) "num of znaupd iters: ", iparam(3)
       
  !! Use ZNEUPD for post processing to find evals && evecs (matrix V). Note that we reuse V. 
  !! NB. 'P' is Schur, we use this because Ritz does not guarantee orthogonality of degenerate vectors 
  !!     -- we use Schur as it is orthogonal.
  CALL ZNEUPD(.TRUE., 'P', selection, evals, v, Ndim, shift, workev, &
            & bmat, Ndim, which, nev, tol, resid, ncv, v, Ndim, &
            & iparam, ipntr, workd, workl, lworkl, rwork, info)

  !! Check if no errors in ZNEUPD
  IF(info .NE. 0) THEN
     WRITE(*,*) "Error code in zneupd. Code: ", info
     STOP
  end IF

  !! Number of converged eigenvalues
  nconv=iparam(5)

  !! Zero everything
  eval=0.0d0; evec = 0.0d0

  !! Extract eval && evec
  eval = evals(1)
  evec = v(:,1)

  CALL largest_eig_timer(3)

  IF(DBG_ARNOLDI .GE. 2) &
         & WRITE(*,'(" ARNOLDI(ARPACK): Loop Timings: ",4(F6.2,X))')  timers

  DEALLOCATE(resid, v, ax, workd, workl, rwork, selection, workev)

  !! Find the element with smallest complex phase, rescale evec by that element (to fix the gauge)
  minPhiLoc = MINLOC(REAL(ATAN(evec)), ABS(REAL(evec)) .GT. tol)
  evec = evec/evec(minPhiLoc(1))

  !! Normalize evec:
  C_norm = SUM(CONJG(evec) * evec)
  evec = evec/sqrt(C_norm)

 CONTAINS

    SUBROUTINE largest_eig_timer(n)
      INTEGER :: n
      CALL cpu_time(ctime)
      timers(n) = timers(n) + ctime-otime
      otime=ctime
    end SUBROUTINE largest_eig_timer
    
 END SUBROUTINE mpo_arnoldi_largest_eig





 !Find the dominant eval && evec using Arnoldi
 SUBROUTINE OLD_mps_arnoldi_largest_eig(app_op_to_vec, eval, evec, which_evec)

  USE arnoldi_utility, ONLY: check_dim_match, update_nev, macheps, tol_factor

  !Apply-Op-To-Vec routine
  INTERFACE
    SUBROUTINE app_op_to_vec(vi, vo, vecFlag)
         COMPLEX(KIND(1.0D0)), INTENT(IN) :: vi(:)
         COMPLEX(KIND(1.0D0)), INTENT(OUT) :: vo(:)
         CHARACTER(LEN=*), INTENT(IN) :: vecFlag
    end SUBROUTINE app_op_to_vec
  end INTERFACE

  !Dominant Eval && Evec to be computed
  COMPLEX(KIND=DP), INTENT(INOUT) :: eval
  COMPLEX(KIND=DP), ALLOCATABLE, INTENT(INOUT) :: evec(:)
  CHARACTER(LEN=*), INTENT(IN) :: which_evec

  !Number of eigenvalues we are trying to find.
  INTEGER :: nev

  !Output evec for testing the residuals
  COMPLEX(KIND=DP), ALLOCATABLE :: evec_out(:)
  
  !evec norm const && location of the min phase element
  COMPLEX(KIND=DP) :: C_norm
  INTEGER :: minPhiLoc(1)

  !Dimension of problem (of matrix), max number singular vals
  INTEGER :: max_nev, ncv, lworkl, nconv, ev, ev_sel
  !If reseting nev, record whether it unchanged.  If we had to shrink nev this time, 
  !we do not let it grow again, on this step, as otherwise we get into tight loops.
  LOGICAL :: nev_unchanged, cease_growing_nev
  INTEGER :: irevcm, info, iparam(11), ipntr(14)
  INTEGER :: Ndim
  CHARACTER :: bmat*1, which*2

  COMPLEX(KIND=DP), ALLOCATABLE :: resid(:), v(:,:), workd(:), ax(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: evals(:), workl(:), workev(:)
  LOGICAL, ALLOCATABLE :: selection(:)
  REAL(KIND=DP), ALLOCATABLE :: rwork(:)
  REAL(KIND=DP) :: tol, shift, norm

  REAL(KIND=DP) :: timers(3), otime, ctime

  EXTERNAL :: znaupd

  !Find size of evec (Our matrix = [Ndim x Ndim], vec = [Ndim])
  max_nev=1
  nev=max_nev    
  cease_growing_nev=.FALSE.

  nevloop: DO 

       !Size of evecs
       Ndim=SIZE(evec)

       !Number of Arnoldi vectors
       ncv=MIN(3*nev+1,Ndim)

       !Work array size
       lworkl=4*ncv*(ncv+5)
       !Tolerance in solver
       tol=1.0D-12 !1.0D-14 !1.0D-08
       !Shift, unused as simple problem
       shift=0.0D0

       !In case we had to re-run after changing nev, need to free this up
       IF(ALLOCATED(resid)) DEALLOCATE(Resid, v, ax, evals, workd, workl, rwork, selection, workev)

       !Work spaces etc. ax, used in multiplication.
       ALLOCATE(resid(Ndim), v(Ndim,ncv), ax(Ndim), evals(nev+1), evec_out(Ndim))
       ALLOCATE(workd(3*Ndim), workl(lworkl), rwork(ncv))
       ALLOCATE(selection(ncv), workev(2*ncv))

       !Initial vec FIXME
       IF(ALLOCATED(evec)) resid = evec

       ! Reverse communication flags etc
       irevcm=0
       info=1 !0          ! Reset residuals. FIXME FIXME FIXME FIXME FIXME CHANGE THIS BACK IF YOU WANNA RANDOM INITIAL VECS FIXME FIXME FIXME FIXME

       iparam(1)=1        ! Use exact shifts
       iparam(3)=5000*nev ! Max iterations
       iparam(7)=1        ! Mode 1, i.e. simple, not generalized.
       bmat='I'           ! Also required for mode1.
       which='LM'         ! Want largest magnitude vals.
       selection=.TRUE.   ! Unreferenced as we ask for all.

       CALL cpu_time(otime); timers=0
       IF (DBG_ARNOLDI .GE. 1) &
            & WRITE(*,'(" ARNOLDI(ARPACK): nev=",I3," Ndim=",I5)') nev, Ndim

       !Run Arnoldi Loop with ZNAUPD && App-Op-To-Vec until convergence
       arnoldi_loop: DO

          !Call ZNAUPD to build Arnoldi basis
          CALL znaupd(irevcm, &
               & bmat, Ndim, which, nev, tol, resid, ncv, v, Ndim,  &
               & iparam, ipntr, workd, workl, lworkl, rwork, info)

          CALL largest_eig_timer(1)

          !Check if no errors in ZNAUPD
          IF (info .NE. 0) THEN
             WRITE(*,*) "Error code in znaupd. Code:", info, " irevcm ", irevcm
             WRITE(*,*) "ARNOLDI(ARPACK): ncv, nev, Ndim", ncv, nev, Ndim
             EXIT arnoldi_loop
          end IF

          !Apply op to vec, or exit, or throw error
          IF(irevcm .EQ. 99) THEN

             !Once irevcm = 99, it means convergence has been achieved -- Exit Arnoldi Loop
             EXIT arnoldi_loop

          ELSEIF(irevcm .EQ. -1 .OR. irevcm .EQ. 1) THEN

             !Apply Linear operator to vector
             CALL app_op_to_vec(workd(ipntr(1):ipntr(1)+Ndim-1), workd(ipntr(2):ipntr(2)+Ndim-1), which_evec)
             CALL largest_eig_timer(2)

          ELSE 
             WRITE(*,*) "Unexpected request in Reverse communication"
          end IF

       end DO arnoldi_loop

       IF(info .EQ. -8) THEN
          !Apparently the root cause of this error code is asking for too many evals.
          CALL update_nev(1, nev, max_nev, -1, nev_unchanged)
          cease_growing_nev=.TRUE.
          
          IF(nev_unchanged) THEN
             WRITE(*,*) "Continued to get error code 7 and cannot decrease nev"
             STOP
          end IF
          CYCLE nevloop
       end IF

       WRITE(*,*) "num of znaupd iters: ", iparam(3) 
       
       !Use ZNEUPD for post processing to find evals && evecs (matrix V). Note that we reuse V. 
       !NB. 'P' is Schur, we use this because Ritz does not guarantee orthogonality of degenerate vectors -- we use Schur as it is orthogonal.
       CALL ZNEUPD(.TRUE., 'P', selection, evals, v, Ndim, shift, workev, &
            & bmat, Ndim, which, nev, tol, resid, ncv, v, Ndim, &
            & iparam, ipntr, workd, workl, lworkl, rwork, info)

       !Check if no errors in ZNEUPD
       IF(info .NE. 0) THEN
          WRITE(*,*) "Error code in zneupd. Code: ", info
          STOP
       end IF

       !Number of converged eigenvalues
       nconv=iparam(5)

       !Check if we have enough evals, if not repeated with incremeneted size.  
       !If we are told to stop growing, we do not do this, and just exit.
       IF((.NOT. cease_growing_nev) .AND. MINVAL(ABS(evals(1:nconv))) .GT. tol_factor*macheps) THEN

          CALL update_nev(1, nev, max_nev, +1, nev_unchanged)

          IF(nev_unchanged) THEN
             WRITE(*,'(" ARNOLDI(ARPACK) WARNING: nev=",I4" nconv=",I4, &
                  &" with smallest eval=",D12.5)') &
                  & nev, nconv, MINVAL(ABS(evals(1:nconv)))
             EXIT nevloop
          end IF
          
          ! Try again with the new value
          CYCLE nevloop
       ELSE
          ! If enough eigenvalues, end the loop.
          EXIT nevloop
       end IF

  END DO nevloop

  !Zero everything
  eval=0.0d0; evec = 0.0d0

  !Extract eval && evec
  eval = evals(1)
  evec = v(:,1)

  !WRITE(*,*)
  !WRITE(*,*) "ARNOLDI(ARPACK): eval = ", eval
  !WRITE(*,*)
  !WRITE(*,*) "ARNOLDI(ARPACK): evec = ", evec
  !WRITE(*,*)

  CALL largest_eig_timer(3)

  IF(DBG_ARNOLDI .GE. 2) &
         & WRITE(*,'(" ARNOLDI(ARPACK): Loop Timings: ",4(F6.2,X))')  timers
    
  !Record a target nev for the next time we start, by choosing whichever is smaller of the current value, or 5 more than
  !required evals.  The +1 is to make sure we always have one more if possible.  
  nev=MIN(nev, 5 + INT(1.1 * COUNT(ABS(evals(1:nconv)) .GT. macheps)))

  IF(DBG_ARNOLDI .GE. 1) &
         & WRITE(*,'(" ARNOLDI(ARPACK): Exiting, non-zero ev:",I3,"/",I3)') &
         & COUNT(ABS(evals(1:nconv)) .GT. macheps), max_nev

  DEALLOCATE(resid, v, ax, workd, workl, rwork, selection, workev)

  !Find the element with smallest complex phase, rescale evec by that element (to fix the gauge)
  minPhiLoc = MINLOC(REAL(ATAN(evec)), ABS(REAL(evec)) .GT. tol)
  evec = evec/evec(minPhiLoc(1))

  !Normalize evec:
  C_norm = SUM(CONJG(evec) * evec)
  evec = evec/sqrt(C_norm)


  !CALL roundTens(evec, 1.0d+12)

  !WRITE(*,*) "ARNOLDI C_norm: ", C_norm
  !WRITE(*,*) "ARNOLDI evec: ", evec

 
  !Check residuals
  !CALL app_op_to_vec(evec, evec_out, which_evec)

  !WRITE(*,*)
  !WRITE(*,*) "arnoldi_calc_largest_eig -- checking residuals: " 
  !WRITE(*,*) MAXVAL(ABS(eval*evec(:) - evec_out(:)))
  !WRITE(*,*)

 CONTAINS

    SUBROUTINE largest_eig_timer(n)
      INTEGER :: n
      CALL cpu_time(ctime)
      timers(n) = timers(n) + ctime-otime
      otime=ctime
    end SUBROUTINE largest_eig_timer
    
 END SUBROUTINE OLD_mps_arnoldi_largest_eig



  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NAME:   arpack_arnoldi_svd
  ! VARIABLES:
  !           app_op_to_vec - Routine to apply matrix to vector
  !           u             - SVD left
  !           sigma         - SVD values.
  !           vh            - SVD right
  !           nev_est       - OPTIONAL: Input: number of ev to start
  !                           by finding.  Output: number actually required.
  ! SYNOPSIS:
  ! Aim to find SVD, i.e. OP = U Sigma VH.  
  ! We can use OP V = U Sigma, so if we know V, Sigma we can find
  ! U.  We want to return U, Sigma, VH.  We also use
  ! OP'= V Sigma UH,  OP' OP V = V Sigma^2
  ! so V is eigenvectors of (OP'OP).  We must at the end invert to
  ! get VH.
  ! nev_est is used as input to give a number of ev (less than storage size)
  ! that we try to find.  On return it is an appropriate value for the nxt
  ! such call of this routine.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE mpo_arnoldi_svd(app_op_to_vec, U, sigma, VH, nev_est)

    USE arnoldi_utility, ONLY: check_dim_match, update_nev, macheps, tol_factor

    INTERFACE
       SUBROUTINE app_op_to_vec(vi, vo, hc)
         COMPLEX(KIND(1.0D0)), INTENT(IN) :: vi(:)
         COMPLEX(KIND(1.0D0)), INTENT(OUT) :: vo(:)
         LOGICAL, INTENT(IN), OPTIONAL :: hc
       end SUBROUTINE app_op_to_vec
    end INTERFACE

    REAL (KIND=DP), INTENT(OUT) :: sigma(:)
    COMPLEX (KIND=DP), INTENT(OUT) :: U(:,:), vh(:,:)
    INTEGER, OPTIONAL, INTENT(INOUT) :: nev_est

    ! Number of eigenvalues we are trying to find.  (Read in from
    ! nev_est if present).
    INTEGER :: nev
    ! Dimension of problem (of matrix), max number singular vals
    INTEGER :: max_nev, ncv, lworkl, nconv, ev, ev_sel
    ! If reseting nev, record whether it unchanged.  If we had
    ! to shrink nev this time, we do not let it grow again,
    ! on this step, as otherwise we get into tight loops.
    LOGICAL :: nev_unchanged, cease_growing_nev
    INTEGER :: irevcm, info, iparam(11), ipntr(14)
    INTEGER:: ldim, rdim, i
    CHARACTER :: bmat*1, which*2

    COMPLEX (KIND=DP), ALLOCATABLE :: resid(:), v(:,:), workd(:), ax(:)
    COMPLEX (KIND=DP), ALLOCATABLE :: evals(:), workl(:), workev(:)
    LOGICAL, ALLOCATABLE :: selection(:)
    REAL (KIND=DP), ALLOCATABLE :: rwork(:)
    REAL (KIND=DP) :: tol, shift, norm

    REAL (KIND=DP) :: timers(4), otime, ctime

    EXTERNAL :: znaupd

    !  Check consistency of all dimensions.  Number of eigen(signular)
    ! values set by size of sigma.
    max_nev=SIZE(sigma)
    CALL CHECK_DIM_MATCH(SIZE(U,2), max_nev, "U/Sigma")
    CALL CHECK_DIM_MATCH(SIZE(VH,1), max_nev, "Sigma/VH")
    

    ! Default to maximum number of evals if no estimate given.  In all
    ! cases make sure we do not exceeed the maximum.
    IF(present(nev_est)) THEN
       nev=MIN(max_nev, nev_est) 
    ELSE
       nev=max_nev
    end IF

    cease_growing_nev=.FALSE.

    nevloop:  DO 
       ldim=SIZE(U,1); rdim=SIZE(VH,2)


       ! Note that since we solve A'A, the relevant dimension comes from
       ! rdim, for what we actually do.  Perhaps it would be smarter to
       ! choose whether we solve A'A or AA' dependent on which is
       ! smaller.  However in practice our matrices are generally square so
       ! it doesn't actually matter.

       ! Number of Arnoldi vectors
       ncv=MIN(3*nev+1,rdim)

       ! Work array size
       lworkl=4*ncv*(ncv+5)
       ! Tolerance in solver
       tol=1.0D-5
       ! Shift, unused as simple problem
       shift=0.0D0


       ! In case we had to re-run after changing nev, need to free this up
       IF(ALLOCATED(resid)) DEALLOCATE(Resid,v,ax,evals,workd,workl,rwork, &
            & selection, workev)

       ! Work spaces etc.  ax, used in multiplication.
       ALLOCATE(resid(rdim), v(rdim,ncv),ax(ldim),evals(nev+1))
       ALLOCATE(workd(3*rdim), workl(lworkl), rwork(ncv))
       ALLOCATE(selection(ncv), workev(2*ncv))

       ! Reverse communication flags etc
       irevcm=0
       info=0           ! Reset residuals.

       iparam(1)=1      ! Use exact shifts
       iparam(3)=10*nev ! Max iterations
       iparam(7)=1      ! Mode 1, i.e. simple, not generalized.
       bmat='I'         ! Also required for mode1.
       which='LM'       ! Want largest magnitude vals.
       selection=.TRUE.  ! Unreferenced as we ask for all.

       CALL cpu_time(otime); timers=0
       IF (DBG_ARNOLDI .GE. 1) &
            & WRITE(*,'(" ARNOLDI(ARPACK): nev=",I3," rdim=",I5)') nev, rdim

       arnoldi_loop: do
          CALL znaupd(irevcm, &
               & bmat, rdim, which, nev, tol, resid, ncv, v, rdim,  &
               & iparam, ipntr, workd, workl, lworkl, rwork, info)

          CALL timer(1)

          IF (info .NE. 0) THEN
             WRITE(*,*) "Error code in znaupd. Code:", info, " irevcm ",irevcm
             EXIT arnoldi_loop
          end IF

          IF(irevcm .EQ. 99) THEN
             EXIT arnoldi_loop
          ELSE IF (irevcm .EQ. -1 .OR. irevcm .EQ. 1) THEN
             ! Use internal subroutines below.  These must be internal
             ! as we need callback to matrices L, R.  Two versions for
             ! normal Hermitian conjugate
             CALL app_op_to_vec(workd(ipntr(1):ipntr(1)+rdim-1), ax)
             CALL timer(2)
             CALL app_op_to_vec(ax, workd(ipntr(2):ipntr(2)+rdim-1),hc=.TRUE.)
             CALL timer(3)
          ELSE 
             WRITE(*,*) "Unexpected request in Reverse communication"
          end IF

       end do arnoldi_loop

       IF(info .EQ. -8) THEN
          ! Apparently the root cause of this error code is asking for
          ! too many evals.
          CALL update_nev(1, nev, max_nev, -1, nev_unchanged)
          cease_growing_nev=.TRUE.
          
          IF (nev_unchanged) THEN
             WRITE(*,*) "Continued to get error code 7 and cannot decrease nev"
             STOP
          end IF
          CYCLE nevloop
       end IF
       

       ! Post processing to find sigma, and v.  Note that we reuse V, and
       ! that now V is the correct order evecs matrix, not VH.  Note that
       ! There are too many columns.  NOte that 'P' is Schur, we use this
       ! because ritz does not guarantee orthogonality of degenerate
       ! vectors.  Schur is correct, as it is orthogonal, and since the
       ! matrix we diagonalise is hermitian (we diagaonalise A' A) then
       ! Schur is the eigenbasis
       CALL ZNEUPD(.TRUE., 'P', selection, evals, v, rdim, shift, workev, &
            & bmat, rdim, which, nev, tol, resid, ncv, v, rdim, &
            & iparam, ipntr, workd, workl, lworkl, rwork, info)

       ! Number of converged eigenvalues
       nconv=iparam(5)


       ! Check if we have enough evals, if not repeated with
       ! incremeneted size.  If we are told to stop growing, we do
       ! not do this, and just exit.
       IF((.NOT. cease_growing_nev) .AND. &
            & MINVAL(ABS(evals(1:nconv))) .GT. tol_factor*macheps) THEN
          CALL update_nev(1, nev, max_nev, +1, nev_unchanged)

          IF (nev_unchanged) THEN
             WRITE(*,'(" ARNOLDI(ARAPACK) WARNING: nev=",I4" nconv=",I4, &
                  &" with smallest eval=",D12.5)') &
                  & nev, nconv, MINVAL(ABS(evals(1:nconv)))
             EXIT nevloop
          end IF
          
          ! Try again with the new value
          CYCLE nevloop
       ELSE
          ! If enough eigenvalues, end the loop.
          EXIT nevloop
       end IF

    end DO nevloop

    ! Zero everything
    U=0.0D0
    VH=0.0D0
    sigma=0.0D0

    ! Number of physical eigenvalues, after discarding negative
    ! evals.
    ev_sel=0
    DO ev=1, nconv
       IF (REAL(evals(ev)) .GE. 0.0D0) THEN
          ! This is a possible eigenvalue, so we record it.
          ev_sel=ev_sel+1
          !!!!!!!!!!!!!!!!!!!!
          ! Record sigma from these
          sigma(ev_sel)=SQRT(REAL(evals(ev)))

          !!!!!!!!!!!!!!!!!!!!
          ! Using V to find U, via: OP_{ij} V_j = U_i Sigma_i.
          ! We use normalisation to fix U scale.
          CALL app_op_to_vec(v(:,ev), U(:, ev_sel))

          norm=SQRT(SUM(ABS(U(:, ev_sel))**2))
          U(:, ev_sel)=U(:, ev_sel)/norm

          !!!!!!!!!!!!!!!!!!!!
          ! Record VH by hermitian conjugate
          vh(ev_sel, :) = CONJG(v(:,ev))
       end IF
       
    end DO

    CALL timer(4)

    IF (DBG_ARNOLDI .GE. 2) &
         & WRITE(*,'(" ARNOLDI(ARPACK): Loop Timings: ",4(F6.2,X))')  timers
    
    ! Record a target nev for the next time we start, by choosing
    ! whichever is smaller of the current value, or 5 more than
    ! required evals.  The +1 is to make sure we always have one more
    ! if possible.  
    nev=MIN(nev, 5 + INT(1.1 * COUNT(ABS(sigma(1:nconv)) .GT. macheps)))

    ! If we passed nev_est, we return the estimate to use for the next
    ! call, so it may be returned again.
    IF (PRESENT(nev_est)) nev_est=nev

    IF(DBG_ARNOLDI .GE. 1) &
         & WRITE(*,'(" ARNOLDI(ARPACK): Exiting, non-zero ev:",I3,"/",I3)') &
         & COUNT(ABS(sigma(1:nconv)) .GT. macheps), max_nev

    DEALLOCATE(resid, v, ax, workd, workl, rwork, selection,workev)

  CONTAINS

    SUBROUTINE timer(n)
      INTEGER :: n
      CALL cpu_time(ctime)
      timers(n) = timers(n) + ctime-otime
      otime=ctime
    end SUBROUTINE timer
    
  END SUBROUTINE mpo_arnoldi_svd

end MODULE arnoldi_wrappers
