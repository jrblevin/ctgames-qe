! model.f90 --- a simple quality ladder model
!
! Copyright (C) 2009-2024 Jason R. Blevins
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation  and/or other materials provided with the distribution.
! 3. Neither the names of the copyright holders nor the names of any
!    contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.

module model
  use osl, wp => osl_wp
  use dataset
  use encoding
  use stdlib_experimental_sparse
  use sparse
  implicit none

  private

  public :: nn, nw, we, mktsize, rho
  public :: np
  public :: maxiter, vf_tol, Q_tol, prmod, nvf
  public :: np1, np2
  public :: theta_lb, theta_ub
  public :: model_init
  public :: model_free
  public :: model_reset
  public :: model_set_theta
  public :: model_print
  public :: model_generate_data
  public :: model_dataset_stats
  public :: model_log_likelihood
  public :: DISCRETE_TIME, DELTA

  ! Model settings
  integer :: nn                  = 5              ! Maximum number of firms
  integer :: nw                  = 7              ! Number of quality levels
  real(wp) :: mktsize            = 2.5_wp         ! Market size
  integer :: we                  = 3              ! Entrant quality
  real(wp) :: rho                = 0.05_wp        ! Discount rate
  integer :: maxiter             = 1000           ! Maximum iterations
  integer :: prmod               = 0              ! Value function verbosity
  real(wp) :: vf_tol             = 1.0e-4_wp      ! Value function tolerance
  real(wp) :: Q_tol              = 1.0e-13_wp     ! Accuracy for expvunif

  ! Model constants
  integer,  parameter :: dim_s   = 2              ! State space dimension
  integer,  parameter :: na_inc  = 3              ! Incumbent actions
  integer,  parameter :: na_ent  = 2              ! Entrant actions
  integer,  parameter :: np      = 6              ! Number of parameters
  integer,  parameter :: np1     = 3              ! First-stage parameters
  integer,  parameter :: np2     = 3              ! Second-stage parameters
  real(wp), parameter :: nr_tol  = 1.0e-6_wp      ! Newton-Raphson tolerance
  real(wp), parameter :: exp_tol = 0.0_wp         ! Expokit tolerance
  real(wp), parameter :: eps     = 1.0e-5_wp      ! Finite-difference delta
  real(wp), parameter :: wstar   = 12.0_wp        ! Constant in g function
  real(wp), parameter :: mc      = 5.0_wp         ! Marginal cost
  real(wp), parameter :: log_max = 700.0_wp       ! Upper bound for log
  integer :: ns = 100                             ! Market structures
  integer :: nk = 0                               ! Number of feasible states
  logical, parameter :: vfi_unif = .false.        ! Apply uniformization to VFI

  ! Upper and lower parameter bounds
  real(wp), parameter, dimension(np) :: theta_lb = [ 0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp ]
  real(wp), parameter, dimension(np) :: theta_ub = [ 5.00_wp, 5.00_wp, 5.00_wp, 5.00_wp, 10.0_wp, 10.00_wp ]

  ! State of the model
  real(wp), dimension(:), allocatable :: v        ! Value function
  real(wp), dimension(:), allocatable :: v_new    ! Value function
  real(wp), dimension(:), allocatable :: pi       ! Flow profit
  real(wp), dimension(:), allocatable :: Emax     ! Emax function
  real(wp), dimension(:,:), allocatable :: ccpi   ! Incumbent choice prob.
  real(wp), dimension(:,:), allocatable :: ccpe   ! Entrant choice prob.
  real(wp), dimension(:,:), allocatable :: ccpi_old ! Previous incumbent choice prob.
  real(wp), dimension(:,:), allocatable :: ccpe_old ! Previous entrant choice prob.

  ! Value function history
  integer :: nvf = 10                             ! Number of v.f. to store
  integer :: ivf                                  ! Index of current v.f.
  real(wp), dimension(:,:), allocatable :: th_old ! Previous parameter values
  real(wp), dimension(:,:), allocatable :: v_old  ! Previous value functions

  ! Internal parameters
  real(wp), dimension(np) :: theta0               ! True parameters
  real(wp) :: lambda_low, lambda_high             ! Named parameters
  real(wp) :: gamma, kappa, eta, fc
  real(wp), parameter :: phi = 0.0_wp             ! Scrap value

  ! State space
  integer, dimension(:,:), allocatable :: is_inv  ! States after investment
  integer, dimension(:,:), allocatable :: is_exit ! States after exit
  integer, dimension(:), allocatable :: is_ent    ! States after entry
  integer, dimension(:), allocatable :: is_dep    ! States after depreciation
  integer, dimension(:,:), allocatable :: lookup  ! Fast decoding table
  integer, dimension(:,:), allocatable :: enc     ! (is,wi) encoding matrix
  integer, dimension(:,:), allocatable :: dec     ! (is,wi) decoding matrix
  integer, dimension(:), allocatable :: k_inv     ! Own investment transitions
  integer, dimension(:), allocatable :: k_dep     ! Depreciation transitions
  integer, dimension(:), allocatable :: k_ent     ! Entrant transitions
  integer, dimension(:,:), allocatable :: k_rinv  ! Rival firm investment
  integer, dimension(:,:), allocatable :: k_rexit ! Rival exit transitions
  integer, dimension(:), allocatable :: k_rent    ! Rival entry transitions
  logical, dimension(:), allocatable :: k_low     ! Low quality indicator
  integer, dimension(:,:), allocatable :: nriv    ! Num. rivals in each state

  ! State Vector Components
  integer, parameter :: k_is = 1                  ! Index of is
  integer, parameter :: k_wi = 2                  ! Index of wi

  ! Named player positions
  integer :: WI_ENTRANT                           ! Potential chain entrant
  integer, parameter :: WI_NATURE = 0             ! Nature
  integer, parameter :: WI_DISCRETE = -1          ! Invalid for discrete-time

  ! Generic actions
  integer, parameter :: A_DISCRETE           = -1 ! Invalid for discrete-time

  ! Dataset
  logical :: DISCRETE_TIME = .false.
  real(wp) :: DELTA = 1.0_wp
  type(dataset_t) :: dat
  type(spmatrix_t) :: dat_counts
  integer, dimension(:), allocatable :: P_state
  type(osl_rng_t), save :: rng_data               ! Random number generator for DGP

  ! Timer variables
  integer :: count0, count1, count_rate
  real(wp) :: time, time0, time1

  ! Sparse matrices for calculating Q
  integer, dimension(:), allocatable :: QCSCp     ! Q matrix in CSC format
  integer, dimension(:), allocatable :: QCSCi     ! Q matrix in CSC format
  real(wp), dimension(:), allocatable :: QCSCx    ! Q matrix in CSC format
  integer, dimension(:), allocatable :: Qi        ! Q matrix in COO format
  integer, dimension(:), allocatable :: Qj        ! Q matrix in COO format
  integer :: Q_nz                                 ! Number of nonzeros
  integer, dimension(:), allocatable :: Q_gamma   ! Q_gamma is kept in vector format
  logical, dimension(:), allocatable :: Q_diag    ! Mask for diagonal elements
  real(wp) :: max_q_ii                            ! Uniformization rate
  type(spmatrix_t) :: Q_ccpi2, Q_ccpi3, Q_ccpe2   ! Q_ccpi2, Q_ccpi3, Q_ccpe2 are CSR

  ! Print human readable state information
  interface print_state
     module procedure print_state_int
     module procedure print_state_vec
  end interface print_state

contains

  !------------------------------------!
  ! Model lifecycle                    !
  !------------------------------------!

  ! Initialize the model and allocate any necessary memory.
  subroutine model_init()
    integer :: k, is, isp, wi, wip, wir
    integer, dimension(nw+1) :: s, sp
    real(wp), dimension(nw) :: profit

    print '(a)', 'Initializing the model...'

    ! nw + 1 support points, maximum value of nn for each point
    print '(a)', 'Setting up the state-space encoding...'
    ns = encoding_init(nw+1, nn+1)

    ! Named player positions.
    WI_ENTRANT     = nw + 1

    ! Generate a lookup table for fast decoding
    print '(a)', 'Generating a fast lookup table for market configurations (s)...'
    allocate(lookup(nw+1,ns))
    do is = 1, ns
       lookup(:,is) = decode(is)
    end do

    ! Find and encode all feasible (s, x, omega) combinations.
    ! Elements of enc equal to zero indicate infeasible states.
    print '(a)', 'Generating an encoding table for all (s, omega) combinations...'
    allocate(enc(ns, nw))
    nk = 0
    enc = -1  ! For trapping invalid values
    do is = 1, ns
       s = lookup(:,is)
       do wi = 1, nw
          if (s(wi) > 0) then
             nk = nk + 1
             enc(is,wi) = nk
          else
             enc(is,wi) = 0
          end if
       end do
    end do
    print '(a,i0)', 'Total states: ', nk

    ! Check for invalid states
    print '(a)', 'Checking for invalid states...'
    if (any(enc < 0) .or. any(enc > nk)) then
       stop 'ERROR: enc contains invalid states [model_init]'
    end if

    ! Build the decoding table
    print '(a)', 'Generating a fast decoding table mapping k to (is, i)...'
    allocate(dec(nk,dim_s))
    dec = -1  ! For trapping invalid values
    do is = 1, ns
       do wi = 1, nw
          k = enc(is,wi)
          if (k > 0) then
             dec(k,k_is) = is
             dec(k,k_wi) = wi
          end if
       end do
    end do

    ! Check for invalid states
    print '(a)', 'Checking for invalid states...'
    if (any(dec < 0) .or. any(dec(:,k_is) > ns) .or. any(dec(:,k_wi) > nw+1)) then
       stop 'ERROR: dec contains invalid states [model_init]'
    end if

    ! Store profits
    print '(a)', 'Precalculating profits...'
    allocate(pi(nk))
    pi = 0.0_wp
    !$OMP PARALLEL DO PRIVATE(profit, is, wi)
    do is = 1, ns
       profit = model_profit(is)
       do wi = 1, nw
          if (enc(is,wi) > 0) then
             pi(enc(is,wi)) = profit(wi)
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    ! Allocate memory and choose an arbitrary starting value for V
    print '(a)', 'Allocating memory for value function storage...'
    allocate(v(nk))
    allocate(v_new(nk), Emax(nk))
    v = pi

    ! Initialize value function history storage
    if (nvf > 0) then
       ! We use NP+1 so as to include PHI and FC
       allocate(v_old(nk, nvf), th_old(np+1, nvf))
       v_old = 0.0_wp
       th_old = 0.0_wp
       v_old(:,1) = pi
       ivf = 1
    end if

    ! Precalculate market state transitions
    print '(a)', 'Precalculating market structure (s) transitions...'
    allocate(is_inv(ns,nw), is_exit(ns,nw), is_ent(ns), is_dep(ns))
    is_inv = 0
    is_exit = 0
    is_ent = 0
    is_dep = 0
    do is = 1, ns
       s = lookup(:,is)

       do wi = 1, nw
          if (s(wi) > 0) then
             ! Investment
             sp = s
             sp(wi) = sp(wi) - 1
             sp(min(wi+1, nw)) = sp(min(wi+1, nw)) + 1
             is_inv(is, wi) = encode(sp)

             ! Exit
             sp = s
             sp(wi) = sp(wi) - 1
             sp(WI_ENTRANT) = sp(WI_ENTRANT) + 1
             is_exit(is, wi) = encode(sp)
          end if
       end do

       ! Entry
       if (s(WI_ENTRANT) > 0) then
          sp = s
          sp(WI_ENTRANT) = sp(WI_ENTRANT) - 1
          sp(we) = sp(we) + 1
          is_ent(is) = encode(sp)
       end if

       ! Market-wide depreciation (must have active firms above level 1).
       ! Each firm moves down one level except for firms at level 1.
       ! If not, the depreciation event still happens but the state is
       ! unchanged.
       if (s(WI_ENTRANT) < nn .and. s(1) < nn) then
          sp = s
          sp(1) = s(1) + s(2)
          do wi = 2, nw - 1
             sp(wi) = s(wi+1)
          end do
          sp(nw) = 0
          sp(nw+1) = s(nw+1)
          is_dep(is) = encode(sp)
       else
          is_dep(is) = is
       end if
    end do

    ! Allocate memory for storing CCPs and for CCP estimation
    print '(a)', 'Allocating memory for CCP storage...'
    allocate(ccpi(na_inc, nk), ccpe(na_ent, ns))
    allocate(ccpi_old(na_inc,nk), ccpe_old(na_ent,ns))

    ! Precalculate state transitions and number of rivals
    print '(a)', 'Precalculating state space (k) transitions...'
    allocate(k_inv(nk), k_dep(nk), k_ent(ns))
    allocate(k_rinv(nk, nw), k_rexit(nk, nw), k_rent(nk), k_low(nk))
    allocate(nriv(nk, nw+1))
    k_inv = -1
    k_rexit = -1
    k_rent = -1
    k_low = .false.
    k_rinv = -1
    k_dep = -1
    k_ent = -1
    nriv = 0
    do k = 1, nk
       ! Decode state K
       is = dec(k,k_is)
       wi = dec(k,k_wi)
       s = lookup(:,is)

       ! Low rate indicator
       k_low(k) = (wi < we)

       ! Number of rivals of each quality for a firm in state K.
       nriv(k,:) = s
       nriv(k,wi) = nriv(k,wi) - 1                ! Omit self
       nriv(k,nw+1) = min(nriv(k,nw+1), 1)        ! Limit to one rival

       ! Own investment: a firm of quality WI invests in state S.
       ! WI increases one step as the firm moves from WI to min(WI+1,NW)
       wip = min(wi+1, nw)
       isp = is_inv(is, wi)
       k_inv(k) = enc(isp, wip)

       ! Investment: a rival firm of quality WIR invests in state S.
       do wir = 1, nw
          if (nriv(k, wir) > 0) then
             ! WI stays the same while a rival moves to min(WIR+1,NW)
             wip = wi
             isp = is_inv(is, wir)
             k_rinv(k, wir) = enc(isp, wip)
          else
             k_rinv(k, wir) = 0
          end if
       end do

       ! Exit: a rival firm of quality I exits in state S.
       do wir = 1, nw
          if (nriv(k, wir) > 0) then
             ! WI stays the same while a rival firm moves from state I to NW+1.
             wip = wi
             isp = is_exit(is, wir)
             k_rexit(k, wir) = enc(isp, wip)
          else
             k_rexit(k, wir) = 0
          end if
       end do

       ! Entry: a potential entrant enters at state WE.
       if (nriv(k, nw + 1) > 0) then
          ! WI stays the same while a single firm moves from NW+1 to WE.
          wip = wi
          isp = is_ent(is)
          k_rent(k) = enc(isp, wip)
       else
          k_rent(k) = 0
       end if

       ! Depreciation: firms moves down one level except for firms at level 1.
       ! If there are no active firms in states 2,...,NW then the state is
       ! unchanged.
       if (s(nw+1) < nn .and. s(1) < nn) then
          wip = max(wi - 1, 1)
          isp = is_dep(is)
          k_dep(k) = enc(isp, wip)
       else
          k_dep(k) = k
       end if
    end do
    do is = 1, ns
       ! Potential entrant
       if (is_ent(is) > 0) then
          k_ent(is) = enc(is_ent(is), we)
       else
          k_ent(is) = 0
       end if
    end do

    ! Precalculate addresses for sparse matrix storage.
    ! Also return exactly the number of nonzeros.
    print '(a)', 'Precalculating Q matrix sparsity pattern...'
    call Q_COO_precalc()
  end subroutine model_init


  ! Resets the model to a fresh initial state, storing the "true"
  ! parameters, and resolving the model.
  subroutine model_reset(theta, verbose)
    real(wp), dimension(np), intent(in) :: theta
    logical, intent(in), optional :: verbose
    logical :: verbose_

    if (present(verbose)) then
      verbose_ = verbose
    else
      verbose_ = .false.
    end if

    ! Store theta for later
    theta0 = theta

    ! Store each component of theta
    call model_set_theta(theta)

    ! Obtain the value function and CCPs
    v = pi ! Simple guess at initial v
    if (verbose_) then
      print '("Obtaining value function...")'
      call tic()
    end if
    call model_vf_vfi()
    if (verbose_) then
        call toc()
    end if
    call model_ccp()

    ! Tests
    call Q_update()
  end subroutine model_reset


  ! Free memory allocated for the model.
  subroutine model_free()
    deallocate(Q_diag)
    call spmatrix_free(Q_ccpi2)
    call spmatrix_free(Q_ccpi3)
    call spmatrix_free(Q_ccpe2)
    call dataset_free(dat)
    deallocate(v_new, Emax)
    deallocate(v, pi)
    deallocate(ccpe_old, ccpi_old)
    deallocate(ccpi, ccpe)
    deallocate(is_inv, is_exit, is_ent, is_dep)
    deallocate(k_inv, k_dep, k_ent)
    deallocate(k_rinv, k_rexit, k_rent, k_low)
    deallocate(lookup)
    deallocate(nriv)
    deallocate(enc)
    deallocate(dec)
    call encoding_free()
  end subroutine model_free


  ! Sets the internal values of all parameters
  subroutine model_set_theta(theta)
    real(wp), dimension(np), intent(in) :: theta

    lambda_low = theta(1) ! Low firm move arrival rate
    lambda_high = theta(2) ! High firm move arrival rate
    gamma = theta(3) ! Market depreciation rate
    kappa = theta(4) ! Investment cost
    eta = theta(5) ! Entry cost
    fc = theta(6) ! Exit cost
  end subroutine model_set_theta


  ! Prints the parameter values with labels
  subroutine model_print()
    print *, 'nn: ', nn
    print *, 'nw: ', nw
    print *, 'ns: ', ns
    print *, 'nk: ', nk
    print *, 'we: ', we
    print *, 'mktsize: ', mktsize
    print *
    print *, 'nr_tol: ', nr_tol
    print *, 'vf_tol: ', vf_tol
    print *, 'Q_tol: ', Q_tol
    print *
    print *, 'lambda_L: ', lambda_low, ' (low move arrival rate)'
    print *, 'lambda_H: ', lambda_high, ' (high move arrival rate)'
    print *, 'gamma:  ', gamma, ' (market depreciation rate)'
    print *, 'kappa:  ', kappa, ' (investment cost)'
    print *, 'eta:    ', eta, ' (entry cost)'
    print *, 'phi:    ', phi, ' (scrap value)'
    print *, 'fc:     ', fc, ' (fixed cost)'
    print *
    if (DISCRETE_TIME) then
       print *, 'Data Sampling: Discrete Time'
       print *, 'Delta: ', DELTA
    else
       print *, 'Data Sampling: Continuous Time'
    end if
    print *, 'Estimation Method: Maximum Likelihood'
  end subroutine model_print


  ! Solve for the value function using value function iteration.
  subroutine model_vf_vfi()
    real(wp), dimension(na_inc) :: vp
    integer :: iter, k, kr, wi, wir, nr, is
    real(wp) :: norm, rate, ratesum, lambda_self, unif_rate

    ! Retrieve value function from a nearby parameter vector
    call model_retrieve_vf()

    ! Update CCPs for initial V
    call model_ccp()

    ! Calculate the uniform rate
    unif_rate = gamma + real(nn, wp) * max(lambda_low, lambda_high) + 0.1_wp

    ! Apply value function iteration until the sup norm of CCPs is
    ! below VF_TOL or more than MAXITER iterations have been performed.
    norm = vf_tol + 1.0_wp
    iter = 0
    if (prmod > 0) then
       print '(a)', 'Performing value function iteration...'
       print '(a8,a17)', 'iter', 'norm'
    end if
    do while ((norm > vf_tol) .and. (iter < maxiter))
       iter = iter + 1

       !$OMP PARALLEL DO PRIVATE(wi, wir, is, kr, nr, rate, ratesum, vp, lambda_self)
       do k = 1, nk
          is = dec(k, k_is)
          wi = dec(k, k_wi)
          ratesum = 0.0_wp

          ! Flow payoff
          if (enc(is,wi) > 0) then
             v_new(k) = pi(k) - fc
          else
             v_new(k) = 0.0_wp
          end if

          ! Market depreciation
          if (k_dep(k) > 0) then
             v_new(k) = v_new(k) + gamma * v(k_dep(k))
             ratesum = ratesum + gamma
          end if

          ! Rival incumbents
          do wir = 1, nw
             nr = nriv(k, wir)
             if (nr > 0) then
                kr = enc(is, wir)                       ! rival state
                ! Total rate at which nr rivals of quality wir make decisions
                if (wir <= we) then
                   rate = real(nr, wp) * lambda_low
                else
                   rate = real(nr, wp) * lambda_high
                end if
                vp(1) = v(k)                            ! nothing
                vp(2) = v(k_rinv(k, wir))               ! invest
                vp(3) = v(k_rexit(k, wir))              ! exit
                v_new(k) = v_new(k) + rate * dot_product(ccpi(:,kr), vp)
                ratesum = ratesum + rate
             end if
          end do

          ! Potential entrants
          if (nriv(k, nw + 1) > 0) then
             v_new(k) = v_new(k) + lambda_low * &
                  (ccpe(1, is) * v(k) + ccpe(2, is) * v(k_rent(k)))
             ratesum = ratesum + lambda_low
          end if

          ! Own moves
          vp(1) = v(k)                              ! nothing
          vp(2) = v(k_inv(k)) - kappa               ! invest
          vp(3) = phi                               ! exit
          Emax(k) = log_sum_exp(vp) + OSL_EULER
          if (wi <= we) then
             lambda_self = lambda_low
          else
             lambda_self = lambda_high
          end if
          v_new(k) = v_new(k) + lambda_self * Emax(k)
          ratesum = ratesum + lambda_self

          ! Uniformization
          if (vfi_unif) then

            if (unif_rate < ratesum) then
               call osl_print('Error: The uniformization rate is too low.')
               call osl_print('unif_rate = ', unif_rate)
               call osl_print('ratesum = ', ratesum)
               call print_k(k)
               stop
            end if

            v_new(k) = v_new(k) + (unif_rate - ratesum) * v(k) ! Self transition
            ratesum = unif_rate
          end if

          ! Normalize
          v_new(k) = v_new(k) / (rho + ratesum)
       end do
       !$OMP END PARALLEL DO

       ! Store previous CCPs
       ccpi_old = ccpi
       ccpe_old = ccpe

       ! Update CCPs for current V
       v = v_new
       call model_ccp()

       ! Convergence criteria
       !$OMP WORKSHARE
       norm = max(maxval(abs(ccpi - ccpi_old)), maxval(abs(ccpe - ccpe_old)))
       !$OMP END WORKSHARE

       ! Print results of iteration
       if (prmod > 0) then
          if (mod(iter, prmod) == 0) then
             print '(i8,g17.5)', iter, norm
          end if
       end if

    end do

    if (iter >= maxiter) then
       call osl_print_heading('Error')
       call model_print()
       call osl_print('iter = ', iter)
       call osl_print('norm = ', norm)
       call osl_print('maxiter = ', maxiter)
       call osl_print('theta0 = ', theta0)
       stop 'model_vf_vfi: value function did not converge!'
    end if

    ! Store this value function and parameter vector
    call model_store_vf()
  end subroutine model_vf_vfi


  ! Store the current value function and parameter vector.
  subroutine model_store_vf()
    real(wp), dimension(np+1) :: theta

    ! We include normalized parameters (phi/fc) here, only so that
    ! switching between normalizations is easier.
    theta = [ lambda_low, lambda_high, gamma, kappa, eta, phi, fc ]
    if (nvf > 0) then
       th_old(:,ivf) = theta
       v_old(:,ivf) = v
       ivf = ivf + 1
       if (ivf == nvf + 1) then
          ivf = 1
       end if
    end if
  end subroutine model_store_vf


  ! Retrieve the stored value function for the nearest parameter vector.
  subroutine model_retrieve_vf()
    real(wp), dimension(nvf) :: d
    real(wp), dimension(np+1) :: theta
    integer :: i, imin

    ! We include normalized parameters (phi/fc) here, only so that
    ! switching between normalizations is easier.
    theta = [ lambda_low, lambda_high, gamma, kappa, eta, phi, fc ]
    if (nvf > 0) then
       ! If value function storage is implemented
       !$OMP PARALLEL DO
       do i = 1, nvf
          d(i) = osl_norm_2(th_old(:,i) - theta)
       end do
       !$OMP END PARALLEL DO
       !$OMP WORKSHARE
       imin = minloc(d, DIM=1)
       v = v_old(:,imin)
       !$OMP END WORKSHARE
    end if
  end subroutine model_retrieve_vf

  !-------------------------------------!
  ! Model implications                  !
  !-------------------------------------!

  ! Calculate and store conditional choice probabilities for incumbent
  ! firms and potential entrants in each state.
  subroutine model_ccp()
    ! Incumbents
    ccpi = transpose(model_ccpi_vec())
    ! Potential entrants
    ccpe = transpose(model_ccpe_vec())
  end subroutine model_ccp

  ! Return incumbent conditional choice probabilities for all states.
  function model_ccpi_vec() result(ccpi)
    real(wp), dimension(nk, na_inc) :: ccpi
    real(wp), dimension(nk) :: expsum

    !$OMP WORKSHARE
    ! Choice-specific value functions
    ccpi(:,1) = v
    ccpi(:,2) = v(k_inv) - kappa
    ccpi(:,3) = phi

    ! Recenter for each state
    ccpi = ccpi - spread(maxval(ccpi, 2), 2, na_inc)

    ! Logit probabilities
    ccpi = exp(ccpi)
    expsum = sum(ccpi, 2)
    ccpi(:, 3) = ccpi(:, 3) / expsum
    ccpi(:, 2) = ccpi(:, 2) / expsum
    ccpi(:, 1) = 1.0_wp - ccpi(:, 3) - ccpi(:, 2)
    !$OMP END WORKSHARE
  end function model_ccpi_vec

  ! Return entrant conditional choice probabilities for all states.
  function model_ccpe_vec() result(ccpe)
    real(wp), dimension(ns, na_ent) :: ccpe

    !$OMP WORKSHARE
    where (k_ent > 0)
       ccpe(:,1) = exp(eta - v(k_ent))
       ccpe(:,1) = ccpe(:,1) / (1 + ccpe(:,1))
       ccpe(:,2) = 1.0_wp - ccpe(:,1)
    elsewhere
       ccpe(:,1) = 0.0_wp
       ccpe(:,2) = 0.0_wp
    end where
    !$OMP END WORKSHARE
  end function model_ccpe_vec

  !------------------------------------!
  ! Likelihood Functions               !
  !------------------------------------!

  ! Full log likelihood function for estimation with continuous time data.
  subroutine model_log_likelihood_continuous(theta, ll)
    real(wp), dimension(:), intent(in) :: theta
    real(wp), intent(out) :: ll
    real(wp), dimension(nk, na_inc) :: ccpi
    real(wp), dimension(ns, na_ent) :: ccpe
    integer :: im, k, n_low, n_high
    real(wp) :: rate
    type(observation_t), pointer :: obs

    call model_set_theta(theta)
    call model_vf_vfi()

    ! Pre-calculate choice probabilities with current parameters
    ccpi = model_ccpi_vec()
    ccpe = model_ccpe_vec()

    ll = 0.0_wp
    !$OMP PARALLEL DO private(obs, n_low, n_high, rate, k) reduction(+:ll)
    do im = 1, dataset_n(dat)
       obs => dataset_first(dat, im)
       do while (associated(obs))

          ! Determine the number of active incumbents and entrants
          call firm_count(obs%is, n_low, n_high)
          rate = real(n_low, wp) * lambda_low + real(n_high, wp) * lambda_high

          ! Accumulate the likelihood in two cases: jumps and moves.
          if (obs%wi == WI_NATURE) then
             ! Jump -----------------------------------------------------------
             ! Expo(gamma)
             ll = ll + blog(osl_exponential_pdf(obs%tau, gamma))
             ! Pr(Tau > tau)
             ll = ll + blog(1.0_wp - osl_exponential_cdf(obs%tau, rate))
          else
             ! Move -----------------------------------------------------------
             ! Expo(rate)
             ll = ll + blog(osl_exponential_pdf(obs%tau, rate))
             ! Pr(Tau > tau)
             ll = ll + blog(1.0_wp - osl_exponential_cdf(obs%tau, gamma))
             ! Pr(a) Choice probability
             if (obs%wi == WI_ENTRANT) then
                ll = ll + blog(ccpe(obs%is, obs%a))   ! bounded log
             else
                k = enc(obs%is, obs%wi)
                ll = ll + blog(ccpi(k, obs%a))   ! bounded log
             end if
          end if

          ! Prepare for another loop
          obs => dataset_next(obs)
       end do
    end do
    !$OMP END PARALLEL DO
    ll = -ll / real(dataset_size(dat), wp)
  end subroutine model_log_likelihood_continuous


  ! Precalculate sparsity pattern for the aggregate intensity matrix
  ! in COO format.
  subroutine Q_COO_precalc
    integer :: is, wi, isp, k
    integer, dimension(nw+1) :: s
    integer :: nzmax, idx
    integer, dimension(:), allocatable :: Qtmp

    type(spmatrix_t) :: Q_ccpi2_coo, Q_ccpi3_coo, Q_ccpe2_coo

    ! Upper bound on number of non-zero Q matrix elements
    nzmax = ns * (3 + nn*3)
    print '(a)', 'Initializing temporary COO-format matrices...'

    allocate(Qi(nzmax), Qj(nzmax), Qtmp(nzmax), Q_gamma(nzmax))
    Q_gamma = 0

    ! NEEDS_WORK: These are using more space than necessary
    call spmatrix_init(Q_ccpi2_coo, nzmax, nk, nzmax, SPMATRIX_TYPE_COO) ! (nzmax x nk) Qx += Q_ccpi2 * ccpi
    call spmatrix_init(Q_ccpi3_coo, nzmax, nk, nzmax, SPMATRIX_TYPE_COO) ! (nzmax x nk) Qx += Q_ccpi3 * ccpi
    call spmatrix_init(Q_ccpe2_coo, nzmax, ns, nzmax, SPMATRIX_TYPE_COO) ! (nzmax x ns) Qx += Q_ccpe2 * ccpe

    idx = 0
    do is = 1, ns
       s = lookup(:,is)

       ! Market-wide depreciation (must have active firms above level 1).
       if (s(nw+1) < nn .and. s(1) < nn) then
          ! depreciation source
          idx = idx + 1
          Qi(idx) = is
          Qj(idx) = is

          ! Store coefficient in Q_gamma
          Q_gamma(idx) = -1

          ! depreciation sink
          idx = idx + 1
          isp = is_dep(is)
          Qi(idx) = is
          Qj(idx) = isp

          ! Store coefficient in Q_gamma
          Q_gamma(idx) = 1
       end if

       ! Incumbent firms
       do wi = 1, nw
          ! If there are no active chain firms in this state, move to the next.
          if (s(wi) == 0) cycle

          ! Encoded rival state
          k = enc(is, wi)

          ! Incumbent chain source
          idx = idx + 1
          Qi(idx) = is
          Qj(idx) = is

          ! Incumbent firms leave at rate
          ! -s(wi) * lambda * (1-ccpi(k,1)) = -s(wi) * lambda * (ccpi(k,2) + ccpi(k,3))
          call spmatrix_set(Q_ccpi2_coo, idx, k, real(-s(wi), wp))
          call spmatrix_set(Q_ccpi3_coo, idx, k, real(-s(wi), wp))

          ! Incumbent chain sink: investment
          idx = idx + 1
          isp = is_inv(is, wi)
          Qi(idx) = is
          Qj(idx) = isp

          ! Incumbent firms invest at rate s(wi) * lambda * ccpi(k,2)
          call spmatrix_set(Q_ccpi2_coo, idx, k, real(s(wi), wp))

          ! Incumbent chain sink: exit
          idx = idx + 1
          isp = is_exit(is, wi)
          Qi(idx) = is
          Qj(idx) = isp

          ! Incumbent firms exit at rate s(wi) * lambda * ccpi(k,3)
          call spmatrix_set(Q_ccpi3_coo, idx, k, real(s(wi), wp))
       end do

       ! Potental entrants
       wi = WI_ENTRANT
       if (s(wi) > 0) then
          ! Potential entrant source
          idx = idx + 1
          Qi(idx) = is
          Qj(idx) = is

          ! Potential entrants leave entrant status at rate -lambda * ccpe(is,2)
          call spmatrix_set(Q_ccpe2_coo, idx, is, -1.0_wp)

          ! Potential entrant sink: entry
          idx = idx + 1
          isp = is_ent(is)
          Qi(idx) = is
          Qj(idx) = isp

          ! Potential entrants enter market at rate lambda * ccpe(is,2)
          call spmatrix_set(Q_ccpe2_coo, idx, is, 1.0_wp)
       end if
    end do

    ! Reduce Q matrix vectors based on actual nonzero count
    Q_nz = idx

    ! Reallocate Qi
    Qtmp = Qi
    deallocate(Qi)
    allocate(Qi(Q_nz))
    Qi = Qtmp(1:Q_nz)

    ! Reallocate Qj
    Qtmp = Qj
    deallocate(Qj)
    allocate(Qj(Q_nz))
    Qj = Qtmp(1:Q_nz)

    ! Reallocate Q_gamma
    Qtmp = Q_gamma
    deallocate(Q_gamma)
    allocate(Q_gamma(Q_nz))
    Q_gamma = Qtmp(1:Q_nz)

    ! Update dimensions of spmatrix_t types before converting to CSR
    Q_ccpi2_coo%rows = Q_nz
    Q_ccpi3_coo%rows = Q_nz
    Q_ccpe2_coo%rows = Q_nz

    ! Calculate diagonal mask for Q elements
    allocate(Q_diag(Q_nz))
    Q_diag = (Qi == Qj)

    ! Convert Q_ccpi2, Q_ccpi3, Q_ccpe2 to CSR format.
    ! This is better for matrix-vector multiplication.
    call spmatrix_coo_to_csr(Q_ccpi2_coo, Q_ccpi2)
    call spmatrix_coo_to_csr(Q_ccpi3_coo, Q_ccpi3)
    call spmatrix_coo_to_csr(Q_ccpe2_coo, Q_ccpe2)
  end subroutine Q_COO_precalc

  ! Calculate the aggregate intensity matrix and store in CSC format.
  subroutine Q_update()
    real(wp), dimension(Q_nz) :: Qx ! Internal COO representation
    real(wp), dimension(nk, na_inc) :: ccpi
    real(wp), dimension(ns, na_ent) :: ccpe
    real(wp), dimension(nk) :: rate

    ! Calculate hazards under current parameters
    where (k_low)
      rate = lambda_low
    elsewhere
      rate = lambda_high
    end where
    ccpi = spread(rate, DIM=2, NCOPIES=na_inc) * model_ccpi_vec()
    ccpe = lambda_low * model_ccpe_vec()

    ! Operate on underlying sparse Q matrix value vector directly
    Qx = Q_gamma * gamma
    Qx = Qx + spmatrix_csr_matvec(Q_ccpi2, ccpi(:,2))
    Qx = Qx + spmatrix_csr_matvec(Q_ccpi3, ccpi(:,3))
    Qx = Qx + spmatrix_csr_matvec(Q_ccpe2, ccpe(:,2))

    ! Calculate maximal diagonal element of Q (in absolute value)
    ! All elements are negative, so we find the minimum instead.
    max_q_ii = -minval(Qx, DIM=1, MASK=Q_diag)

    ! Convert Q to CSC format
    call coo2csc_canonical(Qi, Qj, Qx, QCSCp, QCSCi, QCSCx)
  end subroutine Q_update

  ! Log likelihood function for discretely-sample data.
  ! Uses uniformization algorithm for computing matrix exponential.
  subroutine model_log_likelihood_discrete(theta, ll)
    real(wp), dimension(:), intent(in) :: theta
    real(wp), intent(out) :: ll

    real(wp), dimension(ns,size(P_state,1)) :: basis
    real(wp), dimension(:,:), allocatable :: P
    real(wp) :: count
    integer :: i, j, m_eps
    integer :: is, isp

    ! Set parameters and obtain value function
    call model_set_theta(theta)
    call model_vf_vfi()
    call Q_update()

    ! Allocate memory for transition matrix rows
    allocate(P(ns, size(P_state, 1)))

    ! Precalculate basis vectors
    basis = 0.0_wp
    do i = 1, size(P_state, 1)
       basis(P_state(i), i) = 1.0_wp
    end do

    ! Precalculate truncation term for expvunif
    m_eps = osl_poisson_inv(1.0_wp - Q_tol, max_q_ii)

    ! Calculate matrix exponential rows in parallel
    !$OMP PARALLEL DO
    do i = 1, size(P_state, 1)
       ! Apply uniformization algorithm to calculate is-th row of P
       P(:,i) = expvunif_sparse(basis(:,i), DELTA, max_q_ii, QCSCp, QCSCi, QCSCx, M=m_eps)
    end do
    !$OMP END PARALLEL DO

    ! Calculate the log likelihood in parallel
    ll = 0.0_wp
    !$OMP PARALLEL DO private(is, isp, count, i) reduction(+:ll)
    do j = 1, spmatrix_nz(dat_counts)
       is = dat_counts%row_idx(j)
       isp = dat_counts%col_idx(j)
       count = dat_counts%values(j)

       ! Determine location of row is
       i = 1
       do while (P_state(i) /= is)
          i = i + 1
       end do

       ! Accumulate the log-likelihood
       ll = ll + count * blog(P(isp,i))
    end do
    !$OMP END PARALLEL DO
    ll = -ll / real(dataset_size(dat), wp)
  end subroutine model_log_likelihood_discrete

  ! Log likelihood function
  subroutine model_log_likelihood(theta, ll, ctx)
    real(wp), dimension(:), intent(in) :: theta
    real(wp), intent(out) :: ll
    integer, dimension(:), intent(in), optional :: ctx

    if (present(ctx)) then
       stop 'CTX argument is currently unused'
    end if

    if (DISCRETE_TIME) then
       call model_log_likelihood_discrete(theta, ll)
    else
       call model_log_likelihood_continuous(theta, ll)
    end if
  end subroutine model_log_likelihood


  !------------------------------------!
  ! Profits                            !
  !------------------------------------!

  ! Calculates the profits for each firm for a given market structure
  ! assuming a differentiated products model with a logit demand
  ! system.
  function model_profit(is) result(pi_s)
    integer, intent(in) :: is                     ! market structure index
    real(wp), dimension(nw) :: pi_s               ! profits by quality
    real(wp), dimension(nn) :: pi_w               ! profits by firm
    real(wp), dimension(nn) :: p, p_eps, p_new    ! prices
    real(wp), dimension(nn,nn) :: J               ! Jacobian
    real(wp), dimension(nn) :: y, y_eps           ! foc at p and p_eps
    real(wp), dimension(nn) :: delta              ! descent direction
    real(wp), dimension(nn) :: sigma              ! market shares
    integer, dimension(nw+1) :: s                 ! market structure
    integer, dimension(nn) :: w                   ! quality vector
    real(wp) :: norm
    integer :: iter, i, err

    ! Decode the integer configuration and convert to a state vector
    s = lookup(:,is)
    w = to_w(s)

    ! Use Newton-Raphson to find the Nash equilibrium prices at this w
    ! by solving for the root of the the first order conditions.
    iter = 0
    norm = nr_tol + 1.0_wp
    p = 5.5_wp

    ! repeat until the convergence criteria are met
    do while ((norm > nr_tol) .and. (iter < maxiter))
       iter = iter + 1

       ! evalute the foc at p
       call model_profit_foc(p, w, y, sigma)

       ! calculate a finite difference Jacobian at p
       p_eps = p
       do i = 1, nn
          p_eps(i) = p(i) + eps
          call model_profit_foc(p_eps, w, y_eps, sigma)
          J(:,i) = (y_eps - y) / eps
          p_eps(i) = p(i)
       end do

       ! determine the descent direction
       y = -1.0_wp * y
       delta = solve_linear(J, y, err)
       p_new = p + 0.6_wp * delta

       ! calculate the relative norm and store the new prices
       norm = maxval(abs(p_new - p)) / maxval(abs(p))
       p = p_new

       ! Restart if negative prices are found
       if (any(p < 0.0_wp)) then
          p = real(iter)
          norm = nr_tol + 1.0_wp
       end if
    end do

    if (norm > nr_tol) then
       print *, 'ERROR: Newton-Raphson failed!'
       stop
    end if

    ! Calculate profits
    call model_profit_foc(p, w, y, sigma)
    pi_w = mktsize * (p - mc) * sigma

    ! Return the profits by quality level (excluding nw+1)
    pi_s = 0.0_wp
    do i = 1, nn
       if (w(i) < nw+1) then
          pi_s(w(i)) = pi_w(i)
       end if
    end do
  end function model_profit


  ! Implements the first-order conditions of the profit maximization
  ! problem in the differentiated products demand model.  See equation
  ! (4) of Pakes and McGuire (1994).  The inputs are p, denoting the
  ! vector of prices charged by each firm, and w denoting the market
  ! structure.  The outputs are y, the value of the vector of first
  ! order conditions, and sigma, the vector of market shares.
  subroutine model_profit_foc(p, w, y, sigma)
    real(wp), dimension(nn), intent(in) :: p
    integer, dimension(nn), intent(in) :: w
    real(wp), dimension(nn), intent(out) :: y, sigma
    real(wp), dimension(nn) :: expval

    where (w /= nw+1)
       where (w < wstar)
          expval = exp(w - p)
       elsewhere
          expval = exp(wstar - p) * (2.0_wp - exp(wstar - w))
       end where
    elsewhere
       expval = 0.0_wp
    end where
    sigma = expval / (1.0_wp + sum(expval))
    y = -(p - mc) * (1.0_wp - sigma) + 1.0_wp
  end subroutine model_profit_foc


  !------------------------------------!
  ! Data Generating Process            !
  !------------------------------------!

  ! Draw a player-specific state WI weighted by move arrival rates. The encoded
  ! state is given in IS.  The index of the chosen player among the s(WI)
  ! players in state WI is optionally returned in INDEX.
  function model_draw_position(rng, is, index) result(wi)
    type(osl_rng_t), intent(inout) :: rng
    integer, intent(in) :: is
    integer, intent(out), optional :: index
    integer, dimension(nw+1) :: s
    integer :: wi, i
    real(wp) :: u, cumulative_rate, total_rate
    real(wp), dimension(nw+1) :: rates

    s = lookup(:,is)

    ! Calculate rates for each quality level
    rates = 0.0_wp
    do i = 1, nw
       if (s(i) > 0) then
          if (i <= we) then
             rates(i) = real(s(i), wp) * lambda_low
          else
             rates(i) = real(s(i), wp) * lambda_high
          end if
       end if
    end do

    ! Potential entrants
    if (s(nw+1) > 0) then
       rates(nw+1) = lambda_low  ! Single potential entrant at low rate
    end if

    total_rate = sum(rates)

    if (total_rate > 0.0_wp) then
       u = osl_rng_uniform(rng) * total_rate
       cumulative_rate = 0.0_wp
       wi = 1

       do i = 1, nw+1
          cumulative_rate = cumulative_rate + rates(i)
          if (u <= cumulative_rate .and. s(i) > 0) then
             wi = i
             exit
          end if
       end do

       ! If we selected this quality level, now pick which firm within it
       if (present(index) .and. s(wi) > 0) then
          index = ceiling(real(s(wi), wp) * osl_rng_uniform(rng))
       end if
    else
       ! Fallback if no active players
       wi = nw + 1
       if (present(index)) index = 1
    end if
  end function model_draw_position

  ! Simulate the model until the next event.
  subroutine simulate_forward(rng, is_lag, wi_lag, tau, is, wi, a, index)
    type(osl_rng_t), intent(inout) :: rng
    integer, intent(in) :: is_lag
    integer, intent(out) :: wi_lag
    real(wp), intent(out) :: tau
    integer, intent(out) :: is
    integer, intent(out) :: wi
    integer, intent(out) :: a
    integer, intent(out), optional :: index

    real(wp) :: tau_s, tau_a, u, total_move_rate
    integer :: n_active, k, n_low, n_high

    ! Determine the number of active players by quality level
    call firm_count(is_lag, n_low, n_high)
    n_active = n_low + n_high

    ! Calculate total move arrival rate (competing risks across all players)
    total_move_rate = real(n_low, wp) * lambda_low + real(n_high, wp) * lambda_high

    ! Draw competing times: nature vs. any player move
    tau_s = osl_exponential_rnd(rng, gamma)
    tau_a = osl_exponential_rnd(rng, total_move_rate)

    ! Which happens next, a move or a jump?
    if (tau_s < tau_a) then
       ! Jump -----------------------------------------------------------
       wi_lag = WI_NATURE                   ! nature moved
       wi = WI_NATURE                       ! nature doesn't transition
       a = 0                                ! no move
       is = is_dep(is_lag)                  ! resulting state
       tau = tau_s                          ! store tau

       if (present(index)) then
          index = 0
       end if
    else
       ! Move -----------------------------------------------------------
       tau = tau_a                             ! store tau

       ! Draw which player moves weighted by their arrival rates
       wi_lag = model_draw_position(rng, is_lag, index)

       ! Draw a Uniform(0,1) random variable for simulating choices
       u = osl_rng_uniform(rng)

       if (wi_lag == WI_ENTRANT) then
          ! A potential entrant is moving

          ! Determine the optimal choice given that the state is IS_LAG
          if (u < ccpe(2, is_lag)) then
             ! Entry is optimal
             a = 2
             is = is_ent(is_lag)
             wi = we
          else
             ! Entry is not optimal
             a = 1
             is = is_lag
             wi = nw + 1
          end if
       else
          ! An incumbent is moving
          k = enc(is_lag, wi_lag)

          ! Determine the optimal choice given the state
          if (u < ccpi(1, k)) then
             ! remaining active is optimal
             a = 1
             is = is_lag
             wi = wi_lag
          else if (u < ccpi(1, k) + ccpi(2, k)) then
             ! investing is optimal
             a = 2
             is = is_inv(is_lag, wi_lag)
             wi = min(wi_lag + 1, nw)
          else
             ! exiting is optimal
             a = 3
             is = is_exit(is_lag, wi_lag)
             wi = WI_ENTRANT
          end if
       end if
    end if
  end subroutine simulate_forward

  ! Generate a continuous-time dataset with NM markets each with NT
  ! obervations using the parameters theta.
  subroutine model_generate_data_continuous(nm, nt)
    integer, intent(in) :: nm, nt
    integer :: im, it, a, n
    real(wp) :: t, tau
    integer :: is, is_lag, wi, wi_lag
    type(observation_t) :: obs

    ! For each market
    do im = 1, nm
       t   = 0.0_wp

       ! Generate random initial conditions.  First draw a state and
       ! determine the corresponding actions.
       is_lag = ceiling(ns * osl_rng_uniform(rng_data))

       ! Generate nt observations
       it = 1
       do while (it <= nt)

          ! Determine the current position in the dataset
          n = (im-1)*nt + it

          ! Simulate the next event
          call simulate_forward(rng_data, is_lag, wi_lag, tau, is, wi, a)

          ! Store the observation
          obs%m = im
          obs%t = t
          obs%tau = tau
          obs%is = is_lag
          obs%wi = wi_lag
          obs%a = a
          obs%isp = is
          obs%wip = wi

          ! Increment the time counter
          it = it + 1

          ! Pointer to next observation
          call dataset_append(dat, im, obs)

          ! Prepare for another loop
          t = t + tau
          is_lag = is
          wi_lag = wi
       end do
    end do
  end subroutine model_generate_data_continuous

  ! Using the current parameters, generate a dataset of NM
  ! markets each with NT observations.
  subroutine model_generate_data_discrete(nm, nt)
    integer, intent(in) :: nm, nt
    integer :: im, it, is, m_eps, i, imax
    real(wp) :: u, psum
    integer :: is_d, isp_d
    type(observation_t) :: obs
    real(wp), dimension(ns) :: basis
    real(wp), dimension(ns) :: P_row
    type(observation_t), pointer :: obs_ptr
    integer, dimension(:), allocatable :: P_tmp
    logical :: found

    ! Precalculate truncation term for expvunif
    m_eps = osl_poisson_inv(1.0_wp - Q_tol, max_q_ii)

    ! For each market
    do im = 1, nm
       ! Generate random initial conditions.  First draw a state and
       ! determine the corresponding actions.
       is_d = ceiling(ns * osl_rng_uniform(rng_data))

       ! Generate nt observations
       do it = 1, nt
          ! Basis vector for row is_d of matrix exponential
          basis = 0.0_wp
          basis(is_d) = 1.0_wp

          ! Apply uniformization algorithm to calculate is_d-th row of P
          P_row = expvunif_sparse(basis, DELTA, max_q_ii, QCSCp, QCSCi, QCSCx, M=m_eps)

          ! Take a draw from the discrete distribution given by P_row
          u = osl_rng_uniform(rng_data)
          psum = 0.0_wp
          isp_d = -1
          do is = 1, ns
             psum = psum + P_row(is)
             if (u < psum) then
                isp_d = is
                exit
             end if
          end do

          ! Sanity checks
          if (is_d < 1 .or. is_d > ns) then
             stop 'ERROR: invalid value of is_d in model_generate_data_discrete'
          end if
          if (isp_d < 1 .or. isp_d > ns) then
             stop 'ERROR: invalid value of isp_d in model_generate_data_discrete'
          end if

          ! Store transition from is_d to isp_d
          obs%m = im
          obs%t = it * DELTA
          obs%is = is_d
          obs%isp = isp_d
          obs%tau = DELTA
          obs%a = A_DISCRETE
          obs%wi = WI_DISCRETE
          obs%wip = WI_DISCRETE

          ! Append observation
          call dataset_append(dat, im, obs)

          ! Prepare for another loop
          is_d = isp_d
       end do
    end do

    ! Precalculate P_state vector indicating uniquely observed states
    if (allocated(P_state)) then
       deallocate(P_state)
    end if
    allocate(P_state(min(ns, nt*nm)))

    ! Walk dataset and construct list of all visited states (is) observed
    P_state = 0
    do im = 1, dataset_n(dat)
       obs_ptr => dataset_first(dat, im)
       do while (associated(obs_ptr))
          ! Determine whether this state has been seen before
          found = .false.
          do i = 1, size(P_state)
             if (P_state(i) == 0) then
                ! This is the last previously stored row, stop searching.
                found = .false.
                exit
             else if (P_state(i) == obs_ptr%is) then
                found = .true.
                exit
             end if
          end do

          ! If this is a new state, then add it to the list.
          if (.not. found) then
             if (obs_ptr%is < 1 .or. obs_ptr%is > ns) then
                stop 'ERROR: invalid value of obs_ptr%is in model_generate_data_discrete'
             end if
             P_state(i) = obs_ptr%is
          end if

          obs_ptr => dataset_next(obs_ptr)
       end do
    end do

    ! Find the last calculated row (which precedes location of first
    ! zero, if there are any zeros)
    imax = minloc(P_state, DIM=1)
    if (P_state(imax) == 0) then
       imax = imax - 1
    else
       imax = size(P_state)
    end if

    ! Reallocate P_state with final size
    allocate(P_tmp(imax))
    P_tmp = P_state(1:imax)
    deallocate(P_state)
    allocate(P_state(imax))
    P_state = P_tmp
    deallocate(P_tmp)
  end subroutine model_generate_data_discrete

  ! Using the current parameters (theta), generate a dataset of NM
  ! observations.  Depending on the settings, either discrete-time or
  ! continuous-time data will be generated.  Calling this procedure
  ! triggers a model reset, with a full recalculation of the value
  ! function.
  subroutine model_generate_data(theta, nm, nt, seed)
    real(wp), dimension(np), intent(in) :: theta
    integer, intent(in) :: nm
    integer, intent(in) :: nt
    integer, intent(in) :: seed

    ! Reset the model using the true parameters and given seed
    call model_reset(theta)
    ! Seed the data RNG
    call osl_rng_seed(rng_data, seed)

    ! Initialize storage for the dataset
    call dataset_init(dat, nm)

    ! Call the appropriate data generating procedure
    if (DISCRETE_TIME) then
       ! print '(a,i0,a)', 'Generating discrete-time data (seed = ', seed, ')...'
       call model_generate_data_discrete(nm, nt)
       call model_dataset_count()
    else
       ! print '(a,i0,a)', 'Generating continuous-time data (seed = ', seed, ')...'
       call model_generate_data_continuous(nm, nt)
    end if
  end subroutine model_generate_data


  ! Calculates summary statistics of the current dataset.
  subroutine model_dataset_stats()
    integer :: im, a, nobs
    real(wp) :: t, tau, tau_s, tau_a
    integer :: is, wi_lag
    integer :: n_jump
    integer, dimension(na_inc) :: n_inc
    integer, dimension(na_ent) :: n_ent
    integer :: n_active, w_sum, w_denom, w_max, w_min
    real(wp) :: w_avg, n_avg
    integer, dimension(nw+1) :: s
    integer, dimension(nn) :: w
    type(observation_t), pointer :: obs

    call osl_print_subheading("Sample Summary Statistics")

    tau_s = 0.0_wp
    tau_a = 0.0_wp
    n_jump = 0
    n_inc = 0
    n_ent = 0
    w_avg = 0
    n_avg = 0
    w_denom = 0
    w_max = 0
    w_min = nw+1

    do im = 1, dataset_n(dat)
       t = 0.0_wp
       obs => dataset_first(dat, im)
       do while (associated(obs))
          tau = obs%tau
          wi_lag = obs%wi
          a = obs%a
          is = obs%isp

          s = lookup(:,is)
          w = to_w(s)

          n_active = sum(s(1:nw))
          n_avg = n_avg + real(n_active, wp)
          w_sum = sum(w, MASK=w < nw+1)
          w_denom = w_denom + n_active
          w_avg = w_avg + w_sum
          w_max = max(w_max, maxval(w, MASK=w < nw+1))
          w_min = min(w_min, minval(w, MASK=w < nw+1))

          ! Continuous time data only
          if (.not. DISCRETE_TIME) then
             if (wi_lag == WI_NATURE) then
                ! JUMP -----------------------------------------------------------
                tau_s = tau_s + tau
                n_jump = n_jump + 1
             else
                ! MOVE -----------------------------------------------------------
                tau_a = tau_a + tau
                if (wi_lag == WI_ENTRANT) then
                   ! POTENTIAL ENTRANT -------------------------------------------
                   n_ent(a) = n_ent(a) + 1
                else
                   ! INCUMBENT ---------------------------------------------------
                   n_inc(a) = n_inc(a) + 1
                end if
             end if
          end if

          ! Move to the next observation
          obs => dataset_next(obs)
       end do
    end do

    nobs = dataset_size(dat)
    n_avg = n_avg / real(nobs, wp)
    w_avg = w_avg / real(w_denom, wp)

    call osl_print('Number of observations: ', nobs)
    call osl_print('Avg. num. active firms: ', n_avg)
    call osl_print('Avg. quality:           ', w_avg)
    call osl_print('Min. quality:           ', w_min)
    call osl_print('Max. quality:           ', w_max)
    if (.not. DISCRETE_TIME) then
       call osl_print('Incumbent actions (continue, invest, exit):', n_inc)
       call osl_print('Entrant actions (continue, enter):', n_ent)
       call osl_print('Market depreciation:    ', n_jump)
       call osl_print('Total events:           ', n_jump + sum(n_inc) + sum(n_ent))
    end if
  end subroutine model_dataset_stats

  ! Stores transition counts for discrete time data.
  subroutine model_dataset_count()
    integer :: im, nobs, is, isp
    type(observation_t), pointer :: obs
    type(spmatrix_t) :: counts_coo, temp

    ! This routine works for discrete time data only.  For CT data,
    ! the inter-event times are also important.
    if (.not. DISCRETE_TIME) then
       return
    end if

    ! Allocate sparse matrix in COO format for counting
    nobs = dataset_size(dat)
    call spmatrix_init(counts_coo, ns, ns, nobs, SPMATRIX_TYPE_COO)

    ! Walk dataset and add indicators for state transitions.
    ! Duplicate entries will be summed later, when the matrix is
    ! converted to canonical form.
    do im = 1, dataset_n(dat)
       obs => dataset_first(dat, im)
       do while (associated(obs))
          is = obs%is
          isp = obs%isp
          call spmatrix_set(counts_coo, is, isp, 1.0_wp)
          obs => dataset_next(obs)
       end do
    end do

    ! Convert matrix from COO to CSR to sum duplicates.
    ! Then convert back to COO for easy processing.
    call spmatrix_coo_to_csr(counts_coo, temp)
    call spmatrix_csr_to_coo(temp, dat_counts)
  end subroutine model_dataset_count

  !------------------------------------!
  ! State space                        !
  !------------------------------------!

  ! Converts an s vector, a length (nw+1) vector representing a market
  ! configuration, to a w vector, a vector of length n where each
  ! element gives the state of firm n.  Since we have assumed
  ! anonymity, this is not a one to one mapping.  To resolve the
  ! ambiguity, the returned w vector will be sorted in increasing
  ! order.
  function to_w(s) result(w)
    integer, dimension(nw+1), intent(in) :: s
    integer, dimension(nn) :: w
    integer :: i, off

    off = 0
    do i = 1, nw+1
       w(off+1:off+s(i)) = i
       off = off + s(i)
    end do
  end function to_w

  ! Given a market structure vector, return the total number of active firms.
  ! Includes all incumbents and at most one potential entrant.
  subroutine firm_count(is, n_low, n_high)
    integer, intent(in) :: is
    integer, intent(out) :: n_low, n_high

    n_low = sum(lookup(1:we, is))                 ! Incumbent firms (low)
    n_high = sum(lookup(we+1:nw, is))             ! Incumbent firms (high)
    if (lookup(nw+1, is) > 0) then
       n_low = n_low + 1                          ! Single potential entrant (low)
    end if
   end subroutine firm_count

  !------------------------------------!
  ! Support routines                   !
  !------------------------------------!

  ! Compute the solution to a real linear system of equations
  ! A * X = B where A is an arbitrary N-by-N matrix and B is a
  ! right-hand side vector of length N.  The solution is X a
  ! vector of length N.  The system is solved using the
  ! LU-decomposition via a call LAPACK's DGESV routine.
  !
  ! A nonzero return code ERR indicates an error.  A negative
  ! value ERR = -i indicates that the i-th argument to DGESV
  ! had an illegal value.  A positive value ERR = i indicates
  ! that U(i,i) is exactly zero and, although the decomposition
  ! was completed the U factor is singular and so the system
  ! could not be solved.
  function solve_linear(A, b, err) result(x)
    real(wp), dimension(:,:),       intent(in)  :: A
    real(wp), dimension(size(A,1)), intent(in)  :: b
    real(wp), dimension(size(A,1)) :: x
    integer, intent(out), optional :: err

    real(wp), dimension(size(A,1),size(A,2)) :: A_
    integer, dimension(size(A,1)) :: pivot
    integer :: n, info
    external DGESV

    if (size(A,1) /= size(A,2)) then
       stop 'solve_linear: matrix must be square!'
    end if

    n = size(A,1)
    x = b
    A_ = A
    call DGESV(n, 1, A_, n, pivot, x, n, info)
    if (present(err)) err = info
  end function solve_linear

  !------------------------------------!
  ! Utility Functions                  !
  !------------------------------------!

  ! Return the value of a bounded version of `log` designed to
  ! prevent underflow.
  elemental function blog(x)
    real(wp), intent(in) :: x
    real(wp) :: blog

    if (x < 1.0e-307_wp) then
       blog = log(1.0e-307_wp)
    else
       blog = log(x)
    end if
  end function blog

  ! Given a vector V of length N, this function returns the value
  ! LOG(EXP(V(1)) + ... + EXP(V(N))).  In order to prevent numerical
  ! overflow, it first shifts the vector V by a constant C,
  ! applies EXP and LOG, and then adjusts the result to account
  ! for the shift.
  !
  ! Consider the case when V = (A, B). Note that we can write
  ! EXP(A) + EXP(B) = [ EXP(A-C) + EXP(B-C) ]*EXP(C) for any C.
  ! Furthermore, we have
  ! LOG( ( EXP(A-C) + EXP(B-C) ) * EXP(C) ) = LOG( EXP(A-C) + EXP(B-C) ) + C.
  ! By choosing C appropriately, we can attempt to ensure that EXP is not
  ! evaluated at a number large enough to cause overflow.  Underflow is also
  ! possible as taking LOG of a number close to zero may result in -Inf.
  ! Thus, this function also adjusts for large negative elements of V.  Note that
  ! this function is still not completely robust to very poorly scaled vectors V.
  ! Overflow or underflow is still possible.
  !
  ! For a set of N deterministic values U_1, ..., U_N, if E_1, ..., E_N are iid
  ! and type 1 extreme value distributed, this function returns the expected
  ! value of MAX(U_1 + E_1, ..., U_N + E_N).
  function log_sum_exp(v) result(e)
    real(wp), dimension(:), intent(in) :: v
    real(wp) :: c, e

    ! Choose C to be the element of V that is largest in absolute value.
    if ( maxval(abs(v)) > maxval(v) ) then
       c = minval(v)
    else
       c = maxval(v)
    end if
    e = log(sum(exp(v-c))) + c
  end function log_sum_exp

  ! Start timer
  subroutine tic()
    call system_clock(count0)
    call cpu_time(time0)
  end subroutine tic

  ! Stop timer and report elapsed wall clock and CPU time
  subroutine toc()
    call system_clock(count1, count_rate)
    call cpu_time(time1)
    time = real(count1 - count0) / real(count_rate)
    print '(a,g12.5,a,g12.5,a)', 'Elapsed time: ', &
         time, ' sec. (', time1 - time0, ' sec. cpu time)'
  end subroutine toc

  ! Print a particular state in human readable form.
  subroutine print_k(k)
    integer, intent(in) :: k
    integer :: is, wi
    is = dec(k, k_is)
    wi = dec(k, k_wi)
    print '(a10,i0,a4,i0)', 'State k = ', k, ' of ', nk
    print '(a13,i0)', 'Quality wi = ', wi
    call print_state_int(is)
  end subroutine print_k

  ! Print a particular state in human readable form.
  subroutine print_state_int(is)
    integer, intent(in) :: is
    print '(a11,i0,a4,i0)', 'State is = ', is, ' of ', ns
    call print_state_vec(lookup(:, is))
  end subroutine print_state_int

  ! Print a particular state in human readable form.
  subroutine print_state_vec(s)
    integer, dimension(:), intent(in) :: s
    integer :: wi

    do wi = 1, nw
       print '(a13,i2,2x,i4)', 'Incumbent, wi = ', wi, s(wi)
    end do
    print '(a15,2x,i4)', 'Entrant',  s(WI_ENTRANT)
  end subroutine print_state_vec

end module model
