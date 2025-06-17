! rust_model.f90 -- A single agent renewal model
!
! Copyright (C) 2008-2025 Jason R. Blevins
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

module rust_model
  use osl, wp => osl_wp
  use dataset
  use rust_data
  use matrix
  implicit none

  private

  ! Public interface
  public :: rust_model_dgp
  public :: rust_model_load_data
  public :: rust_model_ll
  public :: rust_model_compare_hazards
  public :: rust_model_dataset_stats
  public :: rust_model_init
  public :: rust_model_free
  public :: rust_model_print
  public :: rust_model_get_ll_scaling
  public :: np_restr, np_full, np_abbe
  public :: LL_PENALTY

  public :: RUST_DATA_CT, RUST_DATA_DT
  public :: rust_model_t, rust_model_ptr_t

  ! Model variants
  public :: RUST_MODEL_HOMOGENEOUS, RUST_MODEL_HETEROGENEOUS, RUST_MODEL_ABBE

  !-----------------------------------------------!
  ! Interfaces                                    !
  !-----------------------------------------------!

  interface divide_maybe
     module procedure divide_maybe_vwp
     module procedure divide_maybe_vi
  end interface divide_maybe

  !-----------------------------------------------!
  ! Global parameters                             !
  !-----------------------------------------------!

  ! Dimensions and model constants
  real(wp), parameter :: rho          = 0.05_wp   ! Discount factor
  integer, parameter :: na            = 2         ! Number of choices
  integer, parameter :: nx            = 90        ! Number of states

  ! Solution method options
  integer, parameter :: nvf = 10                  ! Value function history

  ! Sampling options
  integer, parameter :: RUST_DATA_CT  = 1         ! Observed continuously
  integer, parameter :: RUST_DATA_DT  = 2         ! Observed at fixed intervals

  ! Model variant options
  integer, parameter :: RUST_MODEL_HOMOGENEOUS   = 1  ! Single lambda for all states
  integer, parameter :: RUST_MODEL_HETEROGENEOUS = 2  ! Two lambdas: lambda1, lambda2
  integer, parameter :: RUST_MODEL_ABBE          = 3  ! Fixed lambdas: lambda1 = lambda2 = 1.0

  ! Miscellaneous constants
  integer,  parameter :: MAX_IT       = 10000     ! Maximum no. iterations
  integer,  parameter :: NX_SMALL     = 10        ! State space limit for "small" problem
  real(wp), parameter :: X_SCALE      = 1.0_wp / nx ! Scale X to [ 0, 1 ]
  real(wp), parameter :: VF_TOL       = 1.0e-16_wp ! Value function tolerance
  real(wp), parameter :: LL_PENALTY   = -1.0e99_wp! Penalty value

  ! Structural parameter indices (homogeneous model)
  integer, parameter :: IDX_LAM    = 1
  integer, parameter :: IDX_Q1     = 2
  integer, parameter :: IDX_BETA   = 3
  integer, parameter :: IDX_COST   = 4
  integer, parameter :: np_restr   = 4

  ! Structural parameter indices (heterogeneous model)
  integer, parameter :: IDX_LAM1   = 1
  integer, parameter :: IDX_LAM2   = 2
  integer, parameter :: IDX_Q1_H   = 3
  integer, parameter :: IDX_BETA_H = 4
  integer, parameter :: IDX_COST_H = 5
  integer, parameter :: np_full    = 5

  ! Structural parameter indices (ABBE model - fixed lambdas)
  integer, parameter :: IDX_Q1_A   = 1
  integer, parameter :: IDX_BETA_A = 2
  integer, parameter :: IDX_COST_A = 3
  integer, parameter :: np_abbe    = 3

  !-----------------------------------------------!
  ! Structured data                               !
  !-----------------------------------------------!

  ! Structural parameters
  type :: param_t
     real(wp) :: lambda1  ! Heterogeneous arrival rate (low states)
     real(wp) :: lambda2  ! Heterogeneous arrival rate (high states)
     real(wp) :: q1
     real(wp) :: beta
     real(wp) :: cost
  end type param_t

  !-----------------------------------------------!
  ! Current state of the model                    !
  !-----------------------------------------------!

  type :: rust_model_t
     ! Model settings and defaults
     integer :: model_variant = RUST_MODEL_HETEROGENEOUS ! Model variant
     integer :: np = np_full                      ! Number of parameters
     integer :: data_type = RUST_DATA_CT          ! Sample type
     real(wp) :: delta = 1.0_wp                   ! Observation time interval
     real(wp) :: max_t = 100.0_wp                 ! Observation window
     integer :: nm = 1                            ! Number of "markets"
     integer :: nx_low                            ! Cutoff for heterogeneous model

     ! Debugging
     logical :: debug = .false.

     ! Status flags
     logical :: vf_error = .false.                ! VF calculation worked?

     ! Dataset
     type(dataset_t), allocatable :: data         ! Dataset storage
     logical :: mc_data = .false.                 ! Monte Carlo dataset?

     ! State variables
     real(wp), dimension(nx) :: x

     ! Arrival rates
     real(wp), dimension(nx) :: lambda

     ! Value function iteration
     real(wp), dimension(nx) :: v_prev            ! Previously calculated v.f.

     ! Previous value function storgae
     real(wp), dimension(:,:), allocatable :: theta_old
     real(wp), dimension(:,:), allocatable :: v_old
     integer :: ivf

     ! Parameter bounds
     real(wp), dimension(:), allocatable :: theta_lb, theta_ub

     ! Log likelihood scaling factor
     real(wp) :: ll_scale
  end type rust_model_t

  ! Pointer container for wrapper functions
  type :: rust_model_ptr_t
     type(rust_model_t), pointer :: p => null()
  end type rust_model_ptr_t

contains

  !---------------------------------------------------------------------------!
  ! Module utility routines                                                   !
  !===========================================================================!

  ! Initialize dynamically allocated memory and finalize settings.
  subroutine rust_model_init(self, data_type, nm, max_t, delta, debug, model_variant)
    type(rust_model_t), intent(inout) :: self
    integer, intent(in), optional :: data_type
    integer, intent(in), optional :: nm
    real(wp), intent(in), optional :: max_t
    real(wp), intent(in), optional :: delta
    logical, intent(in), optional :: debug
    integer, intent(in), optional :: model_variant
    integer :: ix

    ! Check settings
    if (present(model_variant)) then
       self%model_variant = model_variant
    end if
    if (present(data_type)) then
       self%data_type = data_type
    end if
    if (present(nm)) then
       self%nm = nm
    end if
    if (present(max_t)) then
       self%max_t = max_t
    end if
    if (present(delta)) then
       self%delta = delta
    end if
    if (present(debug)) then
       self%debug = debug
    else
       self%debug = .false.
    end if

    ! Set up parameter dimensions
    if (self%model_variant == RUST_MODEL_HETEROGENEOUS) then
      self%np = np_full
    else if (self%model_variant == RUST_MODEL_HOMOGENEOUS) then
      self%np = np_restr
    else
      self%np = np_abbe
    end if

    ! Cutoff for low lambda states
    self%nx_low = nx / 2

    ! Set up state variables
    do ix = 1, nx
      self%x(ix) = X_SCALE * real(ix - 1, wp)
    end do

    ! Initialize state-specific arrival rates
    self%lambda = 1.0_wp

    ! Initialize previous value function
    self%v_prev = -1.0_wp

    ! Initialize value function history storage
    if (nvf > 0) then
       allocate(self%v_old(nx, nvf), self%theta_old(self%np, nvf))
       self%v_old = 0.0_wp
       self%theta_old = 0.0_wp
       self%ivf = 1
    end if

    ! Allocate memory for boundary vectors
    allocate(self%theta_lb(self%np))
    allocate(self%theta_ub(self%np))

    ! Set lower and upper bounds based on model variant
    if (self%model_variant == RUST_MODEL_HETEROGENEOUS) then
       self%theta_lb = (/ 1.0e-4_wp, 1.0e-4_wp, 1.0e-4_wp, -50.0_wp, -100.0_wp /)
       self%theta_ub = (/ 5.0_wp,    5.0_wp,    5.0_wp,      0.0_wp,   0.0_wp /)
    else if (self%model_variant == RUST_MODEL_HOMOGENEOUS) then
       self%theta_lb = (/ 1.0e-4_wp, 1.0e-4_wp, -50.0_wp, -100.0_wp /)
       self%theta_ub = (/ 5.0_wp,    5.0_wp,      0.0_wp,   0.0_wp /)
    else
       self%theta_lb = (/ 1.0e-4_wp, -50.0_wp, -100.0_wp /)
       self%theta_ub = (/ 5.0_wp,      0.0_wp,   0.0_wp /)
    end if

    ! Initialize internal dataset
    allocate(self%data)
    call dataset_init(self%data)
  end subroutine rust_model_init

  ! Free dynamically allocated memory.
  subroutine rust_model_free(self)
    type(rust_model_t), intent(inout) :: self

    if (nvf > 0) then
       deallocate(self%v_old)
       deallocate(self%theta_old)
       self%ivf = 1
    end if
    call dataset_free(self%data)
    deallocate(self%data)
    deallocate(self%theta_lb)
    deallocate(self%theta_ub)
  end subroutine rust_model_free


  ! Print parameter values with labels.
  subroutine rust_model_print(self, theta)
    type(rust_model_t), intent(inout) :: self
    real(wp), dimension(:), intent(in) :: theta
    real(wp), dimension(nx,nx) :: Q
    real(wp), dimension(nx) :: pr_a, hazard, period
    real(wp), dimension(nx) :: v
    type(param_t) :: param
    integer :: i

    print '(a)', 'Model Implications'
    print '(a)', '------------------'
    print *

    ! Store individual parameters
    call read_parameters(self, param, theta, self%model_variant)

    select case (self%data_type)
    case (RUST_DATA_CT)
       print '(a,i5,a,f15.2,a)', 'Data: ', self%nm, ' markets observed continuously on [0, ', self%max_t, ']'
    case (RUST_DATA_DT)
       if (self%mc_data) then
         print '(a,i5,a,f15.2,a,f8.2)', 'Data: ', self%nm, ' simulated markets observed discretely on [0, ', &
              self%max_t, '] with Delta =', self%delta
       else
         print '(a,i5,a)', 'Data: ', dataset_size(self%data), ' bus-month observations from Rust (1987) with Delta = 1.0'
       end if
    end select

    print '(a)', 'Value function: Value function iteration'
    print '(a,i5)', 'Number of states: ', nx
    print '(a,i5)', 'Number of actions: ', na
    print *
    print '(a)', 'Structural parameters:'
    print *
    if (self%model_variant == RUST_MODEL_HETEROGENEOUS) then
       print '(100a9)', 'lambda1', 'lambda2', 'q1', 'beta', 'cost'
       print '(100f9.4)', param%lambda1, param%lambda2, param%q1, param%beta, param%cost
    else  ! Homogeneous or ABBE
       print '(100a9)', 'lambda', 'q1', 'beta', 'cost'
       print '(100f9.4)', param%lambda1, param%q1, param%beta, param%cost
    end if
    print *

    ! Calculate the value function and implied hazards
    v = rust_model_vf_vfi(self, param)
    pr_a = rust_model_pr_a(param, v)

    if (nx <= NX_SMALL) then
       ! Report value function
       print '(a14,100f7.2)', 'Value function:', v
       print *

       ! Report replacement CCPs and hazards
       print '(a)', 'Replacement:'
       print '(a14,100f7.2)', 'Hazards:', self%lambda * pr_a
       print *
    end if

    ! Construct overall event hazards
    Q = rust_model_Q(param%q1, self%lambda * pr_a)
    do i = 1, nx
       hazard(i) = -Q(i,i)
    end do
    period = 1.0_wp
    period = divide_maybe(period, hazard - self%lambda * pr_a)

    if (nx <= NX_SMALL) then
       ! Report jump hazards and periods
       print '(a)', 'Jumps:'
       print '(a14,100f7.2)', 'Hazards:', hazard - self%lambda * pr_a
       print *

       ! Report overall event hazards and preiods
       print '(a)', 'Jumps + Replacement:'
       print '(a14,100f7.2)', 'Hazards:', hazard
       print *
    end if
  end subroutine rust_model_print


  !---------------------------------------------------------------------------!
  ! Model implementation                                                      !
  !===========================================================================!

  ! Return the aggregate intensity matrix constructed using the given rates.
  function rust_model_Q(q1, hazard) result(Q)
    real(wp), intent(in) :: q1                    ! One-state jump hazard
    real(wp), dimension(nx), intent(in) :: hazard ! Replacement hazards
    real(wp), dimension(nx,nx) :: Q
    real(wp) :: q_diag
    integer :: i

    ! Nature's intensity matrix
    Q = 0.0_wp
    q_diag = -q1
    do i = 1, nx - 1
       ! Leave state i at rate q1
       Q(i,i) = q_diag

       ! Jump to state i + 1 at rate q1
       Q(i,i+1) = q1
    end do

    ! Agent's intensity matrix
    do i = 1, nx
       ! Exit state i
       Q(i,i) = Q(i,i) - hazard(i)

       ! Arrive back in state 1
       Q(i,1) = Q(i,1) + hazard(i)
    end do
  end function rust_model_Q


  ! Return flow payoffs for each state.
  function rust_model_utility(self, param) result(u)
   type(rust_model_t), intent(in) :: self
    type(param_t), intent(in) :: param
    real(wp), dimension(nx) :: u
    u = param%beta * self%x
  end function rust_model_utility


  ! Solve for the value function using value function iteration.
  function rust_model_vf_vfi(self, param) result(v_new)
    type(rust_model_t), intent(inout) :: self
    type(param_t), intent(in) :: param            ! Structural parameters
    real(wp), dimension(nx, na) :: vm             ! Choice-specific values
    real(wp), dimension(nx) :: v_new              ! Resulting v.f.
    real(wp), dimension(nx) :: v_old              ! Previous v.f.
    real(wp), dimension(nx) :: rate               ! Next event rates
    real(wp), dimension(nx) :: Emax               ! Social surplus function
    real(wp), dimension(nx) :: utility            ! State-specific flow payoffs
    real(wp), dimension(nx) :: lambda             ! State-specific arrival rates
    real(wp) :: q1, beta, cost                    ! Individual parameters
    real(wp) :: norm                              ! Relative v.f. change
    integer :: iter                               ! Iteration counter
    integer :: ix                                 ! State index

    if (self%debug) then
       call osl_print_heading('Value function iteration')
    end if

    ! Store individual parameters for clarity
    lambda = self%lambda
    q1 = param%q1
    beta = param%beta
    cost = param%cost

    ! Calculate flow payoffs
    utility = rust_model_utility(self, param)

    ! Starting value
    if (self%vf_error) then
       v_new = 1.0_wp / (1.0_wp - rho) * utility
    else
       ! Retrieve value function from a nearby parameter vector
       call model_retrieve_vf(self, param)
       v_new = self%v_prev
    end if

    if (self%debug) then
       print '(i10,g17.5,100f8.2)', 0, 0, v_new
    end if

    ! Precalculate value function denominator
    rate(1:nx-1) = rho + lambda(1:nx-1) + q1
    rate(nx) = rho + lambda(nx)

    ! Iterate until relative sup norm is below VF_TOL or MAX_IT is exceeded.
    norm = huge(norm)
    iter = 0
    do while ((norm > VF_TOL) .and. (iter < MAX_IT))
       iter = iter + 1
       v_old = v_new

       ! Flow utility and mileage accumulation
       v_new = utility

       ! Mileage accumulation (except for state ix = nx)
       v_new(1:nx-1) = v_new(1:nx-1) + q1 * v_old(2:nx)

       ! Continuation values
       vm(:,1) = v_old                         ! Continue
       vm(:,2) = v_old(1) + cost               ! Reset
       do ix = 1, nx
          Emax(ix) = log_sum_exp(vm(ix,:)) + OSL_EULER
       end do
       v_new = v_new + lambda * Emax

       ! Denominator
       v_new = v_new / rate

       ! Calculate the relative sup norm of the difference
       norm = maxval(abs(v_new - v_old))
       if (maxval(abs(v_old)) > osl_eps_wp) then
          norm = norm / maxval(abs(v_old))
       end if

       ! Report iteration
       if (self%debug) then
          print '(i10,g17.5,100f8.2)', iter, norm, v_new
       end if

       ! Check for errors
       if (iter == MAX_IT) then
          self%vf_error = .true.
          return
       end if
    end do
    self%vf_error = .false.
    self%v_prev = v_new

    ! Store this value function and parameter vector
    call model_store_vf(self, param, v_new)

    if (self%debug) then
       print '(a27,100f8.2)', 'Value function: ', v_new
    end if
  end function rust_model_vf_vfi

  ! Store the current value function and parameter vector.
  subroutine model_store_vf(self, param, v)
    type(rust_model_t), intent(inout) :: self
    type(param_t), intent(in) :: param
    real(wp), dimension(nx), intent(in) :: v
    real(wp), dimension(self%np) :: theta

    if (nvf > 0) then
       if (self%model_variant == RUST_MODEL_HETEROGENEOUS) then
          theta = [ param%lambda1, param%lambda2, param%q1, param%beta, param%cost ]
       else if (self%model_variant == RUST_MODEL_HOMOGENEOUS) then
          theta = [ param%lambda1, param%q1, param%beta, param%cost ]
       else
          theta = [ param%q1, param%beta, param%cost ]
       end if
       self%theta_old(:,self%ivf) = theta
       self%v_old(:,self%ivf) = v
       self%ivf = self%ivf + 1
       if (self%ivf == nvf + 1) then
          self%ivf = 1
       end if
    end if
  end subroutine model_store_vf

  ! Retrieve the stored value function for the nearest parameter vector.
  ! Store it in self%v_prev.
  subroutine model_retrieve_vf(self, param)
    type(rust_model_t), intent(inout) :: self
    type(param_t), intent(in) :: param
    real(wp), dimension(self%np) :: theta
    real(wp), dimension(nvf) :: d
    integer :: i

    ! If value function storage is implemented
    if (nvf > 0) then
       if (self%model_variant == RUST_MODEL_HETEROGENEOUS) then
          theta = [ param%lambda1, param%lambda2, param%q1, param%beta, param%cost ]
       else if (self%model_variant == RUST_MODEL_HOMOGENEOUS) then
          theta = [ param%lambda1, param%q1, param%beta, param%cost ]
       else
          theta = [ param%q1, param%beta, param%cost ]
       end if
       do i = 1, nvf
          d(i) = osl_norm_2(self%theta_old(:,i) - theta)
       end do
       i = minloc(d, DIM=1)
       self%v_prev = self%v_old(:,i)
    end if
  end subroutine model_retrieve_vf

  ! Return choice probabilities calculated using the given value
  ! function V and the parameters given by PARAM.
  function rust_model_pr_a(param, v) result(pr_a)
    type(param_t), intent(in) :: param          ! Structural parameters
    real(wp), dimension(nx), intent(in) :: v    ! Value function
    real(wp), dimension(nx) :: pr_a             ! Choice probabilities

    pr_a = v - v(1) - param%cost
    pr_a = 1.0_wp / (1.0_wp + bexp(pr_a))
  end function rust_model_pr_a


  !----------------------------------------------------------------------------
  ! Data generating process
  !============================================================================

  ! Simulate the model until the next event using the rates given by
  ! Q1 and HAZARD and starting in the state X_LAG.  Random
  ! numbers are drawn using the provided RNG.  The next realized event
  ! is described by the return values TAU, PLAYER, A, and X.
  subroutine rust_model_simulate_event(q1, hazard, x_lag, rng, tau, player, a, x)
    real(wp), intent(in) :: q1
    real(wp), dimension(nx), intent(in) :: hazard
    integer, intent(in) :: x_lag
    type(osl_rng_t), intent(inout) :: rng
    real(wp), intent(out) :: tau
    integer, intent(out) :: player
    integer, intent(out) :: a
    integer, intent(out) :: x
    real(wp) :: q_diag, tau_x, tau_a
    logical, parameter :: debug = .false.

    ! Draw jump and move arrivals.
    if (x_lag < nx) then
       q_diag = q1
    else
       q_diag = 0.0_wp
    end if
    if (q_diag > 0.0_wp) then
       tau_x = osl_exponential_rnd(rng, q_diag)
    else
       tau_x = huge(tau_x)
    end if
    ! Note that here we simulate arrival times from the replacement
    ! hazard (as opposed to move arrivals at rate lambda followed by
    ! separate continuation and replacement decisions).
    if (hazard(x_lag) > 0.0_wp) then
       tau_a = osl_exponential_rnd(rng, hazard(x_lag))
    else
       tau_a = huge(tau_a)
    end if

    ! Which happens next, a move or a jump?
    if (tau_x < tau_a) then
       ! Jump -----------------------------------------------------------
       player = I_NATURE                       ! Nature moved
       a = A_NATURE                            ! No action
       tau = tau_x                             ! Store tau
       ! Jump forward one state, unless already at nx
       if (x_lag <= nx - 1) then
          x = x_lag + 1
       else
          x = x_lag
       end if
       if (debug) then
          print 2001, tau, 'Nature', 'x+1', x_lag, x
       end if
    else
       ! Move -----------------------------------------------------------
       tau    = tau_a                       ! Store move interval
       player = I_AGENT                     ! The agent moved
       a      = A_REPLACE                   ! Replace
       x      = 1                           ! Resets the state
       if (debug) then
          print 2001, tau, 'Agent', 'Replace', x_lag, x
       end if
    end if
2001 format ("tau = ", g12.5, "  i = ", a12, "  a = ", a12, "  x: ", i3, " -> ", i3)
  end subroutine rust_model_simulate_event


  ! Generate a dataset using parameters THETA with the given number of
  ! markets NM, each observed over an interval of length max_t.  The
  ! dataset will be stored in the internal dataset DATA.  Optionally,
  ! when FILE is present, write the dataset to the file associated
  ! unit FILE.  The actual number of observations in the dataset will
  ! be determined stochastically by the number of events that occur
  ! before time MAX_T.
  subroutine rust_model_dgp(self, theta, seed, file, debug)
    type(rust_model_t), intent(inout) :: self
    real(wp), dimension(:), intent(in) :: theta   ! Structural parameters
    integer, intent(in) :: seed                   ! Random number generator seed
    integer, optional, intent(in) :: file         ! Output file unit number
    logical, optional, intent(in) :: debug        ! Whether to print dataset statistics

    ! Parameters
    real(wp), parameter :: burn_t = 0.0_wp

    ! Local variables
    real(wp), dimension(nx) :: pr_a               ! Choice probabilities
    real(wp), dimension(nx) :: hazard             ! Replacement hazards
    real(wp), dimension(nx) :: v                  ! Value function
    type(osl_rng_t) :: rng                        ! Random number generator
    type(param_t) :: param                        ! Structural parameters
    real(wp) :: t, tau, dt                        ! Jump and move times
    integer :: x, x_lag                           ! State
    integer :: a, player                          ! Move, player
    integer :: obs                                ! Total number of observations
    integer :: i                                  ! Current market observation
    integer :: m                                  ! Current market
    integer :: dx, dx_lag                         ! Discrete time markers
    integer, dimension(nx) :: n_move              ! Number of underlying moves
    real(wp), dimension(nx) :: tau_a              ! Holding times for moves
    integer, dimension(nx) :: n_jump              ! Number of underlying jumps
    real(wp), dimension(nx) :: tau_x              ! Holding times for jumps

    ! Indicate we're generating Monte Carlo data (for reporting purposes)
    self%mc_data = .true.

    ! Read parameters
    call read_parameters(self, param, theta, self%model_variant)

    ! Obtain implied choice probabilities and hazards
    v = rust_model_vf_vfi(self, param)
    pr_a = rust_model_pr_a(param, v)
    hazard = self%lambda * pr_a

    ! Seed the random number generator
    call osl_rng_seed(rng, seed)

    ! Print a header
    if (present(file)) write(file, 1000) 'market', 'time', 'tau', 'player', 'x_lag', 'x', 'a'

    ! Dataset Initializations
    n_move = 0
    n_jump = 0
    obs = 0
    tau_a = 0.0_wp
    tau_x = 0.0_wp

    ! Simulate each market over the interval [0, max_t]
    do m = 1, self%nm

       ! Market Initializations
       t = -burn_t
       tau = 0.0_wp
       x_lag = ceiling(real(nx, wp) * osl_rng_uniform(rng))

       ! Discrete data case
       dt = 0
       dx_lag = x_lag
       dx = x_lag

       ! Simulate an unknown number of events over the interval [0, max_t]
       i = 0
       do
          ! Simulate a single event
          call rust_model_simulate_event(param%q1, hazard, &
               x_lag, rng, tau, player, a, x)

          ! Wait until the burn-in period has passed before recording observations
          if (t >= 0.0_wp) then

             ! Count state jumps and moves by the agent
             if (player == I_NATURE) then
                n_jump(x_lag) = n_jump(x_lag) + 1
                tau_x(x_lag) = tau_x(x_lag) + tau
             else
                n_move(x_lag) = n_move(x_lag) + 1
                tau_a(x_lag) = tau_a(x_lag) + tau
             end if

             ! Write the observation to the specified file handle if requested
             if (present(file)) then
                write(file, 1001) m, t, tau, player, x_lag, x, a
             end if

             ! We store data differently depending on the value of RUST_DATA
             if (self%data_type == RUST_DATA_CT) then

                ! Increment the observation counter
                i = i + 1

                ! In continuous time, stop when time exceeds MAX_T.
                if (t + tau < self%max_t) then
                   ! Store the observation
                   call dataset_append(self%data, M=m, T=(t+tau), TAU=tau, &
                        X=x_lag, I=player, A=a, XP=x)
                else
                   ! Store the final observation for this market
                   obs = obs + i
                   tau = self%max_t - t
                   call dataset_append(self%data, M=m, T=self%max_t, TAU=tau, &
                        X=x_lag, I=I_END, A=A_END, XP=x_lag)
                   exit
                end if

             else if (self%data_type == RUST_DATA_DT) then

                !   ^
                !   |                         *--o
                !   |       *---------o
                !   |                 *-------o
                !   |-------o         :          *------------
                !   |       :         :
                !   +----+----+----+----+----+----+----+---->
                !           t       t+tau

                ! Record potentially several observations since the last event
                do while (dt + self%DELTA < min(t + tau, self%max_t) - osl_eps_wp)
                   ! Update the discrete-time state marker
                   if (dt <= t) then
                      dx = x_lag
                   else if (dt >= t) then
                      dx_lag = x_lag
                   end if

                   ! Increment the observation and time counters
                   i = i + 1
                   dt = dt + self%DELTA

                   ! State stays at x between i and dt
                   call dataset_append(self%data, M=m, T=dt, X=dx_lag, XP=dx, &
                        TAU=self%DELTA, I=I_DISCRETE, A=A_DISCRETE)
                   dx_lag = dx

                end do

                ! Discrete time, stop when discrete counter reaches MAX_T / DELTA.
                if (t + tau > self%max_t - osl_eps_wp) then

                   i = i + 1
                   obs = obs + i
                   dt = dt + self%DELTA
                   call dataset_append(self%data, M=m, T=dt, TAU=self%DELTA, X=dx_lag, &
                        I=I_DISCRETE, A=A_END, XP=x)
                   exit

                end if

             end if

          end if

          ! Prepare to simulate the next event
          x_lag = x
          t = t + tau
       end do

    end do

    ! Report some summary statistics
    if (present(debug)) then
       if (debug) then
          print *
          print '(a)', 'Continuous Time Data Statistics'
          print '(a)', '-------------------------------'
          print *
          print '(a12,i10)', 'DGP seed:', seed
          print *
          print '(a12,3a10)', '', 'Jumps', 'Replace', 'Total'
          print '(a12,3i10)', 'Total', sum(n_jump), sum(n_move), sum(n_jump) + sum(n_move)
          print '(a12,3f10.3)', 'Fraction', &
               real(sum(n_jump), wp) / real(obs, wp), &
               real(sum(n_move), wp) / real(obs, wp), &
               real(sum(n_jump) + sum(n_move), wp) / real(obs, wp)
          print *
          print '(a)', 'Replacement:'
          print '(a14,1000f7.2)', 'Hazards:', divide_maybe(real(n_move, wp), tau_a + tau_x)
          print *
          print '(a)', 'Jumps:'
          print '(a14,1000f7.2)', 'Hazards:', divide_maybe(real(n_jump, wp), tau_a + tau_x)
          print *
          print '(a)', 'Jumps + Replacement:'
          print '(a14,1000f7.2)', 'Hazards:', divide_maybe(real(n_move + n_jump, wp), tau_a + tau_x)
          print *
          if (self%data_type == RUST_DATA_DT) then
             print '(a)', 'Number of events per observation period:'
             print '(a14,1000f7.2)', 'Replacements:', divide_maybe(self%DELTA * n_move, tau_a)
             print '(a14,1000f7.2)', 'Jumps:', divide_maybe(self%DELTA * n_jump, tau_x)
             print '(a14,1000f7.2)', 'All events:', divide_maybe(self%DELTA * (n_move + n_jump), tau_a + tau_x)
             print *
          end if
       end if
    end if

1000 format (7(a8,2x))                            ! Dataset header
1001 format (i8,2x,f8.2,2x,f8.2,2x,i8,2x,i8,2x,i8,2x,i8)  ! Observation format
  end subroutine rust_model_dgp

  subroutine rust_model_load_data(self, groups)
    type(rust_model_t), intent(inout) :: self
    integer, dimension(:) :: groups
    ! Indicate we're using real data (for reporting purposes)
    self%mc_data = .false.
    call osl_print_subheading("Loading Bus Data")
    call rust_data_load_groups(self%data, groups)
  end subroutine rust_model_load_data

  ! Print summary statistics for the current dataset.
  subroutine rust_model_dataset_stats(self)
    type(rust_model_t), intent(inout) :: self

    if (dataset_size(self%data) > 0) then
       if (self%data_type == RUST_DATA_DT) then
          call rust_model_dataset_stats_dt(self)
       else
          call rust_model_dataset_stats_ct(self)
       end if
    end if
  end subroutine rust_model_dataset_stats


  ! Print summary statistics for the current continuous time dataset.
  subroutine rust_model_dataset_stats_ct(self)
    type(rust_model_t), intent(inout) :: self
    real(wp), dimension(nx) :: tau_a
    real(wp), dimension(nx) :: tau_x, tau_all
    integer, dimension(nx) :: nmove
    integer, dimension(nx) :: njump
    type(observation_t), pointer :: obs
    real(wp) :: tau
    integer :: x_lag, x, a

    ! Initialize statistics
    njump = 0
    nmove = 0                                     ! Number of moves in x
    tau_x = 0.0_wp
    tau_a = 0.0_wp

    obs => dataset_first(self%data)
    do while (associated(obs))

       ! Read this observation
       call dataset_get(obs, TAU=tau, A=a, XP=x, X=x_lag)

       ! a = 0 denotes a jump, a = 1 or 2 denotes a move
       if (a == A_NATURE) then
          ! Jump: tally jump intervals and sum interval lengths
          njump(x_lag) = njump(x_lag) + 1
          tau_x(x_lag) = tau_x(x_lag) + tau
       else if (a > 0) then
          ! Move: count choices and move intervals by state
          nmove(x_lag) = nmove(x_lag) + 1
          tau_a(x_lag) = tau_a(x_lag) + tau
       end if

       obs => dataset_next(obs)
    end do

    ! Summary statistics
    tau_all = divide_maybe(tau_x + tau_a, real(njump + nmove, wp))

    ! Report statistics and compare with model predictions
    print *
    print '(a)', 'Observed Dataset Statistics'
    print '(a)', '---------------------------'
    print *
    print '(a14,a)', 'Sampling:', 'continuous time'
    print '(a14,i10)', 'N:', dataset_size(self%data)

    print *
    print '(a14,3a10)', '', 'Jumps', 'Moves', 'Events'
    print '(a14,3i10)', 'Total', sum(njump), sum(nmove), &
         sum(njump) + sum(nmove)
    print '(a12,3f10.3)', 'Fraction', &
         real(sum(njump)) / real(dataset_size(self%data)), &
         real(sum(nmove)) / real(dataset_size(self%data)), &
         real(sum(njump) + sum(nmove)) / real(dataset_size(self%data))
    print *
    print '(a14,100i10)', 'Jumps', njump
    print '(a14,100i10)', 'Moves', nmove
    print *
    print '(a14,100f7.2)', 'Mean tau', tau_all
    print *
    print '(a)', 'Replacement:'
    print '(a14,1000f7.2)', 'Hazards:', divide_maybe(real(nmove, wp), tau_a + tau_x)
    print '(a14,1000f7.2)', 'Periods:', divide_maybe(tau_a + tau_x, real(nmove, wp))
    print *
    print '(a)', 'Jumps:'
    print '(a14,1000f7.2)', 'Hazards:', divide_maybe(real(njump, wp), tau_a + tau_x)
    print '(a14,1000f7.2)', 'Periods:', divide_maybe(tau_a + tau_x, real(njump, wp))
    print *
    print '(a)', 'Jumps + Replacement:'
    print '(a14,1000f7.2)', 'Hazards:', divide_maybe(real(nmove + njump, wp), tau_a + tau_x)
    print '(a14,1000f7.2)', 'Periods:', divide_maybe(tau_a + tau_x, real(nmove + njump, wp))
    print *
  end subroutine rust_model_dataset_stats_ct


  ! Print summary statistics for the current discrete time dataset.
  subroutine rust_model_dataset_stats_dt(self)
    type(rust_model_t), intent(inout) :: self
    type(observation_t), pointer :: obs
    integer :: i, a, x, x_lag, nobs
    real(wp) :: tau, replace
    real(wp), dimension(nx, nx) :: P_hat

    ! Report basic sample information
    print *
    print '(a)', 'Observed Dataset Statistics'
    print '(a)', '---------------------------'
    print *
    print '(a14,a)', 'Sampling: ', 'discrete time'
    print '(a14,i10)', 'N:', dataset_size(self%data)
    print *

    ! Estimate the transition matrix nonparametrically
    P_hat = 0.0_wp
    replace = 0.0_wp
    nobs = 0
    obs => dataset_first(self%data)
    do while (associated(obs))
       ! Read this observation
       call dataset_get(obs, TAU=tau, I=i, A=a, XP=x, X=x_lag)
       nobs = nobs + 1

       ! Sanity checks
       if (i /= I_DISCRETE .and. i /= I_END) then
          call osl_print('ERROR: Unexpected player index in discrete time dataset, i = ', i)
          exit
       end if
       if (a /= A_DISCRETE .and. a /= A_END) then
          call osl_print('ERROR: Unexpected action index in discrete time dataset, a = ', a)
          exit
       end if
       if (abs(tau - self%DELTA) > osl_eps_wp) then
          call osl_print('ERROR: Unexpected observation interval in discrete time dataset, tau = ', tau)
          exit
       end if
       ! Increment transition count
       P_hat(x_lag, x) = P_hat(x_lag, x) + 1.0_wp
       ! Increment replacement count
       if (x == 1 .and. x_lag > 1) then
          replace = replace + 1.0_wp
       end if
       obs => dataset_next(obs)
    end do
    ! Estimate transition matrix
    do i = 1, nx
       if (sum(P_hat(i,:)) > 0.0_wp) then
           P_hat(i,:) = P_hat(i,:) / sum(P_hat(i,:))
       else
           P_hat(i,:) = 0.0_wp
       end if
    end do
    ! Replacement frequency
    replace = replace / real(nobs, wp)

    ! Only print full estimated transition matrix for small state spaces
    if (nx <= NX_SMALL) then
       call osl_print('Est. P(Delta) = ', P_hat, FMT='f7.2')
    else
       call osl_print('Est. P(Delta) (sample) = ', P_hat(1:NX_SMALL, 1:NX_SMALL), FMT='f7.2')
    end if
    print *
    call osl_print('Unconditional replacement frequency = ', replace)
  end subroutine rust_model_dataset_stats_dt


  !---------------------------------------------------------------------------!
  ! Log likelihood                                                            !
  !===========================================================================!

  ! Evaluate the log likelihood function for all parameters THETA.
  ! Obtains V via value function iteration, stores V and CCPs, and
  ! calls the main log-likelihood function.
  subroutine rust_model_ll(theta, ll, ctx)
    real(wp), dimension(:), intent(in) :: theta
    real(wp), intent(out) :: ll
    integer, dimension(:), intent(in), optional :: ctx
    type(rust_model_t), pointer :: self
    type(rust_model_ptr_t) :: model_ptr
    type(observation_t), pointer :: obs
    real(wp), dimension(nx,nx) :: Q, P
    real(wp), dimension(nx) :: v, h
    type(param_t) :: param
    integer :: a, x, x_lag
    real(wp) :: tau

    ! Obtain pointer to model
    if (present(ctx)) then
       model_ptr = transfer(ctx, model_ptr)
       self => model_ptr%p
    else
       stop 'A model pointer must be passed via the CTX argument'
    end if

    ! Store individual parameters
    if (any(theta < self%theta_lb) .or. any(theta > self%theta_ub)) then
       ll = -LL_PENALTY
       return
    else
       call read_parameters(self, param, theta, self%model_variant)
    end if

    ! Value function iteration
    v = rust_model_vf_vfi(self, param)

    if (self%vf_error) then
       ll = -LL_PENALTY
    else
       ! Calculate hazards
       h = self%lambda * rust_model_pr_a(param, v)

       ! Construct the intensity matrix
       Q = rust_model_Q(param%q1, h)

       ! Evaluate the log likelihood for discrete and continuous time separately
       ll = 0.0_wp
       obs => dataset_first(self%data)
       if (self%data_type == RUST_DATA_DT) then

          ! Calculate the aggregate P matrix time interval DELTA
          P = expm(self%DELTA, Q)

          do while (associated(obs))
             call dataset_get(obs, X=x_lag, XP=x)
             ll = ll + blog(P(x_lag, x))
             obs => dataset_next(obs)
          end do

       else

          do while (associated(obs))
             call dataset_get(obs, A=a, TAU=tau, X=x_lag, XP=x)
             if (a == A_END) then
                ll = ll + Q(x_lag, x_lag) * tau
             else
                ll = ll + blog(Q(x_lag, x)) + Q(x_lag, x_lag) * tau
             end if
             obs => dataset_next(obs)
          end do

       end if

       ! Return the rescaled negative log likelihood for minimization
       self%ll_scale = -1.0 / (self%max_t * real(self%data%size, wp))
       ll = self%ll_scale * ll
    end if

  end subroutine rust_model_ll

  ! Returns the scaling factor used for the log likelihood function
  function rust_model_get_ll_scaling(self) result(scale)
    type(rust_model_t), intent(in) :: self
    real(wp) :: scale
    scale = self%ll_scale
  end function rust_model_get_ll_scaling

  !---------------------------------------------------------------------------!
  ! Utility functions                                                         !
  !===========================================================================!

  subroutine rust_model_compare_hazards(self, theta0, theta)
    type(rust_model_t), intent(inout) :: self
    real(wp), dimension(:), intent(in) :: theta0
    real(wp), dimension(:), intent(in) :: theta
    real(wp), dimension(nx) :: model_ccp, model_hazard, model_v
    real(wp), dimension(nx) :: est_ccp, est_hazard, est_v
    type(param_t) :: param

    ! Calculate true model quantities for comparison
    call read_parameters(self, param, theta0, self%model_variant)
    model_v = rust_model_vf_vfi(self, param)
    model_ccp = rust_model_pr_a(param, model_v)
    model_hazard = self%lambda * model_ccp

    ! Calculate estimated model quantities for comparison
    call read_parameters(self, param, theta, self%model_variant)
    est_v = rust_model_vf_vfi(self, param)
    est_ccp = rust_model_pr_a(param, est_v)
    est_hazard = self%lambda * est_ccp

    ! Report several comparisons
    if (nx <= NX_SMALL) then
       print *
       print '(a12,100f12.8)', 'Model hazard:', model_hazard
       print '(a12,100f12.8)', 'Est. hazard:', est_hazard
       print '(a12,100f12.8)', 'Abs. Diff.:', abs(model_hazard - est_hazard)
    else
       print '(a12,100f12.8)', 'Sum Abs. Diff.:', sum(abs(model_hazard - est_hazard))
    end if
  end subroutine rust_model_compare_hazards

  ! Read individual parameters in THETA based on settings and
  ! construct a parameter structure, returned in PARAM.
  subroutine read_parameters(self, param, theta, model_variant)
    type(rust_model_t), intent(inout) :: self
    real(wp), dimension(:), intent(in) :: theta
    type(param_t), intent(out) :: param
    integer, intent(in) :: model_variant

    if (model_variant == RUST_MODEL_HETEROGENEOUS) then
       param%lambda1 = theta(IDX_LAM1)
       param%lambda2 = theta(IDX_LAM2)
       param%q1 = theta(IDX_Q1_H)
       param%beta = theta(IDX_BETA_H)
       param%cost = theta(IDX_COST_H)

       ! Set arrival rates: lambda1 for low states, lambda2 for high states
       self%lambda(1:self%nx_low) = param%lambda1
       self%lambda(self%nx_low+1:nx) = param%lambda2
    else if (model_variant == RUST_MODEL_HOMOGENEOUS) then
       param%lambda1 = theta(IDX_LAM)
       param%lambda2 = theta(IDX_LAM)
       param%q1 = theta(IDX_Q1)
       param%beta = theta(IDX_BETA)
       param%cost = theta(IDX_COST)

       ! Variable, but identical arrival rates
       self%lambda(1:nx) = param%lambda1
    else
       ! Fixed lambda values for ABBE model
       param%lambda1 = 1.0_wp
       param%lambda2 = 1.0_wp
       param%q1 = theta(IDX_Q1_A)
       param%beta = theta(IDX_BETA_A)
       param%cost = theta(IDX_COST_A)

       ! Set arrival rates = 1.0
       self%lambda(1:self%nx_low) = param%lambda1
       self%lambda(self%nx_low+1:nx) = param%lambda2
    end if
  end subroutine read_parameters


  ! Given a vector V of length N, return the value
  ! LOG(EXP(V(1)) + ... + EXP(V(N))).  In order to prevent numerical
  ! overflow, it first shifts the vector V by a constant C, applies
  ! EXP and LOG, and then adjusts the result to account for the shift.
  function log_sum_exp(v) result(e)
    real(wp), dimension(:), intent(in) :: v
    real(wp), dimension(size(v)) :: expv
    real(wp), parameter :: bound = 700.0_wp
    real(wp) :: e, c

    ! Recenter by subtracting C, the largest element of V.
    c = maxval(v)
    expv = v - c

    ! If recentered vector is still poorly scaled, use MAX as approximation.
    if (maxval(expv) > bound) then
       e = c
    else
       expv = bexp(expv)
       e = blog(sum(expv)) + c
    end if
  end function log_sum_exp


  ! Return LOG(X), with bounds imposed to prevent underflow.
  elemental function blog(x) result(l)
    real(wp), intent(in) :: x
    real(wp) :: l
    real(wp), parameter :: bound = 1.0e-307_wp
    real(wp), parameter :: l_bound = log(bound)

    if (x < bound) then
       l = l_bound
    else
       l = log(x)
    end if
  end function blog


  ! Return EXP(X), with bounds imposed to prevent overflow.
  elemental function bexp(x) result(e)
    real(wp), intent(in) :: x
    real(wp) :: e
    real(wp), parameter :: bound = 700.0_wp
    real(wp), parameter :: e_bound = exp(bound)

    if (x > bound) then
       e = e_bound
    else
       e = exp(x)
    end if
  end function bexp


  ! Divide NUMERATOR by DENOMINATOR only when DENOMINATOR is non-zero.
  ! This is the real(wp) version.
  function divide_maybe_vwp(numerator, denominator) result(quotient)
    real(wp), dimension(:), intent(in) :: numerator
    real(wp), dimension(size(numerator)), intent(in) :: denominator
    real(wp), dimension(size(numerator)) :: quotient
    integer :: ix

    do ix = 1, nx
       if (abs(denominator(ix)) > 1.0e-8_wp) then
          quotient(ix) = numerator(ix) / denominator(ix)
       else
          quotient(ix) = 0.0_wp
       end if
    end do
  end function divide_maybe_vwp


  ! Divide NUMERATOR by DENOMINATOR only when DENOMINATOR is non-zero.
  ! This is the integer version.
  function divide_maybe_vi(numerator, denominator) result(quotient)
    integer, dimension(:), intent(in) :: numerator
    integer, dimension(size(numerator)), intent(in) :: denominator
    integer, dimension(size(numerator)) :: quotient
    integer :: ix

    do ix = 1, nx
       if (abs(denominator(ix)) > 0) then
          quotient(ix) = numerator(ix) / denominator(ix)
       else
          quotient(ix) = 0.0_wp
       end if
    end do
  end function divide_maybe_vi

end module rust_model
