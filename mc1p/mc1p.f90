! mc1p.f90 -- single agent Monte Carlo control program
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

program mc1p
  use iso_fortran_env
  use omp_lib
  use osl, wp => osl_wp
  use rust_model
  use lbfgsb

  implicit none

  ! Monte Carlo parameters
  integer :: nmc = 100                            ! Monte Carlo trials
  integer, parameter :: n_start = 3               ! Number of starting values

  ! Model settings
  integer :: data_type
  integer :: nm
  real(wp) :: max_t
  real(wp) :: delta
  integer, parameter :: np = np_full

  ! Model instances
  type(rust_model_t), target :: model
  type(rust_model_ptr_t) :: model_ptr
  integer, dimension(:), allocatable :: model_ctx

  ! Monte Carlo
  integer :: r, k, med_idx
  integer, parameter :: nb = np + 3 ! Number of final statistics stored
  real(wp), dimension(:,:), allocatable :: theta_mc
  real(wp), dimension(:,:), allocatable :: beta_mc
  real(wp), dimension(:), allocatable :: theta_mean, theta_stdev, theta_mean_bias, theta_median_bias, theta_rmse, theta_sort
  real(wp), dimension(:), allocatable :: beta_mean, beta_stdev, beta_mean_bias, beta_median_bias, beta_rmse, beta_sort
  real(wp), dimension(:), allocatable :: ll_theta, ll_theta0
  integer, parameter :: seed_base = 20150403
  character(len=30) :: header, data_str
  real(wp) :: ll_best, ll_current
  real(wp), dimension(np) :: theta_current, theta_best

  ! Parameter vectors
  real(wp), dimension(np) :: theta0 = [ 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp ] ! Placeholder
  real(wp), dimension(np,n_start) :: theta_start ! Starting values array
  real(wp), dimension(nb) :: beta0

  ! Timer variables
  integer :: count0, count1, count_rate
  real(wp) :: time0, time1, time

  ! Command-line arguments
  character(len=100) :: ctrl_file

  call osl_print_heading('CTGAMES SINGLE PLAYER MONTE CARLO EXPERIMENTS')

  print *, 'Compiler: ' // compiler_version()
  print *, 'Compiler options: ' // compiler_options()

  ! Read command-line arguments
  if (command_argument_count() < 1) then
     print *
     print '(a)', 'Usage: mc1p <control_file>'
     stop
  end if
  call get_command_argument(1, ctrl_file)
  call read_control_file(ctrl_file)

  ! Allocate memory
  allocate(theta_mc(np, nmc), beta_mc(nb, nmc))
  allocate(theta_sort(nmc), beta_sort(nmc))
  allocate(theta_mean(np), theta_stdev(np), theta_mean_bias(np), theta_median_bias(np), theta_rmse(np))
  allocate(beta_mean(nb), beta_stdev(nb), beta_mean_bias(nb), beta_median_bias(nb), beta_rmse(nb))
  allocate(ll_theta(nmc), ll_theta0(nmc))

  ! True final statistics
  beta0(1:np) = theta0                            ! Structural parameters
  beta0(np+1) = theta0(5) / theta0(4)             ! Replacement cost ratio
  beta0(np+2) = 0.0_wp                            ! LL of estimates - LL at truth
  beta0(np+3) = 0.0_wp                            ! Clock time

  ! Starting values
  theta_start(:,1) = [ 0.1_wp, 0.2_wp, 1.0_wp, -5.0_wp, -10.0_wp ]
  theta_start(:,2) = [ 0.5_wp, 0.5_wp, 0.5_wp, -5.0_wp, -15.0_wp ]
  theta_start(:,3) = [ 1.0_wp, 1.0_wp, 1.0_wp, -1.0_wp, -10.0_wp ]

  ! Perform nmc Monte Carlo experiments
  call osl_print('Monte Carlo Experiments')
  call osl_print('True parameters:')
  print '(5a12)', 'lambda1     ', 'lambda2     ', 'q1          ', 'beta        ', 'c           '
  print '(5g12.5)', theta0

  !$OMP PARALLEL DO PRIVATE(model, model_ptr, model_ctx, count0, count1, count_rate, &
  !$OMP                     time0, time1, time, theta_current, theta_best, ll_current, &
  !$OMP                     ll_best, k)
  do r = 1, nmc
     ! Initialize the model
     call rust_model_init(model, data_type, nm, max_t, DELTA=delta, &
                          model_variant=RUST_MODEL_HETEROGENEOUS)

     ! Set up a context parameter
     model_ptr%p => model
     model_ctx = transfer(model_ptr, osl_ctx)

     ! Generate data
     call rust_model_dgp(model, THETA=theta0, SEED=seed_base + r * 1708)

     ! Start timer
     call system_clock(count0)
     call cpu_time(time0)

     ! Evaluate log likelihood at truth for reference
     call rust_model_ll(theta0, ll_theta0(r), CTX=model_ctx)

     ! Estimate the model using multiple starting values
     ll_best = HUGE(1.0_wp)   ! Initialize best_ll to a very large number
     theta_best = 0.0_wp      ! Initialize to zero or a default value
     do k = 1, n_start
       ! Perform optimization using L-BFGS-B
       theta_current = theta_start(:,k)
       call lbfgsb_min(rust_model_ll, theta_current, &
                       model%theta_lb, model%theta_ub, ll_current, &
                       EPS=1.0e-13_wp, PGTOL=1.0e-9_wp, CTX=model_ctx)

       ! Check if this run is better (note that we're minimizing -ll)
       if (ll_current < ll_best) then
           ll_best = ll_current
           theta_best = theta_current
       endif
     end do
     ll_theta(r) = ll_best
     theta_mc(:,r) = theta_best

     ! Stop timer
     call system_clock(count1, count_rate)
     call cpu_time(time1)
     time = real(count1 - count0) / real(count_rate)

     ! Calculate reported statistics for trial r
     beta_mc(1:np, r) = theta_mc(:,r)                 ! Structural parameters
     beta_mc(np+1, r) = theta_mc(5,r) / theta_mc(4,r) ! Ratio of replacement cost to mileage cost
     ! Difference of log likelihoods
     ! Since the log likelihood function was negated for minimization, undo that here.
     ! <-------|-----------|--------|--------->
     !   -ll_theta0   -ll_theta     0
     beta_mc(np+2, r) = -(max_t * model%data%size) * (ll_theta(r) - ll_theta0(r))
     beta_mc(np+3, r) = time                          ! Clock time

     ! Clean up
     call rust_model_free(model)

  end do
  !$OMP END PARALLEL DO

  ! Report raw estimates
  call osl_print('Structural Parameter Estimates (theta): ', transpose(theta_mc))
  call osl_print('Auxiliary Estimates (beta): ', transpose(beta_mc))

  ! Report final results
  call osl_print_heading('Results')

  ! Report final second-stage results (theta)
  call osl_sample_stats(theta_mc, MEAN=theta_mean, SD=theta_stdev, &
       TRUTH=theta0, BIAS=theta_mean_bias, RMSE=theta_rmse)

  ! Report final statistics of interest (beta)
  call osl_sample_stats(beta_mc, MEAN=beta_mean, SD=beta_stdev, &
       TRUTH=beta0, BIAS=beta_mean_bias, RMSE=beta_rmse)

  ! Calculate median bias for each parameter
  med_idx = int(nmc / 2.0_wp)
  do k = 1, np
     theta_sort = theta_mc(k,:)
     call osl_quicksort(theta_sort)
     theta_median_bias(k) = theta_sort(med_idx) - theta0(k)
  end do
  do k = 1, nb
     beta_sort = beta_mc(k,:)
     call osl_quicksort(beta_sort)
     beta_median_bias(k) = beta_sort(med_idx) - beta0(k)
  end do

  print '(/,"Final Results:   lam1    lam2    q1      beta        c   c/beta  ll_diff walltime:")'
  print '("True:     ", 100f9.3)', beta0
  print '("Mean:     ", 100f9.3)', beta_mean
  print '("Mean Bias:", 100f9.3)', beta_mean_bias
  print '("Med Bias: ", 100f9.3)', beta_median_bias
  print '("S.D.:     ", 100f9.3)', beta_stdev
  print '("RMSE:     ", 100f9.3)', beta_rmse
  print '("Min:      ", 100f9.3)', minval(beta_mc, DIM=2)
  print '("Max:      ", 100f9.3)', maxval(beta_mc, DIM=2)

  ! Report structural parameter estimates in LaTeX format
  if (data_type == RUST_DATA_DT) then
     write (data_str, '(a,f4.2,a)') '$\Delta = ', delta, '$'
  else
     write (data_str, '(a,f4.2,a)') 'Continuous'
  end if
  print '(/,"Final Structural Parameters (theta):")'
  print '("\begin{table}")'
  print '("  \centering")'
  print '("  \begin{tabular}{rllrrrrrr}")'
  print '("    \hline")'
  print '("    \hline")'


  print '("    ", a8, " & ", a15, " & ", a6, 6("  &", a12), "  \\")', &
     "$M", "Sampling", "", "$\lambda_1$", "$\lambda_2$", "$\gamma$", "$\beta$", "$\mu$", "$\mu/\beta$"
  print '("    \hline")'
  print '("    ", a8, " & ", a15, " & ", a6, 6("  &", f12.3), "  \\")', &
     "$\infty$", "DGP", "True", beta0(1:np+1)
  print '("    \hline")'
  print '("    ", i4, "     & ", a15, " & Mean  ", 6("  &", f12.3), "  \\")', &
     nm, data_str, beta_mean(1:np+1)
  print '("             &                 & S.D.  ", 6("  &", f12.3), "  \\")', beta_stdev(1:np+1)
  print '("    \hline")'
  print '("  \end{tabular} \\")'
  print '("\end{table}")'

  ! Free memory
  deallocate(ll_theta, ll_theta0)
  deallocate(theta_mc, theta_mean, theta_stdev, theta_mean_bias, theta_median_bias, theta_rmse)
  deallocate(theta_sort, beta_sort)
  deallocate(beta_mc, beta_mean, beta_stdev, beta_mean_bias, beta_median_bias, beta_rmse)

contains

  ! Read whitespace-separated control file FILENAME.
  subroutine read_control_file(filename)
    character(len=*), intent(in) :: filename
    integer, parameter :: fh = 15

    ! Input related variables
    character(len=100) :: buffer, label
    integer :: pos
    integer :: ios = 0
    integer :: line = 0

    logical :: file_exists

    ! Check to see if the control file exists
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
       print '(a,a)', 'Error: Control file does not exist: ', trim(filename)
       stop
    end if

    ! Open the control file
    open(fh, file=filename)

    ! ios is negative if an end of record condition is encountered or
    ! if an endfile condition was detected.  It is positive if an
    ! error was detected.  ios is zero otherwise.
    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Act on blank lines as if the were comments
          if (buffer == ' ') buffer = '# '

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, ' ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)

          ! Experiment
          case ('nmc')
             read(buffer, *, iostat=ios) nmc

          ! DGP
          case ('theta0')
             read(buffer, *, iostat=ios) theta0
          case ('sample')
             select case (buffer)
             case ('discrete')
                data_type = RUST_DATA_DT
             case default
                data_type = RUST_DATA_CT
             end select
          case ('markets')
             read(buffer, *, iostat=ios) nm
          case ('max_t')
             read(buffer, *, iostat=ios) max_t
          case ('delta')
             read(buffer, *, iostat=ios) delta

          ! Comments
          case ('#')

          ! Error
          case default
             print *, 'Skipping invalid label at line', line
          end select
       end if
    end do

    close(fh)
  end subroutine read_control_file

end program mc1p
