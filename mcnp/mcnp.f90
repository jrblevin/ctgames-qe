! mcnp.f90 --- control program for quality ladder Monte Carlo experiments
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

program mcnp
  use iso_fortran_env
  use omp_lib
  use osl, wp => osl_wp
  use model
  use lbfgsb
  implicit none

  ! theta: [ lambda_odd, lambda_even, gamma, kappa, eta, fc ]
  real(wp), dimension(np) :: theta0 = [ 0.5, 1.5, 0.2, 0.8, 4.0, 0.1 ]

  ! Monte Carlo
  integer :: nmc = 25                             ! Monte Carlo trials
  integer :: nm  = 40000                          ! number of markets
  integer :: nt  = 1                              ! number of time periods
  integer :: r, k                                 ! trial and starting value indices
  integer :: seed = 274                           ! seed value
  real(wp) :: vf_tol_tmp                          ! Value function tolerance
  real(wp) :: Q_tol_tmp                           ! Matrix exponential tolerance
  integer, parameter :: dgp_seed = 274            ! Initial DGP seed value

  ! L-BFGS-B Optimization
  real(wp) :: eps_lbfgs = 1.0e-8_wp
  real(wp) :: pgtol_lbfgs = 1.0e-9_wp
  integer :: iprint_lbfgs = 1
  integer, parameter :: n_start = 3
  real(wp), dimension(np,n_start) :: theta_start
  real(wp), dimension(np) :: theta_current, theta_best
  real(wp) :: ll_current, ll_best

  ! Monte Carlo Results
  real(wp), dimension(np) :: theta_mean, theta_stdev, theta_bias, theta_rmse, theta_median_bias
  real(wp), dimension(:,:), allocatable :: theta_mc
  real(wp), dimension(:), allocatable :: ll_mc, ll_true
  integer :: start_mc = 0

  ! Timer variables
  integer :: count0, count1, count_rate
  real(wp) :: time

  ! Command-line arguments
  character(len=1000) :: ctrl_file
  character(len=1000) :: restart_file

  call osl_print_heading('CTGAMES MCNP QUALITY LADDER MONTE CARLO EXPERIMENTS')

  print *, 'Compiler: ' // compiler_version()
  print *, 'Compiler options: ' // compiler_options()
  print *, 'OMP Num. Processors: ', omp_get_num_procs()
  !$OMP PARALLEL
  if (omp_get_thread_num() == 0) then
     print *, 'OMP Num. Threads: ', omp_get_num_threads()
  end if
  !$OMP END PARALLEL

  ! Default values
  start_mc = 1
  restart_file = 'restart.ctl'

  ! Read command-line arguments
  if (command_argument_count() < 1) then
     print '(a)', 'Usage: mcnp <control_file>'
     stop
  end if
  call get_command_argument(1, ctrl_file)
  call read_control_file(ctrl_file)

  call osl_print_subheading('Initialization')
  call model_init()

  ! Initialize L-BFGS-B starting values
  theta_start(:,1) = [ 0.3_wp, 0.8_wp, 0.15_wp, 0.5_wp, 3.0_wp, 0.08_wp ]
  theta_start(:,2) = [ 0.7_wp, 1.2_wp, 0.25_wp, 1.0_wp, 5.0_wp, 0.12_wp ]
  theta_start(:,3) = [ 0.5_wp, 1.5_wp, 0.2_wp, 0.8_wp, 4.0_wp, 0.1_wp ]

  call osl_print_subheading('Monte Carlo parameters')
  call osl_print('nmc:', nmc)
  call osl_print('nm:', nm)
  call osl_print('theta0:', theta0)

  if (nmc == 0) then
     ! Solve and report timing only
     call model_reset(theta0)
  else
     ! If restarting previous Monte Carlo experiment
     if (start_mc > 0) then
        call osl_print("Previously completed Monte Carlo experiments: ", start_mc - 1)
        call osl_print('theta_mc:', transpose(theta_mc(:,1:start_mc-1)))
     end if
  end if

  print *
  call model_print()

  ! Monte Carlo experiments
  do r = start_mc, nmc
     print '(/,A18,I4,/,22("="))', 'Monte Carlo trial ', r

     ! Random number generator seed
     seed = dgp_seed + 37 * r
     print '("Seed: ", i8)', seed

     ! Generate data
     ! Use higher tolerances for DGP. Save and restore estimation settings.
     vf_tol_tmp = vf_tol
     vf_tol = 1.0e-12_wp
     Q_tol_tmp = Q_tol
     Q_tol = 1.0e-16_wp
     call model_generate_data(theta0, nm, nt, seed)
     call model_dataset_stats()
     vf_tol = vf_tol_tmp
     Q_tol = Q_tol_tmp

     ! Estimate all parameters jointly
     call osl_print_subheading("Estimation")

     ! Calculate likelihood at true values, noting that the model
     ! returns the negative log likelihood for minimization.
     call model_log_likelihood(theta0, ll_true(r))
     ll_true(r) = -ll_true(r)
     call osl_print('ll(theta0) = ', ll_true(r))

     ! Start timer
     call system_clock(count0)

     ! Optimization by L-BFGS-B with multiple starting values
     ll_best = HUGE(1.0_wp)
     theta_best = 0.0_wp

     do k = 1, n_start
        theta_current = theta_start(:,k)
        call lbfgsb_min(model_log_likelihood, theta_current, &
                        theta_lb, theta_ub, ll_current, &
                        EPS=eps_lbfgs, PGTOL=pgtol_lbfgs, IPRINT=iprint_lbfgs)

        ! Check if this run is better (note that we're minimizing -ll)
        if (ll_current < ll_best) then
           ll_best = ll_current
           theta_best = theta_current
        endif
     end do

     theta_mc(:,r) = theta_best
     ll_mc(r) = -ll_best

     ! Report wall time for current experiment
     call system_clock(count1, count_rate)
     time = real(count1 - count0) / real(count_rate)
     print '(a,i0,a,i0,a)', 'Estimation completed in ', int(time), ' sec.'

     ! Save state in case program is terminated prematurely.
     call write_control_file(restart_file)

     ! Print results
     call osl_print('theta0     = ', theta0)
     call osl_print('theta      = ', theta_mc(:,r))
     call osl_print('ll(theta0) = ', ll_true(r))
     call osl_print('ll(theta)  = ', ll_mc(r))

     ! Check to see if we did at least as good as ll(theta0).
     ! This is an infeasible value, but it is useful as a diagnostic.
     if (ll_mc(r) < ll_true(r)) then
        call osl_print('WARNING: Log likelihood lower at estimates than (infeasible) true parameters.')
     end if

     ! Print results
     if (r > 1) then
        if (r < nmc) then
           call osl_print_heading('Intermediate Monte Carlo Results')
        else
           call osl_print_heading('Final Monte Carlo Results')
           call osl_print('theta_mc:', theta_mc)
        end if
        call osl_sample_stats(theta_mc(:,1:r), &
             MEAN=theta_mean, SD=theta_stdev, &
             TRUTH=theta0, BIAS=theta_bias, RMSE=theta_rmse)
        theta_median_bias = osl_median(theta_mc(:,1:r)) - theta0
        call osl_print('Mean', theta_mean)
        call osl_print('Std. Dev.', theta_stdev)
        call osl_print('Mean Bias', theta_bias)
        call osl_print('Median Bias', theta_median_bias)
        call osl_print('RMSE', theta_rmse)
        call osl_print('Mean LL', sum(ll_mc(1:r)) / real(r, wp))
     end if
  end do

  call model_free()

  deallocate(theta_mc)
  deallocate(ll_mc)
  deallocate(ll_true)

contains

  ! Read whitespace-separated control file FILENAME.
  subroutine read_control_file(filename)
    character(len=*), intent(in) :: filename
    integer, parameter :: fh = 15

    ! Input related variables
    character(len=1000) :: buffer, label
    integer :: pos
    integer :: ios = 0
    integer :: line = 0
    integer :: theta_mc_count = 0
    logical :: file_exists

    ! Check to see if the control file exists
    inquire(FILE=filename, EXIST=file_exists)
    if (.not. file_exists) then
       print '(a,a)', 'Error: Control file does not exist: ', trim(filename)
       stop
    end if

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

          ! Model
          case ('nn')
             read(buffer, *, iostat=ios) nn
          case ('nw')
             read(buffer, *, iostat=ios) nw
          case ('mktsize')
             read(buffer, *, iostat=ios) mktsize
          case ('we')
             read(buffer, *, iostat=ios) we
          case ('rho')
             read(buffer, *, iostat=ios) rho
          case ('maxiter')
             read(buffer, *, iostat=ios) maxiter
          case ('vf_tol')
             read(buffer, *, iostat=ios) vf_tol
          case ('Q_tol')
             read(buffer, *, iostat=ios) Q_tol
          case ('discrete')
             read(buffer, *, iostat=ios) DISCRETE_TIME
          case ('delta')
             read(buffer, *, iostat=ios) DELTA
          case ('nvf')
             read(buffer, *, iostat=ios) nvf
          case ('prmod')
             read(buffer, *, iostat=ios) prmod

          ! Model Parameters
          case ('theta')
             read(buffer, *, iostat=ios) theta0

          ! Monte Carlo Parameters
          case ('nmc')
             read(buffer, *, iostat=ios) nmc
             allocate(theta_mc(np, nmc))
             allocate(ll_mc(nmc))
             allocate(ll_true(nmc))
          case ('nm')
             read(buffer, *, iostat=ios) nm
          case ('nt')
             read(buffer, *, iostat=ios) nt
          case ('start_mc')
             read(buffer, *, iostat=ios) start_mc
          case ('theta_mc')
             ! Read series of START_MC-1 real parameter vectors
             theta_mc_count = theta_mc_count + 1
             read(buffer, *, iostat=ios) theta_mc(:,theta_mc_count)

          ! L-BFGS-B Optimization
          case ('eps_lbfgs')
             read(buffer, *, iostat=ios) eps_lbfgs
          case ('pgtol_lbfgs')
             read(buffer, *, iostat=ios) pgtol_lbfgs
          case ('iprint_lbfgs')
             read(buffer, *, iostat=ios) iprint_lbfgs
          case ('restart')
             read(buffer, *, iostat=ios) restart_file

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

  ! Write a control file suitable for restarting FILENAME.
  subroutine write_control_file(filename)
    character(len=*), intent(in) :: filename
    integer, parameter :: fh = 15
    integer :: j

    open(fh, file=trim(filename))

    call write_comment(fh, trim(filename) // '     -*-conf-*-')
    call write_blank(fh)
    call write_comment(fh, 'Model')
    call write_int(fh, 'nn', nn)
    call write_int(fh, 'nw', nw)
    call write_int(fh, 'we', we)
    call write_rwp(fh, 'mktsize', mktsize)
    call write_rwp(fh, 'rho', rho)
    call write_int(fh, 'maxiter', maxiter)
    call write_rwp(fh, 'vf_tol', vf_tol)
    call write_rwp(fh, 'Q_tol', Q_tol)
    call write_log(fh, 'discrete', DISCRETE_TIME)
    call write_rwp(fh, 'delta', DELTA)
    call write_int(fh, 'nvf', nvf)
    call write_str(fh, 'restart', restart_file)

    call write_blank(fh)
    call write_comment(fh, 'Model Parameters')
    call write_vec(fh, 'theta', theta0)

    call write_blank(fh)
    call write_comment(fh, 'L-BFGS-B Optimization')
    call write_rwp(fh, 'eps_lbfgs', eps_lbfgs)
    call write_rwp(fh, 'pgtol_lbfgs', pgtol_lbfgs)
    call write_int(fh, 'iprint_lbfgs', iprint_lbfgs)

    call write_blank(fh)
    call write_comment(fh, 'Monte Carlo parameters')
    call write_int(fh, 'nmc', nmc)
    call write_int(fh, 'nm', nm)
    call write_int(fh, 'nt', nt)

    call write_blank(fh)
    call write_comment(fh, 'Previous Results')
    call write_int(fh, 'start_mc', r+1)
    do j = 1, r
       call write_vec(fh, 'theta_mc', theta_mc(:,j)) ! Previous Monte Carlo estimates
    end do
    close(fh)
  end subroutine write_control_file

  subroutine write_blank(fh)
    integer, intent(in) :: fh
    write(fh, '("")')
  end subroutine write_blank

  subroutine write_comment(fh, comment)
    integer, intent(in) :: fh
    character(len=*), intent(in) :: comment
    write(fh, '("# ", A)') adjustl(comment)
  end subroutine write_comment

  subroutine write_int(fh, label, n)
    integer, intent(in) :: fh
    character(len=*), intent(in) :: label
    integer, intent(in) :: n
    character(len=20) :: label_
    label_ = label
    write(fh, '(a11,i0)') adjustl(label_), n
  end subroutine write_int

  subroutine write_rwp(fh, label, x)
    integer, intent(in) :: fh
    character(len=*), intent(in) :: label
    real(wp), intent(in) :: x
    character(len=20) :: label_
    label_ = label
    write(fh, '(a11,g22.15)') adjustl(label_), x
  end subroutine write_rwp

  subroutine write_vec(fh, label, v)
    integer, intent(in) :: fh
    character(len=*), intent(in) :: label
    real(wp), dimension(:), intent(in) :: v
    character(len=20) :: label_
    label_ = label
    write(fh, '(a11,999g22.15)') adjustl(label_), v
  end subroutine write_vec

  subroutine write_log(fh, label, b)
    integer, intent(in) :: fh
    character(len=*), intent(in) :: label
    logical, intent(in) :: b
    character(len=20) :: label_
    label_ = label
    write(fh, '(a11,L1)') adjustl(label_), b
  end subroutine write_log

  subroutine write_str(fh, label, str)
    integer, intent(in) :: fh
    character(len=*), intent(in) :: label
    character(len=*), intent(in) :: str
    character(len=20) :: label_
    label_ = label
    write(fh, '(a11,a)') adjustl(label_), trim(str)
  end subroutine write_str

end program mcnp
