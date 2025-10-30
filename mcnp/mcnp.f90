! mcnp.f90 --- control program for quality ladder Monte Carlo experiments
!
! Copyright (C) 2009-2025 Jason R. Blevins
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
  integer :: iprint_lbfgs = -1  ! Quiet mode by default
  integer, parameter :: n_start = 3
  real(wp), dimension(np,n_start) :: theta_start
  real(wp), dimension(np) :: theta_current, theta_best
  real(wp) :: ll_current, ll_best

  ! Monte Carlo Results
  real(wp), dimension(:,:), allocatable :: theta_mc
  real(wp), dimension(:), allocatable :: ll_mc, ll_true
  integer :: start_mc = 0
  integer :: end_mc = 0                           ! End replication (0 = use nmc)

  ! Command-line arguments
  character(len=1000) :: ctrl_file

  ! Environment variable support for parallel execution
  character(len=1000) :: env_val
  character(len=1000) :: results_file = ''
  integer :: total_cores

  call osl_print_heading('CTGAMES MCNP QUALITY LADDER MONTE CARLO EXPERIMENTS')

  call osl_print_subheading('System Information')
  print *, 'Compiler: ' // compiler_version()
  print *, 'Compiler options: ' // compiler_options()
  total_cores = omp_get_num_procs()
  print *, 'OMP Num. Processors: ', total_cores
  !$OMP PARALLEL
  if (omp_get_thread_num() == 0) then
     print *, 'OMP Num. Threads: ', omp_get_num_threads()
  end if
  !$OMP END PARALLEL

  ! Default values
  start_mc = 1

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
     start_mc = 0
     end_mc = -1
     call osl_print_subheading('Model')
     call model_reset(theta0, .true.)
     call model_print()

   else if (nmc > 0) then

      ! Monte Carlo experiments
     call osl_print_subheading('Monte Carlo experiments')

     ! Check for environment variable overrides (for parallel execution)
     call get_environment_variable('MC_START', env_val)
     if (len_trim(env_val) > 0) then
        read(env_val, *) start_mc
        print '(a,i0)', 'Environment override: MC_START = ', start_mc
     end if
     call get_environment_variable('MC_END', env_val)
     if (len_trim(env_val) > 0) then
        read(env_val, *) end_mc
        print '(a,i0)', 'Environment override: MC_END = ', end_mc
     end if
     if (end_mc == 0) end_mc = nmc
     call get_environment_variable('MC_RESULTS_FILE', results_file)

     if (end_mc < nmc) then
        print '(a,i0,a,i0,a,i0)', 'Running a subset of replications: ', start_mc, ' through ', end_mc, ' of ', nmc
     end if
  end if

  do r = start_mc, end_mc
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

     ! Print results
     call osl_print('theta      = ', theta_mc(:,r))
     call osl_print('ll(theta)  = ', ll_mc(r))

  end do

  ! Write results to file if requested (for parallel execution)
  if (len_trim(results_file) > 0) then
     call write_partial_results()
  end if

  call model_free()

  deallocate(theta_mc)
  deallocate(ll_mc)
  deallocate(ll_true)

contains

  ! Write partial Monte Carlo results to file
  subroutine write_partial_results()
    integer :: fh, r_idx

    open(newunit=fh, file=trim(results_file), action='write', status='replace')

    ! Write header
    write(fh, '(a)') '# Partial Monte Carlo Results'
    write(fh, '(a,i0,a,i0)') '# Replications ', start_mc, ' through ', end_mc
    write(fh, '(a)') '# lambda_L lambda_H gamma kappa eta fc'

    ! Write data (one line per replication)
    do r_idx = start_mc, end_mc
       write(fh, '(6(1x,es24.16))') theta_mc(:,r_idx)
    end do

    close(fh)
    print '(a,a)', 'Partial results written to: ', trim(results_file)
  end subroutine write_partial_results

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

end program mcnp
