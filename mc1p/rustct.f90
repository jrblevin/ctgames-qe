! rustct.f90 -- continuous time Rust model estimation
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

program rustct
  use iso_fortran_env
  use omp_lib
  use osl, wp => osl_wp
  use rust_model
  use lbfgsb

  implicit none

  ! Model variant and number of parameters (heterogeneous lambda)
  integer, parameter :: model_variant = RUST_MODEL_HETEROGENEOUS
  integer, parameter :: np_theta = np_full
  ! Use these for the homogeneous model (variable lambda)
  !integer, parameter :: model_variant = RUST_MODEL_HOMOGENEOUS
  !integer, parameter :: np_theta = np_restr
  ! Use these for the ABBE model (fixed lambda)
  !integer, parameter :: model_variant = RUST_MODEL_ABBE
  !integer, parameter :: np_theta = np_abbe

  ! Data groups
  integer, dimension(8), parameter :: groups = [ 1, 2, 3, 4, 5, 6, 7, 8 ]

  ! Model instance
  type(rust_model_t), target :: model
  type(rust_model_ptr_t) :: model_ptr
  integer, dimension(:), allocatable :: model_ctx

  ! Starting values
  integer, parameter :: nstart = 20
  integer :: j, j_best
  real(wp), dimension(np_theta, nstart) :: theta_start

  ! Estimates
  real(wp), dimension(np_theta) :: theta, se, theta_best
  real(wp) :: ll_theta, ll_scale, ll_best
  real(wp), dimension(np_theta, np_theta) :: H, cov

  ! Print heading
  call osl_print_heading('CTGAMES RUST (1987) ESTIMATION')
  print *, 'Compiler: ' // compiler_version()
  print *, 'Compiler options: ' // compiler_options()
  if (model_variant == RUST_MODEL_HETEROGENEOUS) then
     print *, 'Estimating heterogeneous model (5 parameters)'
  else if (model_variant == RUST_MODEL_HOMOGENEOUS) then
     print *, 'Estimating homogeneous model (4 parameters, lambda1 = lambda2)'
  else
     print *, 'Estimating ABBE model (3 parameters, lambda1 = lambda2 = 1.0)'
  end if

  ! Load starting values
  if (model_variant == RUST_MODEL_HETEROGENEOUS) then
     ! 5 parameters: lambda1, lambda2, q1, beta, cost
     theta_start(1:np_full,1)  = [ 0.1_wp, 0.2_wp, 2.0_wp,  -8.0_wp, -20.0_wp ]
     theta_start(1:np_full,2)  = [ 0.2_wp, 0.5_wp, 2.0_wp,  -2.0_wp, -20.0_wp ]
     theta_start(1:np_full,3)  = [ 1.0_wp, 1.0_wp, 1.0_wp,  -1.0_wp, -10.0_wp ]
     theta_start(1:np_full,4)  = [ 0.5_wp, 1.0_wp, 1.0_wp,  -1.0_wp, -20.0_wp ]
     theta_start(1:np_full,5)  = [ 0.1_wp, 0.2_wp, 0.5_wp,  -3.0_wp, -11.0_wp ]
     theta_start(1:np_full,6)  = [ 0.5_wp, 0.5_wp, 1.0_wp,  -5.0_wp, -30.0_wp ]
     theta_start(1:np_full,7)  = [ 2.0_wp, 2.0_wp, 2.0_wp,  -5.0_wp, -10.0_wp ]
     theta_start(1:np_full,8)  = [ 0.1_wp, 1.0_wp, 1.0_wp,  -1.0_wp,  -5.0_wp ]
     theta_start(1:np_full,9)  = [ 0.1_wp, 0.5_wp, 0.5_wp,  -0.5_wp,  -5.0_wp ]
     theta_start(1:np_full,10) = [ 1.0_wp, 2.0_wp, 0.5_wp,  -0.5_wp, -50.0_wp ]
     theta_start(1:np_full,11) = [ 0.3_wp, 0.4_wp, 2.5_wp, -10.0_wp, -25.0_wp ]
     theta_start(1:np_full,12) = [ 0.7_wp, 0.8_wp, 1.5_wp,  -7.0_wp, -30.0_wp ]
     theta_start(1:np_full,13) = [ 1.5_wp, 0.5_wp, 2.0_wp, -12.0_wp, -40.0_wp ]
     theta_start(1:np_full,14) = [ 2.0_wp, 1.5_wp, 0.5_wp, -20.0_wp, -15.0_wp ]
     theta_start(1:np_full,15) = [ 0.2_wp, 0.9_wp, 3.0_wp,  -3.0_wp, -35.0_wp ]
     theta_start(1:np_full,16) = [ 0.8_wp, 1.0_wp, 2.5_wp,  -8.0_wp, -20.0_wp ]
     theta_start(1:np_full,17) = [ 1.2_wp, 2.0_wp, 1.0_wp,  -5.0_wp, -50.0_wp ]
     theta_start(1:np_full,18) = [ 1.0_wp, 1.2_wp, 0.8_wp, -18.0_wp, -70.0_wp ]
     theta_start(1:np_full,19) = [ 0.4_wp, 0.6_wp, 2.2_wp, -15.0_wp, -60.0_wp ]
     theta_start(1:np_full,20) = [ 0.5_wp, 1.8_wp, 1.2_wp, -25.0_wp, -85.0_wp ]
  else if (model_variant == RUST_MODEL_HOMOGENEOUS) then
     ! 4 parameters: lambda, q1, beta, cost
     theta_start(1:np_restr,1)  = [ 0.1_wp, 2.0_wp,  -8.0_wp, -20.0_wp ]
     theta_start(1:np_restr,2)  = [ 0.2_wp, 2.0_wp,  -2.0_wp, -20.0_wp ]
     theta_start(1:np_restr,3)  = [ 1.0_wp, 1.0_wp,  -1.0_wp, -10.0_wp ]
     theta_start(1:np_restr,4)  = [ 0.5_wp, 1.0_wp,  -1.0_wp, -20.0_wp ]
     theta_start(1:np_restr,5)  = [ 0.1_wp, 0.5_wp,  -3.0_wp, -11.0_wp ]
     theta_start(1:np_restr,6)  = [ 0.5_wp, 1.0_wp,  -5.0_wp, -30.0_wp ]
     theta_start(1:np_restr,7)  = [ 2.0_wp, 2.0_wp,  -5.0_wp, -10.0_wp ]
     theta_start(1:np_restr,8)  = [ 0.1_wp, 1.0_wp,  -1.0_wp,  -5.0_wp ]
     theta_start(1:np_restr,9)  = [ 0.1_wp, 0.5_wp,  -0.5_wp,  -5.0_wp ]
     theta_start(1:np_restr,10) = [ 1.0_wp, 0.5_wp,  -0.5_wp, -50.0_wp ]
     theta_start(1:np_restr,11) = [ 0.3_wp, 2.5_wp, -10.0_wp, -25.0_wp ]
     theta_start(1:np_restr,12) = [ 0.7_wp, 1.5_wp,  -7.0_wp, -30.0_wp ]
     theta_start(1:np_restr,13) = [ 1.5_wp, 2.0_wp, -12.0_wp, -40.0_wp ]
     theta_start(1:np_restr,14) = [ 2.0_wp, 0.5_wp, -20.0_wp, -15.0_wp ]
     theta_start(1:np_restr,15) = [ 0.2_wp, 3.0_wp,  -3.0_wp, -35.0_wp ]
     theta_start(1:np_restr,16) = [ 0.8_wp, 2.5_wp,  -8.0_wp, -20.0_wp ]
     theta_start(1:np_restr,17) = [ 1.2_wp, 1.0_wp,  -5.0_wp, -50.0_wp ]
     theta_start(1:np_restr,18) = [ 1.0_wp, 0.8_wp, -18.0_wp, -70.0_wp ]
     theta_start(1:np_restr,19) = [ 0.4_wp, 2.2_wp, -15.0_wp, -60.0_wp ]
     theta_start(1:np_restr,20) = [ 0.5_wp, 1.2_wp, -25.0_wp, -85.0_wp ]
  else
     ! 3 parameters: q1, beta, cost (lambdas fixed at 1.0)
     theta_start(1:np_abbe,1)  = [ 2.0_wp,  -8.0_wp, -20.0_wp ]
     theta_start(1:np_abbe,2)  = [ 2.0_wp,  -2.0_wp, -20.0_wp ]
     theta_start(1:np_abbe,3)  = [ 1.0_wp,  -1.0_wp, -10.0_wp ]
     theta_start(1:np_abbe,4)  = [ 1.0_wp,  -1.0_wp, -20.0_wp ]
     theta_start(1:np_abbe,5)  = [ 0.5_wp,  -3.0_wp, -11.0_wp ]
     theta_start(1:np_abbe,6)  = [ 1.0_wp,  -5.0_wp, -30.0_wp ]
     theta_start(1:np_abbe,7)  = [ 2.0_wp,  -5.0_wp, -10.0_wp ]
     theta_start(1:np_abbe,8)  = [ 1.0_wp,  -1.0_wp,  -5.0_wp ]
     theta_start(1:np_abbe,9)  = [ 0.5_wp,  -0.5_wp,  -5.0_wp ]
     theta_start(1:np_abbe,10) = [ 0.5_wp,  -0.5_wp, -50.0_wp ]
     theta_start(1:np_abbe,11) = [ 2.5_wp, -10.0_wp, -25.0_wp ]
     theta_start(1:np_abbe,12) = [ 1.5_wp,  -7.0_wp, -30.0_wp ]
     theta_start(1:np_abbe,13) = [ 2.0_wp, -12.0_wp, -40.0_wp ]
     theta_start(1:np_abbe,14) = [ 0.5_wp, -20.0_wp, -15.0_wp ]
     theta_start(1:np_abbe,15) = [ 3.0_wp,  -3.0_wp, -35.0_wp ]
     theta_start(1:np_abbe,16) = [ 2.5_wp,  -8.0_wp, -20.0_wp ]
     theta_start(1:np_abbe,17) = [ 1.0_wp,  -5.0_wp, -50.0_wp ]
     theta_start(1:np_abbe,18) = [ 0.8_wp, -18.0_wp, -70.0_wp ]
     theta_start(1:np_abbe,19) = [ 2.2_wp, -15.0_wp, -60.0_wp ]
     theta_start(1:np_abbe,20) = [ 1.2_wp, -25.0_wp, -85.0_wp ]
   end if

  ! Initialize the model
  call rust_model_init(model, DATA_TYPE=RUST_DATA_DT, DELTA=1.0_wp, &
                       model_variant=model_variant)

  ! Set up a context parameter
  model_ptr%p => model
  model_ctx = transfer(model_ptr, osl_ctx)

  ! Load data
  call rust_model_load_data(model, groups)
  call rust_model_dataset_stats(model)

  ll_best = -1.0e99_wp
  j_best = -1

  call osl_print_subheading('Optimization runs:')

  do j = 1, nstart
    ! Maximize likelihood using L-BFGS-B with homogeneous wrapper
    theta = theta_start(:,j) ! j-th starting value
    call lbfgsb_min(rust_model_ll, theta, &
                    model%theta_lb, model%theta_ub, ll_theta, &
                    EPS=1.0e-13_wp, PGTOL=1.0e-9_wp, CTX=model_ctx)

    ! Account for negation due to minimization, and other scaling
    ll_scale = rust_model_get_ll_scaling(model)
    ll_theta = ll_theta / ll_scale
    print '("Starting value j = ", i3, "  LL = ", f20.13)', j, ll_theta
    if ((j == 1) .or. (ll_theta > ll_best)) then
        theta_best = theta
        ll_best = ll_theta
        j_best = j
    end if
  end do

  ! Report best starting value
  call osl_print_subheading('Best run:')
  call osl_print('Best starting value: j = ', j_best)
  call osl_print('LL = ', ll_best, FMT='f20.13')

  ! Calculate standard errors
  call nhessian(rust_model_ll, theta_best, H, EPS=1.0e-8_wp, DATA=model_ctx)
  cov = -ll_scale * inverse(H)
  se = sqrt(diag(cov))

  ! Report final results
  call osl_print_heading('Results')

  ! Report structural parameter estimates in LaTeX format
  print '("\begin{table}")'
  print '("  \centering")'
  if (model_variant == RUST_MODEL_HETEROGENEOUS) then
     print '("  \begin{tabular}{lccccc}")'
     print '("    \hline")'
     print '("    \hline")'
     print '("    Sampling   & $\lambda_1$ & $\lambda_2$ & $q_1$    & $\beta$ & $c$   \\")'
     print '("    \hline")'
     print '("    Estimates  ", 5("  &", f8.3), "  \\")', theta_best
     print '("    S.E.       ", 5("  &", f8.3), "  \\")', se
  else if (model_variant == RUST_MODEL_HOMOGENEOUS) then
     print '("  \begin{tabular}{lcccc}")'
     print '("    \hline")'
     print '("    \hline")'
     print '("    Sampling   & $\lambda$ & $q_1$    & $\beta$ & $c$   \\")'
     print '("    \hline")'
     print '("    Estimates", 4("  &", f8.3), "  \\")', theta_best
     print '("    S.E.     ", 4("  &", f8.3), "  \\")', se
  else
     print '("  \begin{tabular}{lcccc}")'
     print '("    \hline")'
     print '("    \hline")'
     print '("    Sampling   & $\lambda$ & $q_1$    & $\beta$ & $c$   \\")'
     print '("    \hline")'
     print '("    Estimates", "  & 1.000", 3("  &", f8.3), "  \\")', theta_best
     print '("    S.E.     ", "  &  --- ", 3("  &", f8.3), "  \\")', se
  end if
  print '("    \hline")'
  print '("  \end{tabular} \\")'
  print '("\end{table}")'

  ! Clean up
  call rust_model_free(model)

contains

  ! Calculate a finite-difference Hessian of an arbitrary real
  ! function func_n at a given point x.  The Hessian is returned in H.
  ! The stepsize eps is used if given, otherwise, epsilon is used.
  subroutine nhessian(func_n, x, H, eps, data)

    interface
       subroutine func_n(x, y, data)
         use osl, wp => osl_wp
         implicit none
         real(wp), dimension(:), intent(in)  :: x
         real(wp), intent(out) :: y
         integer, dimension(:), intent(in), optional :: data
       end subroutine func_n
    end interface

    real(wp), dimension(:), intent(in) :: x
    real(wp), dimension(size(x),size(x)), intent(out) :: H
    real(wp), optional, intent(in) :: eps
    integer, dimension(:), intent(in), optional :: data

    real(wp), dimension(size(x)) :: hh, xcur, xph, xmh
    real(wp) :: eps_
    real(wp) :: f00, fp0, fm0, fpp, fmm, fmp, fpm
    integer :: i, j

    if (present(eps)) then
       eps_ = eps
    else
       eps_ = 1.0e-8_wp
    end if

    hh = eps_**(1.0_wp / 4.0_wp) * max(abs(x), 1.0_wp)
    xph = x + hh
    xmh = x - hh
    hh = (xph - xmh) / 2.0_wp

    xcur = x
    call func_n(xcur, f00, data)

    do i = 1, size(x)
       ! Evaluate the two unidirectional steps and calculate the diagonal term
       xcur(i) = xph(i)
       call func_n(xcur, fp0, data)
       xcur(i) = xmh(i)
       call func_n(xcur, fm0, data)
       xcur(i) = x(i)

       H(i,i) = ( fp0 - 2.0_wp * f00 + fm0 ) / ( hh(i) * hh(i) )

       do j = i + 1, size(x)
          ! Evaluate the off-diagonal steps using bi-directional steps
          xcur(i) = xph(i)
          xcur(j) = xph(j)
          call func_n(xcur, fpp, data)
          xcur(j) = xmh(j)
          call func_n(xcur, fpm, data)
          xcur(i) = xmh(i)
          call func_n(xcur, fmm, data)
          xcur(j) = xph(j)
          call func_n(xcur, fmp, data)

          H(i,j) = ( fpp + fmm - fmp - fpm ) / ( 4.0_wp * hh(i) * hh(j) )
          H(j,i) = H(i,j)

          ! Return xcur to original point
          xcur(i) = x(i)
          xcur(j) = x(j)
       end do
    end do

  end subroutine nhessian

  ! Return the inverse of a matrix calculated by finding the LU
  ! decomposition.
  function inverse(A) result(Ainv)
    real(wp), dimension(:,:), intent(in) :: A
    real(wp), dimension(size(A,1),size(A,2)) :: Ainv

    real(wp), dimension(size(A,1)) :: work
    integer, dimension(size(A,1)) :: ipiv
    integer :: n, info

    external DGETRF
    external DGETRI

    if (size(A,1) /= size(A,2)) then
       stop 'ERROR in inverse: matrix must be square!'
    end if
    n = size(A,1)

    Ainv = A

    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       print '(a,i8,a)', 'Matrix is numerically singular! (info = ', &
            info, ')'
    end if

    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       print *, 'Matrix inversion failed!'
    end if
  end function inverse

  ! Return the diagonal of a square matrix.
  function diag(A) result(d)
    real(wp), dimension(:,:), intent(in) :: A
    real(wp), dimension(size(A,1)) :: d
    integer :: i

    if (size(A,1) /= size(A,2)) then
       stop 'ERROR in diag: matrix must be square!'
    end if

    do i = 1, size(A,1)
       d(i) = A(i,i)
    end do
  end function diag

end program rustct
