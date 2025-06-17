! lbfgsb.f90 -- wrapper functions for calling L-BFGS-B routines
!
! Copyright (C) 2008-2020 Jason R. Blevins
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

module lbfgsb
  use osl, wp => osl_wp
  implicit none
  private

  public :: lbfgsb_min

contains

  ! Minimize a function FUNC_N using the L-BFGS-B algorithm of Zhu,
  ! Byrd, Lu, and Nocedal.  L-BFGSB-B carries the following condition
  ! for use:
  !
  !     This software is freely available, but we expect that all
  !     publications describing work using this software, or all
  !     commercial products using it, quote at least one of the
  !     references given below.
  !
  !   * R. H. Byrd, P. Lu and J. Nocedal (1995).
  !     A Limited Memory Algorithm for Bound Constrained Optimization.
  !     SIAM Journal on Scientific and Statistical Computing, 16, 1190-1208.
  !
  !   * C. Zhu, R. H. Byrd and J. Nocedal (1997).
  !     L-BFGS-B: Algorithm 778: L-BFGS-B, FORTRAN routines for large scale
  !     bound constrained optimization.
  !     ACM Transactions on Mathematical Software, 23, 550-560.
  subroutine lbfgsb_min(func_n, x, lb, ub, f, eps, pgtol, iprint, ctx)

    interface
       subroutine func_n(x, y, ctx)
         use osl, wp => osl_wp
         implicit none
         real(wp), dimension(:), intent(in)  :: x
         real(wp), intent(out) :: y
         integer, dimension(:), intent(in), optional :: ctx
       end subroutine func_n
    end interface

    ! input variables
    real(wp), dimension(:), intent(inout) :: x
    real(wp), dimension(size(x)), intent(in) :: lb
    real(wp), dimension(size(x)), intent(in) :: ub

    ! output variables
    real(wp), intent(out) :: f

    ! optional variables
    real(wp), intent(in), optional :: eps
    integer, intent(in), optional :: iprint
    real(wp), intent(in), optional :: pgtol

    ! Function context
    integer, dimension(:), intent(in), optional :: ctx

    ! Local L-BFGS-B variables
    real(wp), dimension(size(x)) :: g
    logical :: terminate
    integer, dimension(size(x)) :: nbd

    ! Timer variables
    integer :: count0, count1, count_rate, time
    real(wp) :: time0, time1

    ! L-BFGS-B parameters
    real(wp) :: opt_eps
    integer :: opt_iprint
    logical, parameter :: max = .true.
    integer, parameter :: m = 10
    real(wp), parameter :: factr = 0.0_wp
    real(wp) :: opt_pgtol

    ! workspaces
    character(len=60) :: task, csave
    real(wp), dimension(29) :: dsave
    real(wp), dimension(:), allocatable :: wa
    integer, dimension(3*size(x)) :: iwa
    integer, dimension(44) :: isave
    logical, dimension(4) :: lsave

    external :: SETULB

    ! Start timer
    call system_clock(count0)
    call cpu_time(time0)

    f = 0.0_wp
    g = 0.0_wp
    nbd = 2

    allocate(wa(2*m*size(x)+5*size(x)+11*m*m+8*m))

    ! Optional parameters
    if (present(eps)) then
       opt_eps = eps
    else
       opt_eps = 1.0e-8_wp
    end if
    if (present(iprint)) then
       opt_iprint = iprint
    else
       opt_iprint = -1
    end if
    if (present(pgtol)) then
       opt_pgtol = pgtol
    else
       opt_pgtol = 1.0e-5_wp
    end if

    ! Start L-BFGS-B by calling SETULB with the 'START' task.
    task = 'START'
    terminate = .false.

    do while (.not. terminate)

       ! Call L-BFGS-B
       call SETULB(size(x), m, x, lb, ub, nbd, f, g, factr, opt_pgtol, wa, iwa, &
           task, opt_iprint, csave, lsave, isave, dsave)

       if (task(1:2) .eq. 'FG') then

          ! evaluate the function and approximate its gradient
          call func_n(x, f, ctx)

          call ngrad_central(func_n, x, g, EPS=opt_eps, CTX=ctx)

       else if (task(1:5) .eq. 'NEW_X') then

          ! intermediate output should be reported here and any
          ! alternative stopping criteria should be checked.

          if (opt_iprint > 0) then
             print '("Iteration ", i6, "  F = ", g17.5, "  |G| = ", g17.5)', &
                  isave(30), f, maxval(abs(g))
             print '("X = ", 100g17.5, /)', x
          end if

       else

          terminate = .true.

       endif

    end do

    deallocate(wa)

    ! Stop timer
    call system_clock(count1, count_rate)
    call cpu_time(time1)
    time = int((count1 - count0) / real(count_rate))
    if (opt_iprint > 0) then
       print '(a,i0,a,i0,a)', 'L-BFGS-B completed in ', &
            time, ' sec. (', int(time1 - time0), ' sec. cpu time)'
    end if
  end subroutine lbfgsb_min

  ! Calculate a finite-difference gradient of an arbitrary real
  ! function func_n at a given point x.  The gradient is returned in
  ! df.  The stepsize eps is used if given, otherwise, epsilon is
  ! used.
  subroutine ngrad_central(func_n, x, df, eps, ctx)

    interface
       subroutine func_n(x, y, ctx)
         use osl, wp => osl_wp
         implicit none
         real(wp), dimension(:), intent(in)  :: x
         real(wp), intent(out) :: y
         integer, dimension(:), intent(in), optional :: ctx
       end subroutine func_n
    end interface

    real(wp), dimension(:), intent(in) :: x
    real(wp), dimension(size(x)), intent(out) :: df
    real(wp), intent(in) :: eps
    integer, dimension(:), intent(in), optional :: ctx

    real(wp), dimension(size(x)) :: h, xph, xmh, xcur
    real(wp) :: f_plus, f_minus
    integer :: i

    h = eps**(1.0_wp / 3.0_wp) * max(abs(x), 1.0_wp)
    xph = x + h
    xmh = x - h
    h = xph - xmh
    xcur = x

    do i = 1, size(x)
       xcur(i) = xph(i)
       call func_n(xcur, f_plus, ctx)
       xcur(i) = xmh(i)
       call func_n(xcur, f_minus, ctx)
       df(i) = 0.5_wp * (f_plus - f_minus) / h(i)
       xcur(i) = x(i)
    end do
  end subroutine ngrad_central

end module lbfgsb
