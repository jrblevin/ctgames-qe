! matrix.f90 -- matrix procedures
!
! Copyright (C) 2007-2009 Jason R. Blevins
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

module matrix
  use osl, wp => osl_wp
  implicit none

contains

  ! Set a matrix equal to the identity matrix.
  subroutine identity(I)
    real(wp), dimension(:,:), intent(inout) :: I
    integer :: j

    I = 0.0_wp
    do j = 1, size(I,1)
       I(j,j) = 1.0_wp
    end do
  end subroutine identity


  ! Calculate exp(t*H) for an N-by-N matrix H using Expokit.
  function expm(t, H) result(expH)
    real(wp), intent(in) :: t
    real(wp), dimension(:,:), intent(in) :: H
    real(wp), dimension(size(H,1),size(H,2)) :: expH

    ! Expokit variables
    external :: DGPADM
    integer, parameter :: ideg = 6
    real(wp), dimension(4*size(H,1)*size(H,2) + ideg + 1) :: wsp
    integer, dimension(size(H,1))  :: iwsp
    integer :: iexp, ns, iflag, n

    if (size(H,1) /= size(H,2)) then
       stop 'expm: matrix must be square'
    end if

    n = size(H,1)
    call DGPADM(ideg, n, t, H, n, wsp, size(wsp,1), iwsp, iexp, ns, iflag)
    expH = reshape(wsp(iexp:iexp+n*n-1), shape(expH))
  end function expm


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

end module matrix
