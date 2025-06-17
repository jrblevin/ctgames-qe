! encoding.f90 --- A module for encoding and decoding state vectors
!
! Copyright (C) 2008-2009 Jason R. Blevins
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

! DESCRIPTION
!
! This is a Fortran 2003 implementation of the "probability density"
! encoding algorithm presented in the following article:
!
!   Gowrisankaran, G. (1999). "Efficient representation of state
!   spaces for some dynamic models."  Journal of Economic Dynamics and
!   Control 23 (8), 1077-1098.
!
! This module implements an encoding of the state space $Z(N,M)$
! consisting of all discrete approximations to probability
! distributions on the interval $[0,1]$ with $N-1$ discrete regions
! and $N$ endpoints.  The mass assigned to any of these $N$ endpoints
! may take $M$ discrete values $0/(M-1), 1/(M-1), \dots, (M-1)/(M-1)$
! such that the elements sum to one.
!
! This encoding can be used, for example, to represent the space of
! all probability distributions or to represent the state space in a
! dynamic oligopoly model such as that of Ericson and Pakes (1995).
! In the dynamic oligopoly example, with $n$ firms and $L$ possible
! states, one can represent the state space using $Z(L,n+1)$.  The
! $k$-th element of the state vector then represents the fraction of
! firms in state $k$.

module encoding
  implicit none
  private

  public :: encoding_init, encoding_free
  public :: encode, decode

  integer, dimension(:,:), allocatable :: table   ! encoding table
  integer :: N_                                   ! number of endpoints
  integer :: M_                                   ! sum of values

contains

  ! Initializes the encoding table.
  function encoding_init(N, M) result(n_state)
    integer, intent(in) :: N, M
    integer :: i, j, n_state
    N_ = N
    M_ = M
    allocate(table(N,M))
    table(1,:) = 1
    table(:,1) = 1
    do i = 2, N
       do j = 2, M
          table(i,j) = table(i-1,j) + table(i, j-1)
       end do
    end do

    ! Cardinality of the state space (see Theorem 1, pg. 1084).
    n_state = table(N_, M_)
  end function encoding_init


  ! Frees memory used by the encoding table.
  subroutine encoding_free
    deallocate(table)
  end subroutine encoding_free


  ! Bijectively maps elements of the state space to integers from 1 to
  ! n_state.
  function encode(x) result(enc)
    integer, dimension(N_), intent(in) :: x
    integer :: enc
    integer :: i, xsum, idx

    enc = 0
    xsum = 0
    do i = 1, N_
       idx = N_ - i + 1
       enc = enc + table(idx, M_ - xsum)
       xsum = xsum + x(i)
       enc = enc - table(idx, M_ - xsum)
    end do
    enc = enc + 1
  end function encode


  ! Bijectively maps the integers from 1 to n_state to vectors in the
  ! state space.
  function decode(ix) result(x)
    integer, intent(in) :: ix
    integer, dimension(N_) :: x
    integer :: digit, xsum, i, j, idx, code, enc

    x = 0
    code = ix - 1
    xsum = 0
    do i = 1, N_
       idx = N_ - i + 1
       digit = 0

       ! Search for the largest valid digit that would result in an
       ! encoded value no larger the encoded value ix (see pg. 1088).
       do j = 1, M_ - 1 - xsum
          enc = table(idx, M_ - xsum) - table(idx, M_ - xsum - j)
          if (enc <= code) then
             digit = j
          end if
       end do

       ! Use the remainder of the encoded value to find the next digit.
       code = code - (table(idx, M_ - xsum) - table(idx, M_ - xsum - digit))
       x(i) = digit
       xsum = xsum + x(i)
    end do
  end function decode

end module encoding
