! dataset.f90 -- derived type for storing data and related procedures
!
! Copyright (C) 2008-2013 Jason R. Blevins
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

module dataset
  use osl, wp => osl_wp
  implicit none

  private

  public :: dataset_t
  public :: observation_t

  public :: dataset_init
  public :: dataset_free
  public :: dataset_size
  public :: dataset_first
  public :: dataset_next
  public :: dataset_get
  public :: dataset_append
  public :: dataset_print
  public :: observation_print

  type :: observation_t
     integer :: m                                 ! market
     real(wp) :: t                                ! time
     real(wp) :: tau                              ! time interval
     integer :: i                                 ! player
     integer :: a                                 ! action
     integer :: x                                 ! state just before event
     integer :: xp                                ! new state after event
     type(observation_t), pointer :: next => null()
  end type observation_t

  type :: dataset_t
     type(observation_t), pointer :: head => null()
     type(observation_t), pointer :: tail => null()
     integer :: size = 0
  end type dataset_t

contains

  ! Initialize a dataset.
  subroutine dataset_init(self)
    type(dataset_t), intent(inout) :: self

    if (self%size /= 0) then
       call dataset_free(self)
    end if
    nullify(self%head)
    nullify(self%tail)
    self%size = 0
  end subroutine dataset_init

  ! Print the contents of dataset in tabular form
  subroutine dataset_print(self)
    type(dataset_t), intent(inout) :: self
    type(observation_t), pointer :: obs => null()
    real(wp) :: t, tau
    integer :: m, a, player, x, xp

    print '(7a8)', 'M', 'T', 'TAU', 'A', 'PLAYER', 'X', 'XP'

    obs => dataset_first(self)
    do while (associated(obs))
       call dataset_get(obs, M=m, T=t, TAU=tau, A=a, I=player, X=x, XP=xp)
       print '(i8,1x,f7.2,1x,f7.2,1x,i7,1x,i7,1x,i7,1x,i7)', m, t, tau, a, player, x, xp
       obs => dataset_next(obs)
    end do
  end subroutine dataset_print

  ! Free the list and all allocated observations
  subroutine dataset_free(self)
    type(dataset_t), intent(inout) :: self
    type(observation_t), pointer :: current => null()
    type(observation_t), pointer :: next => null()

    if (associated(self%head)) then
       current => self%head
       do while (associated(current))
          if (associated(current%next)) then
             next => current%next
          else
             nullify(next)
          end if
          deallocate(current)
          nullify(current)
          current => next
       end do
       nullify(self%head)
       nullify(self%tail)
       self%size = 0
    end if
  end subroutine dataset_free

  ! Return the size of the given dataset.
  function dataset_size(self)
    type(dataset_t), intent(inout) :: self
    integer :: dataset_size

    dataset_size = self%size
  end function dataset_size

  ! Return the head observation.
  function dataset_first(self)
    type(dataset_t), intent(inout) :: self
    type(observation_t), pointer :: dataset_first

    nullify(dataset_first)
    if (associated(self%head)) then
       dataset_first => self%head
    end if
  end function dataset_first

  ! Return the next observation.
  function dataset_next(obs)
    type(observation_t), pointer :: obs
    type(observation_t), pointer :: dataset_next

    nullify(dataset_next)
    if (associated(obs%next)) then
       dataset_next => obs%next
    end if
  end function dataset_next

  ! Obtain the elements of the given observation.
  subroutine dataset_get(obs, m, t, tau, i, a, x, xp)
    type(observation_t), pointer :: obs
    integer, intent(out), optional :: m
    real(wp), intent(out), optional :: t
    real(wp), intent(out), optional :: tau
    integer, intent(out), optional :: i
    integer, intent(out), optional :: a
    integer, intent(out), optional :: x
    integer, intent(out), optional :: xp

    if (present(m)) m = obs%m
    if (present(t)) t = obs%t
    if (present(tau)) tau = obs%tau
    if (present(i)) i = obs%i
    if (present(a)) a = obs%a
    if (present(x)) x = obs%x
    if (present(xp)) xp = obs%xp
  end subroutine dataset_get

  ! Append the given observation to the tail.
  subroutine dataset_append(self, m, t, tau, i, a, x, xp)
    type(dataset_t), intent(inout) :: self
    integer, intent(in) :: m
    real(wp), intent(in) :: t
    real(wp), intent(in) :: tau
    integer, intent(in) :: i
    integer, intent(in) :: a
    integer, intent(in) :: x
    integer, intent(in) :: xp
    type(observation_t), pointer :: obs

    allocate(obs)
    obs%m = m
    obs%t = t
    obs%tau = tau
    obs%i = i
    obs%a = a
    obs%x = x
    obs%xp = xp

    if (.not. associated(self%tail)) then
       self%head => obs
       self%tail => obs
    else
       self%tail%next => obs
       self%tail => obs
    end if

    self%size = self%size + 1

  end subroutine dataset_append


  subroutine observation_print(self)
    type(observation_t), intent(in) :: self
    print *, 'Observation:'
    print *, 'm = ', self%m
    print *, 't = ', self%t
    print *, 'tau = ', self%tau
    print *, 'i = ', self%i
    print *, 'a = ', self%a
    print *, 'x = ', self%x
    print *, 'xp = ', self%xp
  end subroutine observation_print

end module dataset
