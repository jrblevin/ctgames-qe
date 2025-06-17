! dataset.f90 -- derived type for storing data and related procedures
!
! Copyright (C) 2008-2021 Jason R. Blevins
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
  public :: series_t
  public :: dataset_init
  public :: dataset_free
  public :: dataset_n
  public :: dataset_size
  public :: dataset_first
  public :: dataset_next
  public :: dataset_append
  public :: dataset_append_known
  public :: dataset_clear_known
  public :: dataset_print_observation

  ! Observation data type
  type :: observation_t
     integer :: m                                 ! Market index
     real(wp) :: t                                ! Time
     real(wp) :: tau                              ! Time interval
     integer :: is                                ! Structure before event
     integer :: wi                                ! Firm's w before event
     integer :: a                                 ! Firm's action
     integer :: isp                               ! New structure after event
     integer :: wip                               ! Firm's w after event

     ! Observed moves in order
     integer, dimension(:,:), allocatable :: known

     ! Next observation pointer
     type(observation_t), pointer :: next => null()
  end type observation_t

  ! Series data type
  type :: series_t
     type(observation_t), pointer :: head => null()
     type(observation_t), pointer :: tail => null()
     integer :: size
  end type series_t

  ! Datset data type
  type :: dataset_t
     integer :: n = 0
     integer :: size = 0
     type(series_t), dimension(:), pointer :: units => null()
  end type dataset_t

contains

  ! Initialize a dataset with N distinct observational units.
  subroutine dataset_init(self, n)
    type(dataset_t), intent(inout) :: self
    integer, intent(in) :: n
    type(series_t), pointer :: series
    integer :: i

    call dataset_free(self)
    self%n = n
    self%size = 0
    allocate(self%units(n))

    do i = 1, n
       series => self%units(i)
       nullify(series%head)
       nullify(series%tail)
       series%size = 0
    end do
  end subroutine dataset_init

  ! Free the dataset and all associated observations.
  subroutine dataset_free(self)
    type(dataset_t), intent(inout) :: self
    type(observation_t), pointer :: current, next
    type(series_t), pointer :: series
    integer :: i

    do i = 1, self%n
       series => self%units(i)
       current => series%head
       do while (associated(current))
          next => current%next
          call dataset_clear_known(current)
          deallocate(current)
          nullify(current)
          current => next
       end do

       nullify(series%head)
       nullify(series%tail)
    end do

    if (associated(self%units)) then
       deallocate(self%units)
    end if
    self%size = 0
    self%n = 0
  end subroutine dataset_free

  ! Return the number of observational units in the dataset.
  function dataset_n(self) result(n)
    type(dataset_t), intent(in) :: self
    integer :: n

    n = self%n
  end function dataset_n

  ! Return the number of observations in the dataset.
  function dataset_size(self) result(size)
    type(dataset_t), intent(in) :: self
    integer :: size

    size = self%size
  end function dataset_size

  ! Return the head observation.
  function dataset_first(self, i) result(first)
    type(dataset_t), intent(in) :: self
    type(observation_t), pointer :: first
    type(series_t), pointer :: series
    integer :: i

    series => self%units(i)
    if (associated(series%head)) then
       first => series%head
    else
       nullify(first)
    end if
  end function dataset_first

  ! Return the next observation.
  function dataset_next(obs) result(next)
    type(observation_t), pointer :: obs
    type(observation_t), pointer :: next

    nullify(next)
    if (associated(obs%next)) then
       next => obs%next
    end if
  end function dataset_next

  ! Append the given observation to the tail.
  subroutine dataset_append(self, i, obs)
    type(dataset_t), intent(inout) :: self
    integer, intent(in) :: i
    type(observation_t) :: obs
    type(observation_t), pointer :: newobs
    type(series_t), pointer :: series

    allocate(newobs)
    newobs = obs
    nullify(newobs%next)

    series => self%units(i)
    if (.not. associated(series%tail)) then
       series%head => newobs
       series%tail => newobs
    else
       series%tail%next => newobs
       series%tail => newobs
    end if

    series%size = series%size + 1
    self%size = self%size + 1
  end subroutine dataset_append

  ! Append a known move to the given observation.
  subroutine dataset_append_known(obs, i, a)
    type(observation_t), intent(inout) :: obs
    integer, intent(in) :: i, a

    integer, dimension(:,:), allocatable :: tmp
    integer :: len

    if (.not. allocated(obs%known)) then
       allocate(obs%known(2,1))
       len = 1
    else
       len = size(obs%known, 2)

       allocate(tmp(2,len))
       tmp = obs%known
       deallocate(obs%known)

       len = len + 1
       allocate(obs%known(2,len))
       obs%known(:,1:len-1) = tmp
       deallocate(tmp)
    end if

    obs%known(1, len) = i
    obs%known(2, len) = a
  end subroutine dataset_append_known

  ! Clear the list of known events
  subroutine dataset_clear_known(obs)
    type(observation_t), intent(inout) :: obs

    if (allocated(obs%known)) then
       deallocate(obs%known)
    end if
  end subroutine dataset_clear_known

  ! Print a particular state in human readable form.
  subroutine dataset_print_observation(obs)
    type(observation_t), intent(in) :: obs

    print '(a15,i12)', 'Market ID', obs%m
    print '(a15,f12.4)', 'Time', obs%t
    print '(a15,f12.4)', 'Tau', obs%tau
    print '(a15,i12)', 'State', obs%is
    print '(a15,i12)', 'Type', obs%wi
    print '(a15,i12)', 'Action', obs%a
    print '(a15,i12)', 'New state', obs%isp
    print '(a15,i12)', 'New type', obs%wip
    if (allocated(obs%known)) then
       print '(a15,100i4)', 'Known types', obs%known(1,:)
       print '(a15,100i4)', 'Known actions', obs%known(2,:)
    end if
  end subroutine dataset_print_observation

end module dataset
