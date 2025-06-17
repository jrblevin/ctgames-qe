! rust_data.f90 -- module for loading data from Rust (1987)
!
! Copyright (C) 2006-2024 Jason R. Blevins
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

module rust_data
  use osl, wp => osl_wp
  use dataset
  implicit none

  private
  public :: rust_data_load_groups

  ! Number of bus groups
  integer, parameter :: NGROUP = 8

  ! Number of busses for each group
  integer, parameter, dimension(NGROUP) :: NOBS = (/ 15, 4, 48, 37, 12, 10, 18, 18 /)

  ! Number periods for each group
  integer, parameter, dimension(NGROUP) :: NPER = (/ 25, 49, 70, 117, 126, 126, 126, 126 /)

  ! Filename containing data
  character(len=12), dimension(NGROUP), save  :: FILENAMES
  data FILENAMES / "g870.asc    ", &              ! Grumman 870
                   "rt50.asc    ", &              ! Chance RT50
                   "t8h203.asc  ", &              ! GMC T8H203 (1979)
                   "a530875.asc ", &              ! GMC A5308 (1975)
                   "a530874.asc ", &              ! GMC A5308 (1974)
                   "a452374.asc ", &              ! GMC A4523 (1974)
                   "a530872.asc ", &              ! GMC A5308 (1972)
                   "a452372.asc " /               ! GMC A4523 (1972)

character(len=17), dimension(NGROUP), save  :: NAMES
data NAMES / "Grumman 870      ", &
             "Chance RT50      ", &
             "GMC T8H203 (1979)", &
             "GMC A5308 (1975) ", &
             "GMC A5308 (1974) ", &
             "GMC A4523 (1974) ", &
             "GMC A5308 (1972) ", &
             "GMC A4523 (1972) " /

  ! Mileage cell size
  integer, parameter :: cell_size = 5000

  ! Dataset constants
  integer, parameter, public :: I_NATURE     = 0  ! Player ID of nature
  integer, parameter, public :: A_NATURE     = 0  ! Indicates moves by nature
  integer, parameter, public :: I_AGENT      = 1  ! Player ID of the agent
  integer, parameter, public :: A_REPLACE    = 2  ! Replacement action
  integer, parameter, public :: I_DISCRETE   = -1 ! No individual moves
  integer, parameter, public :: A_DISCRETE   = -1 ! Indicates discrete data
  integer, parameter, public :: I_END        = -2 ! Indicate end of data
  integer, parameter, public :: A_END        = -2 ! Indicate end of data

contains

  ! Load the specified bus GROUPS.
  subroutine rust_data_load_groups(data, groups)
    type(dataset_t), intent(inout) :: data
    integer, dimension(:), intent(in) :: groups
    integer :: group_idx, i, j, nt, count
    integer, dimension(:), allocatable :: data_x_lag, data_x, data_a

    ! Initialize the dataset
    call dataset_init(data)

    ! Process each group listed in 'groups'
    do i = 1, size(groups)
      group_idx = groups(i)
      nt = NOBS(group_idx) * (NPER(group_idx) - 1)
      count = 0

      allocate(data_x_lag(nt), data_x(nt), data_a(nt))

      ! Load data for the current group
      call rust_data_load_group(group_idx, data_x_lag, data_x, data_a)

      ! Append each observation to the dataset
      do j = 1, nt
        call dataset_append(data, m=group_idx, t=real(j, wp), tau=1.0_wp, i=I_DISCRETE, a=A_DISCRETE, x=data_x_lag(j), xp=data_x(j))
        count = count + 1
      end do

      deallocate(data_x_lag, data_x, data_a)
      print '("* Group ", i2, ": Loaded ", i6, " observations. ", i3, " buses with ", i3, " months per bus")', &
          group_idx, count, NOBS(group_idx), NPER(group_idx) - 1
    end do
    print '("* Loaded ", i6, " observations in total.")', dataset_size(data)
  end subroutine rust_data_load_groups

  ! Load data for a particular bus GROUP and return the choices and
  ! states in DATA_A and DATA_X respectively.  The data files are
  ! assumed to reside in the a subdirectory called data of the current
  ! working directory at runtime.
  subroutine rust_data_load_group(group, data_x_lag, data_x, data_a)
    integer, intent(in) :: group
    integer, dimension(:), intent(out) :: data_x_lag
    integer, dimension(:), intent(out) :: data_x
    integer, dimension(:), intent(out) :: data_a
    character(len=20) :: filename
    integer :: i, j, off, nt
    integer :: bus_number, purchase_month, purchase_year
    integer :: month0, month1, month2, year0, year1, year2
    integer :: miles1, miles2
    integer, dimension(:), allocatable :: mileage
    integer, dimension(:), allocatable :: actions

    nt = NPER(group)
    off = 0

    allocate(mileage(nt), actions(nt))

    ! Open the file for reading
    filename = "data/" // FILENAMES(group)
    open(15, file=filename)

    ! Read the data for each bus
    do i = 1, NOBS(group)
       ! Read the header
       read(15, *) bus_number, purchase_month, purchase_year, month1, &
            year1, miles1, month2, year2, miles2, month0, year0

       ! Read the raw mileage data (all lines at once)
       read(15,*) mileage

       ! Calculate choice vector and discretize mileage
       actions = 0
       if (miles1 > 0.0_wp) where (mileage >= miles1) actions = 1
       if (miles2 > 0.0_wp) where (mileage >= miles2) actions = 2
       where (actions == 1) mileage = mileage - miles1
       where (actions == 2) mileage = mileage - miles2
       actions(2:nt) = actions(2:nt) - actions(1:nt-1)
       actions(1) = 0
       actions = actions + 1

       ! Scale and discretize mileage
       mileage = (mileage / cell_size) + 1

       ! Store and prepare for another loop
       data_x_lag(off+1:off+nt-1) = mileage(1:nt-1)
       data_x(off+1:off+nt-1) = mileage(2:nt)
       data_a(off+1:off+nt-1) = actions(2:nt)

       off = off + nt - 1

    end do

    close(15)
    deallocate(mileage, actions)
  end subroutine rust_data_load_group

end module rust_data
