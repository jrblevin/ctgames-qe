! sparse.f90 -- sparse matrix type and related procedures
!
! Copyright (C) 2021 Jason R. Blevins
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

module sparse
  use osl, wp => osl_wp
  use stdlib_experimental_sparse
  implicit none
  private

  ! Types
  integer, parameter :: SPMATRIX_TYPE_NONE = 0
  integer, parameter :: SPMATRIX_TYPE_COO = 1
  integer, parameter :: SPMATRIX_TYPE_CSR = 2
  integer, parameter :: SPMATRIX_TYPE_CSC = 3

  ! Sparse matrix type
  type :: spmatrix_t
     real(wp), dimension(:), allocatable :: values
     integer, dimension(:), allocatable :: row_idx
     integer, dimension(:), allocatable :: col_idx
     integer, dimension(:), allocatable :: row_ptr
     integer, dimension(:), allocatable :: col_ptr
     integer :: rows = 0
     integer :: cols = 0
     integer :: nzmax = 0
     integer :: nz = 0
     integer :: type = SPMATRIX_TYPE_NONE
  end type spmatrix_t

  public :: spmatrix_t
  public :: SPMATRIX_TYPE_NONE
  public :: SPMATRIX_TYPE_COO
  public :: SPMATRIX_TYPE_CSR
  public :: SPMATRIX_TYPE_CSC
  public :: spmatrix_init
  public :: spmatrix_free
  public :: spmatrix_get
  public :: spmatrix_set
  public :: spmatrix_nz
  public :: spmatrix_coo_to_csr
  public :: spmatrix_coo_to_csc
  public :: spmatrix_csr_to_coo
  public :: spmatrix_csr_matvec
  public :: spmatrix_csc_vecmat
  public :: spmatrix_print_info

  public :: expvunif_sparse

contains

  subroutine spmatrix_init(m, rows, cols, nzmax, type)
    type(spmatrix_t), intent(inout) :: m
    integer, intent(in) :: rows
    integer, intent(in) :: cols
    integer, intent(in) :: nzmax
    integer, intent(in) :: type

    m%rows = rows
    m%cols = cols
    m%nzmax = nzmax
    m%nz = 0
    m%type = type
    allocate(m%values(nzmax))
    if (type == SPMATRIX_TYPE_COO) then
       allocate(m%col_idx(nzmax))
       allocate(m%row_idx(nzmax))
    else if (type == SPMATRIX_TYPE_CSR) then
       allocate(m%col_idx(nzmax))
       allocate(m%row_ptr(rows+1))
    else if (type == SPMATRIX_TYPE_CSC) then
       allocate(m%row_idx(nzmax))
       allocate(m%col_ptr(cols+1))
    end if
  end subroutine spmatrix_init

  subroutine spmatrix_free(m)
    type(spmatrix_t), intent(inout) :: m
    if (allocated(m%values)) then
       deallocate(m%values)
    end if
    if (allocated(m%row_idx)) then
       deallocate(m%row_idx)
    end if
    if (allocated(m%col_idx)) then
       deallocate(m%col_idx)
    end if
    if (allocated(m%row_ptr)) then
       deallocate(m%row_ptr)
    end if
    if (allocated(m%col_ptr)) then
       deallocate(m%col_ptr)
    end if
    m%rows = 0
    m%cols = 0
    m%nzmax = 0
    m%nz = 0
    m%type = SPMATRIX_TYPE_NONE
  end subroutine spmatrix_free

  ! subroutine csc_matrix_init(m, nnz, cols)
  !   type(csc_matrix_t), intent(inout) :: m
  !   integer, intent(in) :: nnz
  !   integer, intent(in) :: rows

  !   allocate(m%values(nnz))
  !   allocate(m%rows(nnz))
  !   allocate(m%col_ptr(cols))
  ! end subroutine csc_matrix_init

  function spmatrix_get(m, i, j) result(x)
    type(spmatrix_t), intent(inout) :: m
    integer, intent(in) :: i, j
    real(wp) :: x

    x = 0.0_wp
    if (m%type == SPMATRIX_TYPE_COO) then
       stop 'ERROR: SPMATRIX_GET NOT IMPLEMENTED FOR COO MATRICES'
    else if (m%type == SPMATRIX_TYPE_CSR) then
       x = csr_getvalue(AP=m%row_ptr, AJ=m%col_idx, AX=m%values, I=i, J=j)
    else if (m%type == SPMATRIX_TYPE_CSC) then
       x = csc_getvalue(AP=m%col_ptr, AI=m%row_idx, Ax=m%values, I=i, J=j)
    else if (m%type == SPMATRIX_TYPE_NONE) then
       stop 'ERROR: SPMATRIX_GET CALLED WITH UNALLOCATED MATRIX'
    end if
  end function spmatrix_get

  subroutine spmatrix_set(m, i, j, x)
    type(spmatrix_t), intent(inout) :: m
    integer, intent(in) :: i, j
    real(wp), intent(in) :: x

    if (m%type == SPMATRIX_TYPE_COO) then
       m%nz = m%nz+1
       m%row_idx(m%nz) = i
       m%col_idx(m%nz) = j
       m%values(m%nz) = x
    else
       stop 'ERROR: SPMATRIX_SET CURRENTLY ONLY IMPLEMENTED FOR COO MATRICES'
    end if
  end subroutine spmatrix_set

  function spmatrix_nz(m) result(nz)
    type(spmatrix_t), intent(in) :: m
    integer :: nz
    nz = m%nz
  end function spmatrix_nz

  subroutine spmatrix_coo_to_csr(m_coo, m_csr)
    type(spmatrix_t), intent(in) :: m_coo
    type(spmatrix_t), intent(inout) :: m_csr

    if (m_csr%nz > 0) then
       call spmatrix_free(m_csr)
    end if
    call spmatrix_init(m_csr, &
         ROWS=m_coo%rows, &
         COLS=m_coo%cols, &
         NZMAX=m_coo%nz, &
         TYPE=SPMATRIX_TYPE_CSR)

    call coo2csr_canonical(AI=m_coo%row_idx(1:m_coo%nz), &
                           AJ=m_coo%col_idx(1:m_coo%nz), &
                           AX=m_coo%values(1:m_coo%nz), &
                           BP=m_csr%row_ptr, &
                           BJ=m_csr%col_idx, &
                           BX=m_csr%values)

    ! Update number of nonzeros
    m_csr%nz = m_coo%nz
  end subroutine spmatrix_coo_to_csr

  subroutine spmatrix_coo_to_csc(m_coo, m_csc)
    type(spmatrix_t), intent(in) :: m_coo
    type(spmatrix_t), intent(inout) :: m_csc

    if (m_csc%nz > 0) then
       call spmatrix_free(m_csc)
    end if
    call spmatrix_init(m_csc, &
         ROWS=m_coo%rows, &
         COLS=m_coo%cols, &
         NZMAX=m_coo%nz, &
         TYPE=SPMATRIX_TYPE_CSC)

    call coo2csc_canonical(AI=m_coo%row_idx(1:m_coo%nz), &
                           AJ=m_coo%col_idx(1:m_coo%nz), &
                           AX=m_coo%values(1:m_coo%nz), &
                           BP=m_csc%col_ptr, &
                           BI=m_csc%row_idx, &
                           BX=m_csc%values)
  end subroutine spmatrix_coo_to_csc

  subroutine spmatrix_csr_to_coo(m_csr, m_coo)
    type(spmatrix_t), intent(in) :: m_csr
    type(spmatrix_t), intent(out) :: m_coo
    integer :: i, j, p
    real(wp) :: x

    if (m_coo%nz > 0) then
       call spmatrix_free(m_coo)
    end if
    call spmatrix_init(m_coo, &
         ROWS=m_csr%rows, &
         COLS=m_csr%cols, &
         NZMAX=m_csr%nz, &
         TYPE=SPMATRIX_TYPE_COO)

    do i = 1, size(m_csr%row_ptr)-1
       do p = m_csr%row_ptr(i), m_csr%row_ptr(i+1)-1
          j = m_csr%col_idx(p)
          x = m_csr%values(p)
          call spmatrix_set(m_coo, i, j, x)
       end do
    end do
  end subroutine spmatrix_csr_to_coo

  function spmatrix_csr_matvec(m, x) result(y)
    type(spmatrix_t), intent(in) :: m
    real(wp), dimension(:) :: x
    real(wp), dimension(m%rows) :: y
    integer :: maxrow

    ! Dimension of csr_matvec result
    maxrow = size(m%row_ptr,1) - 1

    y(1:maxrow) = csr_matvec(m%row_ptr, m%col_idx, m%values, x)
    y(maxrow+1:) = 0.0_wp
  end function spmatrix_csr_matvec

  function spmatrix_csc_vecmat(x, m) result(y)
    real(wp), dimension(:) :: x
    type(spmatrix_t), intent(in) :: m
    real(wp), dimension(m%cols) :: y
    integer :: maxcol

    ! Dimension of csr_matvec result
    maxcol = size(m%col_ptr,1) - 1

    y(1:maxcol) = csc_vecmat(x, m%col_ptr, m%row_idx, m%values)
    y(maxcol+1:) = 0.0_wp
  end function spmatrix_csc_vecmat

  subroutine spmatrix_print_info(m)
    type(spmatrix_t), intent(in) :: m

    call osl_print('spmatrix_t:')
    select case (m%type)
    case (SPMATRIX_TYPE_NONE)
       call osl_print('    type: None')
    case (SPMATRIX_TYPE_COO)
       call osl_print('    type: COO')
    case (SPMATRIX_TYPE_CSR)
       call osl_print('    type: CSR')
    case (SPMATRIX_TYPE_CSC)
       call osl_print('    type: CSC')
    end select
    call osl_print('    rows:', m%rows)
    call osl_print('    cols:', m%cols)
    call osl_print('    nz:', m%nz)
    call osl_print('    nzmax:', m%nzmax)

    call osl_print('    size(values):', size(m%values,1))
    if (allocated(m%row_idx)) then
       call osl_print('    size(row_idx):', size(m%row_idx,1))
    end if
    if (allocated(m%col_idx)) then
       call osl_print('    size(col_idx):', size(m%col_idx,1))
    end if
    if (allocated(m%row_ptr)) then
       call osl_print('    size(row_ptr):', size(m%row_ptr,1))
    end if
    if (allocated(m%col_ptr)) then
       call osl_print('    size(col_ptr):', size(m%col_ptr,1))
    end if
  end subroutine spmatrix_print_info

  !===========================================================================!

  ! Calculate v'*exp(t*H) for an N-by-N matrix H via the
  ! uniformization algorithm of Sherlock (2020). Sparse matrix
  ! version.
  !
  ! Sherlock, Chris (2020).
  ! Direct statistical inference for finite Markov jump processes via
  ! the matrix exponential.
  ! Working Paper, Lancaster University.
  ! https://arxiv.org/abs/1809.07110
  function expvunif_sparse(v, t, rho, Hp, Hi, Hx, m) result(mu)
    real(wp), dimension(:), intent(in) :: v
    real(wp), intent(in) :: t
    real(wp), intent(in) :: rho
    integer, intent(in) :: Hp(:), Hi(:)
    real(wp), intent(in) :: Hx(:)
    real(wp), dimension(size(v,1)) :: vin, vpro, vsum, y, mu
    integer, intent(in), optional :: m

    real(wp), parameter :: BIG = 1.0e100_wp
    real(wp), parameter :: eps = 1.0e-13_wp

    real(wp) :: b, c, fac
    integer :: m_eps, j

    if (present(m)) then
       m_eps = m
    else
       m_eps = osl_poisson_inv(1.0_wp - eps, rho)
    end if

    b = maxval(abs(v))
    vin = v
    c = 0.0_wp
    if (b > BIG) then
       vin = vin / b
       c = c + log(b)
       b = 1.0_wp
    end if

    vpro = vin
    vsum = vpro
    fac = 1.0_wp

    do j = 1, m_eps
       y = csc_vecmat(vpro, Hp, Hi, Hx)
       vpro = (t*y + rho*vpro) / fac

       b = b * rho / fac
       vsum = vsum + vpro

       if (b > BIG) then
          vpro = vpro / b
          vsum = vsum / b
          c = c + log(b)
          b = 1.0_wp
       end if
       fac = fac + 1.0_wp
    end do

    mu = exp(c - rho) * vsum
  end function expvunif_sparse

end module sparse
