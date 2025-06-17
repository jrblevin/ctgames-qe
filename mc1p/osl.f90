! osl.f90
!
! Copyright (C) 2005-2025 Jason R. Blevins
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

module osl
  implicit none

  public

  !=== Precision Parameters ==================================================!

  ! Real precision constants
  integer, parameter :: osl_dp = selected_real_kind(15, 307)  ! double
  integer, parameter :: osl_wp = osl_dp                       ! working

  ! Integer precision constants
  integer, parameter :: osl_int32 = selected_int_kind(9)      ! 32 bit

  ! internal-use private precision constants
  integer, parameter, private :: dp = osl_dp
  integer, parameter, private :: wp = osl_wp

  ! precision constants
  real(dp), parameter :: osl_eps_dp = epsilon(1.0_dp)
  real(wp), parameter :: osl_eps_wp = epsilon(1.0_wp)

  ! A context mold for generic interfaces.
  integer, dimension(1), parameter :: osl_ctx = (/ 0 /)

  !=== Mathematical Constants ================================================!

  real(wp), parameter :: OSL_EULER     = 0.5772156649015328606065120900824024_wp  ! gamma

  !=== Formatted Output ======================================================!

  ! Default edit descriptors
  character(len=*), parameter :: osl_ed_i   = 'i10'
  character(len=*), parameter :: osl_ed_wp  = 'g17.5'
  character(len=*), parameter :: osl_ed_l   = 'l1'
  integer, parameter, private :: fmt_len    = 80

  ! Procedures for printing informative messages
  interface osl_print
     module procedure osl_print_str
     module procedure osl_print_swp
     module procedure osl_print_vwp
     module procedure osl_print_mwp
     module procedure osl_print_si
     module procedure osl_print_vi
     module procedure osl_print_mi
     module procedure osl_print_sl
     module procedure osl_print_vl
     module procedure osl_print_ml
  end interface osl_print

  !=== Norms =================================================================!

  interface osl_norm_2
     module procedure osl_norm_2_vwp
  end interface osl_norm_2

  !=== Array Operations ======================================================!

  ! Swap
  interface osl_swap
     module procedure osl_swap_scalar_wp
  end interface osl_swap

  !=== Pseudo-Random Number Generation =======================================!

  ! Pseudo-random number generator
  public :: osl_rng_t
  public :: osl_rng_seed
  public :: osl_rng_get
  public :: osl_rng_uniform

  ! Continuous distributions
  public :: osl_exponential_rnd

  ! Mersenne Twister parameters
  integer, parameter :: OSL_MT_N = 624
  integer, parameter :: OSL_MT_M = 397
  integer(osl_int32), parameter :: DEFAULT_SEED = 5489_osl_int32
  integer(osl_int32), parameter :: MATRIX_A = -1727483681       ! 0x9908b0df
  integer(osl_int32), parameter :: UPPER_MASK = -2147483647 - 1 ! 0x80000000
  integer(osl_int32), parameter :: LOWER_MASK = 2147483647      ! 0x7fffffff
  integer(osl_int32), parameter :: TMASKB = -1658038656         ! 0x9d2c5680
  integer(osl_int32), parameter :: TMASKC = -272236544          ! 0xefc60000

  integer(osl_int32), dimension(0:1), parameter :: mag01 = (/ 0_osl_int32, MATRIX_A /)

  type :: osl_rng_t
     private
     integer(osl_int32), dimension(OSL_MT_N) :: mt
     integer(osl_int32) :: mti = OSL_MT_N + 2
  end type osl_rng_t

  !=== Sorting ===============================================================!

  public :: osl_quicksort

  interface osl_quicksort
     module procedure osl_quicksort_vwp
  end interface osl_quicksort

  !=== Statistics ============================================================!

  interface osl_mean
     module procedure osl_mean_vwp
     module procedure osl_mean_mwp
  end interface

  interface osl_median
     module procedure osl_median_vwp
  end interface

  interface osl_cov
     module procedure osl_cov_vwp
     module procedure osl_cov_mwp
  end interface

  interface osl_var
     module procedure osl_var_vwp
     module procedure osl_var_mwp
  end interface

contains

  !=== Formatted output =====================================================!

  !--- Informative messages -------------------------------------------------!

  ! Print MSG followed by a REAL(WP) matrix A.  If MSG is the empty string,
  ! then no blank line will be printed.  If FMT is present, it will be used
  ! to format each element of A.  Otherwise, OSL_ED_WP, the default edit
  ! descriptor for REAL(WP) values, will be used.
  subroutine osl_print_mwp(msg, A, fmt)
    character(len=*), intent(in) :: msg
    real(wp), dimension(:,:), intent(in) :: A
    character(len=*), intent(in), optional :: fmt
    character(len=fmt_len) :: row_fmt
    integer :: i

    if (len(msg) > 0) then
       print '(a)', msg
    end if

    if (present(fmt)) then
       write(row_fmt, '("(",i0,"(",a,"))")') size(A,2), fmt
    else
       write(row_fmt, '("(",i0,"(",a,"))")') size(A,2), osl_ed_wp
    end if

    do i = 1, size(A,1)
       print row_fmt, A(i,:)
    end do
  end subroutine osl_print_mwp

  ! Print MSG followed by an INTEGER matrix A.  If MSG is the empty string,
  ! then no blank line will be printed.  If FMT is present, it will be used
  ! to format each element of A.  Otherwise, OSL_ED_I, the default edit
  ! descriptor for INTEGER values, will be used.
  subroutine osl_print_mi(msg, A, fmt)
    character(len=*), intent(in) :: msg
    integer, dimension(:,:), intent(in) :: A
    character(len=*), intent(in), optional :: fmt
    character(len=fmt_len) :: row_fmt
    integer :: i

    if (len(msg) > 0) then
       print '(a)', msg
    end if

    if (present(fmt)) then
       write(row_fmt, '("(",i0,"(",a,"))")') size(A,2), fmt
    else
       write(row_fmt, '("(",i0,"(",a,"))")') size(A,2), osl_ed_i
    end if

    do i = 1, size(A,1)
       print row_fmt, A(i,:)
    end do
  end subroutine osl_print_mi

  ! Print MSG followed by a LOGICAL matrix A.  If MSG is the empty string,
  ! then no blank line will be printed.  If FMT is present, it will be used
  ! to format each element of A.  Otherwise, OSL_ED_L, the default edit
  ! descriptor for LOGICAL values, will be used.
  subroutine osl_print_ml(msg, A, fmt)
    character(len=*), intent(in) :: msg
    logical, dimension(:,:), intent(in) :: A
    character(len=*), intent(in), optional :: fmt
    character(len=fmt_len) :: row_fmt
    integer :: i

    if (len(msg) > 0) then
       print '(a)', msg
    end if

    if (present(fmt)) then
       write(row_fmt, '("(",i0,"(",a,"))")') size(A,2), fmt
    else
       write(row_fmt, '("(",i0,"(",a,"))")') size(A,2), osl_ed_i
    end if

    do i = 1, size(A,1)
       print row_fmt, A(i,:)
    end do
  end subroutine osl_print_ml

  ! Print MSG followed by a REAL(WP) vector VEC.  At most ten elements
  ! will be printed per line.  If MSG is the empty string, then no
  ! blank line will be printed.  If FMT is present, it will be used to
  ! format each element of VEC.  Otherwise, OSL_ED_WP, the default edit
  ! descriptor for REAL(WP) values, will be used.
  subroutine osl_print_vwp(msg, vec, fmt)
    character(len=*), intent(in) :: msg
    real(wp), dimension(:), intent(in) :: vec
    character(len=*), intent(in), optional :: fmt
    character(len=fmt_len) :: fmt_str
    integer :: i, ll, ncols, lines

    if (len(msg) > 0) then
       print '(a)', msg
    end if

    if (present(fmt)) then
       fmt_str = '(10' // fmt // ')'
    else
       fmt_str = '(10' // osl_ed_wp // ')'
    end if

    ncols = size(vec, 1)
    ll = 0
    if (ncols > 10) then
       lines = ncols / 10
       do i = 1, lines
          ll = 10 * (i - 1)
          print fmt_str, vec(1+ll:10+ll)
       end do
       print fmt_str, vec(11+ll:ncols)
    else
       print fmt_str, vec
    end if
  end subroutine osl_print_vwp

  ! Print MSG followed by an INTEGER vector VEC.  At most ten elements
  ! will be printed per line.  If MSG is the empty string, then no
  ! blank line will be printed.  If FMT is present, it will be used to
  ! format each element of VEC.  Otherwise, OSL_ED_I, the default edit
  ! descriptor for INTEGER values, will be used.
  subroutine osl_print_vi(msg, vec, fmt)
    character(len=*), intent(in) :: msg
    integer, dimension(:), intent(in) :: vec
    character(len=*), intent(in), optional :: fmt
    character(len=fmt_len) :: fmt_str
    integer :: i, ll, ncols, lines

    if (len(msg) > 0) then
       print '(a)', msg
    end if

    if (present(fmt)) then
       fmt_str = '(10' // fmt // ')'
    else
       fmt_str = '(10' // osl_ed_i // ')'
    end if

    ncols = size(vec, 1)
    ll = 0
    if (ncols > 10) then
       lines = ncols / 10
       do i = 1, lines
          ll = 10 * (i - 1)
          print fmt_str, vec(1+ll:10+ll)
       end do
       print fmt_str, vec(11+ll:ncols)
    else
       print fmt_str, vec
    end if
  end subroutine osl_print_vi

  ! Print MSG followed by a LOGICAL vector VEC.  At most ten elements
  ! will be printed per line.  If MSG is the empty string, then no
  ! blank line will be printed.  If FMT is present, it will be used to
  ! format each element of VEC.  Otherwise, OSL_ED_L, the default edit
  ! descriptor for LOGICAL values, will be used.
  subroutine osl_print_vl(msg, vec, fmt)
    character(len=*), intent(in) :: msg
    logical, dimension(:), intent(in) :: vec
    character(len=*), intent(in), optional :: fmt
    character(len=fmt_len) :: fmt_str
    integer :: i, ll, ncols, lines

    if (len(msg) > 0) then
       print '(a)', msg
    end if

    if (present(fmt)) then
       fmt_str = '(10' // fmt // ')'
    else
       fmt_str = '(10' // osl_ed_l // ')'
    end if

    ncols = size(vec, 1)
    ll = 0
    if (ncols > 10) then
       lines = ncols / 10
       do i = 1, lines
          ll = 10 * (i - 1)
          print fmt_str, vec(1+ll:10+ll)
       end do
       print fmt_str, vec(11+ll:ncols)
    else
       print fmt_str, vec
    end if
  end subroutine osl_print_vl

  ! Print MSG followed by the REAL(WP) scalar X.  If FMT is present,
  ! it will be used to format X.  Otherwise, OSL_ED_WP, the default
  ! edit descriptor for REAL(WP) values, will be used.
  subroutine osl_print_swp(msg, x, fmt)
    character(len=*), intent(in) :: msg
    real(wp), intent(in) :: x
    character(len=*), intent(in), optional :: fmt

    if (present(fmt)) then
       print '(a,' // fmt // ')', msg, x
    else
       print '(a,' // osl_ed_wp // ')', msg, x
    end if
  end subroutine osl_print_swp

  ! Print MSG followed by the INTEGER scalar X.  If FMT is present,
  ! it will be used to format X.  Otherwise, OSL_ED_I, the default
  ! edit descriptor for INTEGER values, will be used.
  subroutine osl_print_si(msg, x, fmt)
    character(len=*), intent(in) :: msg
    integer, intent(in) :: x
    character(len=*), intent(in), optional :: fmt

    if (present(fmt)) then
       print '(a,' // fmt // ')', msg, x
    else
       print '(a,' // osl_ed_i // ')', msg, x
    end if
  end subroutine osl_print_si

  ! Print MSG followed by the LOGICAL scalar X.  If FMT is present,
  ! it will be used to format X.  Otherwise, OSL_ED_L, the default
  ! edit descriptor for LOGICAL values, will be used.
  subroutine osl_print_sl(msg, x, fmt)
    character(len=*), intent(in) :: msg
    logical, intent(in) :: x
    character(len=*), intent(in), optional :: fmt

    if (present(fmt)) then
       print '(a,' // fmt // ')', msg, x
    else
       print '(a,' // osl_ed_l // ')', msg, x
    end if
  end subroutine osl_print_sl

  ! Prints a single string.
  subroutine osl_print_str(msg)
    character(len=*), intent(in) :: msg
    print '(a)', msg
  end subroutine osl_print_str

  !--- Headings --------------------------------------------------------------!

  ! Print a heading consisting of TEXT underlined by equals signs.
  subroutine osl_print_heading(text)
    character(len=*), intent(in) :: text
    print '(/,a,/,a,/)', text, repeat("=", len_trim(text))
  end subroutine osl_print_heading

  ! Print a subheading consisting of TEXT underlined by hyphens.
  subroutine osl_print_subheading(text)
    character(len=*), intent(in) :: text
    print '(/,a,/,a,/)', text, repeat("-", len_trim(text))
  end subroutine osl_print_subheading

  !=== Norms =================================================================!

  !--- L_2 Norm --------------------------------------------------------------!

  ! Returns the Euclidean norm of a vector defined as
  ! SQRT(DOT_PRODUCT(x,x)).  This function is based on the BLAS
  ! Fortran 77 implementation of DNRM2 which attempts to prevent
  ! overflow.
  function osl_norm_2_vwp(x) result(norm)
    real(wp), dimension(:), intent(in) :: x
    real(wp) :: norm
    real(wp) :: absxi, scale, ssq
    real(wp), parameter :: ONE = 1.0_wp, ZERO = 0.0_wp
    integer :: i, n

    n = size(x)
    if (n == 1) then
       norm = abs(x(1))
    else
       scale = ZERO
       ssq = ONE
       do i = 1, n
          if (x(i) /= ZERO) THEN
             absxi = abs(x(i))
             if (scale < absxi) then
                ssq = ONE + ssq * (scale / absxi)**2
                scale = absxi
             else
                ssq = ssq + (absxi / scale)**2
             end if
          end if
       end do
       norm = scale * sqrt(ssq)
    end if
  end function osl_norm_2_vwp

  !=== Miscellaneous Mathematical Functions ==================================!

  !--- Swap ------------------------------------------------------------------!

  ! Swap the values of two REAL(WP) scalars.
  subroutine osl_swap_scalar_wp(a, b)
    real(wp), intent(inout) :: a, b
    real(wp) :: temp
    temp = a
    a = b
    b = temp
  end subroutine osl_swap_scalar_wp

  !--- Vector-vector operations ----------------------------------------------!

  ! Return the outer product of the vectors A and B.
  function osl_outer_product(a, b) result(C)
    real(wp), dimension(:), intent(in) :: a
    real(wp), dimension(:), intent(in) :: b
    real(wp), dimension(size(a), size(b)) :: C
    C = spread(a, DIM=2, NCOPIES=size(b)) * spread(b, DIM=1, NCOPIES=size(a))
  end function osl_outer_product

  !=== Mersenne Twister ======================================================!

  subroutine osl_rng_seed(self, s)
    type(osl_rng_t), intent(inout) :: self
    integer(osl_int32), intent(in) :: s
    integer :: mti

    self%mt(1) = s
    do mti = 2, OSL_MT_N
       self%mt(mti) = 1812433253 &
            * ieor(self%mt(mti - 1), ishft(self%mt(mti - 1), -30)) &
            + mti - 1
    end do
    self%mti = OSL_MT_N + 1
  end subroutine osl_rng_seed

  function osl_rng_get(self) result(y)
    type(osl_rng_t), intent(inout) :: self
    integer(osl_int32) :: y
    integer :: kk

    ! generate N words at once
    if (self%mti >= OSL_MT_N + 1) then
        if (self%mti == OSL_MT_N + 2) then
           ! osl_rng_seed() has not been called, use default seed
           call osl_rng_seed(self, DEFAULT_SEED)
        end if
        do kk = 1, OSL_MT_N - OSL_MT_M
           y = ior(iand(self%mt(kk), UPPER_MASK), &
                   iand(self%mt(kk+1), LOWER_MASK))
           self%mt(kk) = ieor(ieor(self%mt(kk + OSL_MT_M), ishft(y, -1)), &
                              mag01(iand(y, 1_osl_int32)))
        end do
        do kk = OSL_MT_N - OSL_MT_M + 1, OSL_MT_N - 1
           y = ior(iand(self%mt(kk), UPPER_MASK), &
                   iand(self%mt(kk+1), LOWER_MASK))
           self%mt(kk) = ieor(ieor(self%mt(kk + (OSL_MT_M - OSL_MT_N)), ishft(y, -1)), &
                              mag01(iand(y, 1_osl_int32)))
        end do
        y = ior(iand(self%mt(OSL_MT_N), UPPER_MASK), iand(self%mt(1), LOWER_MASK))
        self%mt(OSL_MT_N) = ieor(ieor(self%mt(OSL_MT_M), ishft(y, -1)), mag01(iand(y, 1_osl_int32)))
        self%mti = 1
    end if

    y = self%mt(self%mti)

    ! Tempering
    y = ieor(y, ishft(y, -11))
    y = ieor(y, iand(ishft(y, 7), TMASKB))
    y = ieor(y, iand(ishft(y, 15), TMASKC))
    y = ieor(y, ishft(y, -18))

    self%mti = self%mti + 1
  end function osl_rng_get

  !--------------------------------------------------------------------------!
  ! Generic RNG Interface                                                    !
  !--------------------------------------------------------------------------!

  ! Returns a uniformly distributed real number on [0,1).
  function osl_rng_uniform(self) result(u)
    type(osl_rng_t), intent(inout) :: self
    real(wp) :: u
    ! 1 / 2^32 = 2.3283064365e-10
    u = 0.5_wp + real(osl_rng_get(self), wp) * 2.3283064365e-10_wp
  end function osl_rng_uniform

  !--------------------------------------------------------------------------!
  ! Exponential Distribution                                                 !
  !--------------------------------------------------------------------------!

  ! Draw from the exponential distribution with rate parameter LAMBDA
  ! using the inverse cdf transformation.  Note that LAMBDA must be
  ! strictly positive, but this condition is not checked for efficiency.
  function osl_exponential_rnd(rng, lambda) result(rnd)
    type(osl_rng_t), intent(inout) :: rng
    real(wp), intent(in) :: lambda
    real(wp) :: rnd

    rnd = -log(osl_rng_uniform(rng)) / lambda
  end function osl_exponential_rnd

  !---------------------------------------------------------------------------!
  ! Quicksort                                                                 !
  !---------------------------------------------------------------------------!

  ! Sort an array of real numbers in ascending order.  This is an
  ! implementation of the Quicksort algorithm of C.A.R. Hoare (1961).
  !
  ! *References*
  !
  !  - Hoare, C. A. R. (1961). Quicksort: Algorithm 64.
  !    _Communications of the ACM_ 4, 321-322.
  subroutine osl_quicksort_vwp(v)
    real(wp), dimension(:), intent(inout) :: v
    integer, parameter :: m = 7
    integer :: n_stack = 50
    integer, dimension(:), allocatable :: s, s_temp
    integer :: i, ir, j, k, jstack, l, n
    real(wp) :: a

    ! initialize stack
    allocate(s(n_stack))

    n = size(v)

    jstack = 0
    l = 1
    ir = n

    do
       if (ir - l < m) then
          do j = l + 1, ir
             a = v(j)
             do i = j - 1 , 1, -1
                if (v(i) <= a) exit
                v(i+1) = v(i)
             end do
             v(i+1) = a
          end do

          if (jstack == 0) return

          ir = s(jstack)
          l = s(jstack-1)
          jstack = jstack - 2

       else

          k = (l+ir)/2
          call osl_swap(v(k),v(l+1))
          if (v(l+1) > v(ir)) call osl_swap(v(l+1), v(ir))
          if (v(l) > v(ir)) call osl_swap(v(l), v(ir))
          if (v(l+1) > v(l)) call osl_swap(v(l+1), v(l))

          i = l + 1
          j = ir
          a = v(l)

          do
             do
                i = i + 1
                if (v(i) >= a) exit
             end do

             do
                j = j - 1
                if (v(j) <= a) exit
             end do

             if (j < i) exit
             call osl_swap(v(i), v(j))
          end do

          v(l) = v(j)
          v(j) = a
          jstack = jstack + 2

          ! double the stack size if necessary
          if (jstack > n_stack) then
             allocate(s_temp(n_stack))
             s_temp = s
             deallocate(s)
             n_stack = 2 * n_stack
             allocate(s(n_stack))
             s(1:size(s_temp,1)) = s_temp
             deallocate(s_temp)
          end if

          if (ir - i + 1 > j - l) then
             s(jstack) = ir
             s(jstack-1) = i
             ir = j-1
          else
             s(jstack) = j-1
             s(jstack-1) = l
             l = i
          end if
       end if
    end do

    ! free stack
    deallocate(s)
  end subroutine osl_quicksort_vwp

  !---------------------------------------------------------------------------!
  ! Sample Statistics                                                         !
  !---------------------------------------------------------------------------!

  !------------------------------------!
  ! Sample mean                        !
  !------------------------------------!

  ! Return the mean of a vector of real numbers V.
  function osl_mean_vwp(v) result(mean)
    real(wp), dimension(:) :: v
    real(wp) :: mean
    mean = sum(v) / real(size(v), wp)
  end function osl_mean_vwp

  ! Return the mean of the columns of a real matrix A.
  function osl_mean_mwp(a) result(mean)
    real(wp), dimension(:,:) :: A
    real(wp), dimension(size(A,1)) :: mean
    mean = sum(A, DIM=2) / real(size(A,2), wp)
  end function osl_mean_mwp

  !------------------------------------!
  ! Sample median                      !
  !------------------------------------!

  ! Return the median of a vector of real numbers V.
  function osl_median_vwp(v) result(median)
    real(wp), dimension(:) :: v
    real(wp), dimension(size(v)) :: vsort
    real(wp) :: median
    integer :: n, imed

    n = size(v, 1)
    vsort = v
    call osl_quicksort(vsort)

    if (mod(n, 2) == 0) then
       imed = n / 2
       median = 0.5_wp * (vsort(imed) + vsort(imed+1))
    else
       imed = n/2 + 1
       median = vsort(imed)
    end if
  end function osl_median_vwp

  !------------------------------------!
  ! Sample covariance                  !
  !------------------------------------!

  ! Return the sample covariance between two vectors X and Y.  If Y
  ! is omitted, return the variance of X.  The provided means MEAN_X
  ! and MEAN_Y will be used if given, otherwise they will be
  ! calculated internally.
  function osl_cov_vwp(x, y, mean_x, mean_y) result(cov)
    real(wp), dimension(:), intent(in) :: x
    real(wp), dimension(size(x)), intent(in) :: y
    real(wp), intent(in), optional :: mean_x, mean_y
    real(wp) :: cov
    real(wp) :: m_x, m_y
    integer :: i

    if (present(mean_x)) then
       m_x = mean_x
    else
       m_x = osl_mean(x)
    end if

    if (present(mean_y)) then
       m_y = mean_y
    else
       m_y = osl_mean(y)
    end if

    cov = 0.0_wp
    do i = 1, size(x)
       cov = cov + (x(i) - m_x) * (y(i) - m_y)
    end do
    cov = cov / real(size(x) - 1, wp)
  end function osl_cov_vwp

  ! Return the sample covariance between the _columns_ of the matrices X
  ! and Y.  The provided means MEAN_X and MEAN_Y will be used if
  ! given, otherwise they will be calculated internally.
  function osl_cov_mwp(x, y, mean_x, mean_y) result(cov)
    real(wp), dimension(:,:), intent(in) :: x
    real(wp), dimension(size(X,1),size(X,2)), intent(in) :: y
    real(wp), dimension(size(X,1)), intent(in), optional :: mean_x, mean_y
    real(wp), dimension(size(X,1),size(X,1)) :: cov
    real(wp), dimension(size(X,1)) :: m_x, m_y
    integer :: i

    if (present(mean_x)) then
       m_x = mean_x
    else
       m_x = osl_mean(x)
    end if

    if (present(mean_y)) then
       m_y = mean_y
    else
       m_y = osl_mean(y)
    end if

    cov = 0.0_wp
    do i = 1, size(X,2)
       cov = cov + osl_outer_product(X(:,i) - m_x, Y(:,i) - m_y)
    end do
    cov = cov / real(size(X,2) - 1, wp)
  end function osl_cov_mwp

  !------------------------------------!
  ! Sample variance                    !
  !------------------------------------!

  ! Return the sample variance of a vector X.  The provided MEAN will
  ! be used if given, otherwise it will be calculated internally.
  function osl_var_vwp(x, mean) result(var)
    real(wp), dimension(:), intent(in) :: x
    real(wp), intent(in), optional :: mean
    real(wp) :: var, m

    if (.not. present(mean)) then
       m = osl_mean(x)
    else
       m = mean
    end if

    var = osl_cov(x, x, m, m)
  end function osl_var_vwp

  !------------------------------------!
  ! Combined sample statistics         !
  !------------------------------------!

  ! Return the sample covariance between the _columns_ of the matrices X
  ! and Y.  If Y is omitted, return the variance of X.  The provided
  ! means will be used if given, otherwise they will be calculated
  ! internally.
  function osl_var_mwp(x, mean) result(var)
    real(wp), dimension(:,:), intent(in) :: x
    real(wp), dimension(size(X,1)), intent(in), optional :: mean
    real(wp), dimension(size(X,1),size(X,1)) :: var
    real(wp), dimension(size(X,1)) :: m

    if (.not. present(mean)) then
       m = osl_mean(x)
    else
       m = mean
    end if

    var = osl_cov(x, x, m, m)
  end function osl_var_mwp

  ! Calculate sample statistics for a sample of vectors.  The
  ! observations must be stored as _columns_ of the variable SAMPLE.  Any
  ! of the following statistics may be calculated (all are optional):
  ! mean (MEAN), standard deviation (SD), variance (VAR), bias
  ! (BIAS), mean squared error (MSE), root mean squared error (RMSE).
  ! Note that TRUTH must be present in order to calculate BIAS, MSE,
  ! or RMSE (zero is returned otherwise).
  subroutine osl_sample_stats(sample, mean, var, sd, truth, bias, mse, rmse)
    real(wp), dimension(:,:), intent(in) :: sample
    real(wp), dimension(size(sample,1)), optional, intent(out) :: mean
    real(wp), dimension(size(sample,1)), optional, intent(out) :: var
    real(wp), dimension(size(sample,1)), optional, intent(out) :: sd
    real(wp), dimension(size(sample,1)), optional, intent(in) :: truth
    real(wp), dimension(size(sample,1)), optional, intent(out) :: bias
    real(wp), dimension(size(sample,1)), optional, intent(out) :: mse
    real(wp), dimension(size(sample,1)), optional, intent(out) :: rmse

    real(wp), dimension(size(sample,1)) :: m
    real(wp), dimension(size(sample,1)) :: v
    integer :: i

    m = osl_mean(sample)
    if (present(mean)) then
       mean = m
    end if

    ! Calculate the component-wise variance, not the covariance matrix
    v = 0.0_wp
    do i = 1, size(sample, 2)
       v = v + (sample(:,i) - m)**2
    end do
    v = v / real(size(sample, 2) - 1, wp)
    if (present(var)) then
       var = v
    end if

    if (present(sd)) then
       sd = sqrt(v)
    end if

    if (present(truth)) then
       if (present(bias)) then
          bias = m - truth
       end if
       if (present(mse)) then
          if (present(bias)) then
             mse = v + bias**2
          else
             mse = v + (m - truth)**2
          end if
       end if
       if (present(rmse)) then
          if (present(mse)) then
             rmse = sqrt(mse)
          else
             if (present(bias)) then
                rmse = sqrt(v + bias**2)
             else
                rmse = sqrt(v + (m - truth)**2)
             end if
          end if
       end if
    else
       if (present(bias)) bias = 0.0_wp
       if (present(mse)) mse = 0.0_wp
       if (present(rmse)) rmse = 0.0_wp
    end if
  end subroutine osl_sample_stats

end module osl
