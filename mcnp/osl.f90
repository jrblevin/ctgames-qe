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

  ! real precision constants
  integer, parameter :: osl_dp = selected_real_kind(15, 307)
  integer, parameter :: osl_wp = osl_dp

  ! integer precision constants
  integer, parameter :: osl_int32 = selected_int_kind(9)

  ! internal-use private precision constants
  integer, parameter, private :: dp = osl_dp
  integer, parameter, private :: wp = osl_wp

  ! precision constants
  real(dp), parameter :: osl_eps_dp = epsilon(1.0_dp)
  real(wp), parameter :: osl_eps_wp = epsilon(1.0_wp)

  ! mathematical constants
  real(wp), parameter :: OSL_EULER = 0.5772156649015328606065120900824024_wp

  ! Default edit descriptors
  character(len=*), parameter :: osl_ed_i   = 'i10'
  character(len=*), parameter :: osl_ed_wp  = 'g17.5'
  integer, parameter, private :: fmt_len    = 80

  ! A context mold for generic interfaces.
  integer, dimension(1), parameter :: osl_ctx = (/ 0 /)

  ! error handling behavior
  integer, parameter :: osl_error_msg_len = 64

  ! status/error codes
  integer, parameter :: OSL_FAILURE  = -1
  integer, parameter :: OSL_SUCCESS  = 0
  integer, parameter :: OSL_EDOM     = 1   ! input domain error
  integer, parameter :: OSL_ERANGE   = 2   ! output range error
  integer, parameter :: OSL_EINVALID = 3   ! invalid argument
  integer, parameter :: OSL_ENOMEM   = 8   ! memory allocation failed
  integer, parameter :: OSL_EMAXITER = 11  ! exceeded maximum iterations
  integer, parameter :: OSL_EMAXEVAL = 12  ! exceeded maximum evaluations
  integer, parameter :: OSL_EUNDRFLW = 15  ! underflow
  integer, parameter :: OSL_EOVRFLW  = 16  ! overflow
  integer, parameter :: OSL_EROUND   = 18  ! roundoff error
  integer, parameter :: OSL_ECONFORM = 19  ! dimensions are non-conforming
  integer, parameter :: OSL_ENOTSQ   = 20  ! matrix is not square
  integer, parameter :: OSL_ESING    = 21  ! singularity detected
  integer, parameter :: OSL_EPOSDEF  = 98  ! matrix not positive definite
  integer, parameter :: OSL_ESYM     = 98  ! matrix not symmetric
  integer, parameter :: OSL_EINIT    = 99  ! object not initialized

  ! Procedures for printing informative messages
  interface osl_print
     module procedure osl_print_str
     module procedure osl_print_swp
     module procedure osl_print_vwp
     module procedure osl_print_mwp
     module procedure osl_print_si
     module procedure osl_print_vi
     module procedure osl_print_mi
  end interface osl_print

  ! Norms
  interface osl_norm_2
     module procedure osl_norm_2_vwp
  end interface osl_norm_2

  interface osl_norm_inf
     module procedure osl_norm_inf_vwp
     module procedure osl_norm_inf_mwp
  end interface osl_norm_inf

  ! IEEE floating point procedures
  interface osl_is_inf
     module procedure osl_is_inf_wp
  end interface osl_is_inf

  ! Swap
  interface osl_swap
     module procedure osl_swap_scalar_wp
     module procedure osl_swap_scalar_i
     module procedure osl_swap_vector_wp
     module procedure osl_swap_vector_i
  end interface osl_swap

  public :: OSL_FACTORIAL_MAX

  public :: osl_factorial

  ! Maximum n such that n! does not overflow
  integer, parameter :: OSL_FACTORIAL_MAX = 170

  ! Lookup table for osl_factorial
  real(wp), dimension(OSL_FACTORIAL_MAX) :: osl_factorial_table
  data osl_factorial_table(1) / 1.0_wp /
  data osl_factorial_table(2) / 2.0_wp /
  data osl_factorial_table(3) / 6.0_wp /
  data osl_factorial_table(4) / 24.0_wp /
  data osl_factorial_table(5) / 120.0_wp /
  data osl_factorial_table(6) / 720.0_wp /
  data osl_factorial_table(7) / 5040.0_wp /
  data osl_factorial_table(8) / 40320.0_wp /
  data osl_factorial_table(9) / 362880.0_wp /
  data osl_factorial_table(10) / 3628800.0_wp /
  data osl_factorial_table(11) / 39916800.0_wp /
  data osl_factorial_table(12) / 479001600.0_wp /
  data osl_factorial_table(13) / 6227020800.0_wp /
  data osl_factorial_table(14) / 87178291200.0_wp /
  data osl_factorial_table(15) / 1307674368000.0_wp /
  data osl_factorial_table(16) / 20922789888000.0_wp /
  data osl_factorial_table(17) / 355687428096000.0_wp /
  data osl_factorial_table(18) / 6402373705728000.0_wp /
  data osl_factorial_table(19) / 121645100408832000.0_wp /
  data osl_factorial_table(20) / 2432902008176640000.0_wp /
  data osl_factorial_table(21) / 51090942171709440000.0_wp /
  data osl_factorial_table(22) / 1124000727777607680000.0_wp /
  data osl_factorial_table(23) / 25852016738884976640000.0_wp /
  data osl_factorial_table(24) / 620448401733239439360000.0_wp /
  data osl_factorial_table(25) / 15511210043330985984000000.0_wp /
  data osl_factorial_table(26) / 403291461126605635584000000.0_wp /
  data osl_factorial_table(27) / 10888869450418352160768000000.0_wp /
  data osl_factorial_table(28) / 304888344611713860501504000000.0_wp /
  data osl_factorial_table(29) / 8841761993739701954543616000000.0_wp /
  data osl_factorial_table(30) / 265252859812191058636308480000000.0_wp /
  data osl_factorial_table(31) / 8222838654177922817725562880000000.0_wp /
  data osl_factorial_table(32) / 263130836933693530167218012160000000.0_wp /
  data osl_factorial_table(33) / 8683317618811886495518194401280000000.0_wp /
  data osl_factorial_table(34) / 2.95232799039604140847618609644e38_wp /
  data osl_factorial_table(35) / 1.03331479663861449296666513375e40_wp /
  data osl_factorial_table(36) / 3.71993326789901217467999448151e41_wp /
  data osl_factorial_table(37) / 1.37637530912263450463159795816e43_wp /
  data osl_factorial_table(38) / 5.23022617466601111760007224100e44_wp /
  data osl_factorial_table(39) / 2.03978820811974433586402817399e46_wp /
  data osl_factorial_table(40) / 8.15915283247897734345611269600e47_wp /
  data osl_factorial_table(41) / 3.34525266131638071081700620534e49_wp /
  data osl_factorial_table(42) / 1.40500611775287989854314260624e51_wp /
  data osl_factorial_table(43) / 6.04152630633738356373551320685e52_wp /
  data osl_factorial_table(44) / 2.65827157478844876804362581101e54_wp /
  data osl_factorial_table(45) / 1.19622220865480194561963161496e56_wp /
  data osl_factorial_table(46) / 5.50262215981208894985030542880e57_wp /
  data osl_factorial_table(47) / 2.58623241511168180642964355154e59_wp /
  data osl_factorial_table(48) / 1.24139155925360726708622890474e61_wp /
  data osl_factorial_table(49) / 6.08281864034267560872252163321e62_wp /
  data osl_factorial_table(50) / 3.04140932017133780436126081661e64_wp /
  data osl_factorial_table(51) / 1.55111875328738228022424301647e66_wp /
  data osl_factorial_table(52) / 8.06581751709438785716606368564e67_wp /
  data osl_factorial_table(53) / 4.27488328406002556429801375339e69_wp /
  data osl_factorial_table(54) / 2.30843697339241380472092742683e71_wp /
  data osl_factorial_table(55) / 1.26964033536582759259651008476e73_wp /
  data osl_factorial_table(56) / 7.10998587804863451854045647464e74_wp /
  data osl_factorial_table(57) / 4.05269195048772167556806019054e76_wp /
  data osl_factorial_table(58) / 2.35056133128287857182947491052e78_wp /
  data osl_factorial_table(59) / 1.38683118545689835737939019720e80_wp /
  data osl_factorial_table(60) / 8.32098711274139014427634118320e81_wp /
  data osl_factorial_table(61) / 5.07580213877224798800856812177e83_wp /
  data osl_factorial_table(62) / 3.14699732603879375256531223550e85_wp /
  data osl_factorial_table(63) / 1.982608315404440064116146708360e87_wp /
  data osl_factorial_table(64) / 1.268869321858841641034333893350e89_wp /
  data osl_factorial_table(65) / 8.247650592082470666723170306800e90_wp /
  data osl_factorial_table(66) / 5.443449390774430640037292402480e92_wp /
  data osl_factorial_table(67) / 3.647111091818868528824985909660e94_wp /
  data osl_factorial_table(68) / 2.480035542436830599600990418570e96_wp /
  data osl_factorial_table(69) / 1.711224524281413113724683388810e98_wp /
  data osl_factorial_table(70) / 1.197857166996989179607278372170e100_wp /
  data osl_factorial_table(71) / 8.504785885678623175211676442400e101_wp /
  data osl_factorial_table(72) / 6.123445837688608686152407038530e103_wp /
  data osl_factorial_table(73) / 4.470115461512684340891257138130e105_wp /
  data osl_factorial_table(74) / 3.307885441519386412259530282210e107_wp /
  data osl_factorial_table(75) / 2.480914081139539809194647711660e109_wp /
  data osl_factorial_table(76) / 1.885494701666050254987932260860e111_wp /
  data osl_factorial_table(77) / 1.451830920282858696340707840860e113_wp /
  data osl_factorial_table(78) / 1.132428117820629783145752115870e115_wp /
  data osl_factorial_table(79) / 8.946182130782975286851441715400e116_wp /
  data osl_factorial_table(80) / 7.156945704626380229481153372320e118_wp /
  data osl_factorial_table(81) / 5.797126020747367985879734231580e120_wp /
  data osl_factorial_table(82) / 4.753643337012841748421382069890e122_wp /
  data osl_factorial_table(83) / 3.945523969720658651189747118010e124_wp /
  data osl_factorial_table(84) / 3.314240134565353266999387579130e126_wp /
  data osl_factorial_table(85) / 2.817104114380550276949479442260e128_wp /
  data osl_factorial_table(86) / 2.422709538367273238176552320340e130_wp /
  data osl_factorial_table(87) / 2.107757298379527717213600518700e132_wp /
  data osl_factorial_table(88) / 1.854826422573984391147968456460e134_wp /
  data osl_factorial_table(89) / 1.650795516090846108121691926250e136_wp /
  data osl_factorial_table(90) / 1.485715964481761497309522733620e138_wp /
  data osl_factorial_table(91) / 1.352001527678402962551665687590e140_wp /
  data osl_factorial_table(92) / 1.243841405464130725547532432590e142_wp /
  data osl_factorial_table(93) / 1.156772507081641574759205162310e144_wp /
  data osl_factorial_table(94) / 1.087366156656743080273652852570e146_wp /
  data osl_factorial_table(95) / 1.032997848823905926259970209940e148_wp /
  data osl_factorial_table(96) / 9.916779348709496892095714015400e149_wp /
  data osl_factorial_table(97) / 9.619275968248211985332842594960e151_wp /
  data osl_factorial_table(98) / 9.426890448883247745626185743100e153_wp /
  data osl_factorial_table(99) / 9.332621544394415268169923885600e155_wp /
  data osl_factorial_table(100) / 9.33262154439441526816992388563e157_wp /
  data osl_factorial_table(101) / 9.42594775983835942085162312450e159_wp /
  data osl_factorial_table(102) / 9.61446671503512660926865558700e161_wp /
  data osl_factorial_table(103) / 9.90290071648618040754671525458e163_wp /
  data osl_factorial_table(104) / 1.02990167451456276238485838648e166_wp /
  data osl_factorial_table(105) / 1.08139675824029090050410130580e168_wp /
  data osl_factorial_table(106) / 1.146280563734708354534347384148e170_wp /
  data osl_factorial_table(107) / 1.226520203196137939351751701040e172_wp /
  data osl_factorial_table(108) / 1.324641819451828974499891837120e174_wp /
  data osl_factorial_table(109) / 1.443859583202493582204882102460e176_wp /
  data osl_factorial_table(110) / 1.588245541522742940425370312710e178_wp /
  data osl_factorial_table(111) / 1.762952551090244663872161047110e180_wp /
  data osl_factorial_table(112) / 1.974506857221074023536820372760e182_wp /
  data osl_factorial_table(113) / 2.231192748659813646596607021220e184_wp /
  data osl_factorial_table(114) / 2.543559733472187557120132004190e186_wp /
  data osl_factorial_table(115) / 2.925093693493015690688151804820e188_wp /
  data osl_factorial_table(116) / 3.393108684451898201198256093590e190_wp /
  data osl_factorial_table(117) / 3.96993716080872089540195962950e192_wp /
  data osl_factorial_table(118) / 4.68452584975429065657431236281e194_wp /
  data osl_factorial_table(119) / 5.57458576120760588132343171174e196_wp /
  data osl_factorial_table(120) / 6.68950291344912705758811805409e198_wp /
  data osl_factorial_table(121) / 8.09429852527344373968162284545e200_wp /
  data osl_factorial_table(122) / 9.87504420083360136241157987140e202_wp /
  data osl_factorial_table(123) / 1.21463043670253296757662432419e205_wp /
  data osl_factorial_table(124) / 1.50614174151114087979501416199e207_wp /
  data osl_factorial_table(125) / 1.88267717688892609974376770249e209_wp /
  data osl_factorial_table(126) / 2.37217324288004688567714730514e211_wp /
  data osl_factorial_table(127) / 3.01266001845765954480997707753e213_wp /
  data osl_factorial_table(128) / 3.85620482362580421735677065923e215_wp /
  data osl_factorial_table(129) / 4.97450422247728744039023415041e217_wp /
  data osl_factorial_table(130) / 6.46685548922047367250730439554e219_wp /
  data osl_factorial_table(131) / 8.47158069087882051098456875820e221_wp /
  data osl_factorial_table(132) / 1.11824865119600430744996307608e224_wp /
  data osl_factorial_table(133) / 1.48727070609068572890845089118e226_wp /
  data osl_factorial_table(134) / 1.99294274616151887673732419418e228_wp /
  data osl_factorial_table(135) / 2.69047270731805048359538766215e230_wp /
  data osl_factorial_table(136) / 3.65904288195254865768972722052e232_wp /
  data osl_factorial_table(137) / 5.01288874827499166103492629211e234_wp /
  data osl_factorial_table(138) / 6.91778647261948849222819828311e236_wp /
  data osl_factorial_table(139) / 9.61572319694108900419719561353e238_wp /
  data osl_factorial_table(140) / 1.34620124757175246058760738589e241_wp /
  data osl_factorial_table(141) / 1.89814375907617096942852641411e243_wp /
  data osl_factorial_table(142) / 2.69536413788816277658850750804e245_wp /
  data osl_factorial_table(143) / 3.85437071718007277052156573649e247_wp /
  data osl_factorial_table(144) / 5.55029383273930478955105466055e249_wp /
  data osl_factorial_table(145) / 8.04792605747199194484902925780e251_wp /
  data osl_factorial_table(146) / 1.17499720439091082394795827164e254_wp /
  data osl_factorial_table(147) / 1.72724589045463891120349865931e256_wp /
  data osl_factorial_table(148) / 2.55632391787286558858117801578e258_wp /
  data osl_factorial_table(149) / 3.80892263763056972698595524351e260_wp /
  data osl_factorial_table(150) / 5.71338395644585459047893286526e262_wp /
  data osl_factorial_table(151) / 8.62720977423324043162318862650e264_wp /
  data osl_factorial_table(152) / 1.31133588568345254560672467123e267_wp /
  data osl_factorial_table(153) / 2.00634390509568239477828874699e269_wp /
  data osl_factorial_table(154) / 3.08976961384735088795856467036e271_wp /
  data osl_factorial_table(155) / 4.78914290146339387633577523906e273_wp /
  data osl_factorial_table(156) / 7.47106292628289444708380937294e275_wp /
  data osl_factorial_table(157) / 1.17295687942641442819215807155e278_wp /
  data osl_factorial_table(158) / 1.85327186949373479654360975305e280_wp /
  data osl_factorial_table(159) / 2.94670227249503832650433950735e282_wp /
  data osl_factorial_table(160) / 4.71472363599206132240694321176e284_wp /
  data osl_factorial_table(161) / 7.59070505394721872907517857094e286_wp /
  data osl_factorial_table(162) / 1.22969421873944943411017892849e289_wp /
  data osl_factorial_table(163) / 2.00440157654530257759959165344e291_wp /
  data osl_factorial_table(164) / 3.28721858553429622726333031164e293_wp /
  data osl_factorial_table(165) / 5.42391066613158877498449501421e295_wp /
  data osl_factorial_table(166) / 9.00369170577843736647426172359e297_wp /
  data osl_factorial_table(167) / 1.50361651486499904020120170784e300_wp /
  data osl_factorial_table(168) / 2.52607574497319838753801886917e302_wp /
  data osl_factorial_table(169) / 4.26906800900470527493925188890e304_wp /
  data osl_factorial_table(170) / 7.25741561530799896739672821113e306_wp /


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

  ! Discrete distributions
  type, public :: osl_discrete_rnd_t
     private
     integer :: n                                 ! size
     real(wp), dimension(:), pointer :: F         ! cutoffs
     integer, dimension(:), pointer :: L          ! alias values
  end type osl_discrete_rnd_t

  public :: osl_quicksort

  interface osl_quicksort
     module procedure osl_quicksort_vwp
  end interface osl_quicksort

  public :: osl_exponential_pdf, osl_exponential_cdf
  public :: osl_poisson_pdf, osl_poisson_inv

  interface osl_mean
     module procedure osl_mean_vwp
     module procedure osl_mean_mwp
     module procedure osl_mean_vi
     module procedure osl_mean_mi
  end interface

  interface osl_median
     module procedure osl_median_vwp
     module procedure osl_median_mwp
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

  !--------------------------------------------------------------------------!
  ! Formatted output                                                         !
  !--------------------------------------------------------------------------!

  !------------------------------------!
  ! Informative messages               !
  !------------------------------------!

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

  ! Prints a single string.
  subroutine osl_print_str(msg)
    character(len=*), intent(in) :: msg
    print '(a)', msg
  end subroutine osl_print_str

  !------------------------------------!
  ! Headings                           !
  !------------------------------------!

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

  !--------------------------------------------------------------------------!
  ! Error Handling                                                           !
  !--------------------------------------------------------------------------!

  ! Return a string describing the error with number ERRNO.
  function osl_error_msg(errno) result(msg)
    integer, intent(in) :: errno
    character(len=osl_error_msg_len) :: msg

    select case (errno)
    case (OSL_FAILURE)
       msg = 'Generic failure'
    case (OSL_SUCCESS)
       msg = 'Generic success'
    case (OSL_EDOM)
       msg = 'Domain (input) error'
    case (OSL_ERANGE)
       msg = 'Range (output) error'
    case (OSL_EINVALID)
       msg = 'Invalid argument'
    case (OSL_ENOMEM)
       msg = 'Memory allocation failed'
    case (OSL_EMAXITER)
       msg = 'Maximum number of iterations exceeded'
    case (OSL_EMAXEVAL)
       msg = 'Maximum number of evaluations exceeded'
    case (OSL_EUNDRFLW)
       msg = 'Numerical underflow'
    case (OSL_EOVRFLW)
       msg = 'Numerical overflow'
    case (OSL_ECONFORM)
       msg = 'Non-conformant dimensions'
    case (OSL_ENOTSQ)
       msg = 'Matrix is not square'
    case (OSL_EINIT)
       msg = 'Object not initialized'
    case default
       msg = 'Unknown error'
    end select
  end function osl_error_msg

  ! If STATUS is present, simply copy ERRNO to STATUS and return.
  ! Otherwise, print an error message and terminate the program.
  subroutine osl_status(errno, status, reason)
    integer, intent(in) :: errno
    integer, intent(out), optional :: status
    character(len=*), intent(in), optional :: reason

    if (present(status)) then
       status = errno
    else if (errno /= OSL_SUCCESS) then
       if (present(reason)) then
          print "('ERROR: ', a, ' (', i2.2, ')')", reason, errno
       else
          print "('ERROR: ', a, ' (', i2.2, ')')", osl_error_msg(errno), errno
       end if
       stop
    end if
  end subroutine osl_status

  !--------------------------------------------------------------------------!
  ! Norms                                                                    !
  !--------------------------------------------------------------------------!

  !------------------------------------!
  ! L_2 Norm                           !
  !------------------------------------!

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

  !------------------------------------!
  ! L_\infty Norm                      !
  !------------------------------------!

  ! Calculate the L_\infty norm of a vector.
  function osl_norm_inf_vwp(x) result(norm)
    real(wp), dimension(:), intent(in) :: x
    real(wp) :: norm
    norm = maxval(abs(x))
  end function osl_norm_inf_vwp

  ! Calculate the L_\infty norm of a matrix.
  function osl_norm_inf_mwp(x) result(norm)
    real(wp), dimension(:,:), intent(in) :: x
    real(wp) :: norm
    norm = maxval(abs(x))
  end function osl_norm_inf_mwp

  !--------------------------------------------------------------------------!
  ! IEEE Floating Point Arithmetic                                           !
  !--------------------------------------------------------------------------!

  ! Return .TRUE. if the given real variable is +Inf or -Inf and
  ! .FALSE. otherwise.
  elemental function osl_is_inf_wp(x)
    real(wp), intent(in):: x
    logical :: osl_is_inf_wp

    osl_is_inf_wp = (abs(x) > huge(x))
  end function osl_is_inf_wp

  !--------------------------------------------------------------------------!
  ! Swap                                                                     !
  !--------------------------------------------------------------------------!

  ! Swap the values of two REAL(WP) scalars.
  subroutine osl_swap_scalar_wp(a, b)
    real(wp), intent(inout) :: a, b
    real(wp) :: temp
    temp = a
    a = b
    b = temp
  end subroutine osl_swap_scalar_wp

  ! Swap the values of two INTEGER scalars.
  subroutine osl_swap_scalar_i(a, b)
    integer, intent(inout) :: a, b
    integer :: temp
    temp = a
    a = b
    b = temp
  end subroutine osl_swap_scalar_i

  ! Swaps the contents of two REAL(WP) vectors.
  subroutine osl_swap_vector_wp(a, b)
    real(wp), dimension(:), intent(inout) :: a, b
    real(wp), dimension(size(a)) :: temp
    temp = a
    a = b
    b = temp
  end subroutine osl_swap_vector_wp

  ! Swaps the contents of two INTEGER vectors.
  subroutine osl_swap_vector_i(a, b)
    integer, dimension(:), intent(inout) :: a, b
    integer, dimension(size(a)) :: temp
    temp = a
    a = b
    b = temp
  end subroutine osl_swap_vector_i

  !--------------------------------------------------------------------------!
  ! Vector-vector operations                                                 !
  !--------------------------------------------------------------------------!

  ! Return the outer product of the vectors A and B.
  function osl_outer_product(a, b) result(C)
    real(wp), dimension(:), intent(in) :: a
    real(wp), dimension(:), intent(in) :: b
    real(wp), dimension(size(a), size(b)) :: C
    C = spread(a, DIM=2, NCOPIES=size(b)) * spread(b, DIM=1, NCOPIES=size(a))
  end function osl_outer_product

  !------------------------------------!
  ! Factorial                          !
  !------------------------------------!

  ! Returns N!, equal to the product of the integers 1, 2, ..., N.
  !
  ! The return value is exact for values of N smaller than 18.
  ! For larger values of N, the error is on the order of
  ! 2 * EPSILON(1.0_WP) * N!.
  function osl_factorial(n, err, status)
    integer, intent(in) :: n
    real(wp), intent(out), optional :: err
    integer, intent(out), optional :: status
    real(wp) :: osl_factorial
    integer :: code
    real(wp) :: error

    error = 0.0_wp
    osl_factorial = 0.0_wp
    code = OSL_SUCCESS

    if (n < 0) then
       code = OSL_EDOM  ! domain is N >= 0
    else if (n == 0) then
       osl_factorial = 1.0_wp
    else if (n < 18) then
       osl_factorial = osl_factorial_table(n)
    else if (n <= osl_factorial_max) then
       osl_factorial = osl_factorial_table(n)
       error = 2.0_wp * epsilon(1.0_wp) * osl_factorial
    else
       code = OSL_EOVRFLW  ! would result in overflow
    end if

    if (present(err)) then
       err = error
    end if

    if (present(status)) then
       status = code
    end if
  end function osl_factorial

  !--------------------------------------------------------------------------!
  ! Mersenne Twister                                                         !
  !--------------------------------------------------------------------------!

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

  !--------------------------------------------------------------------------!
  ! Exponential Distribution                                                 !
  !--------------------------------------------------------------------------!

  function osl_exponential_pdf(x, lambda, status) result(pdf)
    real(wp), intent(in) :: x
    real(wp), intent(in) :: lambda
    integer, intent(out), optional :: status
    real(wp) :: pdf

    if (lambda < osl_eps_wp) then
       pdf = 0.0_wp
       call osl_status(OSL_EINVALID, status, &
            "rate parameter LAMBDA must be positive in OSL_EXPONENTIAL_PDF")
       return
    end if

    if (x >= 0.0_wp) then
       pdf = lambda * exp(-lambda * x)
    else
       pdf = 0.0_wp
    end if
    call osl_status(OSL_SUCCESS, status)
  end function osl_exponential_pdf

  function osl_exponential_cdf(x, lambda, status) result(cdf)
    real(wp), intent(in) :: x
    real(wp), intent(in) :: lambda
    integer, intent(out), optional :: status
    real(wp) :: cdf

    if (lambda < osl_eps_wp) then
       cdf = 0.0_wp
       call osl_status(OSL_EINVALID, status, &
            "rate parameter LAMBDA must be positive in OSL_EXPONENTIAL_CDF")
       return
    end if

    if (x >= 0.0_wp) then
       cdf = 1.0_wp - exp(-lambda * x)
    else
       cdf = 0.0_wp
    end if
    call osl_status(OSL_SUCCESS, status)
  end function osl_exponential_cdf

  !--------------------------------------------------------------------------!
  ! Poisson Distribution                                                     !
  !--------------------------------------------------------------------------!

  ! Returns the pdf of the Poisson distribution.
  function osl_poisson_pdf(k, mu, status) result(pdf)
    integer, intent(in) :: k
    real(wp), intent(in) :: mu
    integer, intent(out), optional :: status
    real(wp) :: pdf, lf

    if (mu < osl_eps_wp) then
       pdf = 0.0_wp
       call osl_status(OSL_EINVALID, status, &
            "invalid rate parameter in OSL_POISSON_CDF")
    else
       lf = log(osl_factorial(k)) ! This can be improved...
       pdf = exp(log(mu) * k - lf - mu)
       call osl_status(OSL_SUCCESS, status)
    end if
  end function osl_poisson_pdf

  ! Returns the inverse cdf of the Poisson distribution.
  !
  ! Using the relationship
  !
  ! k = osl_gamma_inv(p, lambda, 1) - 1
  function osl_poisson_inv(p, mu, status) result(k)
    real(wp), intent(in) :: p, mu
    integer, intent(out), optional :: status
    integer :: k
    integer, parameter :: limit = 20
    real(wp) :: cdf = 0.0_wp

    if (p < 0.0_wp .or. p > 1.0_wp) then
       k = 0
       call osl_status(OSL_EDOM, status, &
            "invalid quantile in OSL_POISSON_INV")
    else if (mu < osl_eps_wp) then
       k = 0
       call osl_status(OSL_EINVALID, status, &
            "invalid rate parameter in OSL_POISSON_INV")
    else
       k = 0
       cdf = 0.0_wp
       do
          cdf = cdf + osl_poisson_pdf(k, mu);
          if (cdf > p) then
             exit
          end if
          k = k + 1
       end do
       call osl_status(OSL_SUCCESS, status)
    end if
  end function osl_poisson_inv


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

  ! Return the mean of a vector of integers V.
  function osl_mean_vi(v) result(mean)
    integer, dimension(:) :: v
    real(wp) :: mean
    mean = real(sum(v), wp) / real(size(v), wp)
  end function osl_mean_vi

  ! Return the mean of the columns of an integer matrix A.
  function osl_mean_mi(A) result(mean)
    integer, dimension(:,:) :: A
    real(wp), dimension(size(A,1)) :: mean
    mean = real(sum(A, DIM=2), wp) / real(size(A,2), wp)
  end function osl_mean_mi

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

  ! Return the component-wise median across the columns of a real
  ! matrix X.  Given an M x N matrix X, the result is an M x 1 column
  ! vector where the elements are the component-wise medians of each
  ! row.
  function osl_median_mwp(x) result(median)
    real(wp), dimension(:,:) :: x
    real(wp), dimension(size(x,1)) :: median
    real(wp), dimension(size(x,2)) :: rowsort
    integer :: rows, cols, i, imed

    rows = size(x, 1)
    cols = size(x, 2)

    do i = 1, rows
       rowsort = x(i,:)
       call osl_quicksort(rowsort)
       if (mod(cols, 2) == 0) then
          imed = cols / 2
          median(i) = 0.5_wp * (rowsort(imed) + rowsort(imed+1))
       else
          imed = cols / 2 + 1
          median(i) = rowsort(imed)
       end if
    end do
  end function osl_median_mwp

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
