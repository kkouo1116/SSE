!------------------------------------------------------------------------------
! SSE reproduction of Fig. 75 (printed page 177) of arXiv:1101.3281,
! A. W. Sandvik, "Computational Studies of Quantum Spin Systems".
!
! Model: 2D columnar-dimerized S=1/2 Heisenberg antiferromagnet
!          H = sum_b J_b S_i(b).S_j(b),   J_1=1, g=J_2/J_1.
! Strong J_2 bonds join (x,y) to (x+1,y) for x=0,2,...,L-2.
! All other nearest-neighbor bonds have J_1. Periodic boundaries are used.
!
! Fig. 75 settings: L=8,16,32,64,128; beta=L (with J_1=1).
!   U_2     = (5/2) [1 - (1/3) <m_s,z^4>/<m_s,z^2>^2]   (Eq. 77, n=3)
!   rho_s,x = <(N_x^+ - N_x^-)^2>/(beta*N)             (Eq. 280)
!   L*rho_s,x = <W_x^2>*L/beta, W_x=(N_x^+ - N_x^-)/L  (Eq. 281)
! No extra 3/2 rotational factor is applied to the finite-size SSE stiffness.
!
! The implementation follows the fixed-length SSE operator string, diagonal
! update [Eqs. (265)-(266) generalized to beta_b=J_b/T], linked-vertex list,
! and operator-loop update in pseudocode {27}-{32}. Diagonal measurements are
! averaged over the n nontrivial propagated states, as in pseudocode {33}.
!
! Build:
!   gfortran -O3 -std=f2008 -Wall -Wextra SSE图75复现.f90 -o sse_fig75
! Full calculation (writes fig75_L*.dat, fig75_plot.gnu, and, if gnuplot is
! available, figure75.pdf):
!   ./sse_fig75
! Quick L=8 test with reduced sampling (run in a separate directory):
!   ./sse_fig75 quick
! Optional quick-mode sampling controls:
!   ./sse_fig75 quick 8 200 500
! where the three numbers are nbins, msteps per bin, and equilibration sweeps.
!
! Production defaults are nbins=20, msteps=2000, isteps=3000 per coupling.
! L=128 is computationally expensive. For publication-quality error bars,
! increase nbins and msteps; do not reduce beta below L for the target figure.
!------------------------------------------------------------------------------

module sse_fig75_module
  use, intrinsic :: iso_fortran_env, only : real64, int64
  implicit none
  private
  public :: sse_config, run_point, get_lattice_sizes, get_g_grid
  public :: write_plot_script

  type :: sse_config
     integer :: ell = 0                 ! spatial linear size L
     integer :: nsite = 0               ! N=L^2
     integer :: nbond = 0               ! N_b=2N
     integer :: n = 0                   ! current expansion order
     integer :: lcut = 0                ! fixed operator-string length
     real(real64) :: beta = 0.0_real64
     real(real64) :: g = 1.0_real64
     integer, allocatable :: spin(:)    ! sigma_i=+/-1, S_i^z=sigma_i/2
     integer, allocatable :: phase(:)   ! (-1)^(x+y)
     integer, allocatable :: op(:)      ! 0, 2*b (diagonal), 2*b+1 (off-diag.)
     integer, allocatable :: link(:)    ! linked-vertex array X(v), v=4*p+l
     integer, allocatable :: vfirst(:), vlast(:)
     integer, allocatable :: site1(:), site2(:)
     integer, allocatable :: bond_dir(:)! 1 for +x bond, 0 for +y bond
     real(real64), allocatable :: beta_b(:) ! J_b/T
  end type sse_config

contains

  function uniform_random() result(r)
    real(real64) :: r
    call random_number(r)
  end function uniform_random

  subroutine initialize_rng(seed, ell, g)
    integer, intent(in) :: seed, ell
    real(real64), intent(in) :: g
    integer :: nseed, i
    integer, allocatable :: seed_array(:)
    integer(int64) :: work

    call random_seed(size=nseed)
    allocate(seed_array(nseed))
    do i = 1, nseed
       work = int(seed, int64) + 104729_int64*int(i, int64) &
            + 1009_int64*int(ell, int64) + int(nint(100000*g), int64)
       seed_array(i) = int(modulo(work, 2147483646_int64) + 1_int64)
    end do
    call random_seed(put=seed_array)
  end subroutine initialize_rng

  integer function site_index(x, y, ell)
    integer, intent(in) :: x, y, ell
    site_index = modulo(x, ell) + ell*modulo(y, ell) + 1
  end function site_index

  subroutine initialize_config(s, ell, g, seed)
    type(sse_config), intent(out) :: s
    integer, intent(in) :: ell, seed
    real(real64), intent(in) :: g
    integer :: x, y, b, i

    if (ell < 2 .or. mod(ell, 2) /= 0) then
       error stop 'L must be an even integer >= 2.'
    end if
    if (g <= 0.0_real64) error stop 'The coupling ratio g=J_2/J_1 must be positive.'

    s%ell = ell
    s%nsite = ell*ell
    s%nbond = 2*s%nsite
    s%beta = real(ell, real64)       ! Fig. 75: beta=L, J_1=1
    s%g = g
    s%n = 0
    s%lcut = max(16, s%nsite/2)

    allocate(s%spin(s%nsite), s%phase(s%nsite))
    allocate(s%vfirst(s%nsite), s%vlast(s%nsite))
    allocate(s%site1(s%nbond), s%site2(s%nbond))
    allocate(s%bond_dir(s%nbond), s%beta_b(s%nbond))
    allocate(s%op(0:s%lcut-1), s%link(0:4*s%lcut-1))
    s%op = 0
    s%link = -1

    ! Strong columnar dimers: N/2 bonds, listed first (b=1,...,N/2).
    b = 0
    do y = 0, ell-1
       do x = 0, ell-2, 2
          b = b + 1
          s%site1(b) = site_index(x, y, ell)
          s%site2(b) = site_index(x+1, y, ell)
          s%bond_dir(b) = 1
          s%beta_b(b) = g*s%beta
       end do
    end do
    ! Weak x bonds: x=1,3,...,L-1, including the periodic seam.
    do y = 0, ell-1
       do x = 1, ell-1, 2
          b = b + 1
          s%site1(b) = site_index(x, y, ell)
          s%site2(b) = site_index(x+1, y, ell)
          s%bond_dir(b) = 1
          s%beta_b(b) = s%beta
       end do
    end do
    ! All y bonds are weak.
    do y = 0, ell-1
       do x = 0, ell-1
          b = b + 1
          s%site1(b) = site_index(x, y, ell)
          s%site2(b) = site_index(x, y+1, ell)
          s%bond_dir(b) = 0
          s%beta_b(b) = s%beta
       end do
    end do
    if (b /= s%nbond) error stop 'Internal bond-count error.'

    do y = 0, ell-1
       do x = 0, ell-1
          i = site_index(x, y, ell)
          if (mod(x+y, 2) == 0) then
             s%phase(i) = 1
          else
             s%phase(i) = -1
          end if
       end do
    end do

    call initialize_rng(seed, ell, g)
    do i = 1, s%nsite
       if (uniform_random() < 0.5_real64) then
          s%spin(i) = -1
       else
          s%spin(i) = 1
       end if
    end do
  end subroutine initialize_config

  subroutine diagonal_update(s)
    type(sse_config), intent(inout) :: s
    integer :: p, b, i, j
    real(real64) :: probability

    do p = 0, s%lcut-1
       if (s%op(p) == 0) then
          b = int(uniform_random()*real(s%nbond, real64)) + 1
          if (b > s%nbond) b = s%nbond
          i = s%site1(b)
          j = s%site2(b)
          if (s%spin(i) == s%spin(j)) cycle
          probability = s%beta_b(b)*real(s%nbond, real64) &
               / (2.0_real64*real(s%lcut-s%n, real64))
          if (uniform_random() < min(1.0_real64, probability)) then
             s%op(p) = 2*b
             s%n = s%n + 1
          end if
       else if (mod(s%op(p), 2) == 0) then
          b = s%op(p)/2
          probability = 2.0_real64*real(s%lcut-s%n+1, real64) &
               / (s%beta_b(b)*real(s%nbond, real64))
          if (uniform_random() < min(1.0_real64, probability)) then
             s%op(p) = 0
             s%n = s%n - 1
          end if
       else
          b = s%op(p)/2
          i = s%site1(b)
          j = s%site2(b)
          s%spin(i) = -s%spin(i)
          s%spin(j) = -s%spin(j)
       end if
    end do
  end subroutine diagonal_update

  subroutine loop_update(s)
    type(sse_config), intent(inout) :: s
    integer :: p, b, i, j, v0, v1, v2, f, l
    integer :: vstart, v, vnext, flag
    logical :: flip_loop

    ! Construction of the linked-vertex list, pseudocode {28}-{29}.
    s%link = -1
    s%vfirst = -1
    s%vlast = -1
    do p = 0, s%lcut-1
       if (s%op(p) == 0) cycle
       v0 = 4*p
       b = s%op(p)/2
       i = s%site1(b)
       j = s%site2(b)
       v1 = s%vlast(i)
       v2 = s%vlast(j)
       if (v1 /= -1) then
          s%link(v1) = v0
          s%link(v0) = v1
       else
          s%vfirst(i) = v0
       end if
       if (v2 /= -1) then
          s%link(v2) = v0+1
          s%link(v0+1) = v2
       else
          s%vfirst(j) = v0+1
       end if
       s%vlast(i) = v0+2
       s%vlast(j) = v0+3
    end do
    do i = 1, s%nsite
       f = s%vfirst(i)
       if (f /= -1) then
          l = s%vlast(i)
          s%link(f) = l
          s%link(l) = f
       end if
    end do

    ! Trace every loop once; flip with probability 1/2, pseudocode {30}-{31}.
    do vstart = 0, 4*s%lcut-1, 2
       if (s%link(vstart) < 0) cycle
       flip_loop = (uniform_random() < 0.5_real64)
       flag = -1
       if (flip_loop) flag = -2
       v = vstart
       do
          p = v/4
          if (flip_loop) s%op(p) = ieor(s%op(p), 1)
          s%link(v) = flag
          v = ieor(v, 1)
          vnext = s%link(v)
          s%link(v) = flag
          v = vnext
          if (v == vstart) exit
          if (v < 0) error stop 'Broken linked-vertex loop.'
       end do
    end do

    ! Update the stored state, including spins with no attached operators.
    do i = 1, s%nsite
       v = s%vfirst(i)
       if (v == -1) then
          if (uniform_random() < 0.5_real64) s%spin(i) = -s%spin(i)
       else
          if (s%link(v) == -2) s%spin(i) = -s%spin(i)
       end if
    end do
  end subroutine loop_update

  subroutine grow_operator_string(s)
    type(sse_config), intent(inout) :: s
    integer :: newcut
    integer, allocatable :: newop(:)

    ! Keep roughly one third of the string as fill-in identities.
    newcut = s%n + s%n/3 + 1
    if (newcut <= s%lcut) return
    allocate(newop(0:newcut-1))
    newop = 0
    newop(0:s%lcut-1) = s%op
    call move_alloc(newop, s%op)
    deallocate(s%link)
    allocate(s%link(0:4*newcut-1))
    s%link = -1
    s%lcut = newcut
  end subroutine grow_operator_string

  subroutine sse_sweep(s)
    type(sse_config), intent(inout) :: s
    call diagonal_update(s)
    call loop_update(s)
    call grow_operator_string(s)
  end subroutine sse_sweep

  subroutine measure(s, m2_avg, m4_avg, rho_l, energy)
    type(sse_config), intent(inout) :: s
    real(real64), intent(out) :: m2_avg, m4_avg, rho_l, energy
    integer :: p, b, i, j, ms, nwind
    integer :: initial_ms
    real(real64) :: m2, sum2, sum4, norm

    ! ms = 2*sum_i phi_i*S_i^z = sum_i phi_i*sigma_i.
    ms = sum(s%phase*s%spin)
    initial_ms = ms
    nwind = 0
    sum2 = 0.0_real64
    sum4 = 0.0_real64
    norm = 1.0_real64/(2.0_real64*real(s%nsite, real64))

    if (s%n == 0) then
       m2 = (real(ms, real64)*norm)**2
       m2_avg = m2
       m4_avg = m2*m2
    else
       do p = 0, s%lcut-1
          if (s%op(p) == 0) cycle
          if (mod(s%op(p), 2) == 1) then
             b = s%op(p)/2
             i = s%site1(b)
             j = s%site2(b)
             ! site1 -> site2 is the positive bond direction. An up spin
             ! initially at site1 is transported in +x on an x-directed bond.
             if (s%bond_dir(b) == 1) then
                if (s%spin(i) > 0) then
                   nwind = nwind + 1   ! N_x^+
                else
                   nwind = nwind - 1   ! N_x^-
                end if
             end if
             ! Both staggered contributions change by -2*phase(i)*spin(i).
             ms = ms - 4*s%phase(i)*s%spin(i)
             s%spin(i) = -s%spin(i)
             s%spin(j) = -s%spin(j)
          end if
          m2 = (real(ms, real64)*norm)**2
          sum2 = sum2 + m2
          sum4 = sum4 + m2*m2
       end do
       m2_avg = sum2/real(s%n, real64)
       m4_avg = sum4/real(s%n, real64)
    end if

    if (ms /= initial_ms) error stop 'Periodic propagation did not restore M_s.'
    ! rho_{s,x}*L = (N_x^+ - N_x^-)^2/(beta*L), since N=L^2.
    rho_l = real(nwind, real64)**2/(s%beta*real(s%ell, real64))
    ! Restore the 1/4 constants subtracted from all bond operators.
    energy = -real(s%n, real64)/(s%beta*real(s%nsite, real64)) &
         + (s%g+3.0_real64)/8.0_real64
  end subroutine measure

  pure real(real64) function binder_from_z(m2, m4)
    real(real64), intent(in) :: m2, m4
    ! For an O(3) order parameter: <m^2>=3<m_z^2>, <m^4>=5<m_z^4>.
    binder_from_z = 2.5_real64*(1.0_real64 - m4/(3.0_real64*m2*m2))
  end function binder_from_z

  subroutine run_point(ell, g, seed, nbins, msteps, isteps, &
       u2, u2_err, rho_l, rho_err, energy, energy_err, n_avg, max_n, final_cut)
    integer, intent(in) :: ell, seed, nbins, msteps, isteps
    real(real64), intent(in) :: g
    real(real64), intent(out) :: u2, u2_err, rho_l, rho_err
    real(real64), intent(out) :: energy, energy_err, n_avg
    integer, intent(out) :: max_n, final_cut
    type(sse_config) :: s
    real(real64), allocatable :: bm2(:), bm4(:), bu2(:), brho(:), be(:)
    real(real64), allocatable :: jack(:)
    real(real64) :: m2, m4, rho, e
    real(real64) :: sm2, sm4, srho, se, sn
    real(real64) :: mean, a, b
    integer :: ib, k

    if (nbins < 2 .or. msteps < 1 .or. isteps < 0) then
       error stop 'Require nbins>=2, msteps>=1, and isteps>=0.'
    end if
    allocate(bm2(nbins), bm4(nbins), bu2(nbins), brho(nbins), be(nbins))
    allocate(jack(nbins))
    call initialize_config(s, ell, g, seed)

    max_n = 0
    do k = 1, isteps
       call sse_sweep(s)
       max_n = max(max_n, s%n)
    end do

    sn = 0.0_real64
    do ib = 1, nbins
       sm2 = 0.0_real64
       sm4 = 0.0_real64
       srho = 0.0_real64
       se = 0.0_real64
       do k = 1, msteps
          call sse_sweep(s)
          call measure(s, m2, m4, rho, e)
          sm2 = sm2 + m2
          sm4 = sm4 + m4
          srho = srho + rho
          se = se + e
          sn = sn + real(s%n, real64)
          max_n = max(max_n, s%n)
       end do
       bm2(ib) = sm2/real(msteps, real64)
       bm4(ib) = sm4/real(msteps, real64)
       bu2(ib) = binder_from_z(bm2(ib), bm4(ib))
       brho(ib) = srho/real(msteps, real64)
       be(ib) = se/real(msteps, real64)
    end do

    ! Central Binder value from the globally averaged moments (not an average
    ! of per-sweep ratios). Jackknife bins supply its statistical error.
    a = sum(bm2)/real(nbins, real64)
    b = sum(bm4)/real(nbins, real64)
    u2 = binder_from_z(a, b)
    do ib = 1, nbins
       m2 = (sum(bm2)-bm2(ib))/real(nbins-1, real64)
       m4 = (sum(bm4)-bm4(ib))/real(nbins-1, real64)
       jack(ib) = binder_from_z(m2, m4)
    end do
    mean = sum(jack)/real(nbins, real64)
    u2_err = sqrt(real(nbins-1, real64)/real(nbins, real64) &
         * sum((jack-mean)**2))

    rho_l = sum(brho)/real(nbins, real64)
    rho_err = sqrt(sum((brho-rho_l)**2)/real(nbins*(nbins-1), real64))
    energy = sum(be)/real(nbins, real64)
    energy_err = sqrt(sum((be-energy)**2)/real(nbins*(nbins-1), real64))
    n_avg = sn/real(nbins*msteps, real64)
    final_cut = s%lcut
    if (max_n >= final_cut) error stop 'Operator-string truncation detected.'
  end subroutine run_point

  subroutine get_lattice_sizes(quick, sizes)
    logical, intent(in) :: quick
    integer, allocatable, intent(out) :: sizes(:)
    if (quick) then
       sizes = [8]
    else
       sizes = [8, 16, 32, 64, 128]
    end if
  end subroutine get_lattice_sizes

  subroutine get_g_grid(ell, g)
    integer, intent(in) :: ell
    real(real64), allocatable, intent(out) :: g(:)
    ! Sampling grids follow the windows and density of points in Fig. 75.
    select case (ell)
    case (8)
       g = [1.80_real64, 1.85_real64, 1.90_real64, 1.95_real64]
    case (16)
       g = [1.80_real64, 1.82_real64, 1.84_real64, 1.86_real64, &
            1.88_real64, 1.90_real64, 1.92_real64, 1.94_real64]
    case (32)
       g = [1.87_real64, 1.88_real64, 1.89_real64, 1.895_real64, &
            1.90_real64, 1.905_real64, 1.91_real64, 1.915_real64, &
            1.92_real64, 1.925_real64, 1.93_real64, 1.94_real64]
    case (64)
       g = [1.900_real64, 1.902_real64, 1.904_real64, 1.906_real64, &
            1.908_real64, 1.910_real64, 1.912_real64, 1.914_real64, &
            1.916_real64, 1.918_real64]
    case (128)
       g = [1.905_real64, 1.907_real64, 1.909_real64, &
            1.911_real64, 1.913_real64, 1.915_real64]
    case default
       error stop 'No Fig. 75 coupling grid is defined for this L.'
    end select
  end subroutine get_g_grid

  subroutine write_plot_script()
    integer :: u
    open(newunit=u, file='fig75_plot.gnu', status='replace', action='write')
    write(u, '(a)') 'set terminal pdfcairo enhanced size 10in,4.5in font "Times,14"'
    write(u, '(a)') 'set output "figure75.pdf"'
    write(u, '(a)') 'set multiplot layout 1,2'
    write(u, '(a)') 'set style line 1 lc rgb "black" lt 1 lw 1.2 pt 6 ps 0.75'
    write(u, '(a)') 'set style line 2 lc rgb "black" lt 1 lw 1.2 pt 7 ps 0.70'
    write(u, '(a)') 'set style line 3 lc rgb "red"   lt 1 lw 1.2 pt 6 ps 0.75'
    write(u, '(a)') 'set style line 4 lc rgb "red"   lt 1 lw 1.2 pt 7 ps 0.70'
    write(u, '(a)') 'set style line 5 lc rgb "black" lt 1 lw 1.2 pt 4 ps 0.70'
    write(u, '(a)') 'set xlabel "g"'
    write(u, '(a)') 'set xrange [1.80:1.975]'
    write(u, '(a)') 'set yrange [0.47:0.79]'
    write(u, '(a)') 'set ylabel "U_2"'
    write(u, '(a)') 'set key at 1.825,0.60 left'
    write(u, '(a)') 'plot "fig75_L8.dat"   u 1:2:3 w yerrorlines ls 1 t "L = 8",   \'
    write(u, '(a)') '     "fig75_L16.dat"  u 1:2:3 w yerrorlines ls 2 t "L = 16",  \'
    write(u, '(a)') '     "fig75_L32.dat"  u 1:2:3 w yerrorlines ls 3 t "L = 32",  \'
    write(u, '(a)') '     "fig75_L64.dat"  u 1:2:3 w yerrorlines ls 4 t "L = 64",  \'
    write(u, '(a)') '     "fig75_L128.dat" u 1:2:3 w yerrorlines ls 5 t "L = 128"'
    write(u, '(a)') 'set xrange [1.88:1.96]'
    write(u, '(a)') 'set yrange [0.80:1.25]'
    write(u, '(a)') 'set ylabel "{/Symbol r}_sL"'
    write(u, '(a)') 'set key at 1.958,1.235 right'
    write(u, '(a)') 'plot "fig75_L8.dat"   u 1:4:5 w yerrorlines ls 1 t "L = 8",   \'
    write(u, '(a)') '     "fig75_L16.dat"  u 1:4:5 w yerrorlines ls 2 t "L = 16",  \'
    write(u, '(a)') '     "fig75_L32.dat"  u 1:4:5 w yerrorlines ls 3 t "L = 32",  \'
    write(u, '(a)') '     "fig75_L64.dat"  u 1:4:5 w yerrorlines ls 4 t "L = 64",  \'
    write(u, '(a)') '     "fig75_L128.dat" u 1:4:5 w yerrorlines ls 5 t "L = 128"'
    write(u, '(a)') 'unset multiplot'
    close(u)
  end subroutine write_plot_script

end module sse_fig75_module

program reproduce_sse_fig75
  use, intrinsic :: iso_fortran_env, only : real64
  use sse_fig75_module
  implicit none
  integer, allocatable :: sizes(:)
  real(real64), allocatable :: g_grid(:)
  integer :: nbins, msteps, isteps, seed
  integer :: il, ig, unit_out, ios
  integer :: max_n, final_cut
  real(real64) :: u2, u2_err, rho_l, rho_err, energy, energy_err, n_avg
  logical :: quick
  character(len=32) :: arg, filename

  quick = .false.
  if (command_argument_count() >= 1) then
     call get_command_argument(1, arg)
     if (trim(arg) == 'quick') quick = .true.
  end if

  seed = 20261106
  nbins = 20
  msteps = 2000
  isteps = 3000
  if (quick) then
     nbins = 8
     msteps = 200
     isteps = 500
     if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg, *, iostat=ios) nbins
        if (ios /= 0) error stop 'Invalid nbins argument.'
     end if
     if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read(arg, *, iostat=ios) msteps
        if (ios /= 0) error stop 'Invalid msteps argument.'
     end if
     if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        read(arg, *, iostat=ios) isteps
        if (ios /= 0) error stop 'Invalid isteps argument.'
     end if
  end if

  call get_lattice_sizes(quick, sizes)
  write(*, '(a)') 'SSE reproduction of Fig. 75: columnar dimers, beta=L, J_1=1'
  write(*, '(a,*(i0,1x))') 'L values: ', sizes
  write(*, '(a,i0,a,i0,a,i0,a)') 'Bins=', nbins, ', sweeps/bin=', msteps, &
       ', equilibration=', isteps

  do il = 1, size(sizes)
     call get_g_grid(sizes(il), g_grid)
     write(filename, '(a,i0,a)') 'fig75_L', sizes(il), '.dat'
     open(newunit=unit_out, file=trim(filename), status='replace', action='write')
     write(unit_out, '(a)') '# Columns: g  U2  U2_jackknife_error  L*rho_sx  rho_error  E/N  E_error  <n>  max_n  final_cutoff'
     write(unit_out, '(a,i0,a,i0,a)') '# L=', sizes(il), ', beta=', sizes(il), ', columnar J_2 dimers along x'
     do ig = 1, size(g_grid)
        call run_point(sizes(il), g_grid(ig), seed, nbins, msteps, isteps, &
             u2, u2_err, rho_l, rho_err, energy, energy_err, n_avg, max_n, final_cut)
        write(unit_out, '(7(es20.10,1x),f20.6,1x,i0,1x,i0)') &
             g_grid(ig), u2, u2_err, rho_l, rho_err, energy, energy_err, &
             n_avg, max_n, final_cut
        flush(unit_out)
        write(*, '(a,i0,a,f8.4,a,f10.6,a,f10.6,a,f10.6,a,f10.6)') &
             'L=', sizes(il), ', g=', g_grid(ig), ', U2=', u2, &
             ' +/- ', u2_err, ', L*rho_sx=', rho_l, ' +/- ', rho_err
     end do
     close(unit_out)
  end do

  if (.not. quick) then
     call write_plot_script()
     call execute_command_line('gnuplot fig75_plot.gnu', exitstat=ios)
     if (ios == 0) then
        write(*, '(a)') 'Created figure75.pdf and fig75_plot.gnu.'
     else
        write(*, '(a)') 'Data complete. Run "gnuplot fig75_plot.gnu" to create the figure.'
     end if
  else
     write(*, '(a)') 'Quick test complete; full mode is required for the complete Fig. 75.'
  end if
end program reproduce_sse_fig75
