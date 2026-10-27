!> \file modfarquhar.f90
!! Farquhar et al. (1980) C3 leaf photosynthesis with Medlyn et al. (2011)
!! stomatal conductance, following CLM5 (Lawrence et al. 2019, CLM5.0 Technical
!! Note, chapter "Stomatal Resistance and Photosynthesis").
!!
!! Used by the multilayer canopy energy balance (modcanopy, lcanopyeb) when
!! ileafphys = 2 in NAMCANOPY.
!!
!! What follows CLM5.0 (src/biogeophys/PhotosynthesisMod.F90, release-clm5.0):
!!  - kinetic constants Kc, Ko, Gamma* and their temperature dependence
!!  - temperature functions (Arrhenius + high-temperature inhibition) for Vcmax,
!!    Jmax, TPU and Rd, with Kattge & Knorr (2007) acclimation of the entropy
!!    terms and of Jmax25/Vcmax25. The 10-day mean temperature of CLM is replaced
!!    by the fixed namelist parameter t_growth.
!!  - Jmax25, TPU25, Rd25 from Vcmax25 (CLM5 without LUNA, satellite-phenology mode)
!!  - electron transport, Rubisco/light/export limited rates and co-limitation
!!  - light inhibition of Rd (factor 0.67)
!!  - canopy decay of Vcmax25 with kn = exp(0.00963 Vcmax25top - 2.43)
!!  - hybrid secant/Brent solver for ci with CLM's tolerances
!! What follows CTSM (main branch, ci_func):
!!  - the Medlyn stomatal conductance without plant hydraulics. In release-clm5.0
!!    the non-hydraulics path always solves Ball-Berry.
!! What differs from CLM:
!!  - leaf boundary layer conductance comes from DALES (leafblc in modcanopy)
!!  - conversion between mol m-2 s-1 and m s-1 uses the local DALES air density
!!  - the ci solver is started from the previous value of ci (warm start)
!!  - optional (lfq_tleaf = .false.): biochemistry at air instead of leaf temperature,
!!    for comparison with A-gs, which uses the air temperature
!!  - no LUNA, no plant hydraulics, C3 only
!!
!!  Internal units are those of CLM: umol m-2 s-1 for fluxes and conductances,
!!  Pa for CO2 and O2 partial pressures, K for temperature.
!!  Conversion to DALES units is done in farquhar_dales.
!!
!!  \author Alexander Platter, with Claude (Anthropic), 2026
!
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
module modfarquhar
  implicit none
  private
  public :: initfarquhar, exitfarquhar, farquhar_set_params, farquhar_dales, farquhar_leaf, &
            vcmax25_can, ivcmax, fq_lim, fq_ncifunc

  save

  ! ---- Namelist NAMFARQUHAR
  integer :: ivcmax      = 1       !< Vcmax25 profile: 1 = vcmax25_top with exponential decay, 2 = from file vcmax.inp.<expnr>
  real    :: vcmax25_top = 41.     !< Vcmax at 25 degC at canopy top [umol m-2 s-1] (ivcmax = 1)
  real    :: kn_vcmax    = -1.     !< Decay coefficient of Vcmax25 with plant area above [-]; < 0: CLM5 formula (ivcmax = 1)
  real    :: t_growth    = 25.     !< Growth temperature for acclimation [degC] (CLM5: 10-day mean air temperature)
  logical :: lrdinhib    = .true.  !< Light inhibition of leaf dark respiration (Rd x 0.67 in light), as CLM5
  real    :: g0_medlyn   = 100.    !< Medlyn intercept [umol H2O m-2 s-1] (CLM5)
  real    :: g1_medlyn   = 4.12    !< Medlyn slope [kPa^0.5] (CLM5, broadleaf evergreen tropical tree)
  real    :: rd_vcmax    = 0.015   !< Rd25 / Vcmax25 [-] (CLM5, C3)
  real    :: tpu_vcmax   = 0.167   !< TPU25 / Vcmax25 [-] (CLM5)
  real    :: rs_max      = 2.e4    !< Maximum stomatal resistance [s m-1] (CLM5)
  logical :: lfq_tleaf   = .true.  !< Temperature of the biochemistry (kinetics, Vcmax, Jmax, TPU, Rd): leaf temperature
                                   !< (.true., CLM5) or air temperature (.false.); the leaf VPD always uses leaf temperature

  ! ---- Constants (CLM5.0, PhotosynthesisMod.F90)
  real, parameter :: rgas      = 8.31447    !< Universal gas constant [J mol-1 K-1]
  real, parameter :: tfrz      = 273.15     !< Freezing point [K]
  real, parameter :: t25       = tfrz + 25. !< Reference temperature [K]
  real, parameter :: o2_frac   = 0.209      !< O2 mole fraction [-]
  real, parameter :: kc25_frac = 404.9e-6   !< Kc at 25 degC [mol mol-1]
  real, parameter :: ko25_frac = 278.4e-3   !< Ko at 25 degC [mol mol-1]
  real, parameter :: cp25_frac = 42.75e-6   !< CO2 compensation point Gamma* at 25 degC [mol mol-1]
  real, parameter :: kcha      = 79430.     !< Activation energy for Kc [J mol-1]
  real, parameter :: koha      = 36380.     !< Activation energy for Ko [J mol-1]
  real, parameter :: cpha      = 37830.     !< Activation energy for Gamma* [J mol-1]
  real, parameter :: vcmaxha   = 72000.     !< Activation energy for Vcmax [J mol-1]
  real, parameter :: jmaxha    = 50000.     !< Activation energy for Jmax [J mol-1]
  real, parameter :: tpuha     = 72000.     !< Activation energy for TPU [J mol-1]
  real, parameter :: lmrha     = 46390.     !< Activation energy for Rd [J mol-1]
  real, parameter :: vcmaxhd   = 200000.    !< Deactivation energy for Vcmax [J mol-1]
  real, parameter :: jmaxhd    = 200000.    !< Deactivation energy for Jmax [J mol-1]
  real, parameter :: tpuhd     = 200000.    !< Deactivation energy for TPU [J mol-1]
  real, parameter :: lmrhd     = 150650.    !< Deactivation energy for Rd [J mol-1]
  real, parameter :: lmrse     = 490.       !< Entropy term for Rd [J mol-1 K-1]
  real, parameter :: fnps      = 0.15       !< Fraction of light absorbed by non-photosynthetic pigments [-]
  real, parameter :: theta_psii = 0.7       !< Curvature of electron transport light response [-]
  real, parameter :: theta_cj  = 0.98       !< Co-limitation Wc - Wj, C3 [-]
  real, parameter :: theta_ip  = 0.95       !< Co-limitation Ai - Wp [-]
  real, parameter :: rdinhib_fac = 0.67     !< Factor on Rd in light (CLM5 light_inhibit)
  real, parameter :: j_to_umol = 4.6        !< PAR energy to photons [umol J-1]
  real, parameter :: vpd_min   = 50.        !< Minimum vapour pressure deficit [Pa]

  ! Absorbed PAR below which a leaf is treated as dark [W m-2 leaf].
  ! modcanopy floors the incoming shortwave at the canopy top at 0.1 W m-2 also at
  ! night, which gives absorbed PAR of a few hundredths of W m-2; this threshold
  ! keeps such leaves dark. 0.1 W m-2 PAR is about 0.5 umol m-2 s-1.
  real, parameter :: apar_min  = 0.1

  ! ---- Derived at initialisation
  real :: tg_clamp        !< t_growth clamped to 11..35 degC as in CLM5
  real :: vcmaxse, jmaxse, tpuse          !< Acclimated entropy terms [J mol-1 K-1]
  real :: vcmaxc, jmaxc, tpuc, lmrc       !< Scaling factors for high-temperature inhibition
  real :: jmax_vcmax                      !< Jmax25 / Vcmax25 [-]
  real, allocatable :: vcmax25_can(:)     !< Vcmax at 25 degC per canopy level [umol m-2 s-1]

  ! ---- Diagnostics of the most recent leaf calculation (read by modcanopy)
  integer :: fq_lim     = 0   !< limiting rate: 0 = dark, 1 = Rubisco (Wc), 2 = light (Wj), 3 = export (Wp)
  integer :: fq_ncifunc = 0   !< number of ci_func evaluations since last reset (reset by the caller)

  ! ---- ci solver settings (CLM5 hybrid and brent)
  integer, parameter :: itmax_hybrid = 40
  integer, parameter :: itmax_brent  = 20
  real,    parameter :: eps_ci  = 1.e-2   !< relative tolerance on ci [-]
  real,    parameter :: tol_f   = 1.e-4   !< tolerance on ci residual [Pa]

  ! Leaf environment shared between the ci solver and ci_func, set in farquhar_leaf
  type leafenv
    real :: vcmax, jmax_pot, tpu, rd   !< temperature-corrected capacities [umol m-2 s-1]
    real :: je                         !< electron transport rate [umol m-2 s-1]
    real :: cp, kc, ko, oair           !< Gamma*, Kc, Ko, O2 [Pa]
    real :: ca, pres                   !< CO2 partial pressure in air, air pressure [Pa]
    real :: gb_mol                     !< boundary layer conductance for H2O [umol m-2 s-1]
    real :: vpd_kpa                    !< vapour pressure deficit leaf - air [kPa]
  end type leafenv

contains

  !> Read NAMFARQUHAR and set up the Vcmax25 profile on the canopy levels.
  !! \param ncanopy number of canopy levels
  !! \param paif    plant area index above the middle of each canopy level [m2 m-2]
  subroutine initfarquhar(ncanopy, paif)
    use modmpi,    only : myid, mpi_logical, mpi_integer, my_real, comm3d, mpierr
    use modglobal, only : ifnamopt, fname_options, ifinput, cexpnr, zf
    implicit none
    integer, intent(in) :: ncanopy
    real,    intent(in) :: paif(ncanopy)

    integer :: ierr, k, n, nlines
    real    :: kn
    real, allocatable :: zin(:), vcin(:), tin(:)
    character(80) :: readstring

    namelist/NAMFARQUHAR/ ivcmax, vcmax25_top, kn_vcmax, t_growth, lrdinhib, &
                          g0_medlyn, g1_medlyn, rd_vcmax, tpu_vcmax, rs_max, lfq_tleaf

    if (myid==0) then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMFARQUHAR,iostat=ierr)
      if (ierr > 0) then
        print *, 'Problem in namoptions NAMFARQUHAR'
        print *, 'iostat error: ', ierr
        stop 'ERROR: Problem in namoptions NAMFARQUHAR'
      endif
      write(6 ,NAMFARQUHAR)
      close(ifnamopt)

      if (ivcmax /= 1 .and. ivcmax /= 2) stop 'ERROR: ivcmax in NAMFARQUHAR must be 1 or 2'
    endif

    call MPI_BCAST(ivcmax      , 1, mpi_integer, 0, comm3d, mpierr)
    call MPI_BCAST(vcmax25_top , 1, my_real    , 0, comm3d, mpierr)
    call MPI_BCAST(kn_vcmax    , 1, my_real    , 0, comm3d, mpierr)
    call MPI_BCAST(t_growth    , 1, my_real    , 0, comm3d, mpierr)
    call MPI_BCAST(lrdinhib    , 1, mpi_logical, 0, comm3d, mpierr)
    call MPI_BCAST(g0_medlyn   , 1, my_real    , 0, comm3d, mpierr)
    call MPI_BCAST(g1_medlyn   , 1, my_real    , 0, comm3d, mpierr)
    call MPI_BCAST(rd_vcmax    , 1, my_real    , 0, comm3d, mpierr)
    call MPI_BCAST(tpu_vcmax   , 1, my_real    , 0, comm3d, mpierr)
    call MPI_BCAST(rs_max      , 1, my_real    , 0, comm3d, mpierr)
    call MPI_BCAST(lfq_tleaf   , 1, mpi_logical, 0, comm3d, mpierr)

    call farquhar_set_params

    allocate(vcmax25_can(ncanopy))

    select case (ivcmax)
    case (1)
      if (kn_vcmax < 0.) then
        kn = exp(0.00963 * vcmax25_top - 2.43)   ! CLM5
      else
        kn = kn_vcmax
      endif
      do k = 1, ncanopy
        vcmax25_can(k) = vcmax25_top * exp(-kn * paif(k))
      end do

    case (2)
      ! File vcmax.inp.<expnr>: lines starting with # are comments; data lines contain
      !   height [m]   Vcmax [umol m-2 s-1]   leaf temperature during measurement [degC]
      ! Heights increasing. Vcmax is converted to 25 degC with the temperature function
      ! of this module, then interpolated linearly to the canopy levels (constant
      ! beyond the lowest and highest height).
      if (myid==0) then
        open(ifinput, file='vcmax.inp.'//cexpnr, status='old', iostat=ierr)
        if (ierr /= 0) stop 'ERROR: ivcmax = 2 but file vcmax.inp.<expnr> not found'
        nlines = 0
        do
          read(ifinput, '(a80)', iostat=ierr) readstring
          if (ierr /= 0) exit
          if (readstring(1:1) /= '#' .and. len_trim(readstring) > 0) nlines = nlines + 1
        end do
        if (nlines < 1) stop 'ERROR: no data in vcmax.inp.<expnr>'
        allocate(zin(nlines), vcin(nlines), tin(nlines))
        rewind(ifinput)
        n = 0
        do
          read(ifinput, '(a80)', iostat=ierr) readstring
          if (ierr /= 0) exit
          if (readstring(1:1) == '#' .or. len_trim(readstring) == 0) cycle
          n = n + 1
          read(readstring, *) zin(n), vcin(n), tin(n)
        end do
        close(ifinput)
        do n = 2, nlines
          if (zin(n) <= zin(n-1)) stop 'ERROR: heights in vcmax.inp.<expnr> must increase'
        end do

        do n = 1, nlines
          vcin(n) = vcin(n) / (ft(tin(n)+tfrz, vcmaxha) * fth(tin(n)+tfrz, vcmaxhd, vcmaxse, vcmaxc))
        end do

        do k = 1, ncanopy
          if (zf(k) <= zin(1)) then
            vcmax25_can(k) = vcin(1)
          else if (zf(k) >= zin(nlines)) then
            vcmax25_can(k) = vcin(nlines)
          else
            do n = 1, nlines-1
              if (zf(k) >= zin(n) .and. zf(k) < zin(n+1)) then
                vcmax25_can(k) = vcin(n) + (vcin(n+1) - vcin(n)) * (zf(k) - zin(n)) / (zin(n+1) - zin(n))
                exit
              endif
            end do
          endif
        end do
        deallocate(zin, vcin, tin)
      endif
      call MPI_BCAST(vcmax25_can, ncanopy, my_real, 0, comm3d, mpierr)
    end select

    if (myid==0) then
      write(6,*) 'Farquhar leaf physiology: Vcmax25 per canopy level [umol m-2 s-1]'
      write(6,'(a6,2a12)') 'k', 'zf [m]', 'Vcmax25'
      do k = 1, ncanopy
        write(6,'(i6,2f12.3)') k, zf(k), vcmax25_can(k)
      end do
    endif

  end subroutine initfarquhar

  !> Derived parameters: acclimation (Kattge & Knorr 2007, as in CLM5) with fixed
  !! growth temperature. Separate from initfarquhar so that the leaf routines can
  !! also be used without the rest of DALES (tests).
  subroutine farquhar_set_params
    implicit none
    tg_clamp   = min(max(t_growth, 11.), 35.)
    vcmaxse    = 668.39 - 1.07 * tg_clamp
    jmaxse     = 659.70 - 0.75 * tg_clamp
    tpuse      = vcmaxse
    jmax_vcmax = 2.59 - 0.035 * tg_clamp
    vcmaxc     = fth25(vcmaxhd, vcmaxse)
    jmaxc      = fth25(jmaxhd , jmaxse )
    tpuc       = fth25(tpuhd  , tpuse  )
    lmrc       = fth25(lmrhd  , lmrse  )
  end subroutine farquhar_set_params

  subroutine exitfarquhar
    implicit none
    if (allocated(vcmax25_can)) deallocate(vcmax25_can)
  end subroutine exitfarquhar

  !> Interface to the DALES canopy: converts DALES units to CLM units, calls
  !! farquhar_leaf and converts the results back.
  subroutine farquhar_dales(co2_ppb, qtair, dens, pres, tleaf, tairk, phi_tot, apar, &   ! in
                            gb_heat, vcmax25,                                         & ! in
                            gcc, fleaf, ci_dales, co2abs, cp_dales, btran)              ! out (ci_dales in/out)
    use modsurfdata, only : phiwp, phifc, MW_CO2, MW_Air, nuco2q
    implicit none
    real, intent(in)    :: co2_ppb   !< CO2 in air [ppb]
    real, intent(in)    :: qtair     !< specific humidity of air [kg kg-1]
    real, intent(in)    :: dens      !< air density [kg m-3]
    real, intent(in)    :: pres      !< air pressure [Pa]
    real, intent(in)    :: tleaf     !< leaf temperature [K]
    real, intent(in)    :: tairk     !< air temperature [K] (used for the biochemistry if .not. lfq_tleaf)
    real, intent(in)    :: phi_tot   !< soil water content [-]
    real, intent(in)    :: apar      !< absorbed PAR [W m-2 leaf]
    real, intent(in)    :: gb_heat   !< leaf boundary layer conductance for heat [m s-1]
    real, intent(in)    :: vcmax25   !< Vcmax at 25 degC of this level [umol m-2 s-1]
    real, intent(out)   :: gcc       !< stomatal conductance for CO2 [m s-1]
    real, intent(out)   :: fleaf     !< net CO2 flux, negative = uptake [mg CO2 m-2 s-1]
    real, intent(inout) :: ci_dales  !< internal CO2 [mg m-3]; in: previous value (warm start)
    real, intent(out)   :: co2abs    !< CO2 in air [mg m-3]
    real, intent(out)   :: cp_dales  !< CO2 compensation point Gamma* [mg m-3]
    real, intent(out)   :: btran     !< soil water stress factor [-]

    real :: cf, ca, eair, gb_mol, ci, an, gs_mol, ag, cp, ppm2mg, tbio

    ! molar density of air [umol m-3], converts m s-1 to umol m-2 s-1. Taken from the
    ! DALES density (CLM: P/(R T)) so that it matches the ppb <-> mg m-3 conversion.
    cf      = dens / (MW_Air * 1.e-3) * 1.e6
    ppm2mg  = (MW_CO2/MW_Air) * dens             ! [mg m-3 per ppm]
    ca      = co2_ppb * 1.e-9 * pres             ! [Pa]
    eair    = qtair * pres / 0.622               ! [Pa], as in f_Ags
    gb_mol  = 1.075 * gb_heat * cf               ! H2O, 1.075 as in modcanopy (Goudriaan and van Laar 1994)
    btran   = max(1.0e-3, min(1.0, (phi_tot-phiwp)/(phifc-phiwp)))   ! as in f_Ags

    ! warm start: previous ci if it is a sensible number (it is not set at the first call)
    ci = ci_dales / ppm2mg * 1.e-6 * pres        ! [Pa]
    if (.not. (ci > 0. .and. ci < 2.*ca)) ci = 0.7 * ca

    if (lfq_tleaf) then
      tbio = tleaf
    else
      tbio = tairk
    endif

    call farquhar_leaf(tbio, tleaf, pres, ca, eair, apar, gb_mol, vcmax25, btran, & ! in
                       ci, an, gs_mol, ag, cp)                                ! out (ci in/out)

    gs_mol   = max(gs_mol, cf / rs_max)                       ! CLM5: rs <= rs_max
    ci       = ca - an * pres * (1.4*gs_mol + 1.6*gb_mol) / (gb_mol*gs_mol)
    ci       = max(ci, 1.e-6)

    gcc      = gs_mol / cf / nuco2q                           ! H2O to CO2 conductance [m s-1]
    fleaf    = -an * MW_CO2 * 1.e-3                           ! [umol m-2 s-1] to [mg m-2 s-1], DALES sign
    ci_dales = ci / pres * 1.e6 * ppm2mg                      ! [mg m-3]
    co2abs   = co2_ppb * 1.e-3 * ppm2mg                       ! [mg m-3]
    cp_dales = cp / pres * 1.e6 * ppm2mg                      ! [mg m-3]

  end subroutine farquhar_dales

  !> Leaf photosynthesis and stomatal conductance (CLM5 units).
  subroutine farquhar_leaf(tbio, tleaf, pres, ca, eair, apar, gb_mol, vcmax25, btran, & ! in
                           ci, an, gs_mol, ag, cp)                                ! out (ci in/out)
    implicit none
    real, intent(in)    :: tbio     !< temperature of the biochemistry [K] (leaf temperature in CLM5)
    real, intent(in)    :: tleaf    !< leaf temperature [K], for the leaf vapour pressure deficit
    real, intent(in)    :: pres     !< air pressure [Pa]
    real, intent(in)    :: ca       !< CO2 partial pressure in air [Pa]
    real, intent(in)    :: eair     !< vapour pressure of air [Pa]
    real, intent(in)    :: apar     !< absorbed PAR [W m-2 leaf]
    real, intent(in)    :: gb_mol   !< boundary layer conductance for H2O [umol m-2 s-1]
    real, intent(in)    :: vcmax25  !< Vcmax at 25 degC [umol m-2 s-1]
    real, intent(in)    :: btran    !< soil water stress factor [-]
    real, intent(inout) :: ci       !< internal CO2 [Pa]; in: first guess
    real, intent(out)   :: an       !< net assimilation [umol CO2 m-2 s-1]
    real, intent(out)   :: gs_mol   !< stomatal conductance for H2O [umol m-2 s-1]
    real, intent(out)   :: ag       !< gross assimilation [umol CO2 m-2 s-1]
    real, intent(out)   :: cp       !< CO2 compensation point Gamma* [Pa]

    type(leafenv) :: le
    real :: qabs, r1, r2, esat, fval

    ! kinetic constants [Pa]
    le%pres = pres
    le%ca   = ca
    le%oair = o2_frac * pres
    le%kc   = kc25_frac * pres * ft(tbio, kcha)
    le%ko   = ko25_frac * pres * ft(tbio, koha)
    le%cp   = cp25_frac * pres * ft(tbio, cpha)
    cp      = le%cp

    ! capacities at leaf temperature; water stress on Vcmax and Rd as in CLM5
    le%vcmax    = vcmax25 * btran * ft(tbio, vcmaxha) * fth(tbio, vcmaxhd, vcmaxse, vcmaxc)
    le%jmax_pot = jmax_vcmax * vcmax25 * ft(tbio, jmaxha) * fth(tbio, jmaxhd, jmaxse, jmaxc)
    le%tpu      = tpu_vcmax  * vcmax25 * ft(tbio, tpuha)  * fth(tbio, tpuhd,  tpuse,  tpuc)
    le%rd       = rd_vcmax   * vcmax25 * btran * ft(tbio, lmrha) * fth(tbio, lmrhd, lmrse, lmrc)

    le%gb_mol = gb_mol
    esat      = 0.611e3 * exp(17.2694 * (tleaf - 273.16) / (tleaf - 35.86))   ! as in f_Ags [Pa]
    le%vpd_kpa = max(esat - min(eair, esat), vpd_min) * 1.e-3

    if (apar <= apar_min) then
      ! dark: no photosynthesis, stomata at minimum conductance
      ag     = 0.
      an     = -le%rd
      gs_mol = g0_medlyn
      fq_lim = 0
      return
    endif

    if (lrdinhib) le%rd = rdinhib_fac * le%rd

    ! electron transport: smaller root of theta J^2 - (I + Jmax) J + I Jmax = 0
    qabs = 0.5 * (1. - fnps) * apar * j_to_umol
    call quadratic(theta_psii, -(qabs + le%jmax_pot), qabs * le%jmax_pot, r1, r2)
    le%je = min(r1, r2)

    call ci_hybrid(le, ci)
    call ci_func(le, ci, fval, an, gs_mol, ag)

  end subroutine farquhar_leaf

  !> Residual of the ci equation ci = ca - (1.4 rb + 1.6 rs) P An, with Medlyn gs.
  !! Follows CTSM ci_func (Medlyn branch).
  subroutine ci_func(le, ci, fval, an, gs_mol, ag)
    implicit none
    type(leafenv), intent(in) :: le
    real, intent(in)  :: ci
    real, intent(out) :: fval, an, gs_mol, ag
    real :: wc, wj, wp, ai, r1, r2, cs, term, aquad, bquad, cquad

    fq_ncifunc = fq_ncifunc + 1

    if (ci > le%cp) then
      wc = le%vcmax * (ci - le%cp) / (ci + le%kc * (1. + le%oair / le%ko))
      wj = le%je    * (ci - le%cp) / (4. * ci + 8. * le%cp)
    else            ! same as CLM's max(ci - cp, 0), but avoids division by zero for ci < 0
      wc = 0.
      wj = 0.
    endif
    wp = 3. * le%tpu

    call quadratic(theta_cj, -(wc + wj), wc * wj, r1, r2)
    ai = min(r1, r2)
    call quadratic(theta_ip, -(ai + wp), ai * wp, r1, r2)
    ag = max(0., min(r1, r2))
    an = ag - le%rd

    if (wp < ai) then
      fq_lim = 3
    else if (wc <= wj) then
      fq_lim = 1
    else
      fq_lim = 2
    endif

    if (an < 0.) then
      ! CTSM: stop iterating, stomata at minimum conductance
      gs_mol = g0_medlyn
      fval   = 0.
      return
    endif

    cs = max(le%ca - 1.4 / le%gb_mol * an * le%pres, 1.e-6)

    ! Medlyn with leaf surface VPD folded into a quadratic for gs [mol m-2 s-1]
    term  = 1.6 * an / (cs / le%pres * 1.e6)
    aquad = 1.0
    bquad = -(2.0 * (g0_medlyn*1.e-6 + term) + (g1_medlyn * term)**2 / (le%gb_mol*1.e-6 * le%vpd_kpa))
    cquad = (g0_medlyn*1.e-6)**2 + (2.0*g0_medlyn*1.e-6 + term * (1.0 - g1_medlyn**2 / le%vpd_kpa)) * term
    call quadratic(aquad, bquad, cquad, r1, r2)
    gs_mol = max(r1, r2) * 1.e6

    fval = ci - le%ca + an * le%pres * (1.4 * gs_mol + 1.6 * le%gb_mol) / (le%gb_mol * gs_mol)

  end subroutine ci_func

  !> Secant iteration on ci, switching to Brent once the root is bracketed (CLM5 hybrid).
  subroutine ci_hybrid(le, x0)
    implicit none
    type(leafenv), intent(in) :: le
    real, intent(inout) :: x0
    real :: x1, f0, f1, x, dx, tol, minx, minf, an, gs, ag
    integer :: iter

    call ci_func(le, x0, f0, an, gs, ag)
    if (f0 == 0.) return

    minx = x0
    minf = f0
    x1   = 0.99 * x0
    call ci_func(le, x1, f1, an, gs, ag)
    if (f1 == 0.) then
      x0 = x1
      return
    endif
    if (abs(f1) < abs(minf)) then
      minx = x1
      minf = f1
    endif

    iter = 0
    do
      iter = iter + 1
      if (f1 == f0) then           ! flat secant: keep best point
        x0 = minx
        exit
      endif
      dx  = -f1 * (x1 - x0) / (f1 - f0)
      x   = x1 + dx
      tol = abs(x) * eps_ci
      if (abs(dx) < tol) then
        x0 = x
        exit
      endif
      x0 = x1
      f0 = f1
      x1 = x
      call ci_func(le, x1, f1, an, gs, ag)
      if (abs(f1) < abs(minf)) then
        minx = x1
        minf = f1
      endif
      if (abs(f1) <= tol_f) then
        x0 = x1
        exit
      endif
      if (f1 * f0 < 0.) then
        call ci_brent(le, x0, x1, f0, f1, tol, x)
        x0 = x
        exit
      endif
      if (iter > itmax_hybrid) then
        x0 = minx
        exit
      endif
    end do

  end subroutine ci_hybrid

  !> Brent's method on a bracketing interval [x1, x2] (as CLM5 brent).
  subroutine ci_brent(le, x1, x2, f1, f2, tol, xout)
    implicit none
    type(leafenv), intent(in) :: le
    real, intent(in)  :: x1, x2, f1, f2, tol
    real, intent(out) :: xout
    real :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm, an, gs, ag
    integer :: iter

    a  = x1
    b  = x2
    fa = f1
    fb = f2
    c  = b
    fc = fb
    d  = b - a
    e  = d
    do iter = 1, itmax_brent
      if ((fb > 0. .and. fc > 0.) .or. (fb < 0. .and. fc < 0.)) then
        c  = a
        fc = fa
        d  = b - a
        e  = d
      endif
      if (abs(fc) < abs(fb)) then
        a  = b
        b  = c
        c  = a
        fa = fb
        fb = fc
        fc = fa
      endif
      tol1 = 2. * eps_ci * abs(b) + 0.5 * tol
      xm   = 0.5 * (c - b)
      if (abs(xm) <= tol1 .or. fb == 0.) exit
      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s = fb / fa
        if (a == c) then
          p = 2. * xm * s
          q = 1. - s
        else
          q = fa / fc
          r = fb / fc
          p = s * (2. * xm * q * (q - r) - (b - a) * (r - 1.))
          q = (q - 1.) * (r - 1.) * (s - 1.)
        endif
        if (p > 0.) q = -q
        p = abs(p)
        if (2. * p < min(3. * xm * q - abs(tol1 * q), abs(e * q))) then
          e = d
          d = p / q
        else
          d = xm
          e = d
        endif
      else
        d = xm
        e = d
      endif
      a  = b
      fa = fb
      if (abs(d) > tol1) then
        b = b + d
      else
        b = b + sign(tol1, xm)
      endif
      call ci_func(le, b, fb, an, gs, ag)
      if (fb == 0.) exit
    end do
    xout = b

  end subroutine ci_brent

  !> Roots of a x^2 + b x + c = 0 (numerically stable form, as CLM quadratic)
  subroutine quadratic(a, b, c, r1, r2)
    implicit none
    real, intent(in)  :: a, b, c
    real, intent(out) :: r1, r2
    real :: disc, q
    disc = b*b - 4.*a*c
    if (disc < 0.) disc = 0.    ! round-off
    q = -0.5 * (b + sign(sqrt(disc), b))
    if (q /= 0.) then
      r1 = q / a
      r2 = c / q
    else
      r1 = 0.
      r2 = 0.
    endif
  end subroutine quadratic

  !> Arrhenius temperature function, 1 at 25 degC
  real function ft(tl, ha)
    implicit none
    real, intent(in) :: tl, ha
    ft = exp(ha / (rgas * t25) * (1. - t25 / tl))
  end function ft

  !> High-temperature inhibition, 1 at 25 degC when scaled with fth25
  real function fth(tl, hd, se, scalef)
    implicit none
    real, intent(in) :: tl, hd, se, scalef
    fth = scalef / (1. + exp((-hd + se * tl) / (rgas * tl)))
  end function fth

  !> Scaling factor for fth
  real function fth25(hd, se)
    implicit none
    real, intent(in) :: hd, se
    fth25 = 1. + exp((-hd + se * t25) / (rgas * t25))
  end function fth25

end module modfarquhar
