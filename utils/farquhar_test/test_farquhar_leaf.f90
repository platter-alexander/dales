! Leaf-level checks of modfarquhar on fixed inputs (no LES).
! Build and run with ./run_test.sh after building DALES in build-debug.
program test_farquhar_leaf
  use modfarquhar, only : farquhar_set_params, farquhar_leaf, farquhar_dales
  implicit none

  real, parameter :: pres = 101325., tfrz = 273.15
  real, parameter :: vcmax25 = 41., btran = 1.
  real :: ca, eair, gb_mol, tleaf, apar, ci, an, gs, ag, cp, esat, resid
  real :: gcc, fleaf, ci_d, co2abs, cp_d, bt, fleaf_check
  integer :: i

  call farquhar_set_params

  ! common conditions: 400 ppm, 25 degC, VPD 1 kPa, gb = 1 mol m-2 s-1
  ca     = 400.e-6 * pres
  tleaf  = 25. + tfrz
  esat   = 0.611e3 * exp(17.2694 * (tleaf - 273.16) / (tleaf - 35.86))
  eair   = esat - 1000.
  gb_mol = 1.e6

  print '(a)', '=== 1. Light response (25 degC, 400 ppm, VPD 1 kPa)'
  print '(a)', '  APAR[W/m2]  An[umol]   gs[mol]   ci/ca    residual[Pa]'
  do i = 0, 8
    apar = 50. * i
    ci   = 0.7 * ca
    call farquhar_leaf(tleaf, pres, ca, eair, apar, gb_mol, vcmax25, btran, ci, an, gs, ag, cp)
    resid = ci - (ca - an * pres * (1.4*gs + 1.6*gb_mol) / (gb_mol*gs))
    print '(f10.1, f10.3, f10.4, f9.3, es14.3)', apar, an, gs*1.e-6, ci/ca, resid
  end do

  print '(a)', '=== 2. Temperature response (APAR 300 W/m2, 400 ppm, VPD 1 kPa)'
  print '(a)', '  Tleaf[C]    An[umol]   gs[mol]   ci/ca'
  do i = 0, 7
    tleaf = 10. + 5.*i + tfrz
    esat  = 0.611e3 * exp(17.2694 * (tleaf - 273.16) / (tleaf - 35.86))
    eair  = esat - 1000.
    ci    = 0.7 * ca
    call farquhar_leaf(tleaf, pres, ca, eair, 300., gb_mol, vcmax25, btran, ci, an, gs, ag, cp)
    print '(f10.1, f10.3, f10.4, f9.3)', tleaf - tfrz, an, gs*1.e-6, ci/ca
  end do

  print '(a)', '=== 3. CO2 response (APAR 300 W/m2, 25 degC, VPD 1 kPa)'
  print '(a)', '  ca[ppm]     An[umol]   gs[mol]   ci[ppm]'
  tleaf = 25. + tfrz
  esat  = 0.611e3 * exp(17.2694 * (tleaf - 273.16) / (tleaf - 35.86))
  eair  = esat - 1000.
  do i = 1, 8
    ca = 100.e-6 * i * pres
    ci = 0.7 * ca
    call farquhar_leaf(tleaf, pres, ca, eair, 300., gb_mol, vcmax25, btran, ci, an, gs, ag, cp)
    print '(f10.0, f10.3, f10.4, f10.1)', ca/pres*1.e6, an, gs*1.e-6, ci/pres*1.e6
  end do

  print '(a)', '=== 4. VPD response (APAR 300 W/m2, 25 degC, 400 ppm)'
  print '(a)', '  VPD[kPa]    An[umol]   gs[mol]   ci/ca'
  ca = 400.e-6 * pres
  do i = 1, 6
    eair = esat - 500. * i
    ci   = 0.7 * ca
    call farquhar_leaf(tleaf, pres, ca, eair, 300., gb_mol, vcmax25, btran, ci, an, gs, ag, cp)
    print '(f10.1, f10.3, f10.4, f9.3)', 0.5*i, an, gs*1.e-6, ci/ca
  end do

  print '(a)', '=== 5. Night (APAR 0) at 20 and 30 degC: expect An = -Rd, gs = g0'
  do i = 0, 1
    tleaf = 20. + 10.*i + tfrz
    ci    = 0.7 * ca
    call farquhar_leaf(tleaf, pres, ca, eair, 0., gb_mol, vcmax25, btran, ci, an, gs, ag, cp)
    print '(a,f5.1,a,f8.4,a,f8.4)', '  T = ', tleaf-tfrz, '  An = ', an, '  gs [mol] = ', gs*1.e-6
  end do

  print '(a)', '=== 6. DALES interface: units and consistency with leafco2f'
  ! DALES-like conditions in the canopy at 99000 Pa
  ci_d = -1.   ! not set yet: interface must fall back to 0.7 ca
  do i = 1, 3
    call farquhar_dales(410.e3, 0.015, 1.15, 99000., 299., 0.3, 100.*i, 0.05, vcmax25, &
                        gcc, fleaf, ci_d, co2abs, cp_d, bt)
    ! flux as DALES computes it in leafco2f (transpiretype = 1)
    fleaf_check = -(co2abs - ci_d) / (1./gcc + 1.4/(1.075*0.05))
    print '(a,f6.0,a,es11.4,a,es11.4,a,f8.4,a,es10.3)', '  APAR=', 100.*i, '  Fleaf[mg/m2/s]=', fleaf, &
          '  leafco2f=', fleaf_check, '  gcc[m/s]=', gcc, '  ci[mg/m3]=', ci_d
  end do
  call farquhar_dales(410.e3, 0.015, 1.15, 99000., 299., 0.3, 0., 0.05, vcmax25, &
                      gcc, fleaf, ci_d, co2abs, cp_d, bt)
  fleaf_check = -(co2abs - ci_d) / (1./gcc + 1.4/(1.075*0.05))
  print '(a,es11.4,a,es11.4,a,f9.6,a,f8.2,a,f8.2)', '  night: Fleaf=', fleaf, '  leafco2f=', fleaf_check, &
        '  gcc=', gcc, '  ci=', ci_d, '  ca=', co2abs

end program test_farquhar_leaf
