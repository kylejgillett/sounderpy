"""
    SOUNDERPY CALC

    Purpose of module:

    Compute sounding parameters from a SounderPy `clean_data` profile: surface-based, mixed-layer,
    and most-unstable parcels (CAPE, CIN, LCL, LFC, EL), entraining CAPE, DCAPE and DCIN, lapse rates,
    the PBL, bulk shear, storm-relative helicity and winds, streamwise vorticity, storm motions,
    the effective inflow layer, and composite indices. `sounding_params(clean_data).calc()` returns
    them in a dict; `vad_params(vad_data).calc()` does the same for VAD wind profiles.

    Docs:
        sounding parameters:   https://sounderpy.app/docs/reference/helpertools.html#returning-a-dictionary-of-profile-parameters
        storm motions:         https://sounderpy.app/docs/reference/plotsettings.html#storm-motions
        parcels:               https://sounderpy.app/docs/reference/plotsettings.html#parcels
        surface modification:  https://sounderpy.app/docs/reference/plotsettings.html#surface-modification

    References:
        Blumberg, W. G., K. T. Halbert, T. A. Supinie, P. T. Marsh, R. L. Thompson, and J. A. Hart, 2017:
            SHARPpy: An open-source sounding analysis toolkit for the atmospheric sciences. Bull. Amer.
            Meteor. Soc., 98, 1625-1636, https://doi.org/10.1175/BAMS-D-15-00309.1 (parcels, indices, and
            storm motions)
        May, R. M., S. C. Arms, P. Marsh, E. Bruning, J. R. Leeman, K. Goebbert, J. E. Thielen, Z. S.
            Bruick, and M. D. Camron, 2023: MetPy: A Python package for meteorological data. Unidata,
            https://doi.org/10.5065/D6WW7G29 (thermodynamic calculations)
        Doswell, C. A., III, and E. N. Rasmussen, 1994: The effect of neglecting the virtual temperature
            correction on CAPE calculations. Wea. Forecasting, 9, 625-629,
            https://doi.org/10.1175/1520-0434(1994)009<0625:TEONTV>2.0.CO;2 (virtual temperature correction
            in CAPE)
        Blanchard, D. O., 1998: Assessing the vertical distribution of convective available potential
            energy. Wea. Forecasting, 13, 870-877,
            https://doi.org/10.1175/1520-0434(1998)013<0870:ATVDOC>2.0.CO;2 (normalized CAPE)
        Peters, J. M., D. R. Chavas, C.-Y. Su, H. Morrison, and B. E. Coffer, 2023: An analytic formula for
            entraining CAPE in midlatitude storm environments. J. Atmos. Sci., 80, 2165-2186,
            https://doi.org/10.1175/JAS-D-23-0003.1; computed with the ecape-parcel package,
            https://github.com/a-urq/ecape-parcel-py (entraining CAPE)
        Emanuel, K. A., 1994: Atmospheric Convection. Oxford University Press, 592 pp. (DCAPE)
        Market, P. S., S. M. Rochette, J. Shewchuk, R. Difani, J. S. Kastman, C. B. Henson, and N. I. Fox,
            2017: Evaluating elevated convection with the downdraft convective inhibition. Atmos. Sci.
            Lett., 18, 76-81, https://doi.org/10.1002/asl.727 (DCIN)
        Gilmore, M. S., and L. J. Wicker, 1998: The influence of midtropospheric dryness on supercell
            morphology and evolution. Mon. Wea. Rev., 126, 943-958,
            https://doi.org/10.1175/1520-0493(1998)126<0943:TIOMDO>2.0.CO;2 (DCAPE)
        Bunkers, M. J., B. A. Klimowski, J. W. Zeitler, R. L. Thompson, and M. L. Weisman, 2000: Predicting
            supercell motion using a new hodograph technique. Wea. Forecasting, 15, 61-79,
            https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2 (Bunkers storm motion)
        Bunkers, M. J., D. A. Barber, R. L. Thompson, R. Edwards, and J. Garner, 2014: Choosing a universal
            mean wind for supercell motion prediction. J. Operational Meteor., 2 (11), 115-129,
            https://doi.org/10.15191/nwajom.2014.0211 (Bunkers storm motion from the effective inflow layer)
        Corfidi, S. F., 2003: Cold pools and MCS propagation: Forecasting the motion of downwind-developing
            MCSs. Wea. Forecasting, 18, 997-1017,
            https://doi.org/10.1175/1520-0434(2003)018<0997:CPAMPF>2.0.CO;2 (Corfidi MCS motion)
        Nixon, C. J., and J. T. Allen, 2021: Anticipating deviant tornado motion using a simple hodograph
            technique. Wea. Forecasting, 36, 219-235, https://doi.org/10.1175/WAF-D-20-0056.1 (deviant
            tornado motion)
        Davies-Jones, R., 1984: Streamwise vorticity: The origin of updraft rotation in supercell storms. J.
            Atmos. Sci., 41, 2991-3006, https://doi.org/10.1175/1520-0469(1984)041<2991:SVTOOU>2.0.CO;2
            (streamwise vorticity)
        Davies-Jones, R., D. W. Burgess, and M. P. Foster, 1990: Test of helicity as a tornado forecast
            parameter. Preprints, 16th Conf. on Severe Local Storms, Kananaskis Park, AB, Canada, Amer.
            Meteor. Soc., 588-592 (SRH)
        Brandt, S., 2022: Storm-Relative Hodograph Plotter (Streamwise_Hodograph_Operational.py). GitHub,
            https://github.com/SamBrandtMeteo/Storm-Relative-Hodograph-Plotter (streamwise vorticity)
        Thompson, R. L., R. Edwards, and C. M. Mead, 2004: An update to the supercell composite and
            significant tornado parameters. Preprints, 22nd Conf. on Severe Local Storms, Hyannis, MA, Amer.
            Meteor. Soc., P8.1, https://www.spc.noaa.gov/publications/thompson/stp_scp.pdf (SCP)
        Thompson, R. L., C. M. Mead, and R. Edwards, 2007: Effective storm-relative helicity and bulk shear
            in supercell thunderstorm environments. Wea. Forecasting, 22, 102-115,
            https://doi.org/10.1175/WAF969.1 (effective inflow layer, effective bulk shear, and effective
            SRH)
        Thompson, R. L., B. T. Smith, J. S. Grams, A. R. Dean, and C. Broyles, 2012: Convective modes for
            significant severe thunderstorms in the contiguous United States. Part II: Supercell and QLCS
            tornado environments. Wea. Forecasting, 27, 1136-1154, https://doi.org/10.1175/WAF-D-11-00116.1
            (STP)
        Hart, J. A., and W. Korotky, 1991: The SHARP Workstation v1.50 Users Guide. NOAA/National Weather
            Service, 30 pp. (EHI)
        Davies, J. M., 1993: Hourly helicity, instability, and EHI in forecasting supercell tornadoes.
            Preprints, 17th Conf. on Severe Local Storms, St. Louis, MO, Amer. Meteor. Soc., 107-111 (EHI)
        Weisman, M. L., and J. B. Klemp, 1982: The dependence of numerically simulated convective storms on
            vertical wind shear and buoyancy. Mon. Wea. Rev., 110, 504-520,
            https://doi.org/10.1175/1520-0493(1982)110<0504:TDONSC>2.0.CO;2 (bulk Richardson number)
        George, J. J., 1960: Weather Forecasting for Aeronautics. Academic Press, 673 pp. (K-index)
        Miller, R. C., 1972: Notes on analysis and severe-storm forecasting procedures of the Air Force
            Global Weather Central. AWS Tech. Rep. 200 (Rev.), Air Weather Service, Scott AFB, IL, 190 pp.
            (Total Totals)
        Storm Prediction Center: Significant Hail Parameter (SHIP).
            https://www.spc.noaa.gov/exper/mesoanalysis/help/help_sigh.html (SHIP)
        Johnson, A. W., and K. E. Sugden, 2014: Evaluation of sounding-derived thermodynamic and
            wind-related parameters associated with large hail events. Electronic J. Severe Storms Meteor.,
            9 (5), 1-42, https://doi.org/10.55599/ejssm.v9i5.57 (SHIP)
        Cohen, A. E., M. C. Coniglio, S. F. Corfidi, and S. J. Corfidi, 2007: Discrimination of mesoscale
            convective system environments using sounding observations. Wea. Forecasting, 22, 1045-1062,
            https://doi.org/10.1175/WAF1040.1 (0-10 km bulk shear and theta-e difference)
        Browning, K. A., and R. Wexler, 1968: The determination of kinematic properties of a wind field
            using Doppler radar. J. Appl. Meteor., 7, 105-113,
            https://doi.org/10.1175/1520-0450(1968)007<0105:TDOKPO>2.0.CO;2 (VAD wind profiles)

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024, 2025, 2026
"""

# sharppy
from .SHARPPYMAIN.sharppy.sharptab.profile import create_profile
from .SHARPPYMAIN.sharppy.sharptab.interp import *
from .SHARPPYMAIN.sharppy.sharptab.winds import *
from .SHARPPYMAIN.sharppy.sharptab.utils import comp2vec, KTS2MS, QC
from .SHARPPYMAIN.sharppy.sharptab.constants import G
from .SHARPPYMAIN.sharppy.sharptab.params import *
from .SHARPPYMAIN.sharppy.sharptab import thermo as sharp_thermo, interp as sharp_interp


try:
    from numpy import trapezoid
except ImportError:
    from numpy import trapz as trapezoid

import numpy as np
import numpy.ma as ma
import warnings
import metpy.calc as mpcalc
import metpy.constants as mpconsts
from metpy.units import units
from ecape_parcel.ecape_calc import calc_ecape

# sounderpy modules
from .utils import modify_surface, find_nearest


# STORM MOTION NAMES
# names accepted for each Bunkers/mean wind storm motion; a custom storm motion is a [direction, speed] list
STORM_MOTION_NAMES = {'rm': ('rm', 'right', 'right_moving', 'right_mover'),
                      'lm': ('lm', 'left', 'left_moving', 'left_mover'),
                      'mw': ('mw', 'mean', 'mean_wind')}

def storm_motion_type(storm_motion):
    # return 'rm', 'lm', 'mw', or 'custom' for a storm motion, or raise an error for an unknown one
    if isinstance(storm_motion, list) and len(storm_motion) == 2:
        return 'custom'
    if isinstance(storm_motion, str):
        for sm_type, names in STORM_MOTION_NAMES.items():
            if storm_motion.casefold() in names:
                return sm_type
    raise ValueError(f"Invalid storm motion: {storm_motion!r}. Use one of {[name for names in STORM_MOTION_NAMES.values() for name in names]}, "
                     f"or a custom storm motion as a [direction, speed] list, ex: [270.0, 25.0]")




# DCAPE/DCIN calculation inherited from MetPy
# https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.downdraft_cape.html
# Based on work by Emanuel, K. A., 1994 & Market et al 2017
# Wet-bulb temperature at every level of a profile (Normand's rule), the same method as MetPy's
# wet_bulb_temperature(): lift each level to its LCL, then follow the moist adiabat (MetPy's moist_lapse()
# equation) back down to the level. All levels are integrated at once with a fixed-step RK4 in ln(p),
# which matches MetPy to within 0.0001 C and is much faster on profiles with many levels
def wet_bulb_profile(pressure, temperature, dewpoint, steps=10):
    lcl_p, lcl_t = mpcalc.lcl(pressure, temperature, dewpoint)
    rd, lv, cp, eps = mpconsts.nounit.Rd, mpconsts.nounit.Lv, mpconsts.nounit.Cp_d, mpconsts.nounit.epsilon

    def dt_dlnp(lnp, t):
        p = np.exp(lnp)
        rs = mpcalc.saturation_mixing_ratio(p * units.Pa, t * units.K).m_as('dimensionless')
        return (rd * t + lv * rs) / (cp + lv * lv * rs * eps / (rd * t ** 2))

    lnp = np.log(lcl_p.m_as('Pa'))
    h = (np.log(pressure.m_as('Pa')) - lnp) / steps
    t = lcl_t.m_as('K')
    for _ in range(steps):
        k1 = dt_dlnp(lnp, t)
        k2 = dt_dlnp(lnp + h / 2, t + h / 2 * k1)
        k3 = dt_dlnp(lnp + h / 2, t + h / 2 * k2)
        k4 = dt_dlnp(lnp + h, t + h * k3)
        t = t + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        lnp = lnp + h
    return (t * units.K).to(temperature.units)



# DCAPE edits and DCIN addition by Kyle Gillett, 2024
def dcape_calc(pressure, temperature, dewpoint):
    # Get layer between 500 and 700 hPa
    p_layer, t_layer, td_layer = mpcalc.get_layer(pressure, temperature, dewpoint,
                                                  bottom=700 * units.hPa,
                                                  depth=200 * units.hPa, interpolate=True)
    theta_e = mpcalc.equivalent_potential_temperature(p_layer, t_layer, td_layer)

    # Find parcel with minimum thetae in the layer
    min_idx = np.argmin(theta_e)
    parcel_start_p = p_layer[min_idx]

    parcel_start_td = td_layer[min_idx]
    parcel_start_wb = mpcalc.wet_bulb_temperature(parcel_start_p, t_layer[min_idx], parcel_start_td)

    # Descend parcel moist adiabatically to surface
    down_pressure = pressure[pressure >= parcel_start_p].to(units.hPa)
    down_parcel_trace = mpcalc.moist_lapse(down_pressure, parcel_start_wb,
                                           reference_pressure=parcel_start_p)

    # Find virtual temperature of parcel and environment
    parcel_virt_temp = mpcalc.virtual_temperature_from_dewpoint(down_pressure, down_parcel_trace,
                                                                down_parcel_trace)
    env_virt_temp = mpcalc.virtual_temperature_from_dewpoint(down_pressure,
                                                             temperature[pressure >= parcel_start_p],
                                                             dewpoint[pressure >= parcel_start_p])

    # calculate total differences (remove units for NumPy)
    total_diff = (env_virt_temp - parcel_virt_temp).to(units.degK).magnitude

    # Where total_diff is zero, there is a level of neutral buoyancy.
    # If a level neutral buoyancy exists (i.e. parcel transitions from
    # cooler than the env. to warmer than the env) then downward momentum
    # transfer via DCAPE is suppressed by DCIN.

    # If there is no value of zero or no sign change in total_diff,
    # then there is no level of neutral buoyancy and the positive area
    # (DCAPE) must extend to the sfc.

    # `locate_bounds()` attempts to identify this

    def locate_bounds(arr):
        # check for value of zero in array
        # if one exists, return its index
        if 0 in arr:
            return np.where(arr == 0)[0][0]

            # else, check for sign changes
        # if one is found, return the index
        # of the value closest to zero
        for i in range(1, len(arr)):
            if np.sign(arr[i - 1]) != np.sign(arr[i]):
                return i if abs(arr[i]) < abs(arr[i - 1]) else i - 1

        # if neither are true, return the appropriate
        # bound by the sign of the array
        if np.all(arr > 0):
            return 0
        elif np.all(arr < 0):
            return -1

    parcel_bound = locate_bounds(total_diff)

    # DCAPE: Integral from top (p_nb) to LNB (index 0)
    cape_diff = (env_virt_temp[parcel_bound:-1] - parcel_virt_temp[parcel_bound:-1]).to(units.degK).magnitude
    cape_lnp = np.log(down_pressure[parcel_bound:-1].magnitude)

    # DCIN: Integral from surface (index 0) to LNB (index p_nb)
    cin_diff = (env_virt_temp[0:parcel_bound] - parcel_virt_temp[0:parcel_bound]).to(units.degK).magnitude
    cin_lnp = np.log(down_pressure[0:parcel_bound].magnitude)

    # Find DCAPE
    dcape = -((287.04749097718457 * units('joule/(kelvin kilogram)'))
              * units.Quantity(trapezoid(cape_diff, cape_lnp), 'K')
              ).to(units('J/kg'))

    # Find DCIN
    dcin = -((287.04749097718457 * units('joule/(kelvin kilogram)'))
             * units.Quantity(trapezoid(cin_diff, cin_lnp), 'K')
             ).to(units('J/kg'))

    return dcape, dcin, down_pressure, down_parcel_trace



# CIN as the negative area between a parcel's starting level and its LFC.
# SHARPpy's parcelx() reports `bminus` as every negative layer below 500 hPa, including the stable air
# above the EL (it saves the CIN at the LFC as `cinh_old`, but the line that would use it is commented
# out). When the EL is below 500 hPa, as with a few J/kg of CAPE, that inflates CIN to -1000 J/kg or more.
# This recomputes it the intended way: dry ascent at constant mixing ratio to the LCL, then moist
# ascent to the LFC, using virtual temperature for the parcel and the environment
def cin_below_lfc(prof, pcl, pres, tmpc, dwpc, dp=-2.0):
    cape = pcl.bplus
    if ma.is_masked(cape) or np.floor(cape) == 0:
        return pcl.bminus  # SHARPpy already sets CIN to 0 without CAPE
    lfc_p = pcl.lfcpres
    if ma.is_masked(lfc_p) or lfc_p is None or not np.isfinite(lfc_p) or lfc_p > pres:
        return pcl.bminus

    lcl_p, lcl_t = sharp_thermo.drylift(pres, tmpc, dwpc)
    mixr = sharp_thermo.mixratio(pres, dwpc)
    theta_k = sharp_thermo.theta(pres, tmpc, 1000.) + 273.15

    p_levels = np.append(np.arange(pres, lfc_p, dp), lfc_p)
    if len(p_levels) < 2:
        return 0.
    tv_parcel = np.empty(len(p_levels))
    for i, p in enumerate(p_levels):
        if p >= lcl_p:
            t = theta_k * (p / 1000.) ** 0.2857 - 273.15
            tv_parcel[i] = sharp_thermo.virtemp(p, t, sharp_thermo.temp_at_mixrat(mixr, p))
        else:
            t = sharp_thermo.wetlift(lcl_p, lcl_t, p)
            tv_parcel[i] = sharp_thermo.virtemp(p, t, t)
    tv_env = np.ma.filled(np.ma.asarray(sharp_interp.vtmp(prof, p_levels), dtype=float), np.nan)
    hght = np.ma.filled(np.ma.asarray(sharp_interp.hght(prof, p_levels), dtype=float), np.nan)

    buoyancy = (tv_parcel - tv_env) / (tv_env + 273.15)
    layer_energy = 9.80665 * (buoyancy[:-1] + buoyancy[1:]) / 2 * np.diff(hght)
    layer_energy = layer_energy[np.isfinite(layer_energy)]
    return float(layer_energy[layer_energy < 0].sum())


# Lift a parcel with SHARPpy, or return an empty one (every value masked) when the profile can't define it.
# Sparse historical soundings (IGRAv2 before ~1970) often have no dewpoints in the lowest levels, or none
# at all; SHARPpy's parcel search then fails (e.g. most_unstable_level finds no levels and raises IndexError)
def lift_parcel(prof, flag, name, pres, tmpc, dwpc):
    start = [pres, tmpc, dwpc]
    if any(ma.is_masked(v) or not np.isfinite(float(v)) for v in start):
        warnings.warn(f"The {name} parcel could not be computed for this sounding (no moisture data at its level)", Warning)
        return Parcel()
    try:
        return parcelx(prof, flag=flag, pres=pres, tmpc=tmpc, dwpc=dwpc)
    except (IndexError, ValueError):
        warnings.warn(f"The {name} parcel could not be computed for this sounding (not enough data)", Warning)
        return Parcel()








# CAPE, CIN and 0-3 km CAPE (J/kg) of a parcel lifted from (pres, tmpc, dwpc), for the stepwise CAPE/CIN
# profiles. These are the same values as SHARPpy's parcelx() `bplus`, `bminus` and `b3km`, but without the
# levels and indices parcelx() also finds (LFC, EL, MPL, lifted indices, cap strength, BRN, and the
# freezing-level CAPEs), so lifting a parcel from every level is about 2-5x faster. Values it can't
# compute are nan
def stepwise_parcel(prof, pres, tmpc, dwpc, dp=-1):
    if prof.pres.compressed().shape[0] < 1:
        return np.nan, np.nan, np.nan
    pbot = prof.pres[prof.sfc]
    ptop = prof.pres[prof.pres.shape[0]-1]
    if pbot > pres:
        pbot = pres

    # dry ascent to the LCL
    pe2, tp2 = sharp_thermo.drylift(pres, tmpc, dwpc)
    if type(pe2) == type(ma.masked) or np.isnan(pe2):
        return np.nan, np.nan, np.nan
    lcl_hght = sharp_interp.to_agl(prof, sharp_interp.hght(prof, pe2))
    theta_parcel = sharp_thermo.theta(pe2, tp2, 1000.)
    mixr = sharp_thermo.mixratio(pres, dwpc)

    # CIN below the LCL, in dp steps, with virtual temperature
    pp = np.arange(pbot, pe2+dp, dp, dtype=type(pbot))
    hh = sharp_interp.hght(prof, pp)
    tv_env = sharp_thermo.virtemp(pp, sharp_thermo.theta(pp, sharp_interp.temp(prof, pp), 1000.), sharp_interp.dwpt(prof, pp))
    tdef = (sharp_thermo.virtemp(pp, theta_parcel, sharp_thermo.temp_at_mixrat(mixr, pp)) - tv_env) / sharp_thermo.ctok(tv_env)
    lyre = G * (tdef[:-1]+tdef[1:]) / 2 * (hh[1:]-hh[:-1])
    totn = lyre[lyre < 0].sum()
    if not totn:
        totn = 0.
    if pbot > pe2:
        pbot = pe2
    if pbot < prof.pres[-1]:
        # the LCL is above the top of the data
        return np.nan, np.nan, np.nan
    lptr = ma.where(pbot >= prof.pres)[0].min()
    uptr = ma.where(ptop <= prof.pres)[0].max()

    # moist ascent from the LCL, one data level at a time
    totp = 0.
    bplus = bminus = b3km = None
    pe1 = pbot
    h1 = sharp_interp.hght(prof, pe1)
    te1 = sharp_interp.vtmp(prof, pe1)
    tp1 = sharp_thermo.wetlift(pe2, tp2, pe1)
    lyre = 0
    h3km = sharp_interp.to_msl(prof, 3000.)
    for i in range(lptr, prof.pres.shape[0]):
        if not QC(prof.tmpc[i]):
            continue
        pe2 = prof.pres[i]
        h2 = prof.hght[i]
        te2 = prof.vtmp[i]
        tp2 = sharp_thermo.wetlift(pe1, tp1, pe2)
        tdef1 = (sharp_thermo.virtemp(pe1, tp1, tp1) - te1) / sharp_thermo.ctok(te1)
        tdef2 = (sharp_thermo.virtemp(pe2, tp2, tp2) - te2) / sharp_thermo.ctok(te2)
        lyre = G * (tdef1 + tdef2) / 2. * (h2 - h1)
        if lyre > 0:
            totp += lyre
        elif pe2 > 500.:
            totn += lyre  # CIN counts negative layers below 500 hPa
        pelast = pe1
        pe1, te1, tp1 = pe2, te2, tp2

        # top of the data
        if i >= uptr and bplus is None:
            pe3, h3, te3, tp3 = pe1, h2, te1, tp1
            if lyre > 0:
                bplus = totp - lyre
                bminus = totn
            else:
                bplus = totp
                bminus = totn + lyre if pe2 > 500. else totn
            pe2 = ptop
            h2 = sharp_interp.hght(prof, pe2)
            te2 = sharp_interp.vtmp(prof, pe2)
            tp2 = sharp_thermo.wetlift(pe3, tp3, pe2)
            tdef3 = (sharp_thermo.virtemp(pe3, tp3, tp3) - te3) / sharp_thermo.ctok(te3)
            tdef2 = (sharp_thermo.virtemp(pe2, tp2, tp2) - te2) / sharp_thermo.ctok(te2)
            lyrf = G * (tdef3 + tdef2) / 2. * (h2 - h3)
            if lyrf > 0:
                bplus += lyrf
            elif pe2 > 500.:
                bminus += lyrf
            if bplus == 0:
                bminus = 0.

        # CAPE up to 3 km AGL
        if lcl_hght < 3000.:
            if sharp_interp.to_agl(prof, h1) <= 3000. and sharp_interp.to_agl(prof, h2) >= 3000. and b3km is None:
                pe3 = pelast
                h3 = sharp_interp.hght(prof, pe3)
                te3 = sharp_interp.vtmp(prof, pe3)
                tp3 = sharp_thermo.wetlift(pe1, tp1, pe3)
                b3km = totp - lyre if lyre > 0 else totp
                pe4 = sharp_interp.pres(prof, h3km)
                if QC(pe2):
                    te2 = sharp_interp.vtmp(prof, pe4)
                    tp2 = sharp_thermo.wetlift(pe3, tp3, pe4)
                    tdef3 = (sharp_thermo.virtemp(pe3, tp3, tp3) - te3) / sharp_thermo.ctok(te3)
                    tdef2 = (sharp_thermo.virtemp(pe4, tp2, tp2) - te2) / sharp_thermo.ctok(te2)
                    lyrf = G * (tdef3 + tdef2) / 2. * (h3km - h3)
                    if lyrf > 0:
                        b3km += lyrf
        else:
            b3km = 0.
        h1 = h2

    if bplus is None:
        bplus = totp
    if np.floor(bplus) == 0:
        bminus = 0.  # no CIN without CAPE
    return bplus, bminus, (np.nan if b3km is None else b3km)



class sounding_params:


    def __init__(self, clean_data, storm_motion='right_moving', modify_sfc=None, include_all_parcels=False, parcel_depth=None):
            self.clean_data = clean_data 
            self.storm_motion = storm_motion
            self.sm_type = storm_motion_type(storm_motion)
            self.modify_sfc = modify_sfc
            self.include_all_parcels = include_all_parcels
            self.parcel_depth = parcel_depth

            
    ######################################################################################################
    ### CALC FUNCTION -- COMPUTE VALUES ###
    ######################################################################################################        
    def calc(self):
        

    

        general = {}
        kinem  = {} 
        thermo = {}
        intrp  = {}

        # declare easy variable names for reuse from `clean_data`
        # first, check for correct units: 
        if self.clean_data['T'].units != 'degree_Celsius':
            self.clean_data['T'] =  self.clean_data['T'].to(units.degC)
        if self.clean_data['Td'].units != 'degree_Celsius':
            self.clean_data['Td'] = self.clean_data['Td'].to(units.degC)





        #################################################################
        ### SET UP THE DATA ###
        #################################################################
        # SFC CORRECTION
        if str(type(self.modify_sfc)) == "<class 'dict'>":
            sounding_data = modify_surface(self.clean_data, self.modify_sfc)
            modify_sfc_txt = True
        else:
            sounding_data = self.clean_data
            modify_sfc_txt = False

        # declare easy variable names for reuse from `clean_data`
        T = sounding_data['T']
        Td = sounding_data['Td']
        p = sounding_data['p']
        z = sounding_data['z']
        u = sounding_data['u']
        v = sounding_data['v']
        wd = mpcalc.wind_direction(u, v)
        ws = mpcalc.wind_speed(u, v)




        # compute useful variables and add them to a new dict of 
        # data that this function will return 
        if sounding_data['site_info']['site-elv'] == 9999:
            general['elevation'] = z[0].m
        else:
            general['elevation'] = int(sounding_data['site_info']['site-elv'])


        general['sfc_pressure']  = (p[0].m*(1-(0.0065*(general['elevation']))/(T[0].m+(0.0065*(general['elevation']))+273.15))**-5.257)*units.hPa
        general['wet_bulb']      = wet_bulb_profile(p, T, Td)
        general['rel_humidity']  = mpcalc.relative_humidity_from_dewpoint(T, Td)*100
        general['spec_humidity'] = (mpcalc.specific_humidity_from_dewpoint(p, Td)*1000)*units.g/units.kg
        general['mix_ratio']     = mpcalc.mixing_ratio_from_specific_humidity(general['spec_humidity'])
        general['virt_temp']     = mpcalc.virtual_temperature(T, general['mix_ratio'])
        general['theta']         = mpcalc.potential_temperature(p, T)
        general['theta_e']       = mpcalc.equivalent_potential_temperature(p, T, Td)
        try:
            general['pwat']      = np.round((mpcalc.precipitable_water(p, Td).to('in')).m, 3)
            if general['pwat'] < 0:
                general['pwat'] = ma.masked
        except:
            general['pwat'] = ma.masked
        ##########################################   



        def interpolate(var,hgts,step):
            levels=np.arange(0,np.max(hgts),step)
            varinterp=np.zeros(len(levels))
            for i in range(0,len(levels)):
                try:
                    lower=np.where(hgts-levels[i]<=0,hgts-levels[i],-np.inf).argmax()
                    varinterp[i]=(((var[lower+1]-var[lower])/(hgts[lower+1]-hgts[lower]))*(levels[i]-hgts[lower])+var[lower])
                except IndexError:
                    # Handle the case where the index is out of bounds
                    varinterp[i] = np.nan 
            return varinterp 

        resolution=100

        rhINTRP     = general['rel_humidity'].m[np.isnan(z)==False]
        mrINTRP     = general['mix_ratio'].m[np.isnan(z)==False]
        thetaINTRP  = general['theta'].m[np.isnan(z)==False]
        thetaeINTRP = general['theta_e'].m[np.isnan(z)==False]
        wbINTRP     = general['wet_bulb'].m[np.isnan(z)==False]
        pINTRP  = p.m[np.isnan(z)==False]
        tINTRP  = T.m[np.isnan(z)==False]
        uINTRP  = u.m[np.isnan(z)==False]
        vINTRP  = v.m[np.isnan(z)==False]
        zINTRP  = z.m[np.isnan(z)==False]
        zINTRP  = zINTRP-zINTRP[0]
        # Interpolation, add to intrp dict
        intrp['rhINTRP'] = interpolate(rhINTRP,zINTRP,resolution)
        intrp['tINTRP']  = interpolate(tINTRP,zINTRP,resolution)
        intrp['mrINTRP'] = interpolate(mrINTRP,zINTRP,resolution)
        intrp['thetaINTRP']  = interpolate(thetaINTRP,zINTRP,resolution)
        intrp['thetaeINTRP'] = interpolate(thetaeINTRP,zINTRP,resolution)
        intrp['wbINTRP'] = interpolate(wbINTRP,zINTRP,resolution)
        intrp['pINTRP'] = interpolate(pINTRP,zINTRP,resolution)
        intrp['uINTRP'] = interpolate(uINTRP,zINTRP,resolution)
        intrp['vINTRP'] = interpolate(vINTRP,zINTRP,resolution)
        intrp['zINTRP'] = interpolate(zINTRP,zINTRP,resolution)

        hgt_lvls = {}
        hgt_var_list = ['h0', 'h05', 'h1', 'h15', 'h2', 'h25', 
                        'h3', 'h35', 'h4', 'h45', 'h5', 'h55', 
                        'h6', 'h65', 'h7', 'h75', 'h8', 'h85', 
                        'h9', 'h10', 'h11', 'h12', 'h13', 'h14']
        hgt_val_list = [0.,500.,1000.,1500.,2000.,2500.,3000.,3500.,
                        4000.,4500.,5000.,5500.,6000.,6500.,7000.,
                        7500.,8000.,8500.,9000.,10000.,11000.,
                        12000.,13000.,14000.]

        for var, val, in zip(hgt_var_list, hgt_val_list):
            try:
                hgt_lvls[var] = np.where(intrp['zINTRP']==val)[0][0]
            except:
                pass
        intrp['hgt_lvls'] = hgt_lvls


        # FREEZING POINT CALCULATION---           
        try: 
            frz_pt_index = intrp['tINTRP'].tolist().index(list(filter(lambda i: i <= 0, intrp['tINTRP'].tolist()))[0])
            general['frz_pt_p'] = intrp['pINTRP'][frz_pt_index]*units.hPa
            general['frz_pt_z'] = intrp['zINTRP'][frz_pt_index]*units.m
        except IndexError:
            general['frz_pt_p'] = ma.masked
            general['frz_pt_z'] = ma.masked
            warnings.warn("This sounding does not have a freezing point (not enough data)", Warning)
            pass

        try: 
            wb_frz_pt_index = intrp['wbINTRP'].tolist().index(list(filter(lambda i: i <= 0, intrp['wbINTRP'].tolist()))[0])
            general['wb_frz_pt_p'] = intrp['pINTRP'][wb_frz_pt_index]*units.hPa
            general['wb_frz_pt_z'] = intrp['zINTRP'][wb_frz_pt_index]*units.m
        except IndexError:
            general['wb_frz_pt_p'] = ma.masked
            general['wb_frz_pt_z'] = ma.masked
            warnings.warn("This sounding does not have a wet bulb freezing point (not enough data)", Warning)
            pass


        ###################################################################
        ### SHARPPY VARIABLES CALC ###
        ###################################################################
        # create sharppy profile object
        prof = create_profile(profile='default',pres=p.m, hght=z.m, tmpc=T.m, dwpc=Td.m, wspd=ws, wdir=wd, 
                                      missing=-9999, strictQC=False)

        ####################################################################
        ### THERMODYNAMIC VALUES ###
        ####################################################################
        #--- SURFACE BASED PARCEL PROPERTIES ---#
        # ---------------------------------------------------------------
        sbpcl = lift_parcel(prof, 1, 'surface-based', p[0].m, T[0].m, Td[0].m)
        thermo['sbT_trace'] = sbpcl.ttrace
        thermo['sbP_trace'] = sbpcl.ptrace
        thermo['sbZ_trace'] = to_agl(prof, hght(prof, thermo['sbP_trace']))  # m AGL, from the profile's own heights
        # compute heights
        thermo['sb_lcl_p'] = sbpcl.lclpres
        thermo['sb_lcl_z'] = sbpcl.lclhght 
        thermo['sb_lfc_p'] = sbpcl.lfcpres
        thermo['sb_lfc_z'] = sbpcl.lfchght
        thermo['sb_el_p']  = sbpcl.elpres  
        thermo['sb_el_z']  = sbpcl.elhght
        thermo['sb_el_T']  = temp(prof, thermo['sb_el_p'])
        thermo['sb_mpl_p'] = sbpcl.mplpres
        thermo['sbcape']   = sbpcl.bplus
        thermo['sbcin']    = cin_below_lfc(prof, sbpcl, p[0].m, T[0].m, Td[0].m)
        thermo['sb3cape']  = sbpcl.b3km
        thermo['sb6cape']  = sbpcl.b6km

        #--- MOST UNSTABLE PARCEL PROPERTIES ---#
        # ---------------------------------------------------------------
        # MU & ML parcels come from the levels with a temperature and dewpoint: MetPy's search picks a level
        # with a missing dewpoint (argmax of nan), and old soundings often have no dewpoint aloft
        has_moist = np.isfinite(p.m) & np.isfinite(T.m) & np.isfinite(Td.m)
        if has_moist.sum() >= 2:
            mupcl_p, mupcl_T, mupcl_Td = mpcalc.most_unstable_parcel(p[has_moist], T[has_moist], Td[has_moist])[0:3]
        else:
            mupcl_p, mupcl_T, mupcl_Td = np.nan*units.hPa, np.nan*units.degC, np.nan*units.degC
        mupcl = lift_parcel(prof, 3, 'most-unstable', mupcl_p.m, mupcl_T.m, mupcl_Td.m)
        thermo['muT_trace'] = mupcl.ttrace
        thermo['muP_trace'] = mupcl.ptrace
        thermo['muZ_trace'] = to_agl(prof, hght(prof, thermo['muP_trace']))  # m AGL, from the profile's own heights
        # compute heights
        thermo['mu_lcl_p'] = mupcl.lclpres 
        thermo['mu_lcl_z'] = mupcl.lclhght 
        thermo['mu_lfc_p'] = mupcl.lfcpres 
        thermo['mu_lfc_z'] = mupcl.lfchght
        thermo['mu_el_p']  = mupcl.elpres
        thermo['mu_el_z']  = mupcl.elhght
        thermo['mu_el_T']  = temp(prof, thermo['mu_el_p'])
        thermo['mu_mpl_p'] = mupcl.mplpres
        thermo['mucape']   = mupcl.bplus
        thermo['mucin']    = cin_below_lfc(prof, mupcl, mupcl_p.m, mupcl_T.m, mupcl_Td.m)
        thermo['mu3cape']  = mupcl.b3km
        thermo['mu6cape']  = mupcl.b6km

        #--- MIXED LAYER PARCEL PROPERTIES ---#
        # ---------------------------------------------------------------
        if has_moist.sum() >= 2:
            mlpcl_p, mlpcl_T, mlpcl_Td = mpcalc.mixed_parcel(p[has_moist], T[has_moist], Td[has_moist], bottom=p[0], depth=100*units.hPa, interpolate=True)[0:3]
        else:
            mlpcl_p, mlpcl_T, mlpcl_Td = np.nan*units.hPa, np.nan*units.degC, np.nan*units.degC
        mlpcl = lift_parcel(prof, 5, 'mixed-layer', mlpcl_p.m, mlpcl_T.m, mlpcl_Td.m)
        thermo['mlT_trace'] = mlpcl.ttrace
        thermo['mlP_trace'] = mlpcl.ptrace
        # compute heights
        thermo['ml_lcl_p'] = mlpcl.lclpres
        thermo['ml_lcl_z'] = mlpcl.lclhght
        thermo['ml_lfc_p'] = mlpcl.lfcpres
        thermo['ml_lfc_z'] = mlpcl.lfchght
        thermo['ml_el_p']  = mlpcl.elpres
        thermo['ml_el_z']  = mlpcl.elhght
        thermo['ml_el_T']  = temp(prof, thermo['ml_el_p'])
        thermo['ml_mpl_p'] = mlpcl.mplpres
        thermo['mlcape']   = mlpcl.bplus
        thermo['mlcin']    = cin_below_lfc(prof, mlpcl, mlpcl_p.m, mlpcl_T.m, mlpcl_Td.m)
        thermo['ml3cape']  = mlpcl.b3km
        thermo['ml6cape']  = mlpcl.b6km

        #--- NORMALIZED CAPE ---#
        # ---------------------------------------------------------------
        thermo['mu_ncape'] = np.round(thermo['mucape'] / ((thermo['mu_el_z'] - thermo['mu_lfc_z'])),1)
        thermo['sb_ncape'] = np.round(thermo['sbcape'] / ((thermo['sb_el_z'] - thermo['sb_lfc_z'])),1)
        thermo['ml_ncape'] = np.round(thermo['mlcape'] / ((thermo['ml_el_z'] - thermo['ml_lfc_z'])),1) 

        #--- DOWNDRAFT CAPE ---#
        # ---------------------------------------------------------------

        if not np.any(np.isnan(sounding_data['Td'])):
            try:
                thermo['dcape'], thermo['dcin'], thermo['dparcel_p'], thermo['dparcel_T'] = dcape_calc(sounding_data['p'],
                                                                                                   sounding_data['T'],
                                                                                                   sounding_data['Td'])
            except ValueError:
                # skip calculation if there is insufficient data (i.e., most ACARS profiles)
                thermo['dcape'] = ma.masked
                thermo['dcin'] = ma.masked
                thermo['dparcel_p'] = ma.masked
                thermo['dparcel_T'] = ma.masked
                warnings.warn("This sounding doesn't have enough data to compute DCAPE & DCIN", Warning)
                pass
        else:
            thermo['dcape'] = ma.masked
            thermo['dcin'] = ma.masked
            thermo['dparcel_p'] = ma.masked
            thermo['dparcel_T'] = ma.masked
            warnings.warn("DCAPE & DCIN could not be calculated for this sounding", Warning)

        #--- PARCEL PROFILES ---#
        # ---------------------------------------------------------------
        if self.include_all_parcels:
            # Calculate parcels. With a `parcel_depth` (m AGL), only levels up to that depth
            # (plus two levels above it, for interpolation) are lifted. Levels above it are nan
            n_parcels = len(p)
            if self.parcel_depth is not None:
                above_depth = np.nonzero((z.m - z.m[0]) >= self.parcel_depth)[0]
                if above_depth.size > 0:
                    n_parcels = min(above_depth[0] + 2, len(p))
            parcels = np.array([stepwise_parcel(prof, p[i].m, T[i].m, Td[i].m) for i in range(n_parcels)], dtype=float).reshape(-1, 3)
            for col, key in enumerate(('cape_profile', 'cin_profile', '3cape_profile')):
                thermo[key] = np.full(len(p), np.nan)
                thermo[key][:n_parcels] = parcels[:, col]
            # Interpolate
            intrp['cape_profileINTRP'] = interpolate(thermo['cape_profile'], zINTRP, resolution)
            intrp['cin_profileINTRP'] = interpolate(thermo['cin_profile'], zINTRP, resolution)
            intrp['3cape_profileINTRP'] = interpolate(thermo['3cape_profile'], zINTRP, resolution)


        # entraining CAPE is computed after the kinematics, so it can use the storm motion
        
        
        

        
        #--- PBL ---#
        # ---------------------------------------------------------------
        thermo['pbl_top'] = pbl_top(prof) 

        # depth (m AGL) and lapse rate (positive when temperature falls with height) of SHARPpy's PBL, from the
        # data levels. Its mean wind is under STORM MOTION
        pbl_z = z.m - z.m[0]
        if not ma.is_masked(thermo['pbl_top']) and np.isfinite(thermo['pbl_top']):
            thermo['pbl_depth'] = float(np.interp(-float(thermo['pbl_top']), -p.m, pbl_z))
        else:
            thermo['pbl_depth'] = ma.masked
        if not ma.is_masked(thermo['pbl_depth']) and thermo['pbl_depth'] > 0:
            thermo['pbl_lapse_rate'] = (T.m[0] - float(np.interp(thermo['pbl_depth'], pbl_z, T.m))) / (thermo['pbl_depth'] / 1000)
        else:
            thermo['pbl_lapse_rate'] = ma.masked

        #--- DGZ ---#
        # ---------------------------------------------------------------
        thermo['dgz'] = dgz(prof)

        #--- HGZ ---#
        # ---------------------------------------------------------------
        thermo['hgz'] = hgz(prof)

        #--- LAPSE RATES ---#
        # ---------------------------------------------------------------
        thermo['lr_36km'] = lapse_rate(prof, 3000, 6000, pres=False)
        thermo['lr_03km'] = lapse_rate(prof, 0, 3000, pres=False)
        thermo['lr_75']   = lapse_rate(prof, 700, 500, pres=True)
        thermo['lr_max']  = max_lapse_rate(prof, lower=0, upper=6000,
                                         interval=250, depth=2000)

        #--- THETA-E DIFFERENCE ---#
        # ---------------------------------------------------------------
        # the minimum theta-e between the level of maximum theta-e in the lowest 3 km and 7 km, minus that
        # maximum (K, negative); larger differences favor severe MCS winds (Cohen et al. 2007)
        low = (intrp['zINTRP'] <= 3000) & np.isfinite(intrp['thetaeINTRP'])
        if np.nanmax(intrp['zINTRP']) >= 3000 and low.any():
            max_idx = np.nanargmax(np.where(low, intrp['thetaeINTRP'], np.nan))
            mid = (intrp['zINTRP'] >= intrp['zINTRP'][max_idx]) & (intrp['zINTRP'] <= 7000)
            thermo['thetae_diff'] = np.nanmin(intrp['thetaeINTRP'][mid]) - intrp['thetaeINTRP'][max_idx]
        else:
            thermo['thetae_diff'] = ma.masked

        #--- K-INDEX & TOTAL TOTALS ---#
        # ---------------------------------------------------------------
        # both use 850 hPa, so they are masked where 850 hPa is below ground
        thermo['k_index']  = k_index(prof)
        thermo['t_totals'] = t_totals(prof)

        #--- BULK RICHARDSON NUMBER ---#
        # ---------------------------------------------------------------
        # SHARPpy's parcelx computes BRN for each parcel: CAPE / (0.5 * (0-6 km mean wind - 0-500 m mean wind)^2)
        # (infinite when there is no shear between the two mean winds, so treat that as missing)
        def finite_or_masked(value):
            return value if not ma.is_masked(value) and np.isfinite(value) else ma.masked
        thermo['sb_brn'] = finite_or_masked(sbpcl.brn)
        thermo['mu_brn'] = finite_or_masked(mupcl.brn)
        # SHARPpy's bulk_rich() returns masked for the ML parcel (an inverted QC check in its non-surface-based
        # branch). The shear term is the same surface-to-6 km term for every parcel here, so reuse the SB one
        thermo['ml_brn'] = finite_or_masked(ma.true_divide(mlpcl.bplus, sbpcl.brnshear))

        #--- SIGNIFICANT HAIL PARAMETER ---#
        # ---------------------------------------------------------------
        # SHARPpy's ship() uses the freezing level height MSL by default; SPC's 2400 m freezing level
        # threshold is AGL, so pass the AGL height
        try:
            thermo['ship'] = ship(prof, mupcl=mupcl, lr75=thermo['lr_75'],
                                  frz_lvl=to_agl(prof, hght(prof, temp_lvl(prof, 0))))
        except Exception:
            thermo['ship'] = ma.masked
        
        #--- TEMPERATURE ADVECTION ---#
        # ---------------------------------------------------------------
        #returns temp_adv (C/hr) array and 2D array of top and bottom bounds of temp advection layer (mb)
        thermo['temp_adv'] = inferred_temp_adv(prof, lat=sounding_data['site_info']['site-latlon'][0])
        ####################################################################
        
        



        ####################################################################
        ### KINEMATIC VALUES ###
        ####################################################################
        #--- EFFECTIVE INFLOW LAYER ---#
        # ---------------------------------------------------------------   
        kinem['eil'] = effective_inflow_layer(prof, ecape=100, ecinh=-150)
        kinem['eil_z'] = [to_agl(prof, hght(prof, kinem['eil'][0])), to_agl(prof, hght(prof, kinem['eil'][1]))]
        kinem['eil_depth'] = kinem['eil_z'][1] - kinem['eil_z'][0]
        
        eil_idx = [find_nearest(intrp['zINTRP'], kinem['eil_z'][0]), find_nearest(intrp['zINTRP'], kinem['eil_z'][1])]



        #--- STORM MOTION ---#
        # ---------------------------------------------------------------   
        bunkers = bunkers_storm_motion(prof)
        sm_rm, sm_lm = bunkers[0:2], bunkers[2:4]
        sm_mw = mean_wind(prof, pbot=850, ptop=250)

        if ma.is_masked(sm_rm[0]) == True:
            try:
                kinem['sm_rm'], kinem['sm_lm'], kinem['sm_mw'] = mpcalc.bunkers_storm_motion(p, u, v, z)
                kinem['sm_rm'], kinem['sm_lm'], kinem['sm_mw'] = kinem['sm_rm'].m, kinem['sm_lm'].m, kinem['sm_mw'].m
            except:
                kinem['sm_rm'], kinem['sm_lm'], kinem['sm_mw'] = ma.masked, ma.masked, ma.masked
                warnings.warn("Bunkers Storm Motion could not be computed for this sounding (not enough data)", Warning)
                pass
        elif np.isnan(sm_rm[0]) == True:
            kinem['sm_rm'], kinem['sm_lm'], kinem['sm_mw'] = ma.masked, ma.masked, ma.masked
            warnings.warn("Bunkers Storm Motion could not be computed for this sounding (not enough data)", Warning)
        else:
            kinem['sm_rm'], kinem['sm_lm'], kinem['sm_mw'] = sm_rm, sm_lm, sm_mw

        #--- PBL MEAN WIND ---#
        # ---------------------------------------------------------------
        # surface to the PBL top (kt), on a 50-point height grid as the PBL can be shallower than the 100 m
        # interpolated grid: the transport wind that boundary-layer moisture, smoke, and blowing snow follow.
        # Only levels with a wind are used (a radiosonde's GPS winds can start after launch); masked when the
        # PBL has none
        kinem['mean_wind_pbl'] = (ma.masked, ma.masked)
        if not ma.is_masked(thermo['pbl_depth']) and thermo['pbl_depth'] > 0:
            u_kt, v_kt = np.asarray(u.to('kt').m, dtype=float), np.asarray(v.to('kt').m, dtype=float)
            has_wind = np.isfinite(u_kt) & np.isfinite(v_kt) & np.isfinite(pbl_z)
            z_pbl = np.linspace(0, thermo['pbl_depth'], 50)
            if has_wind.sum() > 1:
                z_pbl = z_pbl[(z_pbl >= pbl_z[has_wind].min()) & (z_pbl <= pbl_z[has_wind].max())]
            if has_wind.sum() > 1 and z_pbl.size:
                kinem['mean_wind_pbl'] = (float(np.mean(np.interp(z_pbl, pbl_z[has_wind], u_kt[has_wind]))),
                                          float(np.mean(np.interp(z_pbl, pbl_z[has_wind], v_kt[has_wind]))))
            




        #--- USER DEFINED SM PT ---#
        # ---------------------------------------------------------------  
        if self.sm_type == 'custom':
            kinem['sm_u'], kinem['sm_v'] = mpcalc.wind_components(self.storm_motion[1]*units.kts, self.storm_motion[0]*units.deg)
            kinem['sm_u'], kinem['sm_v'] = kinem['sm_u'].m, kinem['sm_v'].m
            
        else:
            
            if ma.is_masked(kinem['sm_rm']) == False:
                kinem['sm_u'] = kinem[f'sm_{self.sm_type}'][0]
                kinem['sm_v'] = kinem[f'sm_{self.sm_type}'][1]

            else:
                kinem['sm_u'] = ma.masked
                kinem['sm_v'] = ma.masked
            



        #--- DEVIANT TORNADO MOTION ---#
        # --------------------------------------------------------------- 
        if ma.is_masked(kinem['sm_u']) == False:
            p300m = pres(prof, to_msl(prof, 300.))
            sfc = prof.pres[prof.sfc]
            mw_u_300, mw_v_300 = mean_wind(prof, pbot=sfc, ptop=p300m)
            kinem['dtm'] = ((kinem['sm_u']+mw_u_300)/2, (kinem['sm_v']+mw_v_300)/2)
        else:
            kinem['dtm'] = ma.masked
          



        #--- CORFIDI MCS MOTION ---#
        #----------------------------------------------------------------
        kinem['mcs'] = corfidi_mcs_motion(prof)





        #--- BULK SHEAR ---#
        # ---------------------------------------------------------------    
        sfc = prof.pres[prof.sfc]
        p500m = pres(prof, to_msl(prof, 500.))
        p1km = pres(prof, to_msl(prof, 1000.))
        p3km = pres(prof, to_msl(prof, 3000.))
        p6km = pres(prof, to_msl(prof, 6000.))
        p9km = pres(prof, to_msl(prof, 9000.))
        p10km = pres(prof, to_msl(prof, 10000.))

        shear_0_to_500  = wind_shear(prof, pbot=sfc, ptop=p500m)
        shear_0_to_1000 = wind_shear(prof, pbot=sfc, ptop=p1km)
        shear_1_to_3000 = wind_shear(prof, pbot=p1km, ptop=p3km)
        shear_3_to_6000 = wind_shear(prof, pbot=p3km, ptop=p6km)
        shear_6_to_9000 = wind_shear(prof, pbot=p6km, ptop=p9km)
        shear_0_to_3000 = wind_shear(prof, pbot=sfc, ptop=p3km)
        shear_0_to_6000 = wind_shear(prof, pbot=sfc, ptop=p6km)
        # 0-10 km shear discriminates weak from derecho-producing MCSs (Cohen et al. 2007)
        shear_0_to_10000 = wind_shear(prof, pbot=sfc, ptop=p10km)
        shear_eil       = wind_shear(prof, pbot=kinem['eil'][0], ptop=kinem['eil'][1])

        kinem['shear_0_to_500']  = comp2vec(shear_0_to_500[0], shear_0_to_500[1])[1]
        kinem['shear_0_to_1000'] = comp2vec(shear_0_to_1000[0], shear_0_to_1000[1])[1]
        kinem['shear_1_to_3000'] = comp2vec(shear_1_to_3000[0], shear_1_to_3000[1])[1]
        kinem['shear_3_to_6000'] = comp2vec(shear_3_to_6000[0], shear_3_to_6000[1])[1]
        kinem['shear_6_to_9000'] = comp2vec(shear_6_to_9000[0], shear_6_to_9000[1])[1]
        kinem['shear_0_to_3000'] = comp2vec(shear_0_to_3000[0], shear_0_to_3000[1])[1]
        kinem['shear_0_to_6000'] = comp2vec(shear_0_to_6000[0], shear_0_to_6000[1])[1]
        kinem['shear_0_to_10000'] = comp2vec(shear_0_to_10000[0], shear_0_to_10000[1])[1]
        kinem['shear_eil'] = comp2vec(shear_eil[0], shear_eil[1])[1]

        # effective bulk wind difference: EIL base to half the MU EL height
        # Thompson, R. L., C. M. Mead, and R. Edwards, 2007
        try:
            ebwd_top = pres(prof, to_msl(prof, kinem['eil_z'][0] + (thermo['mu_el_z'] - kinem['eil_z'][0]) / 2))
            ebwd = wind_shear(prof, pbot=kinem['eil'][0], ptop=ebwd_top)
            kinem['ebwd'] = comp2vec(ebwd[0], ebwd[1])[1]
        except:
            kinem['ebwd'] = ma.masked





        #--- SRH ---#
        # ---------------------------------------------------------------
        if ma.is_masked(kinem['sm_u']) == False:
            
            keys = ['srh_0_to_500', 'srh_0_to_1000', 'srh_1_to_3000', 'srh_3_to_6000', 'srh_6_to_9000',
                    'srh_0_to_3000', 'srh_0_to_6000', 'srh_eil']
            lwrs = [0, 0, 1000, 3000, 6000, 0, 0, kinem['eil_z'][0]]
            uprs = [500, 1000, 3000, 6000, 9000, 3000, 6000, kinem['eil_z'][1]]
            
            for key, lwr, upr in zip(keys, lwrs, uprs):
                try:
                    kinem[key]  = helicity(prof, lwr, upr, stu = kinem['sm_u'], stv = kinem['sm_v'])[0]
                except:
                    kinem[key] = ma.masked
                    pass
        else:
            kinem['srh_0_to_500']  = ma.masked
            kinem['srh_0_to_1000'] = ma.masked
            kinem['srh_1_to_3000'] = ma.masked
            kinem['srh_3_to_6000'] = ma.masked
            kinem['srh_6_to_9000'] = ma.masked
            kinem['srh_6_to_9000'] = ma.masked
            kinem['srh_0_to_3000'] = ma.masked
            kinem['srh_0_to_6000'] = ma.masked
            kinem['srh_eil'] = ma.masked
            warnings.warn("Storm Relative Helicity could not be computed for this sounding (no valid storm motion)", Warning)






        #--- SRW LAYERS ---#
        # --------------------------------------------------------------- 
        if ma.is_masked(kinem['sm_u']) == False:
            def calc_srw_layer(prof, level):
                u, v = components(prof, p=pres(prof, to_msl(prof, level)))
                sru = u - kinem['sm_u']
                srv = v - kinem['sm_v']
                return mpcalc.wind_speed(sru*units.kts, srv*units.kts)
            
            keys = ['srw_0_to_500', 'srw_0_to_1000', 'srw_1_to_3000', 'srw_3_to_6000', 'srw_6_to_9000', 'srw_eil']
            lwrs = [0, 0, 1000, 3000, 6000, kinem['eil_z'][0]]
            uprs = [500, 1000, 3000, 6000, 9000, kinem['eil_z'][1]]
            
            for key, lwr, upr in zip(keys, lwrs, uprs):
                try:
                    kinem[key] = (calc_srw_layer(prof, upr)  + calc_srw_layer(prof, lwr))/2
                except:
                    kinem[key] = ma.masked
                    pass
        else:
            kinem['srw_0_to_500']  = ma.masked
            kinem['srw_0_to_1000'] = ma.masked
            kinem['srw_1_to_3000'] = ma.masked
            kinem['srw_3_to_6000'] = ma.masked
            kinem['srw_6_to_9000'] = ma.masked
            kinem['srw_eil'] = ma.masked
            warnings.warn("Storm Relative Wind could not be computed for this sounding (no valid storm motion)", Warning)





        #--- SRV ---#
        # ---------------------------------------------------------------
        # adopted from Sam Brandt (2022)    
        if ma.is_masked(kinem['sm_u']) == False: 
            # CONVERT TO m/s (uses `sm_u, sm_v` calculated above)
            u_ms = (intrp['uINTRP']/1.94384)
            v_ms = (intrp['vINTRP']/1.94384)
            sm_u_ms = ((intrp['uINTRP'] - kinem['sm_u'])/1.94384)
            sm_v_ms = ((intrp['vINTRP'] - kinem['sm_v'])/1.94384)

            # INTEROPLATED SRW (send back to knots)
            srw = mpcalc.wind_speed(sm_u_ms*units('m/s'), sm_v_ms*units('m/s'))
            kinem['srw'] = (srw.m*1.94384)

            # SHEAR COMPONENTS FOR VORT CALC
            # calc example = change in u over change in z
            dudz = (u_ms[2::]-u_ms[0:-2]) / (intrp['zINTRP'][2::]-intrp['zINTRP'][0:-2])
            dvdz = (v_ms[2::]-v_ms[0:-2]) / (intrp['zINTRP'][2::]-intrp['zINTRP'][0:-2])
            dudz = np.insert(dudz,0,dudz[0])
            dudz = np.insert(dudz,-1,dudz[-1])
            dvdz = np.insert(dvdz,0,dvdz[0])
            dvdz = np.insert(dvdz,-1,dvdz[-1])
            # Shear magnitude, 
            shear=(np.sqrt(dudz**2+dvdz**2)+0.0000001)
            # Vorticity components
            uvort=-dvdz
            vvort=dudz
            # Total horizontal vorticity
            kinem['vort'] = np.sqrt(uvort**2 + vvort**2)
            # Streamwise vorticity: projection of horizontal vorticity onto the storm-relative wind,
            # ω·(V-c)/|V-c|. Positive = streamwise, negative = antistreamwise. Only positive
            # (streamwise) vorticity is kept; antistreamwise levels are set to 0 so they count as
            # "no streamwise vorticity" instead of cancelling streamwise levels in layer means.
            kinem['swv'] = np.clip((sm_u_ms*uvort+sm_v_ms*vvort)/(np.sqrt(sm_u_ms**2+sm_v_ms**2)+0.0000001), 0, None)
            # Streamwiseness fraction: swv/|ω| = cos(angle between ω and the SR wind), as a percent
            # (|ω| equals the shear magnitude). Because swv is clipped, this ranges 0-100%.
            kinem['swv_perc'] = (kinem['swv']/shear)*100

            # layer average streamwiseness and total streamwise vorticity
            kinem['swv_perc_0_to_500']  = np.mean(kinem['swv_perc'][0:5])
            kinem['swv_perc_0_to_1000'] = np.mean(kinem['swv_perc'][0:10])
            kinem['swv_perc_1_to_3000'] = np.mean(kinem['swv_perc'][10:30])
            kinem['swv_perc_3_to_6000'] = np.mean(kinem['swv_perc'][30:60])
            kinem['swv_perc_6_to_9000'] = np.mean(kinem['swv_perc'][60:90])
            kinem['swv_perc_eil']       = np.mean(kinem['swv_perc'][eil_idx[0]:eil_idx[1]])
            kinem['swv_0_to_500']       = np.mean(kinem['swv'][0:5])
            kinem['swv_0_to_1000']      = np.mean(kinem['swv'][0:10])
            kinem['swv_1_to_3000']      = np.mean(kinem['swv'][10:30])
            kinem['swv_3_to_6000']      = np.mean(kinem['swv'][30:60])
            kinem['swv_6_to_9000']      = np.mean(kinem['swv'][60:90])
            kinem['swv_eil']            = np.mean(kinem['swv'][eil_idx[0]:eil_idx[1]])
        else:
            kinem['srw']                = ma.masked
            kinem['swv']                = ma.masked
            kinem['swv_perc']           = ma.masked
            kinem['vort']               = ma.masked
            kinem['swv_perc_0_to_500']  = ma.masked
            kinem['swv_perc_0_to_1000'] = ma.masked
            kinem['swv_perc_1_to_3000'] = ma.masked
            kinem['swv_perc_3_to_6000'] = ma.masked
            kinem['swv_perc_6_to_9000'] = ma.masked
            kinem['swv_perc_eil']       = ma.masked
            kinem['swv_0_to_500']       = ma.masked
            kinem['swv_0_to_1000']      = ma.masked
            kinem['swv_1_to_3000']      = ma.masked
            kinem['swv_3_to_6000']      = ma.masked
            kinem['swv_6_to_9000']      = ma.masked
            kinem['swv_eil']            = ma.masked
            warnings.warn("Streamwise Vorticity could not be computed for this sounding (no valid storm motion)", Warning)




        #--- COMPOSITES ---#
        # ---------------------------------------------------------------
        if ma.is_masked(kinem['sm_u']) == False:
            # effective inflow layer significant tornado parameter
            # Thompson, R. L., B. T. Smith, J. S. Grams, A. R. Dean, and C. Broyles, 2012
            kinem['eil_stp'] = stp_cin(thermo['mlcape'], kinem['srh_eil'], KTS2MS(kinem['ebwd']), thermo['ml_lcl_z'], thermo['mlcin'])
            # effective inflow layer supercell composite
            # Thompson et al. 2004
            kinem['eil_scp'] = scp(thermo['mucape'], kinem['srh_eil'], KTS2MS(kinem['ebwd']))
            # most-unstable energy helicity index
            kinem['ehi_0_to_1000'] = ehi(prof, mupcl, 0, 1000, stu=kinem['sm_u'], stv=kinem['sm_v'])
            kinem['ehi_0_to_3000'] = ehi(prof, mupcl, 0, 3000, stu=kinem['sm_u'], stv=kinem['sm_v'])

        else:
            kinem['eil_stp'] = ma.masked
            kinem['eil_scp'] = ma.masked
            kinem['ehi_0_to_1000'] = ma.masked
            kinem['ehi_0_to_3000'] = ma.masked




        #--- SIMPLE ENTRAINING CAPE ---#
        # ---------------------------------------------------------------
        # ECAPE uses the LFC and EL of the SHARPpy parcels above (as heights MSL), and is masked when the parcel has no CAPE
        for pcl_name, cape_type in (('mu', 'most_unstable'), ('ml', 'mixed_layer'), ('sb', 'surface_based')):
            try:
                cape, lfc, el = (thermo[f'{pcl_name}cape'], thermo[f'{pcl_name}_lfc_z'], thermo[f'{pcl_name}_el_z'])
                if any(ma.is_masked(x) or not np.isfinite(float(x)) for x in (cape, lfc, el)) or float(cape) <= 0:
                    raise ValueError('no CAPE')
                thermo[f'{pcl_name}_ecape'] = calc_ecape(z, p, T, general['spec_humidity'], u, v, cape_type,
                                                       lfc=z[0] + float(lfc)*units.m, el=z[0] + float(el)*units.m,
                                                       storm_motion='user_defined', u_sm=kinem['sm_u']*units.kts, v_sm=kinem['sm_v']*units.kts).m
            except Exception:
                thermo[f'{pcl_name}_ecape'] = ma.masked
                warnings.warn(f"{pcl_name.upper()}-ECAPE could not be computed for this sounding (calculation error)", Warning)
        
        
        
        #--- RH & MIXRAT WITH HEIGHT ---#
        # ---------------------------------------------------------------
        general['rh_0_500']  = np.mean(intrp['rhINTRP'][0:5])
        general['rh_0_1000'] = np.mean(intrp['rhINTRP'][0:10])
        general['rh_1_3000'] = np.mean(intrp['rhINTRP'][10:30])
        general['rh_3_6000'] = np.mean(intrp['rhINTRP'][30:60])
        general['rh_6_9000'] = np.mean(intrp['rhINTRP'][60:90])
        

        general['w_0_500']  = np.mean(intrp['mrINTRP'][0:5]) 
        general['w_0_1000'] = np.mean(intrp['mrINTRP'][0:10]) 
        general['w_1_3000'] = np.mean(intrp['mrINTRP'][10:30]) 
        general['w_3_6000'] = np.mean(intrp['mrINTRP'][30:60]) 
        general['w_6_9000'] = np.mean(intrp['mrINTRP'][60:90]) 
        
        return {"general": general, "thermo": thermo, "kinem": kinem, "intrp": intrp} 

    
    
    ######################################################################################################
    ### PRINT FUNCTION -- PRINT A FEW VALUES ###
    ######################################################################################################  
    def print_vals(self):
        params = self.calc()
        thermo = params["thermo"]
        kinem  = params["kinem"]
        print(' ')
        print(f'    > THERMODYNAMICS --------------------------------------------- ')
        print(f"    --- SBCAPE: {np.round(thermo['sbcape'], 1)} | MUCAPE: {np.round(thermo['mucape'], 1)} | MLCAPE: {np.round(thermo['mlcape'], 1)} | MUECAPE: {np.round(thermo['mu_ecape'], 1)}")
        print(f"    --- MU 0-3: {np.round(thermo['mu3cape'], 1)} | MU 0-6: {np.round(thermo['mu6cape'], 1)} | SB 0-3: {np.round(thermo['sb3cape'], 1)} | SB 0-6: {np.round(thermo['sb6cape'], 1)}")
        print(' ')
        print(f'    > KINEMATICS ------------------------------------------------- ')
        print(f"    --- 0-500 SRW: {np.round(kinem['srw_0_to_500'], 1)} | 0-500 SWV: {np.round(kinem['swv_0_to_500'],3)} | 0-500 SHEAR: {np.round(kinem['shear_0_to_500'], 1)} | 0-500 SRH: {np.round(kinem['srh_0_to_500'], 1)}")
        print(f"    --- 1-3km SRW: {np.round(kinem['srw_1_to_3000'], 1)} | 1-3km SWV: {np.round(kinem['swv_1_to_3000'],3)} | 1-3km SHEAR: {np.round(kinem['shear_1_to_3000'], 1)} | | 1-3km SRH: {np.round(kinem['srh_1_to_3000'], 1)}")
        print(f'    ============================================================== ')
        print(" ")



















        



# VAD Profile Calc functions 

class vad_params:
    
    def __init__(self, vad_data, storm_motion='right_moving', modify_sfc=None):
            self.vad_data = vad_data 
            self.storm_motion = storm_motion
            self.sm_type = storm_motion_type(storm_motion)
            self.modify_sfc = modify_sfc
            
    
    ######################################################################################################
    ### CALC FUNCTION -- COMPUTE VALUES ###
    ######################################################################################################        
    def calc(self):

    
    
        kinem  = {} 
        intrp  = {}

        z   = self.vad_data['z']
        u   = self.vad_data['u']
        v   = self.vad_data['v']
        wd  = mpcalc.wind_direction(u*units.kts, v*units.kts)
        ws  = mpcalc.wind_speed(u*units.kts, v*units.kts)

        # compute useful variables and add them to a new dict of 
        # data that this function will return 
        ##########################################   



        def interpolate(var, hgts, step):
            var = np.asarray(var, dtype=float)
            hgts = np.asarray(hgts, dtype=float)

            # Remove invalid values
            valid = np.isfinite(var) & np.isfinite(hgts)
            var = var[valid]
            hgts = hgts[valid]

            # Sort by height
            order = np.argsort(hgts)
            hgts = hgts[order]
            var = var[order]

            # Highest complete interpolation level available
            max_level = np.floor(np.max(hgts) / step) * step

            # Include the upper endpoint
            levels = np.arange(0, max_level + step, step)

            # Linear interpolation
            varinterp = np.interp(levels, hgts, var)

            return varinterp

        resolution=100
        
        uINTRP  = u[np.isnan(z)==False]
        vINTRP  = v[np.isnan(z)==False]
        zINTRP  = z[np.isnan(z)==False]
        zINTRP  = zINTRP-zINTRP[0]
        
        # Interpolation, add to intrp dict
        intrp['uINTRP'] = interpolate(uINTRP,zINTRP,resolution)
        intrp['vINTRP'] = interpolate(vINTRP,zINTRP,resolution)
        intrp['zINTRP'] = interpolate(zINTRP,zINTRP,resolution)

        hgt_lvls = {}
        hgt_var_list = ['h0', 'h05', 'h1', 'h15', 'h2', 'h25', 
                        'h3', 'h35', 'h4', 'h45', 'h5', 'h55', 
                        'h6', 'h65', 'h7', 'h75', 'h8', 'h85', 
                        'h9']
        
        hgt_val_list = [0.,500.,1000.,1500.,2000.,2500.,3000.,3500.,
                        4000.,4500.,5000.,5500.,6000.,6500.,7000.,
                        7500.,8000.,8500.,9000.]

        for var, val, in zip(hgt_var_list, hgt_val_list):
            try:
                hgt_lvls[var] = np.where(intrp['zINTRP']==val)[0][0]
            except:
                pass
        intrp['hgt_lvls'] = hgt_lvls

        # parameters for layers above the top of the profile are masked
        if intrp['zINTRP'][-1] < 9000:
            warnings.warn(f"This VAD profile only reaches {intrp['zINTRP'][-1]:.0f} m AGL, so parameters for layers above it could not be computed (not enough data)", Warning)

        # whether the profile reaches a given height
        def has_hgt(hgt):
            return hgt in intrp['zINTRP']
        

        ####################################################################
        ### KINEMATIC VALUES ###
        ####################################################################
        
        # THE FOLLOWING FUNCTIONS ARE CONTRUBUTED BY SCOTT THOMAS, NWS GRAND RAPIDS 
        # -------------------------------------------------------------------------  

        # CALC WIND VECTORS
        def calc_vector(u_comp, v_comp):
              mag = np.sqrt(u_comp**2 + v_comp**2)
              dir = np.rad2deg(np.arctan2(u_comp, v_comp)) % 360
              return mag, dir

            
        # CALC LAYER SHEAR
        def calc_shear(u_layer, v_layer, bottom, top, z):
              if not has_hgt(top):
                  return ma.masked
              layer_top = np.where(z == (top))[0][0]
              layer_bot = np.where(z == (bottom))[0][0]
              u_shr = u_layer[layer_top] - u_layer[layer_bot]
              v_shr = v_layer[layer_top] - v_layer[layer_bot]
              shrmag = np.hypot(u_shr, v_shr)
              return shrmag

        
        # CALC MEAN_WIND COMPONENTS 
        def calc_meanwind(u_layer, v_layer, z, top):
              if not has_hgt(top):
                  return ma.masked, ma.masked
              layer_top = np.where(z == (top))[0][0]
              mean_u = np.nanmean(u_layer[:layer_top])
              mean_v = np.nanmean(v_layer[:layer_top])
              return mean_u, mean_v

        
        # CALC BUNKERS STORM MOTION
        def calc_bunkers(u_layer, v_layer, z):
              if not has_hgt(6000):
                  return ma.masked, ma.masked, ma.masked, ma.masked
              layer_top = np.where(z == (6000))[0][0]
              mean_u = np.nanmean(u_layer[:layer_top])
              mean_v = np.nanmean(v_layer[:layer_top])

              layer_top = np.where(z == (6000))[0][0]
              u_shr = u_layer[layer_top] - u_layer[0]
              v_shr = v_layer[layer_top] - v_layer[0]

              dev = 7.5 * 1.94

              dev_amnt = dev / np.hypot(u_shr, v_shr)
              rmu = mean_u + (dev_amnt * v_shr)
              rmv = mean_v - (dev_amnt * u_shr)
              lmu = mean_u - (dev_amnt * v_shr)
              lmv = mean_v + (dev_amnt * u_shr)
              return rmu, rmv, lmu, lmv
            
            
    
        # CALC CORFIDI UP SHEAR & DOWN SHEAR STORM MOTION
        def calc_corfidi(u_layer, v_layer, z, u_mean, v_mean):
              if ma.is_masked(u_mean) or not has_hgt(1500):
                  return ma.masked
              llj_top = np.where(z == (1500))[0][0]
              llj_u = u_layer[:llj_top]
              llj_v = v_layer[:llj_top]

              mag, dir = calc_vector(llj_u, llj_v)
              max=0
              i=0
              for a in mag:
                if mag[i] >= mag[i-1]:
                  max = i

              u_max = llj_u[i]
              v_max = llj_v[i]

              corfidi_up_u = u_mean - u_max
              corfidi_up_v =  v_mean - v_max

              corfidi_down_u = u_mean + corfidi_up_u
              corfidi_down_v = v_mean + corfidi_up_v

              return [corfidi_up_u, corfidi_up_v, corfidi_down_u, corfidi_down_v]

        
        # CALC STORM RELATIVE HELICITY
        def calc_srh(u_comp, v_comp, z, depth, bottom, storm_u, storm_v):
            if not has_hgt((0 if bottom is None else bottom.m) + depth.m):
                return ma.masked
            srh = mpcalc.storm_relative_helicity(height=z, u=u_comp, v=v_comp, depth=depth,
                                           bottom=bottom, storm_u=storm_u, storm_v=storm_v)[0]
            return srh.m
                
            
        # CALC STORM RELATIVE WIND 
        def calc_srw_layer(u, v, bottom, top, z):
            if not has_hgt(top):
                return ma.masked
            layer_top = np.where(z == (top))[0][0]
            layer_bot = np.where(z == (bottom))[0][0]
            
            top_sru = u[layer_top] - kinem['sm_u']
            top_srv = v[layer_top] - kinem['sm_v']
            
            bot_sru = u[layer_bot] - kinem['sm_u']
            bot_srv = v[layer_bot] - kinem['sm_v']
            
            mean_srw = (mpcalc.wind_speed(top_sru*units.kts, top_srv*units.kts) + mpcalc.wind_speed(bot_sru*units.kts, bot_srv*units.kts))/2
            
            return mean_srw.m
        
                
        
            
        #--- BUNKERS STORM MOTION ---#
        # --------------------------------------------------------------- 
        if has_hgt(6000):
            kinem['sm_rm'] = calc_bunkers(intrp['uINTRP'], intrp['vINTRP'], intrp['zINTRP'])[0:2]
            kinem['sm_lm'] = calc_bunkers(intrp['uINTRP'], intrp['vINTRP'], intrp['zINTRP'])[2:4]
            kinem['sm_mw'] = calc_meanwind(intrp['uINTRP'], intrp['vINTRP'], intrp['zINTRP'], 6000)
        else:
            kinem['sm_rm'], kinem['sm_lm'], kinem['sm_mw'] = ma.masked, ma.masked, ma.masked
            warnings.warn("Bunkers Storm Motion could not be computed for this sounding (not enough data)", Warning)
       
        
        #--- USER DEFINED SM PT ---#
        # ---------------------------------------------------------------  
        if self.sm_type == 'custom':
            kinem['sm_u'], kinem['sm_v'] = mpcalc.wind_components(self.storm_motion[1]*units.kts, self.storm_motion[0]*units.deg)
            kinem['sm_u'], kinem['sm_v'] = kinem['sm_u'].m, kinem['sm_v'].m
            
        else:
            
            if ma.is_masked(kinem['sm_rm']) == False:
                kinem['sm_u'] = kinem[f'sm_{self.sm_type}'][0]
                kinem['sm_v'] = kinem[f'sm_{self.sm_type}'][1]

            else:
                kinem['sm_u'] = ma.masked
                kinem['sm_v'] = ma.masked
                
                
                
        #--- DEVIANT TORNADO MOTION ---#
        # --------------------------------------------------------------- 
        if ma.is_masked(kinem['sm_u']) == False:
            mw_u_300, mw_v_300 = calc_meanwind(intrp['uINTRP'],intrp['vINTRP'],intrp['zINTRP'], 300)
            kinem['dtm'] = ((kinem['sm_u']+mw_u_300)/2, (kinem['sm_v']+mw_v_300)/2)
        else:
            kinem['dtm'] = ma.masked
          
        
        
        #--- CORFIDI MCS MOTION ---#
        #----------------------------------------------------------------
        kinem['mcs'] = calc_corfidi(intrp['uINTRP'],intrp['vINTRP'],intrp['zINTRP'],
                                    calc_meanwind(intrp['uINTRP'],intrp['vINTRP'],intrp['zINTRP'], 6000)[0],
                                    calc_meanwind(intrp['uINTRP'],intrp['vINTRP'],intrp['zINTRP'], 6000)[1])


        
        #--- BULK SHEAR ---#
        # ---------------------------------------------------------------  
        kinem['shear_0_to_500']  = calc_shear(intrp['uINTRP'],intrp['vINTRP'], 0, 500, intrp['zINTRP'])
        kinem['shear_0_to_1000'] = calc_shear(intrp['uINTRP'],intrp['vINTRP'], 0, 1000, intrp['zINTRP'])
        kinem['shear_1_to_3000'] = calc_shear(intrp['uINTRP'],intrp['vINTRP'], 1000, 3000, intrp['zINTRP'])
        kinem['shear_3_to_6000'] = calc_shear(intrp['uINTRP'],intrp['vINTRP'], 3000, 6000, intrp['zINTRP'])
        kinem['shear_6_to_9000'] = calc_shear(intrp['uINTRP'],intrp['vINTRP'], 6000, 9000, intrp['zINTRP'])
            
            
            
        #--- SRH ---#
        # ---------------------------------------------------------------
        if ma.is_masked(kinem['sm_u']) == False:
            kinem['srh_0_to_500']  = calc_srh(u_comp=intrp['uINTRP']*units.kts, v_comp=intrp['vINTRP']*units.kts, z=intrp['zINTRP']*units.m, 
                                              bottom=None, depth=500*units.m,  
                                              storm_u=kinem['sm_u']*units.kts, storm_v=kinem['sm_v']*units.kts)
            kinem['srh_0_to_1000'] = calc_srh(u_comp=intrp['uINTRP']*units.kts, v_comp=intrp['vINTRP']*units.kts, z=intrp['zINTRP']*units.m, 
                                              bottom=1000*units.m, depth=2000*units.m, 
                                              storm_u=kinem['sm_u']*units.kts, storm_v=kinem['sm_v']*units.kts)
            kinem['srh_1_to_3000'] = calc_srh(u_comp=intrp['uINTRP']*units.kts, v_comp=intrp['vINTRP']*units.kts, z=intrp['zINTRP']*units.m, 
                                              bottom=1000*units.m, depth=3000*units.m, 
                                              storm_u=kinem['sm_u']*units.kts, storm_v=kinem['sm_v']*units.kts)
            kinem['srh_3_to_6000'] = calc_srh(u_comp=intrp['uINTRP']*units.kts, v_comp=intrp['vINTRP']*units.kts, z=intrp['zINTRP']*units.m, 
                                              bottom=3000*units.m, depth=3000*units.m, 
                                              storm_u=kinem['sm_u']*units.kts, storm_v=kinem['sm_v']*units.kts)
            kinem['srh_6_to_9000'] = calc_srh(u_comp=intrp['uINTRP']*units.kts, v_comp=intrp['vINTRP']*units.kts, z=intrp['zINTRP']*units.m, 
                                              bottom=6000*units.m, depth=3000*units.m, 
                                              storm_u=kinem['sm_u']*units.kts, storm_v=kinem['sm_v']*units.kts)
        else:
            kinem['srh_0_to_500']  = ma.masked
            kinem['srh_0_to_1000'] = ma.masked
            kinem['srh_1_to_3000'] = ma.masked
            kinem['srh_3_to_6000'] = ma.masked
            kinem['srh_6_to_9000'] = ma.masked
            warnings.warn("Storm Relative Helicity could not be computed for this sounding (no valid storm motion)", Warning)


            
        #--- SRW LAYERS ---#
        # --------------------------------------------------------------- 
        if ma.is_masked(kinem['sm_u']) == False:
            kinem['srw_0_to_500']  = calc_srw_layer(intrp['uINTRP'],intrp['vINTRP'], 0, 500, intrp['zINTRP'])
            kinem['srw_0_to_1000'] = calc_srw_layer(intrp['uINTRP'],intrp['vINTRP'], 0, 1000, intrp['zINTRP'])
            kinem['srw_1_to_3000'] = calc_srw_layer(intrp['uINTRP'],intrp['vINTRP'], 1000, 3000, intrp['zINTRP'])
            kinem['srw_3_to_6000'] = calc_srw_layer(intrp['uINTRP'],intrp['vINTRP'], 3000, 6000, intrp['zINTRP'])
            kinem['srw_6_to_9000'] = calc_srw_layer(intrp['uINTRP'],intrp['vINTRP'], 6000, 9000, intrp['zINTRP'])
        else:
            for key in ['srw_0_to_500', 'srw_0_to_1000', 'srw_1_to_3000', 'srw_3_to_6000', 'srw_6_to_9000']:
                kinem[key] = ma.masked

        
        
        #--- SRV ---#
        # ---------------------------------------------------------------
        # adopted from Sam Brandt (2022)    
        if ma.is_masked(kinem['sm_u']) == False: 
            # CONVERT TO m/s (uses `sm_u, sm_v` calculated above)
            u_ms = (intrp['uINTRP']/1.94384)
            v_ms = (intrp['vINTRP']/1.94384)
            sm_u_ms = ((intrp['uINTRP'] - kinem['sm_u'])/1.94384)
            sm_v_ms = ((intrp['vINTRP'] - kinem['sm_v'])/1.94384)

            # INTEROPLATED SRW (send back to knots)
            srw = mpcalc.wind_speed(sm_u_ms*units('m/s'), sm_v_ms*units('m/s'))
            kinem['srw'] = (srw.m*1.94384)

            # SHEAR COMPONENTS FOR VORT CALC
            # calc example = change in u over change in z
            dudz = (u_ms[2::]-u_ms[0:-2]) / (intrp['zINTRP'][2::]-intrp['zINTRP'][0:-2])
            dvdz = (v_ms[2::]-v_ms[0:-2]) / (intrp['zINTRP'][2::]-intrp['zINTRP'][0:-2])
            dudz = np.insert(dudz,0,dudz[0])
            dudz = np.insert(dudz,-1,dudz[-1])
            dvdz = np.insert(dvdz,0,dvdz[0])
            dvdz = np.insert(dvdz,-1,dvdz[-1])
            # Shear magnitude, 
            shear=(np.sqrt(dudz**2+dvdz**2)+0.0000001)
            # Vorticity components
            uvort=-dvdz
            vvort=dudz
            # Total horizontal vorticity
            kinem['vort'] = np.sqrt(uvort**2 + vvort**2)
            # Streamwise vorticity: projection of horizontal vorticity onto the storm-relative wind,
            # ω·(V-c)/|V-c|. Positive = streamwise, negative = antistreamwise. Only positive
            # (streamwise) vorticity is kept; antistreamwise levels are set to 0 so they count as
            # "no streamwise vorticity" instead of cancelling streamwise levels in layer means.
            kinem['swv'] = np.clip((sm_u_ms*uvort+sm_v_ms*vvort)/(np.sqrt(sm_u_ms**2+sm_v_ms**2)+0.0000001), 0, None)
            # Streamwiseness fraction: swv/|ω| = cos(angle between ω and the SR wind), as a percent
            # (|ω| equals the shear magnitude). Because swv is clipped, this ranges 0-100%.
            kinem['swv_perc'] = (kinem['swv']/shear)*100

            # layer average streamwiseness and total streamwise vorticity
            # layer means at 100 m resolution; masked when the profile doesn't reach the layer top
            def layer_mean(values, bottom, top):
                return np.mean(values[bottom // 100:top // 100]) if has_hgt(top) else ma.masked

            kinem['swv_perc_0_to_500']  = layer_mean(kinem['swv_perc'], 0, 500)
            kinem['swv_perc_0_to_1000'] = layer_mean(kinem['swv_perc'], 0, 1000)
            kinem['swv_perc_1_to_3000'] = layer_mean(kinem['swv_perc'], 1000, 3000)
            kinem['swv_perc_3_to_6000'] = layer_mean(kinem['swv_perc'], 3000, 6000)
            kinem['swv_perc_6_to_9000'] = layer_mean(kinem['swv_perc'], 6000, 9000)
            kinem['swv_0_to_500']       = layer_mean(kinem['swv'], 0, 500)
            kinem['swv_0_to_1000']      = layer_mean(kinem['swv'], 0, 1000)
            kinem['swv_1_to_3000']      = layer_mean(kinem['swv'], 1000, 3000)
            kinem['swv_3_to_6000']      = layer_mean(kinem['swv'], 3000, 6000)
            kinem['swv_6_to_9000']      = layer_mean(kinem['swv'], 6000, 9000)
        else:
            kinem['srw']                = ma.masked
            kinem['swv']                = ma.masked
            kinem['swv_perc']           = ma.masked
            kinem['vort']               = ma.masked
            kinem['swv_perc_0_to_500']  = ma.masked
            kinem['swv_perc_0_to_1000'] = ma.masked
            kinem['swv_perc_1_to_3000'] = ma.masked
            kinem['swv_perc_3_to_6000'] = ma.masked
            kinem['swv_perc_6_to_9000'] = ma.masked
            kinem['swv_0_to_500']       = ma.masked
            kinem['swv_0_to_1000']      = ma.masked
            kinem['swv_1_to_3000']      = ma.masked
            kinem['swv_3_to_6000']      = ma.masked
            kinem['swv_6_to_9000']      = ma.masked
            warnings.warn("Streamwise Vorticity could not be computed for this sounding (no valid storm motion)", Warning)
        
        return {'kinem': kinem, 'intrp': intrp} 