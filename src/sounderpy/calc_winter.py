"""
    SOUNDERPY WINTER CALC

    Purpose of module:

    Compute winter-weather sounding parameters: precipitation type, snow growth, and winter
    hazards. `winter_params(clean_data).calc()` runs the standard `sounding_params` calc and adds
    a `winter` dict to what it returns.

    Sources are noted beside each parameter: [peer-reviewed] journal articles, [operational] SPC/NWS
    practice that isn't formally published, [informal] methods without a paper, and [prototype]
    choices made for the winter prototype that still need review.

    (C) KYLE J GILLETT, 2026
"""

# STANDARD PYTHON
import warnings

# NUMPY OPERATIONS
import numpy as np
import numpy.ma as ma

# METPY OPERATIONS
import metpy.calc as mpcalc

# SOUNDERPY MODULES
from .calc import sounding_params
from .utils import modify_surface
from .SHARPPYMAIN.sharppy.sharptab.profile import create_profile
from .SHARPPYMAIN.sharppy.sharptab import watch_type, params as sharp_params, interp as sharp_interp




# CONSTANTS
G = 9.81                # gravity, m/s2
T0 = 273.15             # 0 °C in K
KT_TO_MS = 0.514444
PA_S_TO_UB_S = 10.0     # 1 Pa/s = 10 µb/s, SHARPpy's omega unit




#########################################################################
############################## HELPERS ##################################
#########################################################################

# saturation vapor pressure (hPa) over water and over ice
# Alduchov & Eskridge (1996, J. Appl. Meteor.) [peer-reviewed]
def es_water(T):
    return 6.1094 * np.exp(17.625 * T / (T + 243.04))


def es_ice(T):
    return 6.1121 * np.exp(22.587 * T / (T + 273.86))


# mean / max of values between two heights; masked when the profile doesn't reach the layer top
def layer_mean(values, z, bottom, top):
    sel = (z >= bottom) & (z <= top) & np.isfinite(values)
    return float(np.mean(values[sel])) if sel.sum() > 1 and z.max() >= top else ma.masked


def layer_max(values, z, bottom, top):
    sel = (z >= bottom) & (z <= top) & np.isfinite(values)
    return float(np.max(values[sel])) if sel.any() else ma.masked


# value at one height; masked outside the profile
def value_at(values, z, height):
    return float(np.interp(height, z, values)) if z.min() <= height <= z.max() else ma.masked


# every height (lowest first) where a profile crosses 0 °C, interpolated between levels
def zero_crossings(z, T):
    i = np.where(np.isfinite(T[:-1]) & np.isfinite(T[1:]) & (T[:-1] * T[1:] < 0))[0]
    return list(z[i] + (z[i + 1] - z[i]) * T[i] / (T[i] - T[i + 1]))


# (start, end) indices of each run of True in a boolean array, end inclusive
def true_runs(mask):
    edges = np.diff(np.concatenate(([0], mask.astype(int), [0])))
    return list(zip(np.where(edges == 1)[0], np.where(edges == -1)[0] - 1))


# the dendritic growth zone (-12 to -17 °C) that makes the snow: of every layer in that range (deeper than
# 50 m), the one with the most cloud in it; with no cloud in any, the lowest (SHARPpy's dgz() choice). On a
# 10 m grid. Returns (bottom, top) heights and pressures, or None. A surface colder than -12 °C, or a
# warm layer aloft, can give several such layers, and SHARPpy's dgz() (first -12 and first -17 °C levels)
# can then give an inverted or cloud-free one
def dendritic_growth_zone(z, T, p, cloud_z, cloud_rh_ice, z_top=12000):
    z_f = np.arange(0, min(np.nanmax(z), z_top), 10.0)
    T_f = np.interp(z_f, z, T)
    in_cloud = np.interp(z_f, cloud_z, cloud_rh_ice) >= 90
    layers = [(a, b) for a, b in true_runs((T_f <= -12) & (T_f >= -17)) if z_f[b] - z_f[a] > 50]
    if not layers:
        return None
    a, b = max(layers, key=lambda ab: (in_cloud[ab[0]:ab[1] + 1].sum(), -ab[0]))
    p_f = np.exp(np.interp(z_f[[a, b]], z, np.log(p)))
    return [float(z_f[a]), float(z_f[b])], [float(p_f[0]), float(p_f[1])]


# layers (m AGL, top to bottom) where a mask holds for at least `min_depth`, with their pressures
def find_layers(mask, z, p, min_depth=200):
    return [{"bottom": float(z[a]), "top": float(z[b]), "bottom_p": float(p[a]), "top_p": float(p[b])}
            for a, b in true_runs(mask) if z[b] - z[a] >= min_depth]




#########################################################################
########################## PRECIPITATION TYPE ###########################
#########################################################################

# split the profile below `top` into melting (> 0 °C) and refreezing (< 0 °C) layers, top first, each with its
# energy: g * integral((T - T0) / T0) dz (J/kg), the area between the profile and 0 °C. Bourgouin (2000, Wea.
# Forecasting) [peer-reviewed]. SHARPpy's posneg_temperature() lumps every layer below the first warm layer
# together, so the layers are found here one by one
def melting_refreezing_layers(z, T, top):
    sel = np.isfinite(T) & (z <= top)
    z_full, T_full = z, T
    z, T = z[sel], T[sel]
    if len(z) < 2:
        return []
    # end exactly at `top`
    if z[-1] < top <= z_full.max():
        z, T = np.append(z, top), np.append(T, np.interp(top, z_full, T_full))

    # insert the exact 0 °C crossings so every layer starts and ends at 0 °C
    crossings = zero_crossings(z, T)
    z_all = np.concatenate((z, crossings))
    T_all = np.concatenate((T, np.zeros(len(crossings))))
    order = np.argsort(z_all)
    z, T = z_all[order], T_all[order]

    layers = []
    for i in range(len(z) - 1):
        T_mean = (T[i] + T[i + 1]) / 2
        if T_mean == 0:
            continue
        kind = "warm" if T_mean > 0 else "cold"
        energy = G * T_mean / T0 * (z[i + 1] - z[i])
        extreme = max if kind == "warm" else min
        if layers and layers[-1]["kind"] == kind:
            layers[-1].update(top=z[i + 1], energy=layers[-1]["energy"] + energy,
                              extreme=extreme(layers[-1]["extreme"], T[i + 1]))
        else:
            layers.append({"kind": kind, "bottom": z[i], "top": z[i + 1], "energy": energy, "extreme": extreme(T[i], T[i + 1])})

    # drop slivers (< 50 m with < 1 J/kg)
    layers = [lyr for lyr in layers if lyr["top"] - lyr["bottom"] >= 50 or abs(lyr["energy"]) >= 1]
    for lyr in layers:
        lyr["depth"] = lyr["top"] - lyr["bottom"]
        lyr["surface"] = lyr["bottom"] <= 1.0
    return layers[::-1]


# follow precipitation down from its source, layer by layer, recording each phase change.
# start: from the top temperature of the cloud the source sits in (the source's own when no cloud matches):
#   < -10 °C snow, 0 to -10 °C supercooled drizzle (little ice), > 0 °C rain.
#   Huffman & Norman (1988, MWR); Rauber et al. (2000, J. Appl. Meteor.) [peer-reviewed]
# energies are from the wet-bulb profile, as in the revised Bourgouin method of Birk et al. (2021, Wea. Forecasting,
# 36, 425-438) [peer-reviewed]. Melting layer (energy ME), ice stays frozen / partly melts / melts:
#   surface-based: below 5.6 / 5.6-13.2 / above 13.2 J/kg, Bourgouin (2000) [peer-reviewed]; Birk et al. have no
#     rain-vs-snow cases, so Bourgouin's are the only thresholds set for this
#   elevated: below 9.4 / 9.4-15 / above 15 J/kg, where Birk's ProbSN = 1540 exp(-0.29 ME) is 100 % and 20 %
#   (both start above 2 J/kg, the smallest melting layer in Bourgouin 2000 and Birk et al. 2021)
# refreezing layer under it (energy RE), Birk et al. (2021): ice pellets if RE > 208 + 0.08 ME (freezing rain
#   < 20 %), a mix if RE >= 117 + 0.08 ME, else it stays liquid; partly melted particles refreeze to ice pellets
# liquid reaching a surface at or below 0 °C is freezing rain / drizzle
# [prototype] applied layer by layer (Birk et al. use column totals)
# returns (trace [{"z", "phase"}], surface ptype, why it started as it did)
MELT_SFC, MELT_ALOFT = (5.6, 13.2), (9.4, 15.0)
REFREEZE_MIX, REFREEZE_PL, REFREEZE_SLOPE = 117.0, 208.0, 0.08


def precip_phase_trace(layers, source_T, source_z, cloud_top_T, sfc_T):
    if ma.is_masked(source_T) or ma.is_masked(source_z):
        return [], "None", "no saturated source layer"

    start_T, where = (cloud_top_T, "cloud top") if cloud_top_T is not None else (source_T, "source")
    if start_T > 0:
        phase, why = "RA", f"{start_T:+.0f}°C {where} (liquid)"
    elif start_T > -10:
        phase, why = "DZ", f"{start_T:.0f}°C {where} (little ice)"
    else:
        phase, why = "SN", f"{start_T:.0f}°C {where} (ice)"
    trace = [{"z": float(source_z), "phase": phase}]
    from_ice = phase == "SN"
    melt_pa = 0.0

    for lyr in layers:
        before = phase
        if lyr["kind"] == "warm" and phase in ("SN", "PL"):
            melt_pa = lyr["energy"]
            partial, full = MELT_SFC if lyr["surface"] else MELT_ALOFT
            phase = phase if melt_pa < partial else ("RA/" + phase if melt_pa <= full else "RA")
        elif lyr["kind"] == "cold":
            na = -lyr["energy"]
            if phase.startswith("RA/"):
                phase = "PL"
            elif phase == "RA" and from_ice and melt_pa > 0:
                if na > REFREEZE_PL + REFREEZE_SLOPE * melt_pa:
                    phase = "PL"
                elif na >= REFREEZE_MIX + REFREEZE_SLOPE * melt_pa:
                    phase = "FZRA/PL"
        if phase != before:
            trace.append({"z": float(lyr["bottom"]), "phase": phase})

    if phase in ("RA", "DZ") and sfc_T <= 0:
        phase = "FZ" + phase
        trace.append({"z": 0.0, "phase": phase})
    return trace, phase, why




#########################################################################
############################ WINTER HAZARDS #############################
#########################################################################

# hazard level from a set of ingredient checks: HIGH all met, MOD one short, LOW otherwise.
# `required` must be met for anything above NONE; `needed_for_mod` must be met for MOD
def hazard_level(ok, required, needed_for_mod=None):
    met = sum(ok.values())
    if not ok[required]:
        return "NONE"
    if met == len(ok):
        return "HIGH"
    if met == len(ok) - 1 and (needed_for_mod is None or ok[needed_for_mod]):
        return "MOD"
    return "LOW"


# freezing drizzle POTENTIAL, checked on the lowest cloud (base below 2 km); an ingredients check, not a forecast
#   sfc   surface Tw <= 0 °C: drizzle freezes on contact
#   top   cloud top > -10 °C: little ice, drizzle forms by the warm-rain process (Huffman & Norman 1988; Rauber et al. 2000)
#   dry   mean RHi <= 75 % in the 1 km above the top: no ice falls in from above (McDonough et al. 2004, NCAR FIP);
#         needed for MOD
#   deep  cloud >= 500 m deep: deep enough for drizzle drops to grow
#   lift  mean omega < 0 in the cloud: shown when omega is available, not counted in the level
# [prototype] the 2 km, 75 % and 500 m thresholds
def freezing_drizzle_potential(clouds, z, rh_ice, omega, has_omega, sfc_Tw):
    low = [cloud for cloud in clouds if cloud["bottom"] <= 2000]
    if not low:
        ok = {"sfc": bool(sfc_Tw <= 0), "top": False, "dry": False, "deep": False}
        return {"level": "NONE", "ok": ok, "cloud": None, "rh_above": np.nan, "depth": np.nan, "omega": np.nan, "lift": None}

    cloud = low[0]
    above = (z > cloud["top"]) & (z <= cloud["top"] + 1000)
    rh_above = float(np.mean(rh_ice[above])) if above.sum() > 1 else np.nan
    depth = cloud["top"] - cloud["bottom"]
    ok = {"sfc": bool(sfc_Tw <= 0), "top": bool(cloud["top_T"] > -10),
          "dry": bool(np.isfinite(rh_above) and rh_above <= 75), "deep": bool(depth >= 500)}

    inside = (z >= cloud["bottom"]) & (z <= cloud["top"])
    cloud_omega = float(np.nanmean(omega[inside])) if has_omega and inside.sum() > 1 else np.nan
    lift = bool(cloud_omega < 0) if has_omega else None

    # a cold cloud top makes ice, so it caps the level at LOW; so does moist air above the cloud, where ice from a
    # cloud above can fall in (e.g. KDTX 2023-02-23 00Z: snow from a -33 °C top melts and falls as freezing rain)
    level = hazard_level(ok, "sfc", "dry") if ok["top"] else ("LOW" if ok["sfc"] else "NONE")
    return {"level": level, "ok": ok, "cloud": cloud, "rh_above": rh_above, "depth": depth, "omega": cloud_omega, "lift": lift}


# freezing fog POTENTIAL from four surface ingredients; an ingredients check, not a forecast
#   cold  sfc T < 0 °C: droplets are supercooled and rime on contact
#   sat   sfc RH >= 95 %
#   calm  10 m wind <= 6 kt: mixing doesn't lift or dry out the fog
#   inv   a surface-based inversion traps cold, moist air at the ground
# Gultepe et al. (2007, Pure Appl. Geophys.) [peer-reviewed]. Called freezing fog at any temperature (no separate ice fog)
# [prototype] the 95 % and 6 kt thresholds. Checks use the same rounded values the figure shows
def freezing_fog_potential(sfc_T, sfc_rh, sfc_wind_kt, inv_strength):
    # a missing value (e.g. no surface wind before the radiosonde's GPS locks on) isn't met
    ok = {"cold": bool(sfc_T < 0), "sat": bool(np.isfinite(sfc_rh) and round(sfc_rh) >= 95),
          "calm": bool(np.isfinite(sfc_wind_kt) and round(sfc_wind_kt) <= 6), "inv": bool(inv_strength > 0)}
    return {"level": hazard_level(ok, "cold", "sat"), "ok": ok,
            "values": {"cold": sfc_T, "sat": sfc_rh, "calm": sfc_wind_kt, "inv": inv_strength}}




#########################################################################
############################## MAIN CLASS ###############################
#########################################################################

class winter_params:
    """Winter sounding parameters: everything ``sounding_params`` computes, plus a ``'winter'`` group.

    Takes the same arguments as ``sounding_params``, and ``.calc()`` returns the same dictionary
    (``general``, ``thermo``, ``kinem``, ``intrp``) with ``winter`` added: precipitation type (Bourgouin
    and SPC), melting and refreezing layers, clouds and icing, the DGZ and snow ratio, the mixed layer,
    and winter hazards and winds. Use as ``spy.winter_params(clean_data).calc()``.
    """


    def __init__(self, clean_data, storm_motion='right_moving', modify_sfc=None, include_all_parcels=False, parcel_depth=None):
        # a modified surface goes into both the standard and the winter values
        self.clean_data = modify_surface(clean_data, modify_sfc) if isinstance(modify_sfc, dict) else clean_data
        self.storm_motion = storm_motion
        self.include_all_parcels = include_all_parcels
        self.parcel_depth = parcel_depth


    ######################################################################################################
    ### CALC FUNCTION -- COMPUTE VALUES ###
    ######################################################################################################
    def calc(self):

        # the standard calc first: the reused plot helpers expect its output, and the winter values build on it
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            snd = sounding_params(self.clean_data, self.storm_motion, include_all_parcels=self.include_all_parcels,
                                  parcel_depth=self.parcel_depth).calc()
        general, thermo, intrp = snd['general'], snd['thermo'], snd['intrp']
        winter = {}

        # declare easy variable names for reuse from `clean_data` (data levels)...
        p = self.clean_data['p'].m
        z = self.clean_data['z'].m - self.clean_data['z'].m[0]
        T = self.clean_data['T'].to('degC').m
        Td = self.clean_data['Td'].to('degC').m
        u = self.clean_data['u'].to('kt').m
        v = self.clean_data['v'].to('kt').m
        Tw = np.asarray(getattr(general['wet_bulb'], 'm', general['wet_bulb']), dtype=float)
        has_omega = 'omega' in self.clean_data and np.any(np.asarray(getattr(self.clean_data['omega'], 'm', self.clean_data['omega'])) != 0)
        omega = np.asarray(self.clean_data['omega'].to('Pa/s').m, dtype=float) * PA_S_TO_UB_S if has_omega else None

        # ...and from `intrp` (100 m AGL grid)
        zi = np.asarray(intrp['zINTRP'], dtype=float)
        Ti = np.asarray(intrp['tINTRP'], dtype=float)
        Twi = np.asarray(intrp['wbINTRP'], dtype=float)
        pi = np.asarray(intrp['pINTRP'], dtype=float)
        ui, vi = np.asarray(intrp['uINTRP'], dtype=float), np.asarray(intrp['vINTRP'], dtype=float)
        rhi = np.asarray(intrp['rhINTRP'], dtype=float)
        rhi = rhi * 100 if np.nanmax(rhi) <= 1.5 else rhi
        thetaei = np.asarray(intrp['thetaeINTRP'], dtype=float)
        omegai = np.interp(zi, z, omega) if has_omega else np.full_like(zi, np.nan)

        # RH with respect to ice below 0 °C (water above). Alduchov & Eskridge (1996) [peer-reviewed]
        rh_icei = np.where(Ti < 0, rhi * es_water(Ti) / es_ice(Ti), rhi)

        # profiles for the side panels
        winter.update(z=zi, T=Ti, Tw=Twi, rh=rhi, rh_ice=rh_icei, thetae=thetaei, u=ui, v=vi, omega=omegai, has_omega=has_omega)

        # SHARPpy profile, with omega when there is some (SounderPy's own profile has none)
        wspd = np.hypot(u, v)
        wdir = mpcalc.wind_direction(self.clean_data['u'], self.clean_data['v']).m
        profile_kwargs = dict(profile='default', pres=p, hght=self.clean_data['z'].m, tmpc=T, dwpc=Td,
                              wspd=wspd, wdir=wdir, missing=-9999, strictQC=False)
        if has_omega:
            profile_kwargs['omeg'] = omega
        prof = create_profile(**profile_kwargs)



        ####################################################################
        ### PRECIPITATION TYPE ###
        ####################################################################
        #--- PRECIP SOURCE ---#
        # ---------------------------------------------------------------
        # SHARPpy init_phase(): highest saturated (RH > 80 %) 50 hPa layer below 5 km AGL, limited to rising
        # motion when omega is available. SPC SHARP code (R. Thompson) via SHARPpy [operational]
        src_p, src_phase, src_T, _ = watch_type.init_phase(prof)
        if ma.is_masked(src_p) or src_p is None or src_p <= 0 or src_p == prof.missing:
            winter['source_p'] = winter['source_z'] = winter['source_T'] = ma.masked
        else:
            winter['source_p'] = float(src_p)
            winter['source_z'] = float(sharp_interp.to_agl(prof, sharp_interp.hght(prof, src_p)))
            winter['source_T'] = float(src_T)

        #--- CLOUD & ICING LAYERS ---#
        # ---------------------------------------------------------------
        # [prototype] cloud: RH (over ice below 0 °C) >= 90 %. icing: 0 to -15 °C with RH over water >= 90 %,
        # where supercooled liquid ices aircraft (the ingredients of NCAR's icing products: Bernstein et al. 2005,
        # J. Appl. Meteor.; McDonough et al. 2004, FIP). The 0 to -15 °C range is our choice
        winter['clouds'] = find_layers((rh_icei >= 90) & (zi <= 8000), zi, pi)
        for cloud in winter['clouds']:
            cloud['top_T'] = float(np.interp(cloud['top'], zi, Ti))
        winter['icing_layers'] = find_layers((Ti <= 0) & (Ti >= -15) & (rhi >= 90) & (zi <= 8000), zi, pi)

        # the cloud the source sits in: its top temperature, not the source's, decides whether the cloud makes ice.
        # SHARPpy's source test (RH > 80 % over water) and the cloud test differ, so allow 300 m either way
        winter['source_cloud'] = None
        if not ma.is_masked(winter['source_z']):
            winter['source_cloud'] = next((cloud for cloud in winter['clouds']
                                           if cloud['bottom'] - 300 <= winter['source_z'] <= cloud['top'] + 300), None)

        #--- MELTING & REFREEZING LAYERS ---#
        # ---------------------------------------------------------------
        # energies from the wet-bulb profile: falling precip evaporates and cools unsaturated air toward Tw, so T
        # overstates melting under dry layers (in saturated layers Tw = T). Bourgouin (2000) used T; the wet-bulb
        # version is the revised Bourgouin method of Birk et al. (2021, Wea. Forecasting) [peer-reviewed]
        source_top = winter['source_z'] if not ma.is_masked(winter['source_z']) else 5000.0
        winter['layers'] = melting_refreezing_layers(z, Tw, source_top)
        winter['Tw_levels'] = Tw

        #--- PRECIP TYPE (BOURGOUIN) ---#
        # ---------------------------------------------------------------
        cloud_top_T = winter['source_cloud']['top_T'] if winter['source_cloud'] else None
        winter['phase_trace'], winter['ptype'], winter['ptype_why'] = precip_phase_trace(
            winter['layers'], winter['source_T'], winter['source_z'], cloud_top_T, T[0])
        for step in winter['phase_trace']:
            step['p'] = float(np.interp(step['z'], z, p))

        #--- PRECIP TYPE (SPC BEST GUESS) ---#
        # ---------------------------------------------------------------
        # SHARPpy best_guess_precip(), from SPC SHARP code [operational]. posneg_temperature() needs start=:
        # its default path unpacks init_phase() wrongly and raises ValueError
        try:
            start = winter['source_p'] if not ma.is_masked(winter['source_p']) else 500.0
            pos, neg, _, _ = watch_type.posneg_temperature(prof, start=start)
            winter['sharppy_ptype'] = str(watch_type.best_guess_precip(prof, src_phase, src_p, src_T, pos, neg)).rstrip('.')
        except Exception as ex:
            winter['sharppy_ptype'] = 'NA'
            warnings.warn(f"SHARPpy best guess precip failed: {ex}")

        #--- FREEZING LEVELS & WARM NOSE ---#
        # ---------------------------------------------------------------
        # every 0 °C crossing of T and of Tw; the warm nose is the warmest air temperature in an elevated warm layer
        winter['frz_pts_z'] = zero_crossings(z, T)
        winter['wb_frz_pts_z'] = zero_crossings(z, Tw)
        winter['warm_nose_T'] = winter['warm_nose_p'] = ma.masked
        warm_aloft = [lyr for lyr in melting_refreezing_layers(z, T, source_top) if lyr['kind'] == 'warm' and not lyr['surface']]
        if warm_aloft:
            nose = max(warm_aloft, key=lambda lyr: lyr['extreme'])
            i = np.argmax(np.where((z >= nose['bottom']) & (z <= nose['top']), T, -999))
            winter['warm_nose_T'], winter['warm_nose_p'] = float(T[i]), float(p[i])
        ####################################################################



        ####################################################################
        ### SNOW GROWTH ###
        ####################################################################
        #--- DENDRITIC GROWTH ZONE ---#
        # ---------------------------------------------------------------
        # -12 to -17 °C, the layer in cloud (see dendritic_growth_zone()). Lab basis: Takahashi et al. (1991,
        # J. Meteor. Soc. Japan); Libbrecht (2005, Rep. Prog. Phys.) [peer-reviewed]
        dgz = dendritic_growth_zone(z, T, p, zi, rh_icei)
        has_dgz = dgz is not None
        winter['dgz_z'], winter['dgz_p'] = dgz if has_dgz else ([ma.masked, ma.masked], [ma.masked, ma.masked])
        winter['dgz_depth'] = winter['dgz_z'][1] - winter['dgz_z'][0] if has_dgz else ma.masked
        winter['dgz_rh'] = layer_mean(rhi, zi, *winter['dgz_z']) if has_dgz else ma.masked
        winter['dgz_rh_ice'] = layer_mean(rh_icei, zi, *winter['dgz_z']) if has_dgz else ma.masked
        winter['dgz_omega'] = layer_mean(omegai, zi, *winter['dgz_z']) if has_dgz and has_omega else ma.masked

        #--- OPRH ---#
        # ---------------------------------------------------------------
        # mean omega (µb/s) x precipitable water within the DGZ (in) x mean RH (fraction), all in the DGZ;
        # negative favors snow. As SHARPpy's profile.oprh (SPC winter inset) [operational]
        winter['dgz_pw'] = (float(sharp_params.precip_water(prof, pbot=winter['dgz_p'][0], ptop=winter['dgz_p'][1]))
                            if has_dgz else ma.masked)
        if ma.is_masked(winter['dgz_omega']) or ma.is_masked(winter['dgz_rh']) or ma.is_masked(winter['dgz_pw']):
            winter['oprh'] = ma.masked
        else:
            winter['oprh'] = winter['dgz_omega'] * winter['dgz_pw'] * winter['dgz_rh'] / 100

        #--- SNOW RATIO ---#
        # ---------------------------------------------------------------
        # Kuchera, from the warmest temperature (K) between the surface and 500 hPa: 12 + 2 (271.16 - Tmax) when
        # Tmax > 271.16 K, else 12 + (271.16 - Tmax). E. Kuchera (AFWA, ~2004), unpublished; as documented in the
        # NOAA VLab Forecast Guide [informal]. Not capped; shown only when snow reaches the surface
        Tmax = float(np.max(T[p >= 500])) + T0
        kuchera = max(12 + 2 * (271.16 - Tmax) if Tmax > 271.16 else 12 + (271.16 - Tmax), 0.0)
        winter['slr'] = kuchera if 'SN' in winter['ptype'] else ma.masked
        ####################################################################



        ####################################################################
        ### WINTER HAZARDS ###
        ####################################################################
        #--- SURFACE ---#
        # ---------------------------------------------------------------
        winter['sfc_T_F'] = T[0] * 9 / 5 + 32
        winter['sfc_Tw'] = float(Tw[0])
        winter['sfc_wind_uv'] = (float(u[0]), float(v[0]))
        winter['sfc_wind_kt'] = float(wspd[0])

        #--- WIND CHILL ---#
        # ---------------------------------------------------------------
        # NWS (2001) formula, °F; Osczevski & Bluestein (2005, BAMS) [peer-reviewed]. Defined at <= 50 °F and >= 3 mph
        wind_mph = winter['sfc_wind_kt'] * 1.15078
        if winter['sfc_T_F'] <= 50 and wind_mph >= 3:
            winter['wind_chill'] = (35.74 + 0.6215 * winter['sfc_T_F'] - 35.75 * wind_mph ** 0.16
                                    + 0.4275 * winter['sfc_T_F'] * wind_mph ** 0.16)
        else:
            winter['wind_chill'] = ma.masked

        #--- PBL ---#
        # ---------------------------------------------------------------
        # SHARPpy's PBL depth and lapse rate, from the standard calc; its mean wind is under WINDS
        for key in ('pbl_depth', 'pbl_lapse_rate'):
            winter[key] = thermo[key]

        # gust potential: the strongest wind between the surface and the PBL top, what mixing could bring
        # down. BUFKIT momentum-transfer practice [operational]
        pbl_top_z = float(np.interp(-thermo['pbl_top'], -pi, zi)) if not ma.is_masked(thermo['pbl_top']) else 1000.0
        winter['gust_potential'] = layer_max(np.hypot(ui, vi), zi, 0, max(pbl_top_z, 100))
        winter['max_wind_0_to_2000'] = layer_max(np.hypot(ui, vi), zi, 0, 2000)

        #--- SURFACE INVERSION ---#
        # ---------------------------------------------------------------
        # a surface-based inversion: the warmest point in the lowest 1.5 km when it's > 0.5 °C above the
        # surface (model profiles often cool slightly in the lowest layer first) [prototype]
        i_max = int(np.nanargmax(np.where(zi <= 1500, Ti, -999)))
        inversion = Ti[i_max] - Ti[0] > 0.5
        winter['sfc_inv_depth'] = float(zi[i_max]) if inversion else 0.0
        winter['sfc_inv_strength'] = float(Ti[i_max] - Ti[0]) if inversion else 0.0

        #--- SNOW SQUALL PARAMETER ---#
        # ---------------------------------------------------------------
        # [(0-2 km RH - 60) / 15] x [(4 - 0-2 km theta-e change) / 4] x [0-2 km mean wind / 9 m/s], zero when
        # RH < 60 %, the theta-e change > 4 K, or sfc Tw > 1 °C. Banacos et al. (2014, J. Operational Meteor.) [peer-reviewed]
        rh_0_to_2000 = layer_mean(rhi, zi, 0, 2000)
        wind_0_to_2000 = layer_mean(np.hypot(ui, vi), zi, 0, 2000)
        dthetae = value_at(thetaei, zi, 2000) - thetaei[0] if zi.max() >= 2000 else ma.masked
        if ma.is_masked(rh_0_to_2000) or ma.is_masked(dthetae) or ma.is_masked(wind_0_to_2000):
            winter['snsq'] = ma.masked
        elif rh_0_to_2000 < 60 or dthetae > 4 or Twi[0] > 1:
            winter['snsq'] = 0.0
        else:
            winter['snsq'] = ((rh_0_to_2000 - 60) / 15) * ((4 - dthetae) / 4) * (wind_0_to_2000 * KT_TO_MS / 9)

        #--- FREEZING FOG & FREEZING DRIZZLE ---#
        # ---------------------------------------------------------------
        winter['fzfg'] = freezing_fog_potential(float(T[0]), float(rhi[0]), winter['sfc_wind_kt'],
                                                winter['sfc_inv_strength'])
        winter['fzdz'] = freezing_drizzle_potential(winter['clouds'], zi, rh_icei, omegai, has_omega, float(Twi[0]))
        ####################################################################



        ####################################################################
        ### WINDS ###
        ####################################################################
        #--- MEAN WINDS ---#
        # ---------------------------------------------------------------
        # PBL (from the standard calc): surface to the PBL top, the flow blowing snow and lake-effect bands follow.
        # DGZ [prototype]: the flow where snow grows, a guide to how bands generated aloft drift
        winter['mean_wind_pbl'] = snd['kinem']['mean_wind_pbl']
        winter['mean_wind_dgz'] = ((layer_mean(ui, zi, *winter['dgz_z']), layer_mean(vi, zi, *winter['dgz_z']))
                                   if has_dgz else (ma.masked, ma.masked))

        #--- WINDS AT WINTER LEVELS ---#
        # ---------------------------------------------------------------
        winter['level_winds'] = {level: {'T': float(np.interp(-level, -p, T)), 'u': float(np.interp(-level, -p, u)),
                                         'v': float(np.interp(-level, -p, v))}
                                 for level in (850, 700, 500) if p[0] >= level}
        ####################################################################

        snd['winter'] = winter
        return snd
