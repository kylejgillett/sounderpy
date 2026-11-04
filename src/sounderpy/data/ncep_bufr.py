from datetime import datetime, timedelta, timezone
from urllib.error import HTTPError
from urllib.request import urlopen

import numpy as np




"""
    SOUNDERPY NCEP BUFR SOUNDING DECODER

    Purpose of module:

    Decode NCEP model BUFR sounding files (the 'bufrsnd' station files NCEP
    builds BUFKIT files from) and write them out as BUFKIT text, so that
    bufkit_data.py can read them like any other BUFKIT file.

    NCEP BUFR files carry their own descriptor tables (DX tables) in the
    first few messages, so no outside BUFR tables or libraries are needed.


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
"""



###############
# STANDARD DESCRIPTORS
#########################################################################
'''
- the DX table messages are written with these standard WMO descriptors
- fxy: (mnemonic, units, scale, reference, width in bits)
'''
STD_B = {
    '000001': ('TABLAT', 'CCITT IA5', 0, 0, 24),
    '000002': ('TABLAL1', 'CCITT IA5', 0, 0, 256),
    '000003': ('TABLAL2', 'CCITT IA5', 0, 0, 256),
    '000010': ('FDESC', 'CCITT IA5', 0, 0, 8),
    '000011': ('XDESC', 'CCITT IA5', 0, 0, 16),
    '000012': ('YDESC', 'CCITT IA5', 0, 0, 24),
    '000013': ('ELEMNA1', 'CCITT IA5', 0, 0, 256),
    '000014': ('ELEMNA2', 'CCITT IA5', 0, 0, 256),
    '000015': ('UNITS', 'CCITT IA5', 0, 0, 192),
    '000016': ('SCLSGN', 'CCITT IA5', 0, 0, 8),
    '000017': ('SCALE', 'CCITT IA5', 0, 0, 24),
    '000018': ('REFSGN', 'CCITT IA5', 0, 0, 8),
    '000019': ('REFVAL', 'CCITT IA5', 0, 0, 80),
    '000020': ('WIDTH', 'CCITT IA5', 0, 0, 24),
    '000030': ('DSEQ', 'CCITT IA5', 0, 0, 48),
    '031000': ('DRF1BIT', 'NUMERIC', 0, 0, 1),
    '031001': ('DRF8BIT', 'NUMERIC', 0, 0, 8),
    '031002': ('DRF16BIT', 'NUMERIC', 0, 0, 16),
}
STD_D = {
    '300003': ['000010', '000011', '000012'],
    '300004': ['300003', '000013', '000014', '000015', '000016', '000017',
               '000018', '000019', '000020'],
}

NOMADS_RRFS = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/rrfs/v1.0'
RRFS_SOUNDING_CYCLES = 6      # NCEP only writes RRFS soundings for the 00, 06, 12 & 18Z runs
NOMADS_RRFS_LATENCY = 3       # hours after a run starts before its soundings are usually posted
NOMADS_TIMEOUT = 20




###############
# BUFR MESSAGE DECODING
#########################################################################

class _Bits:
    # read big-endian bit fields from a byte string
    def __init__(self, buf):
        self.buf = buf
        self.pos = 0

    def read(self, width):
        start, end = self.pos, self.pos + width
        if end > 8 * len(self.buf):
            raise ValueError('BUFR data section ended early')
        chunk = int.from_bytes(self.buf[start // 8:(end + 7) // 8], 'big')
        self.pos = end
        return (chunk >> ((8 - end % 8) % 8)) & ((1 << width) - 1)

    def text(self, width):
        return self.read(width).to_bytes(width // 8, 'big').decode('ascii', 'replace')



def split_messages(data):
    # NCEP files wrap each message in Fortran record markers, so find each 'BUFR'
    msgs, i = [], 0
    while True:
        j = data.find(b'BUFR', i)
        if j < 0:
            return msgs
        n = int.from_bytes(data[j + 4:j + 7], 'big')
        if data[j + n - 4:j + n] != b'7777':
            raise ValueError(f'BUFR message at byte {j} is incomplete')
        msgs.append(data[j:j + n])
        i = j + n



def parse_sections(msg):
    # pull the header info, descriptors and data section from one message
    edition = msg[7]
    if edition not in (3, 4):
        raise ValueError(f'BUFR edition {edition} is not supported')
    p = 8
    l1 = int.from_bytes(msg[p:p + 3], 'big')
    if edition == 4:
        category, has_opt = msg[p + 10], msg[p + 9] & 0x80
        date = (int.from_bytes(msg[p + 15:p + 17], 'big'), *msg[p + 17:p + 20])
    else:
        category, has_opt = msg[p + 8], msg[p + 7] & 0x80
        date = (2000 + msg[p + 12], *msg[p + 13:p + 16])
    p += l1
    if has_opt:
        p += int.from_bytes(msg[p:p + 3], 'big')
    l3 = int.from_bytes(msg[p:p + 3], 'big')
    nsub = int.from_bytes(msg[p + 4:p + 6], 'big')
    compressed = bool(msg[p + 6] & 0x40)
    raw = msg[p + 7:p + l3]
    descs = [f'{raw[k] >> 6}{raw[k] & 63:02d}{raw[k + 1]:03d}' for k in range(0, len(raw) - 1, 2)]
    p += l3
    l4 = int.from_bytes(msg[p:p + 3], 'big')
    return {'category': category, 'date': date, 'nsub': nsub,
            'compressed': compressed, 'descs': descs, 'data': msg[p + 4:p + l4]}



class Decoder:
    # decode NCEP BUFR messages using the DX tables found in the file

    def __init__(self):
        self.b = dict(STD_B)
        self.d = {k: list(v) for k, v in STD_D.items()}


    def walk(self, descs, bits, out):
        # expand descriptors in order, appending (mnemonic, value) pairs to out
        i = 0
        width_add, scale_add, next_width = 0, 0, None
        while i < len(descs):
            fxy = descs[i]
            f, x, y = int(fxy[0]), int(fxy[1:3]), int(fxy[3:])
            if f == 0:
                if next_width is not None and fxy not in self.b:
                    bits.read(next_width)
                else:
                    out.append(self._value(fxy, bits, width_add, scale_add))
                next_width = None
                i += 1
            elif f == 1:
                if y == 0:
                    # delayed replication: the count comes from the next descriptor
                    factor = self._value(descs[i + 1], bits, 0, 0)
                    out.append(factor)
                    count, group = factor[1], descs[i + 2:i + 2 + x]
                    i += 2 + x
                else:
                    count, group = y, descs[i + 1:i + 1 + x]
                    i += 1 + x
                for _ in range(count or 0):
                    self.walk(group, bits, out)
            elif f == 2:
                if x == 1:
                    width_add = y - 128 if y else 0
                elif x == 2:
                    scale_add = y - 128 if y else 0
                elif x == 5:
                    out.append(('TEXT', bits.text(8 * y)))
                elif x == 6:
                    next_width = y
                else:
                    raise NotImplementedError(f'BUFR operator {fxy} is not supported')
                i += 1
            else:
                if fxy not in self.d:
                    raise KeyError(f'BUFR sequence {fxy} is not in the file tables')
                # splice the sequence in place: NCEP sequences can be a bare replication
                # operator that applies to the descriptor after the sequence
                descs = descs[:i] + self.d[fxy] + descs[i + 1:]
        return out


    def _value(self, fxy, bits, width_add, scale_add):
        if fxy not in self.b:
            raise KeyError(f'BUFR element {fxy} is not in the file tables')
        mnem, units, scale, ref, width = self.b[fxy]
        if units.startswith('CCITT'):
            return (mnem, bits.text(width).strip())
        width += width_add
        scale += scale_add
        raw = bits.read(width)
        # all bits set means missing, except for replication counts and 1-bit flags
        if width > 1 and raw == (1 << width) - 1 and not mnem.startswith('DRF'):
            return (mnem, None)
        return (mnem, (raw + ref) / 10 ** scale if scale else raw + ref)


    def load_dx(self, sec):
        # read Table A, B and D entries from a DX table message
        vals = [v for _, v in self.walk(sec['descs'], _Bits(sec['data']), [])]
        k = 0
        def nxt():
            nonlocal k
            k += 1
            return vals[k - 1]
        for _ in range(nxt()):
            nxt(), nxt(), nxt()
        for _ in range(nxt()):
            f, x, y = nxt(), nxt(), nxt()
            name, _, units = nxt(), nxt(), nxt()
            ssign, scale, rsign, ref, width = nxt(), nxt(), nxt(), nxt(), nxt()
            scale = int(scale) * (-1 if ssign == '-' else 1)
            ref = int(ref) * (-1 if rsign == '-' else 1)
            self.b[f'{f}{x}{y}'] = (name.split()[0], units, scale, ref, int(width))
        for _ in range(nxt()):
            f, x, y = nxt(), nxt(), nxt()
            nxt()
            self.d[f'{f}{x}{y}'] = [nxt() for _ in range(nxt())]


    def subsets(self, sec):
        # yield the (mnemonic, value) list of each subset in a data message
        if sec['compressed']:
            raise NotImplementedError('compressed BUFR messages are not supported')
        bits = _Bits(sec['data'])
        for _ in range(sec['nsub']):
            start = bits.pos
            vals = self.walk(sec['descs'], bits, [])
            # each NCEP subset starts with its length in bytes; use it to stay aligned
            if vals and vals[0][0] == 'BYTCNT' and vals[0][1]:
                bits.pos = start + 8 * int(vals[0][1])
            yield vals



def decode_soundings(data):
    # decode a model sounding BUFR file into the run time and one dict per forecast hour
    msgs = [parse_sections(m) for m in split_messages(data)]
    dec = Decoder()
    profiles, run_time = [], None
    for sec in msgs:
        if sec['category'] == 11:
            dec.load_dx(sec)
            continue
        run_time = run_time or datetime(*sec['date'])
        for vals in dec.subsets(sec):
            prof = {'levels': []}
            for name, val in vals:
                if name == 'PRES':
                    prof['levels'].append({})
                if name in ('PRES', 'TMDB', 'UWND', 'VWND', 'SPFH', 'OMEG', 'CFRL') and prof['levels']:
                    prof['levels'][-1].setdefault(name, val)
                elif name in ('SMC1', 'STC1'):
                    prof.setdefault(name, []).append(val)
                else:
                    prof.setdefault(name, val)
            profiles.append(prof)
    if run_time is None or not profiles:
        raise ValueError('no soundings were found in the BUFR file')
    return run_time, profiles




###############
# BUFKIT TEXT OUTPUT
#########################################################################

def _fmt(vals):
    return ' '.join('-9999.00' if v is None or not np.isfinite(v) else f'{v:.2f}' for v in vals)



def _arr(levels, key):
    return np.array([np.nan if lvl.get(key) is None else lvl[key] for lvl in levels], dtype=float)



def _dewpoint(q, p_hpa):
    # dewpoint (C) from specific humidity (kg/kg) and pressure (hPa), Bolton (1980)
    e = q * p_hpa / (0.622 + 0.378 * q)
    with np.errstate(divide='ignore', invalid='ignore'):
        lne = np.log(e / 6.112)
        return 243.5 * lne / (17.67 - lne)



def _wetbulb(t, td, p_hpa):
    # wet-bulb temperature (C) from the psychrometric equation, solved with Newton's method
    e = 6.112 * np.exp(17.67 * td / (td + 243.5))
    tw = td + (t - td) / 3
    for _ in range(20):
        es = 6.112 * np.exp(17.67 * tw / (tw + 243.5))
        des = es * 17.67 * 243.5 / (tw + 243.5) ** 2
        gam = 6.6e-4 * p_hpa
        f = es - e - gam * (t - tw)
        tw = tw - f / (des + gam)
    return tw



def _theta_e(t, td, p_hpa):
    # equivalent potential temperature (K), Bolton (1980)
    tk, tdk = t + 273.15, td + 273.15
    e = 6.112 * np.exp(17.67 * td / (td + 243.5))
    r = 0.622 * e / (p_hpa - e)
    tl = 1 / (1 / (tdk - 56) + np.log(tk / tdk) / 800) + 56
    th = tk * (1000 / p_hpa) ** (0.2854 * (1 - 0.28 * r))
    return th * np.exp((3.376 / tl - 0.00254) * 1000 * r * (1 + 0.81 * r))



def _heights(p_hpa, t, q, p_sfc, z_sfc):
    # hydrostatic heights (m) of each level, integrated up from the surface
    tv = (t + 273.15) * (1 + 0.61 * np.nan_to_num(q))
    z = np.empty_like(p_hpa)
    z[0] = z_sfc + 287.04 * tv[0] / 9.80665 * np.log(p_sfc / p_hpa[0])
    for k in range(1, len(p_hpa)):
        z[k] = z[k - 1] + 287.04 * (tv[k - 1] + tv[k]) / 2 / 9.80665 * np.log(p_hpa[k - 1] / p_hpa[k])
    return z



def _k_to_c(v):
    return None if v is None else v - 273.15



def _scale(v, factor):
    return None if v is None else v * factor



def to_bufkit_lines(run_time, profiles):
    # write decoded soundings out as the lines of a BUFKIT text file
    lines = ['',
             'SNPARM = PRES;TMPC;TMWC;DWPC;THTE;DRCT;SKNT;OMEG;CFRL;HGHT',
             'STNPRM = SHOW;LIFT;SWET;KINX;LCLP;PWAT;TOTL;CAPE;LCLT;CINS;EQLV;LFCT;BRCH',
             '']
    sfc_rows = []
    for prof in profiles:
        valid = run_time + timedelta(seconds=prof['FTIM'])
        stamp = valid.strftime('%y%m%d/%H%M')
        lvls = prof['levels']
        p = _arr(lvls, 'PRES') / 100
        t = _arr(lvls, 'TMDB') - 273.15
        q = _arr(lvls, 'SPFH')
        u, v = _arr(lvls, 'UWND'), _arr(lvls, 'VWND')
        td = _dewpoint(q, p)
        drct = np.degrees(np.arctan2(-u, -v)) % 360
        sknt = np.hypot(u, v) * 1.943844
        p_sfc = prof['PRSS'] / 100 if prof.get('PRSS') else p[0]
        hght = _heights(p, t, q, p_sfc, prof['GELV'])

        lines += [f"STID = {prof['RPID']} STNM = {int(prof['STNM'])} TIME = {stamp}",
                  f"SLAT = {prof['CLAT']:.2f} SLON = {prof['CLON']:.2f} SELV = {float(prof['GELV']):.1f}",
                  f"STIM = {int(prof['FTIM'] // 3600)}",
                  '',
                  # the stability indices are left missing; sounderpy computes its own
                  'SHOW = -9999.00 LIFT = -9999.00 SWET = -9999.00 KINX = -9999.00',
                  'LCLP = -9999.00 PWAT = -9999.00 TOTL = -9999.00 CAPE = -9999.00',
                  'LCLT = -9999.00 CINS = -9999.00 EQLV = -9999.00 LFCT = -9999.00',
                  'BRCH = -9999.00',
                  '',
                  'PRES TMPC TMWC DWPC THTE DRCT SKNT OMEG',
                  'CFRL HGHT']
        cols = [p, t, _wetbulb(t, td, p), td, _theta_e(t, td, p), drct, sknt,
                _arr(lvls, 'OMEG'), _arr(lvls, 'CFRL'), hght]
        for k in range(len(p)):
            row = [c[k] for c in cols]
            lines += [_fmt(row[:8]), _fmt(row[8:])]
        lines.append('')

        soil_t = prof.get('STC1') or []
        q2 = prof.get('Q2MS')
        sfc_rows.append([
            prof['STNM'], stamp,
            _scale(prof.get('PMSL'), 0.01), _scale(prof.get('PRSS'), 0.01), _k_to_c(prof.get('TMSK')),
            soil_t[0] if soil_t else None, prof.get('SNFL'), prof.get('WTNS'),
            prof.get('TP01'), prof.get('C01M'), soil_t[1] if len(soil_t) > 1 else None,
            prof.get('LCLD'), prof.get('MCLD'), prof.get('HCLD'), prof.get('SNRA'),
            prof.get('U10M'), prof.get('V10M'), prof.get('R01M'), prof.get('BFGR'),
            _k_to_c(prof.get('T2MS')), _scale(q2, 1000), prof.get('WXTS'), prof.get('WXTP'),
            prof.get('WXTZ'), prof.get('WXTR'), prof.get('USTM'), prof.get('VSTM'),
            prof.get('HLCY'), None, None, _scale(prof.get('CDBP'), 0.01),
            _scale(prof.get('HOVI'), 0.001),
            None if q2 is None or not prof.get('PRSS') else float(_dewpoint(np.array(q2), np.array(prof['PRSS'] / 100)))])

    lines += ['STN YYMMDD/HHMM PMSL PRES SKTC STC1 SNFL WTNS',
              'P01M C01M STC2 LCLD MCLD HCLD',
              'SNRA UWND VWND R01M BFGR T2MS',
              'Q2MS WXTS WXTP WXTZ WXTR USTM',
              'VSTM HLCY SLLH WSYM CDBP VSBK',
              'TD2M']
    for row in sfc_rows:
        vals = row[2:]
        lines += [f'{int(row[0])} {row[1]} ' + _fmt(vals[:6]), _fmt(vals[6:12]), _fmt(vals[12:18]),
                  _fmt(vals[18:24]), _fmt(vals[24:30]), _fmt(vals[30:])]
    return lines




###############
# NOMADS RRFS ACCESS
#########################################################################

def nomads_rrfs_url(stnm, run_dt):
    return (f'{NOMADS_RRFS}/rrfs.{run_dt:%Y%m%d}/{run_dt:%H}/bufrsnd.t{run_dt:%H}z/'
            f'bufr.{int(stnm)}.{run_dt:%Y%m%d%H}')



def latest_expected_rrfs_run(now=None):
    # the newest RRFS run whose soundings should be posted by now
    now = now or datetime.now(timezone.utc).replace(tzinfo=None)
    t = now - timedelta(hours=NOMADS_RRFS_LATENCY)
    return t.replace(hour=t.hour - t.hour % RRFS_SOUNDING_CYCLES, minute=0, second=0, microsecond=0)



def fetch_nomads_rrfs(stnm, run_dt):
    # download one station's RRFS BUFR soundings from NOMADS and return BUFKIT text lines
    with urlopen(nomads_rrfs_url(stnm, run_dt), timeout=NOMADS_TIMEOUT) as resp:
        data = resp.read()
    run_time, profiles = decode_soundings(data)
    return to_bufkit_lines(run_time, profiles)



def newest_nomads_rrfs(stnm, newer_than=None, tries=8):
    # walk back from the newest expected run until NOMADS has one (NOMADS keeps about two days)
    run_dt = latest_expected_rrfs_run()
    for _ in range(tries):
        if newer_than is not None and run_dt <= newer_than:
            return None, None
        try:
            return run_dt, fetch_nomads_rrfs(stnm, run_dt)
        except HTTPError as e:
            if e.code not in (403, 404):
                raise
        except ValueError:
            # the file may still be uploading; try the run before it
            pass
        run_dt -= timedelta(hours=RRFS_SOUNDING_CYCLES)
    return None, None
