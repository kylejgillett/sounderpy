"""
    SOUNDERPY PLOTTING

    Purpose of module:

    The modules that draw SounderPy's figures:

    style.py             : themes, colors, units, layouts, and watermark settings
    components.py        : reusable pieces that figures are built from
    components_winter.py : the winter sounding's own pieces
    operational.py       : full (severe and winter) sounding, hodograph, composite, and VAD figures

    Docs:
        plotting data:  https://sounderpy.app/docs/reference/plottingdata.html

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
"""

from .operational import __full_sounding, __winter_sounding, _resolve_sounding_mode, __full_hodograph, __composite_sounding, __vad_hodograph

