# __init__.py
from .sounderpy import get_model_data 
from .sounderpy import get_obs_data 
from .sounderpy import get_bufkit_data 
from .sounderpy import acars_data
from .sounderpy import build_sounding
from .sounderpy import build_hodograph
from .sounderpy import build_vad_hodograph
from .sounderpy import build_composite
from .sounderpy import print_variables
from .sounderpy import from_cm1
from .sounderpy import from_wrf
from .sounderpy import make_cm1_profile
from .sounderpy import make_wrf_profile
from .utils import interp_data
from .utils import find_nearest
from .utils import get_sfc_index
from .utils import make_sfc_based
from .utils import make_sfc_based_3D
from .utils import get_latlon
from .data.file_io import to_file
from .data.file_io import from_sharppy
from .data.file_io import from_graw
from .data.file_io import from_imet
from .data.file_io import from_windsond
from .profile_utils import merge_profiles
from .profile_utils import smooth_profile
from .calc import sounding_params
from .calc import vad_params
from .calc_winter import winter_params
from .plot.style import set_watermark
from .data.rapruc_grib_parser import parse_rapruc_grib
from .data.rapruc_grib_retriever import get_rapruc_data


# OLD MODULE PATHS
# the data modules moved into sounderpy.data in 3.3.0. The old paths (e.g. `from sounderpy.bufkit_data
# import fetch_bufkit`) still work through these stand-ins, which show a FutureWarning when used.
# Remove this block, and the old paths with it, in a later release.
import sys as _sys
import types as _types
import warnings as _warnings
import importlib as _importlib

class _MovedModule(_types.ModuleType):
    def __init__(self, old_name, new_name):
        super().__init__(old_name)
        self._new_name = new_name

    def __getattr__(self, name):
        if name.startswith('__') and name != '__all__':
            raise AttributeError(name)
        _warnings.warn(f"'{self.__name__}' moved to '{self._new_name}' in SounderPy 3.3.0. The old path "
                       f"will be removed in a later release; import from '{self._new_name}' instead.",
                       FutureWarning, stacklevel=2)
        module = _importlib.import_module(self._new_name)
        if name == '__all__':
            return getattr(module, '__all__', [key for key in vars(module) if not key.startswith('_')])
        return getattr(module, name)

for _name in ('acars_data', 'bufkit_data', 'ncep_bufr', 'obs_data', 'model_data', 'rapruc_grib_retriever',
              'rapruc_grib_parser', 'radar_utils', 'wrf_utils', 'cm1_utils', 'file_io'):
    _moved = _MovedModule(f'{__name__}.{_name}', f'{__name__}.data.{_name}')
    _sys.modules.setdefault(_moved.__name__, _moved)
    # `sounderpy.acars_data` stays the acars_data class
    if _name not in globals():
        globals()[_name] = _moved
del _name, _moved
