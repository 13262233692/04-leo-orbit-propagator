import os
import sys

_ext_path = os.path.dirname(os.path.abspath(__file__))
if _ext_path not in sys.path:
    sys.path.insert(0, _ext_path)

try:
    from ._leo_propagator import (
        Vector3,
        State,
        TLE,
        TLEParser,
        GravityModel,
        RK4Integrator,
        HighPrecisionRK4,
        RK45Integrator,
        OrbitResult,
        OrbitPropagator,
        BatchPropagator,
        KahanSum,
        IntegratorType,
    )
    _CPP_AVAILABLE = True
except ImportError:
    _CPP_AVAILABLE = False

from .scheduler import SatelliteScheduler, generate_mock_tles
from .hdf5_writer import HDF5Writer, save_results_to_hdf5

__version__ = "2.0.0"
__all__ = [
    "Vector3",
    "State",
    "TLE",
    "TLEParser",
    "GravityModel",
    "RK4Integrator",
    "HighPrecisionRK4",
    "RK45Integrator",
    "OrbitResult",
    "OrbitPropagator",
    "BatchPropagator",
    "KahanSum",
    "IntegratorType",
    "SatelliteScheduler",
    "HDF5Writer",
    "save_results_to_hdf5",
    "generate_mock_tles",
]


def is_cpp_available():
    return _CPP_AVAILABLE
