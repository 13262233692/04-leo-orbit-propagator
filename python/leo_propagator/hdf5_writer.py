import h5py
import numpy as np
from typing import List, Optional, Dict, Any
import os
from datetime import datetime


class HDF5Writer:
    def __init__(self, filepath: str, mode: str = "w"):
        self.filepath = filepath
        self.mode = mode
        self._file: Optional[h5py.File] = None

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def open(self):
        dirpath = os.path.dirname(os.path.abspath(self.filepath))
        if dirpath and not os.path.exists(dirpath):
            os.makedirs(dirpath, exist_ok=True)
        self._file = h5py.File(self.filepath, self.mode)

    def close(self):
        if self._file is not None:
            self._file.close()
            self._file = None

    def write_metadata(
        self,
        group_name: str = "/",
        metadata: Optional[Dict[str, Any]] = None,
    ):
        if self._file is None:
            raise RuntimeError("HDF5 file is not open")

        if metadata is None:
            metadata = {
                "created": datetime.now().isoformat(),
                "version": "1.0",
                "coordinate_system": "J2000 ECI",
                "units": "km, km/s, seconds",
            }

        if group_name == "/":
            target = self._file
        else:
            target = self._file.require_group(group_name)

        for key, value in metadata.items():
            if isinstance(value, (int, float, str, np.ndarray)):
                target.attrs[key] = value
            else:
                target.attrs[key] = str(value)

    def write_satellite_result(
        self,
        orbit_result,
        group_prefix: str = "satellites",
        compression: str = "gzip",
        compression_opts: int = 4,
    ):
        if self._file is None:
            raise RuntimeError("HDF5 file is not open")

        sat_id = orbit_result.norad_id if orbit_result.norad_id > 0 else "unknown"
        group_name = f"{group_prefix}/sat_{sat_id:06d}"
        grp = self._file.require_group(group_name)

        times = orbit_result.get_times_array()
        positions = orbit_result.get_positions_array()
        velocities = orbit_result.get_velocities_array()

        grp.create_dataset(
            "time",
            data=times,
            compression=compression,
            compression_opts=compression_opts,
        )
        grp.create_dataset(
            "position",
            data=positions,
            compression=compression,
            compression_opts=compression_opts,
        )
        grp.create_dataset(
            "velocity",
            data=velocities,
            compression=compression,
            compression_opts=compression_opts,
        )

        grp.attrs["norad_id"] = orbit_result.norad_id
        grp.attrs["satellite_name"] = orbit_result.satellite_name
        grp.attrs["num_timesteps"] = len(times)
        grp.attrs["start_time"] = times[0] if len(times) > 0 else 0.0
        grp.attrs["end_time"] = times[-1] if len(times) > 0 else 0.0

    def write_batch_results(
        self,
        orbit_results: List,
        progress_callback: Optional[callable] = None,
    ):
        if self._file is None:
            raise RuntimeError("HDF5 file is not open")

        total = len(orbit_results)
        for idx, result in enumerate(orbit_results):
            self.write_satellite_result(result)
            if progress_callback:
                progress_callback(idx + 1, total)


def save_results_to_hdf5(
    orbit_results: List,
    output_path: str,
    metadata: Optional[Dict[str, Any]] = None,
    progress_callback: Optional[callable] = None,
) -> str:
    with HDF5Writer(output_path, "w") as writer:
        writer.write_metadata("/", metadata)
        writer.write_batch_results(orbit_results, progress_callback)
    return output_path


def load_results_from_hdf5(filepath: str) -> List[Dict[str, Any]]:
    results = []
    with h5py.File(filepath, "r") as f:
        if "satellites" in f:
            sat_group = f["satellites"]
            for sat_name in sat_group:
                grp = sat_group[sat_name]
                result = {
                    "norad_id": grp.attrs.get("norad_id", 0),
                    "satellite_name": grp.attrs.get("satellite_name", ""),
                    "time": grp["time"][:],
                    "position": grp["position"][:],
                    "velocity": grp["velocity"][:],
                }
                results.append(result)
    return results
