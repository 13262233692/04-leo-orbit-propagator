import numpy as np
from typing import List, Callable, Optional, Dict, Any
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from dataclasses import dataclass
import time


@dataclass
class JobConfig:
    job_id: str
    tle_line1: str
    tle_line2: str
    satellite_name: str = ""
    duration_days: float = 7.0
    step_size: float = 1.0
    use_j2: bool = True


@dataclass
class JobResult:
    job_id: str
    norad_id: int
    satellite_name: str
    success: bool
    error_message: str = ""
    num_timesteps: int = 0
    compute_time: float = 0.0


def _process_single_job(config: JobConfig) -> tuple:
    try:
        from ._leo_propagator import TLEParser, OrbitPropagator

        start_time = time.time()
        tle = TLEParser.parse(config.tle_line1, config.tle_line2, config.satellite_name)
        propagator = OrbitPropagator(config.step_size, config.use_j2)
        result = propagator.propagate_tle(tle, config.duration_days)
        compute_time = time.time() - start_time

        return (
            config.job_id,
            True,
            "",
            result,
            compute_time,
        )
    except Exception as e:
        return (
            config.job_id,
            False,
            str(e),
            None,
            0.0,
        )


class SatelliteScheduler:
    def __init__(
        self,
        max_workers: Optional[int] = None,
        step_size: float = 1.0,
        use_j2: bool = True,
    ):
        self.max_workers = max_workers or mp.cpu_count()
        self.step_size = step_size
        self.use_j2 = use_j2
        self._job_queue: List[JobConfig] = []
        self._results_cache: Dict[str, Any] = {}

    def add_job(
        self,
        tle_line1: str,
        tle_line2: str,
        satellite_name: str = "",
        duration_days: float = 7.0,
        job_id: Optional[str] = None,
    ) -> str:
        if job_id is None:
            job_id = f"job_{len(self._job_queue):06d}"

        config = JobConfig(
            job_id=job_id,
            tle_line1=tle_line1,
            tle_line2=tle_line2,
            satellite_name=satellite_name,
            duration_days=duration_days,
            step_size=self.step_size,
            use_j2=self.use_j2,
        )
        self._job_queue.append(config)
        return job_id

    def add_jobs_from_tle_list(self, tle_list: List[tuple]) -> List[str]:
        job_ids = []
        for idx, (line1, line2, name) in enumerate(tle_list):
            jid = self.add_job(line1, line2, name)
            job_ids.append(jid)
        return job_ids

    def run_all(self, progress_callback: Optional[Callable] = None) -> List[JobResult]:
        results = []
        orbit_results = []
        completed = 0
        total = len(self._job_queue)

        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_job = {
                executor.submit(_process_single_job, config): config
                for config in self._job_queue
            }

            for future in as_completed(future_to_job):
                config = future_to_job[future]
                try:
                    job_id, success, error, orbit_result, compute_time = future.result()

                    job_result = JobResult(
                        job_id=job_id,
                        norad_id=0,
                        satellite_name=config.satellite_name,
                        success=success,
                        error_message=error,
                        compute_time=compute_time,
                    )

                    if success and orbit_result is not None:
                        job_result.norad_id = orbit_result.norad_id
                        job_result.num_timesteps = len(orbit_result.times)
                        orbit_results.append(orbit_result)
                        self._results_cache[job_id] = orbit_result

                    results.append(job_result)
                    completed += 1

                    if progress_callback:
                        progress_callback(completed, total, job_result)

                except Exception as e:
                    results.append(
                        JobResult(
                            job_id=config.job_id,
                            norad_id=0,
                            satellite_name=config.satellite_name,
                            success=False,
                            error_message=str(e),
                        )
                    )
                    completed += 1
                    if progress_callback:
                        progress_callback(completed, total, results[-1])

        self._job_queue.clear()
        self._last_orbit_results = orbit_results
        return results

    def get_result(self, job_id: str):
        return self._results_cache.get(job_id)

    def get_all_orbit_results(self) -> List:
        return getattr(self, "_last_orbit_results", [])

    def clear(self):
        self._job_queue.clear()
        self._results_cache.clear()
        if hasattr(self, "_last_orbit_results"):
            del self._last_orbit_results

    @property
    def pending_jobs(self) -> int:
        return len(self._job_queue)


def generate_mock_tles(num_satellites: int, base_altitude_km: float = 550.0) -> List[tuple]:
    tles = []

    for i in range(num_satellites):
        incl = 53.0 + np.random.uniform(-0.5, 0.5)
        raan = (i * 360.0 / num_satellites) % 360.0
        ecc = 0.0001 + np.random.uniform(0, 0.001)
        argp = np.random.uniform(0, 360.0)
        mean_anomaly = np.random.uniform(0, 360.0)

        earth_radius = 6378.137
        mu = 398600.4418
        a = earth_radius + base_altitude_km + np.random.uniform(-10, 10)
        n = np.sqrt(mu / (a ** 3)) * 86400.0 / (2 * np.pi)

        year = 24
        day_of_year = 150.0 + np.random.uniform(0, 10)

        norad_id = 50000 + i

        line1 = (
            f"1 {norad_id:05d}U 24001{(i % 100):02d}  "
            f"{year:02d}{day_of_year:012.8f}  "
            f".00000000  00000-0  10000-3 0  0001"
        )

        line2 = (
            f"2 {norad_id:05d} {incl:8.4f} {raan:8.4f} "
            f"{int(ecc * 10000000):07d} {argp:8.4f} {mean_anomaly:8.4f} "
            f"{n:11.8f}{1:05d}"
        )

        line1 = line1[:68]
        checksum1 = sum(int(c) for c in line1 if c.isdigit()) % 10
        line1 = line1[:68] + str(checksum1)

        line2 = line2[:68]
        checksum2 = sum(int(c) for c in line2 if c.isdigit()) % 10
        line2 = line2[:68] + str(checksum2)

        name = f"SAT-{i:04d}"
        tles.append((line1, line2, name))

    return tles
