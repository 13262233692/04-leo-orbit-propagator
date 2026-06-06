#!/usr/bin/env python3
"""
LEO Orbit Propagator - 批量卫星轨道推演示例

演示如何:
1. 生成模拟的卫星 TLE 数据
2. 使用批量作业调度 API 并行推演 1000 颗卫星
3. 将结果保存到 HDF5 文件
4. 读取并验证结果
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from leo_propagator import (
    SatelliteScheduler,
    generate_mock_tles,
    save_results_to_hdf5,
    load_results_from_hdf5,
    is_cpp_available,
)


def main():
    print("=" * 70)
    print("LEO 卫星星座轨道高精度仿真引擎 - 批量推演示例")
    print("=" * 70)

    print(f"\n[信息] C++ 核心模块状态: {'可用' if is_cpp_available() else '不可用'}")
    if not is_cpp_available():
        print("[错误] 请先运行 build.bat 构建 C++ 扩展模块")
        return 1

    num_satellites = 1000
    duration_days = 7.0
    step_size = 60.0
    output_file = "output/leo_constellation_7days.h5"

    print(f"\n[配置]")
    print(f"  - 卫星数量: {num_satellites}")
    print(f"  - 推演时长: {duration_days} 天")
    print(f"  - 时间步长: {step_size} 秒")
    print(f"  - 输出文件: {output_file}")

    print(f"\n[1/5] 生成模拟 TLE 数据...")
    start_time = time.time()
    tles = generate_mock_tles(num_satellites, base_altitude_km=550.0)
    print(f"  生成 {len(tles)} 颗卫星的 TLE 数据，耗时 {time.time()-start_time:.2f}s")

    print(f"\n[2/5] 初始化批量作业调度器...")
    scheduler = SatelliteScheduler(
        step_size=step_size,
        use_j2=True,
    )

    print(f"\n[3/5] 提交推演任务...")
    job_ids = scheduler.add_jobs_from_tle_list(tles)
    print(f"  已提交 {len(job_ids)} 个推演任务")

    print(f"\n[4/5] 并行执行推演 (7 天 × 1000 颗卫星)...")
    print(f"  预计总时间步数: 每颗卫星约 {int(duration_days * 86400 / step_size)} 步")

    start_time = time.time()

    def progress_callback(completed, total, job_result=None):
        if completed % 50 == 0 or completed == total:
            elapsed = time.time() - start_time
            rate = completed / elapsed if elapsed > 0 else 0
            eta = (total - completed) / rate if rate > 0 else 0
            print(f"  进度: {completed}/{total} ({completed/total*100:.1f}%) | "
                  f"速率: {rate:.1f} 颗/秒 | 预计剩余: {eta:.0f}s")

    results = scheduler.run_all(progress_callback=progress_callback)
    total_time = time.time() - start_time

    success_count = sum(1 for r in results if r.success)
    failed_count = len(results) - success_count

    print(f"\n  推演完成!")
    print(f"  成功: {success_count} 颗, 失败: {failed_count} 颗")
    print(f"  总耗时: {total_time:.2f} 秒")
    print(f"  平均单颗: {total_time/num_satellites*1000:.1f} 毫秒")
    print(f"  吞吐量: {num_satellites/total_time:.2f} 颗/秒")

    if success_count == 0:
        print("\n[错误] 没有成功的推演结果")
        return 1

    print(f"\n[5/5] 序列化结果到 HDF5 文件...")
    orbit_results = scheduler.get_all_orbit_results()

    metadata = {
        "num_satellites": num_satellites,
        "duration_days": duration_days,
        "step_size_seconds": step_size,
        "include_j2": True,
        "coordinate_system": "J2000 ECI (Earth Centered Inertial)",
        "position_units": "km",
        "velocity_units": "km/s",
        "time_units": "seconds from epoch",
    }

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    def hdf5_progress(completed, total):
        if completed % 100 == 0 or completed == total:
            print(f"  写入进度: {completed}/{total}")

    save_path = save_results_to_hdf5(
        orbit_results,
        output_file,
        metadata=metadata,
        progress_callback=hdf5_progress,
    )

    file_size = os.path.getsize(save_path) / (1024 * 1024)
    print(f"\n  HDF5 文件已保存: {save_path}")
    print(f"  文件大小: {file_size:.2f} MB")

    print(f"\n[验证] 读取 HDF5 并验证数据...")
    loaded = load_results_from_hdf5(save_path)
    print(f"  成功读取 {len(loaded)} 颗卫星的数据")

    if loaded:
        sample = loaded[0]
        print(f"\n  示例卫星 (NORAD ID: {sample['norad_id']}):")
        print(f"    名称: {sample['satellite_name']}")
        print(f"    时间步数: {len(sample['time'])}")
        print(f"    时间范围: {sample['time'][0]:.0f}s - {sample['time'][-1]:.0f}s")
        print(f"    初始位置: [{sample['position'][0,0]:.2f}, {sample['position'][0,1]:.2f}, {sample['position'][0,2]:.2f}] km")
        print(f"    初始速度: [{sample['velocity'][0,0]:.4f}, {sample['velocity'][0,1]:.4f}, {sample['velocity'][0,2]:.4f}] km/s")

        altitudes = np.linalg.norm(sample['position'], axis=1) - 6378.137
        print(f"    轨道高度范围: {altitudes.min():.1f} - {altitudes.max():.1f} km")

    print("\n" + "=" * 70)
    print("批量推演任务全部完成!")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    sys.exit(main())
