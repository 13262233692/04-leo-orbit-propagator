#!/usr/bin/env python3
"""
单颗卫星轨道推演示例
使用真实的 ISS (国际空间站) TLE 数据进行测试
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from leo_propagator import (
    TLEParser,
    OrbitPropagator,
    Vector3,
    is_cpp_available,
)


def main():
    print("单颗卫星轨道推演示例")
    print("=" * 50)

    if not is_cpp_available():
        print("错误: C++ 模块未构建，请先运行 build.bat")
        return 1

    line1 = "1 25544U 98067A   24150.59722222  .00016717  00000-0  10270-3 0  9993"
    line2 = "2 25544  51.6416 195.1059 0006763  85.4663 274.7635 15.49559341  1234"

    print("\n解析 TLE 数据...")
    tle = TLEParser.parse(line1, line2, "ISS (ZARYA)")

    print(f"\n卫星信息:")
    print(f"  名称: {tle.name}")
    print(f"  NORAD ID: {tle.norad_id}")
    print(f"  半长轴: {tle.semi_major_axis:.2f} km")
    print(f"  偏心率: {tle.eccentricity:.6f}")
    print(f"  轨道倾角: {tle.inclination * 180 / 3.14159:.2f}°")
    print(f"  轨道周期: {tle.period / 60:.2f} 分钟")
    print(f"  轨道高度: ~{tle.semi_major_axis - 6378.137:.1f} km")

    print(f"\n转换为初始状态向量 (J2000 ECI)...")
    initial_state = TLEParser.to_initial_state(tle)
    print(f"  初始位置: [{initial_state.position.x:.2f}, {initial_state.position.y:.2f}, {initial_state.position.z:.2f}] km")
    print(f"  初始速度: [{initial_state.velocity.x:.4f}, {initial_state.velocity.y:.4f}, {initial_state.velocity.z:.4f}] km/s")

    print(f"\n创建轨道推演器 (RK4, 步长 10s, J2 摄动)...")
    propagator = OrbitPropagator(step_size=10.0, use_j2=True)

    duration_hours = 6
    duration_seconds = duration_hours * 3600
    print(f"\n推演 {duration_hours} 小时的轨道...")

    import time
    start = time.time()
    result = propagator.propagate_state(initial_state, duration_seconds)
    elapsed = time.time() - start

    print(f"\n推演完成!")
    print(f"  时间步数: {len(result.times)}")
    print(f"  推演耗时: {elapsed:.4f} 秒")
    print(f"  每秒推演步数: {len(result.times) / elapsed:.0f}")

    positions = result.get_positions_array()
    velocities = result.get_velocities_array()
    times = result.get_times_array()

    altitudes = np.linalg.norm(positions, axis=1) - 6378.137

    print(f"\n轨道统计:")
    print(f"  近地点高度: {altitudes.min():.2f} km")
    print(f"  远地点高度: {altitudes.max():.2f} km")
    print(f"  平均高度: {altitudes.mean():.2f} km")

    print(f"\n末态 (t = {times[-1]:.0f}s):")
    print(f"  位置: [{positions[-1,0]:.2f}, {positions[-1,1]:.2f}, {positions[-1,2]:.2f}] km")
    print(f"  速度: [{velocities[-1,0]:.4f}, {velocities[-1,1]:.4f}, {velocities[-1,2]:.4f}] km/s")

    print("\n" + "=" * 50)
    print("验证完成!")

    return 0


if __name__ == "__main__":
    sys.exit(main())
