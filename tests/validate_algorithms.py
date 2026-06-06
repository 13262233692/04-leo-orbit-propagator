#!/usr/bin/env python3
"""
核心算法验证脚本 - 不依赖 C++ 编译

验证内容:
1. RK4 积分器对简单谐振子的准确性
2. J2 摄动加速度计算
3. 轨道能量守恒验证
"""

import numpy as np
import math

PI = math.pi
EARTH_MU = 398600.4418
EARTH_RADIUS = 6378.137
EARTH_J2 = 1.082635854e-3


def rk4_step(state, derivative, dt):
    k1 = derivative(state)
    k1_pos = np.array(k1[0]) * dt
    k1_vel = np.array(k1[1]) * dt

    state2 = (state[0] + k1_pos * 0.5, state[1] + k1_vel * 0.5)
    k2 = derivative(state2)
    k2_pos = np.array(k2[0]) * dt
    k2_vel = np.array(k2[1]) * dt

    state3 = (state[0] + k2_pos * 0.5, state[1] + k2_vel * 0.5)
    k3 = derivative(state3)
    k3_pos = np.array(k3[0]) * dt
    k3_vel = np.array(k3[1]) * dt

    state4 = (state[0] + k3_pos, state[1] + k3_vel)
    k4 = derivative(state4)
    k4_pos = np.array(k4[0]) * dt
    k4_vel = np.array(k4[1]) * dt

    new_pos = state[0] + (k1_pos + 2*k2_pos + 2*k3_pos + k4_pos) / 6.0
    new_vel = state[1] + (k1_vel + 2*k2_vel + 2*k3_vel + k4_vel) / 6.0

    return (new_pos, new_vel, state[2] + dt)


def two_body_acceleration(position):
    r = np.linalg.norm(position)
    return -EARTH_MU / (r**3) * position


def j2_acceleration(position):
    r = np.linalg.norm(position)
    r2 = r * r
    r5 = r2 * r2 * r
    z2 = position[2] * position[2]
    ratio = z2 / r2
    factor = -1.5 * EARTH_J2 * EARTH_MU * EARTH_RADIUS**2 / r5

    ax = factor * position[0] * (5.0 * ratio - 1.0)
    ay = factor * position[1] * (5.0 * ratio - 1.0)
    az = factor * position[2] * (5.0 * ratio - 3.0)
    return np.array([ax, ay, az])


def orbit_dynamics(state):
    pos, vel = state[0], state[1]
    acc = two_body_acceleration(pos) + j2_acceleration(pos)
    return (vel, acc)


def orbital_energy(state):
    pos, vel = state[0], state[1]
    r = np.linalg.norm(pos)
    v = np.linalg.norm(vel)
    return 0.5 * v**2 - EARTH_MU / r


def test_harmonic_oscillator():
    print("测试 1: 简谐振子 RK4 积分精度")
    print("-" * 50)

    omega = 1.0
    dt = 0.01
    steps = 1000

    def oscillator_deriv(state):
        x, v = state[0][0], state[1][0]
        return (np.array([v]), np.array([-omega**2 * x]))

    initial_x, initial_v = 1.0, 0.0
    state = (np.array([initial_x]), np.array([initial_v]), 0.0)

    for i in range(steps):
        state = rk4_step(state, oscillator_deriv, dt)

    expected_x = math.cos(omega * state[2])
    expected_v = -omega * math.sin(omega * state[2])

    error_x = abs(state[0][0] - expected_x)
    error_v = abs(state[1][0] - expected_v)

    print(f"  积分时长: {state[2]:.1f}s (周期的 {state[2]/(2*PI):.1f} 倍)")
    print(f"  位置误差: {error_x:.2e}")
    print(f"  速度误差: {error_v:.2e}")

    energy_initial = 0.5 * initial_v**2 + 0.5 * omega**2 * initial_x**2
    energy_final = 0.5 * state[1][0]**2 + 0.5 * omega**2 * state[0][0]**2
    energy_error = abs(energy_final - energy_initial) / energy_initial
    print(f"  能量相对误差: {energy_error:.2e}")

    if error_x < 1e-8 and energy_error < 1e-8:
        print("  ✓ RK4 积分器精度验证通过")
        return True
    else:
        print("  ✗ RK4 积分器精度验证失败")
        return False


def test_j2_acceleration():
    print("\n测试 2: J2 摄动加速度计算")
    print("-" * 50)

    pos_equator = np.array([EARTH_RADIUS + 500.0, 0.0, 0.0])
    acc_j2_eq = j2_acceleration(pos_equator)

    print(f"  赤道位置: r = {np.linalg.norm(pos_equator):.1f} km, z = 0")
    print(f"  J2 加速度: [{acc_j2_eq[0]:.6f}, {acc_j2_eq[1]:.6f}, {acc_j2_eq[2]:.6f}] km/s²")
    print(f"  赤道处 z 分量应为零: {abs(acc_j2_eq[2]) < 1e-10}")

    pos_polar = np.array([0.0, 0.0, EARTH_RADIUS + 500.0])
    acc_j2_pol = j2_acceleration(pos_polar)

    print(f"\n  极地点: r = {np.linalg.norm(pos_polar):.1f} km, z = r")
    print(f"  J2 加速度: [{acc_j2_pol[0]:.6f}, {acc_j2_pol[1]:.6f}, {acc_j2_pol[2]:.6f}] km/s²")
    print(f"  极点处 x,y 分量应为零: {abs(acc_j2_pol[0]) < 1e-10 and abs(acc_j2_pol[1]) < 1e-10}")

    if abs(acc_j2_eq[2]) < 1e-10 and abs(acc_j2_pol[0]) < 1e-10:
        print("  ✓ J2 摄动模型验证通过")
        return True
    else:
        print("  ✗ J2 摄动模型验证失败")
        return False


def test_orbit_energy_conservation():
    print("\n测试 3: 轨道能量守恒 (纯二体, 无 J2)")
    print("-" * 50)

    a = 7000.0
    e = 0.1
    v_circular = math.sqrt(EARTH_MU / a)

    pos = np.array([a * (1 - e), 0.0, 0.0])
    vel = np.array([0.0, v_circular * math.sqrt((1 + e) / (1 - e)), 0.0])
    state = (pos, vel, 0.0)

    def two_body_only(state):
        acc = two_body_acceleration(state[0])
        return (state[1], acc)

    dt = 10.0
    steps = int(2 * math.pi * math.sqrt(a**3 / EARTH_MU) / dt)

    energy0 = orbital_energy(state)

    for i in range(steps):
        state = rk4_step(state, two_body_only, dt)

    energy_final = orbital_energy(state)
    energy_error = abs(energy_final - energy0) / abs(energy0)

    print(f"  轨道半长轴: {a:.0f} km")
    print(f"  偏心率: {e}")
    print(f"  积分步数: {steps} (约 1 个轨道周期)")
    print(f"  初始能量: {energy0:.4f} km²/s²")
    print(f"  末态能量: {energy_final:.4f} km²/s²")
    print(f"  能量相对误差: {energy_error:.2e}")

    if energy_error < 1e-6:
        print("  ✓ 纯二体轨道能量守恒验证通过")
        return True
    else:
        print("  ✗ 纯二体轨道能量守恒验证失败")
        return False


def test_circular_orbit():
    print("\n测试 4: 圆轨道半径验证")
    print("-" * 50)

    r = EARTH_RADIUS + 550.0
    v = math.sqrt(EARTH_MU / r)

    pos = np.array([r, 0.0, 0.0])
    vel = np.array([0.0, v, 0.0])
    state = (pos, vel, 0.0)

    def two_body_only(state):
        acc = two_body_acceleration(state[0])
        return (state[1], acc)

    dt = 1.0
    period = 2 * PI * math.sqrt(r**3 / EARTH_MU)
    steps = int(period / dt)

    radii = []
    for i in range(steps):
        state = rk4_step(state, two_body_only, dt)
        radii.append(np.linalg.norm(state[0]))

    radii = np.array(radii)
    mean_r = radii.mean()
    std_r = radii.std()

    print(f"  标称半径: {r:.2f} km")
    print(f"  平均半径: {mean_r:.2f} km")
    print(f"  半径标准差: {std_r:.2f} km")
    print(f"  相对偏差: {abs(mean_r - r) / r:.2e}")

    if abs(mean_r - r) / r < 1e-6:
        print("  ✓ 圆轨道半径保持验证通过")
        return True
    else:
        print("  ✗ 圆轨道半径保持验证失败")
        return False


def main():
    print("=" * 60)
    print("LEO 轨道仿真引擎 - 核心算法验证")
    print("=" * 60)

    results = []
    results.append(test_harmonic_oscillator())
    results.append(test_j2_acceleration())
    results.append(test_orbit_energy_conservation())
    results.append(test_circular_orbit())

    print("\n" + "=" * 60)
    passed = sum(results)
    total = len(results)
    print(f"验证结果: {passed}/{total} 项通过")

    if passed == total:
        print("✓ 所有核心算法验证通过!")
        print("可以安全构建 C++ 扩展模块。")
    else:
        print("✗ 部分验证失败，请检查算法实现。")

    print("=" * 60)

    return 0 if passed == total else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
