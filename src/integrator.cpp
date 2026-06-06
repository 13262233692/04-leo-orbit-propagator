#include "integrator.h"
#include <cmath>
#include <iostream>

namespace leo_propagator {

RK4Integrator::RK4Integrator(double step_size) : step_size_(step_size) {}

State RK4Integrator::step(const State& state, const DerivativeFunc& derivative) {
    double h = step_size_;

    State k1 = derivative(state);
    k1.position *= h;
    k1.velocity *= h;

    State state2 = state;
    state2.position += k1.position * 0.5;
    state2.velocity += k1.velocity * 0.5;
    State k2 = derivative(state2);
    k2.position *= h;
    k2.velocity *= h;

    State state3 = state;
    state3.position += k2.position * 0.5;
    state3.velocity += k2.velocity * 0.5;
    State k3 = derivative(state3);
    k3.position *= h;
    k3.velocity *= h;

    State state4 = state;
    state4.position += k3.position;
    state4.velocity += k3.velocity;
    State k4 = derivative(state4);
    k4.position *= h;
    k4.velocity *= h;

    Vector3 dp = (k1.position + k2.position * 2.0 + k3.position * 2.0 + k4.position) / 6.0;
    Vector3 dv = (k1.velocity + k2.velocity * 2.0 + k3.velocity * 2.0 + k4.velocity) / 6.0;

    State result = state;
    result.position += dp;
    result.velocity += dv;
    result.time = state.time + h;

    return result;
}

std::vector<State> RK4Integrator::integrate(const State& initial, double duration,
                                             const DerivativeFunc& derivative) {
    int num_steps = static_cast<int>(std::ceil(duration / step_size_)) + 1;
    std::vector<State> states;
    states.reserve(num_steps);

    State current = initial;
    states.push_back(current);

    KahanSum elapsed(0.0);
    KahanSum time_accum(initial.time);

    while (elapsed.value() < duration) {
        double remaining = duration - elapsed.value();
        double current_step = std::min(step_size_, remaining);
        if (current_step <= 0) break;

        double original_step = step_size_;
        step_size_ = current_step;
        current = step(current, derivative);
        step_size_ = original_step;

        time_accum.add(current_step);
        current.time = time_accum.value();

        elapsed.add(current_step);
        states.push_back(current);
    }

    return states;
}

HighPrecisionRK4::HighPrecisionRK4(double step_size) : step_size_(step_size) {}

State HighPrecisionRK4::step(const State& state, const DerivativeFunc& derivative) {
    double h = step_size_;

    State k1 = derivative(state);
    k1.position *= h;
    k1.velocity *= h;

    State state2;
    state2.position = state.position + k1.position * 0.5;
    state2.velocity = state.velocity + k1.velocity * 0.5;
    State k2 = derivative(state2);
    k2.position *= h;
    k2.velocity *= h;

    State state3;
    state3.position = state.position + k2.position * 0.5;
    state3.velocity = state.velocity + k2.velocity * 0.5;
    State k3 = derivative(state3);
    k3.position *= h;
    k3.velocity *= h;

    State state4;
    state4.position = state.position + k3.position;
    state4.velocity = state.velocity + k3.velocity;
    State k4 = derivative(state4);
    k4.position *= h;
    k4.velocity *= h;

    State result;
    result.position = state.position + (k1.position + k2.position * 2.0 + k3.position * 2.0 + k4.position) / 6.0;
    result.velocity = state.velocity + (k1.velocity + k2.velocity * 2.0 + k3.velocity * 2.0 + k4.velocity) / 6.0;
    result.time = state.time + h;

    return result;
}

std::vector<State> HighPrecisionRK4::integrate(const State& initial, double duration,
                                                const DerivativeFunc& derivative) {
    double estimated_steps = duration / step_size_ + 2;
    std::vector<State> states;
    states.reserve(static_cast<size_t>(estimated_steps));

    State current = initial;
    states.push_back(current);

    Vector3Kahan pos_accum;
    Vector3Kahan vel_accum;
    KahanSum time_accum(initial.time);
    KahanSum elapsed(0.0);

    pos_accum.reset(initial.position.x, initial.position.y, initial.position.z);
    vel_accum.reset(initial.velocity.x, initial.velocity.y, initial.velocity.z);

    State prev_state = initial;

    while (elapsed.value() < duration) {
        double remaining = duration - elapsed.value();
        double current_step = std::min(step_size_, remaining);
        if (current_step <= 1e-15) break;

        double original_step = step_size_;
        step_size_ = current_step;
        State next = step(prev_state, derivative);
        step_size_ = original_step;

        Vector3 dp = next.position - prev_state.position;
        Vector3 dv = next.velocity - prev_state.velocity;

        pos_accum.add(dp.x, dp.y, dp.z);
        vel_accum.add(dv.x, dv.y, dv.z);
        time_accum.add(current_step);
        elapsed.add(current_step);

        State corrected;
        corrected.position = Vector3(pos_accum.x(), pos_accum.y(), pos_accum.z());
        corrected.velocity = Vector3(vel_accum.x(), vel_accum.y(), vel_accum.z());
        corrected.time = time_accum.value();

        states.push_back(corrected);
        prev_state = next;
    }

    return states;
}

RK45Integrator::RK45Integrator(double rel_tol, double abs_tol, double max_step, double min_step)
    : rel_tol_(rel_tol), abs_tol_(abs_tol), max_step_(max_step), min_step_(min_step),
      last_step_(std::min(max_step, 1.0)), total_steps_(0), rejected_steps_(0) {}

Vector3 RK45Integrator::computeRKF45(const State& state, const DerivativeFunc& derivative,
                                      double dt, State& y5, State& y4) {
    const double c2 = 1.0/4.0;
    const double c3 = 3.0/8.0;
    const double c4 = 12.0/13.0;
    const double c5 = 1.0;
    const double c6 = 1.0/2.0;

    const double a21 = 1.0/4.0;
    const double a31 = 3.0/32.0;
    const double a32 = 9.0/32.0;
    const double a41 = 1932.0/2197.0;
    const double a42 = -7200.0/2197.0;
    const double a43 = 7296.0/2197.0;
    const double a51 = 439.0/216.0;
    const double a52 = -8.0;
    const double a53 = 3680.0/513.0;
    const double a54 = -845.0/4104.0;
    const double a61 = -8.0/27.0;
    const double a62 = 2.0;
    const double a63 = -3544.0/2565.0;
    const double a64 = 1859.0/4104.0;
    const double a65 = -11.0/40.0;

    const double b1 = 16.0/135.0;
    const double b2 = 0.0;
    const double b3 = 6656.0/12825.0;
    const double b4 = 28561.0/56430.0;
    const double b5 = -9.0/50.0;
    const double b6 = 2.0/55.0;

    const double b1s = 25.0/216.0;
    const double b2s = 0.0;
    const double b3s = 1408.0/2565.0;
    const double b4s = 2197.0/4104.0;
    const double b5s = -1.0/5.0;
    const double b6s = 0.0;

    State k1 = derivative(state);

    State s2;
    s2.position = state.position + k1.position * (a21 * dt);
    s2.velocity = state.velocity + k1.velocity * (a21 * dt);
    State k2 = derivative(s2);

    State s3;
    s3.position = state.position + (k1.position * a31 + k2.position * a32) * dt;
    s3.velocity = state.velocity + (k1.velocity * a31 + k2.velocity * a32) * dt;
    State k3 = derivative(s3);

    State s4;
    s4.position = state.position + (k1.position * a41 + k2.position * a42 + k3.position * a43) * dt;
    s4.velocity = state.velocity + (k1.velocity * a41 + k2.velocity * a42 + k3.velocity * a43) * dt;
    State k4 = derivative(s4);

    State s5;
    s5.position = state.position + (k1.position * a51 + k2.position * a52 + k3.position * a53 + k4.position * a54) * dt;
    s5.velocity = state.velocity + (k1.velocity * a51 + k2.velocity * a52 + k3.velocity * a53 + k4.velocity * a54) * dt;
    State k5 = derivative(s5);

    State s6;
    s6.position = state.position + (k1.position * a61 + k2.position * a62 + k3.position * a63 + k4.position * a64 + k5.position * a65) * dt;
    s6.velocity = state.velocity + (k1.velocity * a61 + k2.velocity * a62 + k3.velocity * a63 + k4.velocity * a64 + k5.velocity * a65) * dt;
    State k6 = derivative(s6);

    y5.position = state.position + (k1.position * b1 + k2.position * b2 + k3.position * b3 + k4.position * b4 + k5.position * b5 + k6.position * b6) * dt;
    y5.velocity = state.velocity + (k1.velocity * b1 + k2.velocity * b2 + k3.velocity * b3 + k4.velocity * b4 + k5.velocity * b5 + k6.velocity * b6) * dt;
    y5.time = state.time + dt;

    y4.position = state.position + (k1.position * b1s + k2.position * b2s + k3.position * b3s + k4.position * b4s + k5.position * b5s + k6.position * b6s) * dt;
    y4.velocity = state.velocity + (k1.velocity * b1s + k2.velocity * b2s + k3.velocity * b3s + k4.velocity * b4s + k5.velocity * b5s + k6.velocity * b6s) * dt;

    Vector3 pos_err = y5.position - y4.position;
    Vector3 vel_err = y5.velocity - y4.velocity;

    return Vector3(
        pos_err.norm() / dt,
        vel_err.norm() / dt,
        0.0
    );
}

State RK45Integrator::step(const State& state, const DerivativeFunc& derivative, double& dt_next) {
    double dt = last_step_;
    const double safety = 0.9;
    const double min_factor = 0.2;
    const double max_factor = 5.0;

    while (true) {
        State y5, y4;
        Vector3 err_vec = computeRKF45(state, derivative, dt, y5, y4);

        Vector3 scale_pos(
            abs_tol_ + rel_tol_ * std::max(std::abs(state.position.x), std::abs(y5.position.x)),
            abs_tol_ + rel_tol_ * std::max(std::abs(state.position.y), std::abs(y5.position.y)),
            abs_tol_ + rel_tol_ * std::max(std::abs(state.position.z), std::abs(y5.position.z))
        );

        Vector3 scale_vel(
            abs_tol_ + rel_tol_ * std::max(std::abs(state.velocity.x), std::abs(y5.velocity.x)),
            abs_tol_ + rel_tol_ * std::max(std::abs(state.velocity.y), std::abs(y5.velocity.y)),
            abs_tol_ + rel_tol_ * std::max(std::abs(state.velocity.z), std::abs(y5.velocity.z))
        );

        Vector3 pos_err = y5.position - y4.position;
        Vector3 vel_err = y5.velocity - y4.velocity;

        double err_pos = std::sqrt(
            (pos_err.x / scale_pos.x) * (pos_err.x / scale_pos.x) +
            (pos_err.y / scale_pos.y) * (pos_err.y / scale_pos.y) +
            (pos_err.z / scale_pos.z) * (pos_err.z / scale_pos.z)
        ) / 3.0;

        double err_vel = std::sqrt(
            (vel_err.x / scale_vel.x) * (vel_err.x / scale_vel.x) +
            (vel_err.y / scale_vel.y) * (vel_err.y / scale_vel.y) +
            (vel_err.z / scale_vel.z) * (vel_err.z / scale_vel.z)
        ) / 3.0;

        double err = std::max(err_pos, err_vel);

        total_steps_++;

        if (err <= 1.0) {
            double factor = (err == 0.0) ? max_factor : safety * std::pow(1.0 / err, 0.2);
            factor = std::min(std::max(factor, min_factor), max_factor);

            dt_next = std::min(dt * factor, max_step_);
            dt_next = std::max(dt_next, min_step_);
            last_step_ = dt;

            return y5;
        } else {
            rejected_steps_++;
            double factor = safety * std::pow(1.0 / err, 0.2);
            factor = std::min(std::max(factor, min_factor), max_factor);
            dt = std::max(dt * factor, min_step_);

            if (dt <= min_step_) {
                dt_next = min_step_;
                last_step_ = min_step_;
                return y5;
            }
        }
    }
}

std::vector<State> RK45Integrator::integrate(const State& initial, double duration,
                                              const DerivativeFunc& derivative) {
    std::vector<State> states;
    states.reserve(static_cast<size_t>(duration / std::min(last_step_, 10.0) + 100));

    State current = initial;
    states.push_back(current);

    KahanSum elapsed(0.0);
    KahanSum time_accum(initial.time);
    Vector3Kahan pos_accum;
    Vector3Kahan vel_accum;

    pos_accum.reset(initial.position.x, initial.position.y, initial.position.z);
    vel_accum.reset(initial.velocity.x, initial.velocity.y, initial.velocity.z);

    double dt_next = last_step_;
    State prev_raw = initial;

    while (elapsed.value() < duration - 1e-12) {
        double remaining = duration - elapsed.value();
        double dt = std::min(dt_next, remaining);
        dt = std::max(dt, min_step_);

        State next_raw = step(prev_raw, derivative, dt_next);
        double actual_dt = next_raw.time - prev_raw.time;

        if (actual_dt <= 0) {
            actual_dt = dt;
        }

        Vector3 dp = next_raw.position - prev_raw.position;
        Vector3 dv = next_raw.velocity - prev_raw.velocity;

        pos_accum.add(dp.x, dp.y, dp.z);
        vel_accum.add(dv.x, dv.y, dv.z);
        time_accum.add(actual_dt);
        elapsed.add(actual_dt);

        State corrected;
        corrected.position = Vector3(pos_accum.x(), pos_accum.y(), pos_accum.z());
        corrected.velocity = Vector3(vel_accum.x(), vel_accum.y(), vel_accum.z());
        corrected.time = time_accum.value();

        states.push_back(corrected);
        prev_raw = next_raw;
    }

    return states;
}

}
