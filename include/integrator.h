#pragma once
#include <vector>
#include <functional>
#include <cmath>
#include "vector3.h"
#include "gravity.h"
#include "numerics.h"

namespace leo_propagator {

using DerivativeFunc = std::function<State(const State&)>;

class HighPrecisionRK4 {
public:
    HighPrecisionRK4(double step_size = 1.0);

    State step(const State& state, const DerivativeFunc& derivative);
    std::vector<State> integrate(const State& initial, double duration,
                                 const DerivativeFunc& derivative);

    void setStepSize(double step) { step_size_ = step; }
    double getStepSize() const { return step_size_; }

private:
    double step_size_;
};

class RK45Integrator {
public:
    RK45Integrator(
        double rel_tol = 1e-12,
        double abs_tol = 1e-15,
        double max_step = 60.0,
        double min_step = 1e-9
    );

    State step(const State& state, const DerivativeFunc& derivative, double& dt_next);
    std::vector<State> integrate(const State& initial, double duration,
                                 const DerivativeFunc& derivative);

    void setRelativeTolerance(double tol) { rel_tol_ = tol; }
    void setAbsoluteTolerance(double tol) { abs_tol_ = tol; }
    void setMaxStep(double max) { max_step_ = max; }
    void setMinStep(double min) { min_step_ = min; }
    double getLastStepSize() const { return last_step_; }
    long getTotalSteps() const { return total_steps_; }
    long getRejectedSteps() const { return rejected_steps_; }

private:
    Vector3 computeRKF45(const State& state, const DerivativeFunc& derivative,
                         double dt, State& y5, State& y4);

    double rel_tol_;
    double abs_tol_;
    double max_step_;
    double min_step_;
    double last_step_;
    long total_steps_;
    long rejected_steps_;
};

class RK4Integrator {
public:
    RK4Integrator(double step_size = 1.0);

    State step(const State& state, const DerivativeFunc& derivative);
    std::vector<State> integrate(const State& initial, double duration,
                                 const DerivativeFunc& derivative);

    void setStepSize(double step) { step_size_ = step; }
    double getStepSize() const { return step_size_; }

private:
    double step_size_;
};

class OrbitDynamics {
public:
    OrbitDynamics(bool use_j2 = true) : use_j2_(use_j2) {}

    State operator()(const State& state) const {
        Vector3 acc = GravityModel::computeAcceleration(state.position, use_j2_);
        return State(state.velocity, acc, 1.0);
    }

private:
    bool use_j2_;
};

}
