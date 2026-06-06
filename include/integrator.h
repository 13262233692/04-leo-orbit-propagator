#pragma once
#include <vector>
#include <functional>
#include "vector3.h"
#include "gravity.h"

namespace leo_propagator {

using DerivativeFunc = std::function<State(const State&)>;

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
