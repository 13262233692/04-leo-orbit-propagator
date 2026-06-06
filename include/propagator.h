#pragma once
#include <vector>
#include <string>
#include "vector3.h"
#include "tle_parser.h"
#include "integrator.h"

namespace leo_propagator {

struct OrbitResult {
    std::vector<Vector3> positions;
    std::vector<Vector3> velocities;
    std::vector<double> times;
    int norad_id;
    std::string satellite_name;
};

class OrbitPropagator {
public:
    OrbitPropagator(double step_size = 1.0, bool use_j2 = true);

    OrbitResult propagateTLE(const TLE& tle, double duration_days = 7.0);
    OrbitResult propagateState(const State& initial_state, double duration_seconds);

    void setStepSize(double step) { integrator_.setStepSize(step); }
    void setUseJ2(bool use) { use_j2_ = use; dynamics_ = OrbitDynamics(use); }

private:
    RK4Integrator integrator_;
    OrbitDynamics dynamics_;
    bool use_j2_;
};

class BatchPropagator {
public:
    BatchPropagator(double step_size = 1.0, bool use_j2 = true);

    std::vector<OrbitResult> propagateBatch(const std::vector<TLE>& tles, 
                                             double duration_days = 7.0);

    void setStepSize(double step) { step_size_ = step; }
    void setUseJ2(bool use) { use_j2_ = use; }

private:
    double step_size_;
    bool use_j2_;
};

}
