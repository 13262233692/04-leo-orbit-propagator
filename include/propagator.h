#pragma once
#include <vector>
#include <string>
#include <memory>
#include "vector3.h"
#include "tle_parser.h"
#include "integrator.h"

namespace leo_propagator {

enum class IntegratorType {
    RK4,
    HIGH_PRECISION_RK4,
    RK45_ADAPTIVE
};

struct OrbitResult {
    std::vector<Vector3> positions;
    std::vector<Vector3> velocities;
    std::vector<double> times;
    int norad_id;
    std::string satellite_name;
    IntegratorType integrator_used;
    long total_steps;
    long rejected_steps;
    double compute_time;
};

class OrbitPropagator {
public:
    OrbitPropagator(
        double step_size = 1.0,
        bool use_j2 = true,
        IntegratorType integrator_type = IntegratorType::HIGH_PRECISION_RK4
    );

    OrbitResult propagateTLE(const TLE& tle, double duration_days = 7.0);
    OrbitResult propagateState(const State& initial_state, double duration_seconds);

    void setStepSize(double step);
    void setUseJ2(bool use) { use_j2_ = use; dynamics_ = OrbitDynamics(use); }
    void setIntegratorType(IntegratorType type);

    IntegratorType getIntegratorType() const { return integrator_type_; }

    void setRKF45Tolerances(double rel_tol, double abs_tol);
    void setRKF45StepLimits(double max_step, double min_step);

private:
    IntegratorType integrator_type_;
    double step_size_;
    bool use_j2_;
    OrbitDynamics dynamics_;

    std::unique_ptr<RK4Integrator> rk4_;
    std::unique_ptr<HighPrecisionRK4> hp_rk4_;
    std::unique_ptr<RK45Integrator> rk45_;

    void ensureIntegrator();
};

class BatchPropagator {
public:
    BatchPropagator(
        double step_size = 1.0,
        bool use_j2 = true,
        IntegratorType integrator_type = IntegratorType::HIGH_PRECISION_RK4
    );

    std::vector<OrbitResult> propagateBatch(const std::vector<TLE>& tles,
                                             double duration_days = 7.0);

    void setStepSize(double step) { step_size_ = step; }
    void setUseJ2(bool use) { use_j2_ = use; }
    void setIntegratorType(IntegratorType type) { integrator_type_ = type; }

private:
    double step_size_;
    bool use_j2_;
    IntegratorType integrator_type_;
};

}
