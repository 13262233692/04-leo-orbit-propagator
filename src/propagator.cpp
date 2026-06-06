#include "propagator.h"
#include "constants.h"
#include <thread>
#include <future>
#include <vector>
#include <chrono>

namespace leo_propagator {

OrbitPropagator::OrbitPropagator(double step_size, bool use_j2, IntegratorType integrator_type)
    : integrator_type_(integrator_type), step_size_(step_size), use_j2_(use_j2),
      dynamics_(use_j2) {
    ensureIntegrator();
}

void OrbitPropagator::ensureIntegrator() {
    switch (integrator_type_) {
        case IntegratorType::RK4:
            if (!rk4_) rk4_ = std::make_unique<RK4Integrator>(step_size_);
            break;
        case IntegratorType::HIGH_PRECISION_RK4:
            if (!hp_rk4_) hp_rk4_ = std::make_unique<HighPrecisionRK4>(step_size_);
            break;
        case IntegratorType::RK45_ADAPTIVE:
            if (!rk45_) rk45_ = std::make_unique<RK45Integrator>();
            break;
    }
}

void OrbitPropagator::setStepSize(double step) {
    step_size_ = step;
    if (rk4_) rk4_->setStepSize(step);
    if (hp_rk4_) hp_rk4_->setStepSize(step);
}

void OrbitPropagator::setIntegratorType(IntegratorType type) {
    integrator_type_ = type;
    ensureIntegrator();
}

void OrbitPropagator::setRKF45Tolerances(double rel_tol, double abs_tol) {
    if (rk45_) {
        rk45_->setRelativeTolerance(rel_tol);
        rk45_->setAbsoluteTolerance(abs_tol);
    }
}

void OrbitPropagator::setRKF45StepLimits(double max_step, double min_step) {
    if (rk45_) {
        rk45_->setMaxStep(max_step);
        rk45_->setMinStep(min_step);
    }
}

OrbitResult OrbitPropagator::propagateTLE(const TLE& tle, double duration_days) {
    State initial = TLEParser::toInitialState(tle);
    double duration_sec = duration_days * SECONDS_PER_DAY;
    OrbitResult result = propagateState(initial, duration_sec);
    result.norad_id = tle.norad_id;
    result.satellite_name = tle.name;
    return result;
}

OrbitResult OrbitPropagator::propagateState(const State& initial_state, double duration_seconds) {
    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<State> states;

    ensureIntegrator();

    long total_steps = 0;
    long rejected_steps = 0;

    switch (integrator_type_) {
        case IntegratorType::RK4:
            states = rk4_->integrate(initial_state, duration_seconds, dynamics_);
            total_steps = states.size() - 1;
            break;
        case IntegratorType::HIGH_PRECISION_RK4:
            states = hp_rk4_->integrate(initial_state, duration_seconds, dynamics_);
            total_steps = states.size() - 1;
            break;
        case IntegratorType::RK45_ADAPTIVE:
            states = rk45_->integrate(initial_state, duration_seconds, dynamics_);
            total_steps = rk45_->getTotalSteps();
            rejected_steps = rk45_->getRejectedSteps();
            break;
    }

    OrbitResult result;
    result.positions.reserve(states.size());
    result.velocities.reserve(states.size());
    result.times.reserve(states.size());

    for (const auto& s : states) {
        result.positions.push_back(s.position);
        result.velocities.push_back(s.velocity);
        result.times.push_back(s.time);
    }

    result.integrator_used = integrator_type_;
    result.total_steps = total_steps;
    result.rejected_steps = rejected_steps;

    auto end_time = std::chrono::high_resolution_clock::now();
    result.compute_time = std::chrono::duration<double>(end_time - start_time).count();

    return result;
}

BatchPropagator::BatchPropagator(double step_size, bool use_j2, IntegratorType integrator_type)
    : step_size_(step_size), use_j2_(use_j2), integrator_type_(integrator_type) {}

std::vector<OrbitResult> BatchPropagator::propagateBatch(const std::vector<TLE>& tles,
                                                          double duration_days) {
    std::vector<OrbitResult> results(tles.size());

    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4;

    size_t chunk_size = (tles.size() + num_threads - 1) / num_threads;

    auto worker = [&](size_t start, size_t end) {
        OrbitPropagator propagator(step_size_, use_j2_, integrator_type_);
        for (size_t i = start; i < end && i < tles.size(); ++i) {
            OrbitResult res = propagator.propagateTLE(tles[i], duration_days);
            res.norad_id = tles[i].norad_id;
            res.satellite_name = tles[i].name;
            results[i] = std::move(res);
        }
    };

    std::vector<std::future<void>> futures;
    for (unsigned int t = 0; t < num_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = start + chunk_size;
        futures.push_back(std::async(std::launch::async, worker, start, end));
    }

    for (auto& f : futures) {
        f.get();
    }

    return results;
}

}
