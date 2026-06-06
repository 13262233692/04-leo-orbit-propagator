#include "propagator.h"
#include "constants.h"
#include <thread>
#include <future>
#include <vector>

namespace leo_propagator {

OrbitPropagator::OrbitPropagator(double step_size, bool use_j2)
    : integrator_(step_size), dynamics_(use_j2), use_j2_(use_j2) {}

OrbitResult OrbitPropagator::propagateTLE(const TLE& tle, double duration_days) {
    State initial = TLEParser::toInitialState(tle);
    double duration_sec = duration_days * SECONDS_PER_DAY;
    return propagateState(initial, duration_sec);
}

OrbitResult OrbitPropagator::propagateState(const State& initial_state, double duration_seconds) {
    std::vector<State> states = integrator_.integrate(initial_state, duration_seconds, dynamics_);

    OrbitResult result;
    result.positions.reserve(states.size());
    result.velocities.reserve(states.size());
    result.times.reserve(states.size());

    for (const auto& s : states) {
        result.positions.push_back(s.position);
        result.velocities.push_back(s.velocity);
        result.times.push_back(s.time);
    }

    return result;
}

BatchPropagator::BatchPropagator(double step_size, bool use_j2)
    : step_size_(step_size), use_j2_(use_j2) {}

std::vector<OrbitResult> BatchPropagator::propagateBatch(const std::vector<TLE>& tles,
                                                          double duration_days) {
    std::vector<OrbitResult> results(tles.size());

    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4;

    size_t chunk_size = (tles.size() + num_threads - 1) / num_threads;

    auto worker = [&](size_t start, size_t end) {
        OrbitPropagator propagator(step_size_, use_j2_);
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
