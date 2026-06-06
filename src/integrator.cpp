#include "integrator.h"
#include <cmath>

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

    State result = state;
    result.position += (k1.position + k2.position * 2.0 + k3.position * 2.0 + k4.position) / 6.0;
    result.velocity += (k1.velocity + k2.velocity * 2.0 + k3.velocity * 2.0 + k4.velocity) / 6.0;
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

    double elapsed = 0.0;
    while (elapsed < duration) {
        double current_step = std::min(step_size_, duration - elapsed);
        if (current_step <= 0) break;

        double original_step = step_size_;
        step_size_ = current_step;
        current = step(current, derivative);
        step_size_ = original_step;

        states.push_back(current);
        elapsed += current_step;
    }

    return states;
}

}
