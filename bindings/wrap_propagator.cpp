#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

#include "vector3.h"
#include "tle_parser.h"
#include "gravity.h"
#include "integrator.h"
#include "propagator.h"
#include "numerics.h"

namespace py = pybind11;
using namespace leo_propagator;

PYBIND11_MODULE(_leo_propagator, m) {
    m.doc() = "LEO Orbit Propagator - High precision satellite orbit simulation engine "
              "with Kahan summation and RK45 adaptive step integrators";

    py::enum_<IntegratorType>(m, "IntegratorType")
        .value("RK4", IntegratorType::RK4)
        .value("HIGH_PRECISION_RK4", IntegratorType::HIGH_PRECISION_RK4)
        .value("RK45_ADAPTIVE", IntegratorType::RK45_ADAPTIVE)
        .export_values();

    py::class_<KahanSum>(m, "KahanSum")
        .def(py::init<>())
        .def(py::init<double>())
        .def("add", &KahanSum::add)
        .def("value", &KahanSum::value)
        .def("reset", py::overload_cast<>(&KahanSum::reset))
        .def("reset", py::overload_cast<double>(&KahanSum::reset))
        .def("__iadd__", &KahanSum::operator+=)
        .def("__float__", &KahanSum::operator double);

    py::class_<Vector3>(m, "Vector3")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def_readwrite("x", &Vector3::x)
        .def_readwrite("y", &Vector3::y)
        .def_readwrite("z", &Vector3::z)
        .def("norm", &Vector3::norm)
        .def("norm_squared", &Vector3::normSquared)
        .def("normalized", &Vector3::normalized)
        .def("dot", &Vector3::dot)
        .def("cross", &Vector3::cross)
        .def("to_array", [](const Vector3& v) {
            return py::array_t<double>(3, &v.x);
        })
        .def("__add__", &Vector3::operator+)
        .def("__sub__", &Vector3::operator-)
        .def("__mul__", &Vector3::operator*)
        .def("__repr__", [](const Vector3& v) {
            return "Vector3(" + std::to_string(v.x) + ", " +
                   std::to_string(v.y) + ", " + std::to_string(v.z) + ")";
        });

    py::class_<State>(m, "State")
        .def(py::init<>())
        .def(py::init<Vector3, Vector3, double>())
        .def_readwrite("position", &State::position)
        .def_readwrite("velocity", &State::velocity)
        .def_readwrite("time", &State::time);

    py::class_<TLE>(m, "TLE")
        .def(py::init<>())
        .def_readwrite("name", &TLE::name)
        .def_readwrite("norad_id", &TLE::norad_id)
        .def_readwrite("line1", &TLE::line1)
        .def_readwrite("line2", &TLE::line2)
        .def_readwrite("epoch_year", &TLE::epoch_year)
        .def_readwrite("epoch_day", &TLE::epoch_day)
        .def_readwrite("inclination", &TLE::inclination)
        .def_readwrite("raan", &TLE::raan)
        .def_readwrite("eccentricity", &TLE::eccentricity)
        .def_readwrite("arg_perigee", &TLE::arg_perigee)
        .def_readwrite("mean_anomaly", &TLE::mean_anomaly)
        .def_readwrite("mean_motion", &TLE::mean_motion)
        .def_readwrite("bstar", &TLE::bstar)
        .def_readwrite("period", &TLE::period)
        .def_readwrite("semi_major_axis", &TLE::semi_major_axis);

    py::class_<TLEParser>(m, "TLEParser")
        .def_static("parse", &TLEParser::parse,
            py::arg("line1"), py::arg("line2"), py::arg("name") = "")
        .def_static("to_initial_state", &TLEParser::toInitialState)
        .def_static("compute_semi_major_axis", &TLEParser::computeSemiMajorAxis)
        .def_static("get_epoch_julian_date", &TLEParser::getEpochJulianDate);

    py::class_<GravityModel>(m, "GravityModel")
        .def_static("compute_acceleration", &GravityModel::computeAcceleration,
            py::arg("position"), py::arg("include_j2") = true)
        .def_static("compute_two_body_acceleration", &GravityModel::computeTwoBodyAcceleration)
        .def_static("compute_j2_acceleration", &GravityModel::computeJ2Acceleration);

    py::class_<RK4Integrator>(m, "RK4Integrator")
        .def(py::init<double>(), py::arg("step_size") = 1.0)
        .def("step", &RK4Integrator::step)
        .def("integrate", &RK4Integrator::integrate)
        .def_property("step_size", &RK4Integrator::getStepSize, &RK4Integrator::setStepSize);

    py::class_<HighPrecisionRK4>(m, "HighPrecisionRK4")
        .def(py::init<double>(), py::arg("step_size") = 1.0)
        .def("step", &HighPrecisionRK4::step)
        .def("integrate", &HighPrecisionRK4::integrate)
        .def_property("step_size", &HighPrecisionRK4::getStepSize, &HighPrecisionRK4::setStepSize)
        .doc() = "RK4 integrator with Kahan compensated summation for high precision";

    py::class_<RK45Integrator>(m, "RK45Integrator")
        .def(py::init<double, double, double, double>(),
            py::arg("rel_tol") = 1e-12,
            py::arg("abs_tol") = 1e-15,
            py::arg("max_step") = 60.0,
            py::arg("min_step") = 1e-9)
        .def("integrate", &RK45Integrator::integrate)
        .def("set_relative_tolerance", &RK45Integrator::setRelativeTolerance)
        .def("set_absolute_tolerance", &RK45Integrator::setAbsoluteTolerance)
        .def("set_max_step", &RK45Integrator::setMaxStep)
        .def("set_min_step", &RK45Integrator::setMinStep)
        .def("get_last_step_size", &RK45Integrator::getLastStepSize)
        .def("get_total_steps", &RK45Integrator::getTotalSteps)
        .def("get_rejected_steps", &RK45Integrator::getRejectedSteps)
        .doc() = "Dormand-Prince RK45 adaptive step size integrator with error control";

    py::class_<OrbitResult>(m, "OrbitResult")
        .def(py::init<>())
        .def_readwrite("positions", &OrbitResult::positions)
        .def_readwrite("velocities", &OrbitResult::velocities)
        .def_readwrite("times", &OrbitResult::times)
        .def_readwrite("norad_id", &OrbitResult::norad_id)
        .def_readwrite("satellite_name", &OrbitResult::satellite_name)
        .def_readwrite("integrator_used", &OrbitResult::integrator_used)
        .def_readwrite("total_steps", &OrbitResult::total_steps)
        .def_readwrite("rejected_steps", &OrbitResult::rejected_steps)
        .def_readwrite("compute_time", &OrbitResult::compute_time)
        .def("get_positions_array", [](const OrbitResult& self) {
            size_t n = self.positions.size();
            py::array_t<double> arr({n, 3});
            auto ptr = arr.mutable_unchecked<2>();
            for (size_t i = 0; i < n; ++i) {
                ptr(i, 0) = self.positions[i].x;
                ptr(i, 1) = self.positions[i].y;
                ptr(i, 2) = self.positions[i].z;
            }
            return arr;
        })
        .def("get_velocities_array", [](const OrbitResult& self) {
            size_t n = self.velocities.size();
            py::array_t<double> arr({n, 3});
            auto ptr = arr.mutable_unchecked<2>();
            for (size_t i = 0; i < n; ++i) {
                ptr(i, 0) = self.velocities[i].x;
                ptr(i, 1) = self.velocities[i].y;
                ptr(i, 2) = self.velocities[i].z;
            }
            return arr;
        })
        .def("get_times_array", [](const OrbitResult& self) {
            size_t n = self.times.size();
            py::array_t<double> arr(n);
            auto ptr = arr.mutable_unchecked<1>();
            for (size_t i = 0; i < n; ++i) {
                ptr(i) = self.times[i];
            }
            return arr;
        });

    py::class_<OrbitPropagator>(m, "OrbitPropagator")
        .def(py::init<double, bool, IntegratorType>(),
            py::arg("step_size") = 1.0,
            py::arg("use_j2") = true,
            py::arg("integrator_type") = IntegratorType::HIGH_PRECISION_RK4)
        .def("propagate_tle", &OrbitPropagator::propagateTLE,
            py::arg("tle"), py::arg("duration_days") = 7.0)
        .def("propagate_state", &OrbitPropagator::propagateState,
            py::arg("initial_state"), py::arg("duration_seconds"))
        .def("set_step_size", &OrbitPropagator::setStepSize)
        .def("set_use_j2", &OrbitPropagator::setUseJ2)
        .def("set_integrator_type", &OrbitPropagator::setIntegratorType)
        .def("get_integrator_type", &OrbitPropagator::getIntegratorType)
        .def("set_rkf45_tolerances", &OrbitPropagator::setRKF45Tolerances,
            py::arg("rel_tol"), py::arg("abs_tol"))
        .def("set_rkf45_step_limits", &OrbitPropagator::setRKF45StepLimits,
            py::arg("max_step"), py::arg("min_step"));

    py::class_<BatchPropagator>(m, "BatchPropagator")
        .def(py::init<double, bool, IntegratorType>(),
            py::arg("step_size") = 1.0,
            py::arg("use_j2") = true,
            py::arg("integrator_type") = IntegratorType::HIGH_PRECISION_RK4)
        .def("propagate_batch", &BatchPropagator::propagateBatch,
            py::arg("tles"), py::arg("duration_days") = 7.0)
        .def("set_step_size", &BatchPropagator::setStepSize)
        .def("set_use_j2", &BatchPropagator::setUseJ2)
        .def("set_integrator_type", &BatchPropagator::setIntegratorType);

    m.attr("__version__") = "2.0.0";
    m.attr("__precision_fix__") = "Kahan summation + RK45 adaptive step";
}
