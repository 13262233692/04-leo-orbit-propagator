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

namespace py = pybind11;
using namespace leo_propagator;

PYBIND11_MODULE(_leo_propagator, m) {
    m.doc() = "LEO Orbit Propagator - High precision satellite orbit simulation engine";

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

    py::class_<OrbitResult>(m, "OrbitResult")
        .def(py::init<>())
        .def_readwrite("positions", &OrbitResult::positions)
        .def_readwrite("velocities", &OrbitResult::velocities)
        .def_readwrite("times", &OrbitResult::times)
        .def_readwrite("norad_id", &OrbitResult::norad_id)
        .def_readwrite("satellite_name", &OrbitResult::satellite_name)
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
        .def(py::init<double, bool>(), 
            py::arg("step_size") = 1.0, py::arg("use_j2") = true)
        .def("propagate_tle", &OrbitPropagator::propagateTLE,
            py::arg("tle"), py::arg("duration_days") = 7.0)
        .def("propagate_state", &OrbitPropagator::propagateState,
            py::arg("initial_state"), py::arg("duration_seconds"))
        .def("set_step_size", &OrbitPropagator::setStepSize)
        .def("set_use_j2", &OrbitPropagator::setUseJ2);

    py::class_<BatchPropagator>(m, "BatchPropagator")
        .def(py::init<double, bool>(),
            py::arg("step_size") = 1.0, py::arg("use_j2") = true)
        .def("propagate_batch", &BatchPropagator::propagateBatch,
            py::arg("tles"), py::arg("duration_days") = 7.0)
        .def("set_step_size", &BatchPropagator::setStepSize)
        .def("set_use_j2", &BatchPropagator::setUseJ2);
}
