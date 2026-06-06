#pragma once
#include <cmath>

namespace leo_propagator {

constexpr double PI = 3.14159265358979323846;
constexpr double DEG2RAD = PI / 180.0;
constexpr double RAD2DEG = 180.0 / PI;

constexpr double EARTH_RADIUS = 6378.137;
constexpr double EARTH_MU = 398600.4418;
constexpr double EARTH_J2 = 1.082635854e-3;
constexpr double EARTH_OMEGA = 7.292115146706979e-5;

constexpr double SECONDS_PER_DAY = 86400.0;
constexpr double MINUTES_PER_DAY = 1440.0;

constexpr double J2000_JD = 2451545.0;

}
