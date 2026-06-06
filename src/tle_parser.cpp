#include "tle_parser.h"
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <sstream>

namespace leo_propagator {

TLE TLEParser::parse(const std::string& line1, const std::string& line2, const std::string& name) {
    TLE tle;
    tle.name = name;
    tle.line1 = line1;
    tle.line2 = line2;

    if (line1.length() < 69 || line2.length() < 69) {
        throw std::invalid_argument("TLE lines must be at least 69 characters");
    }

    if (line1[0] != '1' || line2[0] != '2') {
        throw std::invalid_argument("Invalid TLE format");
    }

    tle.norad_id = std::atoi(line1.substr(2, 5).c_str());
    tle.epoch_year = std::atoi(line1.substr(18, 2).c_str());
    if (tle.epoch_year < 57) {
        tle.epoch_year += 2000;
    } else {
        tle.epoch_year += 1900;
    }
    tle.epoch_day = std::atof(line1.substr(20, 12).c_str());
    tle.bstar = std::atof(line1.substr(53, 8).c_str());

    tle.inclination = readAngle(line2, 8, 8);
    tle.raan = readAngle(line2, 17, 8);
    tle.eccentricity = std::atof(("0." + line2.substr(26, 7)).c_str());
    tle.arg_perigee = readAngle(line2, 34, 8);
    tle.mean_anomaly = readAngle(line2, 43, 8);
    tle.mean_motion = std::atof(line2.substr(52, 11).c_str());

    tle.mean_motion = tle.mean_motion * 2.0 * PI / MINUTES_PER_DAY;
    tle.period = 2.0 * PI / tle.mean_motion;
    tle.semi_major_axis = computeSemiMajorAxis(tle.mean_motion);

    tle.inclination *= DEG2RAD;
    tle.raan *= DEG2RAD;
    tle.arg_perigee *= DEG2RAD;
    tle.mean_anomaly *= DEG2RAD;

    return tle;
}

double TLEParser::readAngle(const std::string& str, int start, int len) {
    return std::atof(str.substr(start, len).c_str());
}

double TLEParser::computeSemiMajorAxis(double mean_motion_rad_per_sec) {
    double mu = EARTH_MU;
    double n = mean_motion_rad_per_sec;
    return std::pow(mu / (n * n), 1.0 / 3.0);
}

double TLEParser::getEpochJulianDate(const TLE& tle) {
    int year = static_cast<int>(tle.epoch_year);
    double day_of_year = tle.epoch_day;

    int a = (year - 1) / 100;
    int b = 2 - a + a / 4;
    double jd = static_cast<int>(365.25 * (year)) + 
                static_cast<int>(30.6001 * 1) + 
                day_of_year + 1720994.5 + b;
    
    return jd;
}

State TLEParser::toInitialState(const TLE& tle) {
    double a = tle.semi_major_axis;
    double e = tle.eccentricity;
    double i = tle.inclination;
    double raan = tle.raan;
    double argp = tle.arg_perigee;
    double M0 = tle.mean_anomaly;
    double mu = EARTH_MU;

    double E = M0;
    for (int iter = 0; iter < 50; ++iter) {
        double E_new = E - (E - e * std::sin(E) - M0) / (1.0 - e * std::cos(E));
        if (std::abs(E_new - E) < 1e-12) {
            E = E_new;
            break;
        }
        E = E_new;
    }

    double nu = 2.0 * std::atan2(
        std::sqrt(1.0 + e) * std::sin(E / 2.0),
        std::sqrt(1.0 - e) * std::cos(E / 2.0)
    );

    double r_norm = a * (1.0 - e * std::cos(E));

    Vector3 r_pqw(
        r_norm * std::cos(nu),
        r_norm * std::sin(nu),
        0.0
    );

    double p = a * (1.0 - e * e);
    double v_mag = std::sqrt(mu / p);

    Vector3 v_pqw(
        -v_mag * std::sin(nu),
        v_mag * (e + std::cos(nu)),
        0.0
    );

    auto rot_z = [](double angle) -> std::array<std::array<double, 3>, 3> {
        double c = std::cos(angle);
        double s = std::sin(angle);
        return {{
            {c, -s, 0},
            {s, c, 0},
            {0, 0, 1}
        }};
    };

    auto rot_x = [](double angle) -> std::array<std::array<double, 3>, 3> {
        double c = std::cos(angle);
        double s = std::sin(angle);
        return {{
            {1, 0, 0},
            {0, c, -s},
            {0, s, c}
        }};
    };

    auto matmul = [](const std::array<std::array<double, 3>, 3>& m, const Vector3& v) -> Vector3 {
        return Vector3(
            m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z,
            m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z,
            m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z
        );
    };

    auto matmul3 = [&](const std::array<std::array<double, 3>, 3>& a,
                       const std::array<std::array<double, 3>, 3>& b,
                       const std::array<std::array<double, 3>, 3>& c,
                       const Vector3& v) -> Vector3 {
        Vector3 res = v;
        res = matmul(c, res);
        res = matmul(b, res);
        res = matmul(a, res);
        return res;
    };

    Vector3 position = matmul3(
        rot_z(raan),
        rot_x(i),
        rot_z(argp),
        r_pqw
    );

    Vector3 velocity = matmul3(
        rot_z(raan),
        rot_x(i),
        rot_z(argp),
        v_pqw
    );

    return State(position, velocity, 0.0);
}

double TLEParser::checkSum(const std::string& line) {
    int sum = 0;
    for (char c : line) {
        if (c >= '0' && c <= '9') {
            sum += c - '0';
        } else if (c == '-') {
            sum += 1;
        }
    }
    return sum % 10;
}

}
