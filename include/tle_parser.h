#pragma once
#include <string>
#include "vector3.h"
#include "constants.h"

namespace leo_propagator {

struct TLE {
    std::string name;
    int norad_id;
    std::string line1;
    std::string line2;

    double epoch_year;
    double epoch_day;
    double inclination;
    double raan;
    double eccentricity;
    double arg_perigee;
    double mean_anomaly;
    double mean_motion;
    double bstar;

    double period;
    double semi_major_axis;
};

class TLEParser {
public:
    static TLE parse(const std::string& line1, const std::string& line2, const std::string& name = "");
    static State toInitialState(const TLE& tle);
    static double computeSemiMajorAxis(double mean_motion_rad_per_sec);
    static double getEpochJulianDate(const TLE& tle);

private:
    static double readAngle(const std::string& str, int start, int len);
    static double checkSum(const std::string& line);
};

}
