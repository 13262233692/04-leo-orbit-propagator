#pragma once
#include <cmath>
#include <limits>

namespace leo_propagator {

class KahanSum {
public:
    KahanSum() : sum_(0.0), compensation_(0.0) {}
    explicit KahanSum(double initial) : sum_(initial), compensation_(0.0) {}

    void add(double value) {
        double y = value - compensation_;
        double t = sum_ + y;
        compensation_ = (t - sum_) - y;
        sum_ = t;
    }

    double value() const { return sum_; }

    void reset() { sum_ = 0.0; compensation_ = 0.0; }
    void reset(double value) { sum_ = value; compensation_ = 0.0; }

    KahanSum& operator+=(double value) {
        add(value);
        return *this;
    }

    operator double() const { return sum_; }

private:
    double sum_;
    double compensation_;
};

class Vector3Kahan {
public:
    Vector3Kahan() : x_(0.0), y_(0.0), z_(0.0) {}

    void add(double dx, double dy, double dz) {
        x_.add(dx);
        y_.add(dy);
        z_.add(dz);
    }

    void reset() { x_.reset(); y_.reset(); z_.reset(); }
    void reset(double x, double y, double z) { x_.reset(x); y_.reset(y); z_.reset(z); }

    double x() const { return x_.value(); }
    double y() const { return y_.value(); }
    double z() const { return z_.value(); }

private:
    KahanSum x_, y_, z_;
};

constexpr double EPSILON_DOUBLE = 2.2204460492503131e-16;
constexpr double MACHINE_EPSILON = std::numeric_limits<double>::epsilon();

}
