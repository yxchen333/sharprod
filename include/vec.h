#pragma once

#include <cmath>

namespace constants
{
    constexpr double pi = 3.1415926535897932384626433832795;
}

class vec
{
public:
    vec() = default;
    vec(double set_x, double set_y) : x(set_x), y(set_y) {}

    vec operator+(const vec& other) const { return vec(x + other.x, y + other.y); }
    vec operator-(const vec& other) const { return vec(x - other.x, y - other.y); }
    vec operator*(double scale) const { return vec(x * scale, y * scale); }
    vec operator/(double scale) const { return vec(x / scale, y / scale); }

    vec& operator+=(const vec& other)
    {
        x += other.x;
        y += other.y;
        return *this;
    }

    double dot(const vec& other) const { return x * other.x + y * other.y; }
    double cross(const vec& other) const { return x * other.y - y * other.x; }
    double mod2() const { return x * x + y * y; }
    double mod() const { return std::sqrt(mod2()); }

    vec normalized() const
    {
        const double length = mod();
        return (length > 1.0e-15) ? *this / length : vec();
    }

    vec rotated(double theta) const
    {
        const double c = std::cos(theta);
        const double s = std::sin(theta);
        return vec(c * x - s * y, s * x + c * y);
    }

    double x = 0.0;
    double y = 0.0;
};

inline vec operator*(double scale, const vec& value)
{
    return value * scale;
}
