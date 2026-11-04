#include "particle.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
void validate_dimensions(double length, double width)
{
    if (length <= 0.0 || width <= 0.0 || width > length)
    {
        throw std::invalid_argument("Particle dimensions require 0 < width <= length.");
    }
}
}

EllipseParticle::EllipseParticle() = default;

EllipseParticle::EllipseParticle(const vec& center, double theta,
    double length, double width)
    : center_(center), theta_(theta), length_(length), width_(width)
{
    validate_dimensions(length_, width_);
}

void EllipseParticle::translate(const vec& new_center) { center_ = new_center; }
void EllipseParticle::rotate(double new_theta) { theta_ = new_theta; }

void EllipseParticle::scale_size(double factor)
{
    if (factor <= 0.0) throw std::invalid_argument("Particle scale must be positive.");
    length_ *= factor;
    width_ *= factor;
}

const vec& EllipseParticle::get_center() const { return center_; }
double EllipseParticle::get_theta() const { return theta_; }
double EllipseParticle::get_length() const { return length_; }
double EllipseParticle::get_width() const { return width_; }
double EllipseParticle::get_area() const { return constants::pi * length_ * width_ / 4.0; }
double EllipseParticle::get_bounding_radius() const { return length_ / 2.0; }

SharpRodParticle::SharpRodParticle() = default;

SharpRodParticle::SharpRodParticle(const vec& center, double theta,
    double length, double width)
    : center_(center), theta_(theta), length_(length), width_(width)
{
    validate_dimensions(length_, width_);
}

void SharpRodParticle::translate(const vec& new_center) { center_ = new_center; }
void SharpRodParticle::rotate(double new_theta) { theta_ = new_theta; }

void SharpRodParticle::scale_size(double factor)
{
    if (factor <= 0.0) throw std::invalid_argument("Particle scale must be positive.");
    length_ *= factor;
    width_ *= factor;
}

const vec& SharpRodParticle::get_center() const { return center_; }
double SharpRodParticle::get_theta() const { return theta_; }
double SharpRodParticle::get_length() const { return length_; }
double SharpRodParticle::get_width() const { return width_; }

double SharpRodParticle::get_circle_radius() const
{
    const double a = length_ / 2.0;
    const double b = width_ / 2.0;
    return (a * a + b * b) / (2.0 * b);
}

double SharpRodParticle::get_circle_offset() const
{
    const double a = length_ / 2.0;
    const double b = width_ / 2.0;
    return (a * a - b * b) / (2.0 * b);
}

double SharpRodParticle::get_area() const
{
    const double a = length_ / 2.0;
    const double radius = get_circle_radius();
    const double offset = get_circle_offset();
    const double ratio = std::clamp(offset / radius, -1.0, 1.0);
    return 2.0 * radius * radius * std::acos(ratio) - 2.0 * a * offset;
}

double SharpRodParticle::get_bounding_radius() const { return length_ / 2.0; }
