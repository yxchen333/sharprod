#pragma once

#include "vec.h"

class EllipseParticle
{
public:
    EllipseParticle();
    EllipseParticle(const vec& center, double theta, double length, double width);

    void translate(const vec& new_center);
    void rotate(double new_theta);
    void scale_size(double factor);

    const vec& get_center() const;
    double get_theta() const;
    double get_length() const;
    double get_width() const;
    double get_area() const;
    double get_bounding_radius() const;

private:
    vec center_;
    double theta_ = 0.0;
    double length_ = 1.0;
    double width_ = 1.0;
};

class SharpRodParticle
{
public:
    SharpRodParticle();
    SharpRodParticle(const vec& center, double theta, double length, double width);

    void translate(const vec& new_center);
    void rotate(double new_theta);
    void scale_size(double factor);

    const vec& get_center() const;
    double get_theta() const;
    double get_length() const;
    double get_width() const;
    double get_circle_radius() const;
    double get_circle_offset() const;
    double get_area() const;
    double get_bounding_radius() const;

private:
    vec center_;
    double theta_ = 0.0;
    double length_ = 1.0;
    double width_ = 1.0;
};
