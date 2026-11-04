#pragma once

#include "particle.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

// Exact curved-particle geometry. Bounding circles are used only for early
// rejection; the final decisions use the analytic ellipse contact function
// or analytic intersections of the circular arcs forming a sharp rod.
namespace hard_overlap
{
namespace detail
{
constexpr double relative_tolerance = 1.0e-12;

struct Matrix2
{
    double xx = 0.0;
    double xy = 0.0;
    double yy = 0.0;
};

inline Matrix2 ellipse_matrix(const EllipseParticle& particle)
{
    const double a2 = 0.25 * particle.get_length() * particle.get_length();
    const double b2 = 0.25 * particle.get_width() * particle.get_width();
    const double c = std::cos(particle.get_theta());
    const double s = std::sin(particle.get_theta());
    return Matrix2{
        a2 * c * c + b2 * s * s,
        (a2 - b2) * c * s,
        a2 * s * s + b2 * c * c
    };
}

inline Matrix2 blend(const Matrix2& a, const Matrix2& b, double lambda)
{
    const double one_minus = 1.0 - lambda;
    return Matrix2{
        one_minus * a.xx + lambda * b.xx,
        one_minus * a.xy + lambda * b.xy,
        one_minus * a.yy + lambda * b.yy
    };
}

inline vec solve(const Matrix2& matrix, const vec& rhs)
{
    const double determinant = matrix.xx * matrix.yy - matrix.xy * matrix.xy;
    return vec(
        (matrix.yy * rhs.x - matrix.xy * rhs.y) / determinant,
        (matrix.xx * rhs.y - matrix.xy * rhs.x) / determinant);
}

inline double ellipse_contact_value(const Matrix2& a, const Matrix2& b,
    const vec& separation, double lambda)
{
    const Matrix2 matrix = blend(a, b, lambda);
    const vec inverse_r = solve(matrix, separation);
    return lambda * (1.0 - lambda) * separation.dot(inverse_r);
}

inline double ellipse_contact_derivative(const Matrix2& a, const Matrix2& b,
    const vec& separation, double lambda)
{
    const Matrix2 matrix = blend(a, b, lambda);
    const vec u = solve(matrix, separation);
    const double q = separation.dot(u);
    const double dxx = b.xx - a.xx;
    const double dxy = b.xy - a.xy;
    const double dyy = b.yy - a.yy;
    const double q_prime = -(dxx * u.x * u.x + 2.0 * dxy * u.x * u.y
        + dyy * u.y * u.y);
    return (1.0 - 2.0 * lambda) * q
        + lambda * (1.0 - lambda) * q_prime;
}

inline vec to_local(const vec& world_point, const SharpRodParticle& particle)
{
    return (world_point - particle.get_center()).rotated(-particle.get_theta());
}

inline vec to_world(const vec& local_point, const SharpRodParticle& particle)
{
    return particle.get_center() + local_point.rotated(particle.get_theta());
}

inline bool point_in_lens_local(const vec& point, double radius, double offset,
    double tolerance)
{
    const vec upper_center(0.0, offset);
    const vec lower_center(0.0, -offset);
    const double limit = radius * radius + tolerance;
    return (point - upper_center).mod2() <= limit
        && (point - lower_center).mod2() <= limit;
}

inline bool point_in_lens(const vec& point, const SharpRodParticle& particle,
    double tolerance)
{
    return point_in_lens_local(to_local(point, particle),
        particle.get_circle_radius(), particle.get_circle_offset(), tolerance);
}

inline int circle_intersections(const vec& c0, double r0, const vec& c1,
    double r1, vec output[2], double tolerance)
{
    const vec delta = c1 - c0;
    const double distance2 = delta.mod2();
    const double distance = std::sqrt(distance2);
    if (distance <= tolerance)
    {
        return 0; // Coincident boundaries are covered by lens arc endpoints.
    }
    if (distance > r0 + r1 + tolerance
        || distance < std::abs(r0 - r1) - tolerance)
    {
        return 0;
    }

    const double along = (r0 * r0 - r1 * r1 + distance2) / (2.0 * distance);
    double height2 = r0 * r0 - along * along;
    if (height2 < -tolerance) return 0;
    height2 = std::max(0.0, height2);

    const vec direction = delta / distance;
    const vec base = c0 + direction * along;
    const vec perpendicular(-direction.y, direction.x);
    const double height = std::sqrt(height2);
    output[0] = base + perpendicular * height;
    if (height <= tolerance) return 1;
    output[1] = base - perpendicular * height;
    return 2;
}

inline double lens_support_local(const SharpRodParticle& particle,
    const vec& direction_world)
{
    const vec direction = direction_world.rotated(-particle.get_theta()).normalized();
    const double a = particle.get_length() / 2.0;
    const double radius = particle.get_circle_radius();
    const double offset = particle.get_circle_offset();
    const double tolerance = relative_tolerance * radius * radius;

    double support = std::max(direction.dot(vec(a, 0.0)),
        direction.dot(vec(-a, 0.0)));
    const std::array<vec, 2> circle_centers = {
        vec(0.0, offset), vec(0.0, -offset)
    };
    for (const vec& circle_center : circle_centers)
    {
        const vec candidate = circle_center + direction * radius;
        if (point_in_lens_local(candidate, radius, offset, tolerance))
        {
            support = std::max(support, direction.dot(candidate));
        }
    }
    return support;
}
}

inline bool overlap(const EllipseParticle& a, const EllipseParticle& b)
{
    const vec separation = b.get_center() - a.get_center();
    const double distance2 = separation.mod2();
    const double bound = a.get_bounding_radius() + b.get_bounding_radius();
    if (distance2 > bound * bound) return false;

    // Each ellipse contains a centered disk whose radius is half its width.
    // Overlap of these two disks guarantees overlap of the ellipses.
    const double inner_bound = 0.5 * (a.get_width() + b.get_width());
    if (distance2 <= inner_bound * inner_bound) return true;

    const detail::Matrix2 matrix_a = detail::ellipse_matrix(a);
    const detail::Matrix2 matrix_b = detail::ellipse_matrix(b);

    // The Perram-Wertheim contact function is concave on [0,1]. Its
    // derivative is positive at 0 and negative at 1, so safeguarded bisection
    // locates the unique maximum without polygonizing either ellipse.
    double lower = 0.0;
    double upper = 1.0;
    for (int iteration = 0; iteration < 80; ++iteration)
    {
        const double middle = 0.5 * (lower + upper);
        if (detail::ellipse_contact_derivative(
                matrix_a, matrix_b, separation, middle) > 0.0)
        {
            lower = middle;
        }
        else
        {
            upper = middle;
        }
    }
    const double lambda = 0.5 * (lower + upper);
    const double contact = detail::ellipse_contact_value(
        matrix_a, matrix_b, separation, lambda);
    return contact <= 1.0 + detail::relative_tolerance;
}

inline bool overlap(const SharpRodParticle& a, const SharpRodParticle& b)
{
    const vec separation = b.get_center() - a.get_center();
    const double distance2 = separation.mod2();
    const double bound = a.get_bounding_radius() + b.get_bounding_radius();
    if (distance2 > bound * bound) return false;

    // The intersection of the two mother circles contains a centered disk
    // whose radius is half the sharp rod's width.
    const double inner_bound = 0.5 * (a.get_width() + b.get_width());
    if (distance2 <= inner_bound * inner_bound) return true;

    const double radius_a = a.get_circle_radius();
    const double radius_b = b.get_circle_radius();
    const double offset_a = a.get_circle_offset();
    const double offset_b = b.get_circle_offset();
    const double scale = std::max(radius_a, radius_b);
    const double tolerance = detail::relative_tolerance * scale * scale;

    if (detail::point_in_lens(a.get_center(), b, tolerance)
        || detail::point_in_lens(b.get_center(), a, tolerance))
    {
        return true;
    }

    // Arc endpoints also cover coincident mother-circle boundaries.
    const std::array<vec, 4> tips_a = {
        detail::to_world(vec(a.get_length() / 2.0, 0.0), a),
        detail::to_world(vec(-a.get_length() / 2.0, 0.0), a),
        detail::to_world(vec(0.0, a.get_width() / 2.0), a),
        detail::to_world(vec(0.0, -a.get_width() / 2.0), a)
    };
    const std::array<vec, 4> tips_b = {
        detail::to_world(vec(b.get_length() / 2.0, 0.0), b),
        detail::to_world(vec(-b.get_length() / 2.0, 0.0), b),
        detail::to_world(vec(0.0, b.get_width() / 2.0), b),
        detail::to_world(vec(0.0, -b.get_width() / 2.0), b)
    };
    for (const vec& point : tips_a)
        if (detail::point_in_lens(point, b, tolerance)) return true;
    for (const vec& point : tips_b)
        if (detail::point_in_lens(point, a, tolerance)) return true;

    const vec normal_a = vec(0.0, 1.0).rotated(a.get_theta());
    const vec normal_b = vec(0.0, 1.0).rotated(b.get_theta());
    const std::array<vec, 2> centers_a = {
        a.get_center() + normal_a * offset_a,
        a.get_center() - normal_a * offset_a
    };
    const std::array<vec, 2> centers_b = {
        b.get_center() + normal_b * offset_b,
        b.get_center() - normal_b * offset_b
    };

    for (const vec& circle_a : centers_a)
    {
        for (const vec& circle_b : centers_b)
        {
            vec intersections[2];
            const int count = detail::circle_intersections(circle_a, radius_a,
                circle_b, radius_b, intersections,
                detail::relative_tolerance * scale);
            for (int index = 0; index < count; ++index)
            {
                if (detail::point_in_lens(intersections[index], a, tolerance)
                    && detail::point_in_lens(intersections[index], b, tolerance))
                {
                    return true;
                }
            }
        }
    }
    return false;
}

inline bool overlap_periodic(const EllipseParticle& a,
    const EllipseParticle& b, double box_length)
{
    EllipseParticle image = b;
    vec delta = b.get_center() - a.get_center();
    delta.x -= std::round(delta.x / box_length) * box_length;
    delta.y -= std::round(delta.y / box_length) * box_length;
    image.translate(a.get_center() + delta);
    return overlap(a, image);
}

inline bool overlap_periodic(const SharpRodParticle& a,
    const SharpRodParticle& b, double box_length)
{
    SharpRodParticle image = b;
    vec delta = b.get_center() - a.get_center();
    delta.x -= std::round(delta.x / box_length) * box_length;
    delta.y -= std::round(delta.y / box_length) * box_length;
    image.translate(a.get_center() + delta);
    return overlap(a, image);
}

inline bool boundary_overlap(const EllipseParticle& particle, double box_length)
{
    const double a = particle.get_length() / 2.0;
    const double b = particle.get_width() / 2.0;
    const double c = std::cos(particle.get_theta());
    const double s = std::sin(particle.get_theta());
    const double extent_x = std::sqrt(a * a * c * c + b * b * s * s);
    const double extent_y = std::sqrt(a * a * s * s + b * b * c * c);
    const double half_box = box_length / 2.0;
    const vec& center = particle.get_center();
    return std::abs(center.x) + extent_x > half_box
        || std::abs(center.y) + extent_y > half_box;
}

inline bool boundary_overlap(const SharpRodParticle& particle, double box_length)
{
    const double extent_x = detail::lens_support_local(particle, vec(1.0, 0.0));
    const double extent_y = detail::lens_support_local(particle, vec(0.0, 1.0));
    const double half_box = box_length / 2.0;
    const vec& center = particle.get_center();
    return std::abs(center.x) + extent_x > half_box
        || std::abs(center.y) + extent_y > half_box;
}
}
