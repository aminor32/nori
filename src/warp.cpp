/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/frame.h>
#include <nori/vector.h>
#include <nori/warp.h>

#include <cmath>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) { return sample; }

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f
                                                                        : 0.0f;
}

float Warp::tentPdf(const float &t) {
    if (t < -1 || 1 < t) {
        return 0;
    } else {
        return 1 - std::abs(t);
    }
}

float Warp::tentCdfInverse(const float &t) {
    if (t < 0 || t > 1) {
        throw NoriException("tent CDF inverse out of range");
    } else if (0 <= t && t <= 0.5) {
        return -1 + std::sqrt(2 * t);
    } else {
        return 1 - std::sqrt(2 - 2 * t);
    }
}

Point2f Warp::squareToTent(const Point2f &sample) {
    return Point2f(tentCdfInverse(sample.x()), tentCdfInverse(sample.y()));
}

float Warp::squareToTentPdf(const Point2f &p) {
    return tentPdf(p.x()) * tentPdf(p.y());
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    // Point2f v = 2 * sample - Point2f(1.f, 1.f);
    float r, theta;

    /* if (v.x() == 0 && v.y() == 0) {
        return Point2f();
    } else if (std::abs(v.x()) > std::abs(v.y())) {
        r = v.x();
        theta = (M_PI / 4) * (v.y() / v.x());
    } else {
        r = v.y();
        theta = M_PI / 2 - (M_PI / 4) * (v.x() / v.y());
    } */

    r = std::sqrt(sample.x());
    theta = 2 * M_PI * sample.y();

    return r * Point2f(std::cos(theta), std::sin(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    if (std::sqrt(p.dot(p)) > 1) {
        return 0;
    } else {
        return 1 / M_PI;
    }
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    const float &x = sample.x(), &y = sample.y();

    return Vector3f(2 * std::sqrt(x - x * x) * std::cos(2 * M_PI * y),
                    2 * std::sqrt(x - x * x) * std::sin(2 * M_PI * y),
                    1 - 2 * x);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return 1 / (4 * M_PI);
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    const float &x = sample.x(), &y = sample.y();

    return Vector3f(std::sqrt(1 - x * x) * std::cos(2 * M_PI * y),
                    std::sqrt(1 - x * x) * std::sin(2 * M_PI * y), x);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    if (v.z() > 0) {
        return 1 / (2 * M_PI);
    } else {
        return 0;
    }
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    Point2f d = squareToUniformDisk(sample);

    return Vector3f(d.x(), d.y(), std::sqrt(1 - d.x() * d.x() - d.y() * d.y()));
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    if (v.z() > 0) {
        return v.z() / M_PI;
    } else {
        return 0;
    }
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    float pi = 2 * M_PI * sample.x();
    float theta = std::atan(std::sqrt(-alpha * alpha * std::log(sample.y())));

    return Vector3f(std::sin(theta) * std::cos(pi),
                    std::sin(theta) * std::sin(pi), std::cos(theta));
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    float alpha2 = alpha * alpha;
    float tanTheta = Frame::tanTheta(m);
    float cosTheta = Frame::cosTheta(m);

    if (cosTheta <= 0) {
        return 0;
    } else {
        return (0.5 * INV_PI) * (2 * std::exp(-tanTheta * tanTheta / alpha2)) /
               (alpha2 * std::pow(cosTheta, 3));
    }
}

NORI_NAMESPACE_END
