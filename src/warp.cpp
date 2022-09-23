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
    } else if (-1 <= t && t < 0) {
        return 1 + t;
    } else {
        return 1 - t;
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
    return Point2f(std::sqrt(sample.x()), 2 * M_PI * sample.y());
}

float Warp::squareToUniformDiskPdf(const Point2f &p) { return 1 / M_PI; }

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
    return 1 / (2 * M_PI);
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    Point2f d = squareToUniformDisk(sample);

    return Vector3f(d.x(), d.y(), std::sqrt(1 - d.x() * d.x() - d.y() * d.y()));
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    return v.z() / M_PI;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");
}

NORI_NAMESPACE_END
