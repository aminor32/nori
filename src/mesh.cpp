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

#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/dpdf.h>
#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/warp.h>

#include <Eigen/Geometry>
#include <cmath>
#include <string>

NORI_NAMESPACE_BEGIN

Mesh::Mesh() {}

Mesh::~Mesh() {
    delete m_bsdf;
    delete m_emitter;
}

void Mesh::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
    }

    // generate discrete probability function
    const uint32_t nTriangle = (uint32_t)m_F.cols();
    dpdf = DiscretePDF(nTriangle);

    for (uint32_t f = 0; f < nTriangle; f++) {
        float area = surfaceArea(f);

        dpdf.append(area);
        m_surface += area;
    }

    dpdf.normalize();
}

float Mesh::surfaceArea(uint32_t index) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v,
                        float &t) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    /* Find vectors for two edges sharing v[0] */
    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = ray.d.cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f) return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = ray.o - p0;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0) return false;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = ray.d.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0) return false;

    /* Ray intersects triangle -> compute t */
    t = edge2.dot(qvec) * inv_det;

    return t >= ray.mint && t <= ray.maxt;
}

BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
    BoundingBox3f result(m_V.col(m_F(0, index)));
    result.expandBy(m_V.col(m_F(1, index)));
    result.expandBy(m_V.col(m_F(2, index)));
    return result;
}

Point3f Mesh::getCentroid(uint32_t index) const {
    return (1.0f / 3.0f) * (m_V.col(m_F(0, index)) + m_V.col(m_F(1, index)) +
                            m_V.col(m_F(2, index)));
}

void Mesh::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EBSDF:
            if (m_bsdf)
                throw NoriException(
                    "Mesh: tried to register multiple BSDF instances!");
            m_bsdf = static_cast<BSDF *>(obj);
            break;

        case EEmitter: {
            Emitter *emitter = static_cast<Emitter *>(obj);
            if (m_emitter != nullptr)
                throw NoriException(
                    "Mesh: tried to register multiple Emitter instances!");
            m_emitter = emitter;
        } break;

        default:
            throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}

Sample Mesh::sampleMesh(Sampler *sampler) const {
    // value for barycentric coordinate
    float zeta1 = sampler->next1D(), zeta2 = sampler->next1D();
    float alpha = 1 - std::sqrt(1 - zeta1), beta = zeta2 * std::sqrt(1 - zeta1);

    // sample face
    uint32_t sampleFace = dpdf.sample(sampler->next1D());
    uint32_t idx0 = m_F(0, sampleFace), idx1 = m_F(1, sampleFace),
             idx2 = m_F(2, sampleFace);

    // sample point on the mesh
    Point3f p0 = m_V.col(idx0), p1 = m_V.col(idx1), p2 = m_V.col(idx2);
    Point3f samplePoint = alpha * p0 + beta * p1 + (1 - alpha - beta) * p2;

    // sample normal
    Normal3f normal;
    if (m_N.size() > 0) {
        // if vertex normal is given
        Normal3f n0 = m_N.col(idx0), n1 = m_N.col(idx1), n2 = m_N.col(idx2);

        normal = alpha * n0 + beta * n1 + (1 - alpha - beta) * n2;
    } else if (m_emitter) {
        // if vertex normal is not given calculate face normal
        normal = (p1 - p0).cross(p2 - p0);
    }

    // get pdf of the sampled face
    float pdf = 1 / dpdf.getSum();

    return Sample(samplePoint, normal.normalized(), pdf);
}

std::string Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  name = \"%s\",\n"
        "  vertexCount = %i,\n"
        "  triangleCount = %i,\n"
        "  bsdf = %s,\n"
        "  emitter = %s\n"
        "]",
        m_name, m_V.cols(), m_F.cols(),
        m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
        m_emitter ? indent(m_emitter->toString()) : std::string("null"));
}

std::string Intersection::toString() const {
    if (!mesh) return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(), t, uv.toString(), indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null"));
}

NORI_NAMESPACE_END
