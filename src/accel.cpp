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

#include <cmath>
#include <set>
#include <nori/accel.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

OctreeNode::OctreeNode(Mesh *inputMesh,
                       uint32_t inputDepth,
                       Point3f *inputMin,
                       std::set<uint32_t> *inputTriangles) : children(new std::set<OctreeNode *>)
{
    mesh = inputMesh;
    depth = inputDepth;
    minPoint = inputMin;

    /* if (depth < 100)
    { */
    buildChildren(inputTriangles);
    /* }
    else
    {
        triangles = inputTriangles;
    } */
};

void OctreeNode::buildChildren(
    std::set<uint32_t> *inputTriangles)
{
    if (inputTriangles == nullptr)
    {
        return;
    }
    else if (inputTriangles->size() < 10)
    {
        triangles = inputTriangles;

        return;
    }
    else
    {
        Point3f min = *minPoint;
        // min에 대한 max의 방향
        Point3f direction = Point3f(1.0, 1.0, 1.0);
        // parent's edge
        float edge = pow(2, -depth);
        BoundingBox3f boundingBox = BoundingBox3f(min, min + direction * edge);

        // center of bounding box
        Point3f mid = boundingBox.getCenter();

        // sub-bounding box에 속하는 triangles를 저장하는 배열
        std::set<uint32_t> childTriangles[8];
        // sub-bounding box의 min을 저장하는 배열
        Point3f mins[8] = {
            min,
            Point3f(min(0, 0), mid(1, 0), min(2, 0)),
            Point3f(mid(0, 0), mid(1, 0), min(2, 0)),
            Point3f(mid(0, 0), min(1, 0), min(2, 0)),
            Point3f(min(0, 0), min(1, 0), mid(2, 0)),
            Point3f(min(0, 0), mid(1, 0), mid(2, 0)),
            mid,
            Point3f(mid(0, 0), min(1, 0), mid(2, 0)),
        };
        BoundingBox3f subBoundingBoxes[8];

        for (int i = 0; i < 8; ++i)
        {
            subBoundingBoxes[i] = BoundingBox3f(mins[i], mins[i] + 0.5 * edge * direction);
        }

        std::set<uint32_t>::iterator it;
        for (it = inputTriangles->begin(); it != inputTriangles->end(); ++it)
        {
            const uint32_t idx = *it;
            const MatrixXu &faces = mesh->getIndices();
            const MatrixXf &vertices = mesh->getVertexPositions();

            uint32_t i0 = faces(0, idx), i1 = faces(1, idx), i2 = faces(2, idx);
            Point3f v0 = vertices.col(i0), v1 = vertices.col(i1), v2 = vertices.col(i2);

            // construct bounding box of triangle
            BoundingBox3f triBoundingBox = BoundingBox3f();
            triBoundingBox.expandBy(v0);
            triBoundingBox.expandBy(v1);
            triBoundingBox.expandBy(v2);

            for (int i = 0; i < 8; ++i)
            {
                // check bounding box of triangle and sub-bounding box overlaps
                if (subBoundingBoxes[i].overlaps(triBoundingBox))
                {
                    childTriangles[i].insert(idx);
                }
            }
        }

        for (int j = 0; j < 8; ++j)
        {
            children->insert(new OctreeNode(mesh, ++depth, &mins[j], &childTriangles[j]));
        }
    }
};

void Accel::addMesh(Mesh *mesh)
{
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

/* returns reference of root of Octree */
OctreeNode *Accel::build()
{
    BoundingBox3f meshBoundingBox = m_mesh->getBoundingBox();
    Point3f min = meshBoundingBox.min;
    uint32_t triangleNum = m_mesh->getTriangleCount();
    std::set<uint32_t> totalTriangles;

    for (uint32_t i = 0; i < triangleNum; ++i)
    {
        totalTriangles.insert(i);
    }

    return new OctreeNode(m_mesh, 0, &min, &totalTriangles);
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const
{
    bool foundIntersection = false; // Was an intersection found so far?
    uint32_t f = (uint32_t)-1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    /* Brute force search through all triangles */
    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx)
    {
        float u, v, t;
        if (m_mesh->rayIntersect(idx, ray, u, v, t))
        {
            /* An intersection was found! Can terminate
               immediately if this is a shadow ray query */
            if (shadowRay)
                return true;
            ray.maxt = its.t = t;
            its.uv = Point2f(u, v);
            its.mesh = m_mesh;
            f = idx;
            foundIntersection = true;
        }
    }

    if (foundIntersection)
    {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.
           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1 - its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh = its.mesh;
        const MatrixXf &V = mesh->getVertexPositions();
        const MatrixXf &N = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                     bary.y() * UV.col(idx1) +
                     bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1 - p0).cross(p2 - p0).normalized());

        if (N.size() > 0)
        {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2))
                    .normalized());
        }
        else
        {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END