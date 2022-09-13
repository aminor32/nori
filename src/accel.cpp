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

#include <set>
#include <nori/accel.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

bool Triangle::operator<(const Triangle &to) const
{
    return false;
}

OctreeNode::OctreeNode(BoundingBox3f inputBoundingBox = TBoundingBox<Point3f>(),
                       std::set<Triangle> inputTriangles = {},
                       std::set<OctreeNode *> inputChildren = {})
{
    boundingBox = inputBoundingBox;
    triangles = inputTriangles;
    children = inputChildren;

    buildChildren();
};

void OctreeNode::buildChildren()
{
    if (triangles.size() < 10)
    {
        std::cout << triangles.size() << std::endl;
        return;
    }
    else
    {
        // min에 대한 max의 방향
        Point3f direction = Point3f(1.0, 1.0, 1.0);
        Point3f min = boundingBox.min;
        Point3f mid = boundingBox.getCenter();
        // length of sub-bounding box's edge
        float halfEdge = (mid(0, 0) - min(0, 0));

        // sub-bounding box에 속하는 triangles를 저장하는 배열
        std::set<Triangle> childTriangles[8];
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
            subBoundingBoxes[i] = BoundingBox3f(mins[i], mins[i] + halfEdge * direction);
        }

        std::set<Triangle>::iterator it;
        for (it = triangles.begin(); it == triangles.end(); ++it)
        {
            Triangle triangle = *it;
            // construct bounding box of triangle
            BoundingBox3f triBoundingBox = BoundingBox3f();
            triBoundingBox.expandBy(triangle.v0);
            triBoundingBox.expandBy(triangle.v1);
            triBoundingBox.expandBy(triangle.v2);

            for (int i = 0; i < 8; ++i)
            {
                // check bounding box of triangle and sub-bounding box overlaps
                if (subBoundingBoxes[i].overlaps(triBoundingBox))
                {
                    childTriangles[i].insert(triangle);
                }
            }
        }

        for (int i = 0; i < 8; ++i)
        {
            children.insert(new OctreeNode(subBoundingBoxes[i], childTriangles[i]));
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
    uint32_t triangleNum = m_mesh->getTriangleCount();
    const MatrixXu &faces = m_mesh->getIndices();
    const MatrixXf &vertices = m_mesh->getVertexPositions();
    std::set<Triangle> totalTriagles;

    for (uint32_t i = 0; i < triangleNum; ++i)
    {
        uint32_t i0 = faces(0, i), i1 = faces(1, i), i2 = faces(2, i);
        const Point3f v0 = vertices.col(i0), v1 = vertices.col(i1), v2 = vertices.col(i2);
        const Triangle currentTriangle = {i, v0, v1, v2};

        totalTriagles.insert(currentTriangle);
    }

    std::cout << totalTriagles.size() << std::endl;

    OctreeNode *root = new OctreeNode(meshBoundingBox, totalTriagles);

    return root;
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