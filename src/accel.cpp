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

#include <nori/accel.h>

#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>

NORI_NAMESPACE_BEGIN

OctreeNode::OctreeNode(Mesh *inputMesh, uint32_t inputDepth,
                       BoundingBox3f &inputBoundingBox,
                       std::vector<uint32_t> &inputTriangles)
    : mesh(inputMesh),
      depth(inputDepth),
      boundingBox(BoundingBox3f(inputBoundingBox.min, inputBoundingBox.max)) {
    if (inputDepth > 10 || inputTriangles.size() < 10) {
        // leaf
        triangles = std::vector<uint32_t>(inputTriangles);

        return;
    } else {
        // not leaf
        const Point3f &min = inputBoundingBox.min;
        const Point3f &center = inputBoundingBox.getCenter();
        const Point3f &extents = inputBoundingBox.getExtents();

        if (std::min(
                {0.5 * extents.x(), 0.5 * extents.y(), 0.5 * extents.z()}) <=
            std::numeric_limits<float>::epsilon()) {
            triangles = std::vector<uint32_t>(inputTriangles);

            return;
        }

        // sub-bounding box의 min을 저장하는 배열
        Point3f mins[8] = {
            min,
            Point3f(min.x(), center.y(), min.z()),
            Point3f(center.x(), center.y(), min.z()),
            Point3f(center.x(), min.y(), min.z()),
            Point3f(min.x(), min.y(), center.z()),
            Point3f(min.x(), center.y(), center.z()),
            center,
            Point3f(center.x(), min.y(), center.z()),
        };
        BoundingBox3f subBoundingBoxes[8];

        for (int i = 0; i < 8; ++i) {
            subBoundingBoxes[i] =
                BoundingBox3f(mins[i], mins[i] + 0.5 * extents);
        }

        // sub-bounding box에 속하는 triangles를 저장하는 배열
        std::array<std::vector<uint32_t>, 8> childTriangles = {{}};

        std::vector<uint32_t>::iterator it;
        for (it = inputTriangles.begin(); it != inputTriangles.end(); ++it) {
            const uint32_t idx = *it;
            const MatrixXu &faces = mesh->getIndices();
            const MatrixXf &vertices = mesh->getVertexPositions();

            uint32_t i0 = faces(0, idx), i1 = faces(1, idx), i2 = faces(2, idx);
            Point3f v0 = vertices.col(i0), v1 = vertices.col(i1),
                    v2 = vertices.col(i2);

            // construct bounding box of triangle
            BoundingBox3f triBoundingBox = BoundingBox3f(v0);
            triBoundingBox.expandBy(v1);
            triBoundingBox.expandBy(v2);

            for (int i = 0; i < 8; ++i) {
                // check bounding box of triangle and sub-bounding box
                // overlaps
                if (subBoundingBoxes[i].overlaps(triBoundingBox)) {
                    childTriangles[i].push_back(idx);
                }
            }
        }

        for (int i = 0; i < 8; ++i) {
            children[i] = new OctreeNode(mesh, depth + 1, subBoundingBoxes[i],
                                         childTriangles[i]);
        }
    }
};

std::string OctreeNode::toString() {
    return tfm::format("OctreeNode: depth=%d, minPoint=%s, trianglesNum=%d",
                       depth, boundingBox.min.toString(), triangles.size());
}

void Accel::addMesh(Mesh *mesh) {
    m_mesh.push_back(mesh);
    m_bbox.expandBy(mesh->getBoundingBox());
}

/* returns reference of root of Octree */
void Accel::build() {
    for (std::vector<Mesh *>::iterator it = m_mesh.begin(); it < m_mesh.end();
         it++) {
        Mesh *mesh = *it;
        BoundingBox3f meshBoundingBox = mesh->getBoundingBox();
        uint32_t triangleNum = mesh->getTriangleCount();
        std::vector<uint32_t> totalTriangles;

        for (uint32_t i = 0; i < triangleNum; ++i) {
            totalTriangles.push_back(i);
        }

        m_root.push_back(
            new OctreeNode(mesh, 0, meshBoundingBox, totalTriangles));
    }

    std::cout << "Octree build done" << std::endl;

    return;
}

void Accel::traverseOctree(OctreeNode *node, Ray3f &ray, Intersection &its,
                           bool shadowRay, bool &foundIntersection,
                           uint32_t &f) const {
    BoundingBox3f &boundingBox = node->boundingBox;
    bool boxIntersection = boundingBox.rayIntersect(ray);

    // not leaf node
    if (boxIntersection && node->children[0] != nullptr) {
        std::array<OctreeNode *, 8> &children = node->children;

        /* std::sort(children.begin(), children.end(),
                  [&ray](const OctreeNode *node1, const OctreeNode *node2) {
                      return node1->boundingBox.distanceTo(ray.o) <
                             node2->boundingBox.distanceTo(ray.o);
                  }); */

        std::array<OctreeNode *, 8>::const_iterator childIt;
        for (childIt = children.begin(); childIt != children.end(); childIt++) {
            const OctreeNode *child = *childIt;

            if (child->boundingBox.distanceTo(ray.o) > ray.maxt) {
                continue;
            }

            traverseOctree(*childIt, ray, its, shadowRay, foundIntersection, f);
        }
    } else if (boxIntersection && node->triangles.size() > 0) {
        for (std::vector<uint32_t>::const_iterator triIt =
                 node->triangles.begin();
             triIt != node->triangles.end(); triIt++) {
            float u, v, t;

            if (node->mesh->rayIntersect(*triIt, ray, u, v, t)) {
                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = node->mesh;
                f = *triIt;
                foundIntersection = true;

                if (shadowRay) {
                    return;
                }
            }
        }
    }
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its,
                         bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t)-1;  // Triangle index of the closest intersection
    Ray3f ray(ray_);  /// Make a copy of the ray (we will need to update its
                      /// '.maxt' value)

    // traverse through meshes
    for (std::vector<OctreeNode *>::const_iterator it = m_root.begin();
         it < m_root.end(); it++) {
        traverseOctree(*it, ray, its, shadowRay, foundIntersection, f);
    }

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.
           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates,
           etc..)
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
            its.uv = bary.x() * UV.col(idx0) + bary.y() * UV.col(idx1) +
                     bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1 - p0).cross(p2 - p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame =
                Frame((bary.x() * N.col(idx0) + bary.y() * N.col(idx1) +
                       bary.z() * N.col(idx2))
                          .normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END