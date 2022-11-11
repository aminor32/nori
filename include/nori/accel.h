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

#pragma once

#include <nori/mesh.h>

#include <array>
#include <vector>

NORI_NAMESPACE_BEGIN

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */

class OctreeNode {
   public:
    // default constructor
    OctreeNode(Mesh *mesh, uint32_t inputDepth, BoundingBox3f &inputBoundingBox,
               std::vector<uint32_t> &inputTriangles);

    // input mesh 저장
    Mesh *mesh;
    // depth of node
    uint32_t depth;
    // bounding box
    BoundingBox3f boundingBox = BoundingBox3f();
    // bounding box에 포함되는 삼각형 저장 (leaf node만)
    std::vector<uint32_t> triangles = {};
    // 자식 노드의 주소를 array로 저장
    std::array<OctreeNode *, 8> children = {nullptr};

    std::string toString();

   private:
    // 자식 노드를 생성하는 함수
    void buildChildren(std::vector<uint32_t> *inputTriangles);
};

class Accel {
   public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure (currently a no-op)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    void traverseOctree(OctreeNode *node, Ray3f &ray, Intersection &its,
                        bool shadowRay, bool &foundIntersection,
                        uint32_t &f) const;
    bool rayIntersect(const Ray3f &ray, Intersection &its,
                      bool shadowRay) const;

   private:
    std::vector<Mesh *> m_mesh = {};  ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;             ///< Bounding box of the entire scene
    std::vector<OctreeNode *> m_root = {};  ///< Root of acceleration data model
};

NORI_NAMESPACE_END