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

#include <nori/bsdf.h>
#include <nori/common.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
   public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        Vector3f n;
        float ni, no;
        float sinTheta_i = Frame::sinTheta(bRec.wi);
        float cosTheta_i = Frame::cosTheta(bRec.wi);

        if (cosTheta_i > 0) {
            // from outside to inside of the mesh
            n = Vector3f(0, 0, 1);
            ni = m_extIOR;
            no = m_intIOR;
        } else {
            // from inside to outside of the mesh
            n = Vector3f(0, 0, -1);
            ni = m_intIOR;
            no = m_extIOR;
        }

        float sinTheta_o = ni * sinTheta_i / no;
        float cosTheta_o = std::sqrt(1 - sinTheta_o * sinTheta_o);
        float fresnel = nori::fresnel(std::abs(cosTheta_i), ni, no);

        if (sample.x() <= fresnel) {
            // reflection (in local frame)
            bRec.wo = Vector3f(-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z());
            bRec.eta = 1.0f;
        } else {
            // refraction (in local frame)
            bRec.wo = sinTheta_o * (-bRec.wi + cosTheta_i * n) / sinTheta_i -
                      cosTheta_o * n;
            bRec.eta = no / ni;
        }

        bRec.measure = EDiscrete;

        return Color3f(1.f);
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }

   private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
