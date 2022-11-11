#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/scene.h>
#include <nori/warp.h>

#include <cmath>
#include <random>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
   public:
    AreaLight(const PropertyList &props) {
        m_radiance = props.getColor("radiance");
    }

    Color3f Le(Normal3f n, Vector3f wo) const {
        return n.dot(wo) > 0 ? m_radiance : Color3f();
    }

    Color3f sampleLe(Sampler *sampler, const Mesh &mesh, Vector3f *wo) const {
        Sample lightSample = mesh.sampleMesh(sampler);

        // sample wo in local frame
        Vector3f localWo = Warp::squareToCosineHemisphere(sampler->next2D());

        Frame lightSampleFrame = Frame(lightSample.n);
        *wo = lightSampleFrame.toWorld(localWo);

        return Le(lightSample.n, *wo);
    }

    std::string toString() const { return "AreaLight[]"; }
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END