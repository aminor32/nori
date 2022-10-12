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
        radiance = props.getColor("radiance");
    }

    Color3f Le(const Mesh &mesh) const {
        /* std::cout << "Le" << std::endl;
        Sample lightSample = mesh.sampleMesh();

        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_real_distribution<float> dist(0, 1);

        Point2f uSquareSample = Point2f(dist(rng), dist(rng));
        Vector3f wo = Warp::squareToCosineHemisphere(uSquareSample);
        Frame lightSampleFrame = Frame(lightSample.n);
        wo = lightSampleFrame.toWorld(wo);

        return lightSample.n.dot(wo) > 0 ? radiance : Color3f(); */

        return Color3f(40, 40, 40);
    }

    std::string toString() const { return "AreaLight[]"; }
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END