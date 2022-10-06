#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/bsdf.h>

#include <string>

NORI_NAMESPACE_BEGIN

class Whitted : public Integrator {
   public:
    Whitted(const PropertyList &props) {
        // no arguments needed
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            return Color3f();
        } else {
            BSDFQueryRecord bsdfQuery = BSDFQueryRecord()
            Color3f fr = its.mesh->getBSDF()->eval()
        }
    }

    std::string toString() const { return "WhittedIntegrator[]"; }
};

NORI_REGISTER_CLASS(Whitted, "whitted");
NORI_NAMESPACE_END