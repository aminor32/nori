#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>

#include <random>
#include <string>
#include <vector>

NORI_NAMESPACE_BEGIN

class Whitted : public Integrator {
   public:
    Whitted(const PropertyList &props) {
        // no arguments needed
    }

    void preprocess(const Scene *scene) {
        const std::vector<Mesh *> &meshes = scene->getMeshes();

        for (std::vector<Mesh *>::const_iterator it = meshes.begin();
             it < meshes.end(); ++it) {
            if ((*it)->isEmitter()) {
                emitters.push_back(*it);
            }
        }
    }

    Color3f sampleIntegralBody(const Scene *scene, Ray3f &ray,
                               Intersection &its) {
        const Mesh &mesh = *(its.mesh);
        Sample lightSample = mesh.sampleMesh();
        Vector3f wi = lightSample.p - its.p;  // not normalized !

        const BSDF &bsdf = *(its.mesh->getBSDF());
        BSDFQueryRecord bsdfQR =
            BSDFQueryRecord(wi.normalized(), -ray.d, EDiscrete);
        Color3f fr = bsdf.eval(bsdfQR);

        Ray3f shadowRay = Ray3f(its.p, wi.normalized());
        shadowRay.maxt = wi.norm();
        float visibility = scene->rayIntersect(shadowRay) ? 0.f : 1.f;
        Color3f geometric =
            visibility *
            std::abs(its.shFrame.cosTheta(its.shFrame.toLocal(wi)) *
                     lightSample.n.dot(wi)) /
            wi.squaredNorm();

        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_int_distribution<float> dist(0, emitters.size());
        const Mesh *emitter = emitters[dist(rng)];

        return fr * geometric * emitter->getEmitter()->Le(*emitter);
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;

        if (!scene->rayIntersect(ray, its)) {
            return Color3f();
        } else {
            Normal3f itsNormal = its.shFrame.n;
            Vector3f reflection = ray.d - 2 * itsNormal.dot(ray.d) * itsNormal;
        }
    }

    std::string toString() const { return "WhittedIntegrator[]"; }

   private:
    std::vector<const Mesh *> emitters = {};
};

NORI_REGISTER_CLASS(Whitted, "whitted");
NORI_NAMESPACE_END