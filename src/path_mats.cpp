#include <nori/bsdf.h>
#include <nori/dpdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <cmath>
#include <vector>

#define MAX_DEPTH 20

NORI_NAMESPACE_BEGIN

class PathMats : public Integrator {
   public:
    PathMats(const PropertyList &props) {
        // no arguments needed
    }

    void preprocess(const Scene *scene) {
        const std::vector<Mesh *> &emitters = scene->getEmitters();

        for (std::vector<Mesh *>::const_iterator it = emitters.begin();
             it < emitters.end(); ++it) {
            const Mesh *emitter = *it;

            emitterDPDF.append(emitter->getDiscretePDF().getSum());
        }

        emitterDPDF.normalize();
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;
        Color3f L(0.f);
        Color3f fr(1.f);
        Ray3f _ray = Ray3f(ray);

        for (int depth = 0; depth < MAX_DEPTH; depth++) {
            // terminate: no intersection between scene and ray
            if (!scene->rayIntersect(_ray, its)) {
                break;
            }

            const Mesh &mesh = *(its.mesh);
            const BSDF &bsdf = *(mesh.getBSDF());

            if (mesh.isEmitter()) {
                L += fr * mesh.getEmitter()->getRadiance();
                break;
            }

            // sample reflected direction on the material
            BSDFQueryRecord bsdfQR =
                BSDFQueryRecord(its.shFrame.toLocal(-_ray.d).normalized());
            fr *= bsdf.sample(bsdfQR, sampler->next2D());

            // update ray
            _ray = Ray3f(its.p, its.toWorld(bsdfQR.wo).normalized());

            // Russian roulette
            if (depth > 3) {
                float probability = std::min<float>(fr.maxCoeff(), 0.99);
                float zeta = sampler->next1D();

                if (zeta > probability) {
                    break;
                } else {
                    fr /= probability;
                }
            }
        }
        return L;
    }

    std::string toString() const { return "PathMatsIntegrator[]"; }

   private:
    DiscretePDF emitterDPDF;
};

NORI_REGISTER_CLASS(PathMats, "path_mats");
NORI_NAMESPACE_END