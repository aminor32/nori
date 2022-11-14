#include <nori/bsdf.h>
#include <nori/dpdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <vector>

#define SAMPLE_NUM 10

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

        if (!scene->rayIntersect(ray, its)) {
            return Color3f();
        } else {
            Color3f Le = Color3f(), Lr = Color3f();
            Ray3f _ray = Ray3f(ray);

            float zeta = sampler->next1D();
            float eta = 1.f;

            for (int i = 0; i < 3 || zeta < eta; i++) {
                // add Le if mesh is emitter
                if (its.mesh->isEmitter()) {
                    Le += its.mesh->getEmitter()->getRadiance();
                    break;
                }

                const Mesh &mesh = *(its.mesh);
                const BSDF &bsdf = *(mesh.getBSDF());

                BSDFQueryRecord bsdfQR =
                    BSDFQueryRecord(its.shFrame.toLocal(-ray.d));
                Color3f samplingWeight = bsdf.sample(bsdfQR, sampler->next2D());
                Color3f fr = bsdf.eval(bsdfQR);

                // generate shadow ray
                Ray3f shadowRay = Ray3f(its.p, bsdfQR.wo);
                Intersection shadowIts;

                if (!scene->rayIntersect(shadowRay, shadowIts)) {
                    break;
                } else {
                    // calculate geometric terms
                    Vector3f dir = (shadowIts.p - its.p).normalized();
                    float d2 = (shadowIts.p - its.p).squaredNorm();

                    Color3f geometric =
                        std::abs((its.shFrame.cosTheta(bsdfQR.wo)) *
                                 shadowIts.shFrame.cosTheta(
                                     shadowIts.toLocal(dir).normalized())) /
                        d2;

                    Lr += (1 / eta) * samplingWeight * fr * geometric;
                }

                // update ray to shadow ray
                _ray = Ray3f(its.p, bsdfQR.wo);
                // update eta
                eta *= bsdfQR.eta;
                // sample zeta
                zeta = sampler->next1D();
            }

            return Le + Lr;
        }
    }

    std::string toString() const { return "PathMatsIntegrator[]"; }

   private:
    DiscretePDF emitterDPDF;
};

NORI_REGISTER_CLASS(PathMats, "path_mats");
NORI_NAMESPACE_END