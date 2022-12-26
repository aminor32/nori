#include <nori/bsdf.h>
#include <nori/dpdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <vector>

#define MAX_DEPTH 20

NORI_NAMESPACE_BEGIN

class PathEms : public Integrator {
   public:
    PathEms(const PropertyList &props) {
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

            // sample emitter
            const std::vector<Mesh *> &emitters = scene->getEmitters();
            int emitterIndex = emitterDPDF.sample(sampler->next1D());
            const Mesh *emitter = emitters[emitterIndex];
            float emitterPdf = emitterDPDF[emitterIndex] || 1;

            // sample a point on the emitter
            Sample emitterSample = emitter->sampleMesh(sampler);
            Vector3f emitterDir = (emitterSample.p - its.p).normalized();
            float d2 = (emitterSample.p - its.p).squaredNorm();

            if (bsdf.isDiffuse()) {
                BSDFQueryRecord bsdfQR = BSDFQueryRecord(
                    its.toLocal(emitterDir).normalized(),
                    its.toLocal(-ray.d).normalized(), ESolidAngle);

                // calculate geometric terms
                Ray3f shadowRay = Ray3f(its.p, emitterDir);
                shadowRay.maxt = std::sqrt(d2) - Epsilon;
                float visibility = scene->rayIntersect(shadowRay) ? 0.f : 1.f;
                Color3f geometric =
                    visibility *
                    std::abs((its.shFrame.cosTheta(
                                 its.toLocal(emitterDir).normalized())) *
                             emitterSample.n.dot(emitterDir)) /
                    d2;

                // calculate Le and update L
                Color3f Le =
                    emitter->getEmitter()->Le(emitterSample.n, -emitterDir);
                L += fr * bsdf.eval(bsdfQR) * geometric * Le /
                     (emitterSample.pdf * emitterPdf);
            }

            // sample reflected direction
            BSDFQueryRecord bsdfQR =
                BSDFQueryRecord(its.toLocal(-_ray.d).normalized());
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

            // next event estimation
            Intersection nextIts;
            if (bsdf.isDiffuse() && scene->rayIntersect(_ray, nextIts) &&
                nextIts.mesh == emitter) {
                // L -= fr * nextIts.mesh->getEmitter()->getRadiance();
                break;
            }
        }

        return L;
    }

    std::string toString() const { return "PathEmsIntegrator[]"; }

   private:
    DiscretePDF emitterDPDF;
};

NORI_REGISTER_CLASS(PathEms, "path_ems");
NORI_NAMESPACE_END