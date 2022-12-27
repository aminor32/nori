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
        Color3f L(0.f), Ld(0.f);
        Color3f throughput(1.f);
        Ray3f _ray(ray);

        for (int depth = 0; depth < MAX_DEPTH; depth++) {
            // terminate: no intersection between scene and ray
            if (!scene->rayIntersect(_ray, its)) {
                break;
            }

            const Mesh &mesh = *(its.mesh);
            const BSDF &bsdf = *(mesh.getBSDF());

            if (depth == 0 && mesh.isEmitter()) {
                L += mesh.getEmitter()->getRadiance();
            }

            // sample emitter
            const std::vector<Mesh *> &emitters = scene->getEmitters();
            int emitterIndex = emitterDPDF.sample(sampler->next1D());
            const Mesh *emitter = emitters[emitterIndex];
            float pEmitter = emitterDPDF[emitterIndex] || 1;

            // sample a point on the emitter
            Sample lightSample = emitter->sampleMesh(sampler);
            Vector3f lightDir = (lightSample.p - its.p).normalized();
            float d2 = (lightSample.p - its.p).squaredNorm();

            if (bsdf.isDiffuse()) {
                BSDFQueryRecord bsdfQR(its.toLocal(lightDir).normalized(),
                                       its.toLocal(-_ray.d).normalized(),
                                       ESolidAngle);

                // calculate geometric term
                Ray3f shadowRay(its.p, lightDir);
                shadowRay.maxt = std::sqrt(d2) - Epsilon;
                float visibility = scene->rayIntersect(shadowRay) ? 0.f : 1.f;
                Color3f geometric =
                    visibility *
                    std::abs((its.shFrame.cosTheta(
                                 its.toLocal(lightDir).normalized())) *
                             lightSample.n.dot(lightDir)) /
                    d2;

                Color3f fr = bsdf.eval(bsdfQR);
                Color3f Le =
                    emitter->getEmitter()->Le(lightSample.n, -lightDir);

                Ld = throughput * fr * geometric * Le /
                     (lightSample.pdf * pEmitter);
                L += Ld;
            }

            // sample reflected direction & update ray
            BSDFQueryRecord bsdfQR(its.toLocal(-_ray.d).normalized());
            throughput *= bsdf.sample(bsdfQR, sampler->next2D());
            _ray = Ray3f(its.p, its.toWorld(bsdfQR.wo).normalized());

            // next event estimation
            Intersection nextIts;
            if (bsdf.isDiffuse() && scene->rayIntersect(_ray, nextIts) &&
                nextIts.mesh == emitter) {
                break;
            }

            // Russian roulette
            if (depth > 3) {
                float probability =
                    std::min<float>(throughput.maxCoeff(), 0.99);
                float zeta = sampler->next1D();

                if (zeta > probability) {
                    break;
                } else {
                    throughput /= probability;
                }
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