#include <nori/bsdf.h>
#include <nori/dpdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <vector>

#define MAX_DEPTH 10

NORI_NAMESPACE_BEGIN

class PathMis : public Integrator {
   public:
    PathMis(const PropertyList &props) {
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
        Color3f Le(0.f), Ld(0.f), Li(0.f);
        Color3f throughput(1.f);
        Ray3f _ray(ray);

        for (int depth = 0; depth < MAX_DEPTH; depth++) {
            // terminate: no intersection between scene and ray
            if (!scene->rayIntersect(_ray, its)) {
                break;
            }

            const Mesh &mesh = *(its.mesh);
            const BSDF &bsdf = *(mesh.getBSDF());

            // sample emitter
            const std::vector<Mesh *> &emitters = scene->getEmitters();
            int emitterIndex = emitterDPDF.sample(sampler->next1D());
            const Mesh *emitter = emitters[emitterIndex];
            float pEmitter = 1.f / emitters.size();

            // sample a point on the emitter
            Sample lightSample = emitter->sampleMesh(sampler);
            Vector3f lightDir = (lightSample.p - its.p).normalized();
            float d2 = (lightSample.p - its.p).squaredNorm();

            if (depth == 0 && mesh.isEmitter()) {
                Le += mesh.getEmitter()->getRadiance();
            }

            if (bsdf.isDiffuse()) {
                BSDFQueryRecord bsdfQR(its.toLocal(lightDir).normalized(),
                                       its.toLocal(-_ray.d).normalized(),
                                       ESolidAngle);

                // change pdf w.r.t solid angle -> area
                float pLight = std::abs(lightSample.n.dot(-_ray.d)) > Epsilon
                                   ? pEmitter * d2 /
                                         (emitter->getDiscretePDF().getSum() *
                                          std::abs(lightSample.n.dot(-_ray.d)))
                                   : Epsilon;
                float pBRDF = bsdf.pdf(bsdfQR);
                float wLight = pLight / (pLight + pBRDF);

                // calculate geometric term
                Ray3f shadowRay(its.p, lightDir);
                shadowRay.maxt = std::sqrt(d2) - Epsilon;
                float visibility = scene->rayIntersect(shadowRay) ? 0.f : 1.f;
                Color3f geometric = visibility *
                                    std::abs(its.shFrame.n.dot(lightDir) *
                                             lightSample.n.dot(lightDir)) /
                                    d2;

                Color3f fr = bsdf.eval(bsdfQR);
                Color3f l = emitter->getEmitter()->Le(lightSample.n, -lightDir);

                Ld += wLight * throughput * fr * geometric * l /
                      (pEmitter * lightSample.pdf);
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

            // sample reflected direction & update ray
            BSDFQueryRecord bsdfQR(its.toLocal(-_ray.d).normalized());
            throughput *= bsdf.sample(bsdfQR, sampler->next2D());
            _ray = Ray3f(its.p, its.toWorld(bsdfQR.wo).normalized());

            // next event estimation
            Intersection nextIts;

            if (!scene->rayIntersect(_ray, nextIts)) {
                break;
            }

            if (nextIts.mesh->isEmitter()) {
                const Mesh &nextMesh = *(nextIts.mesh);
                const BSDF &nextBsdf = *(nextIts.mesh->getBSDF());

                float pBRDF = bsdf.pdf(bsdfQR),
                      pLight = pEmitter * (nextIts.p - its.p).squaredNorm() /
                               (nextMesh.getDiscretePDF().getSum() *
                                std::abs(nextIts.shFrame.n.dot(-_ray.d))),
                      wBRDF = 1.f;

                if (bsdf.isDiffuse()) {
                    wBRDF = pBRDF / (pLight + pBRDF);
                }

                Li += wBRDF * throughput * nextMesh.getEmitter()->getRadiance();
            }
        }

        return Le + Ld + Li;
    }

    std::string toString() const { return "PathMisIntegrator[]"; }

   private:
    DiscretePDF emitterDPDF;
};

NORI_REGISTER_CLASS(PathMis, "path_mis");
NORI_NAMESPACE_END