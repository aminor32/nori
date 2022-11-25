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

    Color3f sampleDirectLight(const Scene *scene, Sampler *sampler,
                              const Ray3f &ray, Intersection &its) const {
        // sample emitter
        const std::vector<Mesh *> &emitters = scene->getEmitters();
        int emitterIndex = emitterDPDF.sample(sampler->next1D());
        const Mesh *emitter = emitters[emitterIndex];
        float emitterPdf = emitterDPDF[emitterIndex] || 1;

        // sample a point on the light source
        Sample lightSample = emitter->sampleMesh(sampler);
        Vector3f lightDir = (lightSample.p - its.p).normalized();
        float d2 = (lightSample.p - its.p).squaredNorm();

        const BSDF &bsdf = *(its.mesh->getBSDF());
        BSDFQueryRecord bsdfQR = BSDFQueryRecord(its.toLocal(lightDir));

        // sample to get measure of BSDF
        bsdf.sample(bsdfQR, sampler->next2D());
        bsdfQR.wo = its.toLocal(-ray.d);
        Color3f fr = bsdf.eval(bsdfQR);

        // calculate geometric terms
        Ray3f shadowRay = Ray3f(its.p, lightDir);
        shadowRay.maxt = std::sqrt(d2) - Epsilon;
        float visibility = scene->rayIntersect(shadowRay) ? 0.f : 1.f;
        Color3f geometric =
            visibility *
            std::abs(
                (its.shFrame.cosTheta(its.toLocal(lightDir).normalized())) *
                lightSample.n.dot(lightDir)) /
            d2;

        // calculate Le
        Color3f Le = emitter->getEmitter()->Le(lightSample.n, -lightDir);

        return fr * geometric * Le / (lightSample.pdf * emitterPdf);
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;

        if (!scene->rayIntersect(ray, its)) {
            return Color3f();
        } else {
            Color3f L = Color3f();
            Ray3f _ray = Ray3f(ray);

            // add Le if first met object is emitter
            if (its.mesh->isEmitter()) {
                L += its.mesh->getEmitter()->getRadiance();
            }

            float throughput = 1.f;
            float eta = 1.f;
            float inf = 1.f;

            for (int i = 0; i < MAX_DEPTH; i++) {
                const Mesh &mesh = *(its.mesh);
                const BSDF &bsdf = *(mesh.getBSDF());
                BSDFQueryRecord bsdfQR =
                    BSDFQueryRecord(its.shFrame.toLocal(-_ray.d).normalized());
                // sample reflected direction on the material
                Color3f samplingWeight = bsdf.sample(bsdfQR, sampler->next2D());
                Color3f fr = bsdf.eval(bsdfQR);

                // generate shadow ray
                Ray3f shadowRay =
                    Ray3f(its.p, its.toWorld(bsdfQR.wo).normalized());
                Intersection shadowIts;

                // terminate: no intersection between ray and scene
                if (!scene->rayIntersect(shadowRay, shadowIts)) {
                    break;
                }

                L += inf * samplingWeight *
                     sampleDirectLight(scene, sampler, _ray, its);

                // update ray to shadow ray to find next intersection
                _ray = Ray3f(its.p, its.toWorld(bsdfQR.wo).normalized());

                // terminate: Russian roulette
                if (i > 3) {
                    float zeta = sampler->next1D();

                    // calculate geometric term to update throughput
                    Vector3f dir = its.toWorld(bsdfQR.wo).normalized();
                    float d2 = (shadowIts.p - its.p).squaredNorm();
                    Color3f geometric =
                        std::abs(its.shFrame.cosTheta(bsdfQR.wo) *
                                 shadowIts.shFrame.cosTheta(
                                     shadowIts.toLocal(dir).normalized())) /
                        d2;

                    // update throughput
                    throughput = std::max<float>(throughput,
                                                 (fr * geometric).maxCoeff());
                    // update eta
                    eta *= bsdfQR.eta;
                    // update inf (infimum of Russian roulette)
                    inf = std::min<float>(throughput * eta * eta, 0.99);

                    if (zeta < inf) {
                        break;
                    }
                }
            }
            return L;
        }
    }

    std::string toString() const { return "PathMatsIntegrator[]"; }

   private:
    DiscretePDF emitterDPDF;
};

NORI_REGISTER_CLASS(PathMats, "path_mats");
NORI_NAMESPACE_END