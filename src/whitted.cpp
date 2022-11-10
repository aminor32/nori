#include <nori/bsdf.h>
#include <nori/dpdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <random>
#include <vector>

#define SAMPLE_NUM 10

NORI_NAMESPACE_BEGIN

class Whitted : public Integrator {
   public:
    Whitted(const PropertyList &props) {
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

    Color3f sampleIntegralBody(const Scene *scene, Sampler *sampler,
                               const Ray3f &ray, Intersection &its) const {
        const std::vector<Mesh *> &emitters = scene->getEmitters();
        int sample = emitterDPDF.sample(sampler->next1D());
        const Mesh *emitter = emitters[sample];

        // sample from light source
        Sample lightSample = emitter->sampleMesh();
        Vector3f lightDir = lightSample.p - its.p;
        float d2 = lightDir.squaredNorm();
        lightDir.normalize();

        const BSDF &bsdf = *(its.mesh->getBSDF());
        BSDFQueryRecord bsdfQR = BSDFQueryRecord(
            its.toLocal(lightDir), its.toLocal(-ray.d), ESolidAngle);
        Color3f fr = bsdf.eval(bsdfQR);

        // calculate geometric terms
        Ray3f shadowRay = Ray3f(its.p, lightDir);
        shadowRay.maxt = std::sqrt(d2);
        float visibility = scene->rayIntersect(shadowRay) ? 0.f : 1.f;
        Color3f geometric =
            visibility *
            std::abs(
                (its.shFrame.cosTheta(its.toLocal(lightDir).normalized())) *
                lightSample.n.dot(lightDir)) /
            d2;

        // calculate Le
        Color3f Le = emitter->getEmitter()->Le(lightSample.n, -lightDir);

        return fr * geometric * Le / lightSample.pdf;
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;

        if (!scene->rayIntersect(ray, its)) {
            return Color3f();
        } else {
            const Mesh &mesh = *(its.mesh);
            const BSDF &bsdf = *(mesh.getBSDF());

            if (bsdf.isDiffuse()) {
                // diffuse
                Color3f Le = Color3f(), Lr = Color3f();

                // add Le if mesh is emitter
                if (mesh.isEmitter()) {
                    Le = mesh.getEmitter()->getRadiance();
                }

                // integral over light sources
                for (int i = 0; i < SAMPLE_NUM; i++) {
                    Lr += sampleIntegralBody(scene, sampler, ray, its);
                }

                return Le + Lr / SAMPLE_NUM;
            } else {
                float zeta = sampler->next1D();

                if (zeta < 0.95) {
                    BSDFQueryRecord bsdfQR =
                        BSDFQueryRecord(its.shFrame.toLocal(-ray.d));
                    Color3f samplingWeight =
                        bsdf.sample(bsdfQR, sampler->next2D());

                    return (1 / 0.95) * samplingWeight *
                           Li(scene, sampler,
                              Ray3f(its.p, its.shFrame.toWorld(bsdfQR.wo)));
                } else {
                    return Color3f();
                }
            }
        }
    }

    std::string toString() const { return "WhittedIntegrator[]"; }

   private:
    DiscretePDF emitterDPDF;
};

NORI_REGISTER_CLASS(Whitted, "whitted");
NORI_NAMESPACE_END