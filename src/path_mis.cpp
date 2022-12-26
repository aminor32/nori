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
        Color3f Lems(0.f), Lmats(0.f);
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
            float pEmitter = emitterDPDF[emitterIndex] || 1;

            // sample a point on the emitter
            Sample lightSample = emitter->sampleMesh(sampler);
            Vector3f lightDir = (lightSample.p - its.p).normalized();
            float d2 = (lightSample.p - its.p).squaredNorm();

            if (mesh.isEmitter()) {
                Lems += mesh.getEmitter()->getRadiance();
                Lmats += throughput * mesh.getEmitter()->getRadiance();
            }

            if (bsdf.isDiffuse()) {
                BSDFQueryRecord bsdfQR(its.toLocal(lightDir).normalized(),
                                       its.toLocal(-ray.d).normalized(),
                                       ESolidAngle);
                Color3f fr = bsdf.eval(bsdfQR);
                float pBRDF = bsdf.pdf(bsdfQR);

                float wLight = pEmitter / (pEmitter + pBRDF);

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

                // calculate Le and update L
                Color3f Le =
                    emitter->getEmitter()->Le(lightSample.n, -lightDir);
                Lems += wLight * throughput * fr * geometric * Le /
                        (lightSample.pdf * pEmitter);
            }

            // sample reflected direction
            BSDFQueryRecord bsdfQR(its.toLocal(-_ray.d).normalized());
            throughput *= bsdf.sample(bsdfQR, sampler->next2D());

            // next event estimation
            Intersection nextIts;
            Ray3f nextRay(its.p, its.toWorld(bsdfQR.wo).normalized());
            if (bsdf.isDiffuse() && scene->rayIntersect(nextRay, nextIts) &&
                nextIts.mesh->isEmitter()) {
                if (nextIts.mesh == emitter) {
                    break;
                }

                float pBRDF = bsdf.pdf(bsdfQR);
                float pLight;

                int i;
                for (i = 0; i < emitters.size(); i++) {
                    if (nextIts.mesh == emitters[i]) {
                        pLight = emitterDPDF[i];
                        break;
                    }
                }

                float wBRDF = pBRDF / (pBRDF + pLight);

                // calculate geometric term
                Color3f geometric =
                    std::abs(its.shFrame.cosTheta(bsdfQR.wo) *
                             nextIts.shFrame.cosTheta(
                                 nextIts.toLocal(-bsdfQR.wo).normalized())) /
                    (nextIts.t * nextIts.t);

                Color3f fr = bsdf.eval(bsdfQR);
                Color3f Le = emitter->getEmitter()->getRadiance();

                Lmats += wBRDF * fr * geometric * Le;
            }

            // update ray
            _ray = Ray3f(its.p, its.toWorld(bsdfQR.wo).normalized());

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

        return Lmats + Lems;
    }

    std::string toString() const { return "PathMisIntegrator[]"; }

   private:
    DiscretePDF emitterDPDF;
};

NORI_REGISTER_CLASS(PathMis, "path_mis");
NORI_NAMESPACE_END