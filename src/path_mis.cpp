#include <nori/bsdf.h>
#include <nori/dpdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <vector>

#define MAX_DEPTH 20

NORI_NAMESPACE_BEGIN

class PathMis : public Integrator {
   public:
    PathMis(const PropertyList &props) {
        // no arguments needed
    }

    void preprocess(const Scene *scene) {
        const std::vector<Mesh *> &emitters = scene->getEmitters();

        for (int i = 0; i < emitters.size(); i++) {
            emitterDPDF.append(1);
        }

        emitterDPDF.normalize();
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;
        Color3f Le(0.f), Ld(0.f), tmp(0.f);
        Color3f throughput(1.f);
        Ray3f _ray(ray);
        bool specular = false;

        for (int depth = 0; depth < MAX_DEPTH; depth++) {
            // terminate: no intersection between scene and ray
            if (!scene->rayIntersect(_ray, its)) {
                break;
            }

            const Mesh &mesh = *(its.mesh);
            const BSDF &bsdf = *(mesh.getBSDF());

            if (mesh.isEmitter()) {
                if (depth == 0) {
                    Le += mesh.getEmitter()->Le(its.shFrame.n, -_ray.d);
                } else if (specular) {
                    Ld += throughput *
                          mesh.getEmitter()->Le(its.shFrame.n, -_ray.d);
                }
            }

            // emitter sampling
            // uniformly sample an emiiter
            const std::vector<Mesh *> &emitters = scene->getEmitters();
            int emitterIndex = emitterDPDF.sample(sampler->next1D());
            const Mesh *emitter = emitters[emitterIndex];
            float pEmitter = 1.f / emitters.size();

            if (bsdf.isDiffuse()) {
                specular = false;

                // sample a point on the emitter
                Sample lightSample = emitter->sampleMesh(sampler);
                Vector3f lightDir = (lightSample.p - its.p).normalized();
                float d2 = (lightSample.p - its.p).squaredNorm();
                BSDFQueryRecord bsdfQR(its.toLocal(lightDir).normalized(),
                                       its.toLocal(-_ray.d).normalized(),
                                       ESolidAngle);

                // calculate w_light
                float pBRDF = bsdf.pdf(bsdfQR);
                float pLight = pEmitter * lightSample.pdf,
                      pLightSolidAngle =
                          pLight * d2 / std::abs(lightSample.n.dot(lightDir));
                float wLight = pLightSolidAngle / (pLightSolidAngle + pBRDF);

                if (std::isnan(wLight)) {
                    wLight = 1.f;
                }

                // calculate geometric term
                Ray3f shadowRay(its.p, lightDir);
                shadowRay.maxt = std::sqrt(d2) - Epsilon;
                int visibility = scene->rayIntersect(shadowRay) ? 0 : 1;
                Color3f geometric = visibility *
                                    std::abs(its.shFrame.n.dot(lightDir) *
                                             lightSample.n.dot(lightDir)) /
                                    d2;

                Color3f fr = bsdf.eval(bsdfQR);
                Color3f l = emitter->getEmitter()->Le(lightSample.n, -lightDir);

                Ld += wLight * throughput * fr * geometric * l / pLight;
            } else {
                specular = true;
            }

            // bsdf sampling
            // sample reflected direction & update ray
            BSDFQueryRecord bsdfQR(its.toLocal(-_ray.d).normalized());
            throughput *= bsdf.sample(bsdfQR, sampler->next2D());
            // Color3f samplingWeight = bsdf.sample(bsdfQR, sampler->next2D());
            _ray = Ray3f(its.p, its.toWorld(bsdfQR.wo).normalized());

            Intersection nextIts;
            if (scene->rayIntersect(_ray, nextIts) &&
                nextIts.mesh->isEmitter()) {
                const Mesh &nextMesh = *(nextIts.mesh);
                Vector3f lightDir = (nextIts.p - its.p).normalized();
                float d2 = (nextIts.p - its.p).squaredNorm();

                // calculate w_brdf
                float pBRDF = bsdf.pdf(bsdfQR);
                float pLight = pEmitter / nextMesh.getSurfaceArea(),
                      pLightSolidAngle =
                          pLight * d2 /
                          std::abs(nextIts.shFrame.n.dot(lightDir));
                float wBRDF = pBRDF / (pLightSolidAngle + pBRDF);

                if (std::isnan(wBRDF)) {
                    wBRDF = 1.f;
                }

                if (bsdf.isDiffuse()) {
                    Ld +=
                        wBRDF * throughput *
                        nextMesh.getEmitter()->Le(nextIts.shFrame.n, -lightDir);
                }
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

        return Le + Ld;
    }

    std::string toString() const { return "PathMisIntegrator[]"; }

   private:
    DiscretePDF emitterDPDF;
};

NORI_REGISTER_CLASS(PathMis, "path_mis");
NORI_NAMESPACE_END