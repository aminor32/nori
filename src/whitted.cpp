#include <nori/bsdf.h>
#include <nori/dpdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>

#include <random>
#include <vector>

#define SAMPLE_NUM 5

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
            Mesh *mesh = *it;

            if (mesh->isEmitter()) {
                emitters.push_back(mesh);
            }
        }
    }

    Color3f sampleIntegralBody(const Scene *scene, const Ray3f &ray,
                               Intersection &its) const {
        // select random emitter
        // TODO: select by pdf
        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_int_distribution<int> dist(0, emitters.size());
        int sample = dpdf.sample(dist(rng));
        const Mesh *emitter = emitters[sample];

        // sample from light source
        Sample lightSample = emitter->sampleMesh();
        Vector3f toLight = lightSample.p - its.p;  // not normalized !

        const BSDF &bsdf = *(its.mesh->getBSDF());
        BSDFQueryRecord bsdfQR =
            BSDFQueryRecord(its.toLocal(toLight).normalized(),
                            its.toLocal(-ray.d).normalized(), ESolidAngle);
        Color3f fr = bsdf.eval(bsdfQR);

        // calculate visibility term
        Ray3f shadowRay = Ray3f(its.p, toLight.normalized());
        shadowRay.maxt = toLight.norm();
        float visibility = scene->rayIntersect(shadowRay) ? 0.f : 1.f;

        // calculate geometric term
        Color3f geometric =
            visibility *
            std::abs(its.shFrame.cosTheta(its.toLocal(toLight)) *
                     lightSample.n.dot(toLight)) /
            toLight.squaredNorm();

        Color3f integralBody = lightSample.pdf * fr * geometric *
                               emitter->getEmitter()->Le(*emitter) /
                               toLight.norm();

        return integralBody;
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
                    Le = mesh.getEmitter()->Le(mesh);
                }

                // integral over light sources
                for (int i = 0; i < SAMPLE_NUM; i++) {
                    Lr += sampleIntegralBody(scene, ray, its);
                }

                return Le + Lr / SAMPLE_NUM;
            } else {
                // specular
                std::random_device rd;
                std::mt19937 rng(rd());
                std::uniform_real_distribution<float> dist(0, 1);

                BSDFQueryRecord bsdfQR = BSDFQueryRecord(-ray.d);
                Color3f result = Color3f();
                for (int i = 0; i < SAMPLE_NUM; i++) {
                    bsdf.sample(bsdfQR, Point2f(dist(rng), dist(rng)));

                    float zeta = dist(rng);

                    if (zeta < 0.95) {
                        result += (1 / 0.95) * (1 / SAMPLE_NUM) *
                                  Li(scene, sampler, Ray3f(its.p, bsdfQR.wo));
                    }
                }
            }
        }
    }

    std::string toString() const { return "WhittedIntegrator[]"; }

   private:
    std::vector<const Mesh *> emitters = {};
    DiscretePDF dpdf;
};

NORI_REGISTER_CLASS(Whitted, "whitted");
NORI_NAMESPACE_END