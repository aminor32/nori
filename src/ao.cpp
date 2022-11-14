#include <nori/frame.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

// Ambient Occlusion Integrator
class AOIntegrator : public Integrator {
   public:
    AOIntegrator(const PropertyList &props) {
        // no properties needed
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;

        if (!scene->rayIntersect(ray, its)) {
            return Color3f();
        } else {
            // position of intersection point
            const Point3f &x = its.p;
            // number of samples per intersection
            int sampleNum = 10;
            // sum
            Color3f l = Color3f();

            for (int i = 0; i < sampleNum; i++) {
                Vector3f sample =
                    Warp::squareToCosineHemisphere(sampler->next2D());
                Vector3f w = its.toWorld(sample);
                float cosTheta = sample.z();
                Ray3f shadowRay = Ray3f(x, w);

                if (!scene->rayIntersect(shadowRay)) {
                    l += INV_PI * cosTheta /
                         Warp::squareToCosineHemispherePdf(sample);
                }
            }

            return Color3f(l / sampleNum);
        }
    }

    std::string toString() const { return "AmbientOcclusionIntegrator[]"; }
};

NORI_REGISTER_CLASS(AOIntegrator, "ao");
NORI_NAMESPACE_END
