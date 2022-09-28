#include <nori/frame.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

#include <algorithm>
#include <cmath>
#include <random>

NORI_NAMESPACE_BEGIN

class NormalIntegrator : public Integrator {
   public:
    NormalIntegrator(const PropertyList &props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its)) return Color3f(0.0f);

        /* Return the component-wise absolute
           value of the shading normal as a color */
        Normal3f n = its.shFrame.n.cwiseAbs();
        return Color3f(n.x(), n.y(), n.z());
    }

    std::string toString() const { return "NormalIntegrator[]"; }
};

class SimpleIntegrator : public Integrator {
   public:
    SimpleIntegrator(const PropertyList &props) {
        lightPosition = props.getPoint("position");
        lightEnergy = props.getColor("energy");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;

        if (!scene->rayIntersect(ray, its)) {
            return Color3f();
        } else {
            // normal at intersection point
            const Normal3f &n = its.shFrame.n;
            // position of intersection point
            const Point3f &x = its.p;
            // direction of light toward intersection point
            Vector3f dir = lightPosition - x;
            // square of distance between ray origin and intersection point
            float cosTheta = n.dot(dir) / (n.norm() * dir.norm());

            Ray3f shadowRay =
                Ray3f(lightPosition, dir.normalized(),
                      std::numeric_limits<float>::epsilon(), dir.norm());
            Intersection shadowIts;
            int v = scene->getAccel()->rayIntersect(shadowRay, shadowIts, true)
                        ? 0
                        : 1;

            if (v == 0) {
                return Color3f(0, 0, 1);
            }

            Color3f l = (lightEnergy * INV_FOURPI / M_PI) *
                        (std::max<float>(0, cosTheta) / dir.squaredNorm()) * v;

            return l;
        }
    }

    std::string toString() const { return "SimpleIntegrator[]"; }

   private:
    Point3f lightPosition;
    Color3f lightEnergy;
};

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

            // set number of samples
            int sampleNum = 20;
            float l = 0;

            for (int i = 0; i < sampleNum; i++) {
                Vector3f sample =
                    Warp::squareToCosineHemisphere(sampler->next2D());
                Vector3f w = its.toWorld(sample);
                float cosTheta = its.shFrame.cosTheta(w - x);

                Ray3f shadowRay = Ray3f(x, w);
                Intersection shadowIts;

                if (!scene->getAccel()->rayIntersect(shadowRay, shadowIts,
                                                     true)) {
                    l += INV_PI * std::max<float>(0, cosTheta);
                }
            }

            std::cout << l << std::endl;

            return Color3f(l);
        }
    }

    std::string toString() const { return "AmbientOcclusionIntegrator[]"; }
};

NORI_REGISTER_CLASS(NormalIntegrator, "normals");
NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_REGISTER_CLASS(AOIntegrator, "ao");
NORI_NAMESPACE_END