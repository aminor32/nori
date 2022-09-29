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
            return Color3f(0.5);
        } else {
            // normal at intersection point
            const Normal3f &n = its.shFrame.n;
            // position of intersection point
            const Point3f &x = its.p;
            // direction of light toward intersection point
            Vector3f dir = lightPosition - x;
            // square of distance between ray origin and intersection point
            float cosTheta = n.dot(dir) / (n.norm() * dir.norm());

            Ray3f shadowRay = Ray3f(lightPosition, dir.normalized());
            shadowRay.maxt = dir.norm();

            if (scene->rayIntersect(shadowRay)) {
                return Color3f();
            } else {
                Color3f l = (lightEnergy * INV_FOURPI / M_PI) *
                            (std::max<float>(0, cosTheta) / dir.squaredNorm());

                return l;
            }
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
            return Color3f(0.5);
        } else {
            // position of intersection point
            const Point3f &x = its.p;
            // number of samples per intersection
            int sampleNum = 5;
            // sum
            Color3f l = Color3f();

            std::random_device rd;
            std::mt19937 rng(rd());
            std::uniform_real_distribution<float> dist(0, 1);

            for (int i = 0; i < sampleNum; i++) {
                Vector3f sample = Warp::squareToCosineHemisphere(
                    Point2f(dist(rng), dist(rng)));
                Vector3f w = its.shFrame.toWorld(sample).normalized();
                float cosTheta = sample.z();
                Ray3f shadowRay = Ray3f(x, w);

                if (!scene->rayIntersect(shadowRay)) {
                    l += INV_PI * std::max<float>(0, cosTheta);
                }
            }

            return Color3f(l * M_PI * M_PI / sampleNum);
        }
    }

    std::string toString() const { return "AmbientOcclusionIntegrator[]"; }
};

NORI_REGISTER_CLASS(NormalIntegrator, "normals");
NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_REGISTER_CLASS(AOIntegrator, "ao");
NORI_NAMESPACE_END