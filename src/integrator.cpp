#include <nori/integrator.h>
#include <nori/scene.h>

#include <algorithm>
#include <cmath>

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
            Vector3f dir = x - lightPosition;
            // square of distance between ray origin and intersection point
            float cosTheta = n.dot((-1) * dir) / (n.norm() * dir.norm());
            float v = 0;

            Ray3f shadowRay =
                Ray3f(lightPosition, dir, 0,
                      dir.norm() - std::numeric_limits<float>::epsilon());
            Intersection shadowIts;

            if (!scene->getAccel()->rayIntersect(shadowRay, shadowIts, true)) {
                v = 1;
            }

            Color3f l = (lightEnergy / (4 * M_PI * M_PI)) *
                        (std::max<float>(0, cosTheta) / dir.squaredNorm()) * v;

            return l;
        }
    }

    std::string toString() const { return "SimpleIntegrator[]"; }

   private:
    Point3f lightPosition;
    Color3f lightEnergy;
};

NORI_REGISTER_CLASS(NormalIntegrator, "normals");
NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END