#include <nori/integrator.h>
#include <nori/scene.h>

#include <cmath>

NORI_NAMESPACE_BEGIN

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
            // direction of light from intersection point
            Vector3f dir = lightPosition - x;
            // square of distance between ray origin and intersection point
            float cosTheta = n.dot(dir) / (n.norm() * dir.norm());

            Ray3f shadowRay = Ray3f(x, dir.normalized());
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

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END