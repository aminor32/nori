#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
   public:
    AreaLight(const PropertyList &props) {
        radiance = props.getColor("radiance");
    }

    Color3f getRadiance(Point3f *lightSample) const { return radiance; }

    std::string toString() const { return "AreaLight[]"; }

   private:
    // radiance of the light source
    Color3f radiance;
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END