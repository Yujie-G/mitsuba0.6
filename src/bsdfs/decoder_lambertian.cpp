/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include "microfacet.h"
#include "ior.h"

#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"
#include <mitsuba/core/warp.h>

MTS_NAMESPACE_BEGIN

class DecoderLambertian : public BSDF {
public:
    DecoderLambertian(const Properties &props)
        : BSDF(props) {
        /* For better compatibility with other models, support both
           'reflectance' and 'diffuseReflectance' as parameter names */
        m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
            props.hasProperty("reflectance") ? "reflectance"
                : "diffuseReflectance", Spectrum(.5f)));
        m_bsdfIndex = props.getInteger("index", 0);
    }

    DecoderLambertian(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {
        m_reflectance = static_cast<Texture *>(manager->getInstance(stream));

        configure();
    }

    void configure() {
        /* Verify the input parameter and fix them if necessary */
        m_reflectance = ensureEnergyConservation(m_reflectance, "reflectance", 1.0f);

        m_components.clear();
        if (m_reflectance->getMaximum().max() > 0)
            m_components.push_back(EDiffuseReflection | EFrontSide
                | (m_reflectance->isConstant() ? 0 : ESpatiallyVarying));
        m_usesRayDifferentials = m_reflectance->usesRayDifferentials();

        BSDF::configure();
    }

	//* [fjh] return whether it is a neural BSDF, false as default *//
	inline bool isNeural() const {
		return true;
	}

    //* [fjh] return index of current bsdf at current intersection, 0 as default *//
    inline int getIndex() const {
        return m_bsdfIndex;
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        return Spectrum(1.0f);
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (measure != ESolidAngle ||
            Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 ||
            ((bRec.component != -1 && bRec.component != 0) ||
             !(bRec.typeMask & EGlossyReflection)))
            return 0.0f;

        return warp::squareToCosineHemispherePdf(bRec.wo);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        SLog(EError, "[DecoderLambertian: sample] not implemented");
        return Spectrum(1.0f);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        
        if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        bRec.wo = warp::squareToCosineHemisphere(sample);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;
        //* [fjh] only use the direction sampling part of Diffuse *//
        //! [fjh] always, pdf = D(w), which is this->pdf(bRec) *//
        // pdf = warp::squareToCosineHemispherePdf(bRec.wo);
        // if (pdf < 1e-10f)
        //     return Spectrum(0.0f);
        // // return eval(bRec, ESolidAngle) / pdf;
        // return Spectrum(1.0f) / pdf;

        /* Side check */
        if (Frame::cosTheta(bRec.wo) <= 0)
            return Spectrum(0.0f);

        /* Jacobian of the half-direction mapping */
        // pdf /= 4.0f * dot(bRec.wo, m);

        pdf = this->pdf(bRec, ESolidAngle);

        if (pdf < 1e-10f)
            return Spectrum(0.0f);

        return Spectrum(1.0f) / pdf;
    }

    void addChild(const std::string &name, ConfigurableObject *child) {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture))
                && (name == "reflectance" || name == "diffuseReflectance")) {
            m_reflectance = static_cast<Texture *>(child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        manager->serialize(stream, m_reflectance.get());
    }

    Float getRoughness(const Intersection &its, int component) const {
        return std::numeric_limits<Float>::infinity();
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "DecoderLambertian[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  reflectance = " << indent(m_reflectance->toString()) << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()

private:
    ref<Texture> m_reflectance;
    int m_bsdfIndex;
};

// ================ Hardware shader implementation ================

class DecoderLambertianShader : public Shader {
public:
    DecoderLambertianShader(Renderer *renderer, const Texture *reflectance)
        : Shader(renderer, EBSDFShader), m_reflectance(reflectance) {
        m_reflectanceShader = renderer->registerShaderForResource(m_reflectance.get());
    }

    bool isComplete() const {
        return m_reflectanceShader.get() != NULL;
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_reflectance.get());
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_reflectanceShader.get());
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
            << "        return vec3(0.0);" << endl
            << "    return " << depNames[0] << "(uv) * inv_pi * cosTheta(wo);" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return " << evalName << "(uv, wi, wo);" << endl
            << "}" << endl;
    }

    MTS_DECLARE_CLASS()
private:
    ref<const Texture> m_reflectance;
    ref<Shader> m_reflectanceShader;
};

Shader *DecoderLambertian::createShader(Renderer *renderer) const {
    return new DecoderLambertianShader(renderer, m_reflectance.get());
}

MTS_IMPLEMENT_CLASS(DecoderLambertianShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(DecoderLambertian, false, BSDF)
MTS_EXPORT_PLUGIN(DecoderLambertian, "DecoderLambertian BRDF");
MTS_NAMESPACE_END
