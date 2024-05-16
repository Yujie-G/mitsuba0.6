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

constexpr float invUINT32_MAX = 1.0f / std::numeric_limits<uint32_t>::max();
static uint32_t x = 114514, y = 1919810, z = 19260817, w = 0;

uint32_t xorshift128(void) {
    uint32_t t = x ^ (x << 11);
    x = y; y = z; z = w;
    return w = w ^ (w >> 19) ^ t ^ (t >> 8);
}

class GaussianDistr {
public:
    GaussianDistr(const Float sigma, const Float amp) : sigma(sigma), amp(amp){};

    inline Float eval(const Vector &m) const {
		return amp * exp(-(m.x * m.x + m.y * m.y) / (2.0f * sigma * sigma));
    }

    inline Float pdf(const Vector &m) const {
        return eval(m);
    }

    /* [fjh] generate Gaussian distribution by Box-Muller */
    inline Normal sample() const {
        Float tmp, x(0.5), y(0.5), r1(0.5), r2(0.5);
        do
        {
            r1 = xorshift128() * invUINT32_MAX;
            r2 = xorshift128() * invUINT32_MAX;
            tmp = std::sqrt(-2 * std::log(r1));
            x = tmp * std::cos(2 * Float(M_PI) * r2) * sigma;
            y = tmp * std::sin(2 * Float(M_PI) * r2) * sigma;
        } while(x * x + y * y >= 1);
        Normal H(x, y, sqrt(1 - x * x - y * y));
        // fprintf(stderr, "r1r2: (%f, %f), H: (%f, %f, %f)\n", r1, r2, x, y, H.z);
        return H;
    }

    Float sigma;
    Float amp;
};

class DecoderGaussian : public BSDF {
public:
    DecoderGaussian(const Properties &props) : BSDF(props) {
        ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();
        
        m_specularReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("specularReflectance", Spectrum(1.0f)));
        std::string filterType = props.getString("filterType", "nearest");
        if (filterType != "nearest")
            m_usesRayDifferentials = true;

        std::string materialName = props.getString("material", "Cu");

        Spectrum intEta, intK;
        if (boost::to_lower_copy(materialName) == "none") {
            intEta = Spectrum(0.0f);
            intK = Spectrum(1.0f);
        } else {
            intEta.fromContinuousSpectrum(InterpolatedSpectrum(
                fResolver->resolve("data/ior/" + materialName + ".eta.spd")));
            intK.fromContinuousSpectrum(InterpolatedSpectrum(
                fResolver->resolve("data/ior/" + materialName + ".k.spd")));
        }

        Float extEta = lookupIOR(props, "extEta", "air");

        m_eta = props.getSpectrum("eta", intEta) / extEta;
        m_k   = props.getSpectrum("k", intK) / extEta;

        MicrofacetDistribution distr(props);
        m_type = distr.getType();
        m_sampleVisible = distr.getSampleVisible();

        m_alphaU = new ConstantFloatTexture(distr.getAlphaU());
        if (distr.getAlphaU() == distr.getAlphaV())
            m_alphaV = m_alphaU;
        else
            m_alphaV = new ConstantFloatTexture(distr.getAlphaV());

        m_sigma = props.getFloat("sigma");
        m_amp = props.getFloat("amp");
        m_pUseDiffuse = props.getFloat("pUseDiffuse", 0.0f);
        m_bsdfIndex = props.getInteger("index", 0);
    }


    DecoderGaussian(Stream *stream, InstanceManager *manager)
     : BSDF(stream, manager) {
        m_type = (MicrofacetDistribution::EType) stream->readUInt();
        m_sampleVisible = stream->readBool();
        m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
        m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
        m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_eta = Spectrum(stream);
        m_k = Spectrum(stream);

        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        BSDF::serialize(stream, manager);

        stream->writeUInt((uint32_t) m_type);
        stream->writeBool(m_sampleVisible);
        manager->serialize(stream, m_alphaU.get());
        manager->serialize(stream, m_alphaV.get());
        manager->serialize(stream, m_specularReflectance.get());
        m_eta.serialize(stream);
        m_k.serialize(stream);
    }

    void configure() {
        unsigned int extraFlags = 0;
        if (m_alphaU != m_alphaV)
            extraFlags |= EAnisotropic;

        if (!m_alphaU->isConstant() || !m_alphaV->isConstant() ||
            !m_specularReflectance->isConstant())
            extraFlags |= ESpatiallyVarying;

        m_components.clear();
        m_components.push_back(EGlossyReflection | EFrontSide | extraFlags);

        /* Verify the input parameters and fix them if necessary */
        m_specularReflectance = ensureEnergyConservation(
            m_specularReflectance, "specularReflectance", 1.0f);

        m_usesRayDifferentials =  m_alphaU->usesRayDifferentials() ||
            m_alphaV->usesRayDifferentials() ||
            m_usesRayDifferentials; //! [fjh] here use the set value directly. *//

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

    void writeEXR(const std::vector<std::vector<float>> &map, std::string filename) const
    {
        std::vector<float> image;
		image.resize(map.size() * map[0].size() * 3, 0.0f);
		for (size_t i = 0; i < map.size(); i++)
		{
			for (size_t j = 0; j < map[0].size(); j++)
			{
				image[(i*map[0].size() + j) * 3 + 0] = map[i][j];
				image[(i*map[0].size() + j) * 3 + 1] = map[i][j];
				image[(i*map[0].size() + j) * 3 + 2] = map[i][j];
			}
		}
		SaveEXRChannel3(image, map.size(), map[0].size(), filename.c_str());
    }

    /// Helper function: reflect \c wi with respect to a given surface normal
    inline Vector reflect(const Vector &wi, const Normal &m) const {
        return 2 * dot(wi, m) * Vector(m) - wi;
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

        /* Calculate the reflection half-vector */
        Vector H = normalize(bRec.wo+bRec.wi);

        /* Construct the gaussian distribution matching the
           sigma values at the current surface position. */
        GaussianDistr distr(m_sigma, m_amp);

        return (1 - m_pUseDiffuse) * distr.pdf(H) / (4 * absDot(bRec.wo, H)) + 
            m_pUseDiffuse * warp::squareToCosineHemispherePdf(bRec.wo);
        // return distr.pdf(H);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        SLog(EError, "[Decoder_gaussian: sample] not implemented");
        return Spectrum(1.0f);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        float useDiffuse = xorshift128() * invUINT32_MAX;

        if (useDiffuse < m_pUseDiffuse)
        { /* [fjh] use diffuse sample */
            if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
                return Spectrum(0.0f);

            bRec.wo = warp::squareToCosineHemisphere(sample);
            bRec.eta = 1.0f;
            bRec.sampledComponent = 0;
            bRec.sampledType = EDiffuseReflection;
            //* [fjh] only use the direction sampling part of Diffuse *//
            //! [fjh] always, pdf = (1 - p) * G(w) + p * D(w), which is this->pdf(bRec) *//
            // pdf = warp::squareToCosineHemispherePdf(bRec.wo);
            // if (pdf < 1e-10f)
            //     return Spectrum(0.0f);
            // // return eval(bRec, ESolidAngle) / pdf;
            // return Spectrum(1.0f) / pdf;
        }
        else
        {
            if (Frame::cosTheta(bRec.wi) < 0 ||
                ((bRec.component != -1 && bRec.component != 0) ||
                !(bRec.typeMask & EGlossyReflection)))
                return Spectrum(0.0f);

            /* Construct the gaussian distribution matching the
            sigma values at the current surface position. */
            GaussianDistr distr(m_sigma, m_amp);

            /* Sample M, the microfacet normal */
            Normal m = distr.sample(); /*[fjh] here only get the G(w) part of pdf, not useful*/

            /* Perfect specular reflection based on the microfacet normal */
            bRec.wo = reflect(bRec.wi, m);
            bRec.eta = 1.0f;
            bRec.sampledComponent = 0;
            bRec.sampledType = EGlossyReflection;
        }

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
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
            if (name == "alpha")
                m_alphaU = m_alphaV = static_cast<Texture *>(child);
            else if (name == "alphaU")
                m_alphaU = static_cast<Texture *>(child);
            else if (name == "alphaV")
                m_alphaV = static_cast<Texture *>(child);
            else if (name == "specularReflectance")
                m_specularReflectance = static_cast<Texture *>(child);
            else
                BSDF::addChild(name, child);
        } else {
            BSDF::addChild(name, child);
        }
    }

    Float getRoughness(const Intersection &its, int component) const {
        return 0.5f * (m_alphaU->eval(its).average()
            + m_alphaV->eval(its).average());
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "DecoderGaussian[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  distribution = " << MicrofacetDistribution::distributionName(m_type) << "," << endl
            << "  sampleVisible = " << m_sampleVisible << "," << endl
            << "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
            << "  alphaV = " << indent(m_alphaV->toString()) << "," << endl
            << "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
            << "  eta = " << m_eta.toString() << "," << endl
            << "  k = " << m_k.toString() << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()

    
private:
    MicrofacetDistribution::EType m_type;
    ref<Texture> m_specularReflectance;
    ref<Texture> m_alphaU, m_alphaV;
    bool m_sampleVisible;
    Spectrum m_eta, m_k;

    //* [fjh] gaussian *//
    Float m_sigma;
    Float m_amp;
    float m_pUseDiffuse;

    int m_bsdfIndex;
};

/**
 * GLSL port of the DecoderGaussian shader. This version is much more
 * approximate -- it only supports the Ashikhmin-Shirley distribution,
 * does everything in RGB, and it uses the Schlick approximation to the
 * Fresnel reflectance of conductors. When the roughness is lower than
 * \alpha < 0.2, the shader clamps it to 0.2 so that it will still perform
 * reasonably well in a VPL-based preview.
 */
class DecoderGaussianShader : public Shader {
public:
    DecoderGaussianShader(Renderer *renderer, const Texture *specularReflectance,
            const Texture *alphaU, const Texture *alphaV, const Spectrum &eta,
            const Spectrum &k) : Shader(renderer, EBSDFShader),
            m_specularReflectance(specularReflectance), m_alphaU(alphaU), m_alphaV(alphaV) {
        m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
        m_alphaUShader = renderer->registerShaderForResource(m_alphaU.get());
        m_alphaVShader = renderer->registerShaderForResource(m_alphaV.get());

        /* Compute the reflectance at perpendicular incidence */
        m_R0 = fresnelConductorExact(1.0f, eta, k);
    }

    bool isComplete() const {
        return m_specularReflectanceShader.get() != NULL &&
               m_alphaUShader.get() != NULL &&
               m_alphaVShader.get() != NULL;
    }

    void putDependencies(std::vector<Shader *> &deps) {
        deps.push_back(m_specularReflectanceShader.get());
        deps.push_back(m_alphaUShader.get());
        deps.push_back(m_alphaVShader.get());
    }

    void cleanup(Renderer *renderer) {
        renderer->unregisterShaderForResource(m_specularReflectance.get());
        renderer->unregisterShaderForResource(m_alphaU.get());
        renderer->unregisterShaderForResource(m_alphaV.get());
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
        parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
        program->setParameter(parameterIDs[0], m_R0);
    }

    void generateCode(std::ostringstream &oss,
            const std::string &evalName,
            const std::vector<std::string> &depNames) const {
        oss << "uniform vec3 " << evalName << "_R0;" << endl
            << endl
            << "float " << evalName << "_D(vec3 m, float alphaU, float alphaV) {" << endl
            << "    float ct = cosTheta(m), ds = 1-ct*ct;" << endl
            << "    if (ds <= 0.0)" << endl
            << "        return 0.0f;" << endl
            << "    alphaU = 2 / (alphaU * alphaU) - 2;" << endl
            << "    alphaV = 2 / (alphaV * alphaV) - 2;" << endl
            << "    float exponent = (alphaU*m.x*m.x + alphaV*m.y*m.y)/ds;" << endl
            << "    return sqrt((alphaU+2) * (alphaV+2)) * 0.15915 * pow(ct, exponent);" << endl
            << "}" << endl
            << endl
            << "float " << evalName << "_G(vec3 m, vec3 wi, vec3 wo) {" << endl
            << "    if ((dot(wi, m) * cosTheta(wi)) <= 0 || " << endl
            << "        (dot(wo, m) * cosTheta(wo)) <= 0)" << endl
            << "        return 0.0;" << endl
            << "    float nDotM = cosTheta(m);" << endl
            << "    return min(1.0, min(" << endl
            << "        abs(2 * nDotM * cosTheta(wo) / dot(wo, m))," << endl
            << "        abs(2 * nDotM * cosTheta(wi) / dot(wi, m))));" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_schlick(float ct) {" << endl
            << "    float ctSqr = ct*ct, ct5 = ctSqr*ctSqr*ct;" << endl
            << "    return " << evalName << "_R0 + (vec3(1.0) - " << evalName << "_R0) * ct5;" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "   if (cosTheta(wi) <= 0 || cosTheta(wo) <= 0)" << endl
            << "        return vec3(0.0);" << endl
            << "   vec3 H = normalize(wi + wo);" << endl
            << "   vec3 reflectance = " << depNames[0] << "(uv);" << endl
            << "   float alphaU = max(0.2, " << depNames[1] << "(uv).r);" << endl
            << "   float alphaV = max(0.2, " << depNames[2] << "(uv).r);" << endl
            << "   float D = " << evalName << "_D(H, alphaU, alphaV)" << ";" << endl
            << "   float G = " << evalName << "_G(H, wi, wo);" << endl
            << "   vec3 F = " << evalName << "_schlick(1-dot(wi, H));" << endl
            << "   return reflectance * F * (D * G / (4*cosTheta(wi)));" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
            << "        return vec3(0.0);" << endl
            << "    return " << evalName << "_R0 * inv_pi * inv_pi * cosTheta(wo);"<< endl
            << "}" << endl;
    }
    MTS_DECLARE_CLASS()
private:
    ref<const Texture> m_specularReflectance;
    ref<const Texture> m_alphaU;
    ref<const Texture> m_alphaV;
    ref<Shader> m_specularReflectanceShader;
    ref<Shader> m_alphaUShader;
    ref<Shader> m_alphaVShader;
    Spectrum m_R0;
};

Shader *DecoderGaussian::createShader(Renderer *renderer) const {
    return new DecoderGaussianShader(renderer,
        m_specularReflectance.get(), m_alphaU.get(), m_alphaV.get(), m_eta, m_k);
}

MTS_IMPLEMENT_CLASS(DecoderGaussianShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(DecoderGaussian, false, BSDF)
MTS_EXPORT_PLUGIN(DecoderGaussian, "DecoderGaussian BRDF");
MTS_NAMESPACE_END
