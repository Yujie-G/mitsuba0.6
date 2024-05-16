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

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>

// [fjh] bsdfList
#include <omp.h>
#include <mitsuba/render/scene.h> // [fjh] to fix python symbol lookup error
#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{diffuse}{Smooth diffuse material}
 * \order{1}
 * \icon{bsdf_diffuse}
 * \parameters{
 *     \parameter{reflectance}{\Spectrum\Or\Texture}{
 *       Specifies the diffuse albedo of the
 *       material \default{0.5}
 *     }
 * }
 *
 * \renderings{
 *     \rendering{Homogeneous reflectance, see \lstref{diffuse-uniform}}
 *         {bsdf_diffuse_plain}
 *     \rendering{Textured reflectance, see \lstref{diffuse-textured}}
 *         {bsdf_diffuse_textured}
 * }
 *
 * The smooth diffuse material (also referred to as ``Lambertian'')
 * represents an ideally diffuse material with a user-specified amount of
 * reflectance. Any received illumination is scattered so that the surface
 * looks the same independently of the direction of observation.
 *
 * Apart from a  homogeneous reflectance value, the plugin can also accept
 * a nested or referenced texture map to be used as the source of reflectance
 * information, which is then mapped onto the shape based on its UV
 * parameterization. When no parameters are specified, the model uses the default
 * of 50% reflectance.
 *
 * Note that this material is one-sided---that is, observed from the
 * back side, it will be completely black. If this is undesirable,
 * consider using the \pluginref{twosided} BRDF adapter plugin.
 * \vspace{4mm}
 *
 * \begin{xml}[caption={A diffuse material, whose reflectance is specified
 *     as an sRGB color}, label=lst:diffuse-uniform]
 * <bsdf type="diffuse">
 *     <srgb name="reflectance" value="#6d7185"/>
 * </bsdf>
 * \end{xml}
 *
 * \begin{xml}[caption=A diffuse material with a texture map,
 *     label=lst:diffuse-textured]
 * <bsdf type="diffuse">
 *     <texture type="bitmap" name="reflectance">
 *         <string name="filename" value="wood.jpg"/>
 *     </texture>
 * </bsdf>
 * \end{xml}
 */
class SmoothDiffuse : public BSDF {
public:
    SmoothDiffuse(const Properties &props)
        : BSDF(props) {
        /* For better compatibility with other models, support both
           'reflectance' and 'diffuseReflectance' as parameter names */
        m_reflectance = new ConstantSpectrumTexture(props.getSpectrum(
            props.hasProperty("reflectance") ? "reflectance"
                : "diffuseReflectance", Spectrum(.5f)));
        m_reflectance_1 = new ConstantSpectrumTexture(props.getSpectrum(
            props.hasProperty("reflectance_1") ? "reflectance_1"
                : "diffuseReflectance_1", Spectrum(.5f)));
        m_decay_to = props.getFloat("decay_to", 0.0f);
    }

    SmoothDiffuse(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager) {
        m_reflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_reflectance_1 = static_cast<Texture *>(manager->getInstance(stream));

        configure();
    }

    void configure() {
        /* Verify the input parameter and fix them if necessary */
        m_reflectance = ensureEnergyConservation(m_reflectance, "reflectance", 1.0f);
        m_reflectance_1 = ensureEnergyConservation(m_reflectance_1, "reflectance_1", 1.0f);

        m_components.clear();
        if (m_reflectance->getMaximum().max() > 0)
            m_components.push_back(EDiffuseReflection | EFrontSide
                | (m_reflectance->isConstant() ? 0 : ESpatiallyVarying));
            m_usesRayDifferentials = m_reflectance->usesRayDifferentials();
        if (m_reflectance_1->getMaximum().max() > 0)
            m_components.push_back(EDiffuseReflection | EFrontSide
                | (m_reflectance_1->isConstant() ? 0 : ESpatiallyVarying));
            m_usesRayDifferentials = m_reflectance_1->usesRayDifferentials();

        BSDF::configure();
    }

    Spectrum getDiffuseReflectance(const Intersection &its, const Vector &wo) const {

        Spectrum eval0 = m_reflectance->eval(its);
        Spectrum eval1 = m_reflectance_1->eval(its);
        float weight = 1.0f - (1.0f - wo.z) * (1.0f - m_decay_to);
        if (weight < 0)
            weight = 0.0f;
        Spectrum evalRes = eval0 * weight + eval1 * (1.0f - weight);

        return evalRes;
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return Spectrum(0.0f);

        return getDiffuseReflectance(bRec.its, bRec.wo) * (INV_PI * Frame::cosTheta(bRec.wo));
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
        if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        return warp::squareToCosineHemispherePdf(bRec.wo);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
        if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        bRec.wo = warp::squareToCosineHemisphere(sample);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;
        return getDiffuseReflectance(bRec.its, bRec.wo);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
        if (!(bRec.typeMask & EDiffuseReflection) || Frame::cosTheta(bRec.wi) <= 0)
            return Spectrum(0.0f);

        bRec.wo = warp::squareToCosineHemisphere(sample);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;
        pdf = warp::squareToCosineHemispherePdf(bRec.wo);
        return getDiffuseReflectance(bRec.its, bRec.wo);
    }
    
	//* [fjh] bsdfList, returning a list of reflectance on top hemisphere *//
	float* bsdfList(int sampleCount, int length, int seed, int highspp=1) const {

		Assert(highspp == 1 && length == 9);
			
		//* [fjh] 9 for wix, wiy, wox, woy, u, v, R, G, B *//
		const int arraySize = sampleCount * length;
		float* result = new float[arraySize];

		//* [fjh] make sampler, generate sample *//
		Properties propos("independent");
		propos.setInteger("seed", seed);
		ref<Sampler> sampler = static_cast<Sampler*> (
			PluginManager::getInstance()->createObject(
				MTS_CLASS(Sampler), 
				propos
			)
		);
		sampler->generate(Point2i(0));

		//* [fjh] OMP stuffs *//
		int tcount = mts_omp_get_max_threads();
		std::vector<Sampler*> samplers(tcount);
		for (size_t i = 0; i < tcount; ++i) {
			ref<Sampler> clonedSampler = sampler->clone();
			clonedSampler->incRef();
			samplers[i] = clonedSampler.get();
		}

		//* [fjh] set to zero *//
#pragma omp parallel for
		for (int i = 0; i < arraySize; i++)
		{
			result[i] = 0.0f;
		}

#if defined(MTS_OPENMP)
	#pragma omp parallel for
#endif
			for (int currentCount = 0; currentCount < sampleCount; currentCount++){
#if defined(MTS_OPENMP)
				int tid = mts_omp_get_thread_num();
#else
				int tid = 0;
#endif

				Sampler* sampler = samplers[tid];
				sampler->setSampleIndex(currentCount);

				//* [fjh] generate wiwo*//
				float random1 = sampler->next1D();
				float random2 = sampler->next1D();
				float random3 = sampler->next1D();
				float random4 = sampler->next1D();
				float random5 = sampler->next1D();
				float random6 = sampler->next1D();

				float wiTheta = random1 * M_PI * 0.5; // [fjh] up to 75 deg
				float wiPhi = random2 * M_PI * 2.0;
				float woTheta = random3 * M_PI * 0.5;
				float woPhi = random4 * M_PI * 2.0;

				Vector wi(
						cos(wiPhi) * sin(wiTheta),
						sin(wiPhi) * sin(wiTheta),
						cos(wiTheta)
				);

				Vector wo(
					cos(woPhi) * sin(woTheta),
					sin(woPhi) * sin(woTheta),
					cos(woTheta)
				);

				//* [fjh] run eval *//
				Intersection its;
				its.uv = Point2(random5, random6);
				BSDFSamplingRecord bRec(its, wi, wo, ERadiance);
				bRec.typeMask = BSDF::EAll;
				ref<Sampler> bRecSampler = static_cast<Sampler*> (
					PluginManager::getInstance()->createObject(
						MTS_CLASS(Sampler), 
						propos
					)
				);
				Spectrum value = eval(bRec, ESolidAngle);

				//* [fjh] assign to result array *//
				int index = currentCount * length;
				result[index + 0] = bRec.wi.x; //* [fjh] keep the interval of [-1, 1] *//
				result[index + 1] = bRec.wi.y;
				result[index + 2] = bRec.wo.x;
				result[index + 3] = bRec.wo.y;
				result[index + 4] = its.uv.x;
				result[index + 5] = its.uv.y;

				result[index + 6] = value[0];
				result[index + 7] = value[1];
				result[index + 8] = value[2];
			} // end-for currentCount

			return result;
		}

	//* [fjh] free bsdfList *//
	void deletebsdfList(float* list)const
	{
		delete list;
		list = NULL;
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
        oss << "SmoothDiffuse[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  reflectance = " << indent(m_reflectance->toString()) << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    ref<Texture> m_reflectance;
    ref<Texture> m_reflectance_1;
    float m_decay_to;
};

// ================ Hardware shader implementation ================

class SmoothDiffuseShader : public Shader {
public:
    SmoothDiffuseShader(Renderer *renderer, const Texture *reflectance)
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

Shader *SmoothDiffuse::createShader(Renderer *renderer) const {
    return new SmoothDiffuseShader(renderer, m_reflectance.get());
}

MTS_IMPLEMENT_CLASS(SmoothDiffuseShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(SmoothDiffuse, false, BSDF)
MTS_EXPORT_PLUGIN(SmoothDiffuse, "Smooth diffuse BRDF")
MTS_NAMESPACE_END
