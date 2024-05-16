#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include "MERL_BRDF.h"
#include <omp.h>
#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN

class MERL : public BSDF
{
public:
    MERL(const Properties &props)
        : BSDF(props)
    {
        m_importance = NULL;
        m_name = props.getString("binary");
    }

    MERL(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager)
    {
        m_name = stream->readString();
        m_importance = static_cast<BSDF *>(manager->getInstance(stream));
        m_importance->incRef();
        configure();
    }

    ~MERL() 
    {
        if (m_importance)
            m_importance->decRef();
    }

    void serialize(Stream *stream, InstanceManager *manager) const
    {
        BSDF::serialize(stream, manager);

        stream->writeString(m_name);
        manager->serialize(stream, m_importance);
    }

    void configure()
    {
        m_components.clear();
        m_components.push_back(EGlossyReflection | EFrontSide);
        m_components.push_back(EDiffuseReflection | EFrontSide);
        m_usesRayDifferentials = false;

        // check if importance sampling-guide has been specified
        if (m_importance == NULL)
            Log(EWarn, "[MERL] Importance sampling guidance not specified! Only evaluation supported.");

        // try to load the MERL BRDF data
        if (!read_brdf(m_name.c_str(), m_data))
            Log(EError, "[MERL] Unable to find \"%s\".", m_name.c_str());

        BSDF::configure();
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const
    {

        // handle diffuse and specular by weighting the evaluation proportional to the relative diffuse and specular component of the importance sampling
        bool hasSpecular = (bRec.typeMask & EGlossyReflection) && (bRec.component == -1 || bRec.component == 0);
        bool hasDiffuse = (bRec.typeMask & EDiffuseReflection) && (bRec.component == -1 || bRec.component == 1);

        // quick bail out
        if ((!hasDiffuse && !hasSpecular) || Frame::cosTheta(bRec.wo) <= 0.0f)
            return Spectrum(0.0f);

        // eval
        double r, g, b;
        double twi = acos(bRec.wi.z);
        double two = acos(bRec.wo.z);
        double pwi = atan2(bRec.wi.y, bRec.wi.x);
        double pwo = atan2(bRec.wo.y, bRec.wo.x);

        lookup_brdf_val(m_data, twi, pwi, two, pwo, r, g, b);

        Spectrum result(0.0f);
        result.fromLinearRGB(r, g, b);

        /* [fjh] I don't know what's going on here */

        // // subtract diffuse if requested
        // if (!hasDiffuse || !hasSpecular)
        // {
        //     BSDFSamplingRecord tRec = bRec;

        //     tRec.typeMask |= EDiffuseReflection;
        //     tRec.typeMask &= ~EGlossyReflection;
        //     Spectrum d = m_importance->eval(tRec, measure);

        //     // compute specular and diffuse component
        //     Spectrum s(0.0f);
        //     for (int i = 0; i < Spectrum::dim; i++)
        //     {
        //         Float diff = result[i] - d[i];
        //         if (diff > 0.0f)
        //             s[i] = diff;
        //         else // s[i] = 0.0f
        //             d[i] = result[i];
        //     }

        //     // copy
        //     if (hasDiffuse)
        //         result = d;
        //     else
        //         result = s;
        // }

        // Done
        return result * Frame::cosTheta(bRec.wo);
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const
    {
        if (m_importance == NULL)
            return 1.0f;
        return m_importance->pdf(bRec, measure);
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const
    {
        if (m_importance == NULL)
            return Spectrum(0.0f);
        m_importance->sample(bRec, pdf, sample);
        if (pdf == 0 || Frame::cosTheta(bRec.wo) <= 0)
            return Spectrum(0.0f);
        else
            return eval(bRec, ESolidAngle) / pdf;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const
    {
        if (m_importance == NULL)
            return Spectrum(0.0f);
        Float pdf;
        return MERL::sample(bRec, pdf, sample);
    }

    Float getRoughness(const Intersection &its, int component) const
    {
        return m_importance->getRoughness(its, component);
    }

    //* [fjh] bsdfList, returning a list of reflectance on top hemisphere *//
    float *bsdfList(int sampleCount, int length, int seed, int highspp = 1) const
    {

        Assert(highspp == 1 && length == 9);

        //* [fjh] 9 for wix, wiy, wox, woy, u, v, R, G, B *//
        const int arraySize = sampleCount * length;
        float *result = new float[arraySize];

        //* [fjh] make sampler, generate sample *//
        Properties propos("independent");
        propos.setInteger("seed", seed);
        ref<Sampler> sampler = static_cast<Sampler *>(
            PluginManager::getInstance()->createObject(
                MTS_CLASS(Sampler),
                propos));
        sampler->generate(Point2i(0));

        //* [fjh] OMP stuffs *//
        int tcount = mts_omp_get_max_threads();
        std::vector<Sampler *> samplers(tcount);
        for (size_t i = 0; i < tcount; ++i)
        {
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
        for (int currentCount = 0; currentCount < sampleCount; currentCount++)
        {
#if defined(MTS_OPENMP)
            int tid = mts_omp_get_thread_num();
#else
            int tid = 0;
#endif

            Sampler *sampler = samplers[tid];
            sampler->setSampleIndex(currentCount);

            //* [fjh] generate wiwo*//
            float random1 = sampler->next1D();
            float random2 = sampler->next1D();
            float random3 = sampler->next1D();
            float random4 = sampler->next1D();
            float random5 = sampler->next1D();
            float random6 = sampler->next1D();

            float wiTheta = random1 * M_PI * 0.5 * 0.83; // [fjh] up to 75 deg
            float wiPhi = random2 * M_PI * 2.0;
            float woTheta = random3 * M_PI * 0.5 * 0.83;
            float woPhi = random4 * M_PI * 2.0;

            Vector wi(
                cos(wiPhi) * sin(wiTheta),
                sin(wiPhi) * sin(wiTheta),
                cos(wiTheta));

            Vector wo(
                cos(woPhi) * sin(woTheta),
                sin(woPhi) * sin(woTheta),
                cos(woTheta));

            //* [fjh] run eval *//
            Intersection its;
            its.uv = Point2(random5, random6);
            BSDFSamplingRecord bRec(its, wi, wo, ERadiance);
            bRec.typeMask = BSDF::EAll;
            ref<Sampler> bRecSampler = static_cast<Sampler *>(
                PluginManager::getInstance()->createObject(
                    MTS_CLASS(Sampler),
                    propos));
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
    void deletebsdfList(float *list) const
    {
        delete list;
        list = NULL;
    }

    void addChild(const std::string &name, ConfigurableObject *child)
    {
        if (child->getClass()->derivesFrom(MTS_CLASS(BSDF)))
        {
            m_importance = static_cast<BSDF *>(child);
            m_importance->incRef();
        }
        else
        {
            BSDF::addChild(name, child);
        }
    }

    std::string toString() const
    {
        std::ostringstream oss;
        oss << "MERL[" << endl;
        oss << "id = \"" << getID() << "\"," << endl;
        oss << "binary = \"" << indent(m_name) << "\"," << endl;
        if (m_importance != NULL)
            oss << "importance = " << indent(m_importance->toString()) << endl;
        else
            oss << "importance = <not yet initialized>" << endl;
        oss << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
private:
    BSDF *m_importance;
    std::vector<double> m_data;
    std::string m_name;
};

// ================ Hardware shader implementation ================

/* MERL shader-- render as a 'black box' */
class MERLShader : public Shader
{
public:
    MERLShader(Renderer *renderer) : Shader(renderer, EBSDFShader)
    {
        m_flags = ETransparent;
    }

    void generateCode(std::ostringstream &oss,
                      const std::string &evalName,
                      const std::vector<std::string> &depNames) const
    {
        oss << "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return vec3(0.0);" << endl
            << "}" << endl;
        oss << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    return vec3(0.0);" << endl
            << "}" << endl;
    }
    MTS_DECLARE_CLASS()
};

Shader *MERL::createShader(Renderer *renderer) const
{
    return new MERLShader(renderer);
}

MTS_IMPLEMENT_CLASS(MERLShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(MERL, false, BSDF)
MTS_EXPORT_PLUGIN(MERL, "MERL BSDF");
MTS_NAMESPACE_END