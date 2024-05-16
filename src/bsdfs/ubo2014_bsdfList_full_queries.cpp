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
#include "ubo2014_btf.h"

// [fjh] bsdfList
#include <omp.h>
#include <mitsuba/render/scene.h> // [fjh] to fix python symbol lookup error
#include <mitsuba/core/plugin.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{diffuse}{Smooth diffuse material}
 */
    class Ubo2014 : public BSDF {
    public:
        Ubo2014(const Properties &props)
            : BSDF(props), btf_model(
                               props.getString("ubo_path").c_str(),
                               props.getInteger("available_texels_u", 400), props.getInteger("available_texels_v", 400),
                               props.getInteger("u_interval", 1), props.getInteger("v_interval", 1),
                               props.getInteger("offset_texels_u", 0), props.getInteger("offset_texels_v", 0),
                               props.getString("interpType", "bilinear").c_str())
        {

            uvscale = props.getFloat("uvscale", 1.f);
            grayScale = props.getBoolean("gray_scale", false);
            singleChannel = props.getBoolean("single_channel", false);
            Assert((grayScale && singleChannel) == false && "grayScale and singleChannel cannot be both true");

            available_texels_u = props.getInteger("available_texels_u", 400);
            available_texels_v = props.getInteger("available_texels_v", 400);
            u_interval = props.getInteger("u_interval", 1); // [fjh] the interval for sparsely sampled BTF. 
            v_interval = props.getInteger("v_interval", 1);

            // [fjh] this line will cause segmentation fault, because btf->Lights is an unique_ptr and therefore cannot be copied
            // please use List-initialization, or you can declare btf_model as a pointer to enable this style. 
            //-- btf_model = UBO2014_btf(props.getString("ubo_path").c_str(), props.getXXX(), ...);
            
            // [fjh] bsdfList
            texel_index_u = props.getInteger("texel_index_u", -1); // [fjh] these indices are real indices (i.e., from 0 to available_texels_u/v * u/v_interval - 1)
            texel_index_v = props.getInteger("texel_index_v", -1);
            Assert(texel_index_u < available_texels_u * u_interval && texel_index_v < available_texels_v * v_interval && "the specified texel index must in [0, UV)");

        }

        Ubo2014(Stream *stream, InstanceManager *manager)
                : BSDF(stream, manager) {


            configure();
        }

        void configure() {

            m_components.clear();

            m_components.push_back(EDiffuseReflection | EFrontSide   |  ESpatiallyVarying);
            m_usesRayDifferentials = true;

            BSDF::configure();
        }

        Spectrum getDiffuseReflectance(const Intersection &its) const {
            return Spectrum(0.f);
        }


        Spectrum eval_model(const BSDFSamplingRecord &bRec) const {

            float u, v;
            if (texel_index_u >= 0)
                u = texel_index_u * 1.0f / (available_texels_u * u_interval);
            else
                u = (bRec.its.uv.x) * uvscale; // [fjh] multiply uvscale
            if (texel_index_v >= 0)
                v = texel_index_v * 1.0f / (available_texels_v * v_interval);
            else
                v = (bRec.its.uv.y) * uvscale; // [fjh] multiply uvscale

            auto result =  btf_model.get_value(u, v, bRec.wo.x, bRec.wo.y, bRec.wo.z, bRec.wi.x, bRec.wi.y, bRec.wi.z);

            Float vals[] = {result.x, result.y, result.z};

            if (grayScale)
            {
                vals[0] = vals[1] = vals[2] = (result.x + result.y + result.z) / 3.0;
            }
            if (singleChannel)
            {
                vals[0] = vals[1] = vals[2] = result.z;
            }

            return Spectrum(vals);
        }

        Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const override {
            if (!(bRec.typeMask & EDiffuseReflection) || measure != ESolidAngle
                || Frame::cosTheta(bRec.wi) <= 0
                || Frame::cosTheta(bRec.wo) <= 0)
                return Spectrum(0.0f);

            return eval_model(bRec) * (INV_PI * Frame::cosTheta(bRec.wo));
        }

        Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
            return 1.0f;
        }

        Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
            return Spectrum(0.0f);
        }

        Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
            return Spectrum(0.0f);
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


        //* [fjh] bsdfList, returning a list of reflectance on top hemisphere *//
        float* bsdfList(int sampleCount, int length, int seed, int highspp=1) const {

            Assert(
                highspp == 1 && texel_index_u >= 0 && texel_index_u < (available_texels_u * u_interval) 
                             && texel_index_v >= 0 && texel_index_v < (available_texels_v * v_interval)
            );
            Assert(sampleCount == 22801 && length == 7);
                
            //* [fjh] 7 for wix, wiy, wox, woy, R, G, B *//
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
            for (int i = 0; i < tcount; ++i) {
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
                    int viewIndex = currentCount / btf_model.btf->LightCount;
                    int lightIndex = currentCount % btf_model.btf->LightCount;
                    auto view = btf_model.btf->Views[viewIndex];
                    auto light = btf_model.btf->Lights[lightIndex];

                    Vector wi(view.x, view.y, view.z);
                    Vector wo(light.x, light.y, light.z);

                    //* [fjh] run eval *//
                    Intersection its;
                    its.uv = Point2(texel_index_u * 1.0f / (available_texels_u * u_interval), texel_index_v * 1.0f / (available_texels_v * v_interval));
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

                    result[index + 4] = value[0];
                    result[index + 5] = value[1];
                    result[index + 6] = value[2];
                } // end-for currentCount

                return result;
            }

        //* [fjh] free bsdfList *//
        void deletebsdfList(float* list)const
        {
            delete list;
            list = NULL;
        }


        Float getRoughness(const Intersection &its, int component) const {
            return std::numeric_limits<Float>::infinity();
        }

        std::string toString() const {
            std::ostringstream oss;
            oss << "Ubo2014[" << endl
                << "  id = \"" << getID() << "\"," << endl
                << "  grayScale = " << grayScale << endl
                << "  singleChannel = " << singleChannel << endl
                << "  available_texels_u = " << available_texels_u  << endl
                << "  available_texels_v = " << available_texels_v << endl
                << "  u_interval = " << u_interval << endl
                << "  v_interval = " << v_interval << endl
                << "  texel_index_u = " << texel_index_u << endl
                << "  texel_index_v = " << texel_index_v << endl
                << "]";
            return oss.str();
        }

        Shader *createShader(Renderer *renderer) const;

        MTS_DECLARE_CLASS()
    private:
        ref<Texture> m_reflectance;
        UBO2014_btf btf_model;
        bool grayScale;
        bool singleChannel;
        float uvscale;
        int available_texels_u, available_texels_v;
        int u_interval, v_interval;

        // [fjh] bsdfList
        int texel_index_u, texel_index_v;
    };

// ================ Hardware shader implementation ================

    class Ubo2014Shader : public Shader {
    public:
        Ubo2014Shader(Renderer *renderer, const Texture *reflectance)
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

    Shader *Ubo2014::createShader(Renderer *renderer) const {
        return new Ubo2014Shader(renderer, m_reflectance.get());
    }

    MTS_IMPLEMENT_CLASS(Ubo2014Shader, false, Shader)
    MTS_IMPLEMENT_CLASS_S(Ubo2014, false, BSDF)
    MTS_EXPORT_PLUGIN(Ubo2014, "Ubo 2014 BRDF")
MTS_NAMESPACE_END
