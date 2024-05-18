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
#include "ubo2003_btf.h" 

#include <omp.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/bitmap.h>

MTS_NAMESPACE_BEGIN
    class Ubo2003 : public BSDF 
    {
    public:
        Ubo2003(const Properties &props)
            : BSDF(props), btf_model(
                               props.getString("ubo_path").c_str(),
                               props.getInteger("available_texels_u", 256), props.getInteger("available_texels_v", 256),
                               props.getInteger("u_interval", 1), props.getInteger("v_interval", 1),
                               props.getInteger("offset_texels_u", 0), props.getInteger("offset_texels_v", 0),
                               props.getString("interpType", "bilinear").c_str())
        {
            uvscale = props.getFloat("uvscale", 1.f);
            grayScale = props.getBoolean("gray_scale", false);
            singleChannel = props.getBoolean("single_channel", false);
            Assert((grayScale && singleChannel) == false && "grayScale and singleChannel cannot be both true");

            // available_texels_u = props.getInteger("available_texels_u", 400);
            // available_texels_v = props.getInteger("available_texels_v", 400);
            u_interval = props.getInteger("u_interval", 1); // [fjh] the interval for sparsely sampled BTF. 
            v_interval = props.getInteger("v_interval", 1);
            
            // [fjh] delaunay triangulation initialization
            if (props.getBoolean("useTriangulation", true))
            {
                int lightCount = btf_model.btf->LightCount;
                SLog(EDebug, "[ubo2003.cpp: Ubo2003()] lightCount: %d, ([%f, %f, %f; ...; %f, %f, %f]\n", lightCount,
                    btf_model.btf->Lights[0].x, btf_model.btf->Lights[0].y, btf_model.btf->Lights[0].z,
                    btf_model.btf->Lights[lightCount-1].x, btf_model.btf->Lights[lightCount-1].y, btf_model.btf->Lights[lightCount-1].z);
                std::vector<dt::Vector3D *> dots;
                for (int i = 0; i < lightCount; i++)
                {
                    auto dir = btf_model.btf->Lights[i];
                    dots.push_back(new dt::Vector3D(dir.x, dir.y, dir.z));
                }
                SLog(EDebug, "[ubo2003.cpp: Ubo2003()] dots.size: %d, ([%f, %f, %f; ...; %f, %f, %f])\n", (int)dots.size(), 
                        dots[0]->X, dots[0]->Y, dots[0]->Z, dots[dots.size()-1]->X, dots[dots.size()-1]->Y, dots[dots.size()-1]->Z);

                dt::DelaunayTriangulation triangulation;
                std::vector<std::tuple<int, int, int> *> mesh = triangulation.GetTriangulationResult(dots);
                SLog(EDebug, "[ubo2003.cpp: Ubo2003()] mesh.size: %d, ([%d, %d, %d; ...; %d, %d, %d])\n", (int)mesh.size(),
                    std::get<0>(*mesh[0]), std::get<1>(*mesh[0]), std::get<2>(*mesh[0]), 
                    std::get<0>(*mesh[mesh.size()-1]), std::get<1>(*mesh[mesh.size()-1]), std::get<2>(*mesh[mesh.size()-1]));

                SLog(EDebug, "[ubo2003.cpp: Ubo2003()] Triangulation statistics: %s\n", triangulation.GetStatistics().c_str());
                btf_model.delaunay_triangles = mesh;

                dt::Vector3D::resetRunningID();
                dt::Triangle::resetRunningID();
            }
        }

        Ubo2003(Stream *stream, InstanceManager *manager)
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
            // if (texel_index_u >= 0) //[fjh] allrandom version cannot specify texel indices
            // u = texel_index_u / 400.0f;
            // else
            u = (bRec.its.uv.x) * uvscale; // [fjh] multiply uvscale
            // if (texel_index_v >= 0)
            // v = texel_index_v / 400.0f;
            // else
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

        Float getRoughness(const Intersection &its, int component) const {
            return std::numeric_limits<Float>::infinity();
        }

        std::string toString() const {
            std::ostringstream oss;
            oss << "Ubo2003[" << endl
                << "  id = \"" << getID() << "\"," << endl
                << "  grayscale = " << grayScale << endl
                << "  singleChannel = " << singleChannel << endl
                << "  u_interval = \"" << u_interval << "\"," << endl
                << "  v_interval = \"" << v_interval << "\"," << endl
                << "]";
            return oss.str();
        }

        Shader *createShader(Renderer *renderer) const;

        MTS_DECLARE_CLASS()
        
    private:
        ref<Texture> m_reflectance;
        UBO2003_btf btf_model;
        bool grayScale;
        bool singleChannel;
        float uvscale;
        int u_interval, v_interval;
    };
// ================ Hardware shader implementation ================

    class Ubo2003Shader : public Shader {
    public:
        Ubo2003Shader(Renderer *renderer, const Texture *reflectance)
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

    Shader *Ubo2003::createShader(Renderer *renderer) const {
        return new Ubo2003Shader(renderer, m_reflectance.get());
    }

    MTS_IMPLEMENT_CLASS(Ubo2003Shader, false, Shader)
    MTS_IMPLEMENT_CLASS_S(Ubo2003, false, BSDF)
    MTS_EXPORT_PLUGIN(Ubo2003, "Ubo 2003 BRDF")

MTS_NAMESPACE_END