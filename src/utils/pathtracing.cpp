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

#include <xercesc/parsers/SAXParser.hpp>
#include <mitsuba/core/appender.h>
#include <mitsuba/render/scenehandler.h>
#include <mitsuba/render/util.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/plugin.h>
#include <boost/algorithm/string.hpp>
#if defined(WIN32)
#include <mitsuba/core/getopt.h>
#endif
#include <mitsuba/render/integrator.h>

#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"

#include <thread> // to calculate the cpu core numbers

using XERCES_CPP_NAMESPACE::SAXParser;
using namespace mitsuba;

MTS_NAMESPACE_BEGIN

const int NUM_CHANNEL = 22;

class PathTracer : public Utility {
public:
	void help() {
		cout <<  "Usage: mtsutil pt [options] <XML file>" << endl;
		cout <<  "Options/Arguments:" << endl;
		cout <<  "   -h          Display this help text" << endl << endl;
		cout <<  "   -D key=val  Define a constant, which can referenced as \"$key\" in the scene" << endl << endl;
		cout <<  "   -o fname    Write the output image to the file denoted by \"fname\"" << endl << endl;
		cout <<  "   -p count    Override the detected number of processors. " << endl << endl;
		cout <<  "   -v          Be more verbose (can be specified twice)" << endl << endl;
		cout <<  "   -L level    Explicitly specify the log level (trace/debug/info/warn/error)" << endl << endl;
		cout <<  "   -q          Quiet mode. (log level = EError)" << endl << endl;
		cout <<  "   -i          Use the spp index as suffix of the output file" << endl << endl;
	}

	void writeOutput(std::vector<float> &m_image, std::string output, const Vector2i &resolution)
	{
		
		const char* channelNames[] = {
			"0_wix", "1_wiy", "2_wox", "3_woy", "4_R", "5_G", "6_B", 
			"7_wix", "8_wiy", "9_wox", "a_woy", "b_R", "c_G", "d_B", 
			"e_isN", "f_u",   "g_v",   "h_ind", "i_dudx", "j_dvdx", "k_dudy", "l_dvdy"
		};

		SaveEXRChannel(NUM_CHANNEL, channelNames, m_image, resolution[0], resolution[1], output.c_str());
		SLog(EInfo, "[PT] output saved into %s", output.c_str());
	}

	inline Float miWeight(Float pdfA, Float pdfB) const {
		pdfA *= pdfA;
		pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	void Li(const RayDifferential &r, RadianceQueryRecord &rRec,
		Vector &outWi_0, Vector &outWo_0, Spectrum &outLighting_0, 
		Vector &outWi_1, Vector &outWo_1, Spectrum &outLighting_1, 
		Float &isNeural, Float &u, Float &v, int &bsdfIndex,
		Float &dudx, Float &dvdx, Float &dudy, Float &dvdy
	) const 
	{
		/* Some aliases and local variables */
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		RayDifferential ray(r);
		bool scattered = false;

		/* Perform the first ray intersection (or ignore if the
		intersection has already been provided). */
		rRec.rayIntersect(ray);
		ray.mint = Epsilon;
		rRec.depth = 1;
		Spectrum throughput(1.0f);
		Float eta = 1.0f;
		bool m_hideEmitters = false;
		Properties pathProps = scene->getIntegrator()->getProperties();
		int m_maxDepth = pathProps.getInteger("maxDepth", 2);
		int m_rrDepth = pathProps.getInteger("rrDepth", 5);
		bool m_strictNormals = pathProps.getBoolean("strictNormals", false);
		int m_integrationType = pathProps.getInteger("integrationType", 2);

		outWi_0 = Vector(0.0f);
		outWo_0 = Vector(0.0f);
		outLighting_0 = Spectrum(0.0f);
		outWi_1 = Vector(0.0f);
		outWo_1 = Vector(0.0f);
		outLighting_1 = Spectrum(0.0f);
		while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {

			if (!its.isValid()) {
				/* If no intersection could be found, potentially return
				radiance from a environment luminaire if it exists */
				if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
					&& (!m_hideEmitters || scattered))
					{
						if (m_integrationType == 0 || m_integrationType == 2)
							outLighting_1 = scene->evalEnvironment(ray);
						else if(m_integrationType == 1)
							outLighting_0 = scene->evalEnvironment(ray);
						isNeural = -1.0f;
						outWi_0 = outWo_0 = outWi_1 = outWo_1 = Vector(0.0f);
					}
				break;
			}

			if ((rRec.depth >= m_maxDepth && m_maxDepth > 0) || (m_strictNormals && dot(ray.d, its.geoFrame.n) * Frame::cosTheta(its.wi) >= 0))
			{

				/* Only continue if:
				   1. The current path length is below the specifed maximum
				   2. If 'strictNormals'=true, when the geometric and shading
				      normals classify the incident direction to the same side */
				break;
			}

			//* [fjh] this operation should be after maxDepth check
			const BSDF *bsdf = its.getBSDF(ray);
			if (bsdf->isNeural())
			{
				isNeural = 1.0f;
				/* [fjh] get bsdf index, uv coordinates and save */
				bsdfIndex = bsdf->getIndex();
				u = its.uv.x;
				v = its.uv.y;
				dudx = its.dudx;
				dvdx = its.dvdx;
				dudy = its.dudy;
				dvdy = its.dvdy;
			}
			else
				isNeural = -1.0f;


			/* ==================================================================== */
			/*                     Direct illumination sampling                     */
			/* ==================================================================== */
			
			/* Estimate the direct illumination if this is requested */
			DirectSamplingRecord dRec(its);

			if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
				(bsdf->getType() & BSDF::ESmooth)) {
				Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
				if (!value.isZero()) {
					const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

					/* Allocate a record for querying the BSDF */
					BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);
					bRec.sampler = rRec.sampler;

					/* Evaluate BSDF * cos(theta) */
					Spectrum bsdfVal(1.0f);
					if (!bsdf->isNeural())
					{
						bsdfVal = bsdf->eval(bRec);
					}

					/* Prevent light leaks due to the use of shading normals */
					if (!bsdfVal.isZero() && (!m_strictNormals
						|| dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

						/* Calculate prob. of having generated that direction
                           using BSDF sampling */
                        Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                            ? bsdf->pdf(bRec) : 0;

                        Spectrum dirLighting(0.0f);

                        if (m_integrationType == 2 || m_integrationType == 4)
                            dirLighting = throughput * value * bsdfVal * miWeight(dRec.pdf, bsdfPdf);
                        else if (m_integrationType == 1)
                            dirLighting = throughput * value * bsdfVal;

						// [fjh] save wiwo & lighting of this part
						if (bsdf->isNormalMap()){
							const Intersection& _its = bRec.its;
							Intersection perturbed(_its);
							perturbed.shFrame = bsdf->getFrame(_its);
							outWi_0 = perturbed.toLocal(_its.toWorld(bRec.wi));
							outWo_0 = perturbed.toLocal(_its.toWorld(bRec.wo));
						}
						else{
							outWi_0 = bRec.wi;
							outWo_0 = bRec.wo;
						}
						outLighting_0 = dirLighting; //* [fjh] Li += dirLighting * eval(wi, wo);
					}
                }
            }

			/* ==================================================================== */
            /*                            BSDF sampling                             */
            /* ==================================================================== */

            /* Sample BSDF * cos(theta) */
            Float bsdfPdf;

            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
			//* [fjh] it should return eval / pdf, but only 1/pdf for now.
            // Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
            Spectrum invBsdfPdf = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
            // if (invBsdfPdf.isZero() || bRec.wo.z < 0.1)
            if (invBsdfPdf.isZero())
                break;

            scattered |= bRec.sampledType != BSDF::ENull;

            /* Prevent light leaks due to the use of shading normals */
            const Vector wo = its.toWorld(bRec.wo);
            Float woDotGeoN = dot(its.geoFrame.n, wo);
            if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                break;

            bool hitEmitter = false;
            Spectrum value(0.0);

            /* Trace a ray in this direction */
            ray = Ray(its.p, wo, ray.time);
            if (scene->rayIntersect(ray, its)) {
                /* Intersected something - check if it was a luminaire */
                if (its.isEmitter()) {
                    value = its.Le(-ray.d);
                    dRec.setQuery(ray, its);
                    hitEmitter = true;
				}
            } else {
				/* Intersected nothing -- perhaps there is an environment map? */
                const Emitter *env = scene->getEnvironmentEmitter();

                if (env) {
                    if (m_hideEmitters && !scattered)
                        break;

                    value = env->evalEnvironment(ray);
                    if (!env->fillDirectSamplingRecord(dRec, ray))
                        break;
                    hitEmitter = true;
                } else {
                    break;
                }
            }

            /* Keep track of the throughput and relative
               refractive index along the path */
            throughput *= invBsdfPdf; //*[fjh] here invBsdfPdf has no eval component
            eta *= bRec.eta;

			/* If a luminaire was hit, estimate the local illumination and
               weight using the power heuristic */
			Spectrum bsdfLighting(0.0f);
            if (hitEmitter &&
                (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
                /* Compute the prob. of generating that direction using the
                   implemented direct illumination sampling technique */
                const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
                    scene->pdfEmitterDirect(dRec) : 0;

				if (m_integrationType == 2 || m_integrationType == 3)
					bsdfLighting = throughput * value * miWeight(bsdfPdf, lumPdf);
				else if (m_integrationType == 0)
					bsdfLighting = throughput * value;

				// [fjh] save wiwo & lighting of this part
				if (bsdf->isNormalMap()){
					const Intersection& _its = bRec.its;
					Intersection perturbed(_its);
					perturbed.shFrame = bsdf->getFrame(_its);
					outWi_1 = perturbed.toLocal(_its.toWorld(bRec.wi));
					outWo_1 = perturbed.toLocal(_its.toWorld(bRec.wo));
				}
				else{
					outWi_1 = bRec.wi;
					outWo_1 = bRec.wo;
				}
				// for (int i = 0; i < 3; i++)
					// bsdfLighting[i] = bsdfLighting[i] < 32768.0f ? bsdfLighting[i] : 32768.0f;
				outLighting_1 = bsdfLighting; //* [fjh] Li += bsdfLighting * eval(wi, wo)
			}

            /* ==================================================================== */
            /*                         Indirect illumination                        */
            /* ==================================================================== */

            /* Set the recursive query type. Stop if no surface was hit by the
               BSDF sample or if indirect illumination was not requested */
            if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                break;
            rRec.type = RadianceQueryRecord::ERadianceNoEmission;

            if (rRec.depth++ >= m_rrDepth) {
                /* Russian roulette: try to keep path weights equal to one,
                   while accounting for the solid angle compression at refractive
                   index boundaries. Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
                if (rRec.nextSample1D() >= q)
                    break;
                throughput /= q;
            }
		}
	}

	double renderSingleSPP(Scene *scene,  const Point2 &offset, int idx_spp, std::string outPath)
	{
		Sampler *_sampler = scene->getSampler();

		size_t tcount = mts_omp_get_max_threads();
		std::vector<Sampler*> samplers(tcount);
		for (size_t i = 0; i < tcount; ++i) {
			ref<Sampler> clonedSampler = _sampler->clone();
			clonedSampler->incRef();
			samplers[i] = clonedSampler.get();
		}

		ref<Sensor> sensor = scene->getSensors()[0];
		const Vector2i resImage = sensor->getFilm()->getSize();

		Float diffScaleFactor = 1.0f /
			std::sqrt((Float)_sampler->getSampleCount());
			
		bool needsApertureSample = sensor->needsApertureSample();
		bool needsTimeSample = sensor->needsTimeSample();
			
		Point2 apertureSample(0.5f);
		Float timeSample = 0.5f;
		uint32_t queryType = RadianceQueryRecord::ESensorRay;

		// create a buffer for all the information
		// woTheta textures are enough, one for wi/wo, and the other for lighting
		std::vector<float> image(resImage.x * resImage.y * NUM_CHANNEL);
		std::vector<float> texture(resImage.x * resImage.y * 3);

		double t0 = omp_get_wtime();

#pragma omp parallel for
		for (int j = 0; j < resImage.y; j++)
		{
			for (int i = 0; i < resImage.x; i++)
			{
#if defined(MTS_OPENMP)
#warning("[pathtracing.cpp] MTS_OPENMP DEFINED.")
				int tid = mts_omp_get_thread_num();
#else
				int tid = 0;
#endif
				RadianceQueryRecord rRec(scene, samplers[tid]);
				Point2i offset(i, j);
				samplers[tid]->generate(offset);
				Point2 samplePos(Point2(offset) + Vector2(rRec.nextSample2D()));
				RayDifferential sensorRay;
				rRec.newQuery(queryType, sensor->getMedium());

				if (needsApertureSample)
					apertureSample = rRec.nextSample2D();
				if (needsTimeSample)
					timeSample = rRec.nextSample1D();
				Spectrum spec = sensor->sampleRayDifferential(
					sensorRay, samplePos, apertureSample, timeSample);
				sensorRay.scaleDifferential(diffScaleFactor);
				
				Vector wi_0(0.0f), wo_0(0.0f), wi_1(0.0f), wo_1(0.0f);
				Spectrum lighting_0, lighting_1;
				Float isNeural = -1.0f;
				Float u = 0.0f, v = 0.0f;
				int bsdfIndex = 0;
				Float dudx = 0.0f;
				Float dvdx = 0.0f;
				Float dudy = 0.0f;
				Float dvdy = 0.0f;
				Li(sensorRay, rRec,
					wi_0, wo_0, lighting_0,
					wi_1, wo_1, lighting_1,
					isNeural, u, v, bsdfIndex,
					dudx, dvdx, dudy, dvdy);
				lighting_0 *= spec;
				lighting_1 *= spec;
				if (wi_0.z < 0 || wo_0.z < 0)
				{
					lighting_0 = Spectrum(0.0f);
				}
				if (wi_1.z < 0 || wo_1.z < 0)
				{
					lighting_1 = Spectrum(0.0f);
				}

				auto it = image.begin() + (j * resImage.x + i) * NUM_CHANNEL;

				*(it++) = wi_0.x;
				*(it++) = wi_0.y;
				*(it++) = wo_0.x;
				*(it++) = wo_0.y;
				*(it++) = lighting_0[0];
				*(it++) = lighting_0[1];
				*(it++) = lighting_0[2];
				*(it++) = wi_1.x;
				*(it++) = wi_1.y;
				*(it++) = wo_1.x;
				*(it++) = wo_1.y;
				*(it++) = lighting_1[0];
				*(it++) = lighting_1[1];
				*(it++) = lighting_1[2];
				*(it++) = isNeural;
				*(it++) = u;
				*(it++) = v;
				*(it++) = bsdfIndex;
				*(it++) = dudx;
				*(it++) = dvdx;
				*(it++) = dudy;
				*(it++) = dvdy;
			}
		}

		double t1 = omp_get_wtime();
		double time = t1 - t0;

		writeOutput(image, outPath, resImage);
		std::vector<float>().swap(image);
		return time;
	}

	void pathtracing(Scene *scene, std::string path, bool useSppIndex)
	{
		ref<Sensor> sensor = scene->getSensors()[0];
		Sampler *sampler = scene->getSampler();
		const int SPP = sampler->getSampleCount();
	
		const int sqrtSPP = sqrt(SPP);
		double time = 0.0f;
		for (float p = 0; p < sqrtSPP; p++) {
			for (float q = 0; q < sqrtSPP; q++) {
				Point2 offset(
					Vector2(1.0 / sqrtSPP * (p + sampler->next1D()), 
					1.0 / sqrtSPP * (q + sampler->next1D()))
				);
				int idx = (int)p*sqrtSPP + (int)q;

				/* deal with output path */
				std::string outPath(path);
				std::string ext = ".exr";
				size_t found = path.find(ext);
				if (!useSppIndex){
					if (found != std::string::npos)
						outPath = outPath;
					else
						outPath = outPath.append(".exr");
				}
				else
				{
					if (found != std::string::npos)
						outPath = outPath.replace(found, ext.length(), std::string("_").append(std::to_string(idx)).append(".exr"));
					else
						outPath = outPath.append(std::string("_").append(std::to_string(idx)).append(".exr"));
				}
					
				time += renderSingleSPP(scene, offset, idx, outPath);
				fprintf(stderr, "[PT] CPU processing: %d/%d spp, time: %f\r", idx, sqrtSPP * sqrtSPP, (float)time);
			}
		}
		cout << endl;
	}

	int run(int argc, char **argv) {

		int optchar;
		char *end_ptr = NULL;

		/* Default settings */
		int nprocs_avail = getCoreCount(), nprocs = nprocs_avail;
		std::string destFile="";
		bool quietMode = false, useSppIndex = false;
		ELogLevel logLevel = EInfo;
		ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver();
		std::map<std::string, std::string, SimpleStringOrdering> parameters;

		if (argc < 2) {
			help();
			return 0;
		}

		optind = 1;
		/* Parse command-line arguments */
		while ((optchar = getopt(argc, argv, "D:o:L:p:vqih")) != -1)
		{
			switch (optchar) {
				case 'D': {
						std::vector<std::string> param = tokenize(optarg, "=");
						if (param.size() != 2)
							SLog(EError, "Invalid parameter specification \"%s\"", optarg);
						parameters[param[0]] = param[1];
					}
					break;
				case 'o':
					destFile = optarg;
					break;
                case 'v':
                    if (logLevel != EDebug)
                        logLevel = EDebug;
                    else
                        logLevel = ETrace;
                    break;
                case 'L': {
                        std::string arg = boost::to_lower_copy(std::string(optarg));
                        if (arg == "trace")
                            logLevel = ETrace;
                        else if (arg == "debug")
                            logLevel = EDebug;
                        else if (arg == "info")
                            logLevel = EInfo;
                        else if (arg == "warn")
                            logLevel = EWarn;
                        else if (arg == "error")
                            logLevel = EError;
                        else
                            SLog(EError, "Invalid log level!");
                    }
                    break;
				case 'p':
					nprocs = strtol(optarg, &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the processor count!");
					break;
                case 'q':
                    quietMode = true;
                    break;
				case 'i':
					useSppIndex = true;
					break;
				case 'h':
				default:
					help();
					return 0;
			}
		}

        /* Initialize OpenMP */
        Thread::initializeOpenMP(nprocs);
		
        /* Prepare for parsing scene descriptions */
        SAXParser* parser = new SAXParser();
        fs::path schemaPath = fileResolver->resolveAbsolute("data/schema/scene.xsd");

        /* Check against the 'scene.xsd' XML Schema */
        parser->setDoSchema(true);
        parser->setValidationSchemaFullChecking(true);
        parser->setValidationScheme(SAXParser::Val_Always);
        parser->setExternalNoNamespaceSchemaLocation(schemaPath.c_str());

        /* Set the handler */
        SceneHandler *handler = new SceneHandler(parameters);
        parser->setDoNamespaces(true);
        parser->setDocumentHandler(handler);
        parser->setErrorHandler(handler);

		std::string lowercase = boost::to_lower_copy(std::string(argv[optind]));
		fs::path
			filename = fileResolver->resolve(argv[optind]),
			filePath = fs::absolute(filename).parent_path(),
			baseName = filename.stem();
		ref<FileResolver> frClone = fileResolver->clone();
		frClone->prependPath(filePath);
		Thread::getThread()->setFileResolver(frClone);

		SLog(EInfo, "Parsing scene description from \"%s\" ..", argv[optind]);

		parser->parse(filename.c_str());
		ref<Scene> scene = handler->getScene();

        /* Configure the logging subsystem */
		Logger *logger = Thread::getThread()->getLogger();
		logger->setLogLevel(logLevel);
        /* Disable the default appenders */
        for (size_t i=0; i<logger->getAppenderCount(); ++i) {
            Appender *appender = logger->getAppender(i);
            if (appender->getClass()->derivesFrom(MTS_CLASS(StreamAppender)))
                logger->removeAppender(appender);
        }
        if (!quietMode)
            logger->addAppender(new StreamAppender(&std::cout));
			
		scene->initialize();		

		pathtracing(scene, (destFile.length() > 0 ? fs::path(destFile) : (filePath / baseName)).string(), useSppIndex);

		return 0;
	}

	MTS_DECLARE_UTILITY()

private:
	//may return 0 when not able to detect
	const uint cpu_count = std::thread::hardware_concurrency() == 0 ? 1: std::thread::hardware_concurrency();
	
};

MTS_EXPORT_UTILITY(PathTracer, "Individual path tracer.")
MTS_NAMESPACE_END
