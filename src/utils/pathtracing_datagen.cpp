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

const int NUM_CHANNEL = 9;

class PathTracerDataGen : public Utility {
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
			"0_wix", "1_wiy", "2_wox", "3_woy", "4_u", "5_v", "6_R", "7_G", "8_B", 
		};

		SaveEXRChannel(NUM_CHANNEL, channelNames, m_image, resolution[0], resolution[1], output.c_str());
		SLog(EInfo, "[PT] output saved into %s", output.c_str());
	}

	void Li(const RayDifferential &r, RadianceQueryRecord &rRec,
		Vector &wi, Vector &wo, Spectrum &bsdfVal,
		Sampler* wo_sampler
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

		while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {

			if (!its.isValid()) {
				/* If no intersection could be found, potentially return
				radiance from a environment luminaire if it exists */
				if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
					&& (!m_hideEmitters || scattered))
					{
						wi = wo = Vector(-1.0f);
						bsdfVal = Spectrum(-1.0f);
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

			/* ==================================================================== */
			/*                     Direct illumination sampling                     */
			/* ==================================================================== */
			
			if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
				(bsdf->getType() & BSDF::ESmooth)) {

				/* [fjh] random generate wo */
				float random1 = wo_sampler->next1D();
				float random2 = wo_sampler->next1D();

				float woTheta = random1 * M_PI * 0.5;// * 0.83; // [fjh] up to 75 deg
				float woPhi = random2 * M_PI * 2.0;

				wo.x = cos(woPhi) * sin(woTheta);
				wo.y = sin(woPhi) * sin(woTheta);
				wo.z = cos(woTheta);
				
				/* Allocate a record for querying the BSDF */
				// BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);
				//! [fjh] we assume this wo is already in shading frame
				BSDFSamplingRecord bRec(its, wo, ERadiance);

				//! [fjh] save wi and wo as a local frame vectors
				wi = bRec.wi;
				wo = bRec.wo;

				/* Evaluate BSDF * cos(theta) */
				bsdfVal = bsdf->eval(bRec);
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
		Properties propos("independent");
		ref<Sampler> _sampler = static_cast<Sampler*> (
			PluginManager::getInstance()->createObject(
				MTS_CLASS(Sampler), 
				propos
			)
		);
		_sampler->generate(Point2i(0));

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
				
				// [fjh] random generate wi
				float random1 = samplers[tid]->next1D();
				float random2 = samplers[tid]->next1D();

				float wiTheta = random1 * M_PI * 0.5 * 0.83; // [fjh] up to 75 deg
				float wiPhi = random2 * M_PI * 2.0;

				Point2 wiSample(
					cos(wiPhi) * sin(wiTheta),
					sin(wiPhi) * sin(wiTheta)
				);
				
				sensor->sampleRayDifferential(sensorRay, samplePos, wiSample, timeSample);
				SLog(EDebug, "o: %s, d: %s", sensorRay.o.toString().c_str(), sensorRay.d.toString().c_str());
				sensorRay.scaleDifferential(diffScaleFactor);
				
				Vector wi(
					wiSample.x,
					wiSample.y,
					sqrt(1.0f - wiSample.x*wiSample.x - wiSample.y*wiSample.y)
				);

				Vector wo(0.0f);
				Spectrum bsdfVal(0.0f);
				Li(sensorRay, rRec, wi, wo, bsdfVal, samplers[tid]);

				auto it = image.begin() + (j * resImage.x + i) * NUM_CHANNEL;

				*(it++) = wi.x;
				*(it++) = wi.y;
				*(it++) = wo.x;
				*(it++) = wo.y;
				*(it++) = samplePos.x / resImage.x;
				*(it++) = samplePos.y / resImage.y;
				*(it++) = bsdfVal[0];
				*(it++) = bsdfVal[1];
				*(it++) = bsdfVal[2];
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
				SLog(EInfo, "[PT] CPU processing: %d/%d spp, time: %f\r", idx, sqrtSPP * sqrtSPP, (float)time);
			}
		}
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

MTS_EXPORT_UTILITY(PathTracerDataGen, "Individual path tracer.")
MTS_NAMESPACE_END
