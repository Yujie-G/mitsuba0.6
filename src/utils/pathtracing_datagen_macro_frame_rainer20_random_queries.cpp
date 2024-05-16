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

#include "ubo2014_btf.h"

#include <thread> // to calculate the cpu core numbers

using XERCES_CPP_NAMESPACE::SAXParser;
using namespace mitsuba;

MTS_NAMESPACE_BEGIN

const int NUM_CHANNEL = 7;

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
	}

	void writeOutput(std::vector<float> &m_image, std::string output, const Vector2i &resolution)
	{
		
		const char* channelNames[] = {
			"0_R", "1_G", "2_B", "3_C", "4_D", "5_E", "6_F", 
		};

		SaveEXRChannel(NUM_CHANNEL, channelNames, m_image, resolution[0], resolution[1], output.c_str());
		SLog(EInfo, "[PT] output saved into %s", output.c_str());
	}

	double render(Scene *scene, std::string outPath)
	{
		ref<Sampler> _sampler = scene->getSensor()->getSampler();
		_sampler->generate(Point2i(0));
		Properties props = scene->getIntegrator()->getProperties();

		/* get a btf and use its directions */
		UBO2014_btf btf_model(props.getString("ubo_path").c_str(), 400, 400, 1, 1, 0, 0, "bilinear");
		const int maxQuery = btf_model.btf->LightCount * btf_model.btf->ViewCount;
		const int numQuery = props.getInteger("numQuery", 400);

		// [fjh] devide rX * rY coordinate space into resoX * resoY texels.
		int resoX = props.getInteger("resoX", 400); 
		int resoY = props.getInteger("resoY", 400);
		float rX = props.getFloat("rX", 5.f);
		float rY = props.getFloat("rY", 5.f);

		size_t tcount = mts_omp_get_max_threads();
		std::vector<Sampler*> samplers(tcount);
		for (size_t i = 0; i < tcount; ++i) {
			ref<Sampler> clonedSampler = _sampler->clone();
			clonedSampler->incRef();
			samplers[i] = clonedSampler.get();
		}

		ref<Sensor> sensor = scene->getSensors()[0];

		double t0 = omp_get_wtime();
		double time = 0.0f;

		#pragma omp parallel for
		for (int curTexel = 0; curTexel < resoX * resoY; curTexel++)
		{
			int uIndex = curTexel % resoY;
			int vIndex = curTexel / resoY;
			float u = uIndex * 1.0f / resoX; // (0 ~ 1)
			float v = vIndex * 1.0f / resoY;

			// create a buffer for all the information
			std::vector<float> image(numQuery * NUM_CHANNEL);
			int tid = mts_omp_get_thread_num();

			for (int i = 0; i < numQuery; i++)
			{
				// [fjh] choose a random query
				int queryIndex = static_cast<int>(samplers[tid]->next1D() * maxQuery);

				// [fjh] generate wi wo
				int viewIndex = queryIndex % btf_model.btf->LightCount;
				int lightIndex = queryIndex / btf_model.btf->LightCount;
				Vector wi(
					btf_model.btf->Views[viewIndex].x,
					btf_model.btf->Views[viewIndex].y,
					btf_model.btf->Views[viewIndex].z
				);
				Vector wo(
					btf_model.btf->Lights[lightIndex].x,
					btf_model.btf->Lights[lightIndex].y,
					btf_model.btf->Lights[lightIndex].z
				);
				float wiPhi = atan2(wi.y, wi.x);
				float woPhi = atan2(wo.y, wo.x);
				float wiTheta = acos(wi.z);
				float woTheta = acos(wo.z);

				// [fjh] make a ray
				RadianceQueryRecord rRec(scene, samplers[tid]);
				rRec.newQuery(RadianceQueryRecord::ESensorRay, sensor->getMedium());
				RayDifferential ray;
				ray.time = 0.0f;
				ray.setOrigin(Point((u * 2 - 1) * rX, (v * 2 - 1)  * rY, 0.0f)); // locate the origin coordinate
				ray.setDirection(wi);
				// SLog(EDebug, "o: %s, d: %s", ray.o.toString().c_str(), ray.d.toString().c_str());
				ray.mint = Epsilon;
				ray.maxt = 1000;
				ray.rxOrigin = ray.o + Point(rX * 2.0f / resoX, 0.0f, 0.0f);
				ray.ryOrigin = ray.o + Point(0.0f, rY * 2.0f / resoY, 0.0f);
				ray.rxDirection = ray.ryDirection = ray.d;
				ray.hasDifferentials = true;

				// [fjh] calculate
				rRec.rayIntersect(ray);
				Intersection &its = rRec.its;
				Spectrum bsdfVal(-1.f);

				if (!its.isValid()){
					wi.x = wi.y = wo.x = wo.y = -1.f;
				}
				else
				{
					const BSDF *bsdf = its.getBSDF(ray);
					/* 
						[fjh] here, wi and wo is actually facing INTO the screen, which is z+ direction.
						and the plane is facing OUT from the screen, which is z-
						so the actual wi and wo for $its should be -wi and -wo
					*/
					Vector shadingWi = its.toLocal(-wi);
					Vector shadingWo = its.toLocal(-wo);
					BSDFSamplingRecord bRec(its, shadingWi, shadingWo, ERadiance); // these vectors passed into bRec should be in shading frame
					bsdfVal = bsdf->eval(bRec);
					SLog(EDebug, "wi: %s, wo: %s -> wi: %s, wo: %s", wi.toString().c_str(), wo.toString().c_str(), bRec.wi.toString().c_str(), bRec.wo.toString().c_str());
				}
				
				// [fjh] assign values
				auto it = image.begin() + i * NUM_CHANNEL;
				*(it++) = bsdfVal[0];
				*(it++) = bsdfVal[1];
				*(it++) = bsdfVal[2];
				*(it++) = wiPhi * INV_PI * 180.0f + (wiPhi < 0 ? 360.0f : 0.0f);
				*(it++) = 90.0f - wiTheta * INV_PI * 180.0f;
				*(it++) = woPhi * INV_PI * 180.0f + (woPhi < 0 ? 360.0f : 0.0f);
				*(it++) = 90.0f - woTheta * INV_PI * 180.0f;
			}

			double t1 = omp_get_wtime();
			time += t1 - t0;

			char idx[12];
			snprintf(idx, 12, "-%06d.exr", curTexel);
			writeOutput(image, outPath + std::string(idx), Vector2i(numQuery, 1));
			std::vector<float>().swap(image);
		}
		return time;
	}

	int run(int argc, char **argv) {

		int optchar;
		char *end_ptr = NULL;

		/* Default settings */
		int nprocs_avail = getCoreCount(), nprocs = nprocs_avail;
		std::string destFile="";
		bool quietMode = false;
		ELogLevel logLevel = EInfo;
		ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver();
		std::map<std::string, std::string, SimpleStringOrdering> parameters;

		if (argc < 2) {
			help();
			return 0;
		}

		optind = 1;
		/* Parse command-line arguments */
		while ((optchar = getopt(argc, argv, "D:o:L:p:vqh")) != -1)
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
		
		/* deal with output path */
		std::string out_path = (destFile.length() > 0 ? fs::path(destFile) : (filePath / baseName)).string();
			
		render(scene, out_path);

		return 0;
	}

	MTS_DECLARE_UTILITY()

private:
	//may return 0 when not able to detect
	const uint cpu_count = std::thread::hardware_concurrency() == 0 ? 1: std::thread::hardware_concurrency();
	
};

MTS_EXPORT_UTILITY(PathTracerDataGen, "Individual path tracer.")
MTS_NAMESPACE_END
