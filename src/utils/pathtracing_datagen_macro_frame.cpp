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
	}

	void writeOutput(std::vector<float> &m_image, std::string output, const Vector2i &resolution)
	{
		
		const char* channelNames[] = {
			"0_wix", "1_wiy", "2_wox", "3_woy", "4_u", "5_v", "6_R", "7_G", "8_B", 
		};

		SaveEXRChannel(NUM_CHANNEL, channelNames, m_image, resolution[0], resolution[1], output.c_str());
		SLog(EInfo, "[PT] output saved into %s", output.c_str());
	}

	double render(Scene *scene, std::string outPath)
	{
		ref<Sampler> _sampler = scene->getSensor()->getSampler();
		_sampler->generate(Point2i(0));
		Properties props = scene->getIntegrator()->getProperties();
		Vector inputWi = normalize(props.getVector("wi", Vector(0.f, 0.f, -1.f)));
		Vector inputWo = normalize(props.getVector("wo", Vector(0.f, 0.f, -1.f)));
		float inputU = props.getFloat("u", -1.f); // [fjh] -1 (default): random uv, -2: uniform uv, 0~1: specified uv
		float inputV = props.getFloat("v", -1.f);
		float rX = props.getFloat("rX", 5.f);
		float rY = props.getFloat("rY", 5.f);
		SLog(EInfo, "wi: %s, wo: %s, u: %f, v: %f", inputWi.toString().c_str(), inputWo.toString().c_str(), inputU, inputV);

		size_t tcount = mts_omp_get_max_threads();
		std::vector<Sampler*> samplers(tcount);
		for (size_t i = 0; i < tcount; ++i) {
			ref<Sampler> clonedSampler = _sampler->clone();
			clonedSampler->incRef();
			samplers[i] = clonedSampler.get();
		}

		ref<Sensor> sensor = scene->getSensors()[0];
		const Vector2i resImage = sensor->getFilm()->getSize();

		// create a buffer for all the information
		std::vector<float> image(resImage.x * resImage.y * NUM_CHANNEL);

		double t0 = omp_get_wtime();

		#pragma omp parallel for
		for (int count = 0; count < resImage.x * resImage.y; count++)
		{
			int tid = mts_omp_get_thread_num();

			// [fjh] random generate wi wo
			float random1 = samplers[tid]->next1D();
			float random2 = samplers[tid]->next1D();
			float random3 = samplers[tid]->next1D();
			float random4 = samplers[tid]->next1D();
			float random5 = samplers[tid]->next1D();
			float random6 = samplers[tid]->next1D();

			float wiTheta = random1 * M_PI * 0.5 * 0.83; // [fjh] up to 75 deg
			float wiPhi = random2 * M_PI * 2.0;
			float woTheta = random3 * M_PI * 0.5 * 0.83; // [fjh] up to 75 deg
			float woPhi = random4 * M_PI * 2.0;

			Vector wi(inputWi);
			Vector wo(inputWo);
			if (inputWi.z < 0)
			{
				wi.x = cos(wiPhi) * sin(wiTheta);
				wi.y = sin(wiPhi) * sin(wiTheta);
				wi.z = cos(wiTheta);
			}
			if (inputWo.z < 0)
			{
				wo.x = cos(woPhi) * sin(woTheta);
				wo.y = sin(woPhi) * sin(woTheta);
				wo.z = cos(woTheta);
			}

			// [fjh] find (fake) uv location in ref plane
			float u = 0.0f, v = 0.0f;
			Point2 samplePos(Point2(count % resImage.x, count / resImage.x) + Vector2(samplers[tid]->next2D())); // film space
			if (inputU == -2.f)
			{
				u = samplePos.x / resImage.x;
			}
			else if (inputU == -1.f)
			{
				u = random5;
			}
			else if (inputU >= 0 && inputU < 1)
			{
				u = inputU * resImage.x;
			}
			else
			{
				SLog(EError, "specified u must be in [0, 1), got %f", inputU);
			}
			if (inputV == -2.f) 
			{
				v = samplePos.y / resImage.y;
			}
			else if (inputV == -1.f)
			{
				v = random6;
			}
			else if (inputV >= 0 && inputV < 1)
			{
				v = inputV * resImage.x;
			}
			else
			{
				SLog(EError, "specified v must be in [0, 1), got %f", inputV);
			}

			// [fjh] make a ray
			RadianceQueryRecord rRec(scene, samplers[tid]);
			rRec.newQuery(RadianceQueryRecord::ESensorRay, sensor->getMedium());
			RayDifferential ray;
			ray.time = 0.0f;
			ray.setOrigin(Point(u * rX * 2 - rX, v * rY * 2 - rY, 0.0f)); // plane size: (-5, 5)
			ray.setDirection(wi);
			// SLog(EDebug, "o: %s, d: %s", ray.o.toString().c_str(), ray.d.toString().c_str());
			ray.mint = Epsilon;
			ray.maxt = 100;
			ray.rxOrigin = ray.o + Point(rX * 2 / resImage.x, 0.0f, 0.0f);
			ray.ryOrigin = ray.o + Point(0.0f, rY * 2 / resImage.y, 0.0f);
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
			auto it = image.begin() + count * NUM_CHANNEL;
			*(it++) = wi.x;
			*(it++) = wi.y;
			*(it++) = wo.x;
			*(it++) = wo.y;
			*(it++) = u; // fake uv
			*(it++) = v;
			*(it++) = bsdfVal[0];
			*(it++) = bsdfVal[1];
			*(it++) = bsdfVal[2];
		}

		double t1 = omp_get_wtime();
		double time = t1 - t0;

		writeOutput(image, outPath, resImage);
		std::vector<float>().swap(image);
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
		size_t found = out_path.find(std::string(".exr"));
		if (found == std::string::npos) // if not found
			out_path = out_path.append(".exr");

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
