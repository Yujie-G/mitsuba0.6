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
#include <mitsuba/render/bsdf.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/spectrum.h>

#include <boost/algorithm/string.hpp>
#if defined(WIN32)
#include <mitsuba/core/getopt.h>
#endif
#include <mitsuba/render/integrator.h>

#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"

// #include "merl.h"

#include <thread> // to calculate the cpu core numbers

using XERCES_CPP_NAMESPACE::SAXParser;
using namespace mitsuba;

MTS_NAMESPACE_BEGIN

const int NUM_CHANNEL = 7;

class MatGen : public Utility {
public:
	void help() {
		cout << "Usage: mtsutil pt [options] <XML file>" << endl;
		cout << "Options/Arguments:" << endl;
		cout << "   -h          Display this help text" << endl << endl;
		cout << "   -D key=val  Define a constant, which can referenced as \"$key\" in the scene" << endl << endl;
		cout << "   -o fname    Write the output image to the file denoted by \"fname\"" << endl << endl;
		cout << "   -p count    Override the detected number of processors. " << endl << endl;
		cout << "   -v          Be more verbose (can be specified twice)" << endl << endl;
		cout << "   -L level    Explicitly specify the log level (trace/debug/info/warn/error)" << endl << endl;
		cout << "   -q          Quiet mode. (log level = EError)" << endl << endl;
		cout << "   -i          Use the spp index as suffix of the output file" << endl << endl;
	}

	void writeOutput(std::vector<float>& m_image, std::string output, const Vector2i& resolution)
	{

		const char* channelNames[] = {
			"0_view_x", "1_view_y", "2_light_x", "3_light_y", "4_R", "5_G", "6_B",
		};

		SaveEXRChannel(NUM_CHANNEL, channelNames, m_image, resolution[0], resolution[1], output.c_str());
		SLog(EInfo, "[MatGen] output saved into %s", output.c_str());
	}

	void parseCommandLine(int argc, char* argv[], std::map<std::string, std::string>& parameters, std::string& destFile, ELogLevel& logLevel, int& nprocs, bool& quietMode, bool& useSppIndex) {
		for (int i = 1; i < argc; ++i) {
			std::string arg = argv[i];
			if (arg[0] == '-') {
				switch (arg[1]) {
				case 'D': {
					if (i + 1 < argc) {
						std::string param = argv[++i];
						auto pos = param.find('=');
						if (pos != std::string::npos && pos != 0) {
							parameters[param.substr(0, pos)] = param.substr(pos + 1);
						}
						else {
							SLog(EError, "Invalid parameter specification \"%s\"", param.c_str());
						}
					}
					break;
				}
				case 'o':
					if (i + 1 < argc) destFile = argv[++i];
					break;
				case 'v':
					if (logLevel != EDebug)
						logLevel = EDebug;
					else
						logLevel = ETrace;
					break;
				case 'L':
					if (i + 1 < argc) {
						std::string arg = boost::to_lower_copy(std::string(argv[++i]));
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
					if (i + 1 < argc) {
						char* end_ptr = nullptr;
						nprocs = strtol(argv[++i], &end_ptr, 10);
						if (*end_ptr != '\0')
							SLog(EError, "Could not parse the processor count!");
					}
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
					return;
				}
			}
		}
	}

	// cross product of two vectors
	void cross_product(const double v1[3], const double v2[3], double(&out)[3]) {
		out[0] = v1[1] * v2[2] - v1[2] * v2[1];
		out[1] = v1[2] * v2[0] - v1[0] * v2[2];
		out[2] = v1[0] * v2[1] - v1[1] * v2[0];
	}

	// normalize vector
	void normalize(double v[3]) {
		// normalize
		double len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
		v[0] = v[0] / len;
		v[1] = v[1] / len;
		v[2] = v[2] / len;
	}

	// rotate vector along one axis
	void rotate_vector(const double vector[3], const double axis[3], double angle,
		double(&out)[3]) {
		double temp;
		double cross[3];
		double cos_ang = std::cos(angle);
		double sin_ang = std::sin(angle);

		out[0] = vector[0] * cos_ang;
		out[1] = vector[1] * cos_ang;
		out[2] = vector[2] * cos_ang;

		temp = axis[0] * vector[0] + axis[1] * vector[1] + axis[2] * vector[2];
		temp = temp * (1.0 - cos_ang);

		out[0] += axis[0] * temp;
		out[1] += axis[1] * temp;
		out[2] += axis[2] * temp;

		cross_product(axis, vector, cross);

		out[0] += cross[0] * sin_ang;
		out[1] += cross[1] * sin_ang;
		out[2] += cross[2] * sin_ang;
	}

	void half_angle_to_std(double theta_h, double theta_d, double phi_d,
		mitsuba::Vector3f& wi_v,
		mitsuba::Vector3f& wo_v) {
		double wi_z = std::cos(theta_d);
		double wi_xy = std::sin(theta_d);
		double wi_x = wi_xy * std::cos(phi_d);
		double wi_y = wi_xy * std::sin(phi_d);
		double wi[] = { wi_x, wi_y, wi_z };
		double wo[] = { -wi_x, -wi_y, wi_z };
		normalize(wi);
		normalize(wo);

		double bi_normal[] = { 0.0, 1.0, 0.0 };
		double normal[] = { 0.0, 0.0, 1.0 };
		double wi_o[3], wo_o[3];

		rotate_vector(wi, bi_normal, theta_h, wi_o);
		rotate_vector(wo, bi_normal, theta_h, wo_o);

		wi_v = Vector3f(wi_o[0], wi_o[1], wi_o[2]);
		wo_v = Vector3f(wo_o[0], wo_o[1], wo_o[2]);
	}

	void save_merl(const std::vector<double>& x, const std::string& path) {
		std::ofstream out_file(path, std::ios::out | std::ios::binary);

		int32_t merl_shape[] = { 90, 90, 180 };

		out_file.write(reinterpret_cast<const char*>(merl_shape),
			sizeof(merl_shape));

		out_file.write(reinterpret_cast<const char*>(x.data()),
			x.size() * sizeof(double));
	}

	void mat_gen(Scene *scene, std::string destFile){
        ref<BSDF> bsdf = scene->getShapes()[0]->getBSDF();
        std::vector<float> m_image;

		Intersection its;
		its.p = Point3(0.0, 0.0, 0.0);
		its.shFrame.n = Normal(0.0, 0.0, 1.0); 
		its.geoFrame.n = its.shFrame.n; 
		its.uv = Point2(0.5, 0.5); 
        Properties props("independent");
        props.setInteger("sampleCount", 16);
        // Instantiate the sampler
        ref<Sampler> sampler = static_cast<Sampler *>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));

		int cnt = 0;
		for (Float theta_i = 0; theta_i <= 90.0; theta_i += 90.0 / 25.0) {
			for (Float phi_i = 0; phi_i < 360.0; phi_i += 360.0 / 25.0) {
				for (Float theta_o = 0; theta_o <= 90.0; theta_o += 90.0 / 25.0) {
					for (Float phi_o = 0; phi_o < 360.0; phi_o += 360.0 / 25.0) {
						Vector wi = sphericalDirection(theta_i, phi_i);
						Vector wo = sphericalDirection(theta_o, phi_o);
                        BSDFSamplingRecord bRec(its, wi, wo);
                        bRec.sampler = sampler;
						Spectrum rgb = bsdf->eval(bRec);

						m_image.push_back(wi[0]);    // view_x
						m_image.push_back(wi[1]);     // view_y
						m_image.push_back(wo[0]);   // light_x
						m_image.push_back(wo[1]);     // light_y
						m_image.push_back(rgb[0]);     // R
						m_image.push_back(rgb[1]);     // G
						m_image.push_back(rgb[2]);     // B
						cnt++;
					}
				}
			}
		}

		std::vector<double> merl(3 * 90 * 90 * 180);
		for (size_t th = 0; th < 90; ++th) {
			double theta_h = static_cast<double>(th) / 90.0;
			theta_h = theta_h * theta_h * M_PI_2;
			for (size_t td = 0; td < 90; ++td) {
				double theta_d = static_cast<double>(td) / 90.0 * M_PI_2;
				for (size_t pd = 0; pd < 180; ++pd) {
					double phi_d = static_cast<double>(pd) / 180.0 * M_PI;

					Vector3f wi, wo;
					half_angle_to_std(theta_h, theta_d, phi_d, wi, wo);

					BSDFSamplingRecord bRec(its, wi, wo);
                    bRec.sampler = sampler;  
					Spectrum rgb = bsdf->eval(bRec);

					constexpr size_t k_stride = 90 * 90 * 180;
					const size_t idx = th * 90 * 180 + td * 180 + pd;

					merl[idx] = rgb[0] / wi[2] * 1500.0;
					merl[idx + k_stride] = rgb[1] / wi[2] * 1500.0 / 1.15;
					merl[idx + k_stride + k_stride] = rgb[2] / wi[2] * 1500.00 / 1.66;
				}
			}
		}

		

		/* deal with output path */
		std::string path = fs::path(destFile).string();
		std::string merl_outPath(path);
		merl_outPath = merl_outPath.append(".merl");
        std::cout<<"save file to "<< merl_outPath<<std::endl;
		save_merl(merl, merl_outPath);

        Vector2i resolution(25*25, 25*25);
		writeOutput(m_image, path.append(".exr"), resolution);
    }
	
	int run(int argc, char** argv) {

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
            std::cout << "DEBUG" <<std::endl;
			switch (optchar) {
				case 'D': {
						std::vector<std::string> param = tokenize(optarg, "=");
						if (param.size() != 2)
							SLog(EError, "Invalid parameter specification \"%s\"", optarg);
						parameters[param[0]] = param[1];
					}
					break;
				case 'o':
                    std::cout << "Option -o with argument: " << optarg << std::endl;
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

		mat_gen(scene, (destFile.length() > 0 ? fs::path(destFile) : (filePath / baseName)).string());
		return 0;
	}

	MTS_DECLARE_UTILITY()

private:
	//may return 0 when not able to detect
	//const uint cpu_count = std::thread::hardware_concurrency() == 0 ? 1 : std::thread::hardware_concurrency();

};

MTS_EXPORT_UTILITY(MatGen, "generate merl sample result")
MTS_NAMESPACE_END
