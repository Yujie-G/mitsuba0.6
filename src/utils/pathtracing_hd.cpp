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

MTS_NAMESPACE_BEGIN

const int NUM_CHANNEL = 22;

class PathTracer : public Utility {
public:
	void help() {
		cout << "Sorry, but you won't get helped." << endl;
		// cout << endl;
		// cout << "Synopsis: Muliple Bounce Glint Computation" << endl;
		// cout << endl;
		// cout << "Usage: mtsutil multibounceglint [options] <Scene XML file or PLY file>" << endl;
		// cout << "Options/Arguments:" << endl;
		// cout << "   -h             Display this help text" << endl << endl;
		// cout << "   -n value       Specify the bounce count" << endl << endl;
		// cout << "   -m true/false  Use hierarchy for pruning or not" << endl << endl;
	}

	void writeOutput(std::vector<float> &m_image, std::string output, const Vector2i &resolution)
	{
		
		const char* channelNames[] = {
			"0_hx", "1_hy", "2_dx", "3_dy", "4_R", "5_G", "6_B", 
			"7_hx", "8_hy", "9_dx", "a_dy", "b_R", "c_G", "d_B", 
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

	// cross product of woTheta vectors
	void cross_product (Float* v1, Float* v2, Float* out)
	{
		out[0] = v1[1]*v2[2] - v1[2]*v2[1];
		out[1] = v1[2]*v2[0] - v1[0]*v2[2];
		out[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}
	
	// normalize vector
	void normalize(Float* v)
	{
		// normalize
		Float len = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
		v[0] = v[0] / len;
		v[1] = v[1] / len;
		v[2] = v[2] / len;
	}

	// rotate vector along one axis
	void rotate_vector(Float* vector, Float* axis, Float angle, Float* out)
	{
		Float temp;
		Float cross[3];
		Float cos_ang = cos(angle);
		Float sin_ang = sin(angle);

		out[0] = vector[0] * cos_ang;
		out[1] = vector[1] * cos_ang;
		out[2] = vector[2] * cos_ang;

		temp = axis[0]*vector[0]+axis[1]*vector[1]+axis[2]*vector[2];
		temp = temp*(1.0-cos_ang);

		out[0] += axis[0] * temp;
		out[1] += axis[1] * temp;
		out[2] += axis[2] * temp;

		cross_product (axis,vector,cross);
		
		out[0] += cross[0] * sin_ang;
		out[1] += cross[1] * sin_ang;
		out[2] += cross[2] * sin_ang;
	}

	// convert standard coordinates to half vector/difference vector coordinates
	void std_coords_to_half_diff_coords(Float theta_in, Float fi_in, Float theta_out, Float fi_out,
									Float& hTheta,Float& hPhi,Float& dTheta,Float& dPhi )
	{

		// compute in vector
		Float in_vec_z = cos(theta_in);
		Float proj_in_vec = sin(theta_in);
		Float in_vec_x = proj_in_vec*cos(fi_in);
		Float in_vec_y = proj_in_vec*sin(fi_in);
		Float in[3]= {in_vec_x,in_vec_y,in_vec_z};
		normalize(in);


		// compute out vector
		Float out_vec_z = cos(theta_out);
		Float proj_out_vec = sin(theta_out);
		Float out_vec_x = proj_out_vec*cos(fi_out);
		Float out_vec_y = proj_out_vec*sin(fi_out);
		Float out[3]= {out_vec_x,out_vec_y,out_vec_z};
		normalize(out);


		// compute halfway vector
		Float half_x = (in_vec_x + out_vec_x)/2.0f;
		Float half_y = (in_vec_y + out_vec_y)/2.0f;
		Float half_z = (in_vec_z + out_vec_z)/2.0f;
		Float half[3] = {half_x,half_y,half_z};
		normalize(half);

		// compute  hTheta, hPhi
		hTheta = acos(half[2]);
		hPhi = atan2(half[1], half[0]);


		Float bi_normal[3] = {0.0, 1.0, 0.0};
		Float normal[3] = { 0.0, 0.0, 1.0 };
		Float temp[3];
		Float diff[3];

		// compute diff vector
		rotate_vector(in, normal , -hPhi, temp);
		rotate_vector(temp, bi_normal, -hTheta, diff);
		
		// compute  dTheta, dPhi	
		dTheta = acos(diff[2]);
		dPhi = atan2(diff[1], diff[0]);

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

	double renderSingleSPP(Scene *scene,  const Point2 &offset, int idx_spp, std::string path, const bool LiLoops)
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

		std::string imageName = path;
		imageName.append("_");
		std::ostringstream ss;
		ss << idx_spp;
		std::string s(ss.str());

		imageName.append(s);
		imageName.append(".exr");

		std::string outPath = path;
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

				// [fjh] convert wiwo xyz to hd xyz
				Float wiTheta_0 = acos(wi_0.z);
				Float woTheta_0 = acos(wo_0.z);
				Float wiPhi_0 = atan2(wi_0.y, wi_0.x);
				Float woPhi_0 = atan2(wo_0.y, wo_0.x);
				double hTheta_0, hPhi_0, dTheta_0, dPhi_0;

				std_coords_to_half_diff_coords(wiTheta_0, wiPhi_0, woTheta_0, woPhi_0, hTheta_0, hPhi_0, dTheta_0, dPhi_0);

				Float wiTheta_1 = acos(wi_1.z);
				Float woTheta_1 = acos(wo_1.z);
				Float wiPhi_1 = atan2(wi_1.y, wi_1.x);
				Float woPhi_1 = atan2(wo_1.y, wo_1.x);
				double hTheta_1, hPhi_1, dTheta_1, dPhi_1;

				std_coords_to_half_diff_coords(wiTheta_1, wiPhi_1, woTheta_1, woPhi_1, hTheta_1, hPhi_1, dTheta_1, dPhi_1);

				auto it = image.begin() + (j * resImage.x + i) * NUM_CHANNEL;

				*(it++) = sin(hTheta_0) * cos(hPhi_0);
				*(it++) = sin(hTheta_0) * sin(hPhi_0);
				*(it++) = sin(dTheta_0) * cos(dPhi_0);
				*(it++) = sin(dTheta_0) * sin(dPhi_0);
				*(it++) = lighting_0[0];
				*(it++) = lighting_0[1];
				*(it++) = lighting_0[2];
				*(it++) = sin(hTheta_1) * cos(hPhi_1);
				*(it++) = sin(hTheta_1) * sin(hPhi_1);
				*(it++) = sin(dTheta_1) * cos(dPhi_1);
				*(it++) = sin(dTheta_1) * sin(dPhi_1);
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

		writeOutput(image, imageName, resImage);
		std::vector<float>().swap(image);
		return time;
	}

	void pathtracing(Scene *scene, const std::string path, const bool LiLoops)
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
				time += renderSingleSPP(scene, offset, idx, path, LiLoops);
				fprintf(stderr, "[PT] CPU processing: %d/%d spp, time: %f\r", idx, sqrtSPP * sqrtSPP, (float)time);
			}
		}
		cout << endl;
	}

	int run(int argc, char **argv) {
		ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver();
		optind = 1;

		// int optchar;
		// bool timeThis = false;

		/* Parse command-line arguments */
		// while ((optchar = getopt(argc, argv, "t")) != -1) {
		// 	switch (optchar) {
		// 		case 't':
		// 			timeThis = true;
		// 			break;
		// 		// case 'n':
		// 		// 	bounceCount = strtol(optarg, &end_ptr,10);
		// 		// 	if (*end_ptr != '\0')
		// 		// 		SLog(EError, "Could not parse the bounce count!");
		// 		// 	break;
		// 	};
		// }

		ref<Scene> scene;

		std::string lowercase = boost::to_lower_copy(std::string(argv[optind]));
		if (boost::ends_with(lowercase, ".xml")) {
			fs::path
				filename = fileResolver->resolve(argv[optind]),
				filePath = fs::absolute(filename).parent_path(),
				baseName = filename.stem();
			ref<FileResolver> frClone = fileResolver->clone();
			frClone->prependPath(filePath);
			Thread::getThread()->setFileResolver(frClone);
			scene = loadScene(argv[optind]);
		} else {
			Log(EError, "The supplied scene filename must end in XML!");
		}

		/* Show some statistics, and make sure it roughly fits in 80cols */
		Logger *logger = Thread::getThread()->getLogger();
		DefaultFormatter *formatter = ((DefaultFormatter *) logger->getFormatter());
		logger->setLogLevel(EDebug);
		formatter->setHaveDate(false);

		scene->initialize();		
		optind++;
		std::string outputPath = std::string(argv[optind]);

		pathtracing(scene, outputPath, true);

		return 0;
	}

	MTS_DECLARE_UTILITY()

private:
	//may return 0 when not able to detect
	const uint cpu_count = std::thread::hardware_concurrency() == 0 ? 1: std::thread::hardware_concurrency();
	
};

MTS_EXPORT_UTILITY(PathTracer, "Individual path tracer.")
MTS_NAMESPACE_END
