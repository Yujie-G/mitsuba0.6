#pragma once

#include <string>
#include <vector>

// #include "powitacq_rgb.h"

void save_merl(const std::vector<double> &x, const std::string &path);
//void half_angle_to_std(double theta_h, double theta_d, double phi_d,
//                       powitacq_rgb::Vector3f &wi_v,
//                       powitacq_rgb::Vector3f &wo_v);

void half_angle_to_std(double theta_h, double theta_d, double phi_d,
    mitsuba::Vector3f& wi_v,
    mitsuba::Vector3f& wo_v);
