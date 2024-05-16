/**
    [fjh] This code is adapted from Alex's NeuMIP.
    Added triangle interpolation for BTF evaluation.
**/

//
// Created by alex on 3/10/21.
//

#ifndef REALBTF_UBO2003_BTF_H
#define REALBTF_UBO2003_BTF_H

#define BTF_UBO2003_IMPLEMENTATION
#include "btf_ubo2003.hh"

// [fjh] Delaunay triangulation
#include <stdexcept>
#include <DelaunayTriangulation/DataStructure.h>
#include <DelaunayTriangulation/Triangulation.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/formatter.h>


inline float square(const float x) {
    return x*x;
}

inline Vector3 cross(const Vector3 &v1, const Vector3 &v2) {
    return Vector3{
        (v1.y * v2.z) - (v1.z * v2.y),
        (v1.z * v2.x) - (v1.x * v2.z),
        (v1.x * v2.y) - (v1.y * v2.x)
    };
}

inline float norm(const Vector3 &v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

inline double determinant(const double matrix[]) {
    // inversed for left handed coordinate system
    double determinant = matrix[2] * matrix[4] * matrix[6]
        + matrix[0] * matrix[5] * matrix[7]
        + matrix[1] * matrix[3] * matrix[8]
        - matrix[0] * matrix[4] * matrix[8]
        - matrix[1] * matrix[5] * matrix[6]
        - matrix[2] * matrix[3] * matrix[7];

    // adjust result based on float number accuracy, otherwise causing deadloop
    return abs(determinant) <= std::numeric_limits<double>::epsilon() ? 0 : determinant;
}

inline double determinant(const Vector3& v0, const Vector3& v1, const Vector3& v2) {
    const double matrix[] = {
        v0.x, v0.y, v0.z,
        v1.x, v1.y, v1.z,
        v2.x, v2.y, v2.z
    };

    return determinant(matrix);
}

class UBO2003_btf {
public:
    BTF *btf = nullptr;
    int available_texels_u, available_texels_v; // [fjh] this parameter works like $trian_uv_size in python codes.
    int u_interval, v_interval;
    int offset_texels_u, offset_texels_v;
    char *interpType = nullptr;
    std::vector<std::tuple<int, int, int> *> delaunay_triangles;

    UBO2003_btf() {}

    UBO2003_btf(
        const char *path,
        int available_texels_u, int available_texels_v,
        int u_interval, int v_interval,
        int offset_texels_u, int offset_texels_v,
        const char *interpType) : available_texels_u(available_texels_u), available_texels_v(available_texels_v),
                                  u_interval(u_interval), v_interval(v_interval),
                                  offset_texels_u(offset_texels_u), offset_texels_v(offset_texels_v)
    {
        btf = LoadBTF(path);
        this->interpType = new char[strlen(interpType) + 1];
        strcpy(this->interpType, interpType); // [fjh] here we must use a strcpy instead of a simple assignment for string members
        assert(
            static_cast<uint32_t>(available_texels_u * u_interval) <= btf->Width && available_texels_u > 0 && u_interval > 0 
            && static_cast<uint32_t>(available_texels_v * v_interval) <= btf->Height && available_texels_v > 0 && v_interval > 0);
    }

    void print_dirs() const{
        for (uint32_t i = 0; i < btf->ViewCount; i++) {
            auto view_value =btf->Lights[i];
            // float theta = atan2(view_value.y, view_value.x);
            // float phi = asin(view_value.z);
            //std::cout << theta*180/M_PI << "\t" << phi*180/M_PI << std::endl;
            std::cout << "["  << view_value.x << ", " << view_value.y << ", " << view_value.z << "]," << std::endl;
            view_value =btf->Views[i];
            std::cout << "["  << view_value.x << ", " << view_value.y << ", " << view_value.z << "]," << std::endl;
        }
    }


    Spectrum get_value(
        float u, float v, 
        float light_x, float light_y, float light_z, 
        float view_x, float view_y, float view_z) const {

        // assert(u >= 0. && u <= 1. && v >= 0. && v <= 1. && "[get_value] uv should between (0, 1)");
        //! use repeat uv
        u = u - floor(u);
        v = v - floor(v);

        int pos_x = floor(u * available_texels_u); // index of texels in this available patch of the BTF
        int pos_y = floor(v * available_texels_v); 

        // make sure that pos < num_texels, so where uv = 1.0 will be num_texels - 1
        if (pos_x == available_texels_u)
            pos_x = available_texels_u -1;
        if (pos_y == available_texels_v)
            pos_y = available_texels_v -1;

        if (!strcmp(interpType, "floor"))
        {
            return get_value_by_integer_xy(pos_x + offset_texels_u, pos_y + offset_texels_v, light_x, light_y, light_z, view_x, view_y, view_z);
        }
        if (!strcmp(interpType, "bilinear"))
        {
            Spectrum res00, res01, res10, res11;
            float r0 = u * available_texels_u - pos_x, r1 = v * available_texels_v - pos_y;

            res00 = get_value_by_integer_xy((pos_x + 0) % available_texels_u + offset_texels_u, (pos_y + 0) % available_texels_v + offset_texels_v, light_x, light_y, light_z, view_x, view_y, view_z);
            res01 = get_value_by_integer_xy((pos_x + 0) % available_texels_u + offset_texels_u, (pos_y + 1) % available_texels_v + offset_texels_v, light_x, light_y, light_z, view_x, view_y, view_z);
            res10 = get_value_by_integer_xy((pos_x + 1) % available_texels_u + offset_texels_u, (pos_y + 0) % available_texels_v + offset_texels_v, light_x, light_y, light_z, view_x, view_y, view_z);
            res11 = get_value_by_integer_xy((pos_x + 1) % available_texels_u + offset_texels_u, (pos_y + 1) % available_texels_v + offset_texels_v, light_x, light_y, light_z, view_x, view_y, view_z);
            // std::cout << "get value at [ubo2003_btf.h line 136]" << std::endl;
            return res00 * (1 - r0) * (1 - r1) + res01 * (1 - r0) * r1 + res10 * r0 * (1 - r1) + res11 * r0 * r1;
        }
        else
        {
            SLog(mitsuba::EError, "interpType is something else!");
            return Spectrum{0.0f, 0.0f, 0.0f};
        }
    }

    /* 
        [fjh] this function takes integer uv coordinates, and returns the corresponding BTF value by potential triangulation interpolation.
        NOTE: before this function, all uv computations are actually on the "shrinked" "skipped-sampled" BTF (i.e, available_texels_u/v). 
        (size: [uv_size, uv_size] = [num_texels // u_interval, num_texels // v_interval])
        in this function, the previously determined integer uv will be multiplied by the $u/v_interval to get the actual uv coordinates on the original BTF. (of course, still integers)
    */
    Spectrum get_value_by_integer_xy(
        uint32_t pos_x, uint32_t pos_y,
        float light_x, float light_y, float light_z,
        float view_x, float view_y, float view_z) const {

        /* delaunay triangle interpolation */

        if (!delaunay_triangles.empty()) // if $useTriangulation is true (by default) in xml
        {
            int light_indices[3] = {0, 0, 0};
            int view_indices[3] = {0, 0, 0};
            float light_weights[3] = {0., 0., 0.};
            float view_weights[3] = {0., 0., 0.};

            if (!find_triangle(light_indices, light_weights, light_x, light_y, light_z))
            {
                // throw std::runtime_error("[ubo2003_btf.h: get_value] cannot find triangle.\n");
                SLog(mitsuba::EDebug, "[ubo2003_btf.h: get_value] cannot find triangle for light: %f, %f, %f\n", light_x, light_y, light_z);
                goto NEAREST;
            }

            if (!find_triangle(view_indices, view_weights, view_x, view_y, view_z))
            {
                // throw std::runtime_error("[ubo2003_btf.h: get_value] cannot find triangle.\n");
                SLog(mitsuba::EDebug, "[ubo2003_btf.h: get_value] cannot find triangle for view: %f, %f, %f\n", view_x, view_y, view_z);
                goto NEAREST;
            }

            uint32_t u_inter = static_cast<uint32_t>(u_interval);
            uint32_t v_inter = static_cast<uint32_t>(v_interval);

            Spectrum res00 = BTFFetchSpectrum(btf, light_indices[0], view_indices[0], pos_x * u_inter, pos_y * v_inter);
            Spectrum res01 = BTFFetchSpectrum(btf, light_indices[0], view_indices[1], pos_x * u_inter, pos_y * v_inter);
            Spectrum res02 = BTFFetchSpectrum(btf, light_indices[0], view_indices[2], pos_x * u_inter, pos_y * v_inter);
            Spectrum res10 = BTFFetchSpectrum(btf, light_indices[1], view_indices[0], pos_x * u_inter, pos_y * v_inter);
            Spectrum res11 = BTFFetchSpectrum(btf, light_indices[1], view_indices[1], pos_x * u_inter, pos_y * v_inter);
            Spectrum res12 = BTFFetchSpectrum(btf, light_indices[1], view_indices[2], pos_x * u_inter, pos_y * v_inter);
            Spectrum res20 = BTFFetchSpectrum(btf, light_indices[2], view_indices[0], pos_x * u_inter, pos_y * v_inter);
            Spectrum res21 = BTFFetchSpectrum(btf, light_indices[2], view_indices[1], pos_x * u_inter, pos_y * v_inter);
            Spectrum res22 = BTFFetchSpectrum(btf, light_indices[2], view_indices[2], pos_x * u_inter, pos_y * v_inter);

            Spectrum res0 = res00 * view_weights[0] + res01 * view_weights[1] + res02 * view_weights[2]; 
            Spectrum res1 = res10 * view_weights[0] + res11 * view_weights[1] + res12 * view_weights[2]; 
            Spectrum res2 = res20 * view_weights[0] + res21 * view_weights[1] + res22 * view_weights[2];

            return res0 * light_weights[0] + res1 * light_weights[1] + res2 * light_weights[2];
        }

        /* nearest */

        else
        {
        NEAREST:
            int light_ind = find_closest_index(light_x, light_y);
            int view_ind = find_closest_index(view_x, view_y);

            return BTFFetchSpectrum(btf, light_ind, view_ind, pos_x * u_interval, pos_y * v_interval);
        }
    }

    /* [fjh] locate the delaunay triangle that contain the given direction (light_x, light_y, light_z) */
    bool find_triangle(int *min_indices, float *weights, const float light_x, const float light_y, const float light_z) const
    {
        for (auto it = delaunay_triangles.begin(); it != delaunay_triangles.end(); it++)
        {
            Vector3 &ver0 = btf->Lights[std::get<0>(**it)];
            Vector3 &ver1 = btf->Lights[std::get<1>(**it)];
            Vector3 &ver2 = btf->Lights[std::get<2>(**it)];
            if (check_triangle(light_x, light_y, light_z, ver0, ver1, ver2))
            {
                min_indices[0] = std::get<0>(**it);
                min_indices[1] = std::get<1>(**it);
                min_indices[2] = std::get<2>(**it);

                intersection_weights(weights, light_x, light_y, light_z, ver0, ver1, ver2);

                return true;
            }
        }
        return false;
    }

    /* [fjh] check if the given vector (light_x, light_y, light_z) pass through the triangle (ver0, ver1, ver2) */
    bool check_triangle(const float light_x, const float light_y, const float light_z, const Vector3& ver0,  const Vector3& ver1, const Vector3& ver2) const {
        if (determinant(ver0, ver1, ver2) < 0) // if vertices are counter-clockwise
            return check_triangle(light_x, light_y, light_z, ver0, ver2, ver1);
        double det[] = { 0, 0, 0 };
        Vector3 point = {light_x, light_y, light_z};
        // [fjh] check if $point is at the same side against all edges
        det[0] = determinant(ver0, ver1, point);
        det[1] = determinant(ver1, ver2, point);
        det[2] = determinant(ver2, ver0, point);
        return (det[0] >= 0 && det[1] >= 0 && det[2] >= 0);
    }

    /* [fjh] Given a vector (light_x, light_y, light_z) that pass through a triangle (ver0, ver1, ver2), find the intersection and the median point coordinate weights */
    void intersection_weights(float* weights, const float light_x, const float light_y, const float light_z, const Vector3& ver0,  const Vector3& ver1, const Vector3& ver2) const {

        Vector3 normal = cross(ver1 - ver0, ver2 - ver0);
        Vector3 point = {light_x, light_y, light_z};
        float b = Dot(normal, point);
        float a = Dot(normal, ver0);
        float r = a / b;

        Vector3 intersectionPoint = r * point;

        float area0 = norm(cross(ver1 - intersectionPoint, ver2 - ver1));
        float area1 = norm(cross(ver2 - intersectionPoint, ver0 - ver2));
        float area2 = norm(cross(ver0 - intersectionPoint, ver1 - ver0));

        weights[0] = area0 / (area0 + area1 + area2);
        weights[1] = area1 / (area0 + area1 + area2);
        weights[2] = area2 / (area0 + area1 + area2);

    }

    int find_closest_index(float light_x, float light_y) const {
        int min_index = 0;
        float min_distance = std::numeric_limits<float>::infinity();

        for (uint32_t cur_index = 0; cur_index < btf->ViewCount; cur_index++) {

            auto &view_value =btf->Lights[cur_index];
            float cur_distance =square(view_value.x - light_x) + square(view_value.y - light_y);

            if (cur_distance < min_distance) {
                min_distance = cur_distance;
                min_index = cur_index;
            }
        }

        //std::cout << min_distance << std::endl;
        return min_index;
    }
    ~UBO2003_btf(){
        DestroyBTF(btf);
        btf = nullptr;
    };


};


#endif //REALBTF_UBO2003_BTF_H
