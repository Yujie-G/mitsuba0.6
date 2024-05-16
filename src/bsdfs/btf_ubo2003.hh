/*   The MIT License
*
*   Free BTF Library for BTF Dataset like Ubo2003
*   Copyright (c) 2022 Hanxiao Sun
*
*   Permission is hereby granted, free of charge, to any person obtaining a copy
*   of this software and associated documentation files (the "Software"), to deal
*   in the Software without restriction, including without limitation the rights
*   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*   copies of the Software, and to permit persons to whom the Software is
*   furnished to do so, subject to the following conditions:
*
*   The above copyright notice and this permission notice shall be included in
*   all copies or substantial portions of the Software.
*
*   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*   THE SOFTWARE.
*/
#ifndef _BTF_UBO2003_HH_
#define _BTF_UBO2003_HH_

#include <memory>
#include <string>
#include <vector>
#include <cmath>
// #include <experimental/filesystem>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include <iostream>
#include <fstream>
#include <sstream>
#include "simpleBMP.h"

#include <assert.h>
#include <cstring>
#include <immintrin.h>

#include <omp.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/formatter.h>

#define BTF_MAX_CHANNEL_COUNT 3
#define PI 3.141592653589793


#define LIGHT_COUNT 81
#define VIEW_COUNT 81
#define RESOLUTION 256


struct Vector2
{
    float x, y;
};

#ifndef _BTF_UBO2003_HH_VECTOR3_
#define _BTF_UBO2003_HH_VECTOR3_
struct Vector3
{
    float x, y, z;
};
#endif
struct Matrix3
{
    Vector3 tangent, binormal, normal;
};

typedef Vector3 Spectrum;
struct BTF;
inline Spectrum BTFFetchSpectrum(const BTF *btf, uint32_t light_vert, uint32_t view_vert, uint32_t x, uint32_t y);
inline Spectrum BTFFetchSpectrum(const BTF *btf, uint32_t lv_idx, uint32_t xy_idx);

BTF *LoadBTF(const char *file_path);
void DestroyBTF(BTF *btf);

struct Image
{
    uint8_t* R;
    uint8_t* G;
    uint8_t* B;
};

#define VIEW_INDEX_STEP 81*256*256
#define LIGHT_INDEX_STEP 256*256

struct BTF
{
    uint32_t ConsineFlag = false;
    uint32_t ChannelCount = 3;
    uint32_t Width = 0,
            Height = 0;
    uint32_t DynamicRangeReduction = false;
    struct
    {
        uint32_t Width = 0;
        uint32_t Height = 0;
    } HeightMapSize;

    uint16_t *HeightMap = nullptr;

    uint32_t ColorModel = 0;
    Vector3 ColorMean;
    Matrix3 ColorTransform;

    uint32_t RowCount = 0,
            ColumnCount = 0,
            DataSize = 0;

    uint8_t *LeftSingularU = nullptr,
            *RightSingularSxV = nullptr;

    uint64_t LeftSingularUSize = 0,
            RightSingularSxVSize = 0;

    Vector3 *Views = nullptr;
    Vector3 *Lights = nullptr;

    uint32_t ViewCount = 0;
    uint32_t LightCount = 0;

    uint32_t *LightIndices = nullptr;
    uint32_t LightTriangleCount = 0;

    uint32_t UElementStride = 0,
            SxVElementStride = 0;

    uint32_t *Offsets = nullptr;
    uint32_t *ComponentCounts = nullptr;
};

uint8_t BTF_data_R[81*81*256*256];
uint8_t BTF_data_G[81*81*256*256];
uint8_t BTF_data_B[81*81*256*256];

struct BTFDeleter
{
    inline void operator()(BTF *btf) { DestroyBTF(btf); }
};

typedef std::unique_ptr<BTF, BTFDeleter> BTFPtr;

inline float Dot(const Vector3 &lhs, const Vector3 &rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

inline Vector3 operator*(const Vector3 &lhs, const Vector3 &rhs) { return Vector3{lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z}; }

inline Vector3 operator*(const Vector3 &lhs, const float &rhs) { return Vector3{lhs.x * rhs, lhs.y * rhs, lhs.z * rhs}; }

inline Vector3 operator*(const float &lhs, const Vector3 &rhs) { return Vector3{lhs * rhs.x, lhs * rhs.y, lhs * rhs.z}; }

inline Vector3 operator+(const Vector3 &lhs, const Vector3 &rhs) { return Vector3{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z}; }

inline Vector3 operator-(const Vector3 &lhs, const Vector3 &rhs) { return Vector3{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z}; }

inline Vector3 operator-(const Vector3 &vec, float scalar) { return Vector3{vec.x - scalar, vec.y - scalar, vec.z - scalar}; }

inline float Dot(const Vector2 &lhs, const Vector2 &rhs) { return lhs.x * rhs.x + lhs.y * rhs.y; }

inline Vector3 ParabolicToCartesianCoordinates(const Vector2 &coordinates)
{
    float len_sq = Dot(coordinates, coordinates);
    return {2.0f * coordinates.x / (len_sq + 1.0f), 2.0f * coordinates.y / (len_sq + 1.0f), (1 - len_sq) / (1 + len_sq)};
}

inline Spectrum BTFFetchSpectrum(BTF *btf, uint32_t light_vert, uint32_t view_vert, uint32_t x, uint32_t y)

{
    assert(light_vert < btf->LightCount && view_vert < btf->LightCount &&
           x < btf->Width && y < btf->Height && "Invalid index");
    uint32_t xy_idx = y * btf->Width + x;

    Spectrum res;

    res.x = BTF_data_R[view_vert*VIEW_INDEX_STEP + light_vert*LIGHT_INDEX_STEP + y*RESOLUTION + x] / 255.0f;
    res.y = BTF_data_G[view_vert*VIEW_INDEX_STEP + light_vert*LIGHT_INDEX_STEP + y*RESOLUTION + x] / 255.0f;
    res.z = BTF_data_B[view_vert*VIEW_INDEX_STEP + light_vert*LIGHT_INDEX_STEP + y*RESOLUTION + x] / 255.0f;
 
    return res;
}

inline float Clampf(float val, float minval, float maxval)
{
    return std::fmax(std::fmin(val, maxval), minval);
}

inline Vector3 Clamp(const Vector3 &val, float smin, float smax)
{
    return Vector3{Clampf(val.x, smin, smax),
                   Clampf(val.y, smin, smax),
                   Clampf(val.z, smin, smax)};
}

inline Vector3 YUVToRGB(const Vector3 &color)
{
    Vector3 vx{1.0f, 0.0f, 1.13983f};
    Vector3 vy{1.0f, -0.39465f, -0.5806f};
    Vector3 vz{1.0f, 2.03211f, 0.0f};

    return Clamp(Vector3{Dot(vx, color), Dot(vy, color), Dot(vz, color)}, 0.0f, 1.0f);
}

inline Vector3 Exp(const Vector3 &val)
{
    return Vector3{expf(val.x), expf(val.y), expf(val.z)};
}

inline Spectrum BTFFetchSpectrum(const BTF *btf, uint32_t lv_idx, uint32_t xy_idx)
{
    uint32_t chan_count = btf->ChannelCount;

    if (chan_count > BTF_MAX_CHANNEL_COUNT)
        return {};

    Spectrum result{};
    float drr_eps = 1e-5f;

    switch (chan_count)
    {
        case 3:
        {
            Vector3 interm_color{};

            switch (btf->ColorModel)
            {
                case 0:
                    break;
                case 11:
                {
                    Vector3 yuv;
                    yuv.x = expf(interm_color.x) - drr_eps;
                    yuv.y = interm_color.y * yuv.x + drr_eps;
                    yuv.z = interm_color.z * yuv.x + drr_eps;
                    interm_color = YUVToRGB(yuv);
                }
                    break;
                default:
                {
                    assert(false && "Unsupported color model");
                }
            }

            if (btf->DynamicRangeReduction)
            {
                interm_color = Exp(interm_color) - drr_eps;
            }
            result = interm_color;
        }
            break;
        default:
        {
            assert(false && "Conversion is unsupported");
        }
            break;
    }

    if (btf->ConsineFlag)
    {
        assert(false && "Stub");
        // TODO
    }

    return result;
}

#endif // _BTF_UBO2003_HH_


#include <xmmintrin.h>
#include <immintrin.h>
#include <fstream>
#include <algorithm>

void DestroyBTF(BTF *btf)
{
    delete[] btf->HeightMap;
    btf->HeightMap = nullptr;
    delete[] btf->LeftSingularU;
    btf->LeftSingularU = nullptr;
    delete[] btf->RightSingularSxV;
    btf->RightSingularSxV = nullptr;
    delete[] btf->Views;
    btf->Views = nullptr;
    delete[] btf->Lights;
    btf->Lights = nullptr;
    delete[] btf->Offsets;
    btf->Offsets = nullptr;
    delete[] btf->ComponentCounts;
    btf->ComponentCounts = nullptr;
    free(btf->LightIndices);
    btf->LightIndices = nullptr;
  

}

#define BTF_MAKE_FOURCC(x, y, z, w) \
    (((x & 0xFFU)) |                \
     ((y & 0xFFU) << 8U) |          \
     ((z & 0xFFU) << 16U) |         \
     ((w & 0xFFU) << 24U))

#define BTF_MAKE_EIGHTCC(c0, c1, c2, c3, c4, c5, c6, c7, c8) \
    (((c0 & 0xFFU))       | \
    ((c1 & 0xFFU) << 8U)  | \
    ((c2 & 0xFFU) << 16U) | \
    ((c3 & 0xFFU) << 24U))| \
    ((c4 & 0xFFU) << 32U))| \
    ((c5 & 0xFFU) << 40U))| \
    ((c6 & 0xFFU) << 48U))| \
    ((c7 & 0xFFU) << 56U))

inline void SinCos_Tan(float omega, float *s, float *c)
{
    float t = tanf(omega * 0.5f);
    float t2 = t * t;
    float div = 1.0f / (1.0f + t2);
    *s = 2.0f * t * div;
    *c = (1.0f - t2) * div;
}

#define FastSinCos SinCos_Tan

inline Vector3 SphereToCartesianCoordinates(const Vector2 &angles)
{
    float cos_theta, sin_theta, cos_phi, sin_phi;
    FastSinCos(angles.x, &sin_theta, &cos_theta);
    FastSinCos(angles.y, &sin_phi, &cos_phi);

    return {cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};
}


BTF btf;
BTF *LoadBTF(const char *file_path)
{
    SLog(mitsuba::EInfo, "load start");
    // initialization
    std::string data_path = file_path;
    
    btf.LightCount = LIGHT_COUNT;
    btf.ViewCount = VIEW_COUNT;
    btf.Lights = new Vector3[LIGHT_COUNT]();
    btf.Views = new Vector3[VIEW_COUNT]();
    // btf->Imgs = new Image*[VIEW_COUNT]();
    btf.Width = RESOLUTION;
    btf.Height = RESOLUTION;

    Vector3 dir;
    char trash;
    int i = 0;

    // [c++11] directory_iterator is in std::tr2::sys namespace
	// [c++17] std::filesystem::directory_iterator
    std::vector<fs::directory_entry> vec_view;
    for (fs::directory_iterator it{data_path}; it != fs::directory_iterator{}; it++)
    {
        vec_view.push_back(*it);
    }

    // for (fs::directory_iterator it_view{data_path}; it_view != fs::directory_iterator{}; it_view++)
    for (auto& it_view : fs::directory_iterator(data_path))
    // #pragma omp parallel for
    // for (auto it_view = vec_view.begin(); it_view < vec_view.end(); it_view++)
    {
        // read tv(theta_view) and pv(phi_view) 
        std::string filename = it_view.path().filename().string();
        SLog(mitsuba::EDebug, "it: %s", filename.c_str());
        std::istringstream iss(filename.c_str());
        int tv;
        int pv;
        if (!filename.compare(0, 2, "tv")) {
            iss >> trash;
            iss >> trash;
            iss >> tv;
        }
        iss >> trash;
        if (!filename.compare(6, 2, "pv")) {
            iss >> trash;
            iss >> trash;
            iss >> pv;
        }
        
        dir.x = std::sin(PI * tv / 180.0) * std::cos(PI * pv / 180.0);
        dir.y = std::sin(PI * tv / 180.0) * std::sin(PI * pv / 180.0);
        dir.z = std::cos(PI * tv / 180.0);

        btf.Lights[i] = dir;
        btf.Views[i] = dir;

        // load image
        int j = 0;
        // for(auto& it_light : std::tr2::sys::directory_iterator(data_path + std::string("//") + filename))
        for(auto& it_light : fs::directory_iterator(data_path + std::string("//") + filename))
        {
            simpleBMP::BMPImg img;
            img.LoadImage(it_light.path().string());

            for(int x = 0; x < RESOLUTION; x++)
                for(int y = 0; y < RESOLUTION; y++)
                {
                    BTF_data_B[i*VIEW_INDEX_STEP + j*LIGHT_INDEX_STEP + x*RESOLUTION + y] = img.At(x,y)[0];
                    BTF_data_G[i*VIEW_INDEX_STEP + j*LIGHT_INDEX_STEP + x*RESOLUTION + y] = img.At(x,y)[1];
                    BTF_data_R[i*VIEW_INDEX_STEP + j*LIGHT_INDEX_STEP + x*RESOLUTION + y] = img.At(x,y)[2];
                }
            ++j;
        }
        ++i;
    }
    SLog(mitsuba::EInfo, "load finish");
    return &btf;
}


