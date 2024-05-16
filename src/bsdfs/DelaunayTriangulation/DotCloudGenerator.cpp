#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
#include <DelaunayTriangulation/DotCloud.h>

#define PI 3.14159265
#define RADIUS 100

using namespace std;
using namespace dt;

vector<Vector3D*> DotCloudGenerator::GetSphericalDots()
{
    int dotCount;
    cout << "Dot Count: ";
    cin >> dotCount;

    ofstream file("Resource\\random_out.txt");
    vector<Vector3D*> dots = vector<Vector3D*>();
    dots.reserve(dotCount);

    srand((unsigned)time(NULL));
    for (int i = 0; i < dotCount; i++)
    {
        Vector3D* dot = GetRandomDotEvenlyDistributed();
        file << "# " << dot->X << " " << dot->Y << " " << dot->Z << " "
            << dot->R << " " << dot->G << " " << dot->B << " " << endl;
        dots.push_back(dot);
    }

    file.close();

    return dots;
}

Vector3D* DotCloudGenerator::GetRandomDot()
{
    // use spherical coordinate
    double phi = (rand() % 360) * PI / 180;
    double theta = (rand() % 360) * PI / 180;

    double x = RADIUS * sin(theta) * cos(phi);
    double y = RADIUS * sin(theta) * sin(phi);
    double z = RADIUS * cos(theta);

    return new Vector3D(x, y, z);
}

Vector3D* DotCloudGenerator::GetRandomDotEvenlyDistributed()
{
    // project random dot in cartesian coordinate to unit sphere
    double x = rand() % 2000 - 1000;
    double y = rand() % 2000 - 1000;
    double z = rand() % 2000 - 1000;

    double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    x = RADIUS * x / r;
    y = RADIUS * y / r;
    z = RADIUS * z / r;

    return new Vector3D(x, y, z);
}