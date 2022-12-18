//
// Created by temp on 3/1/22.
//

#ifndef SPATIALLINKAGES_CUBICSPLINE_H
#define SPATIALLINKAGES_CUBICSPLINE_H

#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

class CubicSpline {
public:
    CubicSpline(){};
    ~CubicSpline(){};

public:
    vector<Vector3d>  ctrPoints;
    int               ctrN;
    vector<Vector3d>  ValuePoints;
    int               valN;

public:
    static double F03(double t);
    static double F13(double t);
    static double F23(double t);
    static double F33(double t);

public:
    vector<Vector3d> GetSplineCurve(int nSize);
    void ValueToControl();

};


#endif //SPATIALLINKAGES_CUBICSPLINE_H
