//
// Created by Yingjie CHNEG
// 20/Dec/2020
//

#ifndef POLYLINEPARSER_BSPLINEFIT_H
#define POLYLINEPARSER_BSPLINEFIT_H

#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

class BSplineFit {
public:
    vector<Vector2d> CtrPoints;
    vector<Vector2d> ValuePoints;

    double sphR;
    Vector3d sphCen;
    double scale_s;

public:
    BSplineFit();
    BSplineFit(vector<Vector2d> _CtrPoints);
    ~BSplineFit();

public:
    static double F03(double t);
    static double F13(double t);
    static double F23(double t);
    static double F33(double t);

    static double DF03(double t);
    static double DF13(double t);
    static double DF23(double t);
    static double DF33(double t);

    static double DDF03(double t);
    static double DDF13(double t);
    static double DDF23(double t);
    static double DDF33(double t);

    static double DDDF03(double t);
    static double DDDF13(double t);
    static double DDDF23(double t);
    static double DDDF33(double t);

public:
    vector<Vector2d> SplineCurve(int nSize);
    Vector2d SplinePoint(double t);
    void ValueToControl();


public:
    Vector2d GetSpCurve(double t);
    Vector2d GetDSpCurve(double t);
    Vector2d GetDDSpCurve(double t);
    Vector2d GetDDDSpCurve(double t);

    Vector3d GetTcurve(double t);
    Vector3d GetDTcurve(double t);
    Vector3d GetDDTcurve(double t);
    Vector3d GetDDDTcurve(double t);

    Vector3d GetSphCurve(double t);
    Vector3d GetDSphCurve(double t);
    Vector3d GetDDSphCurve(double t);
    Vector3d GetDDDSphCurve(double t);

    //todo: GetSphCurve(double t);

private:
    double tin0_2pi(double t);
};


#endif //POLYLINEPARSER_BSPLINEFIT_H
