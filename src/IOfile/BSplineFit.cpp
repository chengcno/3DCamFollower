//
// Created by zz on 2020/12/25.
//

#include "BSplineFit.h"
#include <cmath>
#include <iostream>

BSplineFit::BSplineFit()
{
    scale_s = 0.35;
    // pigeon is 0.4
}

BSplineFit::BSplineFit(vector<Vector2d> _CtrPoints)
{
    CtrPoints = _CtrPoints;
}
BSplineFit::~BSplineFit() {}


/////
/// cubic spline
//////////////////
double BSplineFit::F03(double t)
{
    return pow(1-t,3)/6.0;
}

double BSplineFit::F13(double t)
{
    return (3*t*t*t - 6*t*t + 4)/6.0;
}

double BSplineFit::F23(double t)
{
    return (-3*t*t*t + 3*t*t + 3*t + 1)/6.0;
}

double BSplineFit::F33(double t)
{
    return pow(t,3)/6.0;
}


double BSplineFit::DF03(double t)
{
    return pow(1-t, 2)/2.0 * (-1);
}
double BSplineFit::DF13(double t)
{
    return (3*t*t - 4*t)/2.0;
}
double BSplineFit::DF23(double t)
{
    return (-3*t*t + 2*t + 1)/2.0;
}
double BSplineFit::DF33(double t)
{
    return t*t/2.0;
}

double BSplineFit::DDF03(double t)
{
    return 1-t;
}
double BSplineFit::DDF13(double t)
{
    return 3*t - 2;
}
double BSplineFit::DDF23(double t)
{
    return -3*t + 1;
}
double BSplineFit::DDF33(double t)
{
    return t;
}

double BSplineFit::DDDF03(double t)
{
    return -1;
}
double BSplineFit::DDDF13(double t)
{
    return 3;
}
double BSplineFit::DDDF23(double t)
{
    return -3;
}
double BSplineFit::DDDF33(double t)
{
    return 1;
}

/////
//// Output
/////////
vector<Vector2d> BSplineFit::SplineCurve(int nSize)
{
    vector<Vector2d> splineCurve;
    int Cnum = CtrPoints.size();
    for(int j=0;j<Cnum;j++)
    {
        for (int i = 0; i < nSize; i++)
        {
            double t = i * 1.0 / nSize;
            Vector2d curvePoint;
            curvePoint = CtrPoints[j] * F03(t) + CtrPoints[(j+1)%Cnum] * F13(t) + CtrPoints[(j+2)%Cnum] * F23(t) + CtrPoints[(j+3)%Cnum] * F33(t);
            splineCurve.emplace_back(curvePoint.x(), curvePoint.y());
        }
    }
    return splineCurve;
}

Vector2d BSplineFit::SplinePoint(double t)
{
    Vector2d pnt;
    double theta = 1.0/CtrPoints.size();
    int Num = CtrPoints.size();
    int id = floor(t/theta);
    double reT = (t - id*theta)/theta;
    pnt = CtrPoints[id]*F03(reT) + CtrPoints[(id+1)%Num]*F13(reT) + CtrPoints[(id+2)%Num]*F23(reT) + CtrPoints[(id+3)%Num]*F33(reT);
    return pnt;
}

/////
/// fit (interplote)
////////////
void BSplineFit::ValueToControl()
{
    int nValue = ValuePoints.size();
    SparseMatrix<double> triM(nValue, nValue);
    triM.reserve(3*nValue);

    for(int i=0;i<nValue;i++)
    {
        int bef,nex;
        bef = i-1 > -1 ? i-1 : i-1+nValue;
        nex = i+1 < nValue ? i+1 : i+1 - nValue;
        triM.coeffRef(i,bef) = 1;
        triM.coeffRef(i, i) = 4;
        triM.coeffRef(i, nex) = 1;
    }
    triM.makeCompressed();

    SimplicialCholesky<SparseMatrix<double>> solverTriM(triM);

    VectorXd Vx,Vy,Vz;
    Vx.resize(nValue);
    Vy.resize(nValue);
    Vz.resize(nValue);

    Vx(0) = (ValuePoints[nValue-1].x()*6);
    Vy(0) = (ValuePoints[nValue-1].y()*6);
    for(int i=0;i<nValue-1;i++)
    {
        Vx(i+1) = (ValuePoints[i].x()*6);
        Vy(i+1) = (ValuePoints[i].y()*6);
    }

    ///
    VectorXd Cx,Cy,Cz;
    Cx = solverTriM.solve(Vx);
    Cy = solverTriM.solve(Vy);

    CtrPoints.clear();
    for(int i=0;i<nValue;i++)
        CtrPoints.emplace_back(Cx(i), Cy(i));
}


////////////
//// gradient
//////////

Vector2d BSplineFit::GetSpCurve(double t)
{
    int nSize = CtrPoints.size();
    if(nSize == 0)
        cout<<"Control points are empty!"<<endl;

    t = tin0_2pi(t);

    double theta = 2.0*M_PI/nSize;
    int id = floor(t/theta);
    double rel = (t - id*theta)/theta;
    if(id == nSize) {id = 0;rel = 0;}
    if(id > nSize) cout<< "Get Spline Curve para is wrong!"<<endl;

    Vector2d curvePoint;
    curvePoint = CtrPoints[id] * F03(rel) + CtrPoints[(id+1)%nSize] * F13(rel)
            + CtrPoints[(id+2)%nSize] * F23(rel) + CtrPoints[(id+3)%nSize] * F33(rel);

    return curvePoint;
}

Vector2d BSplineFit::GetDSpCurve(double t)
{
    int nSize = CtrPoints.size();
    if(nSize == 0)
        cout<<"Control points are empty!"<<endl;

    t = tin0_2pi(t);

    double theta = 2.0*M_PI/nSize;
    int id = floor(t/theta);
    double rel = (t - id*theta)/theta;
    if(id == nSize) {id = 0;rel = 0;}
    if(id > nSize) cout<< "Get Spline Curve para is wrong!"<<endl;

    Vector2d DcurvePoint;
    DcurvePoint = CtrPoints[id] * DF03(rel) + CtrPoints[(id+1)%nSize] * DF13(rel)
                 + CtrPoints[(id+2)%nSize] * DF23(rel) + CtrPoints[(id+3)%nSize] * DF33(rel);

    DcurvePoint = DcurvePoint/theta;
    return DcurvePoint;
}

Vector2d BSplineFit::GetDDSpCurve(double t)
{
    int nSize = CtrPoints.size();
    if(nSize == 0)
        cout<<"Control points are empty!"<<endl;

    t= tin0_2pi(t);

    double theta = 2.0*M_PI/nSize;
    int id = floor(t/theta);
    double rel = (t - id*theta)/theta;
    if(id == nSize) {id = 0;rel = 0;}
    if(id > nSize) cout<< "Get Spline Curve para is wrong!"<<endl;

    Vector2d DDcurvePoint;
    DDcurvePoint = CtrPoints[id] * DDF03(rel) + CtrPoints[(id+1)%nSize] * DDF13(rel)
                  + CtrPoints[(id+2)%nSize] * DDF23(rel) + CtrPoints[(id+3)%nSize] * DDF33(rel);

    DDcurvePoint = DDcurvePoint/theta;
    return DDcurvePoint;
}

Vector2d BSplineFit::GetDDDSpCurve(double t)
{
    int nSize = CtrPoints.size();
    if(nSize == 0)
        cout<<"Control points are empty!"<<endl;

    t = tin0_2pi(t);

    double theta = 2.0*M_PI/nSize;
    int id = floor(t/theta);
    double rel = (t - id*theta)/theta;
    if(id == nSize) {id = 0;rel = 0;}
    if(id > nSize) cout<< "Get Spline Curve para is wrong!"<<endl;

    Vector2d DDDcurvePoint;
    DDDcurvePoint = CtrPoints[id] * DDDF03(rel) + CtrPoints[(id+1)%nSize] * DDDF13(rel)
                  + CtrPoints[(id+2)%nSize] * DDDF23(rel) + CtrPoints[(id+3)%nSize] * DDDF33(rel);

    DDDcurvePoint = DDDcurvePoint/theta;
    return DDDcurvePoint;
}


Vector3d BSplineFit::GetTcurve(double t)
{
    Vector2d spPoint = GetSpCurve(t);
    Vector3d Tpoint = Vector3d (0, scale_s* (spPoint.y()), scale_s*spPoint.x());
    return Tpoint;
}

Vector3d BSplineFit::GetDTcurve(double t)
{
    Vector2d Dsp = GetDSpCurve(t);
    Vector3d DTp = Vector3d (0, scale_s*Dsp.y(), scale_s*Dsp.x());
    return DTp;
}

Vector3d BSplineFit::GetDDTcurve(double t)
{
    Vector2d Dsp = GetDDSpCurve(t);
    Vector3d DDTp = Vector3d (0, scale_s*Dsp.y(), scale_s*Dsp.x());
    return DDTp;
}

Vector3d BSplineFit::GetDDDTcurve(double t)
{
    Vector2d Dsp = GetDDDSpCurve(t);
    Vector3d DDDTp = Vector3d (0, scale_s*Dsp.y(), scale_s*Dsp.x());
    return DDDTp;
}

Vector3d BSplineFit::GetSphCurve(double t)
{
    Vector3d sphPoint = GetTcurve(t);
    Vector3d virsphcen = Vector3d (sphCen.x()/abs(sphCen.x()),0,0);
    return (sphPoint - virsphcen).normalized() *sphR + sphCen;
    //return (sphPoint-sphCen).normalized()*sphR + sphCen;
}

Vector3d BSplineFit::GetDSphCurve(double t)
{
    Vector3d sphPoint = GetTcurve(t);
    Vector3d DsphPoint = GetDTcurve(t);
    Vector3d virsphcen = Vector3d (sphCen.x()/abs(sphCen.x()),0,0);

    return (DsphPoint)*sphR/(sphPoint - virsphcen).norm()
    - sphR*(sphPoint - virsphcen)/pow((sphPoint-virsphcen).norm(),3)*(DsphPoint.dot(sphPoint-virsphcen)); ;
}
Vector3d BSplineFit::GetDDSphCurve(double t)
{
    Vector3d sphPoint = GetTcurve(t);
    Vector3d DsphPoint = GetDTcurve(t);
    Vector3d DDsphPoint = GetDDTcurve(t);
    Vector3d virsphcen = Vector3d (sphCen.x()/abs(sphCen.x()),0,0);

    double DPO_1 = 1.0/pow((sphPoint-virsphcen).norm(),3)*(DsphPoint.dot(sphPoint-virsphcen));
    return DDsphPoint*sphR/(sphPoint - virsphcen).norm()
            + 2*DsphPoint*sphR*DPO_1
            + (sphPoint-virsphcen)*sphR*(3.0/pow((sphPoint-virsphcen).norm(),5) * pow(DsphPoint.dot(sphPoint-virsphcen),2)
                                        + (-1)/pow((sphPoint-virsphcen).norm(),3) * (DDsphPoint.dot(sphPoint-virsphcen) + DsphPoint.dot(DsphPoint)));
}
Vector3d BSplineFit::GetDDDSphCurve(double t)
{
    Vector3d sphPoint = GetTcurve(t);
    Vector3d DsphPoint = GetDTcurve(t);
    Vector3d DDsphPoint = GetDDTcurve(t);
    Vector3d DDDsphPoint = GetDDDTcurve(t);
    Vector3d virsphcen = Vector3d (sphCen.x()/abs(sphCen.x()),0,0);

    double DPO_1 = 1.0/pow((sphPoint-virsphcen).norm(),3)*(DsphPoint.dot(sphPoint-virsphcen));
    double DDPO_1 = 3.0/pow((sphPoint-virsphcen).norm(),5) * pow(DsphPoint.dot(sphPoint-virsphcen),2)
    + (-1)/pow((sphPoint-virsphcen).norm(),3) * (DDsphPoint.dot(sphPoint-virsphcen) + DsphPoint.dot(DsphPoint));
    double DPO_3 = -15/pow((sphPoint-virsphcen).norm(),7)*pow(DsphPoint.dot(sphPoint-virsphcen),3)
            + 6/pow((sphPoint-virsphcen).norm(),5)*(DsphPoint.dot(sphPoint-virsphcen))*(DDsphPoint.dot(sphPoint-virsphcen) + DsphPoint.dot(DsphPoint))
            + 3/pow((sphPoint-virsphcen).norm(),5)*(DsphPoint.dot(sphPoint-virsphcen))*(DDsphPoint.dot(sphPoint-virsphcen) + DsphPoint.dot(DsphPoint))
            + (-1)/pow((sphPoint-virsphcen).norm(),3)*(DDDsphPoint.dot(sphPoint-virsphcen) + 3*DDsphPoint.dot(DsphPoint));

    return (DDDsphPoint)*sphR/(sphPoint-virsphcen).norm()
            + 3*DDsphPoint*sphR*DPO_1 + 3*DsphPoint*sphR*DDPO_1 + (sphPoint-virsphcen)*sphR*DPO_3;
}

double BSplineFit::tin0_2pi(double t)
{
    while(t<0)
        t += 2*M_PI;
    while(t>2*M_PI)
        t = t-2*M_PI;
    return t;
}
