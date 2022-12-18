//
// Created by temp on 3/1/22.
//

#include "CubicSpline.h"

double CubicSpline::F03(double t) {
    return pow(1-t,3)/6.0;
}

double CubicSpline::F13(double t) {
    return (3*t*t*t - 6*t*t + 4)/6.0;
}

double CubicSpline::F23(double t) {
    return (-3*t*t*t + 3*t*t + 3*t + 1)/6.0;
}

double CubicSpline::F33(double t) {
    return pow(t,3)/6.0;
}

void CubicSpline::ValueToControl()
{
    valN = ValuePoints.size();
    SparseMatrix<double> triM(valN, valN);
    triM.reserve(3 * valN);

    for(int i=0; i < valN; i++)
    {
        int bef,nex;
        bef = i-1 > -1 ? i-1 : i - 1 + valN;
        nex = i+1 < valN ? i + 1 : i + 1 - valN;
        triM.coeffRef(i,bef) = 1;
        triM.coeffRef(i, i) = 4;
        triM.coeffRef(i, nex) = 1;
    }

    triM.makeCompressed();
    SimplicialCholesky<SparseMatrix<double>> solverTriM(triM);

    ///solverTriM
    VectorXd Vx,Vy,Vz;
    Vx.resize(valN);
    Vy.resize(valN);
    Vz.resize(valN);

    Vx(0) = (ValuePoints[valN - 1].x() * 6);
    Vy(0) = (ValuePoints[valN - 1].y() * 6);
    Vz(0) = (ValuePoints[valN-1].z()*6);
    for(int i=0; i < valN - 1; i++)
    {
        Vx(i+1) = (ValuePoints[i].x()*6);
        Vy(i+1) = (ValuePoints[i].y()*6);
        Vz(i+1) = (ValuePoints[i].z()*6);
    }

    VectorXd Cx,Cy,Cz;
    Cx = solverTriM.solve(Vx);
    Cy = solverTriM.solve(Vy);
    Cz = solverTriM.solve(Vz);

    ctrPoints.clear();
    for(int i=0; i < valN; i++)
        ctrPoints.emplace_back(Cx(i), Cy(i), Cz(i));
    ctrN = valN;
}

vector<Vector3d> CubicSpline::GetSplineCurve(int nSize)
{
    vector<Vector3d> splineCurve;

    for(int j=0;j<ctrN;j++)
    {
        for (int i = 0; i < nSize; i++)
        {
            double t = i * 1.0 / nSize;
            Vector3d curvePoint;
            curvePoint = ctrPoints[j] * F03(t) + ctrPoints[(j+1)%valN] * F13(t) + ctrPoints[(j+2)%valN] * F23(t) + ctrPoints[(j+3)%valN] * F33(t);
            splineCurve.emplace_back(curvePoint.x(), curvePoint.y(), curvePoint.z());
        }
    }
    return splineCurve;
}


