//
// Created by zz on 2020/12/27.
//

#ifndef CAMMECHCREATOR_SPLINEFUNC_H
#define CAMMECHCREATOR_SPLINEFUNC_H

#include <vector>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

class SplineFunc {
public:
    int nSize;
    double theta;
    vector<double> fx;
    vector<double> Ff;
    vector<double> a,b,c,d;
    vector<double> m;

    MatrixXd M_Y;

public:
    SplineFunc(vector<double> _fx);
    SplineFunc(const double* S, int size_s);
    ~SplineFunc();

    double GetFx(double t);
    double GetDFx(double t);
    double GetDDFx(double t);
    double GetDDDFx(double t);

    vector<double> partFx(double t);
    vector<double> partDFx(double t);
    vector<double> partDDFx(double t);
    vector<double> partDDDFx(double t);


public:
    void ComputeCubicSpline();
    bool SameSpline(const double *S);

private:
    void GetIdrel(int &id, double &rel, double t);
    inline int delta(int i,int j);
};


#endif //CAMMECHCREATOR_SPLINEFUNC_H
