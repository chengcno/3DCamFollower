///////////////////////////////////////////////////////////////
//
// LMopt.h
//
//   LM optimizer
//
// by Yingjie Cheng
//
// 01/Dec/2020
//
//
///////////////////////////////////////////////////////////////

#ifndef LMEIGEN_H
#define LMEIGEN_H

#include <Eigen/Eigen>
#include <vector>
#include "IOfile/BSplineFit.h"
#include "IOfile/SplineFunc.h"

using namespace std;
using namespace Eigen;

struct LMFunctor
{
    BSplineFit *bsf;
    SplineFunc *sf;
    Vector3d folCen;
    Vector3d folJoint;
    vector<pair<double, double>> Key;
    int keys;

    // f = f(x) = ( ,fi(x), )
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const;

    // Jco = dfi/dxj
    int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const;

    // Number of data points, i.e. values.
    int m;

    // Returns 'm', the number of values.
    [[nodiscard]] int values() const { return m; }

    // The number of parameters, i.e. inputs.
    int n;

    // Returns 'n', the number of inputs.
    [[nodiscard]] int inputs() const { return n; }

    vector<double> GetKepa(const Eigen::VectorXd &x) const;
};


#endif //LMEIGEN_H
