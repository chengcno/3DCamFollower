//
// Created by temp on 25/11/21.
//

#ifndef SPATIALLINKAGES_LMOPT_H
#define SPATIALLINKAGES_LMOPT_H

#include <vector>
#include <Eigen/Eigen>
using namespace std;
using namespace Eigen;

class LMopt{

public:
    LMopt(int _N);
    ~LMopt();

public:
    int Nx;
    vector<double> res_x, init_x;

    vector<Vector2d> pinPos;
    vector<double> active_x;

public:
    void InitX(vector<double> _x);
    double Solver();

private:
    VectorXd varX, varF;
    MatrixXd Jcob;

private:
    VectorXd GetFx(VectorXd _x);
    MatrixXd GetJcob(VectorXd _x);
    MatrixXd GetHess(VectorXd _x, VectorXd _f);

};


#endif //SPATIALLINKAGES_LMOPT_H
