//
// Created by temp on 25/11/21.
//

#include "LMopt.h"
#include <iostream>
#include "Utility/HelpDefine.h"

LMopt::LMopt(int _N) {
    Nx = _N;
}
LMopt::~LMopt() {}

void LMopt::InitX(vector<double> _x)
{
    init_x = _x;
    Nx = _x.size();
};


double LMopt::Solver()
{
    //start from x0;
    varX.resize(Nx);
    for(int i=0;i<Nx;i++)
        varX(i) = init_x[i];

    //init
    varF = GetFx(varX);
    Jcob = GetJcob(varX);
    MatrixXd Hes;
    VectorXd g;
    Hes = Jcob.transpose()*Jcob;
    g = Jcob.transpose()*varF;
    //L(h) =
    int k_itr = 0;
    double v = 2;
    double tou = 1;
    double u = tou*Jcob.maxCoeff();
    double e2 = 0.000001;
    double e1 = 0.000001;
    bool is_found = g.norm() < e2;

    while(!is_found && k_itr<1000)
    {
        //cout<< " itr= "<< k_itr <<" u=" << u << " v= "<< v
        //    << " Fx= "<< varF.squaredNorm()/2 << endl;
        //cout<< "x=" <<varX.transpose()<<endl;
        //cout<< "f= "<< varF.transpose() <<endl;
        //cout<< "detH=" << Hes.determinant()<<endl;
        k_itr++;
        MatrixXd revM;
        revM.setIdentity(Nx, Nx);
        revM = u*revM + Hes;
        ///
        //revM = GetHess(varX, varF);
        VectorXd h = revM.ldlt().solve(-g);
        //cout<<"g=" <<g.transpose()<<endl;
        //cout<< " |h|= "<<h.norm()<<endl;
        if( h.norm() < e2*(varX.norm() + e2))
            is_found = true;
        else
        {
            VectorXd x_new = varX + h;
            double Fx = varF.squaredNorm()/2.0;
            double L0 = Fx;
            double Lh = L0 + h.dot(g) + (Jcob*h).dot(Jcob*h)/2.0;
            double Fxh = GetFx(x_new).squaredNorm()/2.0;
            double rho = (Fx - Fxh)/(L0 - Lh);
            //cout<< " Fx= "<< Fx << " L0= " << L0 << " Lh= " <<Lh << " Fxh= "<< Fxh << " rho= " <<rho <<endl;

            if(rho > 0)
            {
                varX = x_new;
                varF = GetFx(varX);
                Jcob = GetJcob(varX);
                Hes = Jcob.transpose()*Jcob;
                g = Jcob.transpose()*varF;
                is_found = (g.norm() < e1) || (varF.norm() < e1);
                u = u * _MAX(1/3.0, 1-pow(2*rho-1,3));
                v = 2;
            }
            else
            {
                u = u*v;
                v = 2*v;
            }
        }
    }
    if(is_found)
    {
        //cout<<"find at," <<" itr = " << k_itr << " |F| = "<< varF.norm() << endl;
        //cout<< "x=" <<varX.transpose();
        //cout<< "f= "<< varF.transpose() <<endl;
        //cout<< "detH=" << Hes.determinant()<<endl;
        res_x.clear();
        for(int i=0;i<Nx;i++)
            res_x.push_back(varX(i));
    }
    else
    {
        cout<< " Not found!" << " |F| = " << varF.norm() << endl;
        res_x = init_x;
    }

    return varF.norm();
}

VectorXd LMopt::GetFx(VectorXd _x)
{
    // example 1
    // f1 = x1 - x2/2;
    // f2 = x1 + x2 - 1;
    /*
    VectorXd fx;
    fx.resize(_x.size());
    fx(0) = _x(0) - _x(1)/2;
    fx(1) = _x(0) + _x(1) - 1;
     */

    //example 2
    /*
    VectorXd fx;
    fx.resize(_x.size());
    fx(0) = cos(_x(0)) - 2*sin(_x(1)) + pow(_x(2),2);
    fx(1) = sin(_x(0)) + 2*cos(_x(1)) - _x(2);
    fx(2) = cos(_x(0)) - pow(_x(1),2) + _x(2);
     */

    // four bar linkage
    VectorXd fx;
    fx.resize(_x.size());
    fx(0) = cos(_x(2))*pinPos[0].x() - sin(_x(2))*pinPos[0].y() + _x(0) - pinPos[0].x();
    fx(1) = sin(_x(2))*pinPos[0].x() + cos(_x(2))*pinPos[0].y() + _x(1) - pinPos[0].y();

    fx(2) = cos(_x(2))*pinPos[1].x() - sin(_x(2))*pinPos[1].y() + _x(0)
            - (cos(_x(5))*pinPos[1].x() - sin(_x(5))*pinPos[1].y()+ _x(3) );
    fx(3) = sin(_x(2))*pinPos[1].x() + cos(_x(2))*pinPos[1].y() + _x(1)
            - (sin(_x(5))*pinPos[1].x() + cos(_x(5))*pinPos[1].y()+ _x(4) );

    fx(4) = cos(_x(5))*pinPos[2].x() - sin(_x(5))*pinPos[2].y() + _x(3)
            - (cos(_x(8))*pinPos[2].x() - sin(_x(8))*pinPos[2].y()+ _x(6) );
    fx(5) = sin(_x(5))*pinPos[2].x() + cos(_x(5))*pinPos[2].y() + _x(4)
            - (sin(_x(8))*pinPos[2].x() + cos(_x(8))*pinPos[2].y()+ _x(7) );

    fx(6) = cos(_x(8))*pinPos[3].x() - sin(_x(8))*pinPos[3].y() + _x(6) - pinPos[3].x();
    fx(7) = sin(_x(8))*pinPos[3].x() + cos(_x(8))*pinPos[3].y() + _x(7) - pinPos[3].y();

    fx(8) = _x(2) - active_x[0];
    return fx;
}

MatrixXd LMopt::GetJcob(VectorXd _x)
{
    //example 1
    // J = (1, -0.5 ;
    //      1, 1;)
    /*
    MatrixXd Jx;
    Jx.resize(_x.size(), _x.size());
    Jx(0,0) = 1;
    Jx(0,1) = -0.5;
    Jx(1,0) = 1;
    Jx(1,1) = 1;
*/
    //example 2
    /*
    MatrixXd Jx;
    Jx.resize(_x.size(), _x.size());
    Jx(0,0) = -sin(_x(0));
    Jx(0,1) = -2*cos(_x(1));
    Jx(0,2) = 2*_x(2);
    Jx(1,0) = cos(_x(0));
    Jx(1,1) = -2*sin(_x(1));
    Jx(1,2) = -1;
    Jx(2,0) = -sin(_x(0));
    Jx(2,1) = -2*_x(1);
    Jx(2,2) = 1;
*/
    MatrixXd Jx;
    Jx.resize(_x.size(), _x.size());
    Jx.setZero();
    Jx(0,0) = 1;
    Jx(0, 2) = -sin(_x(2))*pinPos[0].x() - cos(_x(2))*pinPos[0].y();
    Jx(1,1) = 1;
    Jx(1,2) = cos(_x(2))*pinPos[0].x() - sin(_x(2))*pinPos[0].y();

    Jx(2, 0) = 1;
    Jx(2,2) = -sin(_x(2))*pinPos[1].x() - cos(_x(2))*pinPos[1].y();
    Jx(2,3) = -1;
    Jx(2,5) = -(-sin(_x(5))*pinPos[1].x() - cos(_x(5))*pinPos[1].y());
    Jx(3, 1) = 1;
    Jx(3,2) = cos(_x(2))*pinPos[1].x() - sin(_x(2))*pinPos[1].y();
    Jx(3,4) = -1;
    Jx(3,5) =-(cos(_x(5))*pinPos[1].x() - sin(_x(5))*pinPos[1].y());

    Jx(4, 3) = 1;
    Jx(4,5) = -sin(_x(5))*pinPos[2].x() - cos(_x(5))*pinPos[2].y();
    Jx(4,6) = -1;
    Jx(4,8) = -(-sin(_x(8))*pinPos[2].x() - cos(_x(8))*pinPos[2].y());
    Jx(5, 4) = 1;
    Jx(5,5) = cos(_x(5))*pinPos[2].x() - sin(_x(5))*pinPos[2].y();
    Jx(5,7) = -1;
    Jx(5,8) =-(cos(_x(8))*pinPos[2].x() - sin(_x(8))*pinPos[2].y());

    Jx(6,6) = 1;
    Jx(6, 8) = -sin(_x(8))*pinPos[3].x() - cos(_x(8))*pinPos[3].y();
    Jx(7,7) = 1;
    Jx(6,8) = cos(_x(8))*pinPos[3].x() - sin(_x(8))*pinPos[3].y();

    Jx(8,2) = 1;

    return Jx;
}

MatrixXd LMopt::GetHess(VectorXd _x, VectorXd _f)
{
    double f0_22 = -cos(_x(2))*pinPos[0].x() + sin(_x(2))*pinPos[0].y();
    double f1_22 = -sin(_x(2))*pinPos[0].x() - cos(_x(2))*pinPos[0].y();

    double f2_22 = -cos(_x(2))*pinPos[1].x() + sin(_x(2))*pinPos[1].y();
    double f2_55 = cos(_x(5))*pinPos[1].x() - sin(_x(5))*pinPos[1].y();
    double f3_22 = sin(_x(2))*pinPos[1].x() + cos(_x(2))*pinPos[1].y();
    double f3_55 = sin(_x(5))*pinPos[1].x() + cos(_x(5))*pinPos[1].y();

    double f4_55 = -cos(_x(5))*pinPos[2].x() + sin(_x(5))*pinPos[2].y();
    double f4_88 = cos(_x(8))*pinPos[2].x() - sin(_x(8))*pinPos[2].y();
    double f5_55 = -sin(_x(5))*pinPos[2].x() - cos(_x(5))*pinPos[2].y();
    double f5_88 = sin(_x(8))*pinPos[2].x() + cos(_x(8))*pinPos[2].y();

    double f6_88 = -cos(_x(8))*pinPos[3].x() + sin(_x(8))*pinPos[3].y();
    double f7_88 = -sin(_x(8))*pinPos[3].x() - cos(_x(8))*pinPos[3].y();

    MatrixXd Hess;
    Hess.resize(Nx, Nx);
    Hess.setZero();
    Hess(2,2) = f0_22*_f(0) + f1_22*_f(1) + f2_22*_f(2) + f3_22*_f(3);
    Hess(5,5) = f2_55*_f(2) + f3_55*_f(3) + f4_55*_f(4) + f5_55*_f(5);
    Hess(8,8) = f4_88*_f(4) + f5_88*_f(5) + f6_88*_f(6) + f7_88*_f(7);

    Hess = Hess + GetJcob(_x);
    return Hess;
}
