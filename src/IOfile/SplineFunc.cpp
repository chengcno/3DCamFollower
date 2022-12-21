//
// Created by Yingjie CHNEG
// 20/Dec/2020
//

#include "SplineFunc.h"
#include "BSplineFit.h"
#include <iostream>
#include <ctime>


SplineFunc::SplineFunc(vector<double> _fx)
{
    fx = _fx;

    Ff.push_back(0);
    for(auto f : fx)
        Ff.push_back(f);
    Ff.push_back(2.0*M_PI);

    nSize = _fx.size() + 1;
    theta = 2.0*M_PI/nSize;

    DecompressMatrix();
    ComputeCubicSpline();
}

SplineFunc::SplineFunc(const double *S, int size_s)
{
    nSize = size_s+1;
    theta = 2.0*M_PI/nSize;

    Ff.push_back(0);
    for(int i=0;i<size_s;i++)
    {
        fx.push_back(S[i]);
        Ff.push_back(S[i]);
    }
    Ff.push_back(2*M_PI);

    DecompressMatrix();
    ComputeCubicSpline();
}

SplineFunc::~SplineFunc() {}

void SplineFunc::GetIdrel(int &id, double &rel, double t)
{
    id = floor(t/theta);
    rel = (t - id*theta);
    if(id == nSize) {id = 0;}
    if(id > nSize || id < 0 || rel<0) printf("s=s(t) value t is wrong!\n");
}

int SplineFunc::delta(int i, int j)
{
    return i == j ? 1 : 0;
}

double SplineFunc::GetFx(double t)
{
    int id = floor(t/theta);
    double rel = (t - id*theta);
    if(id == nSize) id = 0;
    if(id > nSize || id < 0)
    {
        printf("s=s(t) value t is wrong!\n");
        cout<<"t="<< t<< " id="<<id << " rel="<<rel <<" theta = "<<theta<<endl;
    }

    double SplineY;
    SplineY = a[id] + b[id]*rel + c[id]*pow(rel, 2) +d[id]*pow(rel,3);
    return SplineY;
}

double SplineFunc::GetDFx(double t)
{
    int id = floor(t/theta);
    double rel = (t - id*theta);
    if(id == nSize) {id = 0;}
    if(id > nSize || id < 0 || rel<0) printf("s=s(t) value t is wrong!\n");

    double DSplineY;
    DSplineY = b[id] + 2*c[id]*pow(rel, 1) + 3*d[id]*pow(rel,2);
    return DSplineY;
}

double SplineFunc::GetDDFx(double t)
{
    int id = floor(t/theta);
    double rel = (t - id*theta);
    if(id == nSize) {id = 0;}
    if(id > nSize || id < 0 || rel<0) printf("s=s(t) value t is wrong!\n");

    double DDSplineY;
    DDSplineY = 2*c[id] + 6*d[id]*pow(rel,1);
    return DDSplineY;
}

double SplineFunc::GetDDDFx(double t)
{
    int id = floor(t/theta);
    double rel = (t - id*theta);
    if(id == nSize) {id = 0;}
    if(id > nSize || id < 0 || rel<0) printf("s=s(t) value t is wrong!\n");

    double D3SplineY;
    D3SplineY = 6*d[id];
    return D3SplineY;
}

void SplineFunc::DecompressMatrix()
{
    SparseMatrix<double> triM(nSize + 1, nSize + 1);
    triM.reserve(3 * nSize + 3);

    triM.coeffRef(0, 0) = 2;
    triM.coeffRef(0, 1) = 1;
    triM.coeffRef(0, nSize - 1) = 1;
    triM.coeffRef(0, nSize) = 2;
    for (int i = 1; i < nSize; i++) {
        int bef, nex;
        bef = i - 1;
        nex = i + 1;
        triM.coeffRef(i, bef) = 1;
        triM.coeffRef(i, i) = 4;
        triM.coeffRef(i, nex) = 1;
    }
    triM.coeffRef(nSize, 0) = 1;
    triM.coeffRef(nSize, nSize) = -1;

    triM.makeCompressed();

    solverTriM.compute(triM);
}

void SplineFunc::ComputeCubicSpline() {
    /*
    SparseMatrix<double> triM(nSize + 1, nSize + 1);
    triM.reserve(3 * nSize + 3);

    triM.coeffRef(0, 0) = 2;
    triM.coeffRef(0, 1) = 1;
    triM.coeffRef(0, nSize - 1) = 1;
    triM.coeffRef(0, nSize) = 2;
    for (int i = 1; i < nSize; i++) {
        int bef, nex;
        bef = i - 1;
        nex = i + 1;
        triM.coeffRef(i, bef) = 1;
        triM.coeffRef(i, i) = 4;
        triM.coeffRef(i, nex) = 1;
    }
    triM.coeffRef(nSize, 0) = 1;
    triM.coeffRef(nSize, nSize) = -1;

    triM.makeCompressed();
*/
    /*
    //clock_t st,et;
    //cout<<"start LU "<<endl;
    //st = clock();
    //chrono::steady_clock::time_point time_start=chrono::steady_clock::now();
    //SparseLU<SparseMatrix<double>> solverTriM(triM);
    //et = clock();
    //chrono::steady_clock::time_point time_end=chrono::steady_clock::now();
    //chrono::duration<double> time_used=chrono::duration_cast<chrono::duration<double>>(time_end-time_start);
    //cout<<"finish LU "<<"duration= "<<1000*(et -st)/(double)CLOCKS_PER_SEC<<endl;
    //SimplicialCholesky<SparseMatrix<double>> solverTriM(triM);

    //cout<<"start inverse"<<endl;
    //st = clock();
     */
    ///inverse matrix for gradient
    /*
    SparseMatrix<double> I(nSize+1,nSize+1);
    I.setIdentity();
    MatrixXd M_inv = solverTriM.solve(I);
    //et = clock();
    //cout<<"finish inverse "<<endl;

    /// for Y matrix
    SparseMatrix<double> Ym(nSize+1,nSize+1);
    Ym.reserve(3 * nSize + 3);
    Ym.coeffRef(0,0) = -1;
    Ym.coeffRef(0,1) = 1;
    Ym.coeffRef(0,nSize-1) = 1;
    Ym.coeffRef(0,nSize) = -1;
    for(int i=1;i<nSize;i++)
    {
        Ym.coeffRef(i,i-1) = 1;
        Ym.coeffRef(i,i) = -2;
        Ym.coeffRef(i,i+1) = 1;
    }
    Ym = Ym*6/pow(theta,2);
    Ym.makeCompressed();
    M_Y = M_inv*Ym;
*/
    //cout<<triM<<endl;
    //cout<<M_inv<<endl;
    //cout<<Ym<<endl;
    //cout<<M_Y<<endl;

    VectorXd Vy;
    Vy.resize(nSize + 1);
    Vy(0) = ((Ff[1] - Ff[0]) - (Ff[nSize] - Ff[nSize - 1])) * 6 / pow(theta, 2);

    for (int i = 1; i < nSize; i++) {
        Vy(i) = (Ff[i + 1] + Ff[i - 1] - 2 * Ff[i]) * 6 / pow(theta, 2);
    }
    Vy(nSize) = 0;

    VectorXd Res;
    Res = solverTriM.solve(Vy);

    m.resize(nSize+1);
    a.resize(nSize);
    b.resize(nSize);
    c.resize(nSize);
    d.resize(nSize);
    for (int i = 0; i <= nSize; i++) {
        m[i] = Res(i);
    }
    for (int i = 0; i < nSize; i++) {
        a[i] = Ff[i];
        b[i] = (Ff[i + 1] - Ff[i]) / theta - theta * m[i] / 2 - (m[i + 1] - m[i]) * theta / 6;
        c[i] = m[i] / 2;
        d[i] = (m[i + 1] - m[i]) / 6 / theta;
    }

    //printf("finish solve s=s(t) spline\n");
    //cout<<"res =\n"<<Res<<endl;
    //cout<<theta<<endl;
    //cout << "m0 =" << m[0] << " m_n = " << m[nSize] << endl;
    //cout << "df0=" << b[0] << " df_n= " << b[nSize - 1] + 2 * c[nSize - 1] * theta + 3 * d[nSize - 1] * pow(theta, 2)
    //     << endl;
    //cout << "df0=" << (Ff[1]-Ff[0])/theta - theta*m[0]/3 - theta*m[1]/6 << " df_n= " << (Ff[nSize]-Ff[nSize-1])/theta + theta*m[nSize]/3 + theta*m[nSize-1]/6
    //     << endl;
    //cout<<"df0 = "<<GetDFx(0)<<" fd_n ="<<GetDFx(2*M_PI - 0.00001)<<endl;

    //VectorXd err = triM*Res - Vy;
    //cout<<err<<endl;
}

vector<double> SplineFunc::partFx(double t)
{
    vector<double> partS;
    double partSj;

    int id;
    double rel;
    GetIdrel(id, rel, t);

    for(int j=1;j<nSize;j++)
    {
        partSj = delta(id,j) + ((delta(id+1,j) - delta(id,j))/theta - theta*M_Y(id,j)/2 - theta*(M_Y(id+1,j) - M_Y(id,j))/6)*rel
                + M_Y(id,j)*pow(rel,2)/2 + (M_Y(id+1,j)-M_Y(id,j))/6.0/theta * pow(rel,3);
        partS.push_back(partSj);
    }

    return partS;
}

vector<double> SplineFunc::partDFx(double t)
{
    vector<double> partDS;
    double partDSj;

    int id;
    double rel;
    GetIdrel(id, rel, t);

    for(int j=1;j<nSize;j++)
    {
        partDSj = ((delta(id+1,j) - delta(id,j))/theta - theta*M_Y(id,j)/2 - theta*(M_Y(id+1,j) - M_Y(id,j))/6)
                 + M_Y(id,j)*rel + (M_Y(id+1,j)-M_Y(id,j))/6.0/theta * 3*pow(rel,2);
        partDS.push_back(partDSj);
    }

    return partDS;
}

vector<double> SplineFunc::partDDFx(double t)
{
    vector<double> partDDS;
    double partDDSj;

    int id;
    double rel;
    GetIdrel(id, rel, t);

    for(int j=1;j<nSize;j++)
    {
        partDDSj = M_Y(id,j)+ (M_Y(id+1,j)-M_Y(id,j))/theta * rel;
        partDDS.push_back(partDDSj);
    }

    return partDDS;
}

vector<double> SplineFunc::partDDDFx(double t)
{
    vector<double> partDDDS;
    double partDDDSj;

    int id;
    double rel;
    GetIdrel(id, rel, t);

    for(int j=1;j<nSize;j++)
    {
        partDDDSj = (M_Y(id+1,j)-M_Y(id,j))/theta;
        partDDDS.push_back(partDDDSj);
    }

    return partDDDS;
}


bool SplineFunc::SameSpline(const VectorXd &x)
{
    bool is_same = true;

    for(int i=0;i<nSize-1;i++)
    {
        if(abs(x[i] - Ff[i-1]) > 0.0000001)
        {
            is_same = false;
            //cout<<"update s(t)"<<endl;
            break;
        }
    }

    if(!is_same)
    {
        for(int i=1;i<nSize;i++)
            Ff[i] = x[i-1];
        ComputeCubicSpline();
    }

    return is_same;
}