
#include "LMopt.h"
#include "Utility/HelpFunc.h"

int LMFunctor::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
{
    fvec.resize(m);
    vector<double> kepa = GetKepa(x);
    for(int i=0;i<m-1-keys;i++)
        fvec(i) = kepa[i];
    fvec(m-1-keys) = 0;
    for(int i=1;i<=n+1;i++)
        if(sf->Ff[i] <= sf->Ff[i-1])
            fvec(m-1-keys) += exp(1000*pow(sf->Ff[i-1] - sf->Ff[i], 2))-1;
    for(int i=0; i<keys;i++)
        fvec(i + m-keys) = pow(sf->GetFx(Key[i].first) - Key[i].second ,2)*100;

    return 0;
}

int LMFunctor::df(const VectorXd &x, MatrixXd &fjac) const {

    fjac.resize(m, n);
    fjac.setZero();

    VectorXd epc;
    double tol = 0.00001;
    epc.resize(n);

    for(int itr=0;itr<n;itr++)
    {
        epc.setZero();
        epc(itr) = tol;
        vector<double> Mkepa = GetKepa(x+epc);
        double MF = 0;
        for(int i=1;i<=n+1;i++)
            if(sf->Ff[i] <= sf->Ff[i-1])
                MF += exp(1000*pow(sf->Ff[i-1] - sf->Ff[i], 2))-1;
        vector<double> Mkey(keys);
        for(int i=0; i<keys;i++)
            Mkey[i] = pow(sf->GetFx(Key[i].first) - Key[i].second ,2)*100;

        vector<double> mkepa = GetKepa(x-epc);
        double mf = 0;
        for(int i=1;i<=n+1;i++)
            if(sf->Ff[i] <= sf->Ff[i-1])
                mf += exp(1000*pow(sf->Ff[i-1] - sf->Ff[i], 2))-1;
        vector<double> mkey(keys);
        for(int i=0; i<keys;i++)
            mkey[i] = pow(sf->GetFx(Key[i].first) - Key[i].second ,2)*100;


        for(int i=0;i<m-1-keys;i++)
            fjac(i, itr) = (Mkepa[i] - mkepa[i])/(2*tol);
        fjac(m-1-keys, itr) = (MF - mf)/(2*tol);
        for(int i=0;i<keys;i++)
            fjac(i+m-keys, itr) = (Mkey[i] - mkey[i])/(2*tol);
    }

    return 0;
}

vector<double> LMFunctor::GetKepa(const VectorXd &x) const{

    sf->SameSpline(x);
    vector<Vector3d> pC(m-1-keys);
    for(int i=0;i<m-1-keys;i++) {
        double t = i*2*M_PI/(m-1-keys);
        double s = sf->GetFx(t)/2/M_PI;
        Vector2d tmp = bsf->SplinePoint(s);
        Matrix4d rM = GetRotationMatrix(Vector3d(0,0,t))*
                GetRotationMatrix(tmp.x(), Vector3d(0, 1, 0), folCen) *
                GetRotationMatrix(tmp.y(), Vector3d(0, 0, 1), folCen);
        MultiplyPoint(folJoint, rM, pC[i]);
    }

    vector<double> kepa(m-1-keys);
    for(int i=0;i<m-1-keys;i++)
    {
        int pre = (i-1) >= 0 ? i-1 : m-2-keys;
        int net = (i+1) >= m-1-keys ? 0 : i+1;
        Vector3d A = pC[pre] - pC[i];
        Vector3d B = pC[net] - pC[i];
        double r = (A-B).norm() * A.norm() * B.norm() / (A.cross(B)).norm();
        if(r > 15)
            kepa[i] = 15/r;
        else
            kepa[i] = exp(15/r - 1);
    }

    return kepa;
}