//
// Created by cheng on 27/12/21.
// Modified by CHENG 2022/Aug/08
//

#ifndef SPATIALLINKAGES_READCURVE_H
#define SPATIALLINKAGES_READCURVE_H

#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include <Mesh/Mesh.h>
#include "BSplineFit.h"
#include "Utility/HelpTypedef.h"
#include "CamPart/Cam2R.h"
#include "CamPart/CamC.h"

using namespace std;
using namespace Eigen;

enum SurfaceType {R, T, RR, TT, RT};

class ReadCurve{

public:
    ReadCurve();
    ~ReadCurve(){};

public:
    int pNum;
    MatrixXd posM;

    Vector3d box_min, box_max;
    Matrix4d TranM;

    Mesh     *curMesh;
    Mesh     *supMesh;
    SurfaceType surface_type;

    Vector3d folCen;
    vector<double> folAlphas, folBetas;
    Cam2R* cam2R;
    CamC* camRT;

public:
    void readFile(string stf);
    void funcGene(int type);

    void scaleTo(double _len);   //set origin center and scale
    void rotateTo(Vector3d angleE);     // rotate around 0
    void rotateTo(Matrix4d rM);
    void moveTo(Vector3d startP);       // move to P, and start from P

    void generateMesh(double _radius);

    void readMesh();
    void AppendReadMechs(iglViewer &viewer);
    void ClearViewerList(iglViewer &viewer);
    void AppendOptMechs(iglViewer &viewer);
    void UpdateMotion(iglViewer &viewer, int frame);

    void CreateCamFollower();

private:
    void getMinMax();
    int  closeID();
    void ReplaceOrder();
    void placeCen();

    void ReSample();
    void Smooth();
    void SmileCurve();

    vector<double> Hermite(double x, double y, double dx, double dy, int size);
    vector<Vector3d> HermiteCur(Vector3d X, Vector3d Y, Vector3d dX, Vector3d dY, int size);

    Vector3d Bezier(double t, vector<Vector3d> CtrPts);

private:
    double length;
    int _type;

    BSplineFit* _bsF;
    double radiusR;
};


#endif //SPATIALLINKAGES_READCURVE_H
