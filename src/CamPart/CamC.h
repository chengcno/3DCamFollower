//
// Created by Yingjie CHNEG
// 20/Dec/2020
//

#ifndef SPATIALLINKAGES_CAMC_H
#define SPATIALLINKAGES_CAMC_H


#include <Eigen/Eigen>
#include <vector>
#include "Mesh/Mesh.h"

using namespace Eigen;
using namespace std;

class CamC {
public:
    CamC();
    //CamC(vector<Vector3d> _PitchCurve, Vector3d _follCen);
    ~CamC(){};

public:
    Mesh*               camMesh;
    Mesh*               folMesh;
    double              rCFball, rFol;
    vector<Vector3d>    PitchCurve;
    Vector3d            folCen;
    Vector3d            folAxis;
    vector<double>      folTheta, folTrans;
    Vector3d            pJoint_0, pAnkle_0;

    //Matrix4d            worldMat;
    int                 numN, numM;

private:
    vector<pair<double,Vector3d>>   FollOrient;          // Follower orient
    vector<Quaterniond>             outQuaters;          // follower quats
    vector<double>                  outTrans;            // follower translate

public:
    void CreateCamMesh();
    void CreateFolMesh();
    void ComputeEMechMotion(Matrix4d driMat, Matrix4d &folMat);

private:
    void GetPitchCurve();
    void ComputeFollowOrient();
    // design
    //void GetDesignedCurve();


};


#endif //SPATIALLINKAGES_CAMC_H
