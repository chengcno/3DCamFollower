//
// Created by temp on 7/7/21.
//

#ifndef CONVEXHULL_H
#define CONVEXHULL_H

#include "Mesh/Mesh.h"
#include <igl/copyleft/cgal/convex_hull.h>
#include <Eigen/Eigen>

using namespace Eigen;

class ConVexHull {
public:
    ConVexHull();
    ConVexHull(Mesh* _modelM);
    ~ConVexHull();

public:
    Mesh* modelM;
    Mesh* convexM;

    Mesh* GetConvexHull(Mesh* _modelM);
    Mesh* GetConvexHull(MatrixXd _verM);
    Mesh* GetConvexHull(vector<Vector3d> verList);

    void GetConvexHull_2D(MatrixXd _verM, MatrixXd& re_verM);

};


#endif //CONVEXHULL_H
