# 3D Cam-Follower

A Framework project about 3D Cam-Follower creating, showing and operating.

## Dependencies

The only dependencies are stl, eigen, pugixml, [libigl](http://libigl.github.io/libigl/), CGAL and
the dependencies of the `igl::opengl::glfw::Viewer`.

 libigl and pugixml are added to this git repo as submodule

## Compile

It's better to build it with `release` because of CGAL.

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `3Dcams_bin` binary.

## Run
You would see a ground and 3 axis.

