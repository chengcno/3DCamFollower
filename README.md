# 3D CamFollower Program

[![Status](https://img.shields.io/badge/status-active-success.svg)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](/LICENSE)

![](https://github.com/chengcno/SpatialLinkageCam/blob/main/doc/teaser.png)

This repo is an implementation of [Spatial-Temporal Motion Control via Composite Cam-follower Mechanisms](https://sutd-cgl.github.io/supp/Publication/projects/2021-SIGAsia-3DCamMech/index.html) [Cheng et al. 2021].
If you have any problems when using this code, you can contact me any time through chengyj@mail.ustc.edu.cn

If you make use of this repo in your scientific work, please cite our paper. For your convenience,
you can use the following bibtex snippet:

    @article {Cheng-2021-3DCamMech,
    author  = {Yingjie Cheng and Yucheng Sun and Peng Song and Ligang Liu},
    title   = {Spatial-Temporal Motion Control via Composite Cam-follower Mechanisms},
    journal = {ACM Transactions on Graphics (SIGGRAPH Asia 2021)},
    volume  = {40},
    number  = {6},
    pages   = {270:1 -- 270:15},
    year    = {2021}}

## Table of Contents
- [About](#about)
- [Getting Started](#getting_started)
- [GUI Interface](#usage)
- [Create a Cam-Follower Mechanism](#create_mech)
- [Acknowledgments](#acknowledgement)

## About <a name = "about"></a>
This repo presents a computational approach of designing a cam-follower mechanism that can control their spatial-temporal motions to exactly follow trajectories and timings specified by users, driven by a single actuator.
We implemented our computational design tool in C++ and `libigl` [Jacobson et al. 2018] on a desktop computer with 3.6 GHz 8-Core Intel processor and 16 GB RAM.

## Getting Started <a name = "getting_started"></a>
Our code can be ran on MacOS and Unbuntu (Linux) system. First clone the repository, run CMake to generate Makefiles or CMake/Visual Studio project files, and the rest should just work automatically.


### Compilation

- **MacOS and Ubuntu(Linux)**:

```
cd [current folder path]
mkdir build
cd build
cmake ..
make
```
It's better to build it with `release` because of CGAL.
This should find and build the dependencies and create a `3DCams_Main` binary.

- **Windows**: currently unavailable.


## GUI Interface <a name = "usage"></a>
The control panel is shown below. There are 3 components in the control panel: 
**Operation Control**, **Motion Control**, **Render Control**.
![](https://github.com/chengcno/SpatialLinkageCam/blob/main/doc/UI.png)

- ### Operation Control

  `Read Input Curve`  Read an input curve from folder `/data/inputCurves`

  `Optimization` Optimize a mechanism to realize the input curve.

  `Read CamFollower` Our program can read `mats.dat` files from folder `/data/resultsData`

- ### Operation Control

  `Stop Motion`/`Restart Motion` Interruption or continuation of movement

  `Motion Speed` Adjust the motion velocity

- ### Render Control
  Control the object visualization state.

## Create a Cam-Follower Mechanism <a name = "create_mech"></a>
These instructions give an example to you of how to use our code to generate a cam-follower mechanism by yourself.

### Step 1: import an input curve
Import any input curve file by clicking `Read Input Curve` button.
![](https://github.com/chengcno/SpatialLinkageCam/blob/main/doc/inputCurve.png)

### Step 2: optimize a mechanism 
Click the `Optimization` button. It would cost about ~10sec.
![](https://github.com/chengcno/SpatialLinkageCam/blob/main/doc/optModel.png)

### Step 3: control the movement
Use `Restart Motion`,`Stop Motion`,`Motion Speed` to control the movement. Use mouse to adjust the camera view.
![](https://github.com/chengcno/SpatialLinkageCam/blob/main/doc/motion.png)

### Step 4: Show results demo
Import one results demo from file `data/resultsData/ * /mats.txt` by clicking `Read CamFollower` button.
Here is a view of drawing cats.
![](https://github.com/chengcno/SpatialLinkageCam/blob/main/doc/motion.png)


## Acknowledgements <a name = "acknowledgement"></a>
We thank the reviewers for their valuable comments. 
This work was supported by the SUTD Start-up Research Grant (Number: SRG ISTD 2019 148), 
the National Natural Science Foundation of China (62025207), and Zhejiang Lab (NO. 2019NB0AB03).