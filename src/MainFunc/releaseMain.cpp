
///////////////////////////////////////////////////////////////
//
// main.cpp
//
//   3D Cam Program
//
// by Yingjie Cheng and Peng Song
//
// 06/Aug/2022
//
//
///////////////////////////////////////////////////////////////


#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/per_corner_normals.h>
#include "Utility/HelpTypedef.h"
#include "Utility/HelpDefine.h"
#include "Utility/HelpFunc.h"
#include "Render/libigl_Render.h"
#include "Render/libigl_UI.h"
#include "Mesh/ConVexHull.h"
#include <chrono>
#include "IOfile/ReadCurve.h"
#include "IOfile/IOMechs.h"

///////////////////////////////////////////////////////////////
// Global Variables
///////////////////////////////////////////////////////////////

// the viewer
iglViewer viewer;

// a menu plugin
igl::opengl::glfw::imgui::ImGuiMenu menu;

int windowWidth;
int windowHeight;
static bool isFirstCall;

int AnimateSpeed;

bool visibleMechanism;
bool visibleAxes;
bool visibleGround;

bool is_restart;

libigl_Render myRender;

int frame;

bool is_ReadCurve;
bool is_ReadMechs;
bool is_SaveMechs;
bool is_Optimize;
bool is_FileModel;
bool is_OptModel;

string inputFileName;


///////////////////////////////////////////////////////////////
// Function Declarations
///////////////////////////////////////////////////////////////

void InitSetting();
void InitViewer();


///=========================================================================================///
///                                       Initialization
///=========================================================================================///

void InitSetting()
{
    windowWidth       =  1600;
    windowHeight      =  1000;

    isFirstCall       =  true;

    AnimateSpeed       =  1;

    visibleMechanism  =  true;
    visibleAxes       =  true;
    visibleGround     =  true;

    is_restart = false;
    frame = 0;

    is_ReadCurve = false;
    is_ReadMechs = false;
    is_SaveMechs = false;
    is_Optimize = false;
    is_FileModel = false;
    is_OptModel = false;
}

void InitViewer()
{
    ///set animation
    viewer.core().animation_max_fps = 30.;
    viewer.core().is_animating = false;
    viewer.core().background_color = RowVector4f (1,1,1,0);

    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
}


///=========================================================================================///
///                                       Interaction
///=========================================================================================///

bool key_down(iglViewer &Viewer, unsigned char key, int modifier)
{
    if (key == ' ')
    {
        Viewer.core().is_animating = !Viewer.core().is_animating;
    }
    return false;
}


///=========================================================================================///
///                                       Main Function
///=========================================================================================///

int main(int argc, char *argv[])
{
    setViewerUI(viewer);
    InitViewer();
    InitSetting();

    auto _file = new ReadCurve();
    auto ioMech = new IOMechs();

    vector<igl::opengl::ViewerData> mydataList;

    myRender.RenderScene( viewer, mydataList);

    viewer.callback_pre_draw =
            [&](igl::opengl::glfw::Viewer &)
            {
                myRender.ShowMyModels( viewer, visibleMechanism );
                myRender.ShowGround( viewer, visibleGround );
                myRender.ShowAxes(viewer, visibleAxes );

                if(viewer.core().is_animating)
                {
                    frame += AnimateSpeed;
                    if(is_FileModel) ioMech->UpdateMotion(viewer, frame);
                    else if(is_OptModel)
                    {
                        _file->UpdateMotion(viewer, frame);
                    }
                }

                if(is_ReadCurve)
                {
                    cout<< inputFileName << endl;
                    _file->readFile(inputFileName);
                    _file->ClearViewerList(viewer);
                    _file->AppendReadMechs(viewer);
                    is_ReadCurve = false;
                    is_OptModel = false;
                    is_FileModel = false;
                }

                if(is_Optimize)
                {
                    _file->CreateCamFollower();
                    _file->AppendOptMechs(viewer);
                    is_FileModel = false;
                    is_OptModel = true;
                    is_Optimize = false;
                    frame = 0;
                }

                if(is_ReadMechs)
                {
                    cout<< inputFileName << endl;
                    ioMech->ClearViewerList(viewer);
                    ioMech->ReadInputMechs(inputFileName);
                    ioMech->AppendReadMechs(viewer);
                    is_ReadMechs = false;
                    is_FileModel = true;
                    is_OptModel = false;
                    frame = 1;
                    cout<< "finish read cam-followers" <<endl;
                    ioMech->UpdateMotion(viewer, frame);
                }

                if(isFirstCall)
                {
                    myRender.SetCamera( viewer, 0.1, Vector3f (30, 30, 30) );
                    isFirstCall = false;
                }

                return 0;
            };

    viewer.callback_key_down = &key_down;

    viewer.launch(true,false,"3DCamFollower", windowWidth, windowHeight);

    return 1;
}

