#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <qpushbutton.h>
#include <QtWidgets>
#include <QTimer>
#include "GlInterface/glwidget.h"
#include "mymesh/geo_sur.h"
#include "mymesh/geo_curv.h"
#include "mymesh/contour.h"
#include "mymesh/geo_pc.h"
#include "mymesh/a_multisur.h"


QT_BEGIN_NAMESPACE
class QLabel;
class QMenu;
class QScrollArea;
class QSlider;
class QCheckBox;
class QGroupBox;
class QRadioButton;
QT_END_NAMESPACE
namespace Ui {
class MainWindow;
}






class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void GetArg(int argc,char** argv);
private:
    int argc;
    char** argv;

private:
    Ui::MainWindow *ui;

    QFileDialog *filediag = NULL;


    QTimer timers_forvideo;
    QTimer timers_forcontrol;
    int accTimeout;
    int accTimeoutControl;
    int endTimeoutControl;
    vector<QImage>videoStack;

    //VideoWidgetSurface videoWidget;


    GLWidget *glinterface;

    CrossSections *hidden_crossSection;
    n_rf::Surface *hiddenSurface_forCsView;
    n_rf::Surface *hiddenSurface_forCellsView;
    n_rf::Multi_PointCloud *displayMulti_PointCloud;
    n_rf::Multi_Surfaces *displayMulti_Surfaces;
    //n_rf::Multi_Surfaces *displayMulti_Surfaces2;
    n_rf::Multi_Curves *displayMulti_Curves;

    n_rf::Surface *displaySurface;
    n_rf::Curve *displayCurve;
    n_rf::PointCloud *displayPointCloud;
    n_rf::PointCloud *displayPointCloud_InitNormal;
    n_rf::PointCloud *displayPointCloud_SolveNormal;
    n_rf::Surface *displaySurface_Final;


    n_rf::MultiSurface *displayMulti_Surfaces1;
    n_rf::MultiSurface *displayMulti_Surfaces2;
    n_rf::MultiCellTopo *displayMulti_CellTopo;

    int n_cluster;
    n_rf::Multi_PointCloud *displayMulti_PointCloud_IterInit;


    int GL_state = 0;
    //enum{EMPTY,S_PC,S_Curve,S_MULTI,S_MULTI_M2,S_MULTI_M3,S_MULTI_CS,S_SURF,S_PCurve,S_MultiSurfCut,S_MULTI_USMC};
    enum{EMPTY,S_Curve,S_MULTI_M1,S_MULTI_M2,S_MULTI_CellTopo};
    int currentData;

    vector<int>TopoConstraintInd;
    vector<int>TopoCurInd;
    vector<int>Cell2NTopo;
    int curCellInd;
    vector<int>CellSeq;

    bool isRunDP = false;
    bool GenMode = false;
    bool MultiSurfacecutPlaneMode = false;
    bool specialViewCrossSectionMode  = false;
    bool specialViewCellMode  = false;

    string s_outdir;
    string s_modelname;

private slots:

    void StartTimers_forVideo();
    void StopTimers_forVideo();
    void Timeout_forVideo();

    void StartTimers_forControl();
    void StopTimers_forControl();
    void Timeout_forControl();


    void SaveImage();
    void SaveImage(string fname);
    void renderIntoPixmap(char *pfilename  =NULL);
    void renderIntoPixmap(string filename, int quality = -1);

private slots:

    void readfile();
    void savefile();

    void MainWindowUpdateGL_interface();

    void MainWindowUpdateGL_Multi();
    void MainWindowUpdateGL_Multi2();
    void MainWindowUpdateGL_MMCT();
    void MainWindowUpdateGL_Curve();


    void MainWindowUpdateGL();

    void MainWindowUpdateGL_PartialCurve();

    void MainWindowUpdateGL_MMCT_USMC();

    void MainWindowUpdateGL_PointCloud();

    void MainWindowUpdateGL_Points();
    void MainWindowUpdateGL_Plane();


private slots:
    void SetNMatSlider(int nMat);
    void UpdateCell2Topo(vector<int>&C2NTo);

    void UpdateCellChoices(int pickCell);
    void UpdateTopoChoices(int pickTopo);


    void UpdateVisibilityCellTopos();





    void SetMultiSurfacecutPlaneMode(bool a);

    void RotateCutPlane();

    void UpdataSmoothTransition();
    void ActivatePlayVideo();
    void LoadViewSequence();

    void resizeBallsize();
    void resizeTubesize();
    void resizeDisksize();


private slots:
    void LaplacianSmooth();
    void JuFair();
    void LiepaRefinement();


    void ChangeSmoothingViews();
private slots:
    void runMMmain();
    void InitMMRender(string filename,bool defaultpath);

    void ReadMultiTopoMeta();
    void ReadSuf_Main();
    void ReadContour_Main();

private:


    void ReadSuf1(string filename);
    void ReadSuf2(string filename);

    void ReadCellTopo(string foldername);


private slots:
    void outlinkfunction();



private slots:



    void UpdateViewClusterInits();
    void UpdateVisibilityClusterInits();


private:


    int mode_indicator;
    int n_models;


private slots:



protected:
    //void resizeEvent();
    void resizeEvent(QResizeEvent* event);
    void dropEvent(QDropEvent* event);
    void dragEnterEvent(QDragEnterEvent *ev){ev->accept();}
};

#endif // MAINWINDOW_H
