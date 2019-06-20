#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "ui_dialog.h"
#include "ui_TetManipulation.h"
#include <iostream>
#include <fstream>

using namespace std;

extern double controlBallsize;
extern double controlTubesize;
extern double controlDisksize;



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);


    setAcceptDrops(true);



    GL_state = EMPTY;
    filediag = new QFileDialog;

    glinterface = new GLWidget;

    /******************************************************/

    displayCurve = new n_rf::Curve;
    displaySurface = new n_rf::Surface;
    displayPointCloud = new n_rf::PointCloud;
    displayPointCloud_InitNormal = new n_rf::PointCloud;
    displayPointCloud_SolveNormal = new n_rf::PointCloud;
    displaySurface_Final = new n_rf::Surface;


    displayMulti_PointCloud = new n_rf::Multi_PointCloud;
    displayMulti_Surfaces = new n_rf::Multi_Surfaces;

    displayMulti_Curves = new n_rf::Multi_Curves;

    displayMulti_Surfaces1 = new n_rf::MultiSurface;
    displayMulti_Surfaces2 = new n_rf::MultiSurface;
    displayMulti_CellTopo = new n_rf::MultiCellTopo;

    glinterface->PreSetNum_R_Pointers(6);
    glinterface->PreSetNum_BR_Pointers(3);
    n_cluster = 0;

    /******************************************************/


    hiddenSurface_forCsView = new n_rf::Surface;
    hiddenSurface_forCellsView = new n_rf::Surface;
    hidden_crossSection = new CrossSections();

    /******************************************************/

    connect(&timers_forvideo, SIGNAL(timeout()), this, SLOT(Timeout_forVideo()));
    connect(ui->StartRecordButton,SIGNAL(clicked()),this,SLOT(StartTimers_forVideo()));
    connect(ui->StopRecordButton,SIGNAL(clicked()),this,SLOT(StopTimers_forVideo()));

    connect(ui->PlayVideoButton,SIGNAL(clicked()),this,SLOT(ActivatePlayVideo()));
    connect(ui->LoadViewSequenceButton,SIGNAL(clicked()),this,SLOT(LoadViewSequence()));

    connect(ui->TimerControlStartButton,SIGNAL(clicked()),this,SLOT(StartTimers_forControl()));
    connect(&timers_forcontrol, SIGNAL(timeout()), this, SLOT(Timeout_forControl()));

    /******************************************************/



    /******************************************************/

    /******************************************************/



    /******************************************************/


    connect(ui->MMCrossSectionBox,SIGNAL(toggled(bool)),this,SLOT(MainWindowUpdateGL_interface()));
    connect(ui->MMNonMCurveBox,SIGNAL(toggled(bool)),this,SLOT(MainWindowUpdateGL_interface()));
    connect(ui->MMShowAllCellBox,SIGNAL(toggled(bool)),this,SLOT(MainWindowUpdateGL_interface()));

    connect(ui->FlipMaterialcheckBox,SIGNAL(toggled(bool)),this,SLOT(MainWindowUpdateGL_interface()));

    /******************************************************/

    connect(ui->MMPickSurfSlider,SIGNAL(valueChanged(int)),ui->MMPIckSurfspinBox,SLOT(setValue(int)));
    connect(ui->MMPIckSurfspinBox,SIGNAL(valueChanged(int)),ui->MMPickSurfSlider,SLOT(setValue(int)));

    connect(ui->MMPickSurfSlider,SIGNAL(valueChanged(int)),this,SLOT(MainWindowUpdateGL_interface()));
    /******************************************************/

    connect(ui->MMCellcomboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(UpdateCellChoices(int)));
    connect(ui->MMTopocomboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(UpdateTopoChoices(int)));

    connect(ui->ViewModeComboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(ChangeSmoothingViews()));

    /******************************************************/

    ui->centralWidget->setUpdatesEnabled(true);
    setCentralWidget(ui->centralWidget);
    glinterface->setScreenParameters(ui->mainScreen->width(),ui->mainScreen->height());
    ui->mainScreen->setWidget(glinterface);
    ui->mainScreen->setWidgetResizable(true);

    /******************************************************/


    connect(ui->RenderAmbientSlider,SIGNAL(valueChanged(int)),glinterface,SLOT(SetAmbient(int)));
    connect(ui->RenderDiffuseSlider,SIGNAL(valueChanged(int)),glinterface,SLOT(SetDiffuse(int)));
    connect(ui->RenderSpecularSlider,SIGNAL(valueChanged(int)),glinterface,SLOT(SetSpecular(int)));
    connect(ui->RenderDecadeSlider,SIGNAL(valueChanged(int)),glinterface,SLOT(SetDecade(int)));



    /******************************************************/
    connect(ui->SaveImgButton,SIGNAL(clicked()),this,SLOT(SaveImage()));
    connect(ui->ExportViewButton,SIGNAL(clicked()),glinterface,SLOT(exportViewData()));
    connect(ui->ImportViewButton,SIGNAL(clicked()),glinterface,SLOT(importViewData()));



    connect(ui->RenderSmoothTransitionButton,SIGNAL(clicked()),this,SLOT(UpdataSmoothTransition()));

    connect(ui->SmoothTransitionSlider,SIGNAL(valueChanged(int)),this,SLOT(UpdataSmoothTransition()));

    /******************************************************/



    glinterface->setShowPlane(false);glinterface->setActivationHideHalf(false);glinterface->setInverseHide(false);

    /******************************************************/


    connect(ui->RenderLoadButton_Topo,SIGNAL(clicked()),this,SLOT(ReadMultiTopoMeta()));
    connect(ui->RenderLoadButton_Suf,SIGNAL(clicked()),this,SLOT(ReadSuf_Main()));
    connect(ui->RenderLoadButton_Contour,SIGNAL(clicked()),this,SLOT(ReadContour_Main()));

    /******************************************************/



    connect(ui->RenderSurfaceBox,SIGNAL(toggled(bool)),this,SLOT(MainWindowUpdateGL_interface()));
    connect(ui->RenderNormalBox,SIGNAL(toggled(bool)),this,SLOT(MainWindowUpdateGL_interface()));
    connect(ui->RenderWireBox,SIGNAL(toggled(bool)),this,SLOT(MainWindowUpdateGL_interface()));
    connect(ui->RenderDoubleSide,SIGNAL(toggled(bool)),this,SLOT(MainWindowUpdateGL_interface()));

    connect(ui->RenderNLengthSlider,SIGNAL(valueChanged(int)),this,SLOT(MainWindowUpdateGL_interface()));
    connect(ui->RenderUpNormalSlider,SIGNAL(valueChanged(int)),this,SLOT(MainWindowUpdateGL_interface()));

    //connect(ui->FlipMaterialcheckBox,SIGNAL(toggled(bool)),this,SLOT(MainWindowUpdateGL_Multi()));

    connect(ui->RenderMarkerBox,SIGNAL(toggled(bool)),this,SLOT(MainWindowUpdateGL_interface()));



    connect(ui->RenderLWidthSlider,SIGNAL(valueChanged(int)),glinterface,SLOT(changeLineWidth(int)));


    connect(ui->ballsizeSlider,SIGNAL(valueChanged(int)),this,SLOT(resizeBallsize()));
    connect(ui->tubesizeSlider,SIGNAL(valueChanged(int)),this,SLOT(resizeTubesize()));
    //connect(ui->disksizeSlider,SIGNAL(valueChanged(int)),this,SLOT(resizeBallsize()));

    controlBallsize = ui->ballsizeSlider->value()/50.0;
    controlDisksize = ui->disksizeSlider->value()/50.0;
    controlTubesize = ui->tubesizeSlider->value()/50.0;

    /******************************************************/


    /******************************************************/

    /******************************************************/

    /******************************************************/
    resize(sizeHint());
    glinterface->updateGeometry();
    glinterface->SetAmbient(ui->RenderAmbientSlider->value());
    glinterface->SetDiffuse(ui->RenderDiffuseSlider->value());
    glinterface->SetSpecular(ui->RenderSpecularSlider->value());
    glinterface->SetDecade(ui->RenderDecadeSlider->value());
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::readfile(){
    QString filename =  filediag->getOpenFileName();
    cout<<"test readfile!"<<endl;

}

void MainWindow::InitMMRender(string filename,bool defaultpath){

    cout<<"Define your drop operation."<<endl;
}

void MainWindow::savefile(){
    QString filename =  filediag->getSaveFileName();
    cout<<"test savefile!"<<endl;
    if (filename.length()==0)return;

}
void MainWindow::SaveImage(){
    QString filename =  filediag->getSaveFileName();
    if (filename.length()==0)return;
    renderIntoPixmap(filename.toStdString(),4);

}

void MainWindow::SaveImage(string fname){
    //QString filename =  filediag->getOpenFileName();
    QString filename(fname.data());
    renderIntoPixmap(filename.toStdString(),4);

}

/********************************************************/
void MainWindow::LaplacianSmooth(){

    //displaySurface_multi->SmoothSurfNet_Laplacian(ui->LaplacianspinBox->value());
    //MainWindowUpdateGL_Multi();
}

void MainWindow::JuFair(){
    //displaySurface_multi->SmoothSurfNet_JuFair(ui->JuFairspinBox->value());
    //MainWindowUpdateGL_Multi();
}

void MainWindow::LiepaRefinement(){
   // displaySurface_multi->LiepaRefinement(ui->LiepaSpinBox->value());
    //MainWindowUpdateGL_Multi();
}

void MainWindow::renderIntoPixmap(string filename,int quality){
    QSize size = ui->mainScreen->size();
    //cout<<size.width()<<' '<<size.height()<<endl;
    if (size.isValid()) {
        //QPixmap pixmap = glinterface->renderPixmap(size.width(), size.height(),true);
        QImage pixmap = glinterface->grabFrameBuffer(true);
        //cout<<pixmap.width()<<' '<<pixmap.height()<<endl;
        pixmap.save(QString( filename.c_str()),0,quality);
    }

}
void MainWindow::renderIntoPixmap(char* pfilename)
{
    QSize size = ui->mainScreen->size();
    if (size.isValid()) {
        //QPixmap pixmap = glinterface->renderPixmap(size.width(), size.height(),true);
        QImage pixmap = glinterface->grabFrameBuffer(true);
        cout<<pixmap.width()<<' '<<pixmap.height()<<endl;
        if(pfilename==NULL)pixmap.save("/Users/Research/Geometry/smi_image/a.png",0,4);
        else pixmap.save(pfilename,0,4);
    }

    //    QMediaObject aa(glinterface);
    //    auto recorder = new QMediaRecorder(aa);

    //    QAudioEncoderSettings audioSettings;
    //    audioSettings.setCodec("audio/amr");
    //    audioSettings.setQuality(QMultimedia::HighQuality);

    //    recorder->setAudioSettings(audioSettings);

    //    recorder->setOutputLocation(QUrl::fromLocalFile(fileName));
    //    recorder->record();
}


void MainWindow::resizeEvent(QResizeEvent* event){


    QSize n = event->size();
    QSize o = event->oldSize();

    int w =o.width(),h =o.height();

    n = ui->mainScreen->size();

    cout<<"gl:  "<<n.width()<<' '<<n.height()<<endl;
    glinterface->changeGLwindowsize(n.width(),n.height());

    cout <<"wh: "<<w<<' '<<h<<endl;
}

void MainWindow::dropEvent(QDropEvent* event){

    QList<QUrl> urls = event->mimeData()->urls();

    cout<<"Drop: "<<urls[0].toLocalFile().toStdString()<<endl;

    InitMMRender(urls[0].toLocalFile().toStdString(),false);


}
void MainWindow::MainWindowUpdateGL(){


}


void MainWindow::MainWindowUpdateGL_Multi(){



    cout<<"MainWindowUpdateGL_Multi"<<endl;

    //visibility[ui->leftcomboBox->currentIndex()] = COMPLETELY_VISIBLE;

    displayMulti_Surfaces1->BuildDisplay(infoSurfDisp(true,ui->RenderNormalBox->isChecked(),ui->RenderSurfaceBox->isChecked(),ui->RenderWireBox->isChecked(),
                                                     ui->RenderMarkerBox->isChecked(),ui->FlipMaterialcheckBox->isChecked(),
                                                     ui->RenderNLengthSlider->value(),ui->RenderLWidthSlider->value(),ui->RenderUpNormalSlider->value()));

    vector<bool>b_visibility = displayMulti_Surfaces1->InterpretVisibility(ui->MMCrossSectionBox->isChecked(),ui->MMPickSurfSlider->value()==-1,ui->MMNonMCurveBox->isChecked(),ui->MMPickSurfSlider->value());
    vector<RenderMethod>visibility(b_visibility.size());
    for(int i=0;i<visibility.size();++i)visibility[i] = b_visibility[i]?COMPLETELY_VISIBLE:INVISIBLE;



    glinterface->SetN_BR_Pointers(0, Batch_RenderPointers( displayMulti_Surfaces1->getDisplayVertices(),displayMulti_Surfaces1->getDisplayEdges(),
                                                           displayMulti_Surfaces1->getDisplayFaces(),displayMulti_Surfaces1->getDisplayVerticesNormal(),
                                                           displayMulti_Surfaces1->getDisplayColor(), visibility,
                                                           1 ? COMPLETELY_VISIBLE : INVISIBLE) );



    GL_state = S_MULTI_M1;
    glinterface->DisableExclude(1,0);
    glinterface->update();

}


void MainWindow::MainWindowUpdateGL_Multi2(){


    cout<<"MainWindowUpdateGL_Multi2"<<endl;

    displayMulti_Surfaces2->BuildDisplay(infoSurfDisp(true,ui->RenderNormalBox->isChecked(),ui->RenderSurfaceBox->isChecked(),ui->RenderWireBox->isChecked(),
                                                     ui->RenderMarkerBox->isChecked(),ui->FlipMaterialcheckBox->isChecked(),
                                                     ui->RenderNLengthSlider->value(),ui->RenderLWidthSlider->value(),ui->RenderUpNormalSlider->value()));

    vector<bool>b_visibility = displayMulti_Surfaces2->InterpretVisibility(ui->MMCrossSectionBox->isChecked(),ui->MMPickSurfSlider->value()==-1,ui->MMNonMCurveBox->isChecked(),ui->MMPickSurfSlider->value());
    vector<RenderMethod>visibility(b_visibility.size());
    for(int i=0;i<visibility.size();++i)visibility[i] = b_visibility[i]?COMPLETELY_VISIBLE:INVISIBLE;

    glinterface->SetN_BR_Pointers(1, Batch_RenderPointers( displayMulti_Surfaces2->getDisplayVertices(),displayMulti_Surfaces2->getDisplayEdges(),
                                                           displayMulti_Surfaces2->getDisplayFaces(),displayMulti_Surfaces2->getDisplayVerticesNormal(),
                                                           displayMulti_Surfaces2->getDisplayColor(), visibility,
                                                           1 ? COMPLETELY_VISIBLE : INVISIBLE) );


    GL_state = S_MULTI_M2;
    glinterface->DisableExclude(1,1);
    glinterface->update();



}

void MainWindow::MainWindowUpdateGL_Curve(){

    cout<<"MainWindowUpdateGL_Curve"<<endl;
    displayCurve->BuildDisplay(true,true);

    glinterface->SetN_R_Pointers(0, RenderPointers( displayCurve->getDisplayVertices(),displayCurve->getDisplayEdges(),
                                                    displayCurve->getDisplayFaces(),displayCurve->getDisplayVerticesNormal(),
                                                    displayCurve->getDisplayColor(), COMPLETELY_VISIBLE) );

    GL_state = S_Curve;
    glinterface->updateGL();


}

void MainWindow::MainWindowUpdateGL_PointCloud(){




}

void MainWindow::MainWindowUpdateGL_PartialCurve(){


}

void MainWindow::MainWindowUpdateGL_interface(){

//    if(!GL_state)return;
//    if(GL_state == S_MULTI_M2)MainWindowUpdateGL_MMCT();
//    else if(GL_state == S_MULTI_USMC)MainWindowUpdateGL_MMCT_USMC();
//    else if(GL_state == S_Curve)MainWindowUpdateGL_Curve();
//    else if(GL_state == S_MULTI)MainWindowUpdateGL_Multi();
//    else if(GL_state == S_MULTI_CS)MainWindowUpdateGL_MultiCS();
//    else if(GL_state == S_SURF)MainWindowUpdateGL();
//    else if(GL_state == S_PCurve)MainWindowUpdateGL_PartialCurve();
//    else if(GL_state == S_PC)MainWindowUpdateGL_PointCloud();
//    if(GL_state == S_MULTI_M3)MainWindowUpdateGL_MMCTSeq();

    if(!GL_state)return;
    if(GL_state == S_Curve)MainWindowUpdateGL_Curve();
    else if(GL_state == S_MULTI_M1)MainWindowUpdateGL_Multi();
    else if(GL_state == S_MULTI_M2)MainWindowUpdateGL_Multi2();
    else if(GL_state == S_MULTI_CellTopo)MainWindowUpdateGL_MMCT();



}


void MainWindow::MainWindowUpdateGL_MMCT(){

    if(TopoCurInd.size()<=0)return;
    if(ui->MMPickSurfSlider->value()<-1)return;

    cout<<"MainWindowUpdateGL_MMCT"<<endl;
    displayMulti_CellTopo->BuildDisplay(infoSurfDisp(true,ui->RenderNormalBox->isChecked(),ui->RenderSurfaceBox->isChecked(),ui->RenderWireBox->isChecked(),
                                           ui->RenderMarkerBox->isChecked(),ui->FlipMaterialcheckBox->isChecked(),
                                           ui->RenderNLengthSlider->value(),ui->RenderLWidthSlider->value(),ui->RenderUpNormalSlider->value()),
                              TopoCurInd,ui->MMPickSurfSlider->value());


    UpdateVisibilityCellTopos();
    GL_state = S_MULTI_CellTopo;
    glinterface->DisableExclude(1,2);
    glinterface->update();


}

void MainWindow::MainWindowUpdateGL_MMCT_USMC(){

    cout<<"S_MULTI_USMC"<<endl;

}




void MainWindow::MainWindowUpdateGL_Points(){


}

void MainWindow::MainWindowUpdateGL_Plane(){


}




void MainWindow::RotateCutPlane(){

    glinterface->RotatePlaneOrthogonalToView(1);
}

void MainWindow::UpdateCell2Topo(vector<int>&C2NTo){
    int NCell = C2NTo.size();
    disconnect(ui->MMCellcomboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(UpdateCellChoices(int)));

    Cell2NTopo = C2NTo;
    ui->MMCellcomboBox->setMaxCount(NCell);
    for(int i=0;i<NCell;++i){
        ui->MMCellcomboBox->insertItem(i,QString("Cell_")+QString::number(i));
    }
    ui->MMCellcomboBox->setCurrentIndex(0);
    connect(ui->MMCellcomboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(UpdateCellChoices(int)));

    TopoCurInd.clear();
    TopoCurInd.resize(NCell,0);
    cout<<"NCell: "<<NCell<<endl;

    UpdateCellChoices(0);



}

void MainWindow::UpdateCellChoices(int pickCell){
    if(pickCell<0 || pickCell >= Cell2NTopo.size())return;

    ui->MMTopocomboBox->clear();
    ui->MMTopocomboBox->setEnabled(false);
    curCellInd = pickCell;
    for(int i=0;i<Cell2NTopo[pickCell];++i){
        ui->MMTopocomboBox->insertItem(i,QString("Topo_")+QString::number(i));
    }


    ui->MMTopocomboBox->setEnabled(true);
    if(TopoCurInd.size()<=pickCell){cout<<"Error: UpdateCellChoices: "<< TopoCurInd.size() <<" "<<pickCell<<endl;return;};

    //if(TopoCurInd[pickCell]<0 || TopoCurInd[pickCell]>Cell2NTopo[pickCell]-1 )TopoCurInd[pickCell] = 0;
    ui->MMTopocomboBox->setCurrentIndex(TopoCurInd[pickCell]);

    UpdateTopoChoices(ui->MMTopocomboBox->currentIndex());

    glinterface->SetVisualCell(pickCell);


}

void MainWindow::UpdateTopoChoices(int pickTopo){
    if(TopoCurInd.size()==0)return;
    if(Cell2NTopo.size()==0)return;
    if(!ui->MMTopocomboBox->isEnabled())return;
    if(pickTopo<0 || pickTopo >= Cell2NTopo[ui->MMCellcomboBox->currentIndex()])return;


    int celli = ui->MMCellcomboBox->currentIndex();
    TopoCurInd[celli] = pickTopo;

    MainWindowUpdateGL_MMCT();

}

void MainWindow::SetNMatSlider(int nMat){

    cout<<"nMat: "<<nMat<<endl;
    ui->MMPickSurfSlider->setMaximum(nMat-1);
    ui->MMPIckSurfspinBox->setMaximum(nMat-1);
    ui->MMPickSurfSlider->setMinimum(-1);
    ui->MMPIckSurfspinBox->setMinimum(-1);

}



void MainWindow::SetMultiSurfacecutPlaneMode(bool a){

    MultiSurfacecutPlaneMode = a;
    if(MultiSurfacecutPlaneMode)glinterface->ActivateCutPlane();
    else glinterface->DeActivateCutPlane();
    MainWindowUpdateGL_Multi();

}
void MainWindow::UpdataSmoothTransition(){

    glinterface->ToggleSmoothTransition(ui->SmoothTransitionSlider->value());
}

void MainWindow::ActivatePlayVideo(){


    glinterface->ActivatePlayVideoMode(33);

}
void MainWindow::LoadViewSequence(){

    string path("../view_seq");
    QDir directory(path.data());
    QStringList filters;
    filters << "view_*.dat";
    QStringList filenamelist = directory.entryList(filters, QDir::Files|QDir::NoDotAndDotDot);
    for(auto &a:filenamelist)cout<<a.toStdString()<<endl;
    vector<string>sl(filenamelist.size());
    for(int i=0;i<filenamelist.size();++i)sl[i] = path+string("/")+filenamelist[i].toStdString();
    glinterface->LoadViewSeqs(sl);
}


void MainWindow::resizeBallsize(){
    controlBallsize = ui->ballsizeSlider->value()/50.0;
    MainWindowUpdateGL_interface();
}

void MainWindow::resizeDisksize(){

    controlDisksize = ui->disksizeSlider->value()/50.0;
    MainWindowUpdateGL_interface();

}

void MainWindow::resizeTubesize(){
    controlTubesize = ui->tubesizeSlider->value()/50.0;
    MainWindowUpdateGL_interface();
}



/*******************************************************/
/*******************************************************/
void MainWindow::GetArg(int argc,char** argv){
    this->argc = argc;this->argv = argv;
}

/*******************************************************/
/*******************************************************/

void MainWindow::StartTimers_forVideo(){
    //timers_forvideo.setInterval(ui->RecodIntervalspinBox->value());
    timers_forvideo.setSingleShot(false);
    timers_forvideo.setTimerType(Qt::PreciseTimer);

    accTimeout = 0;

    timers_forvideo.start();

}
void MainWindow::StopTimers_forVideo(){


    timers_forvideo.stop();

    string filename("../videos/test/");
    char ind[50];

    for(int i=0;i<videoStack.size();++i){
        sprintf(ind,"%010d",i);
        string aa = filename + string(ind) + string(".png");
        videoStack[i].save( QString( aa.c_str()),0,100);
    }

    videoStack.clear();


}

void MainWindow::Timeout_forVideo(){


    videoStack.resize(videoStack.size()+1);
    QSize size = ui->mainScreen->size();
    if (size.isValid()) {
        videoStack[videoStack.size()-1] = glinterface->grabFrameBuffer(true);
    }



    //    string filename("../videos/test/");
    //    char ind[50];
    //    sprintf(ind,"%010d",accTimeout);

    //    renderIntoPixmap( filename + string(ind) + string(".png"),100);
    accTimeout++;


}


void MainWindow::StartTimers_forControl(){
    timers_forcontrol.setInterval(1500);
    timers_forcontrol.setSingleShot(false);
    timers_forcontrol.setTimerType(Qt::PreciseTimer);

    accTimeoutControl = 0;
    endTimeoutControl = 6;

    timers_forcontrol.start();


}
void MainWindow::StopTimers_forControl(){



}
void MainWindow::Timeout_forControl(){
    if(accTimeoutControl>=endTimeoutControl){
        timers_forcontrol.stop();
        return;
    }


    vector<int>controlseq({1,3,6,2,5,4});
    glinterface->TriggerShowingMat(controlseq[accTimeoutControl]);
    accTimeoutControl++;
}

/*******************************************************/
/*******************************************************/

void MainWindow::outlinkfunction(){


}

/*******************************************************/
/*******************************************************/

void triggerWarningbox(string s){

    QMessageBox msgBox;
    msgBox.setText(s.data());
    msgBox.exec();

}
/*******************************************************/
/*******************************************************/

