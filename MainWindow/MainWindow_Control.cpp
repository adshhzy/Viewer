#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "ui_dialog.h"
#include <iostream>
#include <fstream>
#include <time.h>

#include "mymesh/readers.h"

using namespace std;

void MainWindow::runMMmain(){


}


void MainWindow::ReadSuf_Main(){

    QString filename =  filediag->getOpenFileName(0,tr("Please select Object files"),"",tr("Multilabel Surface (*.suf);;Object files (*.obj *.off)"));

    if(filename.isEmpty())return;

    ReadSuf1(filename.toStdString());



}

void MainWindow::ReadContour_Main(){

    QString filename =  filediag->getOpenFileName(0,tr("Please select .contour files"),"",tr("Contour (*.contour)"));
    if(filename.isEmpty())return;

    CrossSections CSs;
    CSs.ReadCrossSections(filename.toStdString());
    vector<double>points;
    vector<uint>edges;
    CSs.Stackup(points,edges);
    displayCurve->ImportCurve(points,edges);

    MainWindowUpdateGL_Curve();

}

void MainWindow::ReadMultiTopoMeta(){

    //s_outdir = string("/Users/Research/Geometry/MM/MultiTopo/ConvertFolder/ringc/");

    QString foldername =  filediag->getExistingDirectory();
    if(foldername.isEmpty())return;
    s_outdir = foldername.toStdString();

    ReadSuf1(s_outdir+string("/cross-section/testM.suf"));
    ReadSuf2(s_outdir+string("/suf/outcombine.suf"));
    ReadCellTopo(s_outdir + string("/suf/"));


    ui->ViewModeComboBox->setEnabled(true);
    ChangeSmoothingViews();

}

void MainWindow::ChangeSmoothingViews(){



    cout<<"ChangeSmoothingViews"<<endl;
    if(ui->ViewModeComboBox->currentIndex()==0)MainWindowUpdateGL_MMCT();
    else if(ui->ViewModeComboBox->currentIndex()==1)MainWindowUpdateGL_Multi2();
    else if(ui->ViewModeComboBox->currentIndex()==2)MainWindowUpdateGL_Multi();
//    else if(ui->ViewModeComboBox->currentIndex()==3){
//        MainWindowUpdateGL_MultiCS();
//        ui->MMHighlightLoopSlider->setEnabled(true);
//        ui->MMHighlightLoopspinBox->setEnabled(true);
//    }

}

void MainWindow::ReadCellTopo(string foldername){


    //vector<int>Cell2NTopo;
    readVecFile(foldername + string("/Cell2NTopo_vec.txt"),Cell2NTopo);
    displayMulti_CellTopo->ReadCellTopo(foldername+string("/suf"),Cell2NTopo,false);

    SetNMatSlider(displayMulti_CellTopo->GetNMat());
    UpdateCell2Topo(Cell2NTopo);
    ui->MMPickSurfSlider->setValue(0);

    for(auto a:Cell2NTopo)cout<<a<<' ';cout<<endl;



}

void MainWindow::ReadSuf1(string filename){


    double lscale,pcenter[3];
    displayMulti_Surfaces1->Initialize(filename,false,false,false,0,false,lscale,pcenter,1,NULL);

    SetNMatSlider(displayMulti_Surfaces1->GetNMat());

    MainWindowUpdateGL_Multi();

}

void MainWindow::ReadSuf2(string filename){

    double lscale,pcenter[3];
    displayMulti_Surfaces2->Initialize(filename,false,false,false,0,false,lscale,pcenter,1,NULL);
    //displayMulti_Surfaces2->GetScaleInfo(lscale,pcenter);
    MainWindowUpdateGL_Multi2();
    cout<<"scale: "<< lscale <<endl;

}



void MainWindow::UpdateViewClusterInits(){

    vector<RenderMethod>visibility(n_cluster,INVISIBLE);
    //visibility[ui->RBFIterNormalInitPick->currentIndex()] = COMPLETELY_VISIBLE;
    glinterface->SetN_BR_Pointers_RenderMethods(0,visibility);
    glinterface->updateGL();

}

void MainWindow::UpdateVisibilityClusterInits(){

    //glinterface->SetN_BR_Pointers_SelfVisibility(0,ui->RBFICInits->isChecked()? COMPLETELY_VISIBLE:INVISIBLE);
    glinterface->updateGL();
}

void MainWindow::UpdateVisibilityCellTopos(){


    vector<bool>b_visibility = displayMulti_CellTopo->InterpretVisibility(ui->MMCrossSectionBox->isChecked(),ui->MMShowAllCellBox->isChecked(),ui->MMNonMCurveBox->isChecked(),ui->MMCellcomboBox->currentIndex());
    vector<RenderMethod>visibility(b_visibility.size());
    for(int i=0;i<visibility.size();++i)visibility[i] = b_visibility[i]?COMPLETELY_VISIBLE:INVISIBLE;

    glinterface->SetN_BR_Pointers(2, Batch_RenderPointers( displayMulti_CellTopo->getDisplayVertices(),displayMulti_CellTopo->getDisplayEdges(),
                                                           displayMulti_CellTopo->getDisplayFaces(),displayMulti_CellTopo->getDisplayVerticesNormal(),
                                                           displayMulti_CellTopo->getDisplayColor(), visibility,
                                                           1 ? COMPLETELY_VISIBLE : INVISIBLE) );

    glinterface->update();

}
