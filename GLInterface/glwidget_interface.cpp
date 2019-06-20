
#include "glwidget.h"

#include <QMouseEvent>
#include <QTimer>
#include <QMessageBox>
#include <math.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <GL/glu.h>
#include <GL/glut.h>
#define BUFFER_OFFSET(bytes) ((GLubyte*) NULL + (bytes))



void GLWidget::setScreenParameters(int w, int h){
    s_width = w;
    s_height = h;
    cent_x = w/2;cent_y = h/2;
    mRadius = min(w,h);
    centerSph = QVector3D(cent_x,cent_y,0.0f);
}
void GLWidget::setTextBrowser(QTextBrowser* in_textbrowser){
    textbrowser = in_textbrowser;
    //setText(-1);
}

//void GLWidget::setText(float area){

//    if(!isRead){
//        textbrowser->setPlainText("Please load volumn!");
//        return;
//    }
//    if(!isContour){
//        textbrowser->setPlainText("Please set the threshold and do contouring!");
//        return;
//    }
//    if(isCutMode){
//        textbrowser->setPlainText("Please control the handles and set the cut plane!");
//        return;
//    }
//    if(!isFilled){
//        textbrowser->setPlainText("Please filled holes!");
//        return;
//    }
//    if(area==0){
//        textbrowser->setPlainText("Please sellect a hole by double click");
//        return;
//    }

//    char temp[100];
//    if(area>0){
//        sprintf(temp,"Area of selected hole:\n%f", area);
//    }
//    else{
//        sprintf(temp,"Thickness of the bone:\n%f", -area);
//    }
//    string text(temp);
//    textbrowser->setPlainText(QString::fromStdString(text));
//}

void GLWidget::SetAmbient(int a){
    gl_ambient = float(a)/100.f;
    updateGL();
}
void GLWidget::SetDiffuse(int a){
    gl_diffuse = float(a)/100.f;
    updateGL();

}
void GLWidget::SetSpecular(int a){
    gl_specular = float(a)/100.f;
    updateGL();

}
void GLWidget::SetDecade(int a){
    gl_decade = float(a)/100.f;
    updateGL();

}


void GLWidget::GetRenderingPointers(vector< GLdouble >* vertices, vector<GLuint> *edges,vector<GLuint>* faces,
                          vector< GLdouble >* v_normals, vector<uchar> *v_color, bool isupdate){


    pvertices = vertices;
    pedges = edges;
    pfaces = faces;
    pv_normals = v_normals;
    pv_color = v_color;
    statue = SINGLE;

    if(isupdate)updateGL();

}

void GLWidget::GetPlaneRenderingPointers(vector< GLdouble >* vertices, vector<GLuint> *edges,vector<GLuint>* faces,
                          vector< GLdouble >* v_normals, vector<uchar> *v_color){


    p_planevertices = vertices;
    p_planeedges = edges;
    p_planefaces = faces;
    p_planev_normals = v_normals;
    p_planev_color = v_color;
    //statue = SINGLE;

    updateGL();

}
void GLWidget::GetConstraintRenderingPointers(vector< GLdouble >* vertices, vector<GLuint> *edges,vector<GLuint>* faces,
                          vector< GLdouble >* v_normals, vector<uchar> *v_color){


    p_constraintvertices = vertices;
    p_constraintedges = edges;
    p_constraintfaces = faces;
    p_constraintv_normals = v_normals;
    p_constraintv_color = v_color;

    updateGL();

}
void GLWidget::GetDotsPointers(vector< GLdouble >* points, vector<uchar> *p_color){

    //vector< GLdouble >* ppoints; vector<uchar> *pp_color;
    ppoints = points;

    pp_color = p_color;

    updateGL();

}


void GLWidget::GetRenderingPointers_Multi(vector< vector< GLdouble >* >* vertices, vector<vector<GLuint>*> *edges, vector<vector<GLuint>*>* faces,
                                          vector< vector< GLdouble >* >* v_normals, vector<vector<uchar> *> *v_color){

    pvertices_multi = vertices;
    pedges_multi = edges;
    pfaces_multi = faces;
    pv_normals_multi = v_normals;
    pv_color_multi = v_color;
    statue = MULTI;
    //isCross = isCrossSec;

    //cout<<isCross<<endl;
    updateGL();

}
void GLWidget::GetRenderingPointers_MultiCS(vector< vector< GLdouble >* >* vertices, vector<vector<GLuint>*> *edges, vector<vector<GLuint>*>* faces,
                                          vector< vector< GLdouble >* >* v_normals, vector<vector<uchar> *> *v_color){

    pvertices_multi = vertices;
    pedges_multi = edges;
    pfaces_multi = faces;
    pv_normals_multi = v_normals;
    pv_color_multi = v_color;
    statue = MULTI_CS;
    //isCross = isCrossSec;

    //cout<<isCross<<endl;
    updateGL();

}
void GLWidget::GetRenderingPointers_MultiM2(vector< vector< GLdouble >* >* vertices, vector<vector<GLuint>*> *edges, vector<vector<GLuint>*>* faces,
                                          vector< vector< GLdouble >* >* v_normals, vector<vector<uchar> *> *v_color){

    pvertices_multi = vertices;
    pedges_multi = edges;
    pfaces_multi = faces;
    pv_normals_multi = v_normals;
    pv_color_multi = v_color;
    statue = MULTI_M2;
    //isCross = isCrossSec;

    //cout<<"Run MULTI_M2"<<endl;
    updateGL();

}
void GLWidget::GetRenderingPointers_MultiM3(vector< vector< GLdouble >* >* vertices, vector<vector<GLuint>*> *edges, vector<vector<GLuint>*>* faces,
                                          vector< vector< GLdouble >* >* v_normals, vector<vector<uchar> *> *v_color,
                                            vector<int>&showlist){

    pvertices_multi = vertices;
    pedges_multi = edges;
    pfaces_multi = faces;
    pv_normals_multi = v_normals;
    pv_color_multi = v_color;
    showCellList = showlist;
    statue = MULTI_M3;
    //isCross = isCrossSec;

    //cout<<"Run MULTI_M2"<<endl;
    updateGL();

}

void GLWidget::PreSetNum_R_Pointers(int n){

    R_Pointers.resize(n);
}
void GLWidget::SetN_R_Pointers(int nth, RenderPointers RP){
    if(nth>=R_Pointers.size())return;
    R_Pointers[nth] = RP;

}
void GLWidget::SetN_R_Pointers_RenderMethod(int nth, RenderMethod RM){
    if(nth>=R_Pointers.size())return;
    R_Pointers[nth].Set_RenderMethod(RM);
}

void GLWidget::DisableExclude(int singleOrBatch, int index){

    if(singleOrBatch == 0){
        if(index < 0 || index >= R_Pointers.size())return;
        for(int i=0;i<R_Pointers.size();++i)R_Pointers[i].Render_Method = INVISIBLE;
        R_Pointers[index].Render_Method = COMPLETELY_VISIBLE;

        for(auto &a:BR_Pointers)a.SelfVisibility = INVISIBLE;
    }else if(singleOrBatch == 1){
        if(index < 0 || index >= BR_Pointers.size())return;
        for(int i=0;i<BR_Pointers.size();++i)BR_Pointers[i].SelfVisibility = INVISIBLE;
        BR_Pointers[index].SelfVisibility = COMPLETELY_VISIBLE;

        for(auto &a:R_Pointers)a.Render_Method = INVISIBLE;
    }

    update();

}

void GLWidget::PreSetNum_BR_Pointers(int n){

    BR_Pointers.resize(n);

}

void GLWidget::SetN_BR_Pointers(int nth, Batch_RenderPointers BRP){

    if(nth>=BR_Pointers.size())return;
    BR_Pointers[nth] = BRP;

}

void GLWidget::SetN_BR_Pointers_RenderMethods(int nth, vector<RenderMethod>&RMlist){

    if(nth>=BR_Pointers.size())return;
    BR_Pointers[nth].Set_RenderMethods(RMlist);

}

void GLWidget::SetN_BR_Pointers_SelfVisibility(int nth, RenderMethod RM){
    if(nth>=BR_Pointers.size())return;
    BR_Pointers[nth].Set_SelfVisibility(RM);

}

void GLWidget::GetRayLine(float mouseX,float mouseY){

    GLint aViewport[4];
    GLdouble matMV[16], matProj[16];
    GLdouble wx, wy, wz;  //  temp world x, y, z coords
    //due to unknown reason, I have to handle the coordinate below:
    //mouseX = (mouseX-1)*2;
    //mouseY = (mouseY-1)*2;
    //cout<<"mouse: "<<mouseX<<' '<<mouseY<<endl;
    mouseX*=ratio_width;
    mouseY*=ratio_height;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glScalef(scale,scale,scale);
    glTranslated(m_trans[0],m_trans[1],m_trans[2]);
    glMultMatrixd(m_RotationMatrix);

    glGetDoublev (GL_MODELVIEW_MATRIX, matMV);

    glPopMatrix();

    glGetIntegerv (GL_VIEWPORT, aViewport);
    glGetDoublev (GL_PROJECTION_MATRIX, matProj);
    //  note viewport[3] is height of window in pixels
    mouseY = aViewport[3] - (GLint) mouseY - 1;
    //mouseX = aViewport[2] - (GLint) mouseX - 1;
    //printf ("Coordinates at cursor are (%4d, %4d)\n", mouseX, mouseY);
    //cout<<"vp: "<<aViewport[0]<<' ' <<aViewport[1]<<' '<<aViewport[2] <<' '<<aViewport[3]<<endl;

    gluUnProject ((GLdouble) mouseX, (GLdouble) mouseY, 0,
                  matMV, matProj, aViewport, &wx, &wy, &wz);
    point_near[0] = wx;point_near[1] = wy;point_near[2] = wz;
    gluUnProject ((GLdouble) mouseX, (GLdouble) mouseY, 1,
                  matMV, matProj, aViewport, &wx, &wy, &wz);



    point_far[0] = wx;point_far[1] = wy;point_far[2] = wz;
    vDirection.setX(point_far[0]-point_near[0]);
    vDirection.setY(point_far[1]-point_near[1]);
    vDirection.setZ(point_far[2]-point_near[2]);
    vDirection.normalize();
    for(int i=0;i<3;++i)raydir[i] = point_far[i] - point_near[i];
    //cout<<"p_near: "<<point_near[0] <<' '<<point_near[1]<<' '<<point_near[2] <<endl;
    //cout<<"p_far: "<<point_far[0] <<' '<<point_far[1]<<' '<<point_far[2] <<endl;
}



void GLWidget::ComputeRayIntersectCutplane(){

//    double planePara[4];
//    double point_near[3],point_far[3],raydir[3];
//    double rayIntersectPlane[3];


    double ddd = 0,ddd2 = 0;
    for(int i=0;i<3;++i)ddd+=point_near[i]*planePara[i];
    for(int i=0;i<3;++i)ddd2+=raydir[i]*planePara[i];
    double t = -(planePara[3]+ddd)/ddd2;
    for(int i=0;i<3;++i)rayIntersectPlane[i] = point_near[i]+t*raydir[i];

}




/*////////////////////////////////////////////////*/

void GLWidget::SetVisualSurface(int a){

    pickSurf = a;
    updateGL();



}
void GLWidget::SetCrossSec(bool a){

    isCross = a;
    updateGL();

}

void GLWidget::SetNonMCurve(bool a){

    isNon = a;
    updateGL();

}
void GLWidget::SetVisualCell(int a){
    pickCell = a+1;
    updateGL();

}
void GLWidget::SetIsSingleCell(bool a){
    isShowAllCell = a;
    updateGL();
}


void GLWidget::SetKeyBoardMode(bool a){
    if(isKeyBoardHideMode == false && a){

        showMatLists.clear();
        showMatLists.resize(128,true);

    }

    cout<<"SetKeyBoardMode "<<endl;
    isKeyBoardHideMode = a;
    updateGL();


}
void GLWidget::TriggerShowingMat(int mati){

    cout<<"TriggerShowingMat: "<<mati<<' '<<showMatLists[mati]<<' ';
    if(isKeyBoardHideMode && mati>=0 && mati< 128){

        showMatLists[mati] = !showMatLists[mati] ;
    }
    cout<<showMatLists[mati]<<endl;
    updateGL();


}

/*////////////////////////////////////////////////*/

void GLWidget::ComputeCutPlanePoints(){
    Eigen::Vector3d oriVecE,newVecE,nc;
    for(int j=0;j<2;++j)oriVecE(j) = 0;oriVecE(2) = 1;
    for(int j=0;j<3;++j)newVecE(j) = planePara[j];
    newVecE.normalize();
    nc = oriVecE.cross(newVecE);
    if(nc.norm()<1e-6){nc(0) = 1.0;nc(1) = 0.0;nc(2) = 0.0;}
    nc.normalize();
    double angle = acos(max(-1.,min(1.,oriVecE.dot(newVecE))));
    Eigen::AngleAxisd m_aa;
    m_aa = Eigen::AngleAxisd(angle,nc);


    vector<Eigen::Vector3d>oriP(4);
    oriP[0](0) = -1; oriP[0](1) = -1; oriP[0](2) = 0;
    oriP[1](0) = -1; oriP[1](1) = 1; oriP[1](2) = 0;
    oriP[2](0) = 1; oriP[2](1) = 1; oriP[2](2) = 0;
    oriP[3](0) = 1; oriP[3](1) = -1; oriP[3](2) = 0;

    //planePara[3] = 0.5;
    for(int i=0;i<4;++i)cutplane[i] = m_aa*oriP[i]-planePara[3]*newVecE;

    oriP[0](0) = -1; oriP[0](1) = 0; oriP[0](2) = 0;
    oriP[1](0) = 1; oriP[1](1) = 0; oriP[1](2) = 0;
    oriP[2](0) = 0; oriP[2](1) = -1; oriP[2](2) = 0;
    oriP[3](0) = 0; oriP[3](1) = 1; oriP[3](2) = 0;

    for(int i=0;i<4;++i)crossline[i] = m_aa*oriP[i]-planePara[3]*newVecE;

    for(int i=0;i<4;++i)cout<<planePara[i]<<' ';cout<<endl;

}
void GLWidget::RotatePlane(int xyz,double steps){
    Eigen::Vector3d Raxis,planeNormal,newplaneN;

    if(xyz==0)Raxis = crossline[0] - crossline[1];
    if(xyz==1)Raxis = crossline[2] - crossline[3];
    Raxis.normalize();

    //Raxis = crossline[xy*2+1]-crossline[xy*2];
    //Raxis.normalize();
    //for(int j=0;j<3;++j)Raxis(j) = 0;
    //Raxis(xyz) = 1;

    for(int j=0;j<3;++j)planeNormal(j) = planePara[j];

    Eigen::AngleAxisd m_aa;
    m_aa = Eigen::AngleAxisd(steps,Raxis);
    newplaneN = m_aa*planeNormal;
    newplaneN.normalize();
    for(int j=0;j<3;++j)planePara[j] = newplaneN(j);
    if(xyz==0)planePara[3] = -newplaneN.dot(crossline[0]);
    if(xyz==1)planePara[3] = -newplaneN.dot(crossline[2]);
}
void GLWidget::RotatePlaneOrthogonalToView(int a){
    GetRayLine(s_width/2,s_height/2);

    for(int i=0;i<3;++i)planePara[i] = vDirection[i];

    ComputeCutPlanePoints();

    emit CutPlaneChanges();

    updateGL();






}

void GLWidget::setViewDirection(double *para){



    for(int i=0;i<3;++i)cout<<para[i]<<' ';cout<<endl;
    Eigen::Vector3d zdir,newVecE,nc;
    for(int i=0;i<3;++i)newVecE(i) = para[i];
    for(int i=0;i<3;++i)zdir(i) = 0;
    zdir(2) = 1.;

    newVecE.normalize();
    nc = zdir.cross(newVecE);
    if(nc.norm()<1e-6){nc(0) = 1.0;nc(1) = 0.0;nc(2) = 0.0;}
    nc.normalize();
    double angle = acos(max(-1.,min(1.,zdir.dot(newVecE))));

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glRotated(-angle/3.1415*180, nc(0),nc(1),nc(2));
    //glMultMatrixd(m_RotationMatrix); //accummulate the previous transformations
    glGetDoublev(GL_MODELVIEW_MATRIX, m_RotationMatrix); //update
    glPopMatrix();


    updateGL();


}

void GLWidget::setActivationHideHalf(bool a){

    isHideHalf = a;
    updateGL();

}

void GLWidget::setInverseHide(bool a){

    cout<<"isInverseHide"<<endl;
    isInverseHide = a;
    updateGL();

}
void GLWidget::setShowPlane(bool a){

    isShowPlane = a;
    updateGL();

}
