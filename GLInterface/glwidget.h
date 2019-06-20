/****************************************************************************
**
** Copyright (C) 2013 Digia Plc and/or its subsidiary(-ies).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Digia Plc and its Subsidiary(-ies) nor the names
**     of its contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef GLWIDGET_H
#define GLWIDGET_H
#include <QtWidgets>
#include <QGLWidget>
#include <QGLFunctions>
#include <QTimer>

#include<string>
#include <GL/glu.h>
#include<vector>
#include <eigen3/Eigen/Geometry>

using namespace std;

enum RenderMethod{
    INVISIBLE,
    COMPLETELY_VISIBLE
};

class RenderPointers{
public:
    vector< GLdouble >* pvertices;
    vector<GLuint> *pedges;
    vector<GLuint>* pfaces;
    vector< GLdouble >* pv_normals;
    vector<uchar> *pv_color;
    RenderMethod Render_Method;


    RenderPointers(vector< GLdouble >* vertices, vector<GLuint> *edges,vector<GLuint>* faces,
                              vector< GLdouble >* v_normals, vector<uchar> *v_color, RenderMethod RM):\
        pvertices(vertices),pedges(edges),pfaces(faces),pv_normals(v_normals),pv_color(v_color),Render_Method(RM)
    {}

    RenderPointers(){}

    void Set_RenderMethod(RenderMethod RM){
        Render_Method = RM;
    }
};

class Batch_RenderPointers{

public:
    vector< vector< GLdouble >* >* pvertices_multi;
    vector< vector<GLuint>*> *pedges_multi;
    vector< vector<GLuint>*>* pfaces_multi;
    vector< vector< GLdouble >* >* pv_normals_multi;
    vector< vector<uchar>*>* pv_color_multi;

    vector<RenderMethod>Render_Methods;
    RenderMethod SelfVisibility;

    Batch_RenderPointers(vector< vector< GLdouble >* >* vertices, vector<vector<GLuint>*> *edges, vector<vector<GLuint>*>* faces,
                         vector< vector< GLdouble >* >* v_normals, vector<vector<uchar> *> *v_color,
                         vector<RenderMethod>&RMlist, RenderMethod SelfVisibility):\
        pvertices_multi(vertices),pedges_multi(edges),pfaces_multi(faces),pv_normals_multi(v_normals),pv_color_multi(v_color),Render_Methods(RMlist),SelfVisibility(SelfVisibility)
    {}

    Batch_RenderPointers(){}

    void Set_RenderMethods(vector<RenderMethod>&RMlist){
        Render_Methods = RMlist;
    }

    void Set_SelfVisibility(RenderMethod RM){
        SelfVisibility = RM;
    }

};




class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent = 0);
    ~GLWidget();
    void setScreenParameters(int w, int h);
    int xRotation() const { return xRot; }
    int yRotation() const { return yRot; }
    int zRotation() const { return zRot; }
public slots:


    //void setText(float area);
    void exportViewData(string filename);
    void exportViewData();
    void importViewData(string filename);
    void importViewData();

public slots:
    void changeGLwindowsize(int w, int h);
    void setXRotation(int angle);
    void setYRotation(int angle);
    void setZRotation(int angle);
    void changeLineWidth(int a);
    void changePointSize(int a);




public:
    void pickFaceViaColor(int x, int y, vector<GLdouble> *vertices, vector<GLuint> *faces, vector<uchar> *v_color);

    void setTextBrowser(QTextBrowser* in_textbrowser);

signals:
    void xRotationChanged(int angle);
    void yRotationChanged(int angle);
    void zRotationChanged(int angle);

protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *event);
    void mouseDoubleClickEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
public:
    void paintGL_old();

public:
    void getkeyPressEvent(QKeyEvent *event){keyPressEvent(event);}
private slots:


private:

    enum{VERTEX,FACE,V_NORMAL,F_NORMAL,V_COLOR,NUM};
    enum{LEFT,RIGHT,TOP,BOTTOM};
    GLuint vboIds[NUM];

private:
    int xRot;
    int yRot;
    int zRot;

    QPoint lastPos;
    double m_RotationMatrix[16];
    double m_trans[3];

    double load_m_RotationMatrix[16];
    double load_m_trans[3];
    double load_scale;

    vector<vector<double>> vload_m_RotationMatrix;
    vector<vector<double>> vload_m_trans;
    vector<double> vload_scale;


    float mRadius;
    float cent_x,cent_y;
    QVector3D vDirection;
    QVector3D centerSph;

    float linewidth;
    float pointsize;

    float gl_ambient = 0.5;
    float gl_diffuse = 0.3;
    float gl_specular = 0.0;
    float gl_decade = 0.0;
private:
    GLdouble scale;
    GLdouble zdepth;
    int s_width,s_height;
    int r_width,r_height;
    float ratio_width,ratio_height;
private:
    bool isRead;
    bool isOnConstrain;
private:
    QTextBrowser* textbrowser;
    QTimer timer_forplay;

    int accTimeout;
private:
    GLdouble front_surface[4];
    GLdouble zNear,Zfar,zTrans;
    double point_near[3],point_far[3],raydir[3];
    double rayIntersectPlane[3];
private:
    int pickInd;
//private:
//    n_rf::CurveNet c_frame;
//    n_rf::SurfaceField c_surf;
//    n_rf::TetrahedronField c_vol;
private:
    enum{SINGLE,MULTI,MULTI_M2,MULTI_M3,MULTI_CS};
    int statue;
private:
    bool surf_constrainOn;
private:
    void normalizeAngle(int *angle);
    void drawMesh(vector< GLfloat >* vertices, vector<GLuint>* faces,
                vector< GLfloat >* v_normals = NULL, vector<uchar> *v_color = NULL);
    void drawGrid(vector< GLfloat >* vertices, vector<GLuint>* edges,
                vector< GLfloat >* v_normals = NULL, vector<uchar> *v_color = NULL);
    void drawGrid(vector< GLdouble >* vertices, vector<GLuint> *edges,
                vector< GLdouble >* v_normals = NULL, vector<uchar> *v_color = NULL);
    void drawMeshandGrid(vector< GLdouble >* vertices, vector<GLuint> *edges,vector<GLuint>* faces,
                vector< GLdouble >* v_normals = NULL, vector<uchar> *v_color = NULL);
    void drawPoint(vector<double> *points, vector<uchar> *pcolor);
    void GetRayLine(float mouseX, float mouseY);

//    void drawMeshandGrid_FaceMode(vector< GLdouble >* vertices, vector<GLuint> *edges,vector<GLuint>* faces,
//                vector< GLdouble >* v_normals = NULL, vector<uchar> *f_color = NULL);
/************************************************************************/
private:
    vector< GLdouble >* pvertices; vector<GLuint> *pedges;vector<GLuint>* pfaces;
    vector< GLdouble >* pv_normals ; vector<uchar> *pv_color;
public:
    void GetRenderingPointers(vector< GLdouble >* vertices, vector<GLuint> *edges, vector<GLuint>* faces,
                              vector< GLdouble >* v_normals, vector<uchar> *v_color, bool isupdate = true);
/************************************************************************/
private:
    vector< GLdouble >* p_planevertices; vector<GLuint> *p_planeedges;vector<GLuint>* p_planefaces;
    vector< GLdouble >* p_planev_normals ; vector<uchar> *p_planev_color;
public:
    void GetPlaneRenderingPointers(vector< GLdouble >* vertices, vector<GLuint> *edges,vector<GLuint>* faces,
                              vector< GLdouble >* v_normals, vector<uchar> *v_color);

/************************************************************************/
private:
    vector< GLdouble >* p_constraintvertices; vector<GLuint> *p_constraintedges;vector<GLuint>* p_constraintfaces;
    vector< GLdouble >* p_constraintv_normals ; vector<uchar> *p_constraintv_color;
public:
    void GetConstraintRenderingPointers(vector< GLdouble >* vertices, vector<GLuint> *edges,vector<GLuint>* faces,
                              vector< GLdouble >* v_normals, vector<uchar> *v_color);

/************************************************************************/
private:
    vector< GLdouble >* ppoints = NULL; vector<uchar> *pp_color = NULL;
public:
    void GetDotsPointers(vector< GLdouble >* points, vector<uchar> *p_color);
/************************************************************************/
private:
    vector< vector< GLdouble >* >* pvertices_multi; vector<vector<GLuint>*> *pedges_multi;vector<vector<GLuint>*>* pfaces_multi;
    vector< vector< GLdouble >* >* pv_normals_multi ; vector<vector<uchar>*>* pv_color_multi;
public:
    void GetRenderingPointers_Multi(vector< vector< GLdouble >* >* vertices, vector<vector<GLuint>*> *edges,vector<vector<GLuint>*>* faces,
                              vector< vector< GLdouble >* >* v_normals, vector<vector<uchar> *> *v_color);

    void GetRenderingPointers_MultiM2(vector< vector< GLdouble >* >* vertices, vector<vector<GLuint>*> *edges, vector<vector<GLuint>*>* faces,
                                              vector< vector< GLdouble >* >* v_normals, vector<vector<uchar> *> *v_color);
    void GetRenderingPointers_MultiM3(vector< vector< GLdouble >* >* vertices, vector<vector<GLuint>*> *edges, vector<vector<GLuint>*>* faces,
                                              vector< vector< GLdouble >* >* v_normals, vector<vector<uchar> *> *v_color,
                                                vector<int>&showlist);

    void GetRenderingPointers_MultiCS(vector< vector< GLdouble >* >* vertices, vector<vector<GLuint>*> *edges, vector<vector<GLuint>*>* faces,
                                              vector< vector< GLdouble >* >* v_normals, vector<vector<uchar> *> *v_color);



private:

    vector<RenderPointers>R_Pointers;
    vector<Batch_RenderPointers>BR_Pointers;

public:

    void DisableExclude(int singleOrBatch, int index);

    void PreSetNum_R_Pointers(int n);
    void SetN_R_Pointers(int nth, RenderPointers RP);
    void SetN_R_Pointers_RenderMethod(int nth, RenderMethod RM);

    void PreSetNum_BR_Pointers(int n);
    void SetN_BR_Pointers(int nth, Batch_RenderPointers BRP);
    void SetN_BR_Pointers_RenderMethods(int nth, vector<RenderMethod>&RMlist);
    void SetN_BR_Pointers_SelfVisibility(int nth, RenderMethod RM);

public slots:
    void SetAmbient(int a);
    void SetDiffuse(int a);
    void SetSpecular(int a);
    void SetDecade(int a);




private:
    int pickSurf;
    bool isCross;
    bool isNon;
    bool isShowAllCell;
    bool isKeyBoardHideMode = false;
    int pickCell;
    vector<int>showCellList;
    vector<int>showMatLists;
public slots:
    void SetVisualSurface(int a);
    void SetCrossSec(bool a);
    void SetNonMCurve(bool a);

    void SetVisualCell(int a);
    void SetIsSingleCell(bool a);

    void SetKeyBoardMode(bool a);
    void TriggerShowingMat(int mati);




private:
    double planePara[4];
    bool isplaneOn;
    bool isShowPlane;
    bool isUserAddingConstrainMode;
    bool nextEdges = false;
    Eigen::Vector3d cutplane[4];
    Eigen::Vector3d crossline[4];
    bool isHideHalf;
    bool isInverseHide;
private:
    void ComputeCutPlanePoints();
    void RotatePlane(int xyz, double steps);
    void ComputeRayIntersectCutplane();
public:
    void ActivateCutPlane(){isplaneOn = true;}
    void ActivateUserAddingConstrainMode(){isUserAddingConstrainMode = true;nextEdges = false;isplaneOn = true;}
    bool GetIsUserAddingConstrainMode(){return isUserAddingConstrainMode;}

    double *GetCurrentPlanePara(){return planePara;}
    bool GetNextEdges(){return nextEdges;}

    double *GetRayIntersectPlane(){return rayIntersectPlane;}
signals:

    void addNewPoint();
    void UpdateCurserPoints();
    void CutPlaneChanges();


public slots:

    void DeActivateUserAddingConstrainMode(){isUserAddingConstrainMode = false;isplaneOn = false;}
    void DeActivateCutPlane(){isplaneOn = false;}
    void setActivationHideHalf(bool a);
    void setInverseHide(bool a);
    void setShowPlane(bool a);

    void RotatePlaneOrthogonalToView(int a);

    void ToggleSmoothTransition(int interp_factor_in);
    void ToggleSmoothTransition_forVideo(int loadi, int loadj, double interp_factor_in);

    void LoadViewSeqs(vector<string>&filelists);
    void ActivatePlayVideoMode(int interval);
    void ProcessTimeout();

    void setViewDirection(double *para);
};

#endif // GLWIDGET_H
