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
#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#define BUFFER_OFFSET(bytes) ((GLubyte*) NULL + (bytes))


#include <GLFW/glfw3native.h>
#include <GLFW/glfw3.h>

using namespace std;
GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(parent)
{
    QImage a;

    xRot = 0;
    yRot = 0;
    zRot = 0;

    front_surface[LEFT] = -1;
    front_surface[RIGHT] = 1;
    front_surface[TOP] = 1;
    front_surface[BOTTOM] = -1;
    zNear = 1;
    Zfar = 3;
    zTrans = -2;
    point_near[0] = point_near[1] = point_near[2] = 0;
    point_far[0] = point_far[1] = point_far[2] = 0;
    isRead=false;

    //setText(-1);

    connect(&timer_forplay, SIGNAL(timeout()), this, SLOT(ProcessTimeout()));
    memset(m_RotationMatrix, 0.0, sizeof(double)*16);
    for (unsigned int i=0; i< 4; ++i){
        m_RotationMatrix[i+(4*i)] = 1.0;
    }
    m_trans[0] = 0;m_trans[1] = 0;m_trans[2] = 0;

    memset(load_m_RotationMatrix, 0.0, sizeof(double)*16);
    for (unsigned int i=0; i< 4; ++i){
        load_m_RotationMatrix[i+(4*i)] = 1.0;
    }
    load_m_trans[0] = 0;load_m_trans[1] = 0;load_m_trans[2] = 0;


    statue = SINGLE;

    surf_constrainOn = false;

    pvertices = NULL;
    pedges = NULL;
    pfaces = NULL;
    pv_normals = NULL;
    pv_color = NULL;

    pvertices_multi = NULL;
    pedges_multi = NULL;
    pfaces_multi = NULL;
    pv_normals_multi = NULL;
    pv_color_multi = NULL;


    p_planevertices = NULL;
    p_planeedges = NULL;
    p_planefaces = NULL;
    p_planev_normals = NULL;
    p_planev_color = NULL;

    p_constraintvertices = NULL;
    p_constraintedges = NULL;
    p_constraintfaces = NULL;
    p_constraintv_normals = NULL;
    p_constraintv_color = NULL;





    pickSurf = -1;isCross = true;isNon=true;
    isShowAllCell = true; pickCell = 0;

    for(int i=0;i<4;++i)planePara[i]=0;
    planePara[0] = 1;planePara[3]=0.00;//0.139
    ComputeCutPlanePoints();
    isplaneOn = false;
    isUserAddingConstrainMode = false;

    isHideHalf = false;
    isInverseHide = false;

    setFocusPolicy(Qt::StrongFocus);

    setMouseTracking(true);

    //double aaaa[4] = {0.00732758, 0.257745, -0.966185, -0.04 };0.55946 -0.780027 -0.280289 0.02
     //double aaaa[4] = {0.55946, -0.822394, -0.103308, 0.102 };//0.54893 -0.827268 -0.119597 0.0819508
    //double aaaa[4] = {0.551243, -0.0561991, -0.83245, -0.2};
    //0.00387925 0.947625 0.319361 0.288125 -> liver CS
    //double aaaa[4] = {0.00387925, 0.947625, 0.319361, 0.288125};
    //double aaaa[4] = {-0.788719, 0.585408, 0.187669, 0.12};
    double aaaa[4] = {1, 0, 0, 0.0};


    for(int i=0;i<4;++i)planePara[i] = aaaa[i];
    ComputeCutPlanePoints();

}

GLWidget::~GLWidget()
{
    makeCurrent();

}

void GLWidget::setXRotation(int angle)
{
    normalizeAngle(&angle);
    if (angle != xRot) {
        xRot = angle;
        emit xRotationChanged(angle);
        //updateGL();
    }
}

void GLWidget::setYRotation(int angle)
{
    normalizeAngle(&angle);
    if (angle != yRot) {
        yRot = angle;
        emit yRotationChanged(angle);
        //updateGL();
    }
}

void GLWidget::setZRotation(int angle)
{
    normalizeAngle(&angle);
    if (angle != zRot) {
        zRot = angle;
        emit zRotationChanged(angle);
        //updateGL();
    }
}

void GLWidget::initializeGL()
{
    static const GLfloat lightPos2[4] = { 1.0f, 1.0f, 0.0f, 1.0f };
    static const GLfloat lightPos3[4] = { -1.0f, 1.0f, 0.0f, 1.0f };
    static const GLfloat lightPos4[4] = { 0.0f, 0.0f, -1.0f, 1.0f };
    static const GLfloat lightPos[4] = { 0.0f, 0.0f, 1.0f, 1.0f };
    static const GLfloat reflectance1[4] = { 0.8f, 0.1f, 0.0f, 0.9f };
    static const GLfloat reflectance2[4] = { 0.5f, 0.5f, 0.5f, 0.9f };
    static const GLfloat reflectance3[4] = { 1.f, 1.f, 1.f, 1.f };
    static const GLfloat reflectance4[4] = { 0.0f, 0.8f, 0.2f, 0.9f };

    //glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
    glLightfv(GL_LIGHT1, GL_POSITION, lightPos2);
    glLightfv(GL_LIGHT2, GL_POSITION, lightPos3);
    //glLightfv(GL_LIGHT3, GL_POSITION, lightPos4);

    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,reflectance2);
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,reflectance2);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,reflectance2);
    float a = 10;
    //glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &a);
    glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);
    // glEnable(GL_LIGHT1);
    // glEnable(GL_LIGHT2);
    //glEnable(GL_LIGHT3);

    //glEnable(GL_CULL_FACE);
    //glCullFace(GL_FRONT);

    glEnable(GL_DEPTH_TEST);

    glEnable(GL_COLOR_MATERIAL);

    glEnable(GL_NORMALIZE);
    glEnable(GL_BLEND);
    glEnable(GL_ALPHA_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glGenBuffers(NUM, vboIds);
    scale = 0.5;load_scale = 0.5;
    zdepth = 0;
    linewidth = 2.;
    pointsize = 4.;

    //glPointSize(6.);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_POLYGON_SMOOTH, GL_NICEST);

    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH, GL_NICEST);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_MULTISAMPLE);

    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);



    GLfloat light_ambient[] = { 0.1, 0.1, 0.1, 1.0 };
    GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat light_position[] = { 0.0, 1.0, 1.0, 0.0 };

    GLfloat light_off[] = { .0, .0, .0, 1.0 };

    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glEnable(GL_LIGHTING);

    glEnable(GL_DEPTH_TEST);



    GLfloat light1_position[] = {1., -0., 0.0, 1.0 };

    glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
    //glEnable(GL_LIGHT1);

    GLfloat light2_position[] = { 0., -1.0, 0.0, 1.0 };

    glLightfv(GL_LIGHT2, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT2, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT2, GL_POSITION, light2_position);

    GLfloat light3_position[] = { 1.0, 1.0, 1.0, 0.4 };

    glLightfv(GL_LIGHT3, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT3, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT3, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT3, GL_POSITION, light3_position);

    glEnable(GL_LIGHT0);
//    glEnable(GL_LIGHT1);
//    glEnable(GL_LIGHT2);
//    glEnable(GL_LIGHT3);

    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);


}

void GLWidget::paintGL_old()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    qglClearColor(QColor(255,255,255,255));

    vector<GLfloat> light_ambient(4,gl_ambient);
    vector<GLfloat> light_diffuse(4,gl_diffuse);
    vector<GLfloat> light_specular(4,gl_specular);
    light_ambient[3] = gl_decade;light_diffuse[3] = gl_decade;light_specular[3] = gl_decade;

    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient.data());
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse.data());
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular.data());

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    glScalef(scale,scale,scale);
    glTranslated(m_trans[0],m_trans[1],m_trans[2]);
    glLineWidth(linewidth);

    glMultMatrixd(m_RotationMatrix);




    if((isplaneOn || isUserAddingConstrainMode )&& isHideHalf){
        double ipp[4];
        if(!isInverseHide)for(int i=0;i<4;++i)ipp[i] = planePara[i];
        else{
            for(int i=0;i<4;++i)ipp[i] = -planePara[i];
        }
        if(isplaneOn)ipp[3]+=0.00001;
        glClipPlane(GL_CLIP_PLANE0,ipp);
        glEnable(GL_CLIP_PLANE0);
    }

    if(statue == SINGLE)drawMeshandGrid(pvertices,pedges,pfaces,pv_normals,pv_color);
    else if(statue == MULTI){
        if(pvertices_multi!=NULL){
            //cout<<"MULTI1"<<endl;

            int sss = pvertices_multi->size()-1;
            if(isCross)drawMeshandGrid(pvertices_multi->at(0),pedges_multi->at(0),pfaces_multi->at(0),pv_normals_multi->at(0),pv_color_multi->at(0));
            for(int i =1;i<sss-1;++i){
                //cout<<"multi rendering!"<<endl;
                if(!isKeyBoardHideMode){if(i!=pickSurf+1)continue;}
                else {if(!showMatLists[i-1] || i==1)continue;}
                //if(pickSurf>=0 && i!=pickSurf)continue;
                //if(i==1 || i == sss-2 || i==sss-1)continue;
                //cout<<sss<<endl;
                drawMeshandGrid(pvertices_multi->at(i),pedges_multi->at(i),pfaces_multi->at(i),pv_normals_multi->at(i),pv_color_multi->at(i));
            }
            if(!isKeyBoardHideMode)if(pickSurf==-1)drawMeshandGrid(pvertices_multi->at(sss),pedges_multi->at(sss),pfaces_multi->at(sss),pv_normals_multi->at(sss),pv_color_multi->at(sss));
            if(isNon)drawMeshandGrid(pvertices_multi->at(sss-1),pedges_multi->at(sss-1),pfaces_multi->at(sss-1),pv_normals_multi->at(sss-1),pv_color_multi->at(sss-1));
        }

    }else if(statue == MULTI_CS){
        if(pvertices_multi!=NULL){
            //cout<<"MULTI1"<<endl;

            int sss = pvertices_multi->size()-1;
            if(isCross)drawMeshandGrid(pvertices_multi->at(0),pedges_multi->at(0),pfaces_multi->at(0),pv_normals_multi->at(0),pv_color_multi->at(0));
            for(int i =1;i<sss;++i){
                //cout<<"multi rendering!"<<endl;
                if(pickSurf>=0 && i!=pickSurf+1)continue;
                //if(pickSurf>=0 && i!=pickSurf)continue;
                if(i==1 || i == sss-2 || i==sss-1)continue;
                //cout<<sss<<endl;
                drawMeshandGrid(pvertices_multi->at(i),pedges_multi->at(i),pfaces_multi->at(i),pv_normals_multi->at(i),pv_color_multi->at(i));
            }
            if(isNon)drawMeshandGrid(pvertices_multi->at(sss-1),pedges_multi->at(sss-1),pfaces_multi->at(sss-1),pv_normals_multi->at(sss-1),pv_color_multi->at(sss-1));
        }

        drawMeshandGrid(pvertices,pedges,pfaces,pv_normals,pv_color);

    }else if(statue == MULTI_M2){
        if(pvertices_multi!=NULL){
            int sss = 1+(pvertices_multi->size()-1)/2;
            if(isCross)drawMeshandGrid(pvertices_multi->at(0),pedges_multi->at(0),pfaces_multi->at(0),pv_normals_multi->at(0),pv_color_multi->at(0));
            if(isShowAllCell){

                for(int i =1;i<sss;++i){
                    //if(i!=1)continue;
                    drawMeshandGrid(pvertices_multi->at(i),pedges_multi->at(i),pfaces_multi->at(i),pv_normals_multi->at(i),pv_color_multi->at(i));
                }
                if(isNon)for(int k =1;k<sss;++k){
                    //if(k!=pickCell)continue;
                    int i=k+(pvertices_multi->size()-1)/2;
                    drawMeshandGrid(pvertices_multi->at(i),pedges_multi->at(i),pfaces_multi->at(i),pv_normals_multi->at(i),pv_color_multi->at(i));
                }
            }
            else{

                for(int i =1;i<sss;++i){
                    if(pickCell>-1 && i!=pickCell)continue;
                    drawMeshandGrid(pvertices_multi->at(i),pedges_multi->at(i),pfaces_multi->at(i),pv_normals_multi->at(i),pv_color_multi->at(i));
                }
                if(isNon)for(int k =1;k<sss;++k){
                    if(pickCell>-1 &&k!=pickCell)continue;
                    int i=k+(pvertices_multi->size()-1)/2;
                    drawMeshandGrid(pvertices_multi->at(i),pedges_multi->at(i),pfaces_multi->at(i),pv_normals_multi->at(i),pv_color_multi->at(i));
                }
            }
        }
    }else if(statue == MULTI_M3){
        if(pvertices_multi!=NULL){
            if(isCross)drawMeshandGrid(pvertices_multi->at(0),pedges_multi->at(0),pfaces_multi->at(0),pv_normals_multi->at(0),pv_color_multi->at(0));
            for(int k =0;k<showCellList.size();++k){
                int i=showCellList[k]+1;
                drawMeshandGrid(pvertices_multi->at(i),pedges_multi->at(i),pfaces_multi->at(i),pv_normals_multi->at(i),pv_color_multi->at(i));
            }
        }
    }

    //drawPoint(pvertices);
    if(isplaneOn && !isUserAddingConstrainMode){

        if(isShowPlane){
            glBegin (GL_LINES);
            glColor4f(0.0, 0.0, 0.0, 0.6);
            for(int i = 0;i<4;i++)
                glVertex3d(crossline[i](0),crossline[i](1),crossline[i](2));
            glEnd ();

            glBegin (GL_POLYGON);
            glColor4f(0.8, 0.8, 0.8, 0.3);
            for(int i = 0;i<4;i++)
                glVertex3d(cutplane[i](0),cutplane[i](1),cutplane[i](2));
           glEnd ();

            drawPoint(ppoints,pp_color);
        }

        //drawMeshandGrid(pvertices,pedges,pfaces,pv_normals,pv_color);




    }

    if(isUserAddingConstrainMode){

        glDisable(GL_CLIP_PLANE0);
        if(p_planevertices!=NULL && p_planefaces!=NULL){
            drawMeshandGrid(p_planevertices,p_planeedges,p_planefaces,p_planev_normals,p_planev_color);
        }


        if(p_constraintvertices!=NULL && p_constraintfaces!=NULL){
            //cout<<p_constraintvertices->size()<<endl;
            drawMeshandGrid(p_constraintvertices,p_constraintedges,p_constraintfaces,p_constraintv_normals,p_constraintv_color);
        }
        //drawPoint(ppoints,pp_color);


        //glPointSize(10);

        //glColor4f(0.0, 0.0, 0.0, 0.9);
        //glBegin(GL_POINTS);
        //glVertex3d(rayIntersectPlane[0],rayIntersectPlane[1],rayIntersectPlane[2]);

        glEnd();
        //emit addNewPoint();
    }

    if(isplaneOn || isUserAddingConstrainMode )glDisable(GL_CLIP_PLANE0);


    if(isplaneOn && !isUserAddingConstrainMode)drawMeshandGrid(pvertices,pedges,pfaces,pv_normals,pv_color);


    glPopMatrix();

}


void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    qglClearColor(QColor(255,255,255,255));

    vector<GLfloat> light_ambient(4,gl_ambient);
    vector<GLfloat> light_diffuse(4,gl_diffuse);
    vector<GLfloat> light_specular(4,gl_specular);
    light_ambient[3] = gl_decade;light_diffuse[3] = gl_decade;light_specular[3] = gl_decade;

    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient.data());
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse.data());
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular.data());

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    glScalef(scale,scale,scale);
    glTranslated(m_trans[0],m_trans[1],m_trans[2]);
    glLineWidth(linewidth);

    glMultMatrixd(m_RotationMatrix);

    if(isHideHalf){
        double ipp[4];
        if(!isInverseHide)for(int i=0;i<4;++i)ipp[i] = planePara[i];
        else{
            for(int i=0;i<4;++i)ipp[i] = -planePara[i];
        }
        ipp[3]+=0.00001;
        glClipPlane(GL_CLIP_PLANE0,ipp);
        glEnable(GL_CLIP_PLANE0);
    }



    for(auto &a:R_Pointers){
        if(a.Render_Method==COMPLETELY_VISIBLE)drawMeshandGrid(a.pvertices,a.pedges,a.pfaces,a.pv_normals,a.pv_color);
    }
    for(auto &a:BR_Pointers){
        if(a.SelfVisibility!=COMPLETELY_VISIBLE)continue;
        for(int i =0; i<a.Render_Methods.size(); ++i){
            if(a.Render_Methods[i]==COMPLETELY_VISIBLE)drawMeshandGrid(a.pvertices_multi->at(i),a.pedges_multi->at(i),a.pfaces_multi->at(i),a.pv_normals_multi->at(i),a.pv_color_multi->at(i));
        }
    }

    if(isShowPlane){
        glBegin (GL_LINES);
        glColor4f(0.0, 0.0, 0.0, 0.6);
        for(int i = 0;i<4;i++)
            glVertex3d(crossline[i](0),crossline[i](1),crossline[i](2));
        glEnd ();

        glBegin (GL_POLYGON);
        glColor4f(0.8, 0.8, 0.8, 0.3);
        for(int i = 0;i<4;i++)
            glVertex3d(cutplane[i](0),cutplane[i](1),cutplane[i](2));
        glEnd ();
        drawPoint(ppoints,pp_color);
    }

    if(isHideHalf)glDisable(GL_CLIP_PLANE0);
    glPopMatrix();

}



void GLWidget::resizeGL(int width, int height)
{
    glMatrixMode(GL_VIEWPORT);
    int side = qMin(width, height);
    glViewport(0,0,width,height);

    r_width = width;r_height = height;
    ratio_width = float(r_width)/s_width;ratio_height = float(r_height)/s_height;



    cout << "GLresize: "<<width<<' '<<height<<endl;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //front_surface[LEFT] = -(front_surface[RIGHT] = s_width/s_height);
    //    glFrustum(front_surface[LEFT], front_surface[RIGHT], front_surface[BOTTOM], front_surface[TOP],
    //            zNear, Zfar);
    gluPerspective(30,double(width)/height,zNear,Zfar);
    //gluPerspective(3,s_width/s_height,zNear,Zfar);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(0.0, 0.0, zTrans);
}
void GLWidget::drawPoint(vector<double> *points,vector<uchar>*pcolor){

    if(points==NULL)return;
    //cout<<"Points:  "<<points->size()<<endl;
    glPointSize(10);

    //glColor4f(0.0, 0.0, 0.0, 0.9);
    int nv = points->size()/3;
    double *p_d = points->data();uchar *p_c = pcolor->data();
    glBegin(GL_POINTS);
    //for(int i=0;i<pcolor->size();i++)cout<<int(p_c[i])<<' ';cout<<endl;
    for(int i=0;i<nv;i++){
        auto pdd = p_d+i*3;
        auto pcc = p_c+i*4;
        //glColor4i(133, 133, 133,255);
        glColor4ubv(pcc);
        glVertex3d(pdd[0],pdd[1],pdd[2]);
    }
    //glVertex3d(0,0,0);
    glEnd();
    //glFlush();


}
void GLWidget::drawMesh(vector< GLfloat >* vertices,vector<GLuint>* faces,
                        vector< GLfloat >* v_normals,vector< uchar >* v_color)
{
    glBindBuffer(GL_ARRAY_BUFFER, vboIds[VERTEX]);
    glBufferData(GL_ARRAY_BUFFER, vertices->size()*sizeof(GLfloat), vertices->data(), GL_DYNAMIC_DRAW);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3,GL_FLOAT,0,BUFFER_OFFSET(0));

    if(v_normals!=NULL){
        glBindBuffer(GL_ARRAY_BUFFER, vboIds[V_NORMAL]);
        glBufferData(GL_ARRAY_BUFFER, v_normals->size()*sizeof(GLfloat), v_normals->data(), GL_DYNAMIC_DRAW);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT,0,BUFFER_OFFSET(0));
    }
    if(v_color != NULL){
        glBindBuffer(GL_ARRAY_BUFFER, vboIds[V_COLOR]);
        glBufferData(GL_ARRAY_BUFFER, v_color->size()*sizeof(uchar), v_color->data(), GL_DYNAMIC_DRAW);
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(4,GL_UNSIGNED_BYTE,0,BUFFER_OFFSET(0));

    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vboIds[FACE]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, faces->size()*sizeof(GLuint), faces->data(), GL_DYNAMIC_DRAW);
    // Draw cube geometry using indices from VBO 1
    glDrawElements(GL_TRIANGLES, faces->size(), GL_UNSIGNED_INT, BUFFER_OFFSET(0));
    glDisableClientState(GL_COLOR_ARRAY);


    //glPopMatrix();
}
void GLWidget::drawGrid(vector< GLfloat >* vertices,vector<GLuint>* edges,
                        vector< GLfloat >* v_normals,vector< uchar >* v_color)
{
    glBindBuffer(GL_ARRAY_BUFFER, vboIds[VERTEX]);
    glBufferData(GL_ARRAY_BUFFER, vertices->size()*sizeof(GLfloat), vertices->data(), GL_DYNAMIC_DRAW);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3,GL_FLOAT,0,BUFFER_OFFSET(0));


    if(v_normals!=NULL){
        glBindBuffer(GL_ARRAY_BUFFER, vboIds[V_NORMAL]);
        glBufferData(GL_ARRAY_BUFFER, v_normals->size()*sizeof(GLfloat), v_normals->data(), GL_DYNAMIC_DRAW);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT,0,BUFFER_OFFSET(0));
    }
    if(v_color != NULL){
        glBindBuffer(GL_ARRAY_BUFFER, vboIds[V_COLOR]);
        glBufferData(GL_ARRAY_BUFFER, v_color->size()*sizeof(uchar), v_color->data(), GL_DYNAMIC_DRAW);
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(4,GL_UNSIGNED_BYTE,0,BUFFER_OFFSET(0));

    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vboIds[FACE]);
    //cout<<"render edges: "<<edges->size()<<endl;
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, edges->size()*sizeof(GLuint), edges->data(), GL_DYNAMIC_DRAW);
    // Draw cube geometry using indices from VBO 1
    glDrawElements(GL_LINES, edges->size(), GL_UNSIGNED_INT, BUFFER_OFFSET(0));
    glDisableClientState(GL_COLOR_ARRAY);


    //glPopMatrix();
}

void GLWidget::drawGrid(vector< GLdouble >* vertices,vector<GLuint>* edges,
                        vector< GLdouble >* v_normals,vector< uchar >* v_color)
{
    glBindBuffer(GL_ARRAY_BUFFER, vboIds[VERTEX]);
    glBufferData(GL_ARRAY_BUFFER, vertices->size()*sizeof(GLdouble), vertices->data(), GL_DYNAMIC_DRAW);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3,GL_DOUBLE,0,BUFFER_OFFSET(0));

    if(v_normals!=NULL){
        glBindBuffer(GL_ARRAY_BUFFER, vboIds[V_NORMAL]);
        glBufferData(GL_ARRAY_BUFFER, v_normals->size()*sizeof(GLdouble), v_normals->data(), GL_DYNAMIC_DRAW);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_DOUBLE,0,BUFFER_OFFSET(0));
    }

    if(v_color != NULL){
        glBindBuffer(GL_ARRAY_BUFFER, vboIds[V_COLOR]);
        glBufferData(GL_ARRAY_BUFFER, v_color->size()*sizeof(uchar), v_color->data(), GL_DYNAMIC_DRAW);
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(4,GL_UNSIGNED_BYTE,0,BUFFER_OFFSET(0));

    }

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vboIds[FACE]);
    //cout<<"render edges: "<<edges->size()<<endl;
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, edges->size()*sizeof(GLuint), edges->data(), GL_DYNAMIC_DRAW);
    // Draw cube geometry using indices from VBO 1
    glDrawElements(GL_LINES, edges->size(), GL_UNSIGNED_INT, BUFFER_OFFSET(0));
    glDisableClientState(GL_COLOR_ARRAY);
    //cout<<edges->size()<<' '<<v_color->size()<<' '<<vertices->size()<<endl;


    //glPopMatrix();
}

void GLWidget::drawMeshandGrid(vector< GLdouble >* vertices, vector<GLuint> *edges,vector<GLuint>* faces,
                               vector< GLdouble >* v_normals, vector<uchar> *v_color){

    //glColor3d(1,0,0);
    //cout<<"vertices: "<<vertices<<endl;
    if(vertices==NULL || vertices->size()==0)return;
    glBindBuffer(GL_ARRAY_BUFFER, vboIds[VERTEX]);
    glBufferData(GL_ARRAY_BUFFER, vertices->size()*sizeof(GLdouble), vertices->data(), GL_DYNAMIC_DRAW);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3,GL_DOUBLE,0,BUFFER_OFFSET(0));

    if(v_normals!=NULL && v_normals->size()!=0){
        glBindBuffer(GL_ARRAY_BUFFER, vboIds[V_NORMAL]);
        glBufferData(GL_ARRAY_BUFFER, v_normals->size()*sizeof(GLdouble), v_normals->data(), GL_DYNAMIC_DRAW);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_DOUBLE,0,BUFFER_OFFSET(0));
    }

    if(v_color != NULL && v_color->size()!=0){
        glBindBuffer(GL_ARRAY_BUFFER, vboIds[V_COLOR]);
        glBufferData(GL_ARRAY_BUFFER, v_color->size()*sizeof(uchar), v_color->data(), GL_DYNAMIC_DRAW);
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(4,GL_UNSIGNED_BYTE,0,BUFFER_OFFSET(0));

    }

    if(edges != NULL &&edges->size()!=0){
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vboIds[FACE]);
        //cout<<"render edges: "<<edges->size()<<endl;
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, edges->size()*sizeof(GLuint), edges->data(), GL_DYNAMIC_DRAW);
        // Draw cube geometry using indices from VBO 1
        glDrawElements(GL_LINES, edges->size(), GL_UNSIGNED_INT, BUFFER_OFFSET(0));
    }

    if(faces != NULL && faces->size()!=0){
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vboIds[FACE]);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, faces->size()*sizeof(GLuint), faces->data(), GL_DYNAMIC_DRAW);
        // Draw cube geometry using indices from VBO 1
        glDrawElements(GL_TRIANGLES, faces->size(), GL_UNSIGNED_INT, BUFFER_OFFSET(0));
    }
    glDisableClientState(GL_COLOR_ARRAY);

}

void GLWidget::normalizeAngle(int *angle)
{
    while (*angle < 0)
        *angle += 360 * 16;
    while (*angle > 360 * 16)
        *angle -= 360 * 16;
}

void GLWidget::mouseDoubleClickEvent(QMouseEvent *event){
    lastPos = event->pos();
    cout<<"double click ratio x: " <<float(lastPos.x())/s_width <<"  y: "<<float(lastPos.y())/s_height<<endl;
    cout<<"double click abs x: " <<float(lastPos.x()) <<"  y: "<<float(lastPos.y())<<endl;
    cout<<"double click saved x: " <<s_width <<"  y: "<<s_height<<endl;
    //DrawDblClickLine(float(lastPos.x())/s_width,float(lastPos.y())/s_height);
    GetRayLine(event->x(),event->y());
    if(!isRead)return;
    //float area = c_geometry.UpdateChosenState(point_near,point_far);
    //setText(area);

    updateGL();
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();

    GetRayLine(event->x(),event->y());
    if(isUserAddingConstrainMode){

        if (event->buttons() & Qt::RightButton) {

            if(QApplication::keyboardModifiers () == Qt::ControlModifier){
                GetRayLine(event->x(),event->y());
                ComputeRayIntersectCutplane();
                emit addNewPoint();
                updateGL();
            }
        }
    }
}


void GLWidget::keyPressEvent(QKeyEvent *event){

    if(isplaneOn){
        bool isChangePlane = false;
        double stepi = 0.02;
        if (event->key() == Qt::Key_W){
            planePara[3]+=stepi;
            isChangePlane = true;
        }else if(event->key() == Qt::Key_S){
            planePara[3]-=stepi;
            isChangePlane = true;
        }else if(event->key() == Qt::Key_Up){
            RotatePlane(0,stepi);
            isChangePlane = true;
        }else if(event->key() == Qt::Key_Down){
            RotatePlane(0,-stepi);
            isChangePlane = true;
        }else if(event->key() == Qt::Key_Left){
            RotatePlane(1,stepi);
            isChangePlane = true;
        }else if(event->key() == Qt::Key_Right){
            RotatePlane(1,-stepi);
            isChangePlane = true;
        }
        ComputeCutPlanePoints();
        if(isChangePlane)emit CutPlaneChanges();
        updateGL();
    }
    if(isKeyBoardHideMode){
        int a = event->key() - Qt::Key_0;
        if(a>=0 && a<10) TriggerShowingMat(a);

    }

}


void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    float dx = event->x() - lastPos.x();
    float dy = event->y() - lastPos.y();




    if (event->buttons() & Qt::LeftButton) {

        if(QApplication::keyboardModifiers () == Qt::ShiftModifier){
            m_trans[0] += dx * 0.01f;
            m_trans[1] -= dy * 0.01f;

        }else{
            setXRotation(0.5 * dy);
            setYRotation(0.5 * dx);
            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();
            glRotated(yRot, 0,1.0,0);
            glRotated(xRot, 1.0,0,0);
            glMultMatrixd(m_RotationMatrix); //accummulate the previous transformations
            glGetDoublev(GL_MODELVIEW_MATRIX, m_RotationMatrix); //update
            glPopMatrix();

        }

    }


    lastPos = event->pos();
    if(isUserAddingConstrainMode){

        GetRayLine(event->x(),event->y());
        ComputeRayIntersectCutplane();
        nextEdges = false;
        if (event->buttons() & Qt::RightButton) {
            if(QApplication::keyboardModifiers () == Qt::ControlModifier){
                nextEdges = true;
                emit addNewPoint();
                updateGL();
            }
        }
    }

    updateGL();
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event){

    if(isUserAddingConstrainMode){

        nextEdges = false;
    }


}



void GLWidget::wheelEvent(QWheelEvent *event){
    scale+= event->delta() / 360.0/5.;
    if (scale<0.01)scale = 0.01;
    if (scale>50)scale = 50;
    updateGL();
}

void GLWidget::changeGLwindowsize(int w,int h){
    s_width = w,s_height = h;
    ratio_width = float(r_width)/s_width;ratio_height = float(r_height)/s_height;

}


void GLWidget::changeLineWidth(int a){
    linewidth = a*8./100.;
    updateGL();
}

void GLWidget::changePointSize(int a){
    pointsize = a*16./100.;
    updateGL();
}




void GLWidget::pickFaceViaColor(int x, int y,vector< GLdouble >* vertices, vector<GLuint>* faces, vector<uchar> *v_color){

    unsigned char cc[4];
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glPushAttrib( GL_ALL_ATTRIB_BITS );
    glEnable( GL_DEPTH_TEST );
    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);
    glDisable(GL_ALPHA_TEST);


    qglClearColor(QColor(255,255,255));
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();

    glScalef(scale,scale,scale);
    glTranslated(m_trans[0],m_trans[1],m_trans[2]);

    glMultMatrixd(m_RotationMatrix);



    glBindBuffer(GL_ARRAY_BUFFER, vboIds[VERTEX]);
    glBufferData(GL_ARRAY_BUFFER, vertices->size()*sizeof(GLdouble), vertices->data(), GL_DYNAMIC_DRAW);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3,GL_DOUBLE,0,BUFFER_OFFSET(0));




    glBindBuffer(GL_ARRAY_BUFFER, vboIds[V_COLOR]);
    glBufferData(GL_ARRAY_BUFFER, v_color->size()*sizeof(uchar), v_color->data(), GL_DYNAMIC_DRAW);
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(4,GL_UNSIGNED_BYTE,0,BUFFER_OFFSET(0));





    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vboIds[FACE]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, faces->size()*sizeof(GLuint), faces->data(), GL_DYNAMIC_DRAW);
    // Draw cube geometry using indices from VBO 1
    glDrawElements(GL_TRIANGLES, faces->size(), GL_UNSIGNED_INT, BUFFER_OFFSET(0));



    glPopMatrix();
    glPopAttrib();

    GLint viewport[4];
    glGetIntegerv( GL_VIEWPORT, viewport );
    GLint h = viewport[3];
    GLint w = viewport[2];

    glReadPixels( x, h-y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, cc );

    //    c_surf.CreateNewCluster(c_surf.PickFaceViaColor(cc));

    //pickInd = c_surf.PickFaceViaColor(cc);

    glDisableClientState(GL_COLOR_ARRAY);


}

void GLWidget::exportViewData(){

    //FILE *fp = fopen("/Users/Research/Geometry/RBF/Program/view.dat","w");
    FILE *fp = fopen("../view.dat","w");
    if(fp==NULL)return;

    fwrite(&scale,1,sizeof(double),fp);
    fwrite(m_trans,3,sizeof(double),fp);
    fwrite(m_RotationMatrix,16,sizeof(double),fp);

    fclose(fp);
    cout<<"exportViewData"<<endl;
}

void GLWidget::exportViewData(string filename){
    FILE *fp = fopen(filename.data(),"w");
    if(fp==NULL)return;

    fwrite(&scale,1,sizeof(double),fp);
    fwrite(m_trans,3,sizeof(double),fp);
    fwrite(m_RotationMatrix,16,sizeof(double),fp);

    fclose(fp);
    cout<<"exportViewData"<<endl;
}

void GLWidget::importViewData(string filename){

    FILE *fp = fopen(filename.c_str(),"r");
    //FILE *fp = fopen("../view.dat","r");
    //FILE *fp = fopen("/Users/Research/Geometry/RMFpro/viewdat/biarc.dat","r");
    if(fp==NULL)return;

    fread(&scale,1,sizeof(double),fp);
    fread(m_trans,3,sizeof(double),fp);
    fread(m_RotationMatrix,16,sizeof(double),fp);

    fclose(fp);
    updateGL();
    cout<<"importViewData"<<endl;
}
void GLWidget::importViewData(){

    //FILE *fp = fopen(filename.c_str(),"r");
    FILE *fp = fopen("../view.dat","r");
    //FILE *fp = fopen("/Users/Research/Geometry/RBF/Program/view.dat","r");
    if(fp==NULL)return;

    fread(&scale,1,sizeof(double),fp);
    fread(m_trans,3,sizeof(double),fp);
    fread(m_RotationMatrix,16,sizeof(double),fp);

    fclose(fp);

    load_scale = scale;
    for(int i=0;i<3;++i)load_m_trans[i] = m_trans[i];
    for(int i=0;i<16;++i)load_m_RotationMatrix[i] = m_RotationMatrix[i];


    updateGL();
    cout<<"importViewData"<<endl;
}


void GLWidget::ToggleSmoothTransition(int interp_factor_in){
    glm::dmat4 oriR;
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)oriR[i][j] = load_m_RotationMatrix[i*4+j];

    glm::dmat4 endR(0.0);
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)if(i==j)endR[i][j] = 1.;
    glm::dquat quat_start,quat_end;
    quat_start  = glm::quat_cast(oriR);
    quat_end   = glm::quat_cast(endR);

    // Interpolate half way from original view to the new.
    double interp_factor = double(interp_factor_in)/100.; // 0.0 == original, 1.0 == new

    // First interpolate the rotation
    glm::dquat  quat_interp = glm::slerp (quat_start, quat_end, interp_factor);

    // Then interpolate the translation
    //glm::vec3  pos_interp  = glm::mix   (eye_start,  eye_end,  interp_factor);

    glm::dmat4  view_matrix = glm::mat4_cast (quat_interp); // Setup rotation
    //view_matrix [3]        = glm::vec4 (pos_interp, 1.0);  // Introduce translation

    for(int i=0;i<4;++i)for(int j=0;j<4;++j) m_RotationMatrix[i*4+j] = view_matrix[i][j];


    double endT[3] = {0,1.1,0};

    for(int i=0;i<3;++i)m_trans[i] = (1-interp_factor)*load_m_trans[i] + interp_factor*endT[i];

    double endscale = 0.3;
    scale =  (1-interp_factor)*load_scale + interp_factor*endscale;

    //cout<<"ToggleSmoothTransition"<<endl;
    updateGL();


}

void GLWidget::ToggleSmoothTransition_forVideo(int loadi, int loadj,double interp_factor_in){

    glm::dmat4 oriR;
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)oriR[i][j] = vload_m_RotationMatrix[loadi][i*4+j];

    glm::dmat4 endR;
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)endR[i][j] = vload_m_RotationMatrix[loadj][i*4+j];
    glm::dquat quat_start,quat_end;
    quat_start  = glm::quat_cast(oriR);
    quat_end   = glm::quat_cast(endR);

    // Interpolate half way from original view to the new.
    double interp_factor = double(interp_factor_in); // 0.0 == original, 1.0 == new

    // First interpolate the rotation
    glm::dquat  quat_interp = glm::slerp (quat_start, quat_end, interp_factor);

    // Then interpolate the translation
    //glm::vec3  pos_interp  = glm::mix   (eye_start,  eye_end,  interp_factor);

    glm::dmat4  view_matrix = glm::mat4_cast (quat_interp); // Setup rotation
    //view_matrix [3]        = glm::vec4 (pos_interp, 1.0);  // Introduce translation

    for(int i=0;i<4;++i)for(int j=0;j<4;++j) m_RotationMatrix[i*4+j] = view_matrix[i][j];



    for(int i=0;i<3;++i)m_trans[i] = (1-interp_factor)*vload_m_trans[loadi][i] + interp_factor*vload_m_trans[loadj][i];


    scale =  (1-interp_factor)*vload_scale[loadi] + interp_factor*vload_scale[loadj];

    //cout<<"ToggleSmoothTransition"<<endl;
    updateGL();

}

void GLWidget::LoadViewSeqs(vector<string>&filelists){
    int n = filelists.size();
    vload_m_RotationMatrix.clear();
    vload_m_trans.clear();
    vload_scale.clear();
    vload_m_RotationMatrix.resize(n,vector<double>(16));
    vload_m_trans.resize(n,vector<double>(3));
    vload_scale.resize(n);

    for(int i=0;i<n;++i){

        cout<<filelists[i]<<endl;
        FILE *fp = fopen(filelists[i].data(),"r");
        //FILE *fp = fopen("/Users/Research/Geometry/RMFpro/viewdat/biarc.dat","r");
        if(fp==NULL){
            cout<<"cant open file: "<<filelists[i]<<endl;
            return;
        }


        fread(&(vload_scale[i]),1,sizeof(double),fp);
        fread(vload_m_trans[i].data(),3,sizeof(double),fp);
        fread(vload_m_RotationMatrix[i].data(),16,sizeof(double),fp);

        fclose(fp);

    }


}
void GLWidget::ActivatePlayVideoMode(int interval){








    timer_forplay.setInterval(interval);
    timer_forplay.setSingleShot(false);

    accTimeout = 0;
    timer_forplay.start();

}

void GLWidget::ProcessTimeout(){
    int ival = 50;
    if(accTimeout>=(vload_scale.size()-1)*ival){
        timer_forplay.stop();
        return;
    }

    int fac = accTimeout%ival;
    int n = (accTimeout-fac)/ival;

    //cout<<n<<' '<<n+1<<endl;
    ToggleSmoothTransition_forVideo(n, n+1,double(fac)/double(ival));

    ++accTimeout;



}
