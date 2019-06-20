QMAKE_MAC_SDK = macosx10.14
QT       +=  core gui
QT       += opengl widgets
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TEMPLATE = app
#CONFIG += console
#CONFIG -= app_bundle
#CONFIG -= qt
QT           += opengl

QMAKE_CFLAGS_RELEASE -= -O
QMAKE_CFLAGS_RELEASE -= -O1
QMAKE_CFLAGS_RELEASE -= -O2

QMAKE_CFLAGS_RELEASE *= -O3
CONFIG += c++11
CONFIG +=m64
CONFIG +=-pthread
DEPENDPATH += . $$PWD/../include
INCLUDEPATH += . $$PWD/../include

INCLUDEPATH += $$/opt/local/include
DEPENDPATH += $$/opt/local/include

INCLUDEPATH += /Users/Research/Geometry/RBF/Program/include

HEADERS  += RBFCore/Polygonizer.h \
    mymesh/contour.h \
    mymesh/geo_curv.h \
    mymesh/InfoStruct.h \
    mymesh/my_mesh.h \
    mymesh/readers.h \
    mymesh/UnionFind.h \
    mymesh/utility.h \
    GLInterface/glwidget.h \
    MainWindow/mainwindow.h \
    mymesh/geo_sur.h \
    mymesh/geo_pc.h \
    mymesh/a_multisur.h

SOURCES += main.cpp \
    mymesh/contour.cpp \
    mymesh/contour_fcm.cpp \
    mymesh/contour_geo.cpp \
    mymesh/geo_curv.cpp \
    mymesh/my_mesh.cpp \
    mymesh/readers.cpp \
    mymesh/UnionFind.cpp \
    mymesh/contour_rbf.cpp \
    GLInterface/glwidget_core.cpp \
    GLInterface/glwidget_interface.cpp \
    MainWindow/mainwindow.cpp \
    MainWindow/MainWindow_Control.cpp \
    mymesh/a_geosur_display.cpp \
    mymesh/a_geosur_interface.cpp \
    mymesh/a_geosur.cpp \
    mymesh/geo_curv_display.cpp \
    mymesh/geo_pc.cpp \
    mymesh/geo_multipc.cpp \
    mymesh/a_multisur.cpp


FORMS    += \
    QtForms/mainwindow.ui \



INCLUDEPATH += $$/usr/local/include/

INCLUDEPATH += $$/opt/local/include/eigen3/

