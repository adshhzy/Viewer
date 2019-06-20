#include <iostream>

#include "MainWindow/mainwindow.h"

using namespace std;


int main(int argc, char *argv[])
{


    QApplication a(argc, argv);
    MainWindow w;
    w.GetArg(argc, argv);
    w.show();

    return a.exec();


}
