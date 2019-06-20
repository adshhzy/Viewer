#include "geo_pc.h"
#include "readers.h"
#include <stdio.h>
#include<iostream>
#include<fstream>
#include<limits>
#include <functional>
#include <eigen3/Eigen/Geometry>

namespace n_rf {

int Multi_PointCloud::Import_Multi_PointCloud(vector<vector<double>>* vertices_in, vector<vector<double>>* vn_in){

    numofPC = vertices_in->size();

    PC_buffer.resize(numofPC);

    for(int i=0;i<numofPC;++i){

        PC_buffer[i].ImportPointCloud(vertices_in->at(i),vn_in->at(i));
        if(i!=0)PC_buffer[i].Set_Avedist(PC_buffer[0].Get_Avedist());
    }

    return numofPC;

}

int Multi_PointCloud::Import_Multi_PointCloud(vector<vector<double>>* vertices_in, vector<vector<double>>* vn_in, vector<vector<double>>* vn_cmp){

    numofPC = vertices_in->size();

    PC_buffer.resize(numofPC);

    for(int i=0;i<numofPC;++i){

        PC_buffer[i].ImportPointCloud(vertices_in->at(i),vn_in->at(i), vn_cmp->at(i));
        if(i!=0)PC_buffer[i].Set_Avedist(PC_buffer[0].Get_Avedist());
    }

    return numofPC;

}

int Multi_PointCloud::Import_Multi_PointCloud(vector<string>*pnames,vector<vector<double>>* vertices_in, vector<vector<double>>* vn_in, vector<vector<double>>* vn_cmp){

    pc_names = *pnames;
    pc_nameind.clear();
    for(int i=0;i<pc_names.size();++i)pc_nameind[pc_names[i]] = i;
    Import_Multi_PointCloud(vertices_in,vn_in,vn_cmp);
}

int Multi_PointCloud::Import_Multi_PointCloud(vector<vector<double>>* vertices_in, vector<vector<double>>* vn_in, vector<vector<uint>>* pc_graph_in, vector<vector<double>>* pc_graphw_in){

    numofPC = vertices_in->size();

    PC_buffer.clear();
    PC_buffer.resize(numofPC);

    for(int i=0;i<numofPC;++i){

        PC_buffer[i].ImportPointCloud(&(vertices_in->at(i)),&(vn_in->at(i)), &(pc_graph_in->at(i)), &(pc_graphw_in->at(i)));
        if(i!=0)PC_buffer[i].Set_Avedist(PC_buffer[0].Get_Avedist());
    }

    return numofPC;



}

int Multi_PointCloud::Import_Multi_Graph(vector<vector<uint>>* pc_graph_in, vector<vector<double>>* pc_graphw_in){

    if(numofPC!=pc_graph_in->size())return -1;
    if(numofPC!=pc_graphw_in->size())return -1;

    for(int i=0;i<numofPC;++i){
        PC_buffer[i].ImportGraph(&pc_graph_in->at(i),&pc_graphw_in->at(i));
    }
    return numofPC;

}

int Multi_PointCloud::Get_PCIndex(string fname){
    if(pc_nameind.find(fname)!=pc_nameind.end())return pc_nameind[fname];
    else return -1;
}

void Multi_PointCloud::FlipNormal(){

    for(auto &a:PC_buffer)a.FlipNormal();

}

void Multi_PointCloud::BuildDisplay(infoPCDisp info){

    for(auto &a:PC_buffer)a.BuildDisplay(info);

    display_vertices.resize(PC_buffer.size());
    display_edges.resize(PC_buffer.size());
    display_faces.resize(PC_buffer.size());
    display_normal.resize(PC_buffer.size());
    display_vcolor.resize(PC_buffer.size());


    for(int i=0;i<PC_buffer.size();++i){
        //if(i!=1)continue;
        display_vertices[i] = PC_buffer[i].getDisplayVertices();
        display_edges[i] = PC_buffer[i].getDisplayEdges();
        display_faces[i] = PC_buffer[i].getDisplayFaces();
        display_normal[i] = PC_buffer[i].getDisplayVerticesNormal();
        display_vcolor[i] = PC_buffer[i].getDisplayColor();
    }


}





}//n_rf
