#include "geo_pc.h"
#include "readers.h"
#include <stdio.h>
#include<iostream>
#include<fstream>
#include<limits>
#include <functional>
#include <eigen3/Eigen/Geometry>

extern double controlBallsize;
extern double controlTubesize;
double controlDisksize = 1.0;
void ComputeAngleAxisMatrix(Eigen::AngleAxisd &m_aa,double *oriVec,double *newVec);
namespace n_rf {

bool PointCloud::isload = false;
Mesh PointCloud::sphere;
Mesh PointCloud::cylinder;
Mesh PointCloud::cone;
Mesh PointCloud::disk;

PointCloud::PointCloud(){

    n_vertices=0;
    load();

}


bool PointCloud::load(){
    if(!isload){
        disk.createToy(4);
        sphere.createToy(3);

        cylinder.createToy(2);

        cone.createToy(1);
        cone.ReScale(2.5,2.5,0.5);
        cylinder.ReScale(0.5,0.5,1.0);

        sphere.BuildNeighborTable();
        cylinder.BuildNeighborTable();
        cone.BuildNeighborTable();
        disk.BuildNeighborTable();

        sphere.ComputeFaceNormal(true);
        cylinder.ComputeFaceNormal(true);
        cone.ComputeFaceNormal(true);
        disk.ComputeFaceNormal(true);

        isload = true;
    }
    //sphere.readfile(string("/Users/Research/Geometry/RMFpro/cube.off"));
    return true;
}

void PointCloud::setparameters(){

    n_vertices = vertices.size()/3;

    if(n_vertices==1)avedist = 0.05;
    else{
        avedist = 0;
        for(int i=0;i<n_vertices;++i){
            auto p_v = v_begin(i);
            double mindist = numeric_limits<double>::max();
            for(int j=0;j<n_vertices;++j)if(i!=j){
                mindist = min(mindist,_VerticesDistance(p_v,v_begin(j)));
            }
            avedist+=mindist;
        }
        avedist/=n_vertices;
    }

    //cout<<"avedist: "<<avedist<<endl;
    //cout<<vnormals.size()<<endl;
    int nc00 = 0;
    vnormals_uninorm = vnormals;
    auto p_vno = vnormals_uninorm.data();
    for(int i=0;i<n_vertices;++i){
        auto p_vn = p_vno+i*3;
        double nn = len(p_vn);
        if(nn<1e-16){
            p_vn[0] = 1;p_vn[1] = 0;p_vn[2] = 0;
            ++nc00;
        }else{
            double invsqrtnn = 1/sqrt(nn);
            for(int j=0;j<3;++j)p_vn[j]*=invsqrtnn;
        }

    }
    //cout<<"nc00: "<<nc00<<endl;
    if(vnormals_cmp.size()==vnormals.size()){
        vnormals_uninorm_cmp = vnormals_cmp;
        auto p_vno = vnormals_uninorm_cmp.data();
        for(int i=0;i<n_vertices;++i){
            auto p_vn = p_vno+i*3;
            double nn = len(p_vn);
            if(nn<1e-16){
                p_vn[0] = 1;p_vn[1] = 0;p_vn[2] = 0;
            }else{
                double invsqrtnn = 1/sqrt(nn);
                for(int j=0;j<3;++j)p_vn[j]*=invsqrtnn;
            }
        }
        Compute_CMP_Diff();
    }


}

void PointCloud::Set_Avedist(double avedist_in){
    avedist = avedist_in;
}

double PointCloud::Get_Avedist(){
    return avedist;
}


bool PointCloud::ReadPointCloud(string filename){



	return true;

}

void PointCloud::FlipNormal(){

    for(auto &a:vnormals)a=-a;
    for(auto &a:vnormals_uninorm)a=-a;

}

void PointCloud::Compute_CMP_Diff(){


    if(vnormals_uninorm_cmp.size()!=vnormals_uninorm.size()){
        cout<<"vnormals_uninorm_cmp not set"<<endl;
        return;
    }

    cmp_diff.resize(n_vertices);

    double acc = 0.;
    for(int i=0;i<n_vertices;++i){
        cmp_diff[i] = MyUtility::dot(vnormals_cmp.data()+i*3, vnormals.data()+i*3);
        acc+=cmp_diff[i];
    }
    if(acc>0){
        for(auto &a:cmp_diff)a = acos(a);
    }else{
        for(auto &a:cmp_diff)a = acos(-a);
    }
    for(auto &a:cmp_diff)a = my_PI - a;

    //double maxk =




}


bool PointCloud::ImportPointCloud(vector<double>& vertices_in, vector<double>& vn_in){

    vertices = vertices_in;

    vnormals = vn_in;

    setparameters();
	return true;
}

bool PointCloud::ImportPointCloud(vector<double>* vertices_in, vector<double>* vn_in){

    vertices = *vertices_in;

    vnormals = *vn_in;

    setparameters();
	return true;
}

bool PointCloud::ImportPointCloud(vector<double>* vertices_in, vector<double>* vn_in, vector<double>* vn_in_cmp){

    vertices = *vertices_in;

    vnormals = *vn_in;

    vnormals_cmp = *vn_in_cmp;

    setparameters();
	return true;
}

bool PointCloud::ImportPointCloud(vector<double>& vertices_in, vector<double>& vn_in, vector<double>& vn_in_cmp){
    vertices = vertices_in;

    vnormals = vn_in;

    vnormals_cmp = vn_in_cmp;

    setparameters();
    return true;

}

bool PointCloud::ImportPointCloud(vector<double>* vertices_in, vector<double>* vn_in, vector<uint>* pc_graph_in, vector<double> *pc_graphw_in){
    vertices = *vertices_in;
    vnormals = *vn_in;
    pc_graph = *pc_graph_in;
    pc_graph_weight = *pc_graphw_in;

    setparameters();
	return true;
}

bool PointCloud::ImportPointCloud(vector<double>* vertices_in, vector<double>* vn_in, vector<double>* vn_in_cmp, vector<uint>* pc_graph_in, vector<double>* pc_graphw_in){
    vertices = *vertices_in;
    vnormals = *vn_in;
    pc_graph = *pc_graph_in;
    pc_graph_weight = *pc_graphw_in;
    vnormals_cmp = *vn_in_cmp;
    setparameters();
    return true;
}


bool PointCloud::ImportGraph(vector<uint>* pc_graph_in, vector<double>* pc_graphw_in){
    pc_graph = *pc_graph_in;
    pc_graph_weight = *pc_graphw_in;
	return true;
}

void PointCloud::ClearDisplayBuffer(){

    display_vertices.clear();
    display_vnormals.clear();
    display_edges.clear();
    display_faces.clear();
    display_color.clear();

}



void PointCloud::BuildDisplay(infoPCDisp info){

    BuildDisplay(info.point_method,info.isDoubleside,info.isNormal,info.isUnitNormal,info.nor_len,info.isgraph,info.plotpercentage,info.transparent_cut);
}

void PointCloud::BuildDisplay(int point_method, bool isdoublesidec, bool isshownormal, bool isunitnormal, double norlen, int graphmethod, double plotpercentage, double transparent_cut){

    //cout<<n_vertices<<endl;
    if(n_vertices==0){
        ClearDisplayBuffer();
        return;
    }
    double sphere_scale = 0.80;
    double tube_scale = 0.20;
    double disk_scale = 0.80;


    sphere_scale*=avedist * disscale * controlBallsize;
    disk_scale*=avedist * disscale * controlDisksize;
    tube_scale*=avedist * disscale * controlTubesize;

    if(point_method==0){if(sphere_scale==sphere_scale)sphere.ReScale_uniform(sphere_scale);}
    else{if(disk_scale==disk_scale)disk.ReScale_uniform(disk_scale);}
    //cout<<"isuseball: "<<avedis<<' '<<sphere_scale<<endl;




    display_vertices = vertices;
    //display_edges = edges;
    display_faces.clear();
    //rvector = vertices;

    display_vnormals=display_vertices;



    double *p_spv = NULL, *p_v = NULL,*p_vn = NULL;
    uint *p_sf = NULL;

    //cout<<"ppppp: "<<display_vertices.size()<<endl;

    if(point_method==0){
        auto n_fs = sphere.n_faces*3;
        auto n_vs = sphere.n_vertices*3;
        for(int iv = 0;iv<n_vertices;++iv){

            int offset = display_vertices.size()/3;
            p_v = v_begin(iv);

            for(int j = 0;j<sphere.n_vertices;++j){
                p_spv = sphere.v_begin(j);
                for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+p_spv[k]);
            }
            p_sf = sphere.fv_begin(0);
            for(int j = 0;j<n_fs;++j){
                display_faces.push_back(offset+p_sf[j]);
            }
            p_vn = sphere.vnor_begin(0);
            for(int j = 0;j<n_vs;++j){
                display_vnormals.push_back(p_vn[j]);
            }
        }
    }else {
        double vec[3],nz[3] = {0,0,1};
        Eigen::AngleAxisd t;
        vector<Eigen::Vector3d>ori(disk.n_vertices);
        vector<Eigen::Vector3d>tranformed(disk.n_vertices);
        for(int j = 0;j<disk.n_vertices;++j){
            p_spv = disk.v_begin(j);
            for(int k=0;k<3;++k)ori[j](k) = p_spv[k];
        }
        vector<Eigen::Vector3d>orinor(disk.n_vertices);
        for(int j = 0;j<disk.n_vertices;++j){
            p_spv = disk.vnor_begin(j);
            for(int k=0;k<3;++k)orinor[j](k) = p_spv[k];
        }
        for(int iv = 0;iv<n_vertices;++iv){

            int offset = display_vertices.size()/3;
            p_v = v_begin(iv);
            ComputeAngleAxisMatrix(t,nz,vnu_begin(iv));
            Eigen::Transform<double,3,Eigen::Affine> t_m(t);
            if(point_method==2 && cmp_diff.size()==n_vertices)t_m.scale(Eigen::Vector3d(cmp_diff[iv],cmp_diff[iv],cmp_diff[iv]));


            copyVec(vnu_begin(iv),vec);
            for(int i=0;i<3;++i)vec[i] = -vec[i];

            for(int j = 0;j<disk.n_vertices;++j){
                tranformed[j] = t_m*ori[j];
            }

            for(int i=0;i<3;++i)vec[i] *= 1e-4;
            for(int j = 0;j<disk.n_vertices;++j){
                for(int k=0;k<3;++k)display_vertices.push_back(tranformed[j](k)+p_v[k] - vec[k]);
            }


            for(int j = 0;j<disk.n_vertices;++j){
                for(int k=0;k<3;++k)display_vertices.push_back(tranformed[j](k)+p_v[k] + vec[k]);
            }


            for(int j = 0;j<disk.n_vertices;++j){
                tranformed[j] = t_m*orinor[j];
            }

            for(int j = 0;j<disk.n_vertices;++j){
                for(int k=0;k<3;++k)display_vnormals.push_back(tranformed[j](k));
            }
            for(int j = 0;j<disk.n_vertices;++j){
                for(int k=0;k<3;++k)display_vnormals.push_back(-tranformed[j](k));
            }

            p_sf = disk.fv_begin(0);
            for(int j = 0;j<disk.n_faces*3;++j){
                display_faces.push_back(offset+p_sf[j]);
            }

            offset+=disk.n_vertices;
            for(int j = 0;j<disk.n_faces*3;++j){
                display_faces.push_back(offset+p_sf[j]);
            }

        }


    }

    unsigned char red[4] = {225,0,0,255};
    unsigned char blue[4] = {100,100,225,255};
    unsigned char green[4] = {100,255,100,255};
    unsigned char gray[4] = {215,215,215,255};
    unsigned char yellow[4] = {255,200,0,255};
    unsigned char ccccoooollll[4] = {225,225,225,255};
    unsigned char* pcolor = gray;

    display_color.clear();
    for (int i = 0; i < n_vertices; i++)for(int j =0;j<4;++j)display_color.push_back(red[j]);

    if(point_method==0){
        for(int iv = 0;iv<n_vertices;++iv){
            for(int k =0;k<sphere.n_vertices;++k)for(int j =0;j<4;++j)display_color.push_back(pcolor[j]);
        }

        //for (int i = 0; i <  sphere.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(gray[j]);
        //for (int i = sphere.n_vertices; i < n_vertices * sphere.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(gray[j]);
        //for (int i = 0; i < n_edges * cylinder.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(pcolor[j]);
    }else{
        if(isdoublesidec){
            for(int iv = 0;iv<n_vertices;++iv){
                for(int k =0;k<disk.n_vertices;++k)for(int j =0;j<4;++j)display_color.push_back(green[j]);
                for(int k =0;k<disk.n_vertices;++k)for(int j =0;j<4;++j)display_color.push_back(blue[j]);
            }
        }else {
            for(int iv = 0;iv<n_vertices;++iv){
                for(int k =0;k<disk.n_vertices;++k)for(int j =0;j<4;++j)display_color.push_back(yellow[j]);
                for(int k =0;k<disk.n_vertices;++k)for(int j =0;j<4;++j)display_color.push_back(yellow[j]);
            }
        }
    }

    //cout<<"ADsada"<<endl;
    display_edges.clear();
    if(isshownormal && vnormals.size()==vertices.size()){
        int offset = display_vertices.size()/3;
        double lenscale = 0.005*norlen;
        double newvv[3];
        for(int i=0;i<n_vertices;++i){

            weightedAddVec(1.0,lenscale,v_begin(i),isunitnormal?vnu_begin(i):vn_begin(i),newvv);
            for(int j=0;j<3;++j)display_vertices.push_back(newvv[j]);
            display_edges.push_back(i);
            display_edges.push_back(i+offset);
            for(int j=0;j<3;++j)display_vnormals.push_back(newvv[3]);
            for(int j =0;j<4;++j)display_color.push_back(red[j]);
        }
        if(!isdoublesidec){
            offset = display_vertices.size()/3;
            for(int i=0;i<n_vertices;++i){

                weightedAddVec(1.0,-lenscale,v_begin(i),isunitnormal?vnu_begin(i):vn_begin(i),newvv);
                for(int j=0;j<3;++j)display_vertices.push_back(newvv[j]);
                display_edges.push_back(i);
                display_edges.push_back(i+offset);
                for(int j=0;j<3;++j)display_vnormals.push_back(newvv[3]);
                for(int j =0;j<4;++j)display_color.push_back(red[j]);
            }
        }
    }


    if(graphmethod!=0)if(pc_graph_weight.size()!=0 && pc_graph.size()!=0){

        int ind_offset = display_vertices.size()/3;
        //display_vertices.insert(display_vertices.end(), vertices.begin(),vertices.end());
        //display_vnormals.insert(display_vnormals.end(), vertices.begin(),vertices.end());

        //display_color.clear();
        //for (int i = 0; i < n_vertices; i++)for(int j =0;j<4;++j)display_color.push_back(blue[j]);
        unsigned char t_blue[4] = {0,0,225,255};
        int nplot = fmin(pc_graph_weight.size(),plotpercentage * n_vertices * 2.);
        double cut_thres = pc_graph_weight[0]*transparent_cut;
        //cout<<"cut_thres: "<<cut_thres<<' '<<transparent_cut<<endl;
        if(graphmethod==2)nplot = pc_graph_weight.size();
        for(int i=0;i<nplot;++i){
            auto p_ev = pc_graph.data()+i*2;
            if(graphmethod==3)t_blue[3] = 255;
            else if(graphmethod==1)t_blue[3] = (uchar)fmax(0, fmin(255, pc_graph_weight[i]/cut_thres*255));
            else if(graphmethod==2){
                t_blue[3] = 255;
                t_blue[2] = (uchar)fmax(0, fmin(255, pc_graph_weight[i]/cut_thres*255));
                t_blue[1] = 255 - t_blue[2];
            }

            //cout<<int(t_blue[3])<<' ';
            auto p_v = v_begin(p_ev[0]);
            for(int j=0;j<3;++j)display_vertices.push_back(p_v[j]);
            for(int j=0;j<3;++j)display_vnormals.push_back(p_v[j]);
            p_v = v_begin(p_ev[1]);
            for(int j=0;j<3;++j)display_vertices.push_back(p_v[j]);
            for(int j=0;j<3;++j)display_vnormals.push_back(p_v[j]);


            display_edges.push_back(ind_offset+i*2);
            display_edges.push_back(ind_offset+i*2+1);

            for (int i = 0; i < 2; i++)for(int j =0;j<4;++j)display_color.push_back(t_blue[j]);
        }
    }



}

void PointCloud::BuildDisplay(vector<double>&vcolordegree, int method){





}



}//n_rf
