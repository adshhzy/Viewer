#include "contour.h"
#include "readers.h"
#include "geo_curv.h"
#include<iostream>
#include <assert.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <set>
#include <unordered_map>

template<class BidiIter>
BidiIter random_unique(BidiIter begin, BidiIter end, size_t num_random) {
    size_t left = std::distance(begin, end);
    while (num_random--) {
        BidiIter r = begin;
        std::advance(r, rand()%left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}


void Contour::GetEdgesNormal(vector<double>&enormal, bool isstackup){
    if(!isstackup)enormal.clear();


    for(int i=0;i<m1.size();++i){
        auto p_ev = ev_begin(i);
        double edgeMid[3], edgeVec[3], edgeNor[3];
        MyUtility::_VerticesMidpoint(v_begin(p_ev[0]),v_begin(p_ev[1]),edgeMid);

        MyUtility::minusVec(v_begin(p_ev[1]),v_begin(p_ev[0]),edgeVec);

        MyUtility::cross(plane_para,edgeVec,edgeNor);

        MyUtility::normalize(edgeNor);

        if(m1[i]!=0)MyUtility::inversevec(edgeNor,edgeNor);

        for(int j=0;j<3;++j)enormal.push_back(edgeNor[j]);
    }


}

void Contour::EdgeMidPointToNormal(double eps, vector<double>&re_points,vector<int>&labels){

    //re_points.clear();
    for(int i=0;i<m1.size();++i)if((double)rand() / RAND_MAX > 0.8){
        auto p_ev = ev_begin(i);
        double edgeMid[3], edgeVec[3], edgeNor[3], edgeOut[3];
        MyUtility::_VerticesMidpoint(v_begin(p_ev[0]),v_begin(p_ev[1]),edgeMid);

        MyUtility::minusVec(v_begin(p_ev[1]),v_begin(p_ev[0]),edgeVec);

        MyUtility::cross(plane_para,edgeVec,edgeNor);

        MyUtility::normalize(edgeNor);

        if(m1[i]!=0)MyUtility::inversevec(edgeNor,edgeNor);

        MyUtility::weightedAddVec(1.,eps,edgeMid,edgeNor,edgeOut);
        for(int j=0;j<3;++j)re_points.push_back(edgeOut[j]);
        labels.push_back(-1);

        MyUtility::weightedAddVec(1.,-eps*0.6,edgeMid,edgeNor,edgeOut);
        for(int j=0;j<3;++j)re_points.push_back(edgeOut[j]);
        labels.push_back(1);
    }


}



bool CrossSections::ReadCrossSections_Contourpp(string filename){
    ifstream fin(filename.data());
    if(fin.fail()){
        cout<<"Fail to open input file: "<<filename<<endl;
        return false;
    }

    fin>>n_CrossSec;

    //n_CrossSec = 2;


    CSs.resize(n_CrossSec);
    for(int i=0;i<n_CrossSec;++i){
        CSs[i].ReadContour_Contourpp(fin);
    }
    CalSelfProperty();

    return true;

}

int CrossSections::WriteCrossSections_Contourpp(string filename){
    string outCname = filename + string(".contourpp");
    ofstream outer(outCname.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not open the Contour file " << outCname << endl;
        return -1;
    }

    outer<<std::setprecision(12);
    outer<<CSs.size()<<endl;
    for(auto &b:CSs)b.WriteContour_Contourpp(outer);
    outer.close();
    return 0;

}

void CrossSections::Set_EdgeNormal(vector<uint>&eind, vector<double>&vnormal){
    vector<double>out_Vs;

    vector<uint>tmpedges;
    Stackup(out_Vs,tmpedges);

    n_rf::Curve curveNet;
    curveNet.ImportCurve(out_Vs, tmpedges);
    curveNet.Set_EdgeNormal(eind,vnormal);

    vector<int>e2C;
    vector<int>e2Ce;

    StackUp_EdgeMapper(e2C,e2Ce);

    for(Contour &a:CSs)a.edge_normal.resize(a.n_edges*3);
    int n_totaledges = tmpedges.size()/2;
    auto p_en = curveNet.edge_normal.data();
    for(int i=0;i<n_totaledges;++i){
        copyVec(p_en+i*3, CSs[e2C[i]].edge_normal.data()+e2Ce[i]*3);
    }

}

void CrossSections::StackUp_EdgeMapper(vector<int>&e2C, vector<int>&e2Ce){
    e2C.clear();
    e2Ce.clear();

    for(int i=0;i<n_CrossSec;i++)for(int j=0;j<CSs[i].n_edges;++j){
        e2C.push_back(i);e2Ce.push_back(j);
    }

}

void Contour::Set_EdgeNormal(vector<double>&enormal){



}

bool Contour::ReadContour_Contourpp(ifstream &in){
    for(int i=0;i<4;++i)in>>plane_para[i];
    in>>n_vertices;
    in>>n_edges;

    vertices.resize(n_vertices*3);
    edges.resize(n_edges*2);
    m1.resize(n_edges);
    m2.resize(n_edges);
    edge_normal.resize(n_edges*3);

    for(int i=0;i<vertices.size();++i)in>>vertices[i];
    for(int i=0;i<n_edges;++i){
        in>>edges[i*2];
        in>>edges[i*2+1];
        in>>m1[i];
        in>>m2[i];
        auto p_en = edge_normal.data()+i*3;
        for(int j=0;j<3;++j)in>>p_en[j];
    }

    set<int>mmm;
    for(auto a:m1)mmm.insert(a);
    for(auto a:m2)mmm.insert(a);

    if(0){
        int newind = 0;
        unordered_map<int, int>mapper;
        for(auto a:mmm)mapper[a] = newind++;
        for(auto &a:m1)a = mapper[a];
        for(auto &a:m2)a = mapper[a];
        n_materials = mmm.size();
    }
    else n_materials = *max_element(mmm.begin(),mmm.end())+1;
	
	return true;
}
bool Contour::WriteContour_Contourpp(ofstream &out){

    for(int i=0;i<4;++i)out<<plane_para[i]<<' ';out<<endl;
    out<<n_vertices<<' ';
    out<<n_edges<<endl;
    for(int i=0;i<n_vertices;++i){
        int ind = i*3;
        for(int j=0;j<3;++j)out<<vertices[ind+j]<<' ';
        out<<endl;
    }
    for(int i=0;i<n_edges;++i){
        out<<edges[i*2]<<' ';
        out<<edges[i*2+1]<<' ';
        out<<m1[i]<<' ';
        out<<m2[i]<<' ';
        auto p_en = edge_normal.data()+i*3;
        out<<p_en[0]<<' '<<p_en[1]<<' '<<p_en[2]<<endl;
    }
	return true;
}


void curveNetDownsampling(vector<double>&in_vertices,vector<unsigned int>&in_edges,int nsample, vector<double>&resample_cv,vector<unsigned int>&p2edgeid){



    n_rf::Curve curveNet;


}


/*********************************************************************************************************/
/*********************************************************************************************************/
void curveNetresampling(vector<double>&in_vertices,vector<unsigned int>&in_edges,int nsample, vector<double>&resample_cv,vector<unsigned int>&p2edgeid){



    n_rf::Curve curveNet;

    int n_edges = in_edges.size()/2;
    int n_vertices = in_vertices.size()/3;

    assert(nsample>n_vertices);

    curveNet.ImportCurve(in_vertices, in_edges);

    vector<double>&edge_len = curveNet.edge_len;
    double ave_len = curveNet.ave_elen;
    double total_len = ave_len * n_edges;

    int n_new_vertices = nsample - n_vertices;

    vector<pair<double,uint>>edgelens;


    for(uint i=0;i<n_edges;++i)edgelens.push_back(make_pair(edge_len[i], i));

    sort(edgelens.begin(),edgelens.end(),[](const pair<double,uint> & a, const pair<double,uint> & b)
    {
        return a.first > b.first;
    });


    vector<int>pick_slots(n_edges, 0);

    double pperlen = n_new_vertices/total_len;

    int totalpointfirstround = 0;
    for(uint i=0; i<n_edges; ++i)totalpointfirstround += (pick_slots[i] = int(ceil(pperlen*edge_len[i])));

    int res = totalpointfirstround - n_new_vertices;



    if(res>0){
        vector<unsigned int> indices(n_edges);
        for(uint i=0;i<n_edges;++i)indices[i] = i;
        random_unique(indices.begin(), indices.end(),res);

        for(uint i=0;i<res;++i)--pick_slots[indices[i]];
    }



    vector<double>tmp_resample_cv = in_vertices;
    vector<unsigned int>tmp_p2edgeid(n_vertices);
    for(int i=0;i<n_vertices;++i){
        tmp_p2edgeid[i] = curveNet.ve_begin(i)[0];
    }

    auto sampleseg = [&tmp_resample_cv, &tmp_p2edgeid](double *p1,double *p2,int n, unsigned int eid){

        double vec[3];
        for(int i=0;i<3;++i)vec[i] = (p2[i] - p1[i])/(n+1);

        for(int j=1;j<=n;++j)for(int i=0;i<3;++i)tmp_resample_cv.push_back(p1[i]+vec[i]*j);
        for(int j=1;j<=n;++j)tmp_p2edgeid.push_back(eid);

    };

    for(int i=0;i<n_edges;++i)if(pick_slots[i]>0){

        auto ev = curveNet.ev_begin(i);
        auto p_v1 = curveNet.v_begin(ev[0]);
        auto p_v2 = curveNet.v_begin(ev[1]);

        sampleseg(p_v1,p_v2,pick_slots[i],i);


    }


    vector<int>inds(nsample);

    for(int i=0;i<nsample;++i)inds[i] = i;
    //random_shuffle(inds.begin(), inds.end());


    resample_cv.clear();
    p2edgeid.clear();
    //for(auto a:tmp_resample_cv)resample_cv.push_back(a);
    for(int i=0;i<nsample;++i){
        auto p_v = tmp_resample_cv.data()+inds[i]*3;
        for(int j=0;j<3;++j)resample_cv.push_back(p_v[j]);
        p2edgeid.push_back(tmp_p2edgeid[inds[i]]);

    }


}



