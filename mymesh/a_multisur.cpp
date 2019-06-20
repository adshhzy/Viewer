#include "a_multisur.h"
#include <stdio.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<limits>
#include<cstdlib>
#include <functional>
#include<time.h>
#include<string>
#include <sstream>
#include <iterator>
#include<set>
#include<map>
#include "readers.h"
#include<random>


double randomdouble();


namespace n_rf {

void SufStructure::Construct(const vector<double>&vertices,const vector<unsigned int>&faces2vertices,
                             const vector<int>&facesMat, const vector<int>&CtrEdges,
                             const vector<vector<int>>&in_mappingsToGlobal){

    Vs = vertices;
    Fs = faces2vertices;
    FMs = facesMat;
    Ctrs = CtrEdges;
    mappingsToGlobal = in_mappingsToGlobal;
    nCell = mappingsToGlobal.size();
}


int MultiSurface::Initialize(string filename, bool defaultpath, bool isAutoSmooth, bool isSetnMat, int nMat, bool isCtrl, double lscale, double *pcenter, int colormethod, char* inputcolor){
    clearup();
    bool rere = false;
    if(defaultpath){

        ReadFile("/Users/Research/Geometry/MM/OutputSuf.suf");
        //ReadFile("/Users/Research/Geometry/MM/result/chickenheart/suf_0_5.suf");

        rere = true;
    }
    else{
        rere = ReadFile(filename);
    }
    if(!rere)return -1;

    mixSurf.setparameters();

    mixSurf.BuildNeighborTable_nonmanifold();

    if(ext ==".csl"){

        crossSection.ReScale_Uniform(this->lscale,this->pcenter);
        //cout<<"crossSection: "<<this->lscale<<endl;

        //mixSurf.GetRescaleInfo(this->lscale,this->pcenter);
        //cout<<"mixSurf: "<<this->lscale<<endl;

        mixSurf.ReScale(this->lscale,this->pcenter);

    }
    else if(ext!=".suf")mixSurf.ReScale_uniform(1.0);

    if(ext==".suf"){


        if(!isCtrl){

            GetScaleInfoViaCrossSection(this->lscale,this->pcenter);
            mixSurf.ReScale(this->lscale,this->pcenter);



            //mixSurf.ReScale_uniform(1.0);
        }
        else mixSurf.ReScale(lscale,pcenter);



        if(isAutoSmooth){
            //ComputeSmallMVnormal();
            //mixSurf.Fairing(true,false,false,50,0.5,0.1,&staticV);
            ComputeSmallMVnormal();
            SurfaceFFFairing(0.5,50,&staticV);
        }


        if(isSetnMat)numofSurf = nMat;



    }

    if(ext ==".csl")ReAllocation(false,colormethod,inputcolor,false);
    else ReAllocation(false,colormethod,inputcolor);
    ComputeSmallMVnormal();

    cout<<"number of Materia: "<<numofSurf<<endl;

    return numofSurf;

}

int MultiSurface::InitializeForCutPlane(string filename,bool issmooth){
    clearup();
    bool rere = false;

    //rere = readSufFile(filename,true);
    rere = ReadFile(filename);
    if(!rere)return -1;

    {
        mixSurf.setparameters();



        MergeMaterials();

        mixSurf.BuildNeighborTable_nonmanifold();

        //ChangeMaterials();
        if(issmooth){
            ComputeSmallMVnormal();
            mixSurf.Fairing(true,false,false,20,0.5,0.1,&staticV);
            ComputeSmallMVnormal();
            SurfaceFFFairing(0.5,100,&staticV);
        }
        //ComputeSmallMVnormal();



        //mixSurf.ReScale(1.0,0.5,0.5);
        mixSurf.ReScale_uniform(1.0);


        ReAllocation();
    }


    cout<<"number of Materia: "<<numofSurf<<endl;

    return numofSurf;


}



int MultiSurface::InitializeForCutPlane(string filename,string writefilename,double x,double y,double z){
    clearup();
    bool rere = false;

    //rere = readSufFile(filename,true);
    rere = ReadFile(filename);
    if(!rere)return -1;
    {
        mixSurf.setparameters();

        mixSurf.BuildNeighborTable_nonmanifold();



        mixSurf.ReScale(x,y,z);

        ReAllocation();
    }

    writeObjFile(writefilename,mixSurf.vertices,mixSurf.faces2vertices);

    cout<<"number of Materia: "<<numofSurf<<endl;

    return numofSurf;


}


void MultiSurface::ImportFromSuf(SufStructure &pSuf,bool isAutoSmooth,bool isRescale, bool isReAllocation){



    cout<<"ImportFromSuf"<<endl;

    clearup();


    mixSurf.vertices=pSuf.Vs;

    mixSurf.faces2vertices=pSuf.Fs;

    mixSurf.setparameters();


    surf_group1.resize(mixSurf.n_faces);
    surf_group2.resize(mixSurf.n_faces);
    for(int i=0;i<mixSurf.n_faces;++i){
        surf_group1[i] = pSuf.FMs[i*2];
        surf_group2[i] = pSuf.FMs[i*2+1];

    }
    curvee2v.clear();
    for(auto a:pSuf.Ctrs)curvee2v.push_back(a);

    n_curveedges = curvee2v.size()/2;


    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;


    set<int>matSet;
    for(auto a:surf_group1)matSet.insert(a);
    for(auto a:surf_group2)matSet.insert(a);

    if(1){
        map<int,int>mappM;
        int newMInd = 0;
        for(auto a:matSet)mappM[a] = newMInd++;
        for(auto &a:surf_group1)a = mappM[a];
        for(auto &a:surf_group2)a = mappM[a];
        numofSurf = matSet.size();

    }else {
        numofSurf = *max_element(matSet.begin(),matSet.end())+1;
    }


    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;


    if(isRescale){
        GetScaleInfoViaCrossSection(this->lscale,this->pcenter);
        mixSurf.ReScale(this->lscale,this->pcenter);
    }
    //cout<<"numofSurf: "<<numofSurf<<endl;

    if(isAutoSmooth){
        //ComputeSmallMVnormal();
        //mixSurf.Fairing(true,false,false,50,0.5,0.1,&staticV);
        mixSurf.BuildNeighborTable_nonmanifold();
        ComputeSmallMVnormal();
        SurfaceFFFairing(0.5,100,&staticV);
        //SurfaceFFFairing(0.5,100,&staticV);
    }

    ComputeSmallMVnormal();

    if(isReAllocation)ReAllocation();


}
bool readVector_uint(string filename,vector<uint>&vvv){
    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the file " << filename << endl;
        return false;
    }
    int nnum = 0;
    reader>>nnum;
    nnum*=2;
    vvv.resize(nnum);
    for(int i=0;i<nnum;++i)reader>>vvv[i];
    reader.close();
    return true;
}
int MultiSurface::ReadBatchResult(string inFpre, string insurend, string incontourend, int n_m){

    numofSurf = n_m;

    surf_group1.clear();
    surf_group2.clear();
    mixSurf.clearup();
    curvee2v.clear();
    for(int i=1;i<n_m;++i){
        vector<uint>vvv;
        readVector_uint(inFpre + to_string(i) + incontourend,vvv);
        for(auto &a:vvv)a+=mixSurf.n_vertices;
        for(auto a:vvv)curvee2v.push_back(a);


        Mesh a;
        a.readfile(inFpre + to_string(i) + insurend );
        //a.BuildNeighborTable();
        //a.ReOrientFaces();
        mixSurf.addMesh(a);


        int be = surf_group1.size();
        int nf = a.n_faces;
        surf_group1.resize(be+nf);
        for(int j=be;j<surf_group1.size();++j)surf_group1[j] = i;


    }
    surf_group2.resize(surf_group1.size(),0);
    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;

    mixSurf.ReScale_uniform(1.0);
    mixSurf.BuildNeighborTable_nonmanifold();
    cout<<mixSurf.n_faces<<endl;

    //curvee2v.resize(10);
    //for(int i=0;i<curvee2v.size();++i)curvee2v[i]=i;

    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;
    //for(auto a:curvee2v)cout<<a<<' ';cout<<endl;


    ReAllocation();

    return n_m;

}

int MultiSurface::ReadBatchResult(string inFpre, string insurend, string incontourend, string infaceMat, string inVMat, int n_f,bool ismanifold){

    numofSurf = n_f;

    surf_group1.clear();
    surf_group2.clear();
    mixSurf.clearup();
    curvee2v.clear();

    hidden_group.resize(n_f);
    hidden_groupv.resize(n_f);
    for(int i=0;i<n_f;++i){
        vector<uint>vvv;
        readVector_uint(inFpre + to_string(i) + incontourend,vvv);
        for(auto &a:vvv)a+=mixSurf.n_vertices;
        for(auto a:vvv)curvee2v.push_back(a);


        Mesh a;
        a.readfile(inFpre + to_string(i) + insurend );
        //a.BuildNeighborTable();
        //a.ReOrientFaces();
        mixSurf.addMesh(a);


        int be = surf_group1.size();
        int nf = a.n_faces;
        surf_group1.resize(be+nf);
        for(int j=be;j<surf_group1.size();++j)surf_group1[j] = i;


        readVecFile(inFpre + to_string(i) + infaceMat,hidden_group[i]);
        readVecFile(inFpre + to_string(i) + inVMat,hidden_groupv[i]);



    }
    surf_group2.resize(surf_group1.size(),0);

    mixSurf.ReScale_uniform(1.0);
    mixSurf.BuildNeighborTable_nonmanifold();
    cout<<mixSurf.n_faces<<endl;

    //curvee2v.resize(10);
    //for(int i=0;i<curvee2v.size();++i)curvee2v[i]=i;

    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;
    //for(auto a:curvee2v)cout<<a<<' ';cout<<endl;



    ReAllocation();
    HiddenGroup_SpecialAllocation();




    return n_f;
}

int MultiSurface::ReadBatchResult_specialCellmode(int codeN, string inFpre, string insurend, string inf2Cells, string infaceMat, string inCtrE){



    codeN = 200000;
    mixSurf.readfile(inFpre + to_string(codeN) + insurend );
    mixSurf.ReScale_uniform(1.0);

    mixSurf.BuildUp(false);
    //hiddenF2Mat.clear();hiddenV2Mat.clear();hidden_F2CSs.clear();hidden_F2Cells.clear();hiddencontainer_ctrE.clear();
    vector<int>hiddenF2Mat,verticesMat,faces2Cs,hidden_F2Cells,hiddencontainer_ctrE;
    readVecFile(inFpre + to_string(codeN) + infaceMat,hiddenF2Mat);
    //readVecFile(inFpre + to_string(codeN) + inVMat,hiddenV2Mat);
    //readVecFile(inFpre + to_string(codeN) + inf2Cs,hidden_F2CSs);
    readVecFile(inFpre + to_string(codeN) + inf2Cells,hidden_F2Cells);
    if(!inCtrE.empty())readContourEdgeTxtFile(inFpre + to_string(codeN) + inCtrE,hiddencontainer_ctrE);

    surf_group1 = hiddenF2Mat;surf_group2 = hiddenF2Mat;

    int nMat = *max_element(hiddenF2Mat.begin(),hiddenF2Mat.end())+1;
    colordegree.resize(nMat);
    for(int i=0;i<nMat;++i){
        colordegree[i] = 340*double(i)/double(nMat);
    }
    colordegree[0] = -1;






    numofSurf = *max_element(hidden_F2Cells.begin(),hidden_F2Cells.end())+1;
    vector< vector<uint> >mm_faces(numofSurf);
    vector< vector<uint> >mm_facesmapping(numofSurf);
    vector< vector<uint> >inver_faces_mat(numofSurf);
    for(int i =0;i<mixSurf.n_faces;i++){

        int m1 = hidden_F2Cells[i*2];
        int m2 = hidden_F2Cells[i*2 +1 ];
        auto p_fv = mixSurf.fv_begin(i);


        //if(i>8093)for(int j=0;j<9;++j)cout<<p_fv[j]<<endl;
        if(m1!=-1)for(int j=0;j<3;++j){mm_faces[m1].push_back(p_fv[j]);}



        //if(i>8093)cout<<i<<' '<<m1<<' '<<m2<<endl;
        if(m2!=-1)for(int j=0;j<3;++j){mm_faces[m2].push_back(p_fv[j]);}
        //if(m2!=-1){mm_faces[m2].push_back(p_fv[2]);mm_faces[m2].push_back(p_fv[1]);mm_faces[m2].push_back(p_fv[0]);}

        if(m1!=-1)mm_facesmapping[m1].push_back(i);
        if(m2!=-1)mm_facesmapping[m2].push_back(i);

        if(1){
            if(m1!=-1)inver_faces_mat[m1].push_back(hiddenF2Mat[i]);
            if(m2!=-1)inver_faces_mat[m2].push_back(hiddenF2Mat[i]);
        }

    }

    cout<<"ReadBatchResult_specialCellmode debug 1"<<endl;
    vector< int > newverticeInd (mixSurf.n_vertices,-1);
    inversemap.clear();surfbuffer.clear();Vpos_specialcellmode.clear();CtrE_specialcellmode.clear();
    inversemap.resize(numofSurf);
    surfbuffer.resize(numofSurf);
    Vpos_specialcellmode.resize(numofSurf+1);CtrE_specialcellmode.resize(numofSurf+1);

    int numofPickup = 0;
    vector<uint>&e2v = CtrE_specialcellmode[0];
    vector<double>&vpos = Vpos_specialcellmode[0];

    for(auto &a:newverticeInd)a = -1;
    for(auto a:hiddencontainer_ctrE)newverticeInd[a] = 0;
    for(int j=0;j< newverticeInd.size();++j)if(newverticeInd[j]==0){
        newverticeInd[j] = numofPickup++;
        auto p_v = mixSurf.v_begin(j);
        for(int k=0;k<3;++k)vpos.push_back(p_v[k]);
    }
    for(auto a:hiddencontainer_ctrE)e2v.push_back(newverticeInd[a]);

    for(int i =0;i<numofSurf;i++){

        auto &pinversemap = inversemap[i];
        auto &psurfbuffer = surfbuffer[i];
        auto &pmm_faces = mm_faces[i];
        for(auto &a:newverticeInd)a = -1;
        pinversemap.clear();
        int numofPickup =0;

        for(auto a:pmm_faces)newverticeInd[a] = 0;
        for(auto &a : newverticeInd)if(a==0){
            a = numofPickup;
            ++numofPickup;
        }
        pinversemap.resize(numofPickup);
        for(int j=0;j<newverticeInd.size();++j)if(newverticeInd[j]!=-1){
            pinversemap[newverticeInd[j]] = j;
        }

        for(auto &a:pmm_faces)a = newverticeInd[a];
        psurfbuffer.faces2vertices = pmm_faces;

        psurfbuffer.vertices.resize(numofPickup*3);
        for(int j = 0;j<numofPickup;++j){
            copyVec(mixSurf.v_begin(pinversemap[j]),psurfbuffer.v_begin(j));
        }



        numofPickup = 0;
        vector<uint>&e2v = CtrE_specialcellmode[i+1];
        vector<double>&vpos = Vpos_specialcellmode[i+1];
        for(int j=0,endj = hiddencontainer_ctrE.size()/2;j<endj;++j){
            int a1 = hiddencontainer_ctrE[j*2];
            int a2 = hiddencontainer_ctrE[j*2+1];

            if(newverticeInd[a1]>=0 && newverticeInd[a2]>=0){
                e2v.push_back(a1);
                e2v.push_back(a2);
            }
        }

        for(auto &a:newverticeInd)a = -1;
        for(auto a:e2v)newverticeInd[a] = 0;
        for(int j=0;j< newverticeInd.size();++j)if(newverticeInd[j]!=-1){
            newverticeInd[j] = numofPickup++;
            auto p_v = mixSurf.v_begin(j);
            for(int k=0;k<3;++k)vpos.push_back(p_v[k]);
        }
        for(auto &a:e2v)a = newverticeInd[a];

    }

    for(int i=0;i<numofSurf;++i){
        auto &a = surfbuffer[i];

        a.setparameters();
        a.ReOrientFaces();

        a.colordegree = colordegree[i];
        a.colormethod = 2;

        auto &c_inver_faces_mat = inver_faces_mat[i];
        vector<double>inver_faces_colordegree(c_inver_faces_mat.size());
        for(int j=0;j<c_inver_faces_mat.size();++j)inver_faces_colordegree[j] = colordegree[c_inver_faces_mat[j]];
        a.GetPerFaceColorDegree(inver_faces_colordegree);

        a.BuildUp(false);
        a.reinitflag = true;

        //cout<<a.weighted_fcolor.size()<<endl;
        //a.isbuilddisp = false;
    }

    for(int i =0;i<numofSurf;i++){

        auto &pinversemap = inversemap[i];
        auto &psurfbuffer = surfbuffer[i];
        auto &pfacesmapping = mm_facesmapping[i];


        for(int j=0;j<pfacesmapping.size();++j){
            auto p_fvd = mixSurf.fv_begin(pfacesmapping[j]);
            auto p_fvs = psurfbuffer.fv_begin(j);
            for(int k=0;k<3;++k)p_fvd[k] = pinversemap[p_fvs[k]];
        }
    }


    ComputeSmallMVnormal();
    mixSurf.reinitflag = true;
    mixSurf.colordegree = colordegree[0];
    mixSurf.colormethod = 3;
    vector<double>faces_colordegree(surf_group1.size());
    //for(int j=0;j<surf_group1.size();++j)faces_colordegree[j] = colordegree[min(surf_group1[j],surf_group2[j])==0?max(surf_group1[j],surf_group2[j]):min(surf_group1[j],surf_group2[j])];
    for(int j=0;j<surf_group1.size();++j)faces_colordegree[j] = colordegree[hiddenF2Mat[j]];
    mixSurf.GetPerFaceColorDegree(faces_colordegree);
    for(int j=0;j<surf_group2.size();++j)faces_colordegree[j] = colordegree[hiddenF2Mat[j]];
    mixSurf.GetPerFaceColorDegreeBackside(faces_colordegree);


    return numofSurf;




}






int MultiSurface::GetNMat(){return numofSurf;}



void MultiSurface::SmoothSurfNet_Laplacian(int iter){
    mixSurf.SurfaceLaplacianFairing(0.5,iter,&staticV);

    ReAllocation();
}
void MultiSurface::SmoothSurfNet_JuFair(int iter){

    SurfaceFFFairing(0.5,iter,&staticV);

    ReAllocation();
}

void MultiSurface::SurfaceFFFairing(double lambda,int iter,vector<bool>* staticv){

    vector<bool>isInverseNor(mixSurf.n_faces);
    for(int i=0;i<mixSurf.n_faces;++i){
        if(surf_group1[i]>surf_group2[i])isInverseNor[i] = false;
        else isInverseNor[i] = true;
    }

    mixSurf.SurfaceFFFairing(lambda,iter,&isInverseNor,&staticV);


}

void MultiSurface::LiepaRefinement(double alpha){

    cout<<"LiepaRefinement Disable"<<endl;
//    LuMesh m_meshFair;
//    vector<float> m_verLu(mixSurf.vertices.size());
//    intvector m_faceLu(mixSurf.faces2vertices.size());
//    intvector m_faceMat(surf_group1.size()*2);
//    intvector m_ctrmedgeLu(curvee2v.size());
//    bool m_doSwap = true;

//    auto tmpcurv = curvee2v;

//    auto &mverOri = mixSurf.vertices;
//    auto &mfaceOri = mixSurf.faces2vertices;
//    for(int i=0;i<m_verLu.size();++i)m_verLu[i] = mverOri[i];
//    for(int i=0;i<m_faceLu.size();++i)m_faceLu[i] = mfaceOri[i];
//    for(int i=0;i<m_ctrmedgeLu.size();++i)m_ctrmedgeLu[i] = curvee2v[i];
//    for(int i=0;i<surf_group1.size();++i){m_faceMat[i*2] = surf_group1[i];m_faceMat[i*2+1] = surf_group2[i];}

//    //m_ctrmedgeLu.clear();
//    m_meshFair.InputData(m_verLu,m_faceLu,m_faceMat,m_ctrmedgeLu, m_doSwap);

//    cout<<"Liepa refineMent Start!"<<endl;

//    m_meshFair.LiepaRefine(alpha);

//    cout<<"Liepa refineMent finish!"<<endl;


//    m_meshFair.OutputData(m_verLu,m_faceLu,m_faceMat);

//    clearup();

//    mverOri.resize(m_verLu.size());
//    mfaceOri.resize(m_faceLu.size());
//    for(int i=0;i<m_verLu.size();++i)mverOri[i] = m_verLu[i];
//    for(int i=0;i<m_faceLu.size();++i)mfaceOri[i] = m_faceLu[i];

//    surf_group1.resize(m_faceMat.size()/2);surf_group2.resize(m_faceMat.size()/2);
//    for(int i=0;i<surf_group1.size();++i){ surf_group1[i] = m_faceMat[i*2];surf_group2[i] = m_faceMat[i*2+1];}






//    //writeObjFile(string("../liepatest"),mverOri,mfaceOri);
//    mixSurf.face_Material1 = surf_group1;
//    mixSurf.face_Material2 = surf_group2;

//    mixSurf.setparameters();
//    mixSurf.BuildNeighborTable_nonmanifold();

//    cout<<mixSurf.n_faces<<' '<<surf_group1.size()<<endl;

//    curvee2v = tmpcurv;
//    staticV.clear();
//    staticV.resize(mixSurf.n_vertices,false);
//    for(auto a:curvee2v)staticV[a] = true;



//    ReAllocation();




}
bool MultiSurface::clearup(){

    mixSurf.clearup();

    surf_group1.clear();
    surf_group2.clear();


    surfbuffer.clear();
    inversemap.clear();


    curvee2v.clear();
    curve_inversemap.clear();

    crossSection.reset();
    NonmanifoldNet.reset();

    staticV.clear();

    Vpos_specialcellmode.clear();
    CtrE_specialcellmode.clear();




    display_vertices.clear();
    display_normal.clear();
    display_edges.clear();


    display_field_dot.clear();
    display_vcolor.clear();
    display_faces.clear();



}


bool MultiSurface::ReAllocation(bool ismanifold, int colormethod, char *inputcolors,bool isImportContour){


    cout<<"ReAllocation"<<endl;
    vector< int > newverticeInd (mixSurf.n_vertices,-1);
    for(auto a:curvee2v)newverticeInd[a] = 0;

    int numofPickup = 0;
    for(auto &a : newverticeInd)if(a==0){
        a = numofPickup;
        ++numofPickup;
    }
    curve_inversemap.resize(numofPickup);
    for(int j=0;j<newverticeInd.size();++j)if(newverticeInd[j]!=-1){
        curve_inversemap[newverticeInd[j]] = j;
    }

    vector<uint>e2v;
    for(auto &a:curvee2v)e2v.push_back(newverticeInd[a]);

    vector<double>veveve(numofPickup*3);
    for(int j = 0;j<numofPickup;++j){
        copyVec(mixSurf.v_begin(curve_inversemap[j]),veveve.data()+3*j);
    }

    if(isImportContour)crossSection.ImportCurve(veveve,e2v);

    int nonvind = 0;
    veveve.clear();
    vector<int>nonv(mixSurf.n_vertices,-1);
    vector<uint>nonedges;
    if(0){
        for(int i=0;i<mixSurf.n_vertices;++i)if(mixSurf.vertices_non_manifold[i]){
            nonv[i] = nonvind++;
            auto p_v = mixSurf.v_begin(i);
            for(int j=0;j<3;++j)veveve.push_back(p_v[j]);
        }
        for(int i=0;i<mixSurf.n_edges;++i)if(mixSurf.edge_non_manifold[i]){
            auto p_ev = mixSurf.ev_begin(i);
            nonedges.push_back(nonv[p_ev[0]]);nonedges.push_back(nonv[p_ev[1]]);
        }

        NonmanifoldNet.ImportCurve(veveve,nonedges);
    }




    //cout<<"ReAllocation "<<mixSurf.n_faces<<endl;
    //for(auto a: mixSurf.faces2vertices)cout<<a<<' ';cout<<endl;
    //cout<<mixSurf.faces2vertices.size()/3-mixSurf.n_faces<<endl;

    vector< vector<uint> >mm_faces(numofSurf);
    vector< vector<uint> >inver_faces_mat(numofSurf);
    for(int i =0;i<mixSurf.n_faces;i++){

        int m1 = surf_group1[i];
        int m2 = surf_group2[i];
        auto p_fv = mixSurf.fv_begin(i);


        //if(i>8093)for(int j=0;j<9;++j)cout<<p_fv[j]<<endl;
        for(int j=0;j<3;++j){mm_faces[m1].push_back(p_fv[j]);}



        //if(i>8093)cout<<i<<' '<<m1<<' '<<m2<<endl;
        //for(int j=0;j<3;++j){mm_faces[m2].push_back(p_fv[j]);}
        mm_faces[m2].push_back(p_fv[2]);mm_faces[m2].push_back(p_fv[1]);mm_faces[m2].push_back(p_fv[0]);

        if(1){
            inver_faces_mat[m1].push_back(m2);
            inver_faces_mat[m2].push_back(m1);
        }else{
            inver_faces_mat[m1].push_back(m1==0?m2:(m2==0?m1:m2));
            inver_faces_mat[m2].push_back(m2==0?m1:(m1==0?m2:m1));
        }

    }



    cout<<"reallocation debug 1"<<endl;
    inversemap.resize(numofSurf);
    surfbuffer.resize(numofSurf);
    for(int i =0;i<numofSurf;i++){

        auto &pinversemap = inversemap[i];
        auto &psurfbuffer = surfbuffer[i];
        auto &pmm_faces = mm_faces[i];
        for(auto &a:newverticeInd)a = -1;
        pinversemap.clear();
        int numofPickup =0;

        for(auto a:pmm_faces)newverticeInd[a] = 0;
        for(auto &a : newverticeInd)if(a==0){
            a = numofPickup;
            ++numofPickup;
        }
        pinversemap.resize(numofPickup);
        for(int j=0;j<newverticeInd.size();++j)if(newverticeInd[j]!=-1){
            pinversemap[newverticeInd[j]] = j;
        }

        for(auto &a:pmm_faces)a = newverticeInd[a];
        psurfbuffer.faces2vertices = pmm_faces;

        psurfbuffer.vertices.resize(numofPickup*3);
        for(int j = 0;j<numofPickup;++j){
            copyVec(mixSurf.v_begin(pinversemap[j]),psurfbuffer.v_begin(j));
        }
    }

    colordegree.resize(numofSurf);
    //colormethod = 0;
    if(colormethod ==0){
        for(int i=0;i<numofSurf;++i){
            colordegree[i] = 340*double(i)/double(numofSurf);
        }
    }if(colormethod ==1){
        for(int i=0;i<numofSurf;++i){
            colordegree[i] = 340*double(i)/double(numofSurf-1);
        }
    }
    colordegree[0] = -1;
   // colordegree[1] = 174;


    for(int i=0;i<numofSurf;++i){
        auto &a = surfbuffer[i];




        a.setparameters();

        a.colordegree = colordegree[i];
        a.colormethod = 0;

        auto &c_inver_faces_mat = inver_faces_mat[i];
        vector<double>inver_faces_colordegree(c_inver_faces_mat.size());
        for(int j=0;j<c_inver_faces_mat.size();++j)inver_faces_colordegree[j] = colordegree[c_inver_faces_mat[j]];
        a.GetPerFaceColorDegree(inver_faces_colordegree);

        if(1){a.BuildUp(false);}
        else{
            a.ComputeFaceNormal(false);
            //a.GetRescaleInfo(lscale, pcenter);
        }
        a.reinitflag = true;

        //cout<<a.weighted_fcolor.size()<<endl;
        //a.isbuilddisp = false;
    }
    ComputeSmallMVnormal();
    mixSurf.reinitflag = true;
    mixSurf.colordegree = colordegree[0];
    mixSurf.colormethod = 3;
    vector<double>faces_colordegree(surf_group1.size());
    //for(int j=0;j<surf_group1.size();++j)faces_colordegree[j] = colordegree[min(surf_group1[j],surf_group2[j])==0?max(surf_group1[j],surf_group2[j]):min(surf_group1[j],surf_group2[j])];
    for(int j=0;j<surf_group1.size();++j)faces_colordegree[j] = colordegree[min(surf_group1[j],surf_group2[j])];
    mixSurf.GetPerFaceColorDegree(faces_colordegree);
    //mixSurf.GetPerFaceColorDegreeBackside(faces_colordegree);


    for(int j=0;j<surf_group2.size();++j)faces_colordegree[j] = colordegree[max(surf_group1[j],surf_group2[j])];
    mixSurf.GetPerFaceColorDegreeBackside(faces_colordegree);
    //mixSurf.GetPerFaceColorDegree(faces_colordegree);




    //for(int i=0;i<2;++i)surfbuffer[2].ReOrientFaces();
    //surfbuffer[2].SaveInterface(string("ads"));

    //    for(int i=0;i<numofSurf;++i){
    //        cout<<surfbuffer[i].weighted_fcolor.size()<<endl;
    //    }
    cout<<"reallocation debug end "<<surfbuffer.size()<<endl;

}


bool MultiSurface::HiddenGroup_SpecialAllocation(){


    cout<<"HiddenGroup_SpecialAllocation"<<endl;

    vector<uint>maxM(numofSurf,0);
    for(int i=1;i<numofSurf;++i)maxM[i] = *max_element(hidden_group[i].begin(),hidden_group[i].end());
    int n_m = *max_element(maxM.begin(),maxM.end());
    vector<double>colordegree(n_m);
    for(int i=0;i<n_m;++i){
        colordegree[i] = 360*double(i)/double(n_m+1)+90;
    }
    for(int i=1;i<numofSurf;++i){
        auto &a = surfbuffer[i];

        a.setparameters();

        a.colordegree = 50;
        a.reinitflag = true;


        vector<double>inver_faces_colordegree(a.n_faces);
        for(int j=0;j<inver_faces_colordegree.size();++j)inver_faces_colordegree[j] = colordegree[hidden_group[i][j]];
        a.GetPerFaceColorDegree(inver_faces_colordegree);

        vector<double>inver_vertices_colordegree(a.n_vertices);
        for(int j=0;j<inver_vertices_colordegree.size();++j)inver_vertices_colordegree[j] = colordegree[hidden_groupv[i][j]];
        a.GetPerVertexColorDegree(inver_vertices_colordegree);

        a.hiddenCtr(hidden_groupv[i]);




        //cout<<a.weighted_fcolor.size()<<endl;



        //a.isbuilddisp = false;
    }

}


bool MultiSurface::ReadFile(string filename){



    string filepath;
    string modname;
    string extname;

    SplitFileName(filename,filepath,modname,extname);
    cout<<modname<<' '<<extname<<' '<<filepath<<endl;

    bool issuc = false;


    //issuc = ReadFile(filename);

    if(extname==".suf"){issuc = readSufFile(filename);}
    else if(extname==".aoff"){issuc = readAoffFile(filename,true);}
    else if(extname==".csl"){issuc = readCslFile(filepath+modname,true);}
    else {
        //        if(extname==".obj" || extname==".off"){
        //            surfbuffer.resize(1);
        //            issuc = surfbuffer[0].ReadFile(filename);

        //        }
        issuc = mixSurf.readfile(filename);
        surf_group1.clear();
        surf_group2.clear();
        surf_group1.resize(mixSurf.n_faces,1);
        surf_group2.resize(mixSurf.n_faces,0);
        mixSurf.face_Material1 = surf_group1;
        mixSurf.face_Material2 = surf_group2;

        numofSurf=2;
    }

    if(issuc){
        modelname = modname;
        prepath = filepath;
        ext=extname;
    }

    return issuc;
}


bool MultiSurface::readSufFile(string filename,bool isreArrMat){

    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Suf file " << filename << endl;
        return false;
    }
    mixSurf.clearup();

    cout<<"Reading Suf File"<<endl;



    reader>>mixSurf.n_vertices;
    reader>>mixSurf.n_faces;

    mixSurf.vertices.resize(mixSurf.n_vertices*3);

    mixSurf.faces2vertices.clear();

    surf_group1.clear();
    surf_group2.clear();

    for(int i =0;i<mixSurf.vertices.size();i++){
        reader>>mixSurf.vertices[i];
    }

    int ivalue;
    for(int i =0;i<mixSurf.n_faces;i++){
        for(int j =0;j<3;++j){reader>>ivalue;mixSurf.faces2vertices.push_back(ivalue);}
        reader>>ivalue;surf_group1.push_back(ivalue);
        reader>>ivalue;surf_group2.push_back(ivalue);
    }

    reader>>n_curveedges;
    curvee2v.resize(n_curveedges*2);
    for(int i =0;i<curvee2v.size();i++){
        reader>>curvee2v[i];
    }

    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;


    set<int>matSet;
    for(auto a:surf_group1)matSet.insert(a);
    for(auto a:surf_group2)matSet.insert(a);

    if(isreArrMat){
        map<int,int>mappM;
        int newMInd = 0;
        for(auto a:matSet)mappM[a] = newMInd++;
        for(auto &a:surf_group1)a = mappM[a];
        for(auto &a:surf_group2)a = mappM[a];
        numofSurf = matSet.size();

    }else {
        numofSurf = *max_element(matSet.begin(),matSet.end())+1;
    }


    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;












    cout<<"numofSurf: "<<numofSurf<<endl;
    return true;






}

bool MultiSurface::WriteSuf(string filename){

    vector<int >fmat;
    for(int i=0;i<mixSurf.n_faces;++i){
        fmat.push_back(surf_group1[i]);fmat.push_back(surf_group2[i]);
    }

    vector<int>ctre;
    for(auto a:curvee2v)ctre.push_back(a);

    return writeSufFile(filename,mixSurf.vertices,mixSurf.faces2vertices,fmat,ctre);




}

bool MultiSurface::WriteObj(string filename){




    vector<uint>faces2vertices;
    for(int i=0;i<mixSurf.n_faces;++i){
        auto p_fv = mixSurf.fv_begin(i);
        if(surf_group1[i]<surf_group2[i]){
            for(int j=0;j<3;++j)faces2vertices.push_back(p_fv[j]);

        }else{
            for(int j=2;j>=0;--j)faces2vertices.push_back(p_fv[j]);

        }



    }

    return writeObjFile(filename,mixSurf.vertices,faces2vertices);

}

bool MultiSurface::WritePickMatObj(string filename,int pickmat){



//    vector<double>newvertices;
//    vector<uint>faces2vertices;
//    vector<int>pickvertices(mixSurf.n_vertices,-1);
//    for(int i=0;i<mixSurf.n_faces;++i){

//        if(surf_group1[i]==pickmat || surf_group2[i] == pickmat){
//            auto p_fv = mixSurf.fv_begin(i);
//            for(int j=0;j<3;++j)pickvertices[p_fv[j]] = 0;
//            if(surf_group1[i]<surf_group2[i]){
//                for(int j=0;j<3;++j)faces2vertices.push_back(p_fv[j]);
//            }else{
//                for(int j=2;j>=0;--j)faces2vertices.push_back(p_fv[j]);
//            }
//        }

//    }
//    for(int i=0;i<mixSurf.n_vertices;++i){


//    }

    Surface &picksurf = surfbuffer[pickmat];

    return writeObjFile(filename,picksurf.vertices,picksurf.faces2vertices);

}

bool MultiSurface::readAoffFile(string filename,bool isreArrMat){

    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Aoff file " << filename << endl;
        return false;
    }
    mixSurf.clearup();

    cout<<"Reading Aoff File"<<endl;


    string titl;
    int tmp;


    reader>>titl;

    reader>>mixSurf.n_vertices;
    reader>>mixSurf.n_faces;

    reader>>tmp;
    mixSurf.vertices.resize(mixSurf.n_vertices*3);

    mixSurf.faces2vertices.clear();

    surf_group1.clear();
    surf_group2.clear();

    for(int i =0;i<mixSurf.vertices.size();i++){
        reader>>mixSurf.vertices[i];
    }

    int ivalue;
    for(int i =0;i<mixSurf.n_faces;i++){
        reader>>tmp;
        for(int j =0;j<3;++j){reader>>ivalue;mixSurf.faces2vertices.push_back(ivalue);}

    }
    reader>>titl;
    for(int i =0;i<mixSurf.n_faces;i++){
        reader>>ivalue;surf_group1.push_back(ivalue);
        reader>>ivalue;surf_group2.push_back(ivalue);
    }

    curvee2v.clear();


    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;


    set<int>matSet;
    for(auto a:surf_group1)matSet.insert(a);
    for(auto a:surf_group2)matSet.insert(a);

    if(isreArrMat){
        map<int,int>mappM;
        int newMInd = 0;
        for(auto a:matSet)mappM[a] = newMInd++;
        for(auto &a:surf_group1)a = mappM[a];
        for(auto &a:surf_group2)a = mappM[a];
        numofSurf = matSet.size();

    }else {
        numofSurf = *max_element(matSet.begin(),matSet.end())+1;
    }


    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;










    cout<<"numofSurf: "<<numofSurf<<endl;
    return true;






}



bool MultiSurface::readCslFile(string prefixname,bool isreArrMat){


    string cslfilename = prefixname + string(".csl");
    ifstream reader(cslfilename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Csl file " << cslfilename << endl;
        return false;
    }
    mixSurf.clearup();

    cout<<"Reading Aoff File"<<endl;


    string header;
    reader>>header;
    int nP, nMat;
    reader>>nP; reader>>nMat;
    vector<vector<double>>planeinfo(nP);
    vector<int>plane2nV(nP);
    vector<int>plane2nE(nP);
    vector<vector<double>>plane2vertices(nP);
    vector<vector<unsigned int>>plane2edges(nP);
    vector<vector<int>>plane2edgesMat(nP);


    for(int i=0;i<nP;++i){
        int val;
        int ncomp;
        reader>>val;
        reader>>plane2nV[i];
        reader>>ncomp;
        auto &p_planeinfo = planeinfo[i];
        p_planeinfo.resize(4);
        for(int j=0;j<4;++j){
            reader>>p_planeinfo[j];
        }
        //        reader>>plane2nV[i];
        //        reader>>plane2nE[i];
        auto &p_plane2vertices = plane2vertices[i];
        p_plane2vertices.resize(plane2nV[i]*3);
        for(int nv = plane2nV[i]*3,j=0;j<nv;++j){
            reader>>p_plane2vertices[j];
        }
        reader>>plane2nE[i];
        reader>>val;
        auto &p_plane2edges = plane2edges[i];
        auto &p_plane2edgesMat = plane2edgesMat[i];
        //        p_plane2edges.resize(plane2nE[i]*2);
        //        p_plane2edgesMat.resize(plane2nE[i]*2);
        int val1,val2;
        reader>>val1;
        for(int j=1;j<plane2nE[i];++j){
            reader>>val2;
            p_plane2edges.push_back(val1);p_plane2edges.push_back(val2);
            p_plane2edgesMat.push_back(0);p_plane2edgesMat.push_back(1);

            val1 = val2;
        }

        p_plane2edges.push_back(val1);p_plane2edges.push_back(p_plane2edges[0]);
        p_plane2edgesMat.push_back(0);p_plane2edgesMat.push_back(1);

    }



    /*******************************************************/
    /*******************************************************/


    CrossSections CS;
    CS.ReadCrossSections(planeinfo,plane2vertices,plane2edges,plane2edgesMat);
    vector<double>vvvv;
    vector<uint>e2v;
    CS.Stackup(vvvv,e2v);
    crossSection.ImportCurve(vvvv,e2v);




    mixSurf.readOfffile(prefixname + string(".off"));


    curvee2v.clear();
    staticV.clear();
    staticV.resize(mixSurf.n_vertices,false);
    for(auto a:curvee2v)staticV[a] = true;


    surf_group1.clear();
    surf_group2.clear();
    surf_group1.resize(mixSurf.n_faces,1);
    surf_group2.resize(mixSurf.n_faces,0);
    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;

    numofSurf=2;

    cout<<"numofSurf: "<<numofSurf<<endl;
    return true;

}
void MultiSurface::GetScaleInfo(double &lscale, double *pcenter){

    lscale = this->lscale;
    for(int i=0;i<3;++i)pcenter[i] = this->pcenter[i];

}

void MultiSurface::GetScaleInfoViaCrossSection(double &lscale,double *pcenter){
    for(int i=0;i<3;++i)pcenter[i] = 0;
    double xmin = 999, ymin=9999,zmin=999,xmax=-999,ymax=-999,zmax=-999;
    int nSv = 0;
    for(int i=0;i< mixSurf.n_vertices;++i)if(staticV[i]){
        auto point = mixSurf.v_begin(i);
        xmin = min(xmin,point[0]);xmax = max(xmax,point[0]);
        ymin = min(ymin,point[1]);ymax = max(ymax,point[1]);
        zmin = min(zmin,point[2]);zmax = max(zmax,point[2]);

        for(int j=0;j<3;++j)pcenter[j] += point[j];
        ++nSv;
    }
    for(int j=0;j<3;++j)pcenter[j] /=nSv;
    double largestdis = max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2;
    lscale = 1/largestdis;

    cout<<pcenter[0]<<"  "<<pcenter[1]<<"  "<<pcenter[2]<<"   "<<largestdis<<endl;

}

void MultiSurface::GetScaleInfoViaCrossSection_csl(double &lscale,double *pcenter){
    for(int i=0;i<3;++i)pcenter[i] = 0;
    double xmin = 999, ymin=9999,zmin=999,xmax=-999,ymax=-999,zmax=-999;
    int nSv = 0;
    for(int i=0;i< crossSection.n_vertices;++i){
        auto point = crossSection.v_begin(i);
        xmin = min(xmin,point[0]);xmax = max(xmax,point[0]);
        ymin = min(ymin,point[1]);ymax = max(ymax,point[1]);
        zmin = min(zmin,point[2]);zmax = max(zmax,point[2]);

        for(int j=0;j<3;++j)pcenter[j] += point[j];
        ++nSv;
    }
    for(int j=0;j<3;++j)pcenter[j] /=nSv;
    double largestdis = max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2;
    lscale = 1/largestdis;

    cout<<pcenter[0]<<"  "<<pcenter[1]<<"  "<<pcenter[2]<<"   "<<largestdis<<endl;

}



void MultiSurface::GetColorDegree(vector<double>&out_colordegree){


    out_colordegree = colordegree;

}



void MultiSurface::WriteCutCrossSection(string filename){

    if(Cmode!=S_SURFCUT)return;

    cout<<"write!!!!"<<endl;

    //cutCs.WriteCrossSectionsToYaron(string("../rootinfo/"));

    //cutCs.ManipulateCrossSections();
    //cutCs.WriteCrossSections("../ConvertFolder/root_");

    //cutCs.WriteCrossSections("../ConvertFolder/silicium_");

    //cutCs.WriteCrossSectionsToFCM("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/cutmode/tctrloop.txt",numofSurf);

    //cutCs.WriteCrossSectionsToFCMtCtr("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/cutmode/tctr.txt");

    //cutCs.WriteCrossSectionsToFCM("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/mousebraincutforsurface/tctrloop.txt",numofSurf);

    //cutCs.WriteCrossSectionsToFCMtCtr("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/mousebraincutforsurface/tctr.txt");

    //cutCs.WriteCrossSectionsToFCM("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/newliver/tctrloop.txt",numofSurf);

    //cutCs.WriteCrossSectionsToFCMtCtr("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/newliver/tctr.txt");

    //cutCs.WriteCrossSectionToObj(filename);


    //cutCs.WriteCrossSections("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/mousebraincut/mousebraincut");
    cutCs.WriteCrossSections("/Users/Research/Geometry/FCM/Consistent Contour Networks/Prototyping/Data/newliver/newliver");

    //cutCs.WriteCrossSectionToObj("/Users/Research/Geometry/RBF/data/DataSet/surfacecut/torus.obj");
}
void MultiSurface::DelLastContour(){
    if(Cmode!=S_SURFCUT)return;
    auto &CSs = cutCs.CSs;
    if(CSs.size()>0){
        CSs.resize(CSs.size()-1);
        vector<double>outCtrV;vector<uint>outCtrE;vector<int>outCtrEMat;
        cutCs.Stackup(outCtrV,outCtrE);



        cout<<"CSs.size(): "<<CSs.size()<<endl;
        //for(int i=0;i<4;++i)cout<<para[i]<<' ';cout<<endl;
        crossSection.ImportCurve(outCtrV,outCtrE);

        crossSection.BuildDisplay(false,true);
    }


}
void MultiSurface::DelAllContour(){
    if(Cmode!=S_SURFCUT)return;
    auto &CSs = cutCs.CSs;
    if(CSs.size()>0){
        CSs.resize(0);
        vector<double>outCtrV;vector<uint>outCtrE;vector<int>outCtrEMat;
        cutCs.Stackup(outCtrV,outCtrE);
        cout<<"CSs.size(): "<<CSs.size()<<endl;
        //for(int i=0;i<4;++i)cout<<para[i]<<' ';cout<<endl;
        crossSection.ImportCurve(outCtrV,outCtrE);

        crossSection.BuildDisplay(false,true);
    }


}
void MultiSurface::CutMixSurfaceByPlane(double *para,int sampleStep, int deleteIso){

    if(Cmode!=S_SURFCUT)return;
    auto &CSs = cutCs.CSs;
    CSs.resize(CSs.size()+1);
    Contour &ccCs = CSs[CSs.size()-1];
    vector<double>outCtrV;vector<uint>outCtrE;vector<int>outCtrEMat;
    mixSurf.CutSurfaceByPlane(para,outCtrV,outCtrE,outCtrEMat);

    if(sampleStep>1){
        int n_edges = outCtrE.size()/2;

        vector<int>edgesMat;
        vector<int>edges2Cell(n_edges*2,0);
        for(int i=0;i<n_edges;++i)edges2Cell[i*2]=1;
        for(auto a:outCtrEMat)edgesMat.push_back(a);
        Curve simCurve;
        simCurve.ImportCurve(outCtrV,outCtrE);
        simCurve.CurveNetAnalayze(edgesMat,edges2Cell,2);
        simCurve.CurveNetSimplification(sampleStep,outCtrV,outCtrE,edgesMat,deleteIso);
        outCtrEMat.clear();
        for(auto a:edgesMat)outCtrEMat.push_back(a);
    }

    ccCs.ImportContour(para,outCtrV,outCtrE,outCtrEMat);

    ccCs.CurveSmoothing(200);
    ccCs.ValidationCheck();

    //crossSection.ImportCurve(ccCs.vertices,ccCs.edges);


    cutCs.Stackup(outCtrV,outCtrE);


    cout<<"CSs.size(): "<<CSs.size()<<endl;
    for(int i=0;i<4;++i)cout<<para[i]<<' ';cout<<endl;
    crossSection.ImportCurve(outCtrV,outCtrE);

    crossSection.BuildDisplay(false,true);
}

void MultiSurface::CutMixSurfaceByPlaneTranslation(double *p_trl,double *para,int sampleStep, int deleteIso){

    if(Cmode!=S_SURFCUT)return;
    auto &CSs = cutCs.CSs;
    CSs.resize(CSs.size()+1);
    Contour &ccCs = CSs[CSs.size()-1];
    vector<double>outCtrV;vector<uint>outCtrE;vector<int>outCtrEMat;
    mixSurf.CutSurfaceByPlaneTranslation(p_trl,para,outCtrV,outCtrE,outCtrEMat);

    if(sampleStep>1){
        int n_edges = outCtrE.size()/2;

        vector<int>edgesMat;
        vector<int>edges2Cell(n_edges*2,0);
        for(int i=0;i<n_edges;++i)edges2Cell[i*2]=1;
        for(auto a:outCtrEMat)edgesMat.push_back(a);
        Curve simCurve;
        simCurve.ImportCurve(outCtrV,outCtrE);
        simCurve.CurveNetAnalayze(edgesMat,edges2Cell,2);
        simCurve.CurveNetSimplification(sampleStep,outCtrV,outCtrE,edgesMat,deleteIso);
        outCtrEMat.clear();
        for(auto a:edgesMat)outCtrEMat.push_back(a);
    }

    ccCs.ImportContour(para,outCtrV,outCtrE,outCtrEMat);

    ccCs.CurveSmoothing(500);
    ccCs.ValidationCheck();

    //crossSection.ImportCurve(ccCs.vertices,ccCs.edges);


    cutCs.Stackup(outCtrV,outCtrE);


    cout<<"CSs.size(): "<<CSs.size()<<endl;
    for(int i=0;i<4;++i)cout<<para[i]<<' ';cout<<endl;
    crossSection.ImportCurve(outCtrV,outCtrE);

    crossSection.BuildDisplay(false,true);
}


void MultiSurface::CutMixSurfaceByBatchPlane(vector<vector<double>>&paras){

    if(Cmode!=S_SURFCUT)return;

    vector<double>outCtrV;vector<vector<int>>outCtrE;vector<vector<int>>outCtrEMat;
    vector<double>outCtrVnor;
    mixSurf.CutSurfaceByBatchPlane(paras,outCtrV,outCtrE,outCtrEMat,outCtrVnor);


    vector<uint>out_CtrE;
    vector<int>out_CtrEMat;
    for(auto &a:outCtrE)for(auto b:a)out_CtrE.push_back(b);
    for(auto &a:outCtrEMat)for(auto b:a)out_CtrEMat.push_back(b);
    vector<double>rmf;
    if(true){
        int n_edges = out_CtrE.size()/2;

        vector<int>edgesMat;
        vector<int>edges2Cell(n_edges*2,0);
        for(int i=0;i<n_edges;++i)edges2Cell[i*2]=1;
        for(auto a:out_CtrEMat)edgesMat.push_back(a);
        Curve simCurve;
        simCurve.ImportCurve(outCtrV,out_CtrE);
        simCurve.CurveNetAnalayze(edgesMat,edges2Cell,2);
        simCurve.CurveNormalParalellTransport(outCtrVnor,rmf);
        simCurve.CurveNetSimplification(1,outCtrV,out_CtrE,edgesMat,0);
        //outCtrEMat.clear();
        //for(auto a:edgesMat)out_CtrEMat.push_back(a);
    }


    string filenamessss("../output/ferretCut/ferret_6cut");
    writeCurNetFile(filenamessss + string("_exact"),outCtrV,outCtrE,outCtrEMat,paras,outCtrVnor);
    writeCurNetFile(filenamessss + string("_rmf"),outCtrV,outCtrE,outCtrEMat,paras,rmf);



    crossSection.ImportCurve(outCtrV,out_CtrE);

    crossSection.disscale = 0.3;
    //crossSection.tangent = outCtrVnor;
    //crossSection.tangent.clear();
    crossSection.BuildDisplay(false,true,true);
    cout<<"batch cut with normal"<<endl;
}

vector<bool> MultiSurface::InterpretVisibility(bool isCS, bool isAllMat, bool isNonC, int pickMat){

    vector<bool>re(surfbuffer.size()+3, false);

    if(!isAllMat && (pickMat<0 || pickMat >= surfbuffer.size())){
        cout<<"MultiSurface: pickMat Error, InterpretVisibility failed."<<endl;
        return re;
    }

    int nMat = surfbuffer.size();
    if(isCS)re[0] = true;
    if(isAllMat){

        re[nMat + 2] = true;

    }else{

        re[pickMat+1] = true;

    }

    if(isNonC)re[nMat + 1] = true;
    return re;


}




void MultiSurface::BuildDisplay(infoSurfDisp info, bool rebuildCs){

    //mixSurf.BuildDisplay(info);
    //for(auto &a:surfbuffer)cout<<a.weighted_fcolor.size()<<endl;
    for(auto &a:surfbuffer)a.BuildDisplay(info);
    if(rebuildCs){crossSection.BuildDisplay(false,true,true);}


    mixSurf.BuildDisplay(info);
    //cout<<"aaaaaaaa"<<endl;
    NonmanifoldNet.BuildDisplay(false,true,true);
    display_vertices.resize(surfbuffer.size()+3);
    display_edges.resize(surfbuffer.size()+3);
    display_faces.resize(surfbuffer.size()+3);
    display_normal.resize(surfbuffer.size()+3);
    display_vcolor.resize(surfbuffer.size()+3);

    for(int i=0;i<surfbuffer.size();++i){
        //if(i!=1)continue;
        display_vertices[i+1] = surfbuffer[i].getDisplayVertices();
        display_edges[i+1] = surfbuffer[i].getDisplayEdges();
        display_faces[i+1] = surfbuffer[i].getDisplayFaces();
        display_normal[i+1] = surfbuffer[i].getDisplayVerticesNormal();
        display_vcolor[i+1] = surfbuffer[i].getDisplayColor();
    }
    //    for(int i=0;i<surfbuffer.size();++i){
    //        //if(i!=1)continue;
    //        display_vertices[i+1] = mixSurf.getDisplayVertices();
    //        display_edges[i+1] = mixSurf.getDisplayEdges();
    //        display_faces[i+1] = mixSurf.getDisplayFaces();
    //        display_normal[i+1] = mixSurf.getDisplayVerticesNormal();
    //        display_vcolor[i+1] = mixSurf.getDisplayColor();
    //    }


    //    display_vertices[0] = mixSurf.getDisplayVertices();
    //    display_edges[0] = mixSurf.getDisplayEdges();
    //    display_faces[0] = mixSurf.getDisplayFaces();
    //    display_normal[0] = mixSurf.getDisplayVerticesNormal();
    //    display_vcolor[0] = mixSurf.getDisplayColor();

    display_vertices[0] = crossSection.getDisplayVertices();
    display_edges[0] = crossSection.getDisplayEdges();
    display_faces[0] = crossSection.getDisplayFaces();
    display_normal[0] = crossSection.getDisplayVerticesNormal();
    display_vcolor[0] = crossSection.getDisplayColor();

    int nonind = surfbuffer.size()+1;
    display_vertices[nonind] = NonmanifoldNet.getDisplayVertices();
    display_edges[nonind] = NonmanifoldNet.getDisplayEdges();
    display_faces[nonind] = NonmanifoldNet.getDisplayFaces();
    display_normal[nonind] = NonmanifoldNet.getDisplayVerticesNormal();
    display_vcolor[nonind] = NonmanifoldNet.getDisplayColor();

    int mixind = surfbuffer.size()+2;
    display_vertices[mixind] = mixSurf.getDisplayVertices();
    display_edges[mixind] = mixSurf.getDisplayEdges();
    display_faces[mixind] = mixSurf.getDisplayFaces();
    display_normal[mixind] = mixSurf.getDisplayVerticesNormal();
    display_vcolor[mixind] = mixSurf.getDisplayColor();

    //    display_vertices[mixind] = NonmanifoldNet.getDisplayVertices();
    //    display_edges[mixind] = NonmanifoldNet.getDisplayEdges();
    //    display_faces[mixind] = NonmanifoldNet.getDisplayFaces();
    //    display_normal[mixind] = NonmanifoldNet.getDisplayVerticesNormal();
    //    display_vcolor[mixind] = NonmanifoldNet.getDisplayColor();

    //    mixSurf.getDisplayVertices();
    //    mixSurf.getDisplayEdges();
    //    mixSurf.getDisplayFaces();
    //    mixSurf.getDisplayVerticesNormal();
    //    mixSurf.getDisplayColor();

}


void MultiSurface::BuildDisplay_specialCellmode(int celli){


    //cout<<celli<<endl;

    crossSection.ImportCurve(Vpos_specialcellmode[celli+1],CtrE_specialcellmode[celli+1]);

    crossSection.BuildDisplay(false,true);
    display_vertices[0] = crossSection.getDisplayVertices();
    display_edges[0] = crossSection.getDisplayEdges();
    display_faces[0] = crossSection.getDisplayFaces();
    display_normal[0] = crossSection.getDisplayVerticesNormal();
    display_vcolor[0] = crossSection.getDisplayColor();

}

void MultiSurface::ComputeSmallMVnormal(){

    mixSurf.ComputeFaceNormal(false);
    vector<bool>pisinverN(mixSurf.n_faces);
    for(int i=0;i<mixSurf.n_faces;++i)
        if(surf_group1[i]>surf_group2[i])pisinverN[i]=true;
        else pisinverN[i]=false;
    mixSurf.SetDisplayNormal(pisinverN);
    //mixSurf.ComputeFaceNormal(true,&pisinverN);
    //mixSurf.ComputeFaceNormal(true,NULL);



}
int MultiSurface::ReadandStackContour(string infilename, vector<double>&points, vector<uint> &edges){
    ifstream reader(infilename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Contour file " << infilename << endl;
        return -1;
    }
    int n_c;
    reader>>n_c;

    vector<Contour>mainC(n_c);

    for(int i=0;i<n_c;++i)mainC[i].ReadContour(reader);

    int n_materials = mainC[0].n_materials;




    CrossSections CCC(mainC);
    CCC.Stackup(points,edges);
    return 0;

}

void MultiSurface::MergeMaterials(){



    if(0){
        vector<set<int>>contactness(numofSurf);
        for(int i=0;i<mixSurf.n_faces;++i){
            contactness[surf_group1[i]].insert(surf_group2[i]);
            contactness[surf_group2[i]].insert(surf_group1[i]);
        }

        for(int i=0;i<mixSurf.n_faces;++i){
            cout<<i<<": ";
            for(auto a:contactness[i])cout<<a<<' ';cout<<endl;
        }

        exit(11);
    }
    vector<int>megermap(numofSurf,0);

    if(0){
        for(int i=0;i<megermap.size();++i)megermap[i] = i;



        megermap[1] = 0;megermap[2] = 0;//megermap[5] = 4;
        //    megermap[1] = 0;
        //    megermap[2] = 0;megermap[9] = 0;megermap[12] = 0;

        //    megermap[10] = 7;megermap[17] = 7;megermap[7] = 7;megermap[15] = 7;megermap[16] = 7;megermap[8] = 7;megermap[11] = 7;
        //    megermap[3] = 7;

        //    megermap[4] = 7;megermap[5] = 7;megermap[13] = 7;megermap[14] = 7;


        for(auto &a:surf_group1)a = megermap[a];
        for(auto &a:surf_group2)a = megermap[a];

        vector<bool>restFace(mixSurf.n_faces,false);

        for(int i=0;i<mixSurf.n_faces;++i){
            if(surf_group1[i]!=surf_group2[i])restFace[i]=true;
        }

        vector<int>newVer(mixSurf.n_vertices,-1);
        for(int i=0;i<mixSurf.n_faces;++i)if(restFace[i]){
            auto p_fv= mixSurf.fv_begin(i);
            for(int j=0;j<3;++j)newVer[p_fv[j]] = 0;

        }
        int neInd = 0;
        for(auto &a:newVer)if(a==0)a=neInd++;

        vector<double>verticesPos;
        vector<uint>faces2vertices;

        for(int i=0;i<mixSurf.n_vertices;++i)if(newVer[i]!=-1){
            auto p_v = mixSurf.v_begin(i);
            for(int j=0;j<3;++j)verticesPos.push_back(p_v[j]);
        }
        for(int i=0;i<mixSurf.n_faces;++i)if(restFace[i]){
            auto p_fv= mixSurf.fv_begin(i);
            for(int j=0;j<3;++j)faces2vertices.push_back(newVer[p_fv[j]]);
        }

        vector<int>surf_groupN1,surf_groupN2;
        for(int i=0;i<mixSurf.n_faces;++i)if(restFace[i]){
            surf_groupN1.push_back(surf_group1[i]);
            surf_groupN2.push_back(surf_group2[i]);
        }


        mixSurf.faces2vertices = faces2vertices;
        mixSurf.vertices = verticesPos;
        mixSurf.setparameters();
        surf_group1 = surf_groupN1;
        surf_group2 = surf_groupN2;
    }

    if(1){
        set<int>matSet;
        for(auto a:surf_group1)matSet.insert(a);
        for(auto a:surf_group2)matSet.insert(a);

        map<int,int>mappM;
        int newMInd = 0;
        for(auto a:matSet)mappM[a] = newMInd++;
        for(auto &a:surf_group1)a = mappM[a];
        for(auto &a:surf_group2)a = mappM[a];
        numofSurf = matSet.size();

    }

    mixSurf.face_Material1 = surf_group1;
    mixSurf.face_Material2 = surf_group2;








}




int MultiSurface::splitContour(string infilename, string outfilename, vector<double>&points, vector<uint> &edges){

    ifstream reader(infilename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Contour file " << infilename << endl;
        return -1;
    }
    int n_c;
    reader>>n_c;

    vector<Contour>mainC(n_c);

    for(int i=0;i<n_c;++i)mainC[i].ReadContour(reader);

    int n_materials = mainC[0].n_materials;

    vector< vector<Contour> >mainSC(n_c,vector<Contour>(n_materials));


    CrossSections CCC(mainC);
    CCC.Stackup(points,edges);
    //CCC.MapMaterials();
    //CCC.WriteCrossSections(outfilename);
    //CCC.naiveConstruction();

    for(int i=0;i<n_c;++i)mainC[i].splitContour(mainSC[i]);

    vector<int>m_n_c(n_materials,0);

    for(int i=0;i<n_materials;++i)for(int j=0;j<n_c;++j)if(mainSC[j][i].n_vertices!=0)m_n_c[i]++;

    for(int i=0;i<n_materials;++i){
        string outCname = outfilename + string("_m") + to_string(i) + string(".contour");
        ofstream outer(outCname.data(), ofstream::out);
        if (!outer.good()) {
            cout << "Can not open the Contour file " << outCname << endl;
            return -1;
        }
        outer<<m_n_c[i]<<endl;
        for(int j=0;j<n_c;++j)mainSC[j][i].WriteContour(outer);
        outer.close();
    }

    return n_materials;

}



/**********************************************************************/
vector<bool> MultiCellTopo::InterpretVisibility(bool isCS, bool isAllSurf, bool isNonC, int pickCell){
    vector<bool>re(CellTopo.size()*2+1, false);

    if(pickCell<0 || pickCell >= CellTopo.size()){
        cout<<"pickCell Error, InterpretVisibility failed."<<endl;
        return re;
    }
    int nCell = CellTopo.size();
    if(isCS)re[0] = true;
    if(isAllSurf){

        for(int i=0;i<nCell;++i)re[i+1] = true;
        if(isNonC)for(int i=0;i<nCell;++i)re[i+1+nCell] = true;

    }else{

        re[pickCell+1] = true;
        if(isNonC)re[pickCell+1+nCell] = true;
    }

    return re;

}


void MultiCellTopo::BuildDisplay(infoSurfDisp info,vector<int>&topo,int pickMat){
    //mixSurf.BuildDisplay(info);
    //for(auto &a:surfbuffer)cout<<a.weighted_fcolor.size()<<endl;
    //for(auto &a:CellTopo)for(auto&b : a)b.BuildDisplay(info,false);

    for(int i=0;i<CellTopo.size();++i)if(nTopo[i]>0)for(auto&b : CellTopo[i])b.BuildDisplay(info,false);


    display_vertices.resize(CellTopo.size()*2+1);
    display_edges.resize(CellTopo.size()*2+1);
    display_faces.resize(CellTopo.size()*2+1);
    display_normal.resize(CellTopo.size()*2+1);
    display_vcolor.resize(CellTopo.size()*2+1);

    int pickMatInd = pickMat+1;
    if(pickMat==-1)pickMatInd = nMat + 2;
    //cout<<pickMat<<endl;
    //cout<<pickMatInd<<endl;
    //pickMatInd = 1;
    for(int i=0;i<CellTopo.size();++i)if(nTopo[i]>0){
        //if(i!=1)continue;
        int pickTopo = topo[i];
        display_vertices[i+1] = (CellTopo[i][pickTopo].getDisplayVertices())->at(pickMatInd);
        display_edges[i+1] = (CellTopo[i][pickTopo].getDisplayEdges())->at(pickMatInd);
        display_faces[i+1] = (CellTopo[i][pickTopo].getDisplayFaces())->at(pickMatInd);
        display_normal[i+1] =( CellTopo[i][pickTopo].getDisplayVerticesNormal())->at(pickMatInd);
        display_vcolor[i+1] = (CellTopo[i][pickTopo].getDisplayColor())->at(pickMatInd);
        //cout<<"iii: "<<i<<' '<<(CellTopo[i][pickTopo].getDisplayVertices())->size()<<endl;
    }

    for(int i=0;i<CellTopo.size();++i)if(nTopo[i]>0){
        //if(i!=1)continue;
        int pickTopo = topo[i];
        int nnind = (CellTopo[i][pickTopo].getDisplayVertices())->size()-2;
        display_vertices[i+1+CellTopo.size()] = (CellTopo[i][pickTopo].getDisplayVertices())->at(nnind);
        display_edges[i+1+CellTopo.size()] = (CellTopo[i][pickTopo].getDisplayEdges())->at(nnind);
        display_faces[i+1+CellTopo.size()] = (CellTopo[i][pickTopo].getDisplayFaces())->at(nnind);
        display_normal[i+1+CellTopo.size()] =( CellTopo[i][pickTopo].getDisplayVerticesNormal())->at(nnind);
        display_vcolor[i+1+CellTopo.size()] = (CellTopo[i][pickTopo].getDisplayColor())->at(nnind);
        //cout<<"iii: "<<i<<' '<<(CellTopo[i][pickTopo].getDisplayVertices())->size()<<endl;
    }


    crossS.BuildDisplay(false,true);

    display_vertices[0] = crossS.getDisplayVertices();
    display_edges[0] = crossS.getDisplayEdges();
    display_faces[0] = crossS.getDisplayFaces();
    display_normal[0] = crossS.getDisplayVerticesNormal();
    display_vcolor[0] = crossS.getDisplayColor();



}






void MultiCellTopo::ReadCellTopo(string nameprefix, vector<int>n_Topos, bool isSmooth){

    CellTopo.clear();
    CellTopo.resize(n_Topos.size());
    string ed(".suf");

    Mesh tmp;
    MultiSurface a;

    for(int i=0;i<n_Topos.size();++i)if(n_Topos[i]>0){

        string s_cell = string("_")+to_string(i);
        string s_topo = string("_")+to_string(0);
        //nameprefix+s_cell+s_topo+ed;

        a.readSufFile(nameprefix+s_cell+s_topo+ed);
        tmp.addMesh(a.mixSurf);
        nMat = max(nMat,a.GetNumofSurfaces());
    }

    tmp.GetRescaleInfo(lscale, pcenter);


    //cout<<lscale<<endl;
    //cout<<pcenter[0]<<' '<<pcenter[1]<<' '<<pcenter[2]<<endl;



    nTopo = n_Topos;
    for(int i=0;i<n_Topos.size();++i){
        CellTopo[i].resize(n_Topos[i]);
        string s_cell = string("_")+to_string(i);
        for(int j=0;j<n_Topos[i];++j){
            string s_topo = string("_")+to_string(j);
            CellTopo[i][j].Initialize(nameprefix+s_cell+s_topo+ed,false,true,true,nMat,true,lscale,pcenter);
        }
    }


    for(auto &a:CellTopo)if(a.size()!=0){
        a[0].GetColorDegree(mat_colordegree);
    }
    //CellTopo[0][0].GetColorDegree(mat_colordegree);

    crossS.reset();


    for(int i=0;i<n_Topos.size();++i)if(n_Topos[i]>0){

        crossS.AddCurve(CellTopo[i][0].crossSection);
    }
    crossS.setparameters();
    crossS.BuildEdges(false);


    if(!isSmooth)return;

    GlobalSmoothing();


    //for(auto a:TopoCurInd)cout<<a<<' ';cout<<endl;
    //for(auto &b:pSuf.mappingsToGlobal){for(auto a:b)cout<<a<<' ';cout<<endl;}


}


void MultiCellTopo::ReadCellTopo_picked(string nameprefix, vector<int>&n_Topos, vector<int>&picked_Topos, bool isSmooth){

    CellTopo.clear();
    CellTopo.resize(n_Topos.size());
    string ed(".suf");

    Mesh tmp;
    MultiSurface a;

    for(int i=0;i<n_Topos.size();++i)if(n_Topos[i]>0){

        string s_cell = string("_")+to_string(i);
        string s_topo = string("_")+to_string(0);
        nameprefix+s_cell+s_topo+ed;

        a.readSufFile(nameprefix+s_cell+s_topo+ed);
        tmp.addMesh(a.mixSurf);
        nMat = max(nMat,a.GetNumofSurfaces());
    }

    tmp.GetRescaleInfo(lscale, pcenter);


    //cout<<lscale<<endl;
    //cout<<pcenter[0]<<' '<<pcenter[1]<<' '<<pcenter[2]<<endl;



    nTopo = n_Topos;
    for(int i=0;i<n_Topos.size();++i){
        CellTopo[i].resize(n_Topos[i]);
        string s_cell = string("_")+to_string(i);

        assert(picked_Topos[i]<n_Topos[i]);

        string s_topo = string("_")+to_string(picked_Topos[i]);
        CellTopo[i][picked_Topos[i]].Initialize(nameprefix+s_cell+s_topo+ed,false,false,true,nMat,true,lscale,pcenter);
    }


    for(int i=0;i<n_Topos.size();++i)if(CellTopo[i].size()!=0){
        CellTopo[i][picked_Topos[i]].GetColorDegree(mat_colordegree);
    }
    //CellTopo[0][0].GetColorDegree(mat_colordegree);

    crossS.reset();


    for(int i=0;i<n_Topos.size();++i)if(n_Topos[i]>0){

        crossS.AddCurve(CellTopo[i][picked_Topos[i]].crossSection);
    }
    crossS.setparameters();
    crossS.BuildEdges(false);


    if(!isSmooth)return;

    GlobalSmoothing();


    //for(auto a:TopoCurInd)cout<<a<<' ';cout<<endl;
    //for(auto &b:pSuf.mappingsToGlobal){for(auto a:b)cout<<a<<' ';cout<<endl;}


}

void MultiCellTopo::WriteAllSurface(string nameprefix){


    for(int i=0;i<nTopo.size();++i){

        string s_cell = string("_")+to_string(i);
        for(int j=0;j<nTopo[i];++j){
            string s_topo = string("_")+to_string(j);
            CellTopo[i][j].WriteSuf(nameprefix+s_cell+s_topo);
        }
    }

}

void MultiCellTopo::ReReadSingleCellTopo(string nameprefix, int celli, int c_nTopo){


    cout<<"Function Disable"<<endl;
//    string ed(".suf");

//    nTopo[celli] = c_nTopo;
//    CellTopo[celli].clear();
//    CellTopo[celli].resize(c_nTopo);
//    string s_cell = string("_")+to_string(celli);
//    for(int j=0;j<c_nTopo;++j){
//        string s_topo = string("_")+to_string(j);
//        CellTopo[celli][j].Initialize(nameprefix+s_cell+s_topo+ed,false,false,true,nMat,true,lscale,pcenter);
//    }

//    SufStructure pSuf;
//    MultiSurface multiSuf;
//    vector<vector<double>>TopoVs; vector<vector<uint>>TopoFs;
//    vector<vector<int> > TopoFMs;vector<vector<uint>>TopoCtrs;
//    vector<int>TopoCurInd(nTopo.size(),-1);
//    for(int j=0;j<nTopo.size();++j)if(TopoCurInd[j]+1<nTopo[j]){TopoCurInd[j]+=1;}
//    vector<bool>hasChange(nTopo.size(),false);hasChange[celli] = true;
//    for(int i=0;i<c_nTopo;++i){
//        TopoCurInd[celli]=i;
//        GetAllCellTopos(TopoCurInd, TopoVs, TopoFs,TopoFMs,TopoCtrs);
//        pArr.CreateTotalSurf(TopoVs, TopoFs,TopoFMs,TopoCtrs,false,&pSuf);
//        multiSuf.ImportFromSuf(pSuf,true);
//        UpdateVerticesPosition(multiSuf.mixSurf.vertices,TopoCurInd,pSuf.mappingsToGlobal,hasChange);
//    }


}

void MultiCellTopo::UpdateVerticesPosition(vector<double>&newVpos,vector<int>&pickTopos,vector<vector<int>>&mappingToGlobal,vector<bool>&isChanged){

    for(int i=0;i<pickTopos.size();++i)if(pickTopos[i]>=0)if(isChanged[i]){

        vector<double>&toVs = CellTopo[i][pickTopos[i]].mixSurf.vertices;
        auto &mappingsTog = mappingToGlobal[i];
        auto p_newVd = newVpos.data();
        //toVs.clear();toVs.resize(mappingsTog.size()*3,0);
        auto p_vd = toVs.data();
        for(int j=0;j<mappingsTog.size();++j){

            copyVec(p_newVd+mappingsTog[j]*3,p_vd+j*3);
        }
        CellTopo[i][pickTopos[i]].ComputeSmallMVnormal();
        CellTopo[i][pickTopos[i]].ReAllocation(false);
        //cout<<"CACACA"<<endl;
    }


}
void MultiCellTopo::GlobalSmoothing(){


    cout<<"Function Disable"<<endl;
    exit(12);
//    SufStructure pSuf;
//    MultiSurface multiSuf;

//    vector<vector<double>>TopoVs; vector<vector<uint>>TopoFs;
//    vector<vector<int> > TopoFMs;vector<vector<uint>>TopoCtrs;

//    vector<int>TopoCurInd(nTopo.size(),-1);

//    int maxNCellTopo = *max_element(nTopo.begin(),nTopo.end());
//    for(int i=0;i<maxNCellTopo;++i){
//        vector<bool>hasChange(nTopo.size(),false);
//        for(int j=0;j<nTopo.size();++j)if(TopoCurInd[j]+1<nTopo[j]){TopoCurInd[j]+=1;hasChange[j] = true;}
//        GetAllCellTopos(TopoCurInd, TopoVs, TopoFs,TopoFMs,TopoCtrs);
//        pArr.CreateTotalSurf(TopoVs, TopoFs,TopoFMs,TopoCtrs,false,&pSuf);
//        multiSuf.ImportFromSuf(pSuf,true);
//        UpdateVerticesPosition(multiSuf.mixSurf.vertices,TopoCurInd,pSuf.mappingsToGlobal,hasChange);
//    }

}

void MultiCellTopo::GetAllCellTopos(vector<int>&pickTopos, vector<vector<double>>&TopoVs, vector<vector<uint>>&TopoFs, vector<vector<int> > &TopoFMs,vector<vector<uint>>&TopoCtrs){

    if(pickTopos.size()!=CellTopo.size()){
        cout<<"please input a complete choice of all cells"<<endl;
        return;
    }
    TopoVs.clear();
    TopoFs.clear();
    TopoCtrs.clear();
    TopoFMs.clear();


    TopoVs.resize(CellTopo.size());
    TopoFs.resize(CellTopo.size());
    TopoFMs.resize(CellTopo.size());
    TopoCtrs.resize(CellTopo.size());

    for(int i=0;i<CellTopo.size();i++)if(pickTopos[i]>=0){

        TopoVs[i] = CellTopo[i][pickTopos[i]].mixSurf.vertices;
        TopoFs[i] = CellTopo[i][pickTopos[i]].mixSurf.faces2vertices;
        int nf = CellTopo[i][pickTopos[i]].mixSurf.n_faces;
        auto &sm1 =  CellTopo[i][pickTopos[i]].surf_group1;
        auto &sm2 =  CellTopo[i][pickTopos[i]].surf_group2;
        for(int j=0;j<nf;++j){
            TopoFMs[i].push_back(sm1[j]);
            TopoFMs[i].push_back(sm2[j]);
        }

        TopoCtrs[i] = CellTopo[i][pickTopos[i]].curvee2v;


    }






}


void MultiCellTopo::CutSurfaceByPlane(int celli,int Topoi,double *para,vector<double>&outCtrV,vector<uint>&outCtrE,vector<int>&outCtrEMat){

    CellTopo[celli][Topoi].mixSurf.CutSurfaceByPlane(para,outCtrV,outCtrE,outCtrEMat);

}

void MultiCellTopo::GetRescaleInfo(double &outlscale, double *outpcenter){

    outlscale = lscale;
    for(int i=0;i<3;++i)outpcenter[i] = pcenter[i];


}

void MultiCellTopo::GetColorDegree(vector<double>&out_colordegree){


    out_colordegree = mat_colordegree;

}
double MultiCellTopo::GetLabel2Colordegree(int label){

    return mat_colordegree[label];
}
}//n_rf
