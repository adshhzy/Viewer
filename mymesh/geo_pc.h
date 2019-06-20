#ifndef GEO_PC_H
#define GEO_PC_H



#include<vector>
#include<math.h>
#include"my_mesh.h"
#include"InfoStruct.h"
#include<unordered_map>


using namespace std;
using namespace MyUtility;


namespace n_rf {


class PointCloud{


public:

    int n_vertices;

    vector<double>vertices;
    vector<double>vnormals;

    vector<double>vnormals_uninorm;


    vector<double>vnormals_cmp;
    vector<double>vnormals_uninorm_cmp;

    vector<double>cmp_diff;
public:
    vector<uint>pc_graph;
    vector<double>pc_graph_weight;

public:
    double avedist;

public:




    vector<double>display_vertices;
    vector<double>display_vnormals;
    vector<uint>display_edges;
    vector<uint>display_faces;
    vector<unsigned char>display_color;


public:
    double disscale = 1.0;

    static bool isload;
    static Mesh sphere;
    static Mesh cylinder;
    static Mesh cone;
    static Mesh disk;
    int n_ver_f1;
    int n_col_f1;
    int n_face_f1;
    int n_vnor_f1;

private:
    bool load();


public:
    PointCloud();

    void setparameters();

    bool ReadPointCloud(string filename);

    bool ImportPointCloud(vector<double>& vertices_in, vector<double>& vn_in);
    bool ImportPointCloud(vector<double>* vertices_in, vector<double>* vn_in);

    bool ImportPointCloud(vector<double>& vertices_in, vector<double>& vn_in, vector<double>& vn_in_cmp);
    bool ImportPointCloud(vector<double>* vertices_in, vector<double>* vn_in, vector<double>* vn_in_cmp);

    bool ImportPointCloud(vector<double>* vertices_in, vector<double>* vn_in, vector<uint>* pc_graph_in, vector<double>* pc_graphw_in);

    bool ImportPointCloud(vector<double>* vertices_in, vector<double>* vn_in, vector<double>* vn_in_cmp, vector<uint>* pc_graph_in, vector<double>* pc_graphw_in);

    bool ImportGraph(vector<uint>* pc_graph_in, vector<double>* pc_graphw_in);

    void Set_Avedist(double avedist_in);
    double Get_Avedist();
    void Compute_CMP_Diff();

    void FlipNormal();
public:

    void BuildDisplay(int point_method, bool isdoublesidec, bool isshownormal, bool isunitnormal, double norlen, int graphmethod = 0, double plotpercentage = 0., double transparent_cut = 1.);
    void BuildDisplay(vector<double>&vcolordegree, int method = 0);

    void BuildDisplay(infoPCDisp info);
    void ClearDisplayBuffer();


    vector<double>* getDisplayVertices(){return &display_vertices;}
    vector<double>* getDisplayVerticesNormal(){return &display_vnormals;}
    vector<uint>* getDisplayEdges(){return &display_edges;}
    vector<uint>* getDisplayFaces(){return &display_faces;}
    vector<unsigned char>* getDisplayColor(){return &display_color;}


public:

    inline double* v_begin(int v_ind){ return &(vertices[v_ind * 3]); }
    inline double* v_end(int v_ind){ return &(vertices[(v_ind+1) * 3]); }

    inline double* vn_begin(int v_ind){ return &(vnormals[v_ind * 3]); }
    inline double* vn_end(int v_ind){ return &(vnormals[(v_ind+1) * 3]); }

    inline double* vnu_begin(int v_ind){ return &(vnormals_uninorm[v_ind * 3]); }
    inline double* vnu_end(int v_ind){ return &(vnormals_uninorm[(v_ind+1) * 3]); }
};



class Multi_PointCloud{

public:
    int numofPC;
    vector<PointCloud>PC_buffer;

    unordered_map<string,int>pc_nameind;
    vector<string>pc_names;

    int Import_Multi_PointCloud(vector<vector<double>>* vertices_in, vector<vector<double>>* vn_in);
    int Import_Multi_PointCloud(vector<vector<double>>* vertices_in, vector<vector<double>>* vn_in, vector<vector<double>>* vn_cmp);
    int Import_Multi_PointCloud(vector<string>*pnames,vector<vector<double>>* vertices_in, vector<vector<double>>* vn_in, vector<vector<double>>* vn_cmp);
    int Import_Multi_PointCloud(vector<vector<double>>* vertices_in, vector<vector<double>>* vn_in, vector<vector<uint>>* pc_graph_in, vector<vector<double>>* pc_graphw_in);
    int Import_Multi_Graph(vector<vector<uint>>* pc_graph_in, vector<vector<double>>* pc_graphw_in);

    void BuildDisplay(infoPCDisp info);

    int Get_PCIndex(string fname);

    void FlipNormal();

public:
    vector< vector<double>* > display_vertices;
    vector< vector<double>* >display_normal;
    vector< vector<uint>* >display_edges;

    vector< vector<uint>* >display_field_dot;
    vector< vector<unsigned char>* >display_vcolor;
    vector< vector<uint>* >display_faces;


public:
    vector<vector<double>*>* getDisplayVertices(){return &display_vertices;}
    vector<vector<double>*>* getDisplayVerticesNormal(){return &display_normal;}
    vector<vector<uint>*>* getDisplayEdges(){return &display_edges;}
    vector<vector<uint>*>* getDisplayFaces(){return &display_faces;}
    vector<vector<uchar>*>* getDisplayColor(){return &display_vcolor;}

};





} // n_rf





#endif // GEO_PC_H
