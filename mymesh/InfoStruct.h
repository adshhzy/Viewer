#ifndef INFOSTRUCT_H
#define INFOSTRUCT_H



class infoSet{

public:
    bool isrescale;
    bool isvmf;
    bool isrmf;
    bool isnnormal;
    bool isbinormal;
    bool istangent;
    bool isrvector;
    bool isovector;
    bool isVector;
    bool isSurface;
    double veclen;
    double vecthickness;
    double boxlen;
    double boxthickness;

    bool isshowall;
    int pickid;

    bool isSurfColor;
    bool isColorFirst;

    infoSet(){}
    infoSet(bool isrescale,bool isvmf, bool isrmf, bool isnnormal,bool isbinormal, bool istangent,
            bool isrvector,bool isovector,
            bool isVector, bool isSurface,
            double veclen,double vecthickness,
            double boxlen,double boxthickness,
            bool isshowall,int pickid,
            bool isSurfColor, bool isColorFirst):
        isrescale(isrescale),
        isvmf(isvmf),isrmf(isrmf),
        isnnormal(isnnormal),isbinormal(isbinormal),istangent(istangent),
        isrvector(isrvector),isovector(isovector),
        isVector(isVector),isSurface(isSurface),
        veclen(veclen),vecthickness(vecthickness),
        boxlen(boxlen),boxthickness(boxthickness),
        isshowall(isshowall),pickid(pickid),
        isSurfColor(isSurfColor),isColorFirst(isColorFirst)
    {}







};

class infoFrame{
public:
    int init_div;
    int subdiv;
    int iter;
    bool isrefine;
    bool isrotate;
    bool isopt;
    int mode;
    double alpha;
    bool isInverseInit;
    bool isRandomInit;
    bool isReInit;

    infoFrame(){}

    infoFrame(int init_div, int subdiv,int iter,bool isrefine,bool isrotate,bool isopt,int mode,double alpha,
              bool isReInit = false,bool isRandomInit = false,bool isInverseInit = false):
        init_div(init_div),subdiv(subdiv),iter(iter),
        isrefine(isrefine),isrotate(isrotate),isopt(isopt),
        mode(mode),alpha(alpha),isInverseInit(isInverseInit),
        isRandomInit(isRandomInit),isReInit(isReInit)
    {}


};



class infoSurfDisp{
public:
    bool isField;
    bool isNormal;
    bool isSurface;
    bool isWire;
    bool isSingularity;
    bool ismark;
    int length;
    int width;
    int upnormal;

    infoSurfDisp():length(30),upnormal(30){}

    infoSurfDisp(bool isField,bool isNormal,bool isSurface,bool isWire,bool isSingularity,bool ismark,int length,int width,int upnormal):
        isField(isField),isNormal(isNormal),isSurface(isSurface),isWire(isWire),isSingularity(isSingularity),ismark(ismark),length(length),width(width),upnormal(upnormal)
    {}


};

class infoSurfCompute{
public:
    bool isGS;

    infoSurfCompute():isGS(false){}

    infoSurfCompute(bool isGS):
        isGS(isGS)
    {}


};

class infoSurfOutput{
public:
    int n_vertices;
    int n_faces;
    int n_edges;
    double OptTime;
    double FinalEnergy;
    int numOfSingularity;

};

class infoVolDisp{
public:
    bool isField;
    bool isNormal;
    bool isSurface;
    bool isWire;
    bool isSingularity;
    bool isIso;
    int length;
    int upnormal;
    bool isSlice;
    int xyz;
    int slicenum;
    double thickness;
    bool isLineMode;
    //int upnormal;

    infoVolDisp():length(30),isIso(false){}

    infoVolDisp(bool isField,bool isNormal,bool isSurface,bool isWire,bool isSingularity,
                bool isIso,int length,int upnormal,
                bool isSlice,int xyz,int slicenum,double thickness,bool isLineMode):
        isField(isField),isNormal(isNormal),isSurface(isSurface),isWire(isWire),isSingularity(isSingularity),
        isIso(isIso),length(length),upnormal(upnormal),
        isSlice(isSlice),xyz(xyz),slicenum(slicenum),thickness(thickness),isLineMode(isLineMode)
    {}


};



class infoPCDisp{
public:
    int point_method;
    //bool isUseball;
    bool isDoubleside;
    bool isNormal;
    bool isUnitNormal;
    double nor_len;
    int isgraph;
    double plotpercentage;
    double transparent_cut;

    infoPCDisp():point_method(0),isNormal(true){}
    infoPCDisp(int point_method, bool isDoubleside, bool isNormal, bool isUnitNormal, double nor_len, int isgraph, double plotpercentage,double transparent_cut):
        point_method(point_method),isDoubleside(isDoubleside),isNormal(isNormal),isUnitNormal(isUnitNormal),nor_len(nor_len),isgraph(isgraph), plotpercentage(plotpercentage),transparent_cut(transparent_cut){}
};

#endif // INFOSTRUCT_H
