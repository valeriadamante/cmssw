#ifndef LCTDEBUG_CLASSES
#define LCTDEBUG_CLASSES
#include <vector>

class LCTDebugobject{
    
    public:
        LCTDebugobject();
        void SetbendingangleLayer1(int eighthStripDiff_);
        void SetbendingangleLayer2(int eighthStripDiff_);
        void SetGEMClusterKeyStripLayer1(int cl_es_);
        void SetGEMClusterKeyStripLayer2(int cl_es_);
        void SetresidualLayer1(int residualwithalignment);
        void SetresidualLayer2(int residualwithalignment);
        void SetresidualnoalignmentcorrectionLayer1(int residualwithoutalignment_);
        void SetresidualnoalignmentcorrectionLayer2(int residualwithoutalignment_);
        void SetbendinganglenoalignmentcorrectionLayer1(int bendinganglenoalignmentcorrection_);
        void SetbendinganglenoalignmentcorrectionLayer2(int bendinganglenoalignmentcorrection_);
        void Setidentifiers(int KeyWG_, int bx_, int Bend_,int KeyStrip_, int slope_);
        void SetClusterRoll1(int roll_);
        void SetClusterRoll2(int roll_);
        void SetClusterBxLayer1(int bx_);
        void SetClusterBxLayer2(int bx_);
        int GetbendingangleLayer1();
        int GetbendingangleLayer2();
        int GetbendinganglenoalignmentcorrectionLayer1();
        int GetbendinganglenoalignmentcorrectionLayer2();
        void Setlayer2bool(bool isLayer2);
        void SetLayer1Match(bool layer1match);
        void SetLayer2Match(bool layer2match);
        bool Getlayer2bool();
        bool GetLayer1Match();
        bool GetLayer2Match();
        std::vector<int> Getidentifiers();
        int GetGEMClusterKeyStripLayer1();
        int GetGEMClusterKeyStripLayer2();
        int GetCLCTKeyStrip();
        int GetresidualLayer1();
        int GetresidualLayer2();
        int GetresidualnoalignmentcorrectionLayer1();
        int GetresidualnoalignmentcorrectionLayer2();
        int GetClusterRoll1();
        int GetClusterRoll2();
        int GetClusterBxLayer1();
        int GetClusterBxLayer2();
    private:    
        int bendinganglelayer1;
        int bendinganglelayer2;
        int bendinganglenoalignmentcorrectionlayer1;
        int bendinganglenoalignmentcorrectionlayer2;
        int KeyWG;
        int bx;
        int Bend;
        int KeyStrip;
        int slope;
        bool layer2bool;
        bool layer1match;
        bool layer2match;
        int cl_eslayer1;
        int cl_eslayer2;
        int residualwithalignmentlayer1;
        int residualwithalignmentlayer2;
        int residualwithoutalignmentlayer1;
        int residualwithoutalignmentlayer2;
        int cluster_rolllayer1;
        int cluster_rolllayer2;
        int cluster_bxlayer1;
        int cluster_bxlayer2;
};











#endif