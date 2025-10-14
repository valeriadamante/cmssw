#include "DataFormats/LCTDebug/interface/LCTDebug.h"


LCTDebugobject::LCTDebugobject(){
    // bendingangle=-99;
    // bendinganglenoalignmentcorrection = -99;
    KeyWG=0;
    bx=0;
    Bend=0;
    KeyStrip=0;
    slope=0;
}

void LCTDebugobject::SetbendingangleLayer1(int eighthStripDiff_){
    bendinganglelayer1 = eighthStripDiff_;
}

void LCTDebugobject::SetbendingangleLayer2(int eighthStripDiff_){
    bendinganglelayer2 = eighthStripDiff_;
}

void LCTDebugobject::SetbendinganglenoalignmentcorrectionLayer1(int bendinganglenoalignmentcorrection_){
    bendinganglenoalignmentcorrectionlayer1 = bendinganglenoalignmentcorrection_;
}

void LCTDebugobject::SetbendinganglenoalignmentcorrectionLayer2(int bendinganglenoalignmentcorrection_){
    bendinganglenoalignmentcorrectionlayer2 = bendinganglenoalignmentcorrection_;
}

void LCTDebugobject::Setidentifiers(int KeyWG_, int bx_, int Bend_,int KeyStrip_, int slope_){
    KeyWG = KeyWG_;
    bx = bx_;
    Bend = Bend_;
    KeyStrip = KeyStrip_;
    slope = slope_;
}

void LCTDebugobject::SetGEMClusterKeyStripLayer1(int cl_es_){
    cl_eslayer1 = cl_es_;
}

void LCTDebugobject::SetGEMClusterKeyStripLayer2(int cl_es_){
    cl_eslayer2 = cl_es_;
}

void LCTDebugobject::SetresidualLayer1(int residualwithalignment_){
    residualwithalignmentlayer1 = residualwithalignment_;
}

void LCTDebugobject::SetresidualLayer2(int residualwithalignment_){
    residualwithalignmentlayer2 = residualwithalignment_;
}

void LCTDebugobject::SetresidualnoalignmentcorrectionLayer1(int residualwithoutalignment_){
    residualwithoutalignmentlayer1 = residualwithoutalignment_;
}

void LCTDebugobject::SetresidualnoalignmentcorrectionLayer2(int residualwithoutalignment_){
    residualwithoutalignmentlayer2 = residualwithoutalignment_;
}

void LCTDebugobject::SetClusterRoll1(int roll_){
    cluster_rolllayer1 = roll_;
}

void LCTDebugobject::SetClusterRoll2(int roll_){
    cluster_rolllayer2 = roll_;
}

void LCTDebugobject::SetClusterBxLayer1(int bx_){
    cluster_bxlayer1 = bx_;
}

void LCTDebugobject::SetClusterBxLayer2(int bx_){
    cluster_bxlayer2 = bx_;
}

void LCTDebugobject::SetLayer1Match(bool match){
    layer1match = match;
}

void LCTDebugobject::SetLayer2Match(bool match){
    layer2match = match;
}

int LCTDebugobject::GetbendingangleLayer1(){
    return bendinganglelayer1;
}

int LCTDebugobject::GetbendingangleLayer2(){
    return bendinganglelayer2;
}

int LCTDebugobject::GetbendinganglenoalignmentcorrectionLayer1(){
    return bendinganglenoalignmentcorrectionlayer1;
}

int LCTDebugobject::GetbendinganglenoalignmentcorrectionLayer2(){
    return bendinganglenoalignmentcorrectionlayer2;
}

std::vector<int> LCTDebugobject::Getidentifiers(){
    std::vector<int> output;
    output.push_back(KeyWG);
    output.push_back(bx);
    output.push_back(Bend);
    output.push_back(KeyStrip);
    output.push_back(slope);
    return output;
}

void LCTDebugobject::Setlayer2bool(bool isLayer2){
    layer2bool = isLayer2;
}

bool LCTDebugobject::Getlayer2bool(){
    return layer2bool;
}

int LCTDebugobject::GetGEMClusterKeyStripLayer1(){
    return cl_eslayer1;
}

int LCTDebugobject::GetGEMClusterKeyStripLayer2(){
    return cl_eslayer2;
}

int LCTDebugobject::GetresidualLayer1(){
    return residualwithalignmentlayer1;
}

int LCTDebugobject::GetresidualLayer2(){
    return residualwithalignmentlayer2;
}

int LCTDebugobject::GetresidualnoalignmentcorrectionLayer1(){
    return residualwithoutalignmentlayer1;
}

int LCTDebugobject::GetresidualnoalignmentcorrectionLayer2(){
    return residualwithoutalignmentlayer2;
}

int LCTDebugobject::GetClusterRoll1(){
    return cluster_rolllayer1;
}

int LCTDebugobject::GetClusterRoll2(){
    return cluster_rolllayer2;
}

int LCTDebugobject::GetClusterBxLayer1(){
    return cluster_bxlayer1;
}

int LCTDebugobject::GetClusterBxLayer2(){
    return cluster_bxlayer2;
}

bool LCTDebugobject::GetLayer1Match(){
    return layer1match;
}

bool LCTDebugobject::GetLayer2Match(){
    return layer2match;
}