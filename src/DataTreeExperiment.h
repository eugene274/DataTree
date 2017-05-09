#ifndef DataTreeExperiment_H
#define DataTreeExperiment_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"
#include "DataTreeDetector.h"

class DataTreeExperiment : public TObject
{
    
public:
  
    DataTreeExperiment();
    ~DataTreeExperiment();
    
    void Init()
    {
	CBM_Reco_Vertex -> Init();
	CBM_Reco_Tracks -> Init();
	CBM_Reco_V0_Tracks_RecoPid -> Init();
	CBM_Reco_V0_Tracks_MCPid -> Init();
	CBM_Reco_PSD -> Init();
	CBM_Reco_TOF -> Init();
	CBM_MC_Vertex -> Init();
	CBM_MC_Tracks -> Init();
	
	NA61_Reco_Main_Vertex -> Init();
	NA61_Reco_Main_Vertex_Tracks -> Init();
	NA61_Reco_Primary_Vertex -> Init();
	NA61_Reco_Primary_Vertex_Tracks -> Init();
	NA61_Reco_PSD -> Init();
	NA61_Reco_BPD -> Init();
	NA61_Reco_Triggers -> Init();
	NA61_Reco_WFA -> Init();
    }
    
    void ClearEvent()
    {
	CBM_Reco_Vertex -> ClearEvent();
	CBM_Reco_Tracks -> ClearEvent();
	CBM_Reco_V0_Tracks_RecoPid -> ClearEvent();
	CBM_Reco_V0_Tracks_MCPid -> ClearEvent();
	CBM_Reco_PSD -> ClearEvent();
	CBM_Reco_TOF -> ClearEvent();
	CBM_MC_Vertex -> ClearEvent();
	CBM_MC_Tracks -> ClearEvent();
	
	NA61_Reco_Main_Vertex -> ClearEvent();
	NA61_Reco_Main_Vertex_Tracks -> ClearEvent();
	NA61_Reco_Primary_Vertex -> ClearEvent();
	NA61_Reco_Primary_Vertex_Tracks -> ClearEvent();
	NA61_Reco_PSD -> ClearEvent();
	NA61_Reco_BPD -> ClearEvent();
	NA61_Reco_Triggers -> ClearEvent();
	NA61_Reco_WFA -> ClearEvent();
	
	RunId = -999;
	EventId = -999;
	EventTimestamp = -999;
	RPAngle = -999;
	Bx = -999;
	By = -999;
    }
    int GetRunId(){return RunId;}
    int GetEventId(){return EventId;}
    double GetEventTimestamp(){return EventTimestamp;}
    double GetRPAngle(){return RPAngle;}
    double GetBx(){return Bx;}
    double GetBy(){return By;}
    
    void SetRunId(int pValue){ RunId = pValue; }
    void SetEventId(int pValue){ EventId = pValue; }
    void SetEventTimestamp(double pValue){ EventTimestamp = pValue; }
    void SetRPAngle(double pValue){ RPAngle = pValue; }
    void SetImpactParameter(double pX, double pY){ Bx = pX; By = pY; }
    void SetBx(double pX){ Bx = pX; }
    void SetBy(double pY){ By = pY; }
    
    DataTreeDetector* Get_CBM_Reco_Vertex(){ return CBM_Reco_Vertex; }
    DataTreeDetector* Get_CBM_Reco_Tracks(){ return CBM_Reco_Tracks; }
    DataTreeDetector* Get_CBM_Reco_V0_Tracks_RecoPid(){ return CBM_Reco_V0_Tracks_RecoPid; }
    DataTreeDetector* Get_CBM_Reco_V0_Tracks_MCPid(){ return CBM_Reco_V0_Tracks_MCPid; }
    DataTreeDetector* Get_CBM_Reco_PSD(){ return CBM_Reco_PSD; }
    DataTreeDetector* Get_CBM_Reco_TOF(){ return CBM_Reco_TOF; }
    DataTreeDetector* Get_CBM_MC_Vertex(){ return CBM_MC_Vertex; }
    DataTreeDetector* Get_CBM_MC_Tracks(){ return CBM_MC_Tracks; }
    DataTreeDetector* Get_NA61_Reco_Main_Vertex(){ return NA61_Reco_Main_Vertex; }
    DataTreeDetector* Get_NA61_Reco_Main_Vertex_Tracks(){ return NA61_Reco_Main_Vertex_Tracks; }
    DataTreeDetector* Get_NA61_Reco_Primary_Vertex(){ return NA61_Reco_Primary_Vertex; }
    DataTreeDetector* Get_NA61_Reco_Primary_Vertex_Tracks(){ return NA61_Reco_Primary_Vertex_Tracks; }
    DataTreeDetector* Get_NA61_Reco_PSD(){ return NA61_Reco_PSD; }
    DataTreeDetector* Get_NA61_Reco_BPD(){ return NA61_Reco_BPD; }
    DataTreeDetector* Get_NA61_Reco_Triggers(){ return NA61_Reco_Triggers; }
    DataTreeDetector* Get_NA61_Reco_WFA(){ return NA61_Reco_WFA; }
    
private:

    int RunId;
    int EventId;
    double EventTimestamp;
    double RPAngle;			//Reaction plane angle
    double Bx;		//Impact parameter
    double By;		//Impact parameter
  
    DataTreeDetector* CBM_Reco_Vertex;
    DataTreeDetector* CBM_Reco_Tracks;
    DataTreeDetector* CBM_Reco_V0_Tracks_RecoPid;
    DataTreeDetector* CBM_Reco_V0_Tracks_MCPid;
    DataTreeDetector* CBM_Reco_PSD;
    DataTreeDetector* CBM_Reco_TOF;
    DataTreeDetector* CBM_MC_Vertex;
    DataTreeDetector* CBM_MC_Tracks;
    
    DataTreeDetector* NA61_Reco_Main_Vertex;
    DataTreeDetector* NA61_Reco_Main_Vertex_Tracks;
    DataTreeDetector* NA61_Reco_Primary_Vertex;
    DataTreeDetector* NA61_Reco_Primary_Vertex_Tracks;
    DataTreeDetector* NA61_Reco_PSD;
    DataTreeDetector* NA61_Reco_BPD;
    DataTreeDetector* NA61_Reco_Triggers;
    DataTreeDetector* NA61_Reco_WFA;
    
    ClassDefNV(DataTreeExperiment, 1)
};

#endif