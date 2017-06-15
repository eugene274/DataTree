#ifndef DataTreeConverterNA61_H
#define DataTreeConverterNA61_H 1

#include <iostream>

#include "TClonesArray.h"

#include "TFile.h"

#include "TString.h"
#include "TChain.h"
#include "TTree.h"
#include "DataTreeEvent.h"
#include "DataTreeTrack.h"

class TClonesArray;

class TDirectory;



class DataTreeConverterNA61
{
    
public:
  
    DataTreeConverterNA61();
    ~DataTreeConverterNA61();
    
    void SetInputFile(TString filename) { sInputFileName = filename; }
    void SetOutputFile(TString filename) { sOutputFileName = filename; }
    void SetNEvents(int fValue) { nEvents = fValue; }
    void Body();
      virtual void Init();
	void Chain_Init();
	void OutputTree_Init();
	void DataTreeEvent_Init();
      virtual void Exec(int event);
	void Clear_Event();
	void Read_Event();
      virtual void Finish();
      
private:
    
    DataTreeEvent* DTEvent;
  
    int fCurEvent;
    int nEvents;
    int nRealEvents;
    
    static const int nInfinity = 10000000;
    static const int nUndefinedValue = -999;
    
    static const int nRawReco = 2;
    static const int nTriggers_Simple = 6;
    static const int nTriggers = 6;
    static const int nBPD = 3;
    static const int nBPDcomponents = 3;
    static const int nPSD_Modules = 45;
    static const int nPSD_Sections = 10;
    static const int nTPC_Tracks = 2000;
    static const int nMaxWFAsignals = 2000;

	static const int nTriggers_Combined = 4;
    
    TString sInputFileName;
    TChain* fChain;
    
    TString sOutputFileName;
    TFile* fTreeFile;
    TTree* fTreeQA;
    
    int Run_Id;
    int Event_Id;
    int Event_Timestamp;
    float RPAngle;
    float ImpactParameter;
    
    bool T1;
    bool T2;
    bool T4;
    
    float triggersADC[nRawReco][nTriggers];
    float BPD_Position[nRawReco][nBPD][nBPDcomponents];
    
    //Main Vertex
    float Main_Vertex_X;
    float Main_Vertex_Y;
    float Main_Vertex_Z;

    //Primary Vertex
    float Primary_Vertex_X;
    float Primary_Vertex_Y;
    float Primary_Vertex_Z;
    bool Has_Primary_Vertex;
    int Primary_Vertex_Q;
    
//     float FitVertexX;
//     float FitVertexY;
//     float FitVertexZ;
//     int FitVertexQ;
    
    float PSD_module_Energy[nPSD_Modules];
    float PSD_section_Energy[nPSD_Modules][nPSD_Sections];
    float PSD_module_X[nPSD_Modules];
    float PSD_module_Y[nPSD_Modules];
    float PSD_module_Z[nPSD_Modules];

    float TPC_main_vtx_pT[nTPC_Tracks];
	float TPC_main_vtx_px[nTPC_Tracks];
	float TPC_main_vtx_py[nTPC_Tracks];
	float TPC_main_vtx_pz[nTPC_Tracks];
    float TPC_main_vtx_eta[nTPC_Tracks];
    float TPC_main_vtx_phi[nTPC_Tracks];
	float TPC_main_vtx_impact_point_X[nTPC_Tracks];
	float TPC_main_vtx_impact_point_Y[nTPC_Tracks];
	float TPC_main_vtx_impact_point_Z[nTPC_Tracks];
    float TPC_main_vtx_track_pT[nTPC_Tracks];
    float TPC_main_vtx_track_eta[nTPC_Tracks];
    float TPC_main_vtx_track_phi[nTPC_Tracks];
	float TPC_main_vtx_track_px[nTPC_Tracks];
	float TPC_main_vtx_track_py[nTPC_Tracks];
	float TPC_main_vtx_track_pz[nTPC_Tracks];
	float TPC_main_vtx_track_point_x[nTPC_Tracks];
	float TPC_main_vtx_track_point_y[nTPC_Tracks];
	float TPC_main_vtx_track_point_z[nTPC_Tracks];
	float TPC_main_vtx_track_first_point_x[nTPC_Tracks];
	float TPC_main_vtx_track_first_point_y[nTPC_Tracks];
	float TPC_main_vtx_track_first_point_z[nTPC_Tracks];
	float TPC_main_vtx_track_last_point_x[nTPC_Tracks];
	float TPC_main_vtx_track_last_point_y[nTPC_Tracks];
	float TPC_main_vtx_track_last_point_z[nTPC_Tracks];

	float TPC_main_vtx_ext_track_px[nTPC_Tracks];
	float TPC_main_vtx_ext_track_py[nTPC_Tracks];
	float TPC_main_vtx_ext_track_pz[nTPC_Tracks];
	float TPC_main_vtx_ext_first_track_px[nTPC_Tracks];
	float TPC_main_vtx_ext_first_track_py[nTPC_Tracks];
	float TPC_main_vtx_ext_first_track_pz[nTPC_Tracks];
	float TPC_main_vtx_ext_last_track_px[nTPC_Tracks];
	float TPC_main_vtx_ext_last_track_py[nTPC_Tracks];
	float TPC_main_vtx_ext_last_track_pz[nTPC_Tracks];
	float TPC_main_vtx_ext_point_x[nTPC_Tracks];
	float TPC_main_vtx_ext_point_y[nTPC_Tracks];
	float TPC_main_vtx_ext_point_z[nTPC_Tracks];
	float TPC_main_vtx_ext_first_point_x[nTPC_Tracks];
	float TPC_main_vtx_ext_first_point_y[nTPC_Tracks];
	float TPC_main_vtx_ext_first_point_z[nTPC_Tracks];
	float TPC_main_vtx_ext_last_point_x[nTPC_Tracks];
	float TPC_main_vtx_ext_last_point_y[nTPC_Tracks];
	float TPC_main_vtx_ext_last_point_z[nTPC_Tracks];

	float TPC_main_test_ext_track_px[nTPC_Tracks];
	float TPC_main_test_ext_track_py[nTPC_Tracks];
	float TPC_main_test_ext_track_pz[nTPC_Tracks];
	float TPC_main_test_ext_first_track_px[nTPC_Tracks];
	float TPC_main_test_ext_first_track_py[nTPC_Tracks];
	float TPC_main_test_ext_first_track_pz[nTPC_Tracks];
	float TPC_main_test_ext_last_track_px[nTPC_Tracks];
	float TPC_main_test_ext_last_track_py[nTPC_Tracks];
	float TPC_main_test_ext_last_track_pz[nTPC_Tracks];
	float TPC_main_test_ext_point_x[nTPC_Tracks];
	float TPC_main_test_ext_point_y[nTPC_Tracks];
	float TPC_main_test_ext_point_z[nTPC_Tracks];
	float TPC_main_test_ext_first_point_x[nTPC_Tracks];
	float TPC_main_test_ext_first_point_y[nTPC_Tracks];
	float TPC_main_test_ext_first_point_z[nTPC_Tracks];
	float TPC_main_test_ext_last_point_x[nTPC_Tracks];
	float TPC_main_test_ext_last_point_y[nTPC_Tracks];
	float TPC_main_test_ext_last_point_z[nTPC_Tracks];

	float TPC_main_PSD_ext_track_px[nTPC_Tracks];
	float TPC_main_PSD_ext_track_py[nTPC_Tracks];
	float TPC_main_PSD_ext_track_pz[nTPC_Tracks];
	float TPC_main_PSD_ext_point_x[nTPC_Tracks];
	float TPC_main_PSD_ext_point_y[nTPC_Tracks];
	float TPC_main_PSD_ext_point_z[nTPC_Tracks];
    
    int TPC_main_vtx_track_NofHits[nTPC_Tracks];
    int TPC_main_vtx_track_PidHypo[nTPC_Tracks];
    int TPC_main_vtx_track_Flag[nTPC_Tracks];
    float TPC_main_vtx_track_chi2[nTPC_Tracks];
    int TPC_main_vtx_track_ndf[nTPC_Tracks];
    
    int TPC_main_vtx_track_charge[nTPC_Tracks];
    
    int TPC_main_vtx_track_nClusters_Total[nTPC_Tracks];
    int TPC_main_vtx_track_nClusters_TPCV1[nTPC_Tracks];
    int TPC_main_vtx_track_nClusters_TPCV2[nTPC_Tracks];
    int TPC_main_vtx_track_nClusters_TPCVmain[nTPC_Tracks];
    int TPC_main_vtx_track_nClusters_TPCVgap[nTPC_Tracks];
    
    int TPC_main_vtx_track_nClustersPotential_Total[nTPC_Tracks];
    int TPC_main_vtx_track_nClustersPotential_TPCV1[nTPC_Tracks];
    int TPC_main_vtx_track_nClustersPotential_TPCV2[nTPC_Tracks];
    int TPC_main_vtx_track_nClustersPotential_TPCVmain[nTPC_Tracks];
    int TPC_main_vtx_track_nClustersPotential_TPCVgap[nTPC_Tracks];
    
    double TPC_main_vtx_track_EnergyClusters_Total[nTPC_Tracks];
    double TPC_main_vtx_track_EnergyClusters_TPCV1[nTPC_Tracks];
    double TPC_main_vtx_track_EnergyClusters_TPCV2[nTPC_Tracks];
    double TPC_main_vtx_track_EnergyClusters_TPCVmain[nTPC_Tracks];
    double TPC_main_vtx_track_EnergyClusters_TPCVgap[nTPC_Tracks];
        

    float TPC_primary_vtx_pT[nTPC_Tracks];
    float TPC_primary_vtx_eta[nTPC_Tracks];
    float TPC_primary_vtx_phi[nTPC_Tracks];
    float TPC_primary_vtx_track_pT[nTPC_Tracks];
    float TPC_primary_vtx_track_eta[nTPC_Tracks];
    float TPC_primary_vtx_track_phi[nTPC_Tracks];
    
    int TPC_primary_vtx_track_NofHits[nTPC_Tracks];
    int TPC_primary_vtx_track_PidHypo[nTPC_Tracks];
    int TPC_primary_vtx_track_Flag[nTPC_Tracks];
    float TPC_primary_vtx_track_ChiSq[nTPC_Tracks];
    int TPC_primary_vtx_track_NDF[nTPC_Tracks];
    
    float TPC_primary_vtx_impact_point_X[nTPC_Tracks];
    float TPC_primary_vtx_impact_point_Y[nTPC_Tracks];
    float TPC_primary_vtx_impact_point_Z[nTPC_Tracks];
    
    int TPC_primary_vtx_track_charge[nTPC_Tracks];
    
    int TPC_primary_vtx_track_nClusters_Total[nTPC_Tracks];
    int TPC_primary_vtx_track_nClusters_TPCV1[nTPC_Tracks];
    int TPC_primary_vtx_track_nClusters_TPCV2[nTPC_Tracks];
    int TPC_primary_vtx_track_nClusters_TPCVmain[nTPC_Tracks];
    int TPC_primary_vtx_track_nClusters_TPCVgap[nTPC_Tracks];
    
    int TPC_primary_vtx_track_nClustersPotential_Total[nTPC_Tracks];
    int TPC_primary_vtx_track_nClustersPotential_TPCV1[nTPC_Tracks];
    int TPC_primary_vtx_track_nClustersPotential_TPCV2[nTPC_Tracks];
    int TPC_primary_vtx_track_nClustersPotential_TPCVmain[nTPC_Tracks];
    int TPC_primary_vtx_track_nClustersPotential_TPCVgap[nTPC_Tracks];
    
    double TPC_primary_vtx_track_EnergyClusters_Total[nTPC_Tracks];
    double TPC_primary_vtx_track_EnergyClusters_TPCV1[nTPC_Tracks];
    double TPC_primary_vtx_track_EnergyClusters_TPCV2[nTPC_Tracks];
    double TPC_primary_vtx_track_EnergyClusters_TPCVmain[nTPC_Tracks];
    double TPC_primary_vtx_track_EnergyClusters_TPCVgap[nTPC_Tracks];

	bool isTriggers_Combined[nRawReco][nTriggers_Combined];

    float WFA_TimeStructure[nTriggers_Simple][nMaxWFAsignals];
    int    WFA_NumberOfSignalHits[nTriggers_Simple];
    
//     float TPC_track_pT[nTPC_Tracks];
//     float TPC_track_eta[nTPC_Tracks];
//     float TPC_track_phi[nTPC_Tracks];
//     
//     int TPC_track_NofHits[nTPC_Tracks];
//     int TPC_track_PidHypo[nTPC_Tracks];
//     int TPC_track_Flag[nTPC_Tracks];
//     double TPC_track_ChiSq[nTPC_Tracks];
//     int TPC_track_NDF[nTPC_Tracks];
//     
//     float TPC_track_DCAtoVertex_X[nTPC_Tracks];
//     float TPC_track_DCAtoVertex_Y[nTPC_Tracks];
//     
//     int TPC_track_nClusters_Total[nTPC_Tracks];
//     int TPC_track_nClusters_TPCV1[nTPC_Tracks];
//     int TPC_track_nClusters_TPCV2[nTPC_Tracks];
//     int TPC_track_nClusters_TPCVmain[nTPC_Tracks];
//     int TPC_track_nClusters_TPCVgap[nTPC_Tracks];
//     
//     int TPC_track_nClustersPotential_Total[nTPC_Tracks];
//     int TPC_track_nClustersPotential_TPCV1[nTPC_Tracks];
//     int TPC_track_nClustersPotential_TPCV2[nTPC_Tracks];
//     int TPC_track_nClustersPotential_TPCVmain[nTPC_Tracks];
//     int TPC_track_nClustersPotential_TPCVgap[nTPC_Tracks];
    
    ClassDefNV(DataTreeConverterNA61, 1)
//     ClassDef(DataTreeConverterNA61,1);
};

#endif
