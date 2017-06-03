#include "DataTreeConverterNA61.h"
#include <iostream>
#include <fstream>
#include <vector>
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

#include "TDirectory.h"

//=================================================================> MAIN <===============================================================
DataTreeConverterNA61::DataTreeConverterNA61()
:
  nEvents (1),
  fCurEvent (0)
{
      DTEvent = new DataTreeEvent();
}
//--------------------------------------------------------------------------------------------------
DataTreeConverterNA61::~DataTreeConverterNA61()
{
    
}
//--------------------------------------------------------------------------------------------------
void DataTreeConverterNA61::Body()
{
//     cout << "Init" << endl;
    Init();
    nRealEvents = fChain->GetEntries();
//     cout << nRealEvents << " " << nEvents << endl;//TEST
    if (fChain->GetEntries() > nEvents) nRealEvents = nEvents;
//     cout << nRealEvents << " " << nEvents << endl;//TEST
//     cout << "Exec" << endl;
    for (int i=0;i<nRealEvents;i++)
    {
	Exec(i);
    }
//     cout << "Finish" << endl;
    Finish();
}
//=================================================================> INIT <===============================================================
void DataTreeConverterNA61::Init()
{
//     cout << "Chain_Init" << endl;
    Chain_Init();
//     cout << "OutputTree_Init" << endl;
    OutputTree_Init();
//     cout << "DataTreeEvent_Init" << endl;
    DataTreeEvent_Init();
}
//--------------------------------------------------------------------------------------------------
void DataTreeConverterNA61::DataTreeEvent_Init()
{
    for (int i=0;i<nPSD_Modules;i++)
    {
	DTEvent -> AddPSDModule(10);
    }
    TString sTriggerNames[8] = {"S1","S2","S3","V1","PSD","T1","T2","T4"};	
    for (int i=0;i<8;i++)
    {
	DTEvent -> AddTrigger(sTriggerNames[i]);
    }
    for (int i=0;i<nBPD;i++)
    {
	DTEvent -> AddBPD();
    }
	for (int i=0;i<nTriggers_Simple;i++)
	{
	DTEvent -> AddWFA();
	}
}
//--------------------------------------------------------------------------------------------------
void DataTreeConverterNA61::OutputTree_Init()
{
    fTreeFile = new TFile(sOutputFileName, "RECREATE");
    fTreeFile -> cd();
    fTreeQA = new TTree("fDataTree","fDataTree");
    fTreeQA -> SetMaxTreeSize(90000000);
    
    fTreeQA -> Branch("DTEvent", &DTEvent, 256000, 3);    
}
//--------------------------------------------------------------------------------------------------
void DataTreeConverterNA61::Chain_Init()
{
    fChain = new TChain("fTreeQA","fTreeQA");
    fChain -> Add(sInputFileName);
    
    fChain->SetBranchAddress("Run_Id",&Run_Id);
    fChain->SetBranchAddress("Event_Id",&Event_Id);
//     fChain->SetBranchAddress("RPAngle",&RPAngle);
//     fChain->SetBranchAddress("ImpactParameter",&ImpactParameter);
    
    fChain->SetBranchAddress("Event_Timestamp",&Event_Timestamp);
    
//     fChain->SetBranchAddress("T1",&T1);
//     fChain->SetBranchAddress("T2",&T2);
//     fChain->SetBranchAddress("T4",&T4);
    
    fChain->SetBranchAddress("triggersADC",&triggersADC);
    
    fChain->SetBranchAddress("BPD_Position",&BPD_Position);
    
    fChain->SetBranchAddress("Main_Vertex_X",&Main_Vertex_X);
    fChain->SetBranchAddress("Main_Vertex_Y",&Main_Vertex_Y);
    fChain->SetBranchAddress("Main_Vertex_Z",&Main_Vertex_Z);
    fChain->SetBranchAddress("Primary_Vertex_X",&Primary_Vertex_X);
    fChain->SetBranchAddress("Primary_Vertex_Y",&Primary_Vertex_Y);
    fChain->SetBranchAddress("Primary_Vertex_Z",&Primary_Vertex_Z);
    fChain->SetBranchAddress("Has_Primary_Vertex",&Has_Primary_Vertex);
    fChain->SetBranchAddress("Primary_Vertex_Q",&Primary_Vertex_Q);
    
//     fChain->SetBranchAddress("FitVertexX",&FitVertexX);
//     fChain->SetBranchAddress("FitVertexY",&FitVertexY);
//     fChain->SetBranchAddress("FitVertexZ",&FitVertexZ);
//     fChain->SetBranchAddress("FitVertexQ",&FitVertexQ);
    
    fChain->SetBranchAddress("PSD_module_Energy",&PSD_module_Energy);
    fChain->SetBranchAddress("PSD_section_Energy",&PSD_section_Energy);
    fChain->SetBranchAddress("PSD_module_X",&PSD_module_X);
    fChain->SetBranchAddress("PSD_module_Y",&PSD_module_Y);
    fChain->SetBranchAddress("PSD_module_Z",&PSD_module_Z);

    fChain->SetBranchAddress("TPC_main_vtx_pT",&TPC_main_vtx_pT);
	fChain->SetBranchAddress("TPC_main_vtx_px",&TPC_main_vtx_px);
	fChain->SetBranchAddress("TPC_main_vtx_py",&TPC_main_vtx_py);
	fChain->SetBranchAddress("TPC_main_vtx_pz",&TPC_main_vtx_pz);
    fChain->SetBranchAddress("TPC_main_vtx_eta",&TPC_main_vtx_eta);
    fChain->SetBranchAddress("TPC_main_vtx_phi",&TPC_main_vtx_phi);
    fChain->SetBranchAddress("TPC_main_vtx_track_pT",&TPC_main_vtx_track_pT);
    fChain->SetBranchAddress("TPC_main_vtx_track_eta",&TPC_main_vtx_track_eta);
    fChain->SetBranchAddress("TPC_main_vtx_track_phi",&TPC_main_vtx_track_phi);
	fChain->SetBranchAddress("TPC_main_vtx_track_px",&TPC_main_vtx_track_px);
	fChain->SetBranchAddress("TPC_main_vtx_track_py",&TPC_main_vtx_track_py);
	fChain->SetBranchAddress("TPC_main_vtx_track_pz",&TPC_main_vtx_track_pz);
	fChain->SetBranchAddress("TPC_main_vtx_track_point_x",&TPC_main_vtx_track_point_x);
	fChain->SetBranchAddress("TPC_main_vtx_track_point_y",&TPC_main_vtx_track_point_y);
	fChain->SetBranchAddress("TPC_main_vtx_track_point_z",&TPC_main_vtx_track_point_z);
	fChain->SetBranchAddress("TPC_main_vtx_track_first_point_x",&TPC_main_vtx_track_first_point_x);
	fChain->SetBranchAddress("TPC_main_vtx_track_first_point_y",&TPC_main_vtx_track_first_point_y);
	fChain->SetBranchAddress("TPC_main_vtx_track_first_point_z",&TPC_main_vtx_track_first_point_z);
	fChain->SetBranchAddress("TPC_main_vtx_track_last_point_x",&TPC_main_vtx_track_last_point_x);
	fChain->SetBranchAddress("TPC_main_vtx_track_last_point_y",&TPC_main_vtx_track_last_point_y);
	fChain->SetBranchAddress("TPC_main_vtx_track_last_point_z",&TPC_main_vtx_track_last_point_z);
    fChain->SetBranchAddress("TPC_main_vtx_impact_point_X",&TPC_main_vtx_impact_point_X);
    fChain->SetBranchAddress("TPC_main_vtx_impact_point_Y",&TPC_main_vtx_impact_point_Y);
    fChain->SetBranchAddress("TPC_main_vtx_impact_point_Z",&TPC_main_vtx_impact_point_Z);
   
	fChain->SetBranchAddress("TPC_main_vtx_ext_track_px", &TPC_main_vtx_ext_track_px);
    fChain->SetBranchAddress("TPC_main_vtx_ext_track_py", &TPC_main_vtx_ext_track_py);
    fChain->SetBranchAddress("TPC_main_vtx_ext_track_pz", &TPC_main_vtx_ext_track_pz);
    fChain->SetBranchAddress("TPC_main_vtx_ext_first_track_px", &TPC_main_vtx_ext_first_track_px);
    fChain->SetBranchAddress("TPC_main_vtx_ext_first_track_py", &TPC_main_vtx_ext_first_track_py);
    fChain->SetBranchAddress("TPC_main_vtx_ext_first_track_pz", &TPC_main_vtx_ext_first_track_pz);
    fChain->SetBranchAddress("TPC_main_vtx_ext_last_track_px", &TPC_main_vtx_ext_last_track_px);
    fChain->SetBranchAddress("TPC_main_vtx_ext_last_track_py", &TPC_main_vtx_ext_last_track_py);
    fChain->SetBranchAddress("TPC_main_vtx_ext_last_track_pz", &TPC_main_vtx_ext_last_track_pz);
    fChain->SetBranchAddress("TPC_main_vtx_ext_point_x", &TPC_main_vtx_ext_point_x);
    fChain->SetBranchAddress("TPC_main_vtx_ext_point_y", &TPC_main_vtx_ext_point_y);
    fChain->SetBranchAddress("TPC_main_vtx_ext_point_z", &TPC_main_vtx_ext_point_z);
    fChain->SetBranchAddress("TPC_main_vtx_ext_first_point_x", &TPC_main_vtx_ext_first_point_x);
    fChain->SetBranchAddress("TPC_main_vtx_ext_first_point_y", &TPC_main_vtx_ext_first_point_y);
    fChain->SetBranchAddress("TPC_main_vtx_ext_first_point_z", &TPC_main_vtx_ext_first_point_z);
    fChain->SetBranchAddress("TPC_main_vtx_ext_last_point_x", &TPC_main_vtx_ext_last_point_x);
    fChain->SetBranchAddress("TPC_main_vtx_ext_last_point_y", &TPC_main_vtx_ext_last_point_y);
    fChain->SetBranchAddress("TPC_main_vtx_ext_last_point_z", &TPC_main_vtx_ext_last_point_z);
    fChain->SetBranchAddress("TPC_main_test_ext_track_px", &TPC_main_test_ext_track_px);
    fChain->SetBranchAddress("TPC_main_test_ext_track_py", &TPC_main_test_ext_track_py);
    fChain->SetBranchAddress("TPC_main_test_ext_track_pz", &TPC_main_test_ext_track_pz);
    fChain->SetBranchAddress("TPC_main_test_ext_first_track_px", &TPC_main_test_ext_first_track_px);
    fChain->SetBranchAddress("TPC_main_test_ext_first_track_py", &TPC_main_test_ext_first_track_py);
    fChain->SetBranchAddress("TPC_main_test_ext_first_track_pz", &TPC_main_test_ext_first_track_pz);
    fChain->SetBranchAddress("TPC_main_test_ext_last_track_px", &TPC_main_test_ext_last_track_px);
    fChain->SetBranchAddress("TPC_main_test_ext_last_track_py", &TPC_main_test_ext_last_track_py);
    fChain->SetBranchAddress("TPC_main_test_ext_last_track_pz", &TPC_main_test_ext_last_track_pz);
    fChain->SetBranchAddress("TPC_main_test_ext_point_x", &TPC_main_test_ext_point_x);
    fChain->SetBranchAddress("TPC_main_test_ext_point_y", &TPC_main_test_ext_point_y);
    fChain->SetBranchAddress("TPC_main_test_ext_point_z", &TPC_main_test_ext_point_z);
    fChain->SetBranchAddress("TPC_main_test_ext_first_point_x", &TPC_main_test_ext_first_point_x);
    fChain->SetBranchAddress("TPC_main_test_ext_first_point_y", &TPC_main_test_ext_first_point_y);
    fChain->SetBranchAddress("TPC_main_test_ext_first_point_z", &TPC_main_test_ext_first_point_z);
    fChain->SetBranchAddress("TPC_main_test_ext_last_point_x", &TPC_main_test_ext_last_point_x);
    fChain->SetBranchAddress("TPC_main_test_ext_last_point_y", &TPC_main_test_ext_last_point_y);
    fChain->SetBranchAddress("TPC_main_test_ext_last_point_z", &TPC_main_test_ext_last_point_z);

    fChain->SetBranchAddress("TPC_main_PSD_ext_track_px", &TPC_main_PSD_ext_track_px);
    fChain->SetBranchAddress("TPC_main_PSD_ext_track_py", &TPC_main_PSD_ext_track_py);
    fChain->SetBranchAddress("TPC_main_PSD_ext_track_pz", &TPC_main_PSD_ext_track_pz);
    fChain->SetBranchAddress("TPC_main_PSD_ext_point_x", &TPC_main_PSD_ext_point_x);
    fChain->SetBranchAddress("TPC_main_PSD_ext_point_y", &TPC_main_PSD_ext_point_y);
    fChain->SetBranchAddress("TPC_main_PSD_ext_point_z", &TPC_main_PSD_ext_point_z);

	fChain->SetBranchAddress("TPC_main_vtx_track_charge", &TPC_main_vtx_track_charge);
    fChain->SetBranchAddress("TPC_main_vtx_track_nClusters_Total",&TPC_main_vtx_track_nClusters_Total);
    fChain->SetBranchAddress("TPC_main_vtx_track_nClusters_TPCV1",&TPC_main_vtx_track_nClusters_TPCV1);
    fChain->SetBranchAddress("TPC_main_vtx_track_nClusters_TPCV2",&TPC_main_vtx_track_nClusters_TPCV2);
    fChain->SetBranchAddress("TPC_main_vtx_track_nClusters_TPCVmain",&TPC_main_vtx_track_nClusters_TPCVmain);
	fChain->SetBranchAddress("TPC_main_vtx_track_chi2",&TPC_main_vtx_track_chi2);
	fChain->SetBranchAddress("TPC_main_vtx_track_ndf",&TPC_main_vtx_track_ndf);
//     fChain->SetBranchAddress("TPC_main_vtx_track_nClusters_TPCVgap",&TPC_main_vtx_track_nClusters_TPCVgap);
    
//     fChain->SetBranchAddress("TPC_main_vtx_track_nClustersPotential_Total",&TPC_main_vtx_track_nClustersPotential_Total);
//     fChain->SetBranchAddress("TPC_main_vtx_track_nClustersPotential_TPCV1",&TPC_main_vtx_track_nClustersPotential_TPCV1);
//     fChain->SetBranchAddress("TPC_main_vtx_track_nClustersPotential_TPCV2",&TPC_main_vtx_track_nClustersPotential_TPCV2);
//     fChain->SetBranchAddress("TPC_main_vtx_track_nClustersPotential_TPCVmain",&TPC_main_vtx_track_nClustersPotential_TPCVmain);
//     fChain->SetBranchAddress("TPC_main_vtx_track_nClustersPotential_TPCVgap",&TPC_main_vtx_track_nClustersPotential_TPCVgap);
    
    fChain->SetBranchAddress("TPC_main_vtx_track_EnergyClusters_Total",&TPC_main_vtx_track_EnergyClusters_Total);
    fChain->SetBranchAddress("TPC_main_vtx_track_EnergyClusters_TPCV1",&TPC_main_vtx_track_EnergyClusters_TPCV1);
    fChain->SetBranchAddress("TPC_main_vtx_track_EnergyClusters_TPCV2",&TPC_main_vtx_track_EnergyClusters_TPCV2);
    fChain->SetBranchAddress("TPC_main_vtx_track_EnergyClusters_TPCVmain",&TPC_main_vtx_track_EnergyClusters_TPCVmain);
//     fChain->SetBranchAddress("TPC_main_vtx_track_EnergyClusters_TPCVgap",&TPC_main_vtx_track_EnergyClusters_TPCVgap);

    fChain->SetBranchAddress("TPC_primary_vtx_pT",&TPC_primary_vtx_pT);
    fChain->SetBranchAddress("TPC_primary_vtx_eta",&TPC_primary_vtx_eta);
    fChain->SetBranchAddress("TPC_primary_vtx_phi",&TPC_primary_vtx_phi);
    fChain->SetBranchAddress("TPC_primary_vtx_track_pT",&TPC_primary_vtx_track_pT);
    fChain->SetBranchAddress("TPC_primary_vtx_track_eta",&TPC_primary_vtx_track_eta);
    fChain->SetBranchAddress("TPC_primary_vtx_track_phi",&TPC_primary_vtx_track_phi);
    fChain->SetBranchAddress("TPC_primary_vtx_impact_point_X",&TPC_primary_vtx_impact_point_X);
    fChain->SetBranchAddress("TPC_primary_vtx_impact_point_Y",&TPC_primary_vtx_impact_point_Y);
    fChain->SetBranchAddress("TPC_primary_vtx_impact_point_Z",&TPC_primary_vtx_impact_point_Z);
    
    fChain->SetBranchAddress("TPC_primary_vtx_track_nClusters_Total",&TPC_primary_vtx_track_nClusters_Total);
    fChain->SetBranchAddress("TPC_primary_vtx_track_nClusters_TPCV1",&TPC_primary_vtx_track_nClusters_TPCV1);
    fChain->SetBranchAddress("TPC_primary_vtx_track_nClusters_TPCV2",&TPC_primary_vtx_track_nClusters_TPCV2);
    fChain->SetBranchAddress("TPC_primary_vtx_track_nClusters_TPCVmain",&TPC_primary_vtx_track_nClusters_TPCVmain);
//     fChain->SetBranchAddress("TPC_primary_vtx_track_nClusters_TPCVgap",&TPC_primary_vtx_track_nClusters_TPCVgap);
    
//     fChain->SetBranchAddress("TPC_primary_vtx_track_nClustersPotential_Total",&TPC_primary_vtx_track_nClustersPotential_Total);
//     fChain->SetBranchAddress("TPC_primary_vtx_track_nClustersPotential_TPCV1",&TPC_primary_vtx_track_nClustersPotential_TPCV1);
//     fChain->SetBranchAddress("TPC_primary_vtx_track_nClustersPotential_TPCV2",&TPC_primary_vtx_track_nClustersPotential_TPCV2);
//     fChain->SetBranchAddress("TPC_primary_vtx_track_nClustersPotential_TPCVmain",&TPC_primary_vtx_track_nClustersPotential_TPCVmain);
//     fChain->SetBranchAddress("TPC_primary_vtx_track_nClustersPotential_TPCVgap",&TPC_primary_vtx_track_nClustersPotential_TPCVgap);
    
    fChain->SetBranchAddress("TPC_primary_vtx_track_EnergyClusters_Total",&TPC_primary_vtx_track_EnergyClusters_Total);
    fChain->SetBranchAddress("TPC_primary_vtx_track_EnergyClusters_TPCV1",&TPC_primary_vtx_track_EnergyClusters_TPCV1);
    fChain->SetBranchAddress("TPC_primary_vtx_track_EnergyClusters_TPCV2",&TPC_primary_vtx_track_EnergyClusters_TPCV2);
    fChain->SetBranchAddress("TPC_primary_vtx_track_EnergyClusters_TPCVmain",&TPC_primary_vtx_track_EnergyClusters_TPCVmain);
//     fChain->SetBranchAddress("TPC_primary_vtx_track_EnergyClusters_TPCVgap",&TPC_primary_vtx_track_EnergyClusters_TPCVgap);


    
//     fChain->SetBranchAddress("TPC_track_pT",&TPC_track_pT);
//     fChain->SetBranchAddress("TPC_track_eta",&TPC_track_eta);
//     fChain->SetBranchAddress("TPC_track_phi",&TPC_track_phi);
//     fChain->SetBranchAddress("TPC_track_DCAtoVertex_X",&TPC_track_DCAtoVertex_X);
//     fChain->SetBranchAddress("TPC_track_DCAtoVertex_Y",&TPC_track_DCAtoVertex_Y);
    
//     fChain->SetBranchAddress("TPC_track_nClusters_Total",&TPC_track_nClusters_Total);
//     fChain->SetBranchAddress("TPC_track_nClusters_TPCV1",&TPC_track_nClusters_TPCV1);
//     fChain->SetBranchAddress("TPC_track_nClusters_TPCV2",&TPC_track_nClusters_TPCV2);
//     fChain->SetBranchAddress("TPC_track_nClusters_TPCVmain",&TPC_track_nClusters_TPCVmain);
//     fChain->SetBranchAddress("TPC_track_nClusters_TPCVgap",&TPC_track_nClusters_TPCVgap);
//     
//     fChain->SetBranchAddress("TPC_track_nClustersPotential_Total",&TPC_track_nClustersPotential_Total);
//     fChain->SetBranchAddress("TPC_track_nClustersPotential_TPCV1",&TPC_track_nClustersPotential_TPCV1);
//     fChain->SetBranchAddress("TPC_track_nClustersPotential_TPCV2",&TPC_track_nClustersPotential_TPCV2);
//     fChain->SetBranchAddress("TPC_track_nClustersPotential_TPCVmain",&TPC_track_nClustersPotential_TPCVmain);
//     fChain->SetBranchAddress("TPC_track_nClustersPotential_TPCVgap",&TPC_track_nClustersPotential_TPCVgap);
//     
//     fChain->SetBranchAddress("TPC_track_EnergyClusters_Total",&TPC_track_EnergyClusters_Total);
//     fChain->SetBranchAddress("TPC_track_EnergyClusters_TPCV1",&TPC_track_EnergyClusters_TPCV1);
//     fChain->SetBranchAddress("TPC_track_EnergyClusters_TPCV2",&TPC_track_EnergyClusters_TPCV2);
//     fChain->SetBranchAddress("TPC_track_EnergyClusters_TPCVmain",&TPC_track_EnergyClusters_TPCVmain);
//     fChain->SetBranchAddress("TPC_track_EnergyClusters_TPCVgap",&TPC_track_EnergyClusters_TPCVgap);

	fChain->SetBranchAddress("isTriggers_Combined",&isTriggers_Combined);

	fChain->SetBranchAddress("WFA_TimeStructure",&WFA_TimeStructure);
	fChain->SetBranchAddress("WFA_NumberOfSignalHits",&WFA_NumberOfSignalHits);
}

//=================================================================> EXEC <===============================================================
void DataTreeConverterNA61::Exec(int event)
{
    fChain -> GetEntry(event);
    Read_Event();
    DTEvent -> Process();
    fTreeQA -> Fill();
    Clear_Event();
}
//--------------------------------------------------------------------------------------------------
void DataTreeConverterNA61::Clear_Event()
{
    DTEvent -> ClearEvent();
}
//--------------------------------------------------------------------------------------------------
void DataTreeConverterNA61::Read_Event()
{
    DTEvent -> SetRunId(Run_Id);
    DTEvent -> SetEventId(Event_Id);
    DTEvent -> SetRPAngle(RPAngle);
    DTEvent -> SetImpactParameter(ImpactParameter);
    DTEvent -> SetEventTimestamp(Event_Timestamp);
    
    DTEvent -> GetTrigger(0) -> SetSignal(triggersADC[0][0]);	//S1
    DTEvent -> GetTrigger(1) -> SetSignal(triggersADC[0][1]);	//S2
    DTEvent -> GetTrigger(2) -> SetSignal(triggersADC[0][4]);	//V1p -> S3
    DTEvent -> GetTrigger(3) -> SetSignal(triggersADC[0][3]);	//V1
    DTEvent -> GetTrigger(4) -> SetSignal(triggersADC[0][5]);	//PSD
    DTEvent -> GetTrigger(5) -> SetBool(isTriggers_Combined[0][0]); //T1
    DTEvent -> GetTrigger(6) -> SetBool(isTriggers_Combined[0][1]); //T2
    DTEvent -> GetTrigger(7) -> SetBool(isTriggers_Combined[0][3]); //T4
    
    for (int i=0;i<nBPD;i++)
    {
	DTEvent -> GetBPD(i) -> SetPosition(BPD_Position[1][i][0],BPD_Position[1][i][1],BPD_Position[1][i][2]);
    }

	for (int i=0;i<nTriggers_Simple;i++)
    {
		DTEvent -> GetWFA(i) -> SetNHits(i,WFA_NumberOfSignalHits[i]);
		for (int j=0;j<nMaxWFAsignals;j++){
			DTEvent -> GetWFA(i) -> SetTime(i,j,WFA_TimeStructure[i][j]);
			// TESTING WFA signal
			if (i == 0 && j == 10) std::cout << "WFA from TreeNA61: " << WFA_TimeStructure[i][j] << "; WFA from DataTree: " << DTEvent->GetWFA(i)->GetTime(i,j) << std::endl;;
		}
    }
    
    DTEvent -> SetVertexPosition(Main_Vertex_X,Main_Vertex_Y,Main_Vertex_Z,0);
    DTEvent -> SetVertexPosition(Primary_Vertex_X,Primary_Vertex_Y,Primary_Vertex_Z,1);
    DTEvent -> SetHasVertex(Has_Primary_Vertex,1);
//     DTEvent -> SetVertexPosition(FitVertexX,FitVertexY,FitVertexZ);
    DTEvent -> SetVertexQuality(Primary_Vertex_Q,1);
    
    for (int i=0;i<nPSD_Modules;i++)
    {
	DataTreePSDModule* mod = DTEvent -> GetPSDModule(i);
	mod -> SetPosition(PSD_module_X[i],PSD_module_Y[i],PSD_module_Z[i]);
	for (int j=0;j<nPSD_Sections;j++)
	{
	    DataTreePSDSection* sec = mod -> GetSection(j);
	    sec -> AddEnergy(PSD_section_Energy[i][j]);
	}
    }
    for (int i=0;i<nTPC_Tracks;i++)
    {
	if (TPC_main_vtx_pT[i] < 0) break;
	DTEvent -> AddTrack();
	DataTreeTrack* track = DTEvent -> GetLastTrack();
	
	track -> SetType(0);//Main vertex track
	track -> SetPt(0,TPC_main_vtx_pT[i]);
	track -> SetPx(0,TPC_main_vtx_px[i]);
	track -> SetPy(0,TPC_main_vtx_py[i]);
	track -> SetPz(0,TPC_main_vtx_pz[i]);
	track -> SetEta(0,TPC_main_vtx_eta[i]);
	track -> SetPhi(0,TPC_main_vtx_phi[i]);
	track -> SetDCA(0,TPC_main_vtx_impact_point_X[i]-Main_Vertex_X,TPC_main_vtx_impact_point_Y[i]-Main_Vertex_Y,TPC_main_vtx_impact_point_Z[i]-Main_Vertex_Z);

	track -> SetTrackPointComponent(0,0,TPC_main_vtx_impact_point_X[i]);
	track -> SetTrackPointComponent(0,1,TPC_main_vtx_impact_point_Y[i]);
	track -> SetTrackPointComponent(0,2,TPC_main_vtx_impact_point_Z[i]);
	track -> SetTrackPointComponent(1,0,TPC_main_vtx_track_point_x[i]);
	track -> SetTrackPointComponent(1,1,TPC_main_vtx_track_point_y[i]);
	track -> SetTrackPointComponent(1,2,TPC_main_vtx_track_point_z[i]);
	track -> SetTrackPointComponent(2,0,TPC_main_vtx_track_last_point_x[i]);
	track -> SetTrackPointComponent(2,1,TPC_main_vtx_track_last_point_y[i]);
	track -> SetTrackPointComponent(2,2,TPC_main_vtx_track_last_point_z[i]);

	track -> SetExtPointComponent(0,0,TPC_main_vtx_ext_point_x[i]);
	track -> SetExtPointComponent(0,1,TPC_main_vtx_ext_point_y[i]);
	track -> SetExtPointComponent(0,2,TPC_main_vtx_ext_point_z[i]);
	track -> SetExtPointComponent(1,0,TPC_main_vtx_ext_last_point_x[i]);
	track -> SetExtPointComponent(1,1,TPC_main_vtx_ext_last_point_y[i]);
	track -> SetExtPointComponent(1,2,TPC_main_vtx_ext_last_point_z[i]);	
	track -> SetExtPointComponent(2,0,TPC_main_test_ext_point_x[i]);
	track -> SetExtPointComponent(2,1,TPC_main_test_ext_point_y[i]);
	track -> SetExtPointComponent(2,2,TPC_main_test_ext_point_z[i]);
	track -> SetExtPointComponent(3,0,TPC_main_PSD_ext_point_x[i]);
	track -> SetExtPointComponent(3,1,TPC_main_PSD_ext_point_y[i]);
	track -> SetExtPointComponent(3,2,TPC_main_PSD_ext_point_z[i]);

	track -> SetExtPComponent(0,0,TPC_main_vtx_ext_track_px[i]);
	track -> SetExtPComponent(0,1,TPC_main_vtx_ext_track_py[i]);
	track -> SetExtPComponent(0,2,TPC_main_vtx_ext_track_pz[i]);
	track -> SetExtPComponent(1,0,TPC_main_vtx_ext_last_track_px[i]);
	track -> SetExtPComponent(1,1,TPC_main_vtx_ext_last_track_py[i]);
	track -> SetExtPComponent(1,2,TPC_main_vtx_ext_last_track_pz[i]);
	track -> SetExtPComponent(2,0,TPC_main_test_ext_track_px[i]);
	track -> SetExtPComponent(2,1,TPC_main_test_ext_track_py[i]);
	track -> SetExtPComponent(2,2,TPC_main_test_ext_track_pz[i]);
	track -> SetExtPComponent(3,0,TPC_main_PSD_ext_track_px[i]);
	track -> SetExtPComponent(3,1,TPC_main_PSD_ext_track_py[i]);
	track -> SetExtPComponent(3,2,TPC_main_PSD_ext_track_pz[i]);

	track-> SetNofHits(0,TPC_main_vtx_track_nClusters_Total[i],0);
	track-> SetNofHits(0,TPC_main_vtx_track_nClusters_TPCV1[i],1);
	track-> SetNofHits(0,TPC_main_vtx_track_nClusters_TPCV2[i],2);
	track-> SetNofHits(0,TPC_main_vtx_track_nClusters_TPCVmain[i],3);
	track-> SetNofHits(0,TPC_main_vtx_track_nClusters_TPCVgap[i],4);
	
	track-> SetNofHitsPotential(0,TPC_main_vtx_track_nClustersPotential_Total[i],0);
	track-> SetNofHitsPotential(0,TPC_main_vtx_track_nClustersPotential_TPCV1[i],1);
	track-> SetNofHitsPotential(0,TPC_main_vtx_track_nClustersPotential_TPCV2[i],2);
	track-> SetNofHitsPotential(0,TPC_main_vtx_track_nClustersPotential_TPCVmain[i],3);
	track-> SetNofHitsPotential(0,TPC_main_vtx_track_nClustersPotential_TPCVgap[i],4);
	
	track-> SetdEdx(0,TPC_main_vtx_track_EnergyClusters_Total[i],0);
	track-> SetdEdx(0,TPC_main_vtx_track_EnergyClusters_TPCV1[i],1);
	track-> SetdEdx(0,TPC_main_vtx_track_EnergyClusters_TPCV2[i],2);
	track-> SetdEdx(0,TPC_main_vtx_track_EnergyClusters_TPCVmain[i],3);
	track-> SetdEdx(0,TPC_main_vtx_track_EnergyClusters_TPCVgap[i],4);
	
	track-> SetCharge(0,TPC_main_vtx_track_charge[i]);
	track-> SetChiSq(0,TPC_main_vtx_track_chi2[i]);
	track-> SetNDF(0,TPC_main_vtx_track_ndf[i]);
	
	track -> SetPt(1,TPC_main_vtx_track_pT[i]);
	track -> SetEta(1,TPC_main_vtx_track_eta[i]);
	track -> SetPhi(1,TPC_main_vtx_track_phi[i]);
	track -> SetDCA(1,TPC_main_vtx_impact_point_X[i]-Primary_Vertex_X,TPC_main_vtx_impact_point_Y[i]-Primary_Vertex_Y,TPC_main_vtx_impact_point_Z[i]-Primary_Vertex_Z);
	
	track-> SetNofHits(1,TPC_main_vtx_track_nClusters_Total[i],0);
	track-> SetNofHits(1,TPC_main_vtx_track_nClusters_TPCV1[i],1);
	track-> SetNofHits(1,TPC_main_vtx_track_nClusters_TPCV2[i],2);
	track-> SetNofHits(1,TPC_main_vtx_track_nClusters_TPCVmain[i],3);
	track-> SetNofHits(1,TPC_main_vtx_track_nClusters_TPCVgap[i],4);
	
	track-> SetNofHitsPotential(1,TPC_main_vtx_track_nClustersPotential_Total[i],0);
	track-> SetNofHitsPotential(1,TPC_main_vtx_track_nClustersPotential_TPCV1[i],1);
	track-> SetNofHitsPotential(1,TPC_main_vtx_track_nClustersPotential_TPCV2[i],2);
	track-> SetNofHitsPotential(1,TPC_main_vtx_track_nClustersPotential_TPCVmain[i],3);
	track-> SetNofHitsPotential(1,TPC_main_vtx_track_nClustersPotential_TPCVgap[i],4);
	
	track-> SetdEdx(1,TPC_main_vtx_track_EnergyClusters_Total[i],0);
	track-> SetdEdx(1,TPC_main_vtx_track_EnergyClusters_TPCV1[i],1);
	track-> SetdEdx(1,TPC_main_vtx_track_EnergyClusters_TPCV2[i],2);
	track-> SetdEdx(1,TPC_main_vtx_track_EnergyClusters_TPCVmain[i],3);
	track-> SetdEdx(1,TPC_main_vtx_track_EnergyClusters_TPCVgap[i],4);
	
	track-> SetCharge(1,TPC_main_vtx_track_charge[i]);
	track-> SetChiSq(1,TPC_main_vtx_track_chi2[i]);
	track-> SetNDF(1,TPC_main_vtx_track_ndf[i]);
    }
    for (int i=0;i<nTPC_Tracks;i++)
    {
	if (TPC_primary_vtx_pT[i] < 0) break;
	DTEvent -> AddTrack();
	DataTreeTrack* track = DTEvent -> GetLastTrack();
	
	track -> SetType(1);//Primary vertex track
	
	track -> SetPt(0,TPC_primary_vtx_pT[i]);
	track -> SetEta(0,TPC_primary_vtx_eta[i]);
	track -> SetPhi(0,TPC_primary_vtx_phi[i]);
	track -> SetDCA(0,TPC_primary_vtx_impact_point_X[i]-Primary_Vertex_X,TPC_primary_vtx_impact_point_Y[i]-Primary_Vertex_Y,TPC_primary_vtx_impact_point_Z[i]-Primary_Vertex_Z);
	
	track-> SetNofHits(0,TPC_primary_vtx_track_nClusters_Total[i],0);
	track-> SetNofHits(0,TPC_primary_vtx_track_nClusters_TPCV1[i],1);
	track-> SetNofHits(0,TPC_primary_vtx_track_nClusters_TPCV2[i],2);
	track-> SetNofHits(0,TPC_primary_vtx_track_nClusters_TPCVmain[i],3);
	track-> SetNofHits(0,TPC_primary_vtx_track_nClusters_TPCVgap[i],4);
	
	track-> SetNofHitsPotential(0,TPC_primary_vtx_track_nClustersPotential_Total[i],0);
	track-> SetNofHitsPotential(0,TPC_primary_vtx_track_nClustersPotential_TPCV1[i],1);
	track-> SetNofHitsPotential(0,TPC_primary_vtx_track_nClustersPotential_TPCV2[i],2);
	track-> SetNofHitsPotential(0,TPC_primary_vtx_track_nClustersPotential_TPCVmain[i],3);
	track-> SetNofHitsPotential(0,TPC_primary_vtx_track_nClustersPotential_TPCVgap[i],4);
	
	track-> SetdEdx(0,TPC_primary_vtx_track_EnergyClusters_Total[i],0);
	track-> SetdEdx(0,TPC_primary_vtx_track_EnergyClusters_TPCV1[i],1);
	track-> SetdEdx(0,TPC_primary_vtx_track_EnergyClusters_TPCV2[i],2);
	track-> SetdEdx(0,TPC_primary_vtx_track_EnergyClusters_TPCVmain[i],3);
	track-> SetdEdx(0,TPC_primary_vtx_track_EnergyClusters_TPCVgap[i],4);
	
	track-> SetCharge(0,TPC_primary_vtx_track_charge[i]);
	track-> SetChiSq(0,TPC_primary_vtx_track_ChiSq[i]);
	track-> SetNDF(0,TPC_primary_vtx_track_NDF[i]);
	
	track -> SetPt(1,TPC_primary_vtx_track_pT[i]);
	track -> SetEta(1,TPC_primary_vtx_track_eta[i]);
	track -> SetPhi(1,TPC_primary_vtx_track_phi[i]);
	track -> SetDCA(1,TPC_primary_vtx_impact_point_X[i]-Main_Vertex_X,TPC_primary_vtx_impact_point_Y[i]-Main_Vertex_Y,TPC_primary_vtx_impact_point_Z[i]-Main_Vertex_Z);
	
	track-> SetNofHits(1,TPC_primary_vtx_track_nClusters_Total[i],0);
	track-> SetNofHits(1,TPC_primary_vtx_track_nClusters_TPCV1[i],1);
	track-> SetNofHits(1,TPC_primary_vtx_track_nClusters_TPCV2[i],2);
	track-> SetNofHits(1,TPC_primary_vtx_track_nClusters_TPCVmain[i],3);
	track-> SetNofHits(1,TPC_primary_vtx_track_nClusters_TPCVgap[i],4);
	
	track-> SetNofHitsPotential(1,TPC_primary_vtx_track_nClustersPotential_Total[i],0);
	track-> SetNofHitsPotential(1,TPC_primary_vtx_track_nClustersPotential_TPCV1[i],1);
	track-> SetNofHitsPotential(1,TPC_primary_vtx_track_nClustersPotential_TPCV2[i],2);
	track-> SetNofHitsPotential(1,TPC_primary_vtx_track_nClustersPotential_TPCVmain[i],3);
	track-> SetNofHitsPotential(1,TPC_primary_vtx_track_nClustersPotential_TPCVgap[i],4);
	
	track-> SetdEdx(1,TPC_primary_vtx_track_EnergyClusters_Total[i],0);
	track-> SetdEdx(1,TPC_primary_vtx_track_EnergyClusters_TPCV1[i],1);
	track-> SetdEdx(1,TPC_primary_vtx_track_EnergyClusters_TPCV2[i],2);
	track-> SetdEdx(1,TPC_primary_vtx_track_EnergyClusters_TPCVmain[i],3);
	track-> SetdEdx(1,TPC_primary_vtx_track_EnergyClusters_TPCVgap[i],4);
	
	track-> SetCharge(1,TPC_primary_vtx_track_charge[i]);
	track-> SetChiSq(1,TPC_primary_vtx_track_ChiSq[i]);
	track-> SetNDF(1,TPC_primary_vtx_track_NDF[i]);
    }
//     for (int i=0;i<nTPC_Tracks;i++)
//     {
// 	if (TPC_track_pT[i] <= nUndefinedValue) break;
// 	DTEvent -> AddTrack();
// 	DataTreeTrack* track = DTEvent -> GetTrack(i);
// 	
// 	track -> SetPt(0,TPC_track_pT[i]);
// 	track -> SetEta(0,1/2.*TPC_track_eta[i]);//NA61 code specifics (bug fixing)
// 	track -> SetPhi(0,TPC_track_phi[i]);
// 
// 	track -> SetDCA(0,TPC_track_DCAtoVertex_X[i],TPC_track_DCAtoVertex_Y[i],nUndefinedValue);
// 	
// 	track-> SetNofHits(0,TPC_track_nClusters_Total[i],0);
// 	track-> SetNofHits(0,TPC_track_nClusters_TPCV1[i],1);
// 	track-> SetNofHits(0,TPC_track_nClusters_TPCV2[i],2);
// 	track-> SetNofHits(0,TPC_track_nClusters_TPCVmain[i],3);
// 	track-> SetNofHits(0,TPC_track_nClusters_TPCVgap[i],4);
// 	
// 	track-> SetNofHitsPotential(0,TPC_track_nClustersPotential_Total[i],0);
// 	track-> SetNofHitsPotential(0,TPC_track_nClustersPotential_TPCV1[i],1);
// 	track-> SetNofHitsPotential(0,TPC_track_nClustersPotential_TPCV2[i],2);
// 	track-> SetNofHitsPotential(0,TPC_track_nClustersPotential_TPCVmain[i],3);
// 	track-> SetNofHitsPotential(0,TPC_track_nClustersPotential_TPCVgap[i],4);
//     }
}
//================================================================> FINISH <==============================================================
void DataTreeConverterNA61::Finish()
{
    fTreeQA -> Write();
    fTreeFile -> Write();
    fTreeFile -> Close();
    std::cout << "Output in file " << sOutputFileName << std::endl;
}

ClassImp(DataTreeConverterNA61)
