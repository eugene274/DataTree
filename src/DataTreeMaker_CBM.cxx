//TODO runid, eventid, vertex, fitter!, match in STS, constants

#include "DataTreeMaker_CBM.h"
#include <iostream>
#include <fstream>
#include <vector>
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

#include "TDirectory.h"
#include "CbmPsdHit.h"
#include "CbmPsdDigi.h"
#include "CbmMCTrack.h"
#include "CbmStsTrack.h"
#include "CbmGlobalTrack.h"
#include "CbmTrackMatchNew.h"
#include "CbmTofHit.h"
// #include "CbmTrackMatch.h"

#include "FairMCEventHeader.h"
#include "UEvent.h"
#include "UParticle.h"

#include "CbmL1PFFitter.h"
#include "CbmKFVertex.h"
#include "L1Field.h"
#include "KFParticleTopoReconstructor.h"

#include "DataTreeEvent.h"
#include "DataTreeTrack.h"

//=================================================================> MAIN <===============================================================
DataTreeMaker_CBM::DataTreeMaker_CBM()
  : FairTask("DataTreeMaker_CBM",1),
  uEvent (new UEvent()),
  fTreeEvents (new TChain("events","events"))/*,
  DTEvent (new DataTreeEvent())*/
  
//     l1 (new CbmL1())
{
      DTEvent = new DataTreeEvent();
}
DataTreeMaker_CBM::~DataTreeMaker_CBM()
{
    
}
//=================================================================> INIT <===============================================================
//----------------------------------
InitStatus DataTreeMaker_CBM::Init()
{
    Init_Input();
    Init_Output();
}
//----------------------------------
void DataTreeMaker_CBM::Init_Input()
{
    FairRootManager* ioman = FairRootManager::Instance();
    fPrimVtx = (CbmVertex*) ioman->GetObject("PrimaryVertex");
    fHeader = (FairMCEventHeader*) ioman->GetObject("MCEventHeader.");
    //BEGIN temp for bad phi data
    fTreeEvents -> Add(sInputFileName);
    fTreeEvents -> SetBranchAddress("event",&uEvent);
    //END temp for bad phi data
//     fHeader = (CbmMCEventHeader*) ioman->GetObject("MCEventHeader.");
    flistPSDhit = (TClonesArray*) ioman->GetObject("PsdHit");
    flistPSDdigit = (TClonesArray*) ioman->GetObject("PsdDigi");
    flistMCtrack = (TClonesArray*) ioman->GetObject("MCTrack");
    flistSTSRECOtrack = (TClonesArray*) ioman->GetObject("StsTrack");
    flistSTStrackMATCH = (TClonesArray*) ioman->GetObject("StsTrackMatch");
    fGlobalTrackArray = (TClonesArray*) ioman->GetObject("GlobalTrack");
    fTofHitArray = (TClonesArray*) ioman->GetObject("TofHit");
    fTofHitMatchArray = (TClonesArray*) ioman->GetObject("TofHitMatch");
    
    Init_PSD();
}

//--------------------------------------------------------------------------------------------------
void DataTreeMaker_CBM::Init_PSD()
{
    int idx = 0;
    fPSD_X[idx++] = -39;
    fPSD_X[idx++] = -19;
    fPSD_X[idx++] =   1;
    fPSD_X[idx++] =  21;
    fPSD_X[idx++] =  41;
    fPSD_X[idx++] =  61;
    fPSD_X[idx++] = -59;
    fPSD_X[idx++] = -39;
    fPSD_X[idx++] = -19;
    fPSD_X[idx++] =   1;
    fPSD_X[idx++] =  21;
    fPSD_X[idx++] =  41;
    fPSD_X[idx++] =  61;
    fPSD_X[idx++] =  81;
    fPSD_X[idx++] = -59;
    fPSD_X[idx++] = -39;
    fPSD_X[idx++] = -19;
    fPSD_X[idx++] =   1;
    fPSD_X[idx++] =  21;
    fPSD_X[idx++] =  41;
    fPSD_X[idx++] =  61;
    fPSD_X[idx++] =  81;
    fPSD_X[idx++] = -59;
    fPSD_X[idx++] = -39;
    fPSD_X[idx++] = -19;
    fPSD_X[idx++] =   1;
    fPSD_X[idx++] =  21;
    fPSD_X[idx++] =  41;
    fPSD_X[idx++] =  61;
    fPSD_X[idx++] =  81;
    fPSD_X[idx++] = -59;
    fPSD_X[idx++] = -39;
    fPSD_X[idx++] = -19;
    fPSD_X[idx++] =   1;
    fPSD_X[idx++] =  21;
    fPSD_X[idx++] =  41;
    fPSD_X[idx++] =  61;
    fPSD_X[idx++] =  81;
    fPSD_X[idx++] = -39;
    fPSD_X[idx++] = -19;
    fPSD_X[idx++] =   1;
    fPSD_X[idx++] =  21;
    fPSD_X[idx++] =  41;
    fPSD_X[idx++] =  61;

    idx = 0;
    fPSD_Y[idx++] = -50;
    fPSD_Y[idx++] = -50;
    fPSD_Y[idx++] = -50;
    fPSD_Y[idx++] = -50;
    fPSD_Y[idx++] = -50;
    fPSD_Y[idx++] = -50;
    fPSD_Y[idx++] = -30;
    fPSD_Y[idx++] = -30;
    fPSD_Y[idx++] = -30;
    fPSD_Y[idx++] = -30;
    fPSD_Y[idx++] = -30;
    fPSD_Y[idx++] = -30;
    fPSD_Y[idx++] = -30;
    fPSD_Y[idx++] = -30;
    fPSD_Y[idx++] = -10;
    fPSD_Y[idx++] = -10;
    fPSD_Y[idx++] = -10;
    fPSD_Y[idx++] = -10;
    fPSD_Y[idx++] = -10;
    fPSD_Y[idx++] = -10;
    fPSD_Y[idx++] = -10;
    fPSD_Y[idx++] = -10;
    fPSD_Y[idx++] =  10;
    fPSD_Y[idx++] =  10;
    fPSD_Y[idx++] =  10;
    fPSD_Y[idx++] =  10;
    fPSD_Y[idx++] =  10;
    fPSD_Y[idx++] =  10;
    fPSD_Y[idx++] =  10;
    fPSD_Y[idx++] =  10;
    fPSD_Y[idx++] =  30;
    fPSD_Y[idx++] =  30;
    fPSD_Y[idx++] =  30;
    fPSD_Y[idx++] =  30;
    fPSD_Y[idx++] =  30;
    fPSD_Y[idx++] =  30;
    fPSD_Y[idx++] =  30;
    fPSD_Y[idx++] =  30;
    fPSD_Y[idx++] =  50;
    fPSD_Y[idx++] =  50;
    fPSD_Y[idx++] =  50;
    fPSD_Y[idx++] =  50;
    fPSD_Y[idx++] =  50;
    fPSD_Y[idx++] =  50;
}
  
//----------------------------------
void DataTreeMaker_CBM::Init_Output()
{
    Init_DataTreeEvent();
    Init_OutputTree();
}

//--------------------------------------------------------------------------------------------------
void DataTreeMaker_CBM::Init_DataTreeEvent()
{
    for (int i=0;i<nPSD_Modules;i++)
    {
	DTEvent -> AddPSDModule(10);
    }
}

//--------------------------------------------------------------------------------------------------
void DataTreeMaker_CBM::Init_OutputTree()
{
    fTreeFile = new TFile(sOutputFileName, "RECREATE");
    fTreeFile -> cd();
    fDataTree = new TTree("fDataTree","fDataTree");
//     fDataTree -> SetMaxTreeSize(9000000);
    
//     fDataTree -> Branch("DTEvent", &DTEvent, 256000, 3);    
    fDataTree -> Branch("DTEvent", &DTEvent);
}

//=================================================================> EXEC <===============================================================
void DataTreeMaker_CBM::Exec(Option_t* opt)
{
    Clear_Event();
    
    Read_Event();
    Read_PSD();
    Read_STS();
    Read_MC_Tracks();
    Link_STS();
    Read_V0_Candidate(0);
    Read_V0_Candidate(1);
    Read_TOF_Hits();
    
    DTEvent -> Process();
    fDataTree -> Fill();
}
//--------------------------------------------------------------------------------------------------
void DataTreeMaker_CBM::Clear_Event()
{
    MCTrackIDs.clear();
    nMCTrackIDs = 0;
    TrackIDs.clear();
    nTrackIDs = 0;
    DTEvent -> ClearEvent();
}

//--------------------------------------------------------------------------------------------------
void DataTreeMaker_CBM::Read_Event()
{
    if (!fHeader) cout << "No fHeader!" << endl;
    else
    {
	DTEvent -> SetRPAngle(fHeader -> GetRotZ());
	DTEvent -> SetImpactParameter(fHeader -> GetB());
	DTEvent -> SetRunId(fHeader -> GetRunID());
	DTEvent -> SetEventId(fHeader -> GetEventID());
	DTEvent -> SetMCVertexPosition(fHeader->GetX(),fHeader->GetY(),fHeader->GetZ());
    }
//     //BEGIN temp for bad phi data
//     fTreeEvents->GetEntry(fCurEvent-1);
//     if (Generator == 1){DTEvent -> SetRPAngle(uEvent->GetPhi());}//used only for DCM-QGSM
// //     cout << "RPAngle = " << DTEvent -> GetRPAngle() << endl;
//     //END temp for bad phi data
}

//--------------------------------------------------------------------------------------------------
void DataTreeMaker_CBM::Read_PSD()
{
    for (int i=0; i<nPSD_Modules; i++)
    {
	DTEvent -> GetPSDModule(i) -> SetPosition(fPSD_X[i], fPSD_Y[i], nUndefinedValue);
    }
    
    int nPSDdigits = flistPSDdigit->GetEntriesFast();
    CbmPsdDigi* digit = NULL;
    for (int i=0; i<nPSDdigits; i++)
    {
	digit  = (CbmPsdDigi*) flistPSDdigit -> At(i);
	if (!digit) continue;
	DTEvent -> GetPSDModule(digit->GetModuleID()-1) -> GetSection(digit->GetSectionID()-1) -> AddEnergy(digit->GetEdep());
    }
}

//--------------------------------------------------------------------------------------------------
int DataTreeMaker_CBM::GetMCTrackMatch(int idx)
{
    int result = nUndefinedValue;
    for (int i=0;i<DTEvent->GetNTracks();i++)
    {
	if (DTEvent->GetTrack(i)->GetMCTrackId() == i) result = i;
    }
    return result;
}

//--------------------------------------------------------------------------------------------------
void DataTreeMaker_CBM::Read_MC_Tracks()
{
    CbmMCTrack* mctrack;
    Double_t p, px, py, pz, energy, mass, charge, pT, phi, eta, MC_px, MC_py, MC_pz, MC_p, MC_pT, MC_phi, MC_eta, MC_E, MC_y, chi2;
    int nTracks = flistMCtrack->GetEntries();
    for (int i=0;i<nTracks;i++)
    {
	mctrack = (CbmMCTrack*) flistMCtrack->At(i);
	int motherid = mctrack->GetMotherId();
	if (motherid != -1) continue;
	int type = mctrack->GetPdgCode();
	if (type < 1000000000)
	{
	    charge = mctrack->GetCharge();
	    mass = mctrack->GetMass();
	}
	else
	{
	    //pdg = 1000000000 + 10*1000*z + 10*a + i;
	    charge = TMath::Floor( ( type - 1000000000 ) / 10000 );
	    mass = TMath::Floor( ( type - 1000000000 -  10000 * charge ) / 10 );
	}
	MC_px = mctrack->GetPx();
	MC_py = mctrack->GetPy();
	MC_pz = mctrack->GetPz();
	MC_p = TMath::Sqrt(MC_px*MC_px + MC_py*MC_py + MC_pz*MC_pz);
	MC_pT = TMath::Sqrt(MC_px*MC_px + MC_py*MC_py);
	MC_phi = TMath::ATan2(MC_py,MC_px);
	if (MC_phi<-TMath::Pi()) MC_phi+=2*TMath::Pi();
	if (MC_phi>TMath::Pi()) MC_phi-=2*TMath::Pi();
	MC_eta = TMath::Log((MC_p+MC_pz)/(MC_p-MC_pz))/2.;
	MC_E = TMath::Sqrt(mass*mass+MC_p*MC_p);
	MC_y = TMath::Log((MC_E+MC_pz)/(MC_E-MC_pz))/2.;
	
	MCTrackIDs.push_back(i);
	DTEvent -> AddMCTrack();
	DataTreeMCTrack* DTMCTrack = DTEvent -> GetLastMCTrack();
	DTMCTrack->SetPx(MC_px);
	DTMCTrack->SetPy(MC_py);
	DTMCTrack->SetPz(MC_pz);
	DTMCTrack->SetP(MC_p);
	DTMCTrack->SetPt(MC_pT);
	DTMCTrack->SetEta(MC_eta);
	DTMCTrack->SetPhi(MC_phi);
	DTMCTrack->SetMass(mass);
	DTMCTrack->SetCharge(charge);
	DTMCTrack->SetEnergy(MC_E);
	DTMCTrack->SetY(MC_y);
	DTMCTrack->SetPdgId(type);
	DTMCTrack->SetMotherId(motherid);
	DTMCTrack->SetTrackId(GetMCTrackMatch(i));
    }
}

//--------------------------------------------------------------------------------------------------
void DataTreeMaker_CBM::Read_STS()
{
    CbmStsTrack* track;
    CbmTrackMatchNew* match;
    CbmMCTrack* mctrack;
	
    Double_t p, px, py, pz, energy, mass, pT, phi, eta, MC_px, MC_py, MC_pz, MC_p, MC_pT, MC_phi, MC_eta, chi2, Energy, y;
    Int_t type, trackID, motherid;
    bool matched = false;
    
    const FairTrackParam *trackParam;
    TVector3 momRec;

    Int_t nSTStracks = flistSTSRECOtrack->GetEntries();

    // Extrapolation track parameters back to primary vertex
    vector<CbmStsTrack> vRTracks;
    vector<CbmStsTrack> vRTracks_old;
    vRTracks.resize(nSTStracks);   
    
    //BEGIN fitter
    CbmL1PFFitter fitter;
    vector<float> vChiToPrimVtx;
    CbmKFVertex kfVertex;
    if(fPrimVtx) kfVertex = CbmKFVertex(*fPrimVtx);
    vector<L1FieldRegion> vField;
    vector<int> Pdg;
    Pdg.resize(nSTStracks);
    
    for (Int_t i=0; i<nSTStracks; i++)
    {
	matched = false;
	vRTracks[i] = *( (CbmStsTrack*) flistSTSRECOtrack->At(i));
	
	if(!track)
	{
	    cout << "ERROR: empty track!";
	}
        else
	{
	    match = (CbmTrackMatchNew*) flistSTStrackMATCH->At(i);
	    if (match != NULL)
	    {
		trackID = match->GetMatchedLink().GetIndex();
		if (trackID >= 0)
		{
		    mctrack = (CbmMCTrack*) flistMCtrack->At(trackID);
		    if (mctrack)
		    {
			matched = true;
			mass = mctrack->GetMass();
		    }
		}
	    }
	    track = &vRTracks[i];
	    
	    TrackIDs.push_back(i);
	    DTEvent->AddTrack();
	    //BEGIN first point
	    trackParam = track->GetParamFirst();	
	    trackParam->Momentum(momRec);
		
	    px = momRec.X();
	    py = momRec.Y();
	    pz = momRec.Z();
	    pT = TMath::Sqrt(px*px+py*py);
	    p = TMath::Sqrt(px*px+py*py+pz*pz);
	    phi = TMath::ATan2(py,px);//+(TMath::Pi()/2);
	    if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
	    if (phi>TMath::Pi()) phi-=2*TMath::Pi();
	    eta = TMath::Log((p+pz)/(p-pz))/2.;
	    
	    DataTreeTrack* DTTrack = DTEvent->GetTrack(i);
	    
	    DTTrack->SetPt(DataTreeTrack::kFirstPointTrack,pT);
	    DTTrack->SetPhi(DataTreeTrack::kFirstPointTrack,phi);
	    DTTrack->SetEta(DataTreeTrack::kFirstPointTrack,eta);
	    DTTrack->SetPx(DataTreeTrack::kFirstPointTrack,px);
	    DTTrack->SetPy(DataTreeTrack::kFirstPointTrack,py);
	    DTTrack->SetPz(DataTreeTrack::kFirstPointTrack,pz);
	    DTTrack->SetP(DataTreeTrack::kFirstPointTrack,p);
	    if (matched)
	    {
		Energy = TMath::Sqrt(mass*mass+p*p);
		y = TMath::Log((Energy+pz)/(Energy-pz))/2.;
		DTTrack->SetEnergy(DataTreeTrack::kFirstPointTrack,Energy);
		DTTrack->SetY(DataTreeTrack::kFirstPointTrack,y);
	    }
	    
	    DTTrack->SetNofHits(DataTreeTrack::kFirstPointTrack,track->GetNofHits());
	    DTTrack->SetFlag(DataTreeTrack::kFirstPointTrack,track->GetFlag());
	    DTTrack->SetChiSq(DataTreeTrack::kFirstPointTrack,track->GetChiSq());
	    DTTrack->SetNDF(DataTreeTrack::kFirstPointTrack,track->GetNDF());
	    
	    //END first point
	
	    //BEGIN last point
	    trackParam = track->GetParamLast();	
	    trackParam->Momentum(momRec);
		
	    px = momRec.X();
	    py = momRec.Y();
	    pz = momRec.Z();
	    pT = TMath::Sqrt(px*px+py*py);
	    p = TMath::Sqrt(px*px+py*py+pz*pz);
	    phi = TMath::ATan2(py,px);//+(TMath::Pi()/2);
	    if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
	    if (phi>TMath::Pi()) phi-=2*TMath::Pi();
	    eta = TMath::Log((p+pz)/(p-pz))/2.;
		    
	    DTTrack->SetPt(DataTreeTrack::kLastPointTrack,pT);
	    DTTrack->SetPhi(DataTreeTrack::kLastPointTrack,phi);
	    DTTrack->SetEta(DataTreeTrack::kLastPointTrack,eta);
	    DTTrack->SetPx(DataTreeTrack::kLastPointTrack,px);
	    DTTrack->SetPy(DataTreeTrack::kLastPointTrack,py);
	    DTTrack->SetPz(DataTreeTrack::kLastPointTrack,pz);
	    DTTrack->SetP(DataTreeTrack::kLastPointTrack,p);
	    if (matched)
	    {
		Energy = TMath::Sqrt(mass*mass+p*p);
		y = TMath::Log((Energy+pz)/(Energy-pz))/2.;
		DTTrack->SetEnergy(DataTreeTrack::kLastPointTrack,Energy);
		DTTrack->SetY(DataTreeTrack::kLastPointTrack,y);
	    }
	    
	    DTTrack->SetNofHits(DataTreeTrack::kLastPointTrack,track->GetNofHits());
	    DTTrack->SetFlag(DataTreeTrack::kLastPointTrack,track->GetFlag());
	    DTTrack->SetChiSq(DataTreeTrack::kLastPointTrack,track->GetChiSq());
	    DTTrack->SetNDF(DataTreeTrack::kLastPointTrack,track->GetNDF());
	    
	    //END last point
	    
	    Pdg[i] = 211;
	}
    }
    
    fitter.Fit(vRTracks, Pdg);
    fitter.GetChiToVertex(vRTracks, vField, vChiToPrimVtx, kfVertex, 1.e9f);//tracks array, field, dca over error
    //END fitter
    for (Int_t i=0; i<nSTStracks; i++)
    {
	matched = false;
	trackID = nUndefinedValue;
	track = &vRTracks[i];
	if(!track)
	{
	    cout << "ERROR: empty track!";
	}
	else
	{
	    match = (CbmTrackMatchNew*) flistSTStrackMATCH->At(i);
	    if (match != NULL)
	    {
		trackID = match->GetMatchedLink().GetIndex();
		if (trackID >= 0)
		{
		    mctrack = (CbmMCTrack*) flistMCtrack->At(trackID);
		    if (mctrack)
		    {
			matched = true;
			mass = mctrack->GetMass();
		    }
		}
	    }
	    //BEGIN vertex point
	    
	    trackParam = track->GetParamFirst();	
	    trackParam->Momentum(momRec);
		
	    px = momRec.X();
	    py = momRec.Y();
	    pz = momRec.Z();
	    pT = TMath::Sqrt(px*px+py*py);
	    p = TMath::Sqrt(px*px+py*py+pz*pz);
	    phi = TMath::ATan2(py,px);//+(TMath::Pi()/2);
	    if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
	    if (phi>TMath::Pi()) phi-=2*TMath::Pi();
	    eta = TMath::Log((p+pz)/(p-pz))/2.;
	    
	    
	    DataTreeTrack* DTTrack = DTEvent -> GetTrack(i);
	    DTTrack->SetPt(DataTreeTrack::kVtxPointTrack,pT);
	    DTTrack->SetPhi(DataTreeTrack::kVtxPointTrack,phi);
	    DTTrack->SetEta(DataTreeTrack::kVtxPointTrack,eta);
	    DTTrack->SetPx(DataTreeTrack::kVtxPointTrack,px);
	    DTTrack->SetPy(DataTreeTrack::kVtxPointTrack,py);
	    DTTrack->SetPz(DataTreeTrack::kVtxPointTrack,pz);
	    DTTrack->SetP(DataTreeTrack::kVtxPointTrack,p);
	    if (matched)
	    {
		Energy = TMath::Sqrt(mass*mass+p*p);
		y = TMath::Log((Energy+pz)/(Energy-pz))/2.;
		DTTrack->SetEnergy(DataTreeTrack::kVtxPointTrack,Energy);
		DTTrack->SetY(DataTreeTrack::kVtxPointTrack,y);
	    }
	    
	    DTTrack->SetNofHits(DataTreeTrack::kVtxPointTrack,track->GetNofHits());
	    DTTrack->SetFlag(DataTreeTrack::kVtxPointTrack,track->GetFlag());
	    DTTrack->SetChiSq(DataTreeTrack::kVtxPointTrack,track->GetChiSq());
	    DTTrack->SetVtxChiSq(DataTreeTrack::kVtxPointTrack,vChiToPrimVtx[i]);
	    DTTrack->SetNDF(DataTreeTrack::kVtxPointTrack,track->GetNDF());
	    
	    //END vertex point
	    if (matched) DTTrack -> SetMCTrackId(trackID);
	}
    }
}

//--------------------------------------------------------------------------------------------------
void DataTreeMaker_CBM::Link_STS()
{
    bool found = false;
//     std::cout << "==================== MC" << std::endl;
    for (int i=0;i<DTEvent->GetNTracks();i++)
    {
	found = false;
	for (int j=0;j<DTEvent->GetNMCTracks();j++)
	{
// 	    std::cout<<j<<" "<<MCTrackIDs.at(j)<<" " <<DTEvent->GetTrack(i)->GetMCTrackId()<<" "<<DTEvent->GetMCTrack(j)->GetId() <<std::endl;
	    if (MCTrackIDs.at(j) == DTEvent->GetTrack(i)->GetMCTrackId())
	    {
// 		std::cout<<"track id: "<<i<<"; before: " << DTEvent->GetTrack(i)->GetMCTrackId()<<"; after: "<<j << std::endl;
		found = true;
		DTEvent->GetTrack(i)->SetMCTrackId(j);
	    }
	}
	if (!found)
	{
// 	    std::cout<<"track id: "<<i<<"; before: " << DTEvent->GetTrack(i)->GetMCTrackId()<<"; after: "<<nUndefinedValue << std::endl;
	    DTEvent->GetTrack(i)->SetMCTrackId(nUndefinedValue);
	}
    }
    for (int j=0;j<DTEvent->GetNMCTracks();j++)
    {
	found = false;
	for (int i=0;i<DTEvent->GetNTracks();i++)
	{
	    if (TrackIDs.at(i) == DTEvent->GetMCTrack(j)->GetTrackId())
	    {
		found = true;
		DTEvent->GetMCTrack(j)->SetTrackId(i);
	    }
	}
	if (!found) DTEvent->GetMCTrack(j)->SetTrackId(nUndefinedValue);
    }
}

//--------------------------------------------------------------------------------------------------
void DataTreeMaker_CBM::Read_TOF_Hits()
{
    for (Int_t igt = 0; igt < fGlobalTrackArray->GetEntries(); igt++)
    {  
	const CbmGlobalTrack* globalTrack = static_cast<const CbmGlobalTrack*>(fGlobalTrackArray->At(igt));

	Int_t stsTrackIndex = globalTrack->GetStsTrackIndex();
	if( stsTrackIndex<0 ) continue;

	CbmStsTrack* cbmStsTrack = (CbmStsTrack*) flistSTSRECOtrack->At(stsTrackIndex);
	
	const FairTrackParam *stsPar = cbmStsTrack->GetParamLast(); 
	TVector3 mom;
	stsPar->Momentum(mom);

	Double_t p = mom.Mag();
	Int_t q = stsPar->GetQp() > 0 ? 1 : -1;
	
	Double_t l = globalTrack->GetLength();// l is calculated by global tracking
	
	Double_t time;
	Int_t tofHitIndex = globalTrack->GetTofHitIndex();
	if (tofHitIndex >= 0) 
	{
	    const CbmTofHit* tofHit = static_cast<const CbmTofHit*>(fTofHitArray->At(tofHitIndex));
	    if(!tofHit) continue;
	    time = tofHit->GetTime();
	
	    Double_t m2 = p*p*(1./((l/time/SpeedOfLight)*(l/time/SpeedOfLight))-1.);
	    
	    DataTreeTOFHit* DTTOFHit = DTEvent->AddTOFHit();
	    DTTOFHit -> SetPosition(tofHit->GetX(),tofHit->GetY(),tofHit->GetZ());
	    DTTOFHit -> SetTime(time);
	    DTTOFHit -> SetMass2(m2);
	    DTTOFHit -> SetLength(l);
	    DTTOFHit -> SetBeta(l/time/SpeedOfLight);
	    DTTOFHit -> SetP(p);
	    DTTOFHit -> SetCharge(q);
	    
	    DataTreeTrack* track = DTEvent->GetTrack(stsTrackIndex);
	    track -> SetTOFHitId(DTEvent->GetNTOFHits()-1);
	}
	else
	    continue;
    }
}

//--------------------------------------------------------------------------------------------------
void DataTreeMaker_CBM::Read_V0_Candidate(int UseMCpid = 0)
{
    const KFParticleTopoReconstructor* tr;
    if (!UseMCpid){tr = fCbmKFParticleFinder_TOF -> GetTopoReconstructor();}
    if (UseMCpid){tr = fCbmKFParticleFinder_MC -> GetTopoReconstructor();}
    if (!tr) cout << "DataTreeMaker_CBM::Read_V0_Candidate_TOF: ERROR: no KFParticleTopoReconstructor!" << endl;
//     printf("Particles: %d\n",tr->GetParticles().size());
    const int ConstNV0Types = DTEvent -> GetNV0Types();
//     cout << "DataTreeMaker_CBM::Read_V0_Candidate: ConstNV0Types = " << ConstNV0Types << endl;
    int V0Mult[ConstNV0Types];
    for (int i=0;i<ConstNV0Types;i++){V0Mult[i]=0;}
	
    for(unsigned int iP=0; iP<tr->GetParticles().size(); iP++)
    {
    //  printf("PDG: %d\n",tr->GetParticles()[iP].GetPDG());
	bool accept_V0 = false;
	for (int i=0;i<ConstNV0Types;i++)
	{
	    if (tr->GetParticles()[iP].GetPDG() == DTEvent -> GetNV0Pdg(i)){accept_V0 = true; V0Mult[i]++;}
	}
	if (!accept_V0) continue;
	
	const KFParticle& V0 = tr->GetParticles()[iP];
	DataTreeV0Candidate* V0Candidate;
	if (!UseMCpid){V0Candidate = DTEvent -> AddV0CandidateTOFpid();}
	if (UseMCpid){V0Candidate = DTEvent -> AddV0CandidateMCpid();}
	if (!V0Candidate) cout << "DataTreeMaker_CBM::Read_V0_Candidate_TOF: ERROR: no V0Candidate!" << endl;
	
	V0Candidate -> SetPx(V0.GetPx());
	V0Candidate -> SetPy(V0.GetPy());
	V0Candidate -> SetPz(V0.GetPz());
	V0Candidate -> SetPt(V0.GetPt());
	V0Candidate -> SetEta(V0.GetEta());
	V0Candidate -> SetPhi(V0.GetPhi());
	V0Candidate -> SetP(V0.GetP());
	V0Candidate -> SetPdgId(V0.GetPDG());
	V0Candidate -> SetMass(V0.GetMass());
	V0Candidate -> SetChiSq(V0.GetChi2());
	V0Candidate -> SetCharge((int) V0.GetQ());
	

// 	printf("V0 mother p: %f, %f, %f\n",V0.GetPx(),V0.GetPy(),V0.GetPz());

// 	if(V0id>nV0){ printf("Number of V0 candidates exceed maximum (%d)",nV0); continue;}

	if(V0.NDaughters()!=nV0Daughters){ printf("Number of daughters not %d (%d)!\n",nV0Daughters, V0.NDaughters()); continue;}

	for(int iDaughter=0; iDaughter<V0.NDaughters(); iDaughter++)
	{
	    int daugherIndex = V0.DaughterIds()[iDaughter];
	    const KFParticle& daughter = tr->GetParticles()[daugherIndex];
	    
	    V0Candidate -> AddDaughter();
	    DataTreeV0Candidate* Daughter = V0Candidate -> GetDaughter(iDaughter);
	    Daughter -> SetPx(daughter.GetPx());
	    Daughter -> SetPy(daughter.GetPy());
	    Daughter -> SetPz(daughter.GetPz());
	    Daughter -> SetPt(daughter.GetPt());
	    Daughter -> SetEta(daughter.GetEta());
	    Daughter -> SetPhi(daughter.GetPhi());
	    Daughter -> SetP(daughter.GetP());
	    Daughter -> SetPdgId(daughter.GetPDG());
	    Daughter -> SetMass(daughter.GetMass());
	    Daughter -> SetChiSq(daughter.GetChi2());
	    if( daughter.NDaughters()==1 ){Daughter -> SetTrackId(daughter.DaughterIds()[0]);}
	}
    }

    int V0Mult_All = 0;
    for (int i=0;i<ConstNV0Types;i++)
    {
	if (!UseMCpid){DTEvent -> SetNV0SpecificCandidatesTOFpid(i,V0Mult[i]);}
	if (UseMCpid){DTEvent -> SetNV0SpecificCandidatesMCpid(i,V0Mult[i]);}
	V0Mult_All+=V0Mult[i];
    }
    if (!UseMCpid){DTEvent -> SetNV0CandidatesTOFpid(V0Mult_All);}
    if (UseMCpid){DTEvent -> SetNV0CandidatesMCpid(V0Mult_All);}
// printf("Multiplicities: %d, %d, %d\n",V0Mult[0],V0Mult[1],V0Mult[2]);
}

//================================================================> FINISH <==============================================================
void DataTreeMaker_CBM::Finish()
{
    cout << "DataTreeMaker_CBM::Finish" << endl;
    fDataTree -> Write();
    fTreeFile -> Write();
    fTreeFile -> Close();
}

ClassImp(DataTreeMaker_CBM)