#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "DataTreeEvent.h"
#include "DataTreeConverterNA61.h"
#include "TProfile.h"



using namespace std;

bool SaveOutput = true;
const int nMaxWFAsignals = 2000;
const int n_trigg = 2;
TString fInputDir;
TString fOutputDir;
TChain* fChain;
DataTreeEvent* DTEvent;
DataTreeMCTrack* mctrack;
DataTreeTrack* track;
TH1D* hMreco;
TH1D* hMsim;
TH1D* hE;
TH2D* hMEcorr;
TH2D* hMMcorr;
TH2D* hMBcorr;
TH2D* hEBcorr;
TH1D* hNTriggers;
TH1D* hNBPDs;
TH1D* hNdf;
TH1D* hChi2;
TH1D* hNhits;
TH1D* hNhitPot;
TH2D* hMCTrackPtEta;
TH2D* hMCTrackPtPhi;
TH2D* hMCTrackEtaPhi;
TH2D* hMCTrackProtonPtY;
TH2D* hMCTrackPionPtY;
TH2D* hMCTrackKaonPtY;
TH2D* hTrackPtEta;
TH2D* hTrackPtPhi;
TH2D* hTrackEtaPhi;
TH2D* hMatchTrackPt;
TH2D* hMatchTrackEta;
TH2D* hMatchTrackPhi;
TH1D* hTrackDCA[3];
TH1D* hCutsTrackDCA[3][3];
TH1D* hVertexZ[3];
TH1D* hVtxAllZ;
TH1D* hTrackdEdx;
TH1D* hTrackChiSq;
TH1D* hTrackCharge;
TH2D* hPSDGeometry;
TH2D* hTofPos;
TH2D* hTofNeg;
TH2D* hTofPosBeta;
TH2D* hTofNegBeta;
TH2D* hTofProtons;
TH2D* hTofPions;
TH2D* hTofKaons;
TH1D* hVertexQ;
TH2D* hVertexZQ;
////////////////////////
TH1D* hEnergyPSD;
TH2D* hPSDvsTPC;
TH1D* hTriggerT4;
TH1D* hVertexZcut;
TH2D* hVertexXYcut;
TH1D* hWFA;
TH1D* hWFAT1;
TH1D* hWFAT4;
TH2D* hBPD1xy;
TH2D* hBPD2xy;
TH2D* hBPD3xy;
TProfile* hBPD1xvstime;
TProfile* hBPD1yvstime;
TH1D* hNClusters[5];
TH1D* hpT;
TH1D* hphi;
TH1D* heta;
TH2D* hetaVSphi;
TH1D* hChi2NDF;
TH1D* hCharge;
TH1D* hp[4]; //px,py,pz,p;
TH2D* hPdEdx;
bool IsGoodEvent();
bool IsGoodTriggerT4();
bool IsGoodVertex();
bool IsGoodBPD();
int  WFAcoincidence();
int EventId;
namespace EnumTrigger{
    enum eTriggerId{
        kS1 = 0,
        kS2 = 1,
        kS3 = 2,
        kV1 = 3,
        kPSD= 4,
        kT1 = 5,
        kT2 = 6,
        kT4 = 7    
    };
}

int TriggersIdConvert[8]={-1,-1,-1,-1,-1,0,-1,1};
//TH1D* hqpVSdEdx;
////////////////////////
void Init_Histograms()
{
    TString s_proj[3] = {TString("x"), TString("y"), TString("z")};
    Int_t nbins = 1000;
    
    Double_t pt_low = 0.;    
    Double_t pt_high = 3.5;
    Double_t eta_low = -2.;
    Double_t eta_high = 8.;
    Double_t phi_low = -3.2;
    Double_t phi_high = 3.2;
    Double_t y_low = -2.;
    Double_t y_high = 5.;
    Double_t p_low = 0.;
    Double_t p_high = 10.;
    
    
    hMreco = new TH1D("hMreco","Mreco", 500, 0, 500);
    hMsim = new TH1D("hMsim", "Msim", 1500, 0, 1500);
    hE = new TH1D("hE", "E", 1000, 0, 10000);
    hMEcorr = new TH2D("hMEcorr", "M : E", 500, 0, 500, 1000, 0, 10000);
    hMMcorr = new TH2D("hMMcorr", "Mreco : Msim", 500, 0, 500, 1500, 0, 1500);
    hMBcorr = new TH2D("hMBcorr", "Mreco : B", 500, 0, 500, 200, 0, 20);
    hEBcorr = new TH2D("hEBcorr", "Epsd : B", 500, 0, 100, 200, 0, 20);
    hNTriggers = new TH1D("hNTriggers", "NTriggers", 10, 0., 10.);
    hNBPDs =     new TH1D("hNBPDs"    , "NBPDs"    , 5 , 0.,  5.);
    
    hNdf =     new TH1D("hNdf", "hNdf", 20, 0, 20);
    hChi2 =    new TH1D("hChi2", "hChi2", 500, 0, 20);
    hNhits =   new TH1D("hNhits", "hNhits", 250, 0, 250);
    hNhitPot = new TH1D("hNhitPot", "hNhitPot", 12, 0, 12);
    hMCTrackPtEta = new TH2D("hMCTrackPtEta", "MC pT : eta", nbins, pt_low, pt_high, nbins, eta_low, eta_high);
    hMCTrackPtPhi = new TH2D("hMCTrackPtPhi", "MC pT : phi", nbins, pt_low, pt_high, nbins, phi_low, phi_high);
    hMCTrackEtaPhi = new TH2D("hMCTrackEtaPhi", "MC eta : phi", nbins, eta_low, eta_high, nbins, phi_low, phi_high);
    hMCTrackProtonPtY = new TH2D("hMCTrackProtonPtY", "MC proton Y : pT", nbins, y_low, y_high, nbins, pt_low, pt_high);
    hMCTrackPionPtY = new TH2D("hMCTrackPionPtY", "MC pion Y : pT", nbins, y_low, y_high, nbins, pt_low, pt_high);
    hMCTrackKaonPtY = new TH2D("hMCTrackKaonPtY", "MC kaon Y : pT", nbins, y_low, y_high, nbins, pt_low, pt_high);
    
    
    hTrackPtEta = new TH2D("hTrackPtEta", "pT : eta", nbins, eta_low, eta_high, nbins, pt_low, pt_high);
    hTrackPtPhi = new TH2D("hTrackPtPhi", "pT : phi", nbins, pt_low, pt_high, nbins, phi_low, phi_high);
    hTrackEtaPhi = new TH2D("hTrackEtaPhi", "eta : phi", nbins, phi_low, phi_high, nbins, eta_low, eta_high);
    hMatchTrackPt = new TH2D("hMatchTrackPt", "pT reco : sim", nbins, pt_low, pt_high, nbins, pt_low, pt_high);
    hMatchTrackEta = new TH2D("hMatchTrackEta", "eta reco : sim", nbins, eta_low, eta_high, nbins, eta_low, eta_high);
    hMatchTrackPhi = new TH2D("hMatchTrackPhi", "phi reco : sim", nbins, phi_low, phi_high, nbins, phi_low, phi_high);
    for (Int_t i=0;i<3;++i){
        hTrackDCA[i] = new TH1D(Form("hTrackDCA%i",i),Form("TrackDCA %s", s_proj[i].Data()),1000,-50.,50.);
        for (Int_t j=0;j<3;++j){
            hCutsTrackDCA[i][j] = new TH1D(Form("hCutsTrackDCA%i%i",i,j),Form("TrackDCA %s with Q-cut %i", s_proj[i].Data(),j),1000,-50.,50.);
        }
        hVertexZ[i] = new TH1D(Form("hVertexZ%i",i),Form("Vertex Z %i",i), 2000,-1000,1000);
    }
    hVtxAllZ = new TH1D("hVtxAllZ","Vertex without quality cut",2000,-1000,1000);
    hTrackdEdx =   new TH1D("hTrackdEdx"  ,"TrackdEdx"  , 100,0. ,10.);
    hTrackChiSq =  new TH1D("hTrackChiSq" ,"TrackChiSq" , 500,0. ,50.);
    hTrackCharge = new TH1D("hTrackCharge","TrackCharge", 100,-2.,2. );
    hPSDGeometry = new TH2D("hPSDGeometry", "PSD Geometry", 2000, -100, 100, 2000, -100, 100);
    hTofPos = new TH2D("hTofPos", "hTofPos", nbins, p_low, p_high, nbins, -0.5, 2);
    hTofNeg = new TH2D("hTofNeg", "hTofNeg", nbins, p_low, p_high, nbins, -0.5, 2);
    hTofPosBeta = new TH2D("hTofPosBeta", "hTofPosBeta", nbins, p_low, p_high, nbins, 0, 1.1);
    hTofNegBeta = new TH2D("hTofNegBeta", "hTofNegBeta", nbins, p_low, p_high, nbins, 0, 1.1);
    hTofProtons = new TH2D("hTofProtons", "hTofProtons", nbins, p_low, p_high, nbins, -0.5, 2);
    hTofPions = new TH2D("hTofPions", "hTofPions", nbins, p_low, p_high, nbins, -0.5, 2);
    hTofKaons = new TH2D("hTofKaons", "hTofKaons", nbins, p_low, p_high, nbins, -0.5, 2);
    hVertexQ  = new TH1D("hVertexQ","Vertex quality",3,0.,3.);
    hVertexZQ = new TH2D("hVertexZQ","Vertex Z vs quality",3,0.,3.,2000,-1000.,1000.);
    hEnergyPSD    = new TH1D("hEnergyPSD","dEdx in the PSD modules",nbins,0.,1000.);
    hPSDvsTPC     = new TH2D("hPSDvsTPC","dEdx in PSD vs Multiplicity in TPC",nbins,0.,1000.,500,0.,500.);
    hTriggerT4    = new TH1D("hTriggerT4","Trigger T4",2,0.,2.);
    hVertexZcut   = new TH1D("hVertexZcut","Vertex Z with quality cut",2000,-1000.,1000.);
    hVertexXYcut  = new TH2D("hVertexXYcut","Vertex XY with quality cut",2000,-100.,100.,2000,-100.,100.);
    
    hWFA = new TH1D("hWFA","WFA beam (combined)",200,-20000,20000);
    hWFAT1 = new TH1D("hWFAT1","WFA beam (T1)",200,-20000,20000);
    hWFAT4 = new TH1D("hWFAT4","WFA beam (T4)",200,-20000,20000);
    hBPD1xy = new TH2D("hBPD1xy","BPD1:xy plane",100,-1.2,1.2,100,-1.2,1.2);
    hBPD2xy = new TH2D("hBPD2xy","BPD2:xy plane",100,-1.,0.2,100,-1.,1.);
    hBPD3xy = new TH2D("hBPD3xy","BPD3:xy plane",100,-1.2,0.,100,-1.,0.6);
    hBPD1xvstime = new TProfile("hBPD1xvstime","BPD1x vs time (Profile)",100,-10000,10000,-1,1);
    hBPD1yvstime = new TProfile("hBPD1yvstime","BPD1y vs time (Profile)",100,-10000,10000,-1,1);
    
    for (int i=0;i<5;i++){
        hNClusters[i] = new TH1D(Form("hNClusters%i",i),Form("N Clusters %i",i),nbins,0.,100.);
    }
    hpT       = new TH1D("hpT","p_{T}",nbins,pt_low,pt_high);
    hphi      = new TH1D("hphi","#phi",nbins,phi_low,phi_high);
    heta      = new TH1D("heta","#eta",nbins,eta_low,eta_high);
    hetaVSphi = new TH2D("hetaVSphi","#eta : #phi",nbins,eta_low,eta_high,nbins,phi_low,phi_high);
    hChi2NDF  = new TH1D("hChi2NDF","#chi^{2}/NDF",500,0.,1.);
    hCharge   = new TH1D("hCharge","charge",20,-10.,10.);
    
    for (int i=0;i<4;i++){
        hp[i] = new TH1D(Form("hp%i",i),Form("Momentum %i",i),nbins,p_low,p_high); //px,py,pz,p;
    }
    //hqpVSdEdx;
    hPdEdx = new TH2D("hPdEdx", "dE/dx: qp", 500, -10, 10, 500, 0, 4);
    
}
void Init_Chain(TChain* InChain)
{
    fChain = InChain;
    fChain -> SetBranchAddress("DTEvent",&DTEvent);
}
void Fill_Histograms(int type=0)
{
    //rewrite with for loop
    if (DTEvent->GetVertexQuality(1)==0) hVertexZ[0]->Fill(DTEvent->GetVertexPositionComponent(2,1));
    if (DTEvent->GetVertexQuality(1)==1) hVertexZ[1]->Fill(DTEvent->GetVertexPositionComponent(2,1));
    if (DTEvent->GetVertexQuality(1)==2) hVertexZ[2]->Fill(DTEvent->GetVertexPositionComponent(2,1));
    hVertexQ->Fill(DTEvent->GetVertexQuality(1));
    hVertexZQ->Fill(DTEvent->GetVertexQuality(1),DTEvent->GetVertexPositionComponent(2,1));
    hVtxAllZ->Fill(DTEvent->GetVertexPositionComponent(2,1));  
    
    double Mreco = DTEvent -> GetNTracks();
    double Msim = DTEvent -> GetNMCTracks();
    double E = DTEvent -> GetPSDEnergy();
    double B = DTEvent -> GetImpactParameter();
    double Ntrig = DTEvent -> GetNTriggers();
    double NBPD  = DTEvent -> GetNBPDs();
    
    // square cuts and BPD fill
    
    //if((DTEvent->GetBPD(0)->GetPositionComponent(0)) > -0.5 && (DTEvent->GetBPD(0)->GetPositionComponent(0)) < 0.3 && (DTEvent->GetBPD(0)->GetPositionComponent(1))>-0.7 && (DTEvent->GetBPD(0)->GetPositionComponent(1))<1.1){
    hBPD1xy->Fill(DTEvent->GetBPD(0)->GetPositionComponent(0),DTEvent->GetBPD(0)->GetPositionComponent(1));
    //  }
    
    //if((DTEvent->GetBPD(1)->GetPositionComponent(0))>-0.7 && (DTEvent->GetBPD(1)->GetPositionComponent(0))<-0.1 && (DTEvent->GetBPD(1)->GetPositionComponent(1))>-0.4 && (DTEvent->GetBPD(1)->GetPositionComponent(1))<0.3){
    hBPD2xy->Fill(DTEvent->GetBPD(1)->GetPositionComponent(0),DTEvent->GetBPD(1)->GetPositionComponent(1));
    //  }
    
    //if((DTEvent->GetBPD(2)->GetPositionComponent(0))>-1.05 && (DTEvent->GetBPD(2)->GetPositionComponent(0))<-0.35 && (DTEvent->GetBPD(2)->GetPositionComponent(1))>-0.55 && (DTEvent->GetBPD(2)->GetPositionComponent(1))<0.05){
    hBPD3xy->Fill(DTEvent->GetBPD(2)->GetPositionComponent(0),DTEvent->GetBPD(2)->GetPositionComponent(1));
    //    }
    
    //BPD time dependence
    
    hBPD1xvstime->Fill((DTEvent -> GetEventTimestamp()),(DTEvent->GetBPD(0)->GetPositionComponent(0)));
    hBPD1yvstime->Fill((DTEvent -> GetEventTimestamp()),(DTEvent->GetBPD(0)->GetPositionComponent(1)));
    
    // WFA fill'n'cuts
    // WFA T1(S1_1) Trigger (i_trigg=0)
    // WFA T4 Trigger (i_trigg=1)
    
    for (Int_t j=0;j<nMaxWFAsignals;j++){
            hWFAT1->Fill(DTEvent->GetWFA(0)->GetTime(0,j));
            hWFAT4->Fill(DTEvent->GetWFA(1)->GetTime(1,j));
    }
    
    //hMreco -> Fill(Mreco);
    hMsim -> Fill(Msim);
    hE -> Fill(E);
    //hMEcorr -> Fill(Mreco, E);
    hMMcorr -> Fill(Mreco, Msim);
    hMBcorr -> Fill(Mreco, B);
    hEBcorr -> Fill(E, B);
    hNTriggers->Fill(Ntrig);
    hNBPDs    ->Fill(NBPD);
    
    double pT;
    double eta;
    double phi;
    double mc_pT;
    double mc_eta;
    double mc_phi;
    double mc_y;
    int    flag;
    int    Nstations;
    int    charge;
    double NhitsTotal;
    double NhitsVTPC1 ;
    double NhitsVTPC2 ;
    double NhitsMTPC ;
    double Chi2  ;
    double dEdx ;
    int    Mult = 0;
    
    for (int i=0;i<DTEvent->GetNTracks();i++)
    {
        track = DTEvent -> GetTrack(i);
        if (track->GetType()!=type) continue;
        NhitsTotal = track->GetNofHits(0,0);
        NhitsVTPC1 = track->GetNofHits(0,1);
        NhitsVTPC2 = track->GetNofHits(0,2);
        NhitsMTPC = track->GetNofHits(0,3);
        
        Chi2 = track->GetChiSq(0);
        charge = track->GetCharge(0);
        dEdx = track->GetdEdx(0);
        
        pT = track->GetPt(0);
        eta = track->GetEta(0);
        phi = track->GetPhi(0);
        Double_t dcax = track->GetDCAComponent(0,0);
        Double_t dcay = track->GetDCAComponent(0,1);
        Bool_t cut =  ( NhitsMTPC < 90 ) && ( NhitsMTPC == 0 || NhitsMTPC > 10  )  ;
        cut = cut &&  ( NhitsVTPC1 == 0 || NhitsVTPC1 > 10  ) &&  ( NhitsVTPC2 == 0 || NhitsVTPC2 > 10  ) ;
        cut = cut &&  ( NhitsVTPC2 < 72 )  ;
        cut = cut &&  ( Chi2/track->GetNDF(0) < 10 ) &&  ( NhitsVTPC2 == 0 || NhitsVTPC2 > 10  ) ;
        cut = cut &&  ( dcax*dcax + dcay*dcay < 9 ) ;
        if (!cut) continue;
        hNdf      -> Fill ((Double_t)track->GetNDF(0));  
        if (track->GetNDF(0)!=0) hChi2     -> Fill (track->GetChiSq(0)/track->GetNDF(0));   
        hNhits    -> Fill (track->GetNofHits(0,0));  
        hNhitPot  -> Fill (track->GetNofHitsPotential(0,0));
        flag = track->GetFlag(0);
        Nstations = track->GetNStations(0);
        hTrackPtEta -> Fill(eta,pT);
        hTrackPtPhi -> Fill(pT,phi);
        hTrackEtaPhi -> Fill(phi,eta);
        for (Int_t j=0;j<3;++j){
            hTrackDCA[j]->Fill(track->GetDCAComponent(0,j));
            if (DTEvent->GetVertexQuality(1)==0) hCutsTrackDCA[j][0]->Fill(track->GetDCAComponent(0,j));
            if (DTEvent->GetVertexQuality(1)==1) hCutsTrackDCA[j][1]->Fill(track->GetDCAComponent(0,j));
            if (DTEvent->GetVertexQuality(1)==2) hCutsTrackDCA[j][2]->Fill(track->GetDCAComponent(0,j));
        }
        hTrackdEdx  ->Fill(track->GetdEdx(0,0));
        hTrackChiSq ->Fill(track->GetChiSq(0));
        hTrackCharge->Fill(track->GetCharge(0));        
        
        hPdEdx -> Fill ( charge*TMath::Log10(pT*TMath::CosH(eta)), dEdx );
        Mult++;
        
    }
    hMreco -> Fill(Mult);
    hMEcorr -> Fill(Mult, E);
    
}


void Write_Histograms(TString outFileName)
{
    TFile *MyFile = new TFile(outFileName.Data(),"RECREATE"); 
    if ( MyFile->IsOpen() ) cout << endl << "File opened successfully" << endl;
    else cout << endl << "Could not open the file! Exit..." << endl; 
    hMBcorr -> Write();
    hEBcorr -> Write();
    
    hTrackPtEta -> Write();
    hTrackPtPhi -> Write();
    hTrackEtaPhi -> Write();
    hMCTrackPtEta -> Write();
    hMCTrackPtPhi -> Write();
    hMCTrackEtaPhi -> Write();
    hMatchTrackPt -> Write();
    hMatchTrackEta -> Write();
    hMatchTrackPhi -> Write(); 
    hMreco -> Write();
    hMsim -> Write();
    hE -> Write();
    hMEcorr -> Write();
    hMMcorr -> Write();
    hNTriggers->Write();
    hNBPDs    ->Write();
    hVertexQ ->Write();
    hVertexZQ->Write();
    
    for (Int_t i=0;i<3;++i){
        hTrackDCA[i]->Write();
        for (Int_t j=0;j<3;++j){
            hCutsTrackDCA[i][j]->Write();
        }
    }
    for (Int_t i=0;i<3;++i){
        hVertexZ[i]->Write();
    }
    hVtxAllZ    ->Write();
    hTrackdEdx  ->Write();
    hTrackChiSq ->Write();
    hTrackCharge->Write();
    hNdf     ->Write();
    hChi2    ->Write();
    hNhits   ->Write();
    hNhitPot ->Write();
    
    hTofPos  ->Write();
    hTofNeg  ->Write();  
    hTofPosBeta  ->Write();
    hTofNegBeta  ->Write();  
    
    hTofProtons ->Write();  
    hTofPions ->Write();  
    hTofKaons ->Write();  
    
    hMCTrackProtonPtY -> Write();
    hMCTrackPionPtY ->   Write();
    hMCTrackKaonPtY ->   Write();
    hPdEdx -> Write();
    
    hBPD1xy->Write();
    hBPD2xy->Write();
    hBPD3xy->Write();
    hWFA->Write();
    hWFAT1->Write();
    hWFAT4->Write();
    hBPD1xvstime->Write();
    hBPD1yvstime->Write();
    
    MyFile->Close();     
    
    
}

bool IsGoodBPD(){
    // BPD cuts
    //rewrite in a form independent if statement plus introduce variable which is BPDgoodevent
    if((DTEvent->GetBPD(0)->GetPositionComponent(0))<-0.5 || DTEvent->GetBPD(0)->GetPositionComponent(0) > 0.3 || (DTEvent->GetBPD(0)->GetPositionComponent(1))<-0.7 || (DTEvent->GetBPD(0)->GetPositionComponent(1))>1.1)
        return false;
    
    
    if((DTEvent->GetBPD(1)->GetPositionComponent(0))<-0.7 || (DTEvent->GetBPD(1)->GetPositionComponent(0))>-0.1 || (DTEvent->GetBPD(1)->GetPositionComponent(1))<-0.4 || (DTEvent->GetBPD(1)->GetPositionComponent(1))>0.3)
        return false;
    
    
    if((DTEvent->GetBPD(2)->GetPositionComponent(0))<-1.05 || (DTEvent->GetBPD(2)->GetPositionComponent(0))>-0.35 || (DTEvent->GetBPD(2)->GetPositionComponent(1))<-0.55 || (DTEvent->GetBPD(2)->GetPositionComponent(1))>0.05)
        return false;
    return true;
}  
//WFA cuts
int WFAcoincidence(int TriggerId){
    int WFAcounter = 0;
    bool WFAjunk=false;
    float WFAsignal; 
    
        for (Int_t j=0;j<nMaxWFAsignals;j++)
        {
            WFAsignal = DTEvent->GetWFA(TriggersIdConvert[TriggerId])->GetTime(TriggersIdConvert[TriggerId],j);
            if(WFAsignal==-999) WFAjunk=true;
             if(TMath::Abs(WFAsignal)<1e-20 && WFAsignal!=0)  WFAjunk=true;
            if(WFAjunk==false) {
                WFAcounter++;
//                 std::cout<<EventId<<" WFAcounter= "<<WFAsignal<<"WFAcounter="<<WFAcounter<<std::endl;
            }
        }
    if(WFAcounter==0 && DTEvent->GetTrigger(TriggerId)->GetBool()==true){
        std::cout<<"WFAcounter=0, check DataTreeConverterNA61"<<std::endl;
        exit(0);
    }
    return WFAcounter;
}

bool IsGoodEvent(){
    if(!IsGoodBPD()) return false;
   if(WFAcoincidence(EnumTrigger::kT1)!=1) return false;
   if(WFAcoincidence(EnumTrigger::kT4)!=1) return false;
        if(!IsGoodVertex()) return false;
       if(!IsGoodTriggerT4()) return false;
    return true;
    
}

bool IsGoodVertex(){
    int fittedVertexID=1;
    if ( DTEvent->GetVertexQuality(fittedVertexID) != 0 ) return false;    //0 means good quality
    if ( DTEvent->GetVertexPositionComponent(2, fittedVertexID) > -590 || DTEvent->GetVertexPositionComponent(2, fittedVertexID) < -594 ) return false;    
    if ( DTEvent->GetVertexPositionComponent(1, fittedVertexID) > 0.4  || DTEvent->GetVertexPositionComponent(1, fittedVertexID) < -0.8 ) return false;    
    if ( DTEvent->GetVertexPositionComponent(0, fittedVertexID) > 0.0  || DTEvent->GetVertexPositionComponent(0, fittedVertexID) < -1.1 ) return false;
    return true;
}

bool IsGoodTriggerT4(){
    if ( DTEvent->GetTrigger(7)->GetBool() == false) return false;
    return true;
}

void DataTreeEventQA(int nEvents=-1,int type=0, TChain* inChain=0, TString outputFileName="./hist_out.root")
{
    Init_Histograms();
    Init_Chain(inChain);
    if(nEvents==-1) nEvents = fChain->GetEntries();
    std::cout << "Entries = " << nEvents << std::endl;
    for (int i=0;i<nEvents;i++)
    {
        EventId=i;
        if ( i%1000==0) std::cout << '\r' << i << "/" << nEvents << std::flush;
        fChain->GetEntry(i);
        
        if(IsGoodEvent()){
            Fill_Histograms(type);
        }
        
    }
    Write_Histograms(outputFileName.Data());
}
