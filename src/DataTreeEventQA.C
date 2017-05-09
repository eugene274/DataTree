#include <iostream>
#include <fstream>
//#include "/lustre/nyx/hades/user/parfenov/Soft/CBMRoot/src/analysis/flow/DataTree/DataTreeTrack.h"
{
bool SaveOutput = true;

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
//TH1D* hWFA;

TH1D* hNClusters[5];
TH1D* hpT;
TH1D* hphi;
TH1D* heta;
TH2D* hetaVSphi;
TH1D* hChi2NDF;
TH1D* hCharge;
TH1D* hp[4] //px,py,pz,p;
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
    hE = new TH1D("hE", "E", 3000, 0, 3000);
    hMEcorr = new TH2D("hMEcorr", "M : E", 500, 0, 500, 500, 0, 100);
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
    
    
    hTrackPtEta = new TH2D("hTrackPtEta", "pT : eta", nbins, pt_low, pt_high, nbins, eta_low, eta_high);
    hTrackPtPhi = new TH2D("hTrackPtPhi", "pT : phi", nbins, pt_low, pt_high, nbins, phi_low, phi_high);
    hTrackEtaPhi = new TH2D("hTrackEtaPhi", "eta : phi", nbins, eta_low, eta_high, nbins, phi_low, phi_high);
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
	hPSDvsTPC     = new TH1D("hPSDvsTPC","dEdx in PSD vs Multiplicity in TPC",nbins,0.,1000.,500,0.,500.);
	hTriggerT4    = new TH1D("hTriggerT4","Trigger T4",2,0.,2.);
	hVertexZcut   = new TH1D("hVertexZcut","Vertex Z with quality cut",2000,-1000.,1000.);
	hVertexXYcut  = new TH1D("hVertexXYcut","Vertex XY with quality cut",2000,-100.,100.,2000,-100.,100.);
	//hWFA;

	for (int i=0;i<5;i++){
		hNClusters[i] = new TH1D(Form("hNClusters%i",i),Form("N Clusters %i",i),nbins,0.,100.);
	}
	hpT       = new TH1D("hpT","p_{T}",nbins,pt_low,pt_high);
	hphi      = new TH1D("hphi","#phi",nbins,phi_low,phi_high);
	heta      = new TH1D("heta","#eta",nbins,eta_low,eta_high);
	hetaVSphi = new TH1D("hetaVSphi","#eta : #phi",nbins,eta_low,eta_high,nbins,phi_low,phi_high);
	hChi2NDF  = new TH1D("hChi2NDF","#chi^{2}/NDF",500,0.,1.);
	hCharge   = new TH1D("hCharge","charge",20,-10.,10.);
	
	for (int i=0;i<4;i++){
		hp[4] = new TH1D(Form("hp%i",i),Form("Momentum %i",i),nbins,p_low,p_high); //px,py,pz,p;
	}
	//hqpVSdEdx;


    
}

/*void Init_Chain(Int_t num = 0, Int_t step = 10)
{
    fChain = new TChain("fDataTree","fDataTree");
    TString dir;
    for (UInt_t i=num*step; i<(num+1)*step; i++)
    {
        dir = ""; dir += i;
        fChain->Add( dir + "/tree.root");
    }    
    fChain -> SetBranchAddress("DTEvent",&DTEvent);
}
*/
void Init_Chain(TChain* InChain)
{
	fChain = InChain;
	fChain -> SetBranchAddress("DTEvent",&DTEvent);

}

/*void Init_Chain(TString InputFileName)
{
	fChain = new TChain("fDataTree");
	fChain ->Add(InputFileName.Data());
	fChain -> SetBranchAddress("DTEvent",&DTEvent);

}*/

void Fill_Histograms(int type=0)
{
    double Mreco = DTEvent -> GetNTracks();
    double Msim = DTEvent -> GetNMCTracks();
    double E = DTEvent -> GetPSDEnergy();
    double B = DTEvent -> GetImpactParameter();
	double Ntrig = DTEvent -> GetNTriggers();
	double NBPD  = DTEvent -> GetNBPDs();

//     if (B<14) return;
    
    hMreco -> Fill(Mreco);
    hMsim -> Fill(Msim);
    hE -> Fill(E);
    hMEcorr -> Fill(Mreco, E);
    hMMcorr -> Fill(Mreco, Msim);
    hMBcorr -> Fill(Mreco, B);
    hEBcorr -> Fill(E, B);

	hNTriggers->Fill(Ntrig);
	hNBPDs    ->Fill(NBPD);
    
	hVertexQ->Fill(DTEvent->GetVertexQuality(1));
	hVertexZQ->Fill(DTEvent->GetVertexQuality(1),DTEvent->GetVertexPositionComponent(0,2));

	if (DTEvent->GetVertexQuality(1)==0) hVertexZ[0]->Fill(DTEvent->GetVertexPositionComponent(0,2));
	if (DTEvent->GetVertexQuality(1)==1) hVertexZ[1]->Fill(DTEvent->GetVertexPositionComponent(0,2));
	if (DTEvent->GetVertexQuality(1)==2) hVertexZ[2]->Fill(DTEvent->GetVertexPositionComponent(0,2));
	hVtxAllZ->Fill(DTEvent->GetVertexPositionComponent(0,2));

    double pT;
    double eta;
    double phi;
    double mc_pT;
    double mc_eta;
    double mc_phi;
    double mc_y;
	int    flag;
	int    Nstations;
    
    for (int i=0;i<DTEvent->GetNTracks();i++)
    {
	track = DTEvent -> GetTrack(i);
	if (track->GetType()!=type) continue;

        hNdf      -> Fill (track->GetNDF(0));  
        if (track->GetNDF(0)!=0) hChi2     -> Fill (track->GetChiSq(0)/track->GetNDF(0));   
        hNhits    -> Fill (track->GetNofHits(0,0));  
        hNhitPot  -> Fill (track->GetNofHitsPotential(0,0));

        //if (track->GetMCTrackId() < 0) continue;
	//mctrack = DTEvent -> GetMCTrack(track->GetMCTrackId());

        //Bool_t cut = ( track->GetNofHits(0,0) > 3 ) && ( track->GetChiSq(0)/track->GetNDF(0) < 10 ) && ( TMath::Abs(track->GetVtxChiSq(0)) < 3)  ;
        
        //if (!cut) continue;
        
        pT = track->GetPt(0);
	eta = track->GetEta(0);
	phi = track->GetPhi(0);
	flag = track->GetFlag(0);
	Nstations = track->GetNStations(0);
	//mc_pT = mctrack->GetPt();
	//mc_eta = mctrack->GetEta();
	//mc_phi = mctrack->GetPhi();
        //mc_p = mctrack->GetP();
	hTrackPtEta -> Fill(pT,eta);
	hTrackPtPhi -> Fill(pT,phi);
	hTrackEtaPhi -> Fill(eta,phi);
// 	hMCTrackPtEta -> Fill(mc_pT,mc_eta);
// 	hMCTrackPtPhi -> Fill(mc_pT,mc_phi);
// 	hMCTrackEtaPhi -> Fill(mc_eta,mc_phi);
	//hMatchTrackPt -> Fill(pT,mc_pT);
	//hMatchTrackEta -> Fill(eta,mc_eta);
	//hMatchTrackPhi -> Fill(phi,mc_phi);

	for (Int_t j=0;j<3;++j){
		hTrackDCA[j]->Fill(track->GetDCAComponent(0,j));
		if (DTEvent->GetVertexQuality(1)==0) hCutsTrackDCA[j][0]->Fill(track->GetDCAComponent(0,j));
		if (DTEvent->GetVertexQuality(1)==1) hCutsTrackDCA[j][1]->Fill(track->GetDCAComponent(0,j));
		if (DTEvent->GetVertexQuality(1)==2) hCutsTrackDCA[j][2]->Fill(track->GetDCAComponent(0,j));
	}
	hTrackdEdx  ->Fill(track->GetdEdx(0,0));
	hTrackChiSq ->Fill(track->GetChiSq(0));
	hTrackCharge->Fill(track->GetCharge(0));
        
//         if (track->GetTOFHitId() > 0 && mctrack->GetPdgId() == 2212 )
        
        //if (track->GetTOFHitId() < 0) continue;
        //if ( mctrack->GetMotherId() != -1) continue;
        
        //DataTreeTOFHit *tof = DTEvent->GetTOFHit(track->GetTOFHitId());
        
        //if ( tof->GetBeta () > 0.93 + 0.1*tof->GetP() ) continue;
        
        //if (mctrack->GetPdgId() == 2212)  hTofProtons->Fill(tof->GetP(), tof->GetMass2 () );
        //if (mctrack->GetPdgId() == 211)   hTofPions->Fill(tof->GetP(), tof->GetMass2 () );
        //if (mctrack->GetPdgId() == 321)  hTofKaons->Fill(tof->GetP(), tof->GetMass2 () );
        
//        if (tof->GetCharge() > 0)
//        {
//            hTofPos->Fill(tof->GetP(), tof->GetMass2 () );
//            hTofPosBeta->Fill(tof->GetP(), tof->GetBeta () );
            
//        }
//        if (tof->GetCharge() < 0)
//        {
//            hTofNeg->Fill(tof->GetP(), tof->GetMass2 () );
//            hTofNegBeta->Fill(tof->GetP(), tof->GetBeta () );
//        }
        
        
    }
}
/*
void Write_Histograms(Int_t num = 0)
{
    TCanvas* cEventCorrelations = new TCanvas("cEventCorrelations", "Event variables", 900, 900);
    cEventCorrelations -> Divide(3,2);
    cEventCorrelations -> cd(1);
    hMreco -> Draw("");
    cEventCorrelations -> cd(2);
    hMsim -> Draw("");
    cEventCorrelations -> cd(3);
    hE -> Draw("");
    cEventCorrelations -> cd(4);
    hMEcorr -> Draw("colz");
    cEventCorrelations -> cd(5);
    hMMcorr -> Draw("colz");
    if (SaveOutput)
    {
//	cEventCorrelations -> SaveAs(fOutputDir+"cEventCorrelations.pdf");
//	cEventCorrelations -> SaveAs(fOutputDir+"cEventCorrelations.C");
    }
    
    TCanvas* cTrackCorrelations = new TCanvas("cTrackCorrelations", "Track Correlations", 900, 900);
    cTrackCorrelations -> Divide(3,3);
    cTrackCorrelations -> cd(1);
    hTrackPtEta -> Draw("colz");
    cTrackCorrelations -> cd(2);
    hTrackPtPhi -> Draw("colz");
    cTrackCorrelations -> cd(3);
    hTrackEtaPhi -> Draw("colz");
    cTrackCorrelations -> cd(4);
    hMCTrackPtEta -> Draw("colz");
    cTrackCorrelations -> cd(5);
    hMCTrackPtPhi -> Draw("colz");
    cTrackCorrelations -> cd(6);
    hMCTrackEtaPhi -> Draw("colz");
    cTrackCorrelations -> cd(7);
    hMatchTrackPt -> Draw("colz");
    cTrackCorrelations -> cd(8);
    hMatchTrackEta -> Draw("colz");
    cTrackCorrelations -> cd(9);
    hMatchTrackPhi -> Draw("colz");
    
    if (SaveOutput)
    {
//	cTrackCorrelations -> SaveAs(fOutputDir+"cTrackCorrelations.pdf");
//	cTrackCorrelations -> SaveAs(fOutputDir+"cTrackCorrelations.C");
    }
    
    TCanvas* cPSDGeometry = new TCanvas("cPSDGeometry", "PSD geometry", 600, 600);
    cPSDGeometry -> cd();
    for (int i=0;i<DTEvent->GetNPSDModules();i++)
    {
	DataTreePSDModule* mod = DTEvent->GetPSDModule(i);
	hPSDGeometry -> Fill(mod->GetPositionComponent(0),mod->GetPositionComponent(1),i+1);
    }
    hPSDGeometry -> SetLineWidth(2);
    hPSDGeometry -> SetMarkerSize(2);
    hPSDGeometry -> Draw("TEXT");
    if (SaveOutput)
    {
//	cPSDGeometry -> SaveAs(fOutputDir+"cPSDGeometry.pdf");
//	cPSDGeometry -> SaveAs(fOutputDir+"cPSDGeometry.C");
    }
///////////////////////////////////////////////////////////////////////////
    TFile *MyFile = new TFile(Form("hist_%d.root", num),"RECREATE"); 
    if ( MyFile->IsOpen() ) cout << "File opened successfully" << endl; 

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
    
    MyFile->Close();     
    
    
}
*/
void Write_Histograms(TString outFileName)
{
	TFile *MyFile = new TFile(outFileName.Data(),"RECREATE"); 
    if ( MyFile->IsOpen() ) cout << "File opened successfully" << endl; 

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
    
    MyFile->Close();     
    
    
}


/*
void DataTreeEventQA(Int_t num=0,Int_t step=10, int type, TString inputDir = "", TString outputDir = "QA")
{
    fInputDir = inputDir;
    fOutputDir = outputDir;
    Init_Histograms();
    Init_Chain(num,step);
    int nevents = fChain->GetEntries();
    int outputstep = 100;
    std::cout << "Entries = " << nevents << std::endl;
    for (int i=0;i<nevents;i++)
    {
	if ( i+1 % outputstep == 0) std::cout << i << "/" << nevents << "\r" << std::flush;
	fChain->GetEntry(i);
	Fill_Histograms(type);
    }
    Write_Histograms(num);
}

void DataTreeEventQA(int type, TString inputFileName = "", TString outputFileName="./hist_out.root")
{
    Init_Histograms();
    Init_Chain(inputFileName.Data());
    int nevents = fChain->GetEntries();
    int outputstep = 1;
    std::cout << "Entries = " << nevents << std::endl;
    for (int i=0;i<nevents;i++)
    {
	if ( i+1 % outputstep == 0) std::cout << i << "/" << nevents << "\r" << std::flush;
	fChain->GetEntry(i);
	Fill_Histograms(type);
    }
    Write_Histograms(outputFileName.Data());
}
*/

void DataTreeEventQA(int type, TChain* inChain, TString outputFileName="./hist_out.root")
{
    Init_Histograms();
    Init_Chain(inChain);
    int nevents = fChain->GetEntries();
    int outputstep = 100;
    std::cout << "Entries = " << nevents << std::endl;
    for (int i=0;i<nevents;i++)
    {
	if ( i+1 % outputstep == 0) std::cout << i << "/" << nevents << "\r" << std::flush;
	fChain->GetEntry(i);
	Fill_Histograms(type);
    }
    Write_Histograms(outputFileName.Data());
}

}
