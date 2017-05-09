#ifndef DataTreeMaker_CBM_H
#define DataTreeMaker_CBM_H 1

#include "CbmMCEventData.h"
#include "CbmPsdEventData.h"
#include "CbmStsEventData.h"

#include "CbmKFPartEfficiencies.h"
#include "CbmKFParticleFinder.h"

#include "FairTask.h"
#include "CbmVertex.h"
#include <vector>
#include "TLorentzVector.h"
#include <map>
#include <cstring>

#include "UEvent.h"
#include "TGraphErrors.h"
#include <iostream>
#include "TClonesArray.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "FairMCEventHeader.h"
#include "UEvent.h"

#include "CbmKFVertex.h"
#include "CbmTrackMatch.h"

#include "DataTreeEvent.h"
#include "DataTreeTrack.h"

class TClonesArray;
class CbmVertex;
class TDirectory;
class TH1F;
class TProfile;
class TH2F;


class DataTreeMaker_CBM : public FairTask
{
    
const double SpeedOfLight = 29.9792458;

public:
  
    DataTreeMaker_CBM();
    ~DataTreeMaker_CBM();
    
    virtual InitStatus Init();
    virtual void Exec(Option_t* opt);
    virtual void Finish();
    
    void SetPSDCoordinatesFileName(TString fileName) { sPSDCoordinatesFileName = fileName; } // File containing PSD module (x, y) coordinates in LAB
    void SetOutputFile(TString filename) { sOutputFileName = filename; }
    void SetInputFile(TString filename) { sInputFileName = filename; }
    void SetCbmKFParticleFinder_MC(CbmKFParticleFinder* kf) {fCbmKFParticleFinder_MC=kf;}
    void SetCbmKFParticleFinder_TOF(CbmKFParticleFinder* kf) {fCbmKFParticleFinder_TOF=kf;}
    void SetGenerator(int fValue){Generator = fValue;}
    
    
private:
    
    DataTreeEvent* DTEvent;
  
    TString sPSDCoordinatesFileName;
    TString sOutputFileName;

    CbmKFParticleFinder* fCbmKFParticleFinder_TOF;
    CbmKFParticleFinder* fCbmKFParticleFinder_MC;

    int Generator = 0; //0 -- URQMD, 1 -- DCM_QGSM
    //===================================================
    
    int fCurEvent;
//     CbmMCEventHeader* fHeader;
    CbmVertex* fPrimVtx;
    FairMCEventHeader* fHeader;
    TClonesArray* flistPSDhit;
    TClonesArray* flistPSDdigit;
    TClonesArray* flistMCtrack;
    TClonesArray* flistSTSRECOtrack;
    TClonesArray* flistSTStrackMATCH;
    TClonesArray *fGlobalTrackArray;
    TClonesArray *fTofHitArray;
    TClonesArray *fTofHitMatchArray;
    
    void Init_Input();
    void Init_PSD();
    void Init_Output();
    void Init_OutputTree();
    void Init_DataTreeEvent();
    void Clear_Event();
    void Read_Event();
    void Read_PSD();
    void Read_STS();
    void Link_STS();
    void Read_TOF_Hits();
    void Read_MC_Tracks();
      int GetMCTrackMatch(int idx);
    void Read_V0_Candidate(int UseMCpid);
    
    const int nInfinity = 10000000;
    const int nUndefinedValue = -999;
    static const int nV0Daughters = 2;
    
    static const int nPSD_Modules = 44;
    static const int nTPC_Tracks = 500;
    
    std::vector<int> MCTrackIDs;
    int nMCTrackIDs;
    std::vector<int> TrackIDs;
    int nTrackIDs;
    
    double fPSD_X[nPSD_Modules];
    double fPSD_Y[nPSD_Modules];
    double fPSD_Z[nPSD_Modules];
    
    TFile* fTreeFile;
    TTree* fDataTree;
    TChain* fTreeEvents;//temp for bad phi data
    UEvent* uEvent;//temp for bad phi data
    TString sInputFileName;//temp for bad phi data
    
    ClassDefNV(DataTreeMaker_CBM, 1)
//     ClassDef(DataTreeMaker_CBM,1);
};

#endif