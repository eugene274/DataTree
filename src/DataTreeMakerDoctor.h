#ifndef DataTreeMakerDoctor_H
#define DataTreeMakerDoctor_H 1

#include "CbmStsEventData.h"


#include "FairTask.h"
#include <vector>
#include "TLorentzVector.h"
#include <map>
#include <cstring>

#include <iostream>
#include "TClonesArray.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "CbmTrackMatch.h"

#include "DataTreeEvent.h"
#include "DataTreeTrack.h"

class TClonesArray;
class CbmVertex;
class TDirectory;
class TH1F;
class TProfile;
class TH2F;

class DataTreeMakerDoctor : public FairTask
{
    
const double SpeedOfLight = 29.9792458;

public:
  
    DataTreeMakerDoctor();
    ~DataTreeMakerDoctor();
    
    virtual InitStatus Init();
    virtual void Exec(Option_t* opt);
    virtual void Finish();
    
    void SetOutputFile(TString filename) { sOutputFileName = filename; }
    void SetInputFile(TString filename) { sInputFileName = filename; }
    void SetGenerator(int fValue){Generator = fValue;}
    
    
private:
    
    DataTreeEvent* DTEvent;
   
    int Generator = 0; //0 -- URQMD, 1 -- DCM_QGSM
    
    int fCurEvent;
    TClonesArray* flistSTSRECOtrack;
    TClonesArray* flistSTStrackMATCH;
    
    TChain* fInputChain;
    TString sInputFileName;
    
    void InputChain_Init();
    void OutputTree_Init();
    void Cure_Tracks();
    
    const int nInfinity = 10000000;
    const int nUndefinedValue = -999;
    
    static const int nTPC_Tracks = 500;
    
    TFile* fTreeFile;
    TTree* fDataTree;
    TString sOutputFileName;
    
    ClassDefNV(DataTreeMakerDoctor, 1)
//     ClassDef(DataTreeMakerDoctor,1);
};

#endif