//TODO runid, eventid, vertex, fitter!, match in STS, constants

#include "DataTreeMakerDoctor.h"
#include <iostream>
#include <fstream>
#include <vector>
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

#include "TDirectory.h"
#include "CbmStsTrack.h"
#include "CbmTrackMatchNew.h"

#include "DataTreeEvent.h"
#include "DataTreeTrack.h"

//=================================================================> MAIN <===============================================================
DataTreeMakerDoctor::DataTreeMakerDoctor()
  : FairTask("DataTreeMakerDoctor",1)
{
      DTEvent = new DataTreeEvent();
      fCurEvent = 0;
}

DataTreeMakerDoctor::~DataTreeMakerDoctor()
{
    
}
//=================================================================> INIT <===============================================================
InitStatus DataTreeMakerDoctor::Init()
{
    fCurEvent = 0;
    FairRootManager* ioman = FairRootManager::Instance();
    flistSTSRECOtrack = (TClonesArray*) ioman->GetObject("StsTrack");
    flistSTStrackMATCH = (TClonesArray*) ioman->GetObject("StsTrackMatch");
    
    InputChain_Init();
    OutputTree_Init();
}

//--------------------------------------------------------------------------------------------------
void DataTreeMakerDoctor::InputChain_Init()
{
    std::cout << sInputFileName << std::endl;
    fInputChain = new TChain("fDataTree","fDataTree");
    fInputChain -> Add(sInputFileName);
    fInputChain -> SetBranchAddress("DTEvent",&DTEvent);
}

//--------------------------------------------------------------------------------------------------
void DataTreeMakerDoctor::OutputTree_Init()
{
    fTreeFile = new TFile(sOutputFileName, "RECREATE");
    fTreeFile -> cd();
    fDataTree = new TTree("fDataTree","fDataTree");
    fDataTree -> SetMaxTreeSize(9000000);
    
    fDataTree -> Branch("DTEvent", &DTEvent, 256000, 3);    
}

//=================================================================> EXEC <===============================================================
void DataTreeMakerDoctor::Exec(Option_t* opt)
{    
    fInputChain -> GetEntry(fCurEvent);
    Cure_Tracks();
    fDataTree -> Fill();
    fCurEvent++;
}

//--------------------------------------------------------------------------------------------------
void DataTreeMakerDoctor::Cure_Tracks()
{
    CbmStsTrack* track;
    CbmTrackMatchNew* match;
	
    int trackID;
    int nSTStracks = flistSTSRECOtrack->GetEntries();
    int tdx = 0;
    int mdx = 0;
    
    for (Int_t i=0; i<nSTStracks; i++)
    {
	track = (CbmStsTrack*) flistSTSRECOtrack->At(i);
	if(!track) continue;

	match = (CbmTrackMatchNew*) flistSTStrackMATCH->At(i);
	if (match != NULL)
	{
	    trackID = match->GetMatchedLink().GetIndex();
	    if (trackID < 0) std::cout << "ERROR: TrackID < 0!" << std::endl;
	    if (trackID >= 0)
	    {
		DTEvent->GetTrack(tdx)->SetMCTrackId(mdx);
		mdx++;
	    }
	}
	else
	{
	    DTEvent->GetTrack(tdx)->SetMCTrackId(nUndefinedValue);
	}
	tdx++;
    }
}

//================================================================> FINISH <==============================================================
void DataTreeMakerDoctor::Finish()
{
    cout << "DataTreeMakerDoctor::Finish" << endl;
    fDataTree -> Write();
    fTreeFile -> Write();
    fTreeFile -> Close();
}

ClassImp(DataTreeMakerDoctor)