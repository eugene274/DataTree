#include "DataTreeExperiment.h"
#include <iostream>
#include <vector>
#include "TObject.h"
#include "DataTreeTrack.h"
#include "DataTreePSDModule.h"
#include "DataTreeTOFHit.h"
#include "DataTreeV0Candidate.h"
#include "DataTreeMCTrack.h"
#include "DataTreeTrigger.h"
#include "DataTreeBPD.h"

DataTreeExperiment::DataTreeExperiment() : TObject(),
    CBM_Reco_Vertex (new DataTreeDetector()),
    CBM_Reco_Tracks (new DataTreeDetector()),
    CBM_Reco_V0_Tracks_RecoPid (new DataTreeDetector()),
    CBM_Reco_V0_Tracks_MCPid (new DataTreeDetector()),
    CBM_Reco_PSD (new DataTreeDetector()),
    CBM_Reco_TOF (new DataTreeDetector()),
    CBM_MC_Vertex (new DataTreeDetector()),
    CBM_MC_Tracks (new DataTreeDetector()),
    
    NA61_Reco_Main_Vertex (new DataTreeDetector()),
    NA61_Reco_Main_Vertex_Tracks (new DataTreeDetector()),
    NA61_Reco_Primary_Vertex (new DataTreeDetector()),
    NA61_Reco_Primary_Vertex_Tracks (new DataTreeDetector()),
    NA61_Reco_PSD (new DataTreeDetector()),
    NA61_Reco_BPD (new DataTreeDetector()),
    NA61_Reco_Triggers (new DataTreeDetector()),
    NA61_Reco_WFA (new DataTreeDetector())
{
//     std::cout << "Constructor" << std::endl;
}
DataTreeExperiment::~DataTreeExperiment()
{

}


ClassImp(DataTreeExperiment)

