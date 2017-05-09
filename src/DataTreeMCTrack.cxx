#include "DataTreeMCTrack.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeMCTrack::DataTreeMCTrack(int idx) : TObject()
{
    SetId(idx);
    SetUndefinedValues();
}

DataTreeMCTrack::~DataTreeMCTrack()
{
    
}

void DataTreeMCTrack::SetUndefinedValues()
{
    double fValue = -999;//TEST
    ifPrimary = false;
    px = fValue;
    py = fValue;
    pz = fValue;
    p = fValue;
    pT = fValue;
    eta = fValue;
    phi = fValue;
    Mass = fValue;
    Rapidity = fValue;
    Energy = fValue;
    Charge = fValue;
    MotherId = fValue;
    TOFPosition[0] = fValue;
    TOFPosition[1] = fValue;
    TOFPosition[2] = fValue;
    TOF_pT = fValue;
    TOF_eta = fValue;
    TOF_phi = fValue;
    PSDPosition[0] = fValue;
    PSDPosition[1] = fValue;
    PSDPosition[2] = fValue;
    PSD_pT = fValue;
    PSD_phi = fValue;
    PSD_eta = fValue;
    PSDEnergy = fValue;
    TrackId = fValue;
    PSDSectionId = fValue;
    TOFSegmentId = fValue;
}
    
    


ClassImp(DataTreeMCTrack)

