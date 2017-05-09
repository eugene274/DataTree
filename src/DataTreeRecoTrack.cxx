#include "DataTreeRecoTrack.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeRecoTrack::DataTreeRecoTrack(int idx) : TObject()
{
    SetId(idx);
}
DataTreeRecoTrack::~DataTreeRecoTrack()
{
    
}

ClassImp(DataTreeRecoTrack)