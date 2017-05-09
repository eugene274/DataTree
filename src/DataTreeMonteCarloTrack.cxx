#include "DataTreeMonteCarloTrack.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeMonteCarloTrack::DataTreeMonteCarloTrack(int idx) : TObject()
{
    SetId(idx);
}
DataTreeMonteCarloTrack::~DataTreeMonteCarloTrack()
{
    
}

ClassImp(DataTreeMonteCarloTrack)