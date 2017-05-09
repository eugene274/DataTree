#include "DataTreeMontCarloTrack.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeMontCarloTrack::DataTreeMontCarloTrack(int idx) : TObject()
{
    SetId(idx);
}
DataTreeMontCarloTrack::~DataTreeMontCarloTrack()
{
    
}

ClassImp(DataTreeMontCarloTrack)