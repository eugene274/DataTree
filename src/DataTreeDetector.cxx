#include "DataTreeDetector.h"
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

DataTreeDetector::DataTreeDetector() : TObject()
{
//     std::cout << "Constructor" << std::endl;
}
DataTreeDetector::~DataTreeDetector()
{

}


ClassImp(DataTreeDetector)

