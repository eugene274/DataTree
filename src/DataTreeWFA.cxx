#include "DataTreeWFA.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeWFA::DataTreeWFA(int idx) : TObject()
{
    SetId(idx);
}
DataTreeWFA::~DataTreeWFA()
{
    
}

ClassImp(DataTreeWFA)
