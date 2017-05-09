#include "DataTreeModule.h"
#include <iostream>
#include <vector>
#include "TObject.h"

DataTreeModule::DataTreeModule(int idx) : TObject(),
nSections(0)
{
    SetId(idx);
}

DataTreeModule::DataTreeModule(int idx, int pNsections) : TObject(),
nSections(0)
{
    SetId(idx);
    nSections = pNsections;
}

DataTreeModule::~DataTreeModule()
{
    
}

ClassImp(DataTreeModule)

