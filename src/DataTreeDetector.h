#ifndef DataTreeDetector_H
#define DataTreeDetector_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"
#include "DataTreeRecoTrack.h"
#include "DataTreeMonteCarloTrack.h"
#include "DataTreeModule.h"
#include "DataTreeTOFHit.h"

class DataTreeDetector : public TObject
{
    
public:
  
    DataTreeDetector();
    ~DataTreeDetector();
  
//Main
//--------------------------------------------------------  
void Init()
{
    for (int i=0;i<RecoTracks.size();i++)
    {
	RecoTracks.at(i) -> SetUndefinedValues();
    }
    for (int i=0;i<MCTracks.size();i++)
    {
	MCTracks.at(i) -> SetUndefinedValues();
    }
    for (int i=0;i<Modules.size();i++)
    {
	Modules.at(i) -> ClearPosition();
	Modules.at(i) -> ClearEnergy();
    }
    for (int i=0;i<Hits.size();i++)
    {
	Hits.at(i) -> ClearEvent();
    }
}
//--------------------------------------------------------   
void ClearEvent()
{
    for (int i=0;i<RecoTracks.size();i++)
    {
	RecoTracks.at(i) -> Init();
    }
    for (int i=0;i<MCTracks.size();i++)
    {
	MCTracks.at(i) -> Init();
    }
    for (int i=0;i<Modules.size();i++)
    {
	Modules.at(i) -> Init();
    }
}
//Setters
//--------------------------------------------------------
void AddRecoTrack(DataTreeRecoTrack* track){ RecoTracks.push_back(track); }
void AddMCTrack(DataTreeMonteCarloTrack* track){ MCTracks.push_back(track); }
void AddModule(DataTreeModule* module){ Modules.push_back(module); }
void AddTOFHit(DataTreeTOFHit* hit){ Hits.push_back(hit); }
//Getters
//--------------------------------------------------------
DataTreeRecoTrack* GetRecoTrack(int idx){ if (idx < RecoTracks.size()) return RecoTracks.at(idx); else { std::cout << "Wrong index! " << idx << "/" << RecoTracks.size() << std::endl; return 0x0;} }
DataTreeMonteCarloTrack* GetMonteCarloTrack(int idx){ if (idx < MCTracks.size()) return MCTracks.at(idx); else { std::cout << "Wrong index! " << idx << "/" << MCTracks.size() << std::endl; return 0x0;} }
DataTreeModule* GetModule(int idx){ if (idx < Modules.size()) return Modules.at(idx); else { std::cout << "Wrong index! " << idx << "/" << Modules.size() << std::endl; return 0x0;} }
DataTreeTOFHit* GetTOFHit(int idx){ if (idx < Hits.size()) return Hits.at(idx); else { std::cout << "Wrong index! " << idx << "/" << Hits.size() << std::endl; return 0x0;} }
    
private:
    
    std::vector <DataTreeRecoTrack*> RecoTracks;
    std::vector <DataTreeMonteCarloTrack*> MCTracks;
    std::vector <DataTreeModule*> Modules;
    std::vector <DataTreeTOFHit*> Hits;
    
    ClassDefNV(DataTreeDetector, 1)
};

#endif