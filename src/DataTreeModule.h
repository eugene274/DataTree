#ifndef DataTreeModule_H
#define DataTreeModule_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"
#include "DataTreePSDSection.h"
#include "TMath.h"

class DataTreeModule : public TObject
{
    
public:
  
    DataTreeModule(int idx = 0);
    DataTreeModule(int idx, int pNsections);
    ~DataTreeModule();
    
int GetId(){return id;}

//========================================================
//--------------------------------------------------------
void Init()
{
    for (int i=0;i<nSections;i++)
    {
	SectionEnergy.push_back(0);
    }
}
//--------------------------------------------------------
void ClearPosition()
{
    Position[0]=-999;
    Position[1]=-999;
    Position[2]=-999;
}
//--------------------------------------------------------
void ClearEnergy()
{
    for (int i=0;i<SectionEnergy.size();i++)
    {
	SectionEnergy.at(i) = 0;
    }
}

//========================================================
//Position
//--------------------------------------------------------
void SetPosition(double fX, double fY, double fZ)
{
    Position[0]=fX;
    Position[1]=fY;
    Position[2]=fZ;
}
//--------------------------------------------------------
void SetPositionComponent(int idx, double fValue)
{
    Position[idx] = fValue;
}
//--------------------------------------------------------
double GetPositionComponent(int idx)
{
    return Position[idx];
}
//--------------------------------------------------------
double GetPhi()
{
    return TMath::ATan2(Position[1],Position[0]);
}
//========================================================
//Number of sections
//--------------------------------------------------------
void SetNSections(int pNsections)
{
    nSections = pNsections;
}
//--------------------------------------------------------
int GetNSections()
{
    return nSections;
}
//========================================================
//Energy
//--------------------------------------------------------
void AddEnergy(int section_id, double energy)
{
    if (section_id < SectionEnergy.size()) SectionEnergy.at(section_id) += energy;
    else std::cout << "Wrong section id: " << section_id << " of max " << SectionEnergy.size() << std::endl;
}
//--------------------------------------------------------
void SetEnergy(int section_id, double energy)
{
    if (section_id < SectionEnergy.size()) SectionEnergy.at(section_id) = energy;
    else std::cout << "Wrong section id: " << section_id << " of max " << SectionEnergy.size() << std::endl;
}
//--------------------------------------------------------
double GetEnergy(int section_id)
{
    if (section_id < SectionEnergy.size()) return SectionEnergy.at(section_id);
    else
    {
	std::cout << "Wrong section id: " << section_id << " of max " << SectionEnergy.size() << std::endl;
	return -999;
    }
}
//--------------------------------------------------------
double GetEnergy()
{
    double energy = 0;
    for (int i=0;i<SectionEnergy.size();i++)
    {
	energy += SectionEnergy.at(i);
    }
    return energy;
}
    
private:    
    void SetId(int idx){id = idx;}
      
    int nSections;			//number of section in the module
    
    int id;
    double Position[3];			//Position of the module in lab frame
    double Energy;			//energy deposit in the module
    int nFiredSections;			//number of sections where the energy was deposited
    std::vector<double> SectionEnergy;		//sections in the module
    
    ClassDefNV(DataTreeModule, 1)
};

#endif