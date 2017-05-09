#ifndef DataTreeRecoTrack_H
#define DataTreeRecoTrack_H 1

#include <vector>
#include <iostream>
#include <fstream>
#include "TClonesArray.h"
#include "TObject.h"
#include "TMath.h"

// #include "DataTreeConstants.h"

class DataTreeRecoTrack : public TObject
{
    
public:
    static const int nSubDetectors = 5;
    static const int nIndices = 5;
    
    DataTreeRecoTrack(int idx = 0);
    ~DataTreeRecoTrack();
    
    int GetId(){return id;}

//========================================================
//--------------------------------------------------------
void Init()
{
    for (int i=0;i<nPoints;i++)
    {
	pT.push_back(-999);
	eta.push_back(-999);
	phi.push_back(-999);
	px.push_back(-999);
	py.push_back(-999);
	pz.push_back(-999);
	p.push_back(-999);
	PdgId.push_back(-999);
	Mass.push_back(-999);
	Rapidity.push_back(-999);
	Charge.push_back(-999);
	for (int j=0;j<nSubDetectors;j++)
	{
	    NofHits[j].push_back(-999);
	    NofHitsPotential[j].push_back(-999);
	    dEdx[j].push_back(-999);
	}
	Flag.push_back(-999);
	Chi2.push_back(-999);
	NDF.push_back(-999);
	for (int j=0;j<3;j++)
	{
	    Impact_Point[j].push_back(-999);
	}
	Length.push_back(-999);
	for (int j=0;j<nIndices;j++)
	{
	    Index[j].push_back(-999);
	}
    }
    Daughters.clear();
}
//--------------------------------------------------------
void SetNpoints(int idx)
{
    nPoints = idx;
}
//Setters
//--------------------------------------------------------
void SetUndefinedValues()
{
    for (int i=0;i<nPoints;i++)
    {
	pT.at(i) = -999;
	eta.at(i) = -999;
	phi.at(i) = -999;
	px.at(i) = -999;
	py.at(i) = -999;
	pz.at(i) = -999;
	p.at(i) = -999;
	Rapidity.at(i) = -999;
	PdgId.at(i) = -999;
	Mass.at(i) = -999;
	Charge.at(i) = -999;
	for (int j=0;j<nSubDetectors;j++)
	{
	    NofHits[j].at(i) = -999;
	    NofHitsPotential[j].at(i) = -999;
	    dEdx[j].at(i) = -999;
	}
	Flag.at(i) = -999;
	Chi2.at(i) = -999;
	NDF.at(i) = -999;
	for (int j=0;j<3;j++)
	{
	    Impact_Point[j].at(i) = -999;
	}
	Length.at(i) = -999;
	for (int j=0;j<nIndices;j++)
	{
	    Index[j].at(i) = -999;
	}
    }
    Daughters.clear();
}
//--------------------------------------------------------
void SetPt(int idx, double pValue){ if (idx < pT.size()) pT.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << pT.size() << std::endl;}
void SetEta(int idx, double pValue){ if (idx < eta.size()) eta.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << eta.size() << std::endl;}
void SetPhi(int idx, double pValue){ if (idx < phi.size()) phi.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << phi.size() << std::endl;}
void SetPx(int idx, double pValue){ if (idx < px.size()) px.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << px.size() << std::endl;}
void SetPy(int idx, double pValue){ if (idx < py.size()) py.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << py.size() << std::endl;}
void SetPz(int idx, double pValue){ if (idx < pz.size()) pz.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << pz.size() << std::endl;}
void SetP(int idx, double pValue){ if (idx < p.size()) p.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << p.size() << std::endl;}
void SetPdgId(int idx, double pValue){ if (idx < PdgId.size()) PdgId.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << PdgId.size() << std::endl;}
void SetMass(int idx, double pValue){ if (idx < Mass.size()) Mass.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << Mass.size() << std::endl;}
void SetRapidity(int idx, double pValue){ if (idx < Rapidity.size()) Rapidity.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << Rapidity.size() << std::endl;}
void SetCharge(int idx, double pValue){ if (idx < Charge.size()) Charge.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << Charge.size() << std::endl;}
void SetFlag(int idx, double pValue){ if (idx < Flag.size()) Flag.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << Flag.size() << std::endl;}
void SetChi2(int idx, double pValue){ if (idx < Chi2.size()) Chi2.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << Chi2.size() << std::endl;}
void SetNDF(int idx, double pValue){ if (idx < NDF.size()) NDF.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << NDF.size() << std::endl;}
void SetLength(int idx, double pValue){ if (idx < Length.size()) Length.at(idx) = pValue; else std::cout << "Wrong index! " << idx << "/" << Length.size() << std::endl;}
//--------------------------------------------------------
void SetNofHits(int idx, int jdx, double pValue)
{
    if (jdx<nSubDetectors)
    {
	if (idx < NofHits[jdx].size()) NofHits[jdx].at(idx) = pValue;
	else std::cout << "Wrong index! " << idx << "/" << NofHits[jdx].size() << std::endl;
    }
    else std::cout << "Wrong index! " << jdx << "/" << nSubDetectors << std::endl;
}
//--------------------------------------------------------
void SetNofHitsPotential(int idx, int jdx, double pValue)
{
    if (jdx<nSubDetectors)
    {
	if (idx < NofHitsPotential[jdx].size()) NofHitsPotential[jdx].at(idx) = pValue;
	else std::cout << "Wrong index! " << idx << "/" << NofHitsPotential[jdx].size() << std::endl;
    }
    else std::cout << "Wrong index! " << jdx << "/" << nSubDetectors << std::endl;
}
//--------------------------------------------------------
void SetdEdx(int idx, int jdx, double pValue)
{
    if (jdx<nSubDetectors)
    {
	if (idx < dEdx[jdx].size()) dEdx[jdx].at(idx) = pValue;
	else std::cout << "Wrong index! " << idx << "/" << dEdx[jdx].size() << std::endl;
    }
    else std::cout << "Wrong index! " << jdx << "/" << nSubDetectors << std::endl;
}
//--------------------------------------------------------
void SetImpactPointComponent(int idx, int jdx, double pValue)
{
    if (jdx<3)
    {
	if (idx < Impact_Point[jdx].size()) Impact_Point[jdx].at(idx) = pValue;
	else std::cout << "Wrong index! " << idx << "/" << Impact_Point[jdx].size() << std::endl;
    }
    else std::cout << "Wrong index! " << jdx << "/" << 3 << std::endl;
}
//--------------------------------------------------------
void SetImpactPoint(int idx, double pX, double pY, double pZ)
{
    if (idx < Impact_Point[0].size()) { Impact_Point[0].at(idx) = pX; Impact_Point[1].at(idx) = pY; Impact_Point[2].at(idx) = pZ; }
    else std::cout << "Wrong index! " << idx << "/" << Impact_Point[0].size() << std::endl;
}
//--------------------------------------------------------
void SetIndex(int idx, int jdx, int pValue)
{
    if (jdx<nIndices)
    {
	if (idx < Index[jdx].size()) Index[jdx].at(idx) = pValue;
	else std::cout << "Wrong index! " << idx << "/" << Index[jdx].size() << std::endl;
    }
    else std::cout << "Wrong index! " << jdx << "/" << nIndices << std::endl;
}
//--------------------------------------------------------
//Getters
double GetPt(int idx){ if (idx < pT.size()) return  pT.at(idx); else std::cout << "Wrong index! " << idx << "/" << pT.size() << std::endl;}
double GetEta(int idx){ if (idx < eta.size()) return  eta.at(idx); else std::cout << "Wrong index! " << idx << "/" << eta.size() << std::endl;}
double GetPhi(int idx){ if (idx < phi.size()) return  phi.at(idx); else std::cout << "Wrong index! " << idx << "/" << phi.size() << std::endl;}
double GetPx(int idx){ if (idx < px.size()) return  px.at(idx); else std::cout << "Wrong index! " << idx << "/" << px.size() << std::endl;}
double GetPy(int idx){ if (idx < py.size()) return  py.at(idx); else std::cout << "Wrong index! " << idx << "/" << py.size() << std::endl;}
double GetPz(int idx){ if (idx < pz.size()) return  pz.at(idx); else std::cout << "Wrong index! " << idx << "/" << pz.size() << std::endl;}
double GetP(int idx){ if (idx < p.size()) return  p.at(idx); else std::cout << "Wrong index! " << idx << "/" << p.size() << std::endl;}
double GetPdgId(int idx){ if (idx < PdgId.size()) return  PdgId.at(idx); else std::cout << "Wrong index! " << idx << "/" << PdgId.size() << std::endl;}
double GetMass(int idx){ if (idx < Mass.size()) return  Mass.at(idx); else std::cout << "Wrong index! " << idx << "/" << Mass.size() << std::endl;}
double GetRapidity(int idx){ if (idx < Rapidity.size()) return  Rapidity.at(idx); else std::cout << "Wrong index! " << idx << "/" << Rapidity.size() << std::endl;}
double GetCharge(int idx){ if (idx < Charge.size()) return  Charge.at(idx); else std::cout << "Wrong index! " << idx << "/" << Charge.size() << std::endl;}
double GetFlag(int idx){ if (idx < Flag.size()) return  Flag.at(idx); else std::cout << "Wrong index! " << idx << "/" << Flag.size() << std::endl;}
double GetChi2(int idx){ if (idx < Chi2.size()) return  Chi2.at(idx); else std::cout << "Wrong index! " << idx << "/" << Chi2.size() << std::endl;}
double GetNDF(int idx){ if (idx < NDF.size()) return  NDF.at(idx); else std::cout << "Wrong index! " << idx << "/" << NDF.size() << std::endl;}
double GetLength(int idx){ if (idx < Length.size()) return  Length.at(idx); else std::cout << "Wrong index! " << idx << "/" << Length.size() << std::endl;}
//--------------------------------------------------------
double GetNofHits(int idx, int jdx)
{
    if (jdx<nSubDetectors)
    {
	if (idx < NofHits[jdx].size()) return  NofHits[jdx].at(idx);
	else std::cout << "Wrong index! " << idx << "/" << NofHits[jdx].size() << std::endl;
    }
    else std::cout << "Wrong index! " << jdx << "/" << nSubDetectors << std::endl;
}
//--------------------------------------------------------
double GetNofHitsPotential(int idx, int jdx)
{
    if (jdx<nSubDetectors)
    {
	if (idx < NofHitsPotential[jdx].size()) return  NofHitsPotential[jdx].at(idx);
	else std::cout << "Wrong index! " << idx << "/" << NofHitsPotential[jdx].size() << std::endl;
    }
    else std::cout << "Wrong index! " << jdx << "/" << nSubDetectors << std::endl;
}
//--------------------------------------------------------
double GetdEdx(int idx, int jdx)
{
    if (jdx<nSubDetectors)
    {
	if (idx < dEdx[jdx].size()) return  dEdx[jdx].at(idx);
	else std::cout << "Wrong index! " << idx << "/" << dEdx[jdx].size() << std::endl;
    }
    else std::cout << "Wrong index! " << jdx << "/" << nSubDetectors << std::endl;
}
//--------------------------------------------------------
double GetImpactPointComponent(int idx, int jdx)
{
    if (jdx<3)
    {
	if (idx < Impact_Point[jdx].size()) return  Impact_Point[jdx].at(idx);
	else std::cout << "Wrong index! " << idx << "/" << Impact_Point[jdx].size() << std::endl;
    }
    else std::cout << "Wrong index! " << jdx << "/" << 3 << std::endl;
}
//--------------------------------------------------------
int GetIndex(int idx, int jdx)
{
    if (jdx<nIndices)
    {
	if (idx < Index[jdx].size()) return  Index[jdx].at(idx);
	else std::cout << "Wrong index! " << idx << "/" << Index[jdx].size() << std::endl;
    }
    else std::cout << "Wrong index! " << jdx << "/" << nIndices << std::endl;
}
//--------------------------------------------------------
//Daughters

void AddDaughterTrack(DataTreeRecoTrack* track){ Daughters.push_back(track); }
DataTreeRecoTrack* GetDaughterTrack(int idx){ if (idx < Daughters.size()) return Daughters.at(idx); else std::cout << "Wrong index! " << idx << "/" << Daughters.size() << std::endl; }


private:    
    void SetId(int idx){id = idx;}
  
    int id;
    
    int nPoints;
    
    std::vector <double> pT;
    std::vector <double> eta;
    std::vector <double> phi;
    std::vector <double> px;
    std::vector <double> py;
    std::vector <double> pz;
    std::vector <double> p;
    std::vector <double> PdgId;
    std::vector <double> Mass;
    std::vector <double> Rapidity;
    std::vector <double> Charge;
    std::vector <double> NofHits[nSubDetectors];
    std::vector <double> NofHitsPotential[nSubDetectors];
    std::vector <double> dEdx[nSubDetectors];
    std::vector <double> Flag;
    std::vector <double> Chi2;
    std::vector <double> NDF;
    std::vector <double> Impact_Point[3];
    std::vector <double> Length;
    std::vector <int> Index[nIndices];

    std::vector <DataTreeRecoTrack*> Daughters;
    ClassDefNV(DataTreeRecoTrack, 1)
};

#endif