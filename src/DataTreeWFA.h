#ifndef DataTreeWFA_H
#define DataTreeWFA_H 1

#include <vector>
#include <iostream>
#include "TClonesArray.h"
#include "TObject.h"

class DataTreeWFA : public TObject
{
    
public:
  
    DataTreeWFA(int idx = 0);
    ~DataTreeWFA();
        
    int GetId(){return id;}
    
    int GetNHits(int idx){return NHits[idx];}
    double GetTime(int idx, int jdx){return TimeWFA[idx][jdx];}

    void SetTime(int idx, int jdx, double fValue){TimeWFA[idx][jdx] = fValue;}
    void SetNHits(int idx, int fValue){NHits[idx] = fValue;}
    
    
private:    
    void SetId(int idx){id = idx;}
    
    double TimeWFA[6][2000];
    int    NHits[6];
    int id;
    
    ClassDefNV(DataTreeWFA, 1)
};

#endif