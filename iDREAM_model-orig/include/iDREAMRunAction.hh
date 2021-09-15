#ifndef iDREAMRunAction_h
#define iDREAMRunAction_h 1

class iDREAMHistoManager;

//#include <fstream>
//using namespace std;

#include "G4UserRunAction.hh"
#include "globals.hh"
class G4Run;

class iDREAMRunAction: public G4UserRunAction
{
  public:
    iDREAMRunAction();
   ~iDREAMRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
    
    void DisplayProgress(G4int);
    void FillEmitedHist(G4double);
    void FillEmitedEnergy(G4double);
    void FillDetectedHist(G4double);
    void FillArrivalHist(G4double);

  private:
    G4int eventsNumber;
    G4int printModulo;
    iDREAMHistoManager* histEmited;
    iDREAMHistoManager* histEmitedeeeee;
    iDREAMHistoManager* histDetected;
    iDREAMHistoManager* histArrival;
    //ofstream Emited;
    //ofstream Detected;
    //ofstream Arrival;
    int k;
    unsigned int start_time;
    unsigned int end_time;
    unsigned int search_time;
};

#endif
