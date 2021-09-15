#ifndef iDREAMPMTSD_h
#define iDREAMPMTSD_h 1

#include "G4VSensitiveDetector.hh"
class G4Step;
class iDREAMRunAction;

class iDREAMPMTSD: public G4VSensitiveDetector 
{
  public:
    iDREAMPMTSD(G4String);
    ~iDREAMPMTSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);
    
  private:
    iDREAMRunAction* runAction;
    G4double totalPhotons;
    G4int lambda;
};

#endif