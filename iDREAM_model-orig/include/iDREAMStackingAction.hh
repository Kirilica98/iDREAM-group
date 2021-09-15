#ifndef iDREAMStackingAction_H
#define iDREAMStackingAction_H 1

#include "G4UserStackingAction.hh"
#include "globals.hh"
class G4Track;

class iDREAMStackingAction: public G4UserStackingAction
{
  public:
    iDREAMStackingAction();
    ~iDREAMStackingAction();

  public:
    void PrepareNewEvent();
    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);
    
    G4double emitedPhotons;
    G4double gammaEnergy;
};

#endif

