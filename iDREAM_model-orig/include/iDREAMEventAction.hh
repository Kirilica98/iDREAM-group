#ifndef iDREAMEventAction_h
#define iDREAMEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
class G4Event;

class iDREAMEventAction: public G4UserEventAction
{
  public:
    iDREAMEventAction();
   ~iDREAMEventAction();

    void EndOfEventAction(const G4Event*);
};

#endif

