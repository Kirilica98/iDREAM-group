#include "iDREAMStackingAction.hh"

#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

iDREAMStackingAction::iDREAMStackingAction() {}

iDREAMStackingAction::~iDREAMStackingAction() {}

void iDREAMStackingAction::PrepareNewEvent()
{
  emitedPhotons = 0;
  gammaEnergy = 0;
}

G4ClassificationOfNewTrack iDREAMStackingAction::ClassifyNewTrack(const G4Track* track)
{
  G4String particleName = track->GetDefinition()->GetParticleName();

  if (particleName == "opticalphoton") {
//    emitedPhotons++;
	emitedPhotons+=(track->GetTotalEnergy())*MeV;
}

  if (particleName == "Gd158")
  	G4cout << "Gd158" << G4endl;

  if (particleName == "Gd156")
  	G4cout << "Gd156" << G4endl;

  if (particleName == "deuteron")
  	G4cout << "deuteron" << G4endl;

  if (particleName == "gamma") {
  	gammaEnergy += (track->GetTotalEnergy())*MeV;
  	G4cout << "Gamma energy = " << (track->GetTotalEnergy())*MeV << G4endl;
  	}

//G4cout << "iDREAMStackingAction::ClassifyNewTrack" << G4endl;
  
  return fWaiting;
}
 

