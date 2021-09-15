#include "iDREAMEventAction.hh"
#include "iDREAMRunAction.hh"
#include "iDREAMStackingAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

iDREAMEventAction::iDREAMEventAction() {}

iDREAMEventAction::~iDREAMEventAction() {}

void iDREAMEventAction::EndOfEventAction(const G4Event* event)
{
  // в конце каждого события
  iDREAMRunAction* runAction = (iDREAMRunAction*) G4RunManager::GetRunManager()->GetUserRunAction();
  iDREAMStackingAction* stackingAction = (iDREAMStackingAction*) G4RunManager::GetRunManager()->GetUserStackingAction();
  
  // набираем гистограмму испущенных фотонов
  runAction->FillEmitedHist(stackingAction->emitedPhotons);
  // и гамма квантов
  runAction->FillEmitedEnergy(stackingAction->gammaEnergy);
  // отображаем прогресс моделирования
  runAction->DisplayProgress(event->GetEventID());

  //G4cout << "iDREAMEventAction::EndOfEventAction" << G4endl;
}

