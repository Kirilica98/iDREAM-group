#include "iDREAMRunAction.hh"
#include "iDREAMHistoManager.h"

#include "G4Run.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

//#include <ctime>

iDREAMRunAction::iDREAMRunAction() {}

iDREAMRunAction::~iDREAMRunAction() {}

void iDREAMRunAction::BeginOfRunAction(const G4Run* run)
{
  start_time =  clock(); // начальное время

  eventsNumber = run->GetNumberOfEventToBeProcessed();
  printModulo = eventsNumber/100;
  if (printModulo == 0)
    printModulo = 1;

  // создаем гистограмму
  // от 0 до 1000, с 1000 каналов
  histEmited = new iDREAMHistoManager(10000);   //0 - 4 294 967 295
  histEmitedeeeee = new iDREAMHistoManager(10000);
  histDetected = new iDREAMHistoManager(10000);
  histArrival = new iDREAMHistoManager(50000000);
  k=0;
  G4cout << "iDREAMRunAction::BeginOfRunAction" << G4endl;
}

void iDREAMRunAction::DisplayProgress(G4int eventID)
{
  //if (eventID%printModulo == 0)
    //G4cout << "Progress: " << 100*eventID/eventsNumber << "%\r" << std::flush;
}

void iDREAMRunAction::FillEmitedHist(G4double photons)
{
  if (photons > 0)
    histEmited->fill(photons);
  k++;
  G4cout << "NomberOfEvent = " << k << G4endl;
}

void iDREAMRunAction::FillEmitedEnergy(G4double ph)
{
    histEmitedeeeee->fill(ph);
  //k++;
  //G4cout << "NomberOfEvent = " << k << G4endl;
}

void iDREAMRunAction::FillDetectedHist(G4double photons)
{
  if (photons > 0)
    histDetected->fill(photons);
}

void iDREAMRunAction::FillArrivalHist(G4double time)
{
  histArrival->fill(time/ns);
}

void iDREAMRunAction::EndOfRunAction(const G4Run*)
{
  histEmited->save("spectrum-emited.dat", "\"#photons\", N");
  histEmitedeeeee->save("spectrum-Eeeeeeeeeeenergy.dat", "\"#photons\", N");
  histDetected->save("spectrum-detected.dat", "\"#photons\", N");
  histArrival->save("spectrum-arrival.dat", "\"time, ns\", N");
  G4cout << "iDREAMRunAction::EndOfRunAction" << G4endl;
  end_time = clock(); // конечное время
  search_time = end_time - start_time; // искомое время
  G4cout << "search_time = " << search_time/1000000.0 << G4endl;
}

