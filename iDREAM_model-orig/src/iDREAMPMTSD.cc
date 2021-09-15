#include "iDREAMPMTSD.hh"
#include "iDREAMRunAction.hh"

#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

iDREAMPMTSD::iDREAMPMTSD(G4String name): G4VSensitiveDetector(name)
{
  // получаем указатель на класс iDREAMRunAction
  // мы будев вызывать его метод iDREAMRunAction::FillHist
  // для заполнения гистограммы спектра поглощенной энергии
  runAction = (iDREAMRunAction*) G4RunManager::GetRunManager()->GetUserRunAction();
  //G4cout << "iDREAMPMTSD::iDREAMPMTSD---------------------iDREAMPMTSD::iDREAMPMTSD------------------------iDREAMPMTSD::iDREAMPMTSD-------" << G4endl;
}

iDREAMPMTSD::~iDREAMPMTSD() {}

void iDREAMPMTSD::Initialize(G4HCofThisEvent*)
{
  // в начале события сбрасываем энергию, поглощенную детектором
  totalPhotons = 0;
  //G4cout << "iDREAMPMTSD::Initialize(G4HCofThisEvent*)---------------------iDREAMPMTSD::Initialize(G4HCofThisEvent*)-------------" << G4endl;
}

G4bool iDREAMPMTSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  G4double Energy[] = {210., 220., 230., 240., 250., 260., 270., 280., 290., 300., 310., 320., 330., 340., 350., 360., 370., 380., 390., 400.,
                       410., 420., 430., 440., 450., 460., 470., 480., 490., 500., 510., 520., 530., 540., 550., 560., 570., 580., 590., 600.,
                       610., 620., 630., 640., 650., 660., 670., 680., 690., 700.};
  G4double Probability[] = {0., 0.00169, 0.00291, 0.00094, 0.0007, 0.00108, 0.00434, 0.0193, 0.058, 0.112, 0.176, 0.223, 0.256, 0.278, 0.289,
                            0.312, 0.3085, 0.305, 0.294, 0.283, 0.274, 0.265, 0.254, 0.243, 0.2275, 0.212, 0.198, 0.184, 0.1735, 0.163, 0.1400,
                            0.117, 0.10655, 0.0661, 0.0544, 0.0427, 0.0336, 0.0245, 0.01775, 0.011, 0.00722, 0.00343, 0.00172, 0.00077, 0.00048,
                            0.00019, 0.000155, 0.00012, 0.00010, 0.00008};

  G4Track* track = step->GetTrack();
  G4String particleName = track->GetDefinition()->GetParticleName();
  //G4cout << "iDREAMPMTSD::ProcessHits----------------------------------------------------------------------------------------------------" << totalPhotons << G4endl;
  
  if (particleName == "opticalphoton") {
    //G4double fff=(h_Planck*c_light)/(track->GetTotalEnergy())*1000000;
    //G4cout << "WaveLangth = " << fff << G4endl;
    for (G4int i=0;i<37;i++) {
      if ((Energy[i]<=((h_Planck*c_light)/(track->GetTotalEnergy())*1000000))&&(((h_Planck*c_light)/(track->GetTotalEnergy())*1000000)<Energy[i+1])) {
        lambda=i+1;
        //G4cout << "lambda = " << lambda << G4endl;
        break;
      }
    }

    if (((h_Planck*c_light)/(track->GetTotalEnergy())*1000000)<Energy[0]) {
      lambda=0;
    }
    if (((h_Planck*c_light)/(track->GetTotalEnergy())*1000000)>Energy[36]) {
      lambda=0;
    }

    G4double lambda_rand = G4UniformRand();
    if (lambda_rand<Probability[lambda]) {
      totalPhotons++;
    }
    
    G4double timea = step->GetPreStepPoint()->GetGlobalTime();
    runAction->FillArrivalHist(timea);
    //G4cout << "Time: " << timea << G4endl;
    //G4cout << "totalPhotons: " << totalPhotons << G4endl;
  
    // уничтожаем частицу
    track->SetTrackStatus(fStopAndKill);
  }

  return true;
}

void iDREAMPMTSD::EndOfEvent(G4HCofThisEvent*)
{
  // сохраняем энергию, накопленную за событие в детекторе
  // в густограмму
  runAction->FillDetectedHist(totalPhotons);
  //G4cout << "totalPhotons: " << totalPhotons << G4endl;
  //G4cout << "iDREAMPMTSD::EndOfEvent" << G4endl;
}
