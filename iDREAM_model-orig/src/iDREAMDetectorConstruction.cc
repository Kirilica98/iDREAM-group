//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B1DetectorConstruction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "iDREAMDetectorConstruction.hh"
#include "iDREAMPMTSD.hh"

#include "G4RunManager.hh"
#include "G4RunManager.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "globals.hh"
#include "G4UImanager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>

#include <iostream>
#include <fstream>
using namespace std;
#define SIZE(x) sizeof(x)/sizeof(*x)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

iDREAMDetectorConstruction::iDREAMDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

iDREAMDetectorConstruction::~iDREAMDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//void iDREAMDDetectorConstruction::ConstructMaterials(){
  
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* iDREAMDetectorConstruction::Construct()
{  
  // Get nist material manager
  //G4NistManager* nist = G4NistManager::Instance();
  G4double  density;  
  G4int   ncomponents, natoms;
  G4double  fractionmass;
  G4double mc;
  G4double  temperature = 290.0;
  
  G4NistManager* nistManager = G4NistManager::Instance();

  //***Elements
  G4Element* fH = nistManager->FindOrBuildElement(1);
  G4Element* fC = nistManager->FindOrBuildElement(6);
  G4Element* fO = nistManager->FindOrBuildElement(8);
  G4Element* fN = nistManager->FindOrBuildElement(7);
  G4Element* fFe = nistManager->FindOrBuildElement(26);
  G4Element* fSi = nistManager->FindOrBuildElement(14);
  G4Element* fCr = nistManager->FindOrBuildElement(24);
  G4Element* fMn = nistManager->FindOrBuildElement(12);
  G4Element* fNi = nistManager->FindOrBuildElement(28);  
  G4Element* fGd = nistManager->FindOrBuildElement(64);

  //***Materials
  //Vacuum
  G4Material* fVacuum = nistManager->FindOrBuildMaterial("G4_Galactic");
  //Air
  G4Material* fAir = nistManager->FindOrBuildMaterial("G4_AIR");
  //N2
  density = 27.*mg/cm3;
  G4double pressure = 5.*atmosphere;
  //G4double temperature = 325.*kelvin;
  G4Material* fN2 = new G4Material("N2", density, ncomponents=1, kStateGas, temperature, pressure);
  fN2->AddElement(fN, natoms=2);
  //Glass
  G4Material* fGlass = new G4Material("Glass", density=1.032*g/cm3, 2);
  fGlass->AddElement(fC, 91.533*perCent);
  fGlass->AddElement(fH, 8.467*perCent);
  //StainlessSteel
  G4Material* fStainlessSteel = new G4Material("StainlessSteel", density=8.06*g/cm3, 6);
  fStainlessSteel->AddElement(fC, 0.001);
  fStainlessSteel->AddElement(fSi, 0.007);
  fStainlessSteel->AddElement(fCr, 0.18);
  fStainlessSteel->AddElement(fMn, 0.01);
  fStainlessSteel->AddElement(fFe, 0.712);
  fStainlessSteel->AddElement(fNi, 0.09);
  //Acrylic
  G4Material* fAcrylic = new G4Material("Acrylic", density=1.19*g/cm3, ncomponents=3, kStateSolid, temperature);
  fAcrylic->AddElement(fH, natoms=6);
  fAcrylic->AddElement(fC, natoms=4);
  fAcrylic->AddElement(fO, natoms=2);  
  //LAB
  //mol. weight = 120.19
  G4Material* fLAB = new G4Material("LAB", density=0.86*g/cm3, ncomponents=2, kStateLiquid, 290.0);
  fLAB->AddElement(fC, natoms=8);
  fLAB->AddElement(fH, natoms=9);
  //PPO
  // Mol. Weight = 221.3
  G4Material* fPPO = new G4Material("PPO", density=0.3*g/cm3, ncomponents=4, kStateSolid, 290.0);
  fPPO->AddElement(fC, natoms=15);
  fPPO->AddElement(fH, natoms=11);
  fPPO->AddElement(fN, natoms=1);
  fPPO->AddElement(fO, natoms=1);
  //bisMSB
  // Mol. Weight = 310.4
  G4Material* fbisMSB = new G4Material("bisMSB", density=0.45*g/cm3, ncomponents=2, kStateSolid, 290.0);
  fbisMSB->AddElement(fC, natoms=24);
  fbisMSB->AddElement(fH, natoms=22);
  //GadoliniumDPM
  G4Material* fGadoliniumDPM = new G4Material("GadoliniumDPM", density=1.3*g/cm3, 4, kStateSolid, temperature);
  fGadoliniumDPM->AddElement(fGd,natoms=1);
  fGadoliniumDPM->AddElement(fH,natoms=57);
  fGadoliniumDPM->AddElement(fO,natoms=6);
  fGadoliniumDPM->AddElement(fC,natoms=33);

     
  //***Material properties tables

  const G4int num = 3;

  G4double glass_Energy[num] = {7.0*eV, 7.07*eV, 7.14*eV};
  G4double glass_RIND[num] = {1.49, 1.49, 1.49};
  G4double glass_AbsLength[num]={420.*cm, 420.*cm, 420.*cm};
  G4MaterialPropertiesTable *glass_mt = new G4MaterialPropertiesTable();
  glass_mt->AddProperty("ABSLENGTH", glass_Energy, glass_AbsLength, num);
  glass_mt->AddProperty("RINDEX", glass_Energy, glass_RIND, num);
  fGlass->SetMaterialPropertiesTable(glass_mt);

  G4double vacuum_Energy[num]={2.0*eV, 7.0*eV, 7.14*eV};
  G4double vacuum_RIND[num]={1., 1., 1.};
  G4MaterialPropertiesTable *vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND, num);
  fVacuum->SetMaterialPropertiesTable(vacuum_mt);
  fAir->SetMaterialPropertiesTable(vacuum_mt);//Give air the same rindex

  /*G4MaterialPropertiesTable *Acrylic_mt = new G4MaterialPropertiesTable();
  FillMaterialPropertiesTable("Acrylic",Acrylic_mt,"ABSLENGTH",1.0);
  FillMaterialPropertiesTable("Acrylic",Acrylic_mt,"RINDEX",1.0);
  fAcrylic->SetMaterialPropertiesTable(Acrylic_mt);
  *//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //LAB_mt
  G4double LAB_MolW = (fH->GetA()*9 + fC->GetA()*8)/g/mole;
  G4double bisMSB_MolW = (fH->GetA()*22 + fC->GetA()*24)/g/mole;
  G4double PPO_MolW = (fH->GetA()*11 + fC->GetA()*15 + fN->GetA()*1 + fO->GetA()*1)/g/mole;
  G4double GadoliniumDPM_MolW = (fH->GetA()*57 + fC->GetA()*33 + fO->GetA()*6 + fGd->GetA())/g/mole;
  const G4int LABALnum = 261;
  G4MaterialPropertiesTable *LAB_mt = new G4MaterialPropertiesTable();
  G4double LABAL_Energy[LABALnum] = {340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362,
                                     363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385,
                                     386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408,
                                     409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431,
                                     432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454,
                                     455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477,
                                     478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500,
                                     501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523,
                                     524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546,
                                     547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569,
                                     570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592,
                                     593, 594, 595, 596, 597, 598, 599, 600}; //lambda

  G4double LABAbsLength[LABALnum] = {0.698223, 0.734847, 0.74621,  0.738596, 0.711958, 0.676471, 0.645311, 0.616021, 0.587679, 0.563287,
                                     0.553947, 0.566964, 0.60235,  0.656034, 0.723824, 0.791065, 0.84329,  0.873832, 0.882712, 0.868589,
                                     0.828806, 0.76192,  0.686089, 0.623091, 0.577519, 0.535505, 0.486877, 0.440015, 0.410486, 0.409326,
                                     0.440015, 0.499764, 0.584515, 0.686089, 0.798335, 0.89916,  0.978141, 1.02428,  1.03898,  1.01234,
                                     0.94207,  0.849891, 0.765951, 0.710793, 0.677527, 0.634933, 0.569193, 0.492956, 0.439125, 0.425779,
                                     0.458117, 0.54219,  0.681781, 0.873832, 1.11358,  1.39197,  1.70312,  2.05827,  2.41275,  2.71434,
                                     2.93442,  3.12442,  3.39293,  3.77647,  4.34468,  5.08542,  5.93055,  6.7795,   7.75387,  8.8451,
                                     10.1899,  11.7377,  13.6829,  16.0138,  18.9648,  22.2373,  26.0057,  30.3702,  36.4341,  41.6789,
                                     46.2507,  52.9627,  57.9832,  61.341,   69.7102,  73.6092,  78.1105,  81.6343,  83.3579,  88.6315,
                                     92.0115,  95.4493,  97.8141,  96.9407,  102.67,   110.789,  109.948,  108.846,  114.288,  109.948,
                                     118.015,  122.336,  116.122,  124.084,  124.084,  129.254,  135.294,  137.871,  134.874,  149.242,
                                     148.223,  148.223,  150.275,  154.005,  163.885,  160.256,  158.502,  160.256,  172.339,  175.828,
                                     177.99,   177.263,  180.205,  181.713,  182.477,  203.894,  200.136,  197.407,  171.658,  183.247,
                                     182.477,  192.166,  196.513,  188.824,  227.379,  222.715,  228.576,  243.986,  243.986,  268.083,
                                     260.057,  308.01,   297.462,  303.702,  301.593,  312.442,  336.662,  358.921,  355.979,  353.085,
                                     324.1,    353.085,  368.046,  324.1,    368.046,  368.046,  413.614,  334.073,  314.706,  324.1,
                                     368.046,  368.046,  371.192,  405.883,  398.435,  384.331,  452.39,   549.74,   457.152,  493.516,
                                     405.883,  504.994,  510.935,  564.019,  542.868,  620.421,  668.145,  579.059,  748.784,  700.475,
                                     748.784,  761.92,   761.92,   700.475,  611.682,  700.475,  658.022,  523.246,  499.189,  457.152,
                                     398.435,  387.763,  314.706,  314.706,  374.392,  429.995,  477.247,  804.249,  1142.88,  1173.77,
                                     1608.5,   1670.36,  2895.3,   2895.3,   3102.1,   1206.37,  1447.65,  1085.74,  904.78,   987.033,
                                     904.78,   804.249,  603.187,  668.145,  549.74,   447.726,  425.779,  380.96,   321.7,    358.921,
                                     334.073,  329.011,  347.436,  339.293,  377.647,  339.293,  308.01,   355.979,  297.462,  231.008,
                                     226.195,  219.341,  266.438,  329.011,  380.96,   447.726,  482.549,  564.019,  1034.03,  1240.84,
                                     21714.7,  2068.07,  21714.7,  43429.4,  43429.4,  43429.4,  43429.4,  43429.4,  43429.4,  8685.89,
                                     8685.89,  2171.47,  2068.07,  851.558,  775.526,  443.158,  611.682,  421.645,  334.073,  291.473,
                                     227.379};
  
  for(int i=0; i < LABALnum; ++i) {
    LABAL_Energy[i]=(h_Planck*c_light)/(LABAL_Energy[i]*nm);//////////QUESTION
    LABAbsLength[i]=LABAbsLength[i]*m;
  }
  
  const G4int LABRnum = 2;
  G4double LABR_Energy[LABRnum] = {200, 700};
  G4double LABRindex[LABRnum] = {1.5048, 1.5048};

  for(int i=0; i < LABRnum; ++i) {
     LABR_Energy[i]=(h_Planck*c_light)/(LABR_Energy[i]*nm);
  }

  LAB_mt->AddProperty("ABSLENGTH", LABAL_Energy, LABAbsLength, SIZE(LABAL_Energy));
  LAB_mt->AddProperty("RINDEX", LABR_Energy, LABRindex, SIZE(LABR_Energy));
  LAB_mt->AddConstProperty("LAB_MolW",LAB_MolW);
  fLAB->SetMaterialPropertiesTable(LAB_mt);
  
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Target Scintiilator  
  //-------------------------------------------------------------------
  G4double  LABFractionVolume = 1.0;    // 1.0 volume 
  G4double  PPOConcetration = 3E-3*g/cm3;   // 3 g/L
  G4double  bisMSBConcetration = 2.0E-5*g/cm3;  // 20 mg/L
  G4double  GadoliniumConcetration = 1.0E-3*g/cm3;  // 1 g/L (metallic Gd)
  G4double  GadoliniumDPMConcetration = GadoliniumConcetration*GadoliniumDPM_MolW/(fGd->GetA()/g/mole);
  
  density = LABFractionVolume*fLAB->GetDensity() + PPOConcetration + bisMSBConcetration + GadoliniumDPMConcetration;

  G4Material* fScint_LAB_Gd = new G4Material("Scint_LAB_Gd", density, ncomponents=4, kStateLiquid, temperature);
  G4MaterialPropertiesTable *Scint_LAB_Gd_mt = new G4MaterialPropertiesTable();

  fScint_LAB_Gd->AddMaterial(fLAB, fractionmass = LABFractionVolume*fLAB->GetDensity()/density);  
  mc = (density*fractionmass*1.0E+6/LAB_MolW)/g/mole;
  Scint_LAB_Gd_mt->AddConstProperty("LAB_MC",mc);
  
  fScint_LAB_Gd->AddMaterial(fPPO, fractionmass=PPOConcetration/density);
  mc = (density*fractionmass*1.0E+6/PPO_MolW)/g/mole;
  Scint_LAB_Gd_mt->AddConstProperty("PPO_MC",mc);
  
  fScint_LAB_Gd->AddMaterial(fbisMSB, fractionmass=bisMSBConcetration/density);
  mc = (density*fractionmass*1.0E+6/bisMSB_MolW)/g/mole;
  Scint_LAB_Gd_mt->AddConstProperty("bisMSB_MC",mc); 

  fScint_LAB_Gd->AddMaterial(fGadoliniumDPM, fractionmass=GadoliniumDPMConcetration/density);
  mc = (density*fractionmass*1.0E+6/GadoliniumDPM_MolW)/g/mole;
  Scint_LAB_Gd_mt->AddConstProperty("GadoliniumDPM_MC",mc); 

  
  G4double ScintGdAbsLength[LABALnum] = {0.693014, 0.729365, 0.740644, 0.733087, 0.706648, 0.671426, 0.6405, 0.611429, 0.5833, 0.559092, 0.549823, 0.562745, 0.597869, 
                                         0.651155, 0.718441, 0.785182, 0.837019, 0.867335, 0.87615, 0.862133, 0.822648, 0.756261, 0.680996, 0.618468, 0.573236, 0.531536,
                                         0.483272, 0.43676, 0.407452, 0.406302, 0.436773, 0.496083, 0.580211, 0.681039, 0.792464, 0.892559, 0.970979, 1.01681, 1.03144,
                                         1.00506, 0.935377, 0.843961, 0.760736, 0.706089, 0.673183, 0.631028, 0.565912, 0.490392, 0.437121, 0.424035, 0.456287, 0.539905, 
                                         0.678654, 0.869422, 1.10757, 1.38405, 1.69311, 2.04575, 2.39782, 2.69732, 2.91591, 3.1047, 3.37129, 3.75224, 4.31647, 5.05182, 
                                         5.89086, 6.73382, 7.70121, 8.78454, 10.1196, 11.6561, 13.5872, 15.9012, 18.8305, 22.0789, 25.8195, 30.1518, 36.1708, 41.3769, 45.9148, 
                                         52.5771, 57.5607, 60.8935, 69.2004, 73.071, 77.5386, 81.0366, 82.7477, 87.9823, 91.3376, 94.7495, 97.0969, 96.2304, 101.917, 109.976, 
                                         109.142, 108.048, 113.45, 109.142, 117.149, 121.439, 115.271, 123.174, 123.175, 128.306, 134.302, 136.86, 133.885, 148.146, 147.137, 
                                         147.137, 149.173, 152.877, 162.683, 159.081, 157.341, 159.083, 171.075, 174.539, 176.685, 175.964, 178.886, 180.385, 181.143, 202.4, 
                                         198.67, 195.965, 170.41, 181.912, 181.146, 190.767, 195.088, 187.453, 225.723, 221.095, 226.919, 242.206, 242.214, 266.131, 258.166,
                                         305.774, 295.328, 301.508, 299.427, 310.179, 334.226, 356.314, 353.412, 350.534, 321.794, 350.533, 365.389, 321.815, 365.442, 365.966, 
                                         410.763, 333.54, 312.477, 321.981, 365.598, 367.261, 370.381, 404.815, 397.426, 383.427, 450.979, 547.6, 455.704, 491.797, 404.817, 
                                         503.189, 509.087, 561.776, 540.779, 617.758, 665.122, 576.704, 745.162, 697.217, 745.163, 758.203, 758.204, 697.217, 609.086, 697.22, 
                                         655.084, 521.313, 497.437, 455.712, 397.434, 386.847, 314.333, 314.332, 373.571, 428.762, 475.657, 800.221, 1136.33, 1166.99, 1598.48, 
                                         1659.87, 2875.68, 2875.68, 3080.94, 1199.35, 1438.83, 1079.62, 900.009, 981.65, 900.011, 800.233, 600.674, 665.142, 547.622, 446.373, 
                                         424.582, 380.103, 321.283, 358.233, 333.581, 328.537, 346.825, 338.743, 376.815, 338.742, 307.703, 355.31, 297.229, 231.275, 226.498, 
                                         219.693, 266.436, 328.537, 380.112, 446.384, 480.941, 561.803, 1028.31, 1233.57, 21554.7, 2054.63, 21554.7, 43107.4, 43107.4, 43107.4, 
                                         43107.4, 43107.4, 43107.4, 8623.07, 8623.08, 2157.27, 2054.63, 847.206, 771.739, 441.837, 609.104, 420.488, 333.578, 291.292, 227.672};

  G4double ScintGdES_Energy[LABALnum] = {339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362,
                                     363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385,
                                     386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408,
                                     409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431,
                                     432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454,
                                     455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477,
                                     478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500,
                                     501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523,
                                     524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546,
                                     547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569,
                                     570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587, 588, 589, 590, 591, 592,
                                     593, 594, 595, 596, 597, 598, 599}; //lambda

  G4double ScintGdEmissionSpectra[LABALnum] = {0, 607586, 740041, 857636, 942503, 1.01562e+06, 1.05831e+06, 1.08089e+06, 1.1029e+06, 1.11071e+06, 1.11319e+06, 1.13253e+06, 
                                               1.17763e+06, 1.22859e+06, 1.30412e+06, 1.39074e+06, 1.48388e+06, 1.57634e+06, 1.66983e+06, 1.76751e+06, 1.86507e+06, 1.93489e+06, 
                                               1.99217e+06, 2.02235e+06, 2.04047e+06, 2.02822e+06, 2.02252e+06, 1.99239e+06, 1.98439e+06, 1.9292e+06, 1.91338e+06, 1.89833e+06, 
                                               1.85941e+06, 1.85233e+06, 1.83476e+06, 1.81543e+06, 1.81372e+06, 1.80177e+06, 1.78167e+06, 1.77201e+06, 1.75171e+06, 1.73136e+06, 
                                               1.70003e+06, 1.67753e+06, 1.63343e+06, 1.60316e+06, 1.55068e+06, 1.50379e+06, 1.46437e+06, 1.42094e+06, 1.37559e+06, 1.32982e+06, 
                                               1.28625e+06, 1.23931e+06, 1.20216e+06, 1.15769e+06, 1.11433e+06, 1.08245e+06, 1.04395e+06, 1.00837e+06, 972317, 937816, 910793, 
                                               890576, 865925, 848339, 824220, 795546, 774618, 755067, 734443, 708016, 680963, 661675, 633376, 617302, 590476, 577408, 554976, 
                                               528823, 504544, 476371, 456860, 434192, 417048, 403378, 388564, 371550, 358736, 345707, 334269, 322723, 309707, 294853, 282567, 
                                               271564, 260188, 247916, 236694, 224614, 214970, 201935, 193197, 183785, 173908, 167760, 159448, 153222, 146024, 138101, 134631, 
                                               125070, 121322, 114858, 110659, 109035, 103356, 97221.2, 95534.9, 91247.8, 86738.4, 81834.6, 76895, 71944.6, 70422.5, 67401, 64309.8,
                                               62108.5, 58729.8, 57596.1, 53235.1, 51029.2, 49852.3, 47677.8, 45249.7, 42486.4, 41541, 39209.7, 36883.4, 35937.9, 32990.5, 33171.1, 
                                               31344.7, 29670.6, 27467.8, 26933.9, 25045.7, 23249.7, 23020.6, 21699.1, 21531.4, 20015.2, 19929.6, 18086.1, 17720.2, 16364.3, 15692.4, 
                                               14553.7, 14556.3, 13390.1, 13227.6, 13271.6, 12903.1, 11811.8, 10476.5, 10712.9, 9696.32, 9043.16, 9943.89, 9286.7, 8637.96, 7940.41, 
                                               0.298013, 0.291391, 0.284768, 0.271523, 0.264901, 0.258278, 0.251656, 0.245033, 0.238411, 0.225166, 0.218543, 0.211921, 0.205298, 
                                               0.198675, 0.192053, 0.18543, 0.178808, 0.172185, 0.165563, 0.15894, 0.152318, 0.145695, 0.139073, 0.13245, 0.125828, 0.125828, 
                                               0.119205, 0.112583, 0.112583, 0.10596, 0.0993377, 0.0993377, 0.0927152, 0.0927152, 0.0860927, 0.0860927, 0.0794702, 0.0794702, 
                                               0.0794702, 0.0728477, 0.0728477, 0.0728477, 0.0662252, 0.0662252, 0.0662252, 0.0662252, 0.0662252, 0.0596026, 0.0596026, 0.0596026, 
                                               0.0596026, 0.0596026, 0.0529801, 0.0529801, 0.0529801, 0.0529801, 0.0463576, 0.0463576, 0.0463576, 0.0463576, 0.0397351, 0.0397351, 
                                               0.0397351, 0.0397351, 0.0397351, 0.0397351, 0.0397351, 0.0331126, 0.0331126, 0.0331126, 0.0331126, 0.0331126, 0.0331126, 0.0331126, 
                                               0.0331126, 0.0264901, 0.0264901, 0.0264901, 0.0264901, 0.0264901, 0.0264901, 0.0264901, 0.0264901, 0.0264901, 0.0264901, 0.0264901, 
                                               0.0264901, 0.0264901, 0.0264901};
  
  for(int i=0; i < LABALnum; ++i) {
    ScintGdES_Energy[i]=(h_Planck*c_light)/(ScintGdES_Energy[i]*nm);//////////QUESTION
    ScintGdAbsLength[i]=ScintGdAbsLength[i]*m;
  }
    
  Scint_LAB_Gd_mt->AddProperty("ABSLENGTH", LABAL_Energy, ScintGdAbsLength, SIZE(LABAL_Energy));
  Scint_LAB_Gd_mt->AddProperty("RINDEX", LABR_Energy, LABRindex, SIZE(LABR_Energy));
  Scint_LAB_Gd_mt->AddProperty("FASTCOMPONENT", ScintGdES_Energy, ScintGdEmissionSpectra, LABALnum);
//  Scint_LAB_Gd_mt->AddProperty("SLOWCOMPONENT", Scnt_PP, Scnt_SLOW, LABALnum);
  Scint_LAB_Gd_mt->AddConstProperty("SCINTILLATIONYIELD", 11500./MeV);
  Scint_LAB_Gd_mt->AddConstProperty("RESOLUTIONSCALE", 2.0);
  Scint_LAB_Gd_mt->AddConstProperty("FASTTIMECONSTANT", 1.6*ns);
//  Scint_LAB_Gd_mt->AddConstProperty("SLOWTIMECONSTANT", 10.*ns);
  Scint_LAB_Gd_mt->AddConstProperty("YIELDRATIO", 1.0);

  fScint_LAB_Gd->SetMaterialPropertiesTable(Scint_LAB_Gd_mt); 
  
  LABFractionVolume = 1.0;  // 1.0 volume 
  PPOConcetration = 3E-3*g/cm3;   // 3 g/L
  bisMSBConcetration = 2.0E-5*g/cm3;  // 20 mg/L

  density = LABFractionVolume*fLAB->GetDensity() + PPOConcetration + bisMSBConcetration;
    
  G4Material* fScint_LAB = new G4Material("Scint_LAB", density,ncomponents=3, kStateLiquid, temperature);
  G4MaterialPropertiesTable *Scint_LAB_mt = new G4MaterialPropertiesTable();

  fScint_LAB->AddMaterial(fLAB, fractionmass = LABFractionVolume*fLAB->GetDensity()/density);
  mc = (density*fractionmass*1.0E+6/LAB_MolW)/g/mole;
  Scint_LAB_mt->AddConstProperty("LAB_MC",mc);
  
  fScint_LAB->AddMaterial(fPPO, fractionmass=PPOConcetration/density);
  mc = (density*fractionmass*1.0E+6/PPO_MolW)/g/mole;
  Scint_LAB_mt->AddConstProperty("PPO_MC",mc);
  
  fScint_LAB->AddMaterial(fbisMSB, fractionmass=bisMSBConcetration/density);
  mc = (density*fractionmass*1.0E+6/bisMSB_MolW)/g/mole;
  Scint_LAB_mt->AddConstProperty("bisMSB_MC",mc);


  G4double ScintAbsLength[LABALnum] = {0.696121, 0.732634, 0.743963, 0.736372, 0.709814, 0.674434, 0.643368, 0.614166, 0.58591, 0.561591, 0.552279, 0.565257, 0.600536, 0.654059, 
                                       0.721645, 0.788683, 0.840751, 0.871201, 0.880055, 0.865974, 0.826311, 0.759627, 0.684024, 0.621216, 0.575782, 0.533895, 0.485414, 0.438695, 
                                       0.409256, 0.408101, 0.438708, 0.498284, 0.582788, 0.684069, 0.795993, 0.896537, 0.975308, 1.02134, 1.03604, 1.00954, 0.939546, 0.84772, 
                                       0.764121, 0.70923, 0.676176, 0.633831, 0.568423, 0.492563, 0.439052, 0.425907, 0.458301, 0.542291, 0.681658, 0.873278, 1.11249, 1.3902, 
                                       1.70064, 2.05485, 2.40848, 2.70931, 2.92885, 3.11846, 3.38621, 3.76882, 4.33552, 5.07409, 5.9168, 6.76343, 7.73501, 8.82304, 10.164, 11.7071, 
                                       13.6467, 15.9708, 18.913, 22.1758, 25.9329, 30.2844, 36.3302, 41.5593, 46.1174, 52.8094, 57.815, 61.1628, 69.5068, 73.3944, 77.882, 81.3956, 
                                       83.1142, 88.3722, 91.7422, 95.1696, 97.5273, 96.6568, 102.369, 110.464, 109.625, 108.527, 113.953, 109.626, 117.669, 121.977, 115.782, 123.72, 
                                       123.72, 128.875, 134.897, 137.467, 134.478, 148.803, 147.788, 147.788, 149.834, 153.553, 163.403, 159.786, 158.037, 159.786, 171.832, 175.311, 
                                       177.467, 176.743, 179.676, 181.18, 181.941, 203.294, 199.547, 196.826, 171.155, 182.71, 181.941, 191.602, 195.937, 188.27, 226.709, 222.059, 
                                       227.904, 243.266, 243.267, 267.292, 259.291, 307.099, 296.583, 302.805, 300.703, 311.52, 335.667, 357.86, 354.927, 352.041, 323.144, 352.042, 
                                       366.958, 323.143, 366.959, 366.958, 412.389, 333.088, 313.78, 323.147, 366.959, 366.961, 370.095, 404.683, 397.261, 383.199, 451.053, 548.108, 
                                       455.8, 492.055, 404.685, 503.498, 509.422, 562.347, 541.257, 618.581, 666.157, 577.343, 746.556, 698.395, 746.557, 759.656, 759.656, 698.396, 
                                       609.869, 698.399, 656.073, 521.703, 497.72, 455.808, 397.269, 386.634, 313.796, 313.795, 373.298, 428.738, 475.843, 801.862, 1139.48, 1170.27, 
                                       1603.69, 1665.37, 2886.62, 2886.63, 3092.81, 1202.78, 1443.33, 1082.52, 902.097, 984.104, 902.099, 801.874, 601.42, 666.177, 548.13, 446.427, 
                                       424.538, 379.86, 320.777, 357.892, 333.129, 328.063, 346.433, 338.315, 376.557, 338.313, 307.135, 354.956, 296.614, 230.365, 225.566, 218.731, 
                                       265.684, 328.063, 379.869, 446.438, 481.15, 562.374, 1030.98, 1237.16, 21649.4, 2061.89, 21649.4, 43298.7, 43298.7, 43298.7, 43298.7, 43298.7, 
                                       43298.7, 8659.78, 8659.79, 2164.99, 2061.89, 849.057, 773.25, 441.87, 609.888, 420.426, 333.126, 290.65, 226.746};

  for(int i=0; i < LABALnum; ++i) {
    ScintAbsLength[i]=ScintAbsLength[i]*m;
  }

  Scint_LAB_mt->AddProperty("ABSLENGTH", LABAL_Energy, ScintAbsLength, SIZE(LABAL_Energy));
  Scint_LAB_mt->AddProperty("RINDEX", LABR_Energy, LABRindex, SIZE(LABR_Energy));
  Scint_LAB_mt->AddProperty("FASTCOMPONENT", ScintGdES_Energy, ScintGdEmissionSpectra, LABALnum);
//  Scint_LAB_Gd_mt->AddProperty("SLOWCOMPONENT", Scnt_PP, Scnt_SLOW, LABALnum);
  Scint_LAB_mt->AddConstProperty("SCINTILLATIONYIELD", 11500./MeV);
  Scint_LAB_mt->AddConstProperty("RESOLUTIONSCALE", 2.0);
  Scint_LAB_mt->AddConstProperty("FASTTIMECONSTANT", 1.6*ns);
//  Scint_LAB_Gd_mt->AddConstProperty("SLOWTIMECONSTANT", 10.*ns);
  Scint_LAB_mt->AddConstProperty("YIELDRATIO", 1.0);
  
  fScint_LAB->SetMaterialPropertiesTable(Scint_LAB_mt);


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the Birks Constant for the scintillator
  fScint_LAB_Gd->GetIonisation()->SetBirksConstant(0.11*mm/MeV);
  fScint_LAB->GetIonisation()->SetBirksConstant(0.11*mm/MeV);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //G4Material* fBufferMaterial = fLAB;
  //G4Material* fCatcherMaterial = fScint_LAB;
  //G4Material* fTargetMaterial = fScint_LAB_Gd;
  
  // Envelope parameters
  //
  //G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //

  G4double thick = 3*mm; // толщина стенок

  G4double world_sizeX = 5.*m;
  G4double world_sizeY = 5.*m;
  G4double world_sizeZ = 5.*m;
  G4Material* world_mat = fVacuum;
  //G4Material* world_mat = fAir;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Внешний бак детектора
  // 
  G4double pRMin_ex = 0*mm; //это диаметры
  G4double pRMax_ex = 1871*mm;
  G4double pDz_ex = 1301*mm;
  G4double pSPhi_ex = 0*deg;
  G4double pDPhi_ex = 360*deg;
  G4Material* external_mat = fStainlessSteel;

  G4Tubs* solidExternal =    
    new G4Tubs("External",                    //its name
        0.5*pRMin_ex, 0.5*(pRMax_ex + 2*thick), 0.5*(pDz_ex + 2*thick), pSPhi_ex, pDPhi_ex); //its size
      
  G4LogicalVolume* logicExternal =                         
    new G4LogicalVolume(solidExternal,       //its solid
                        external_mat,        //its material
                        "External");         //its name
               
  G4PVPlacement *physExternal = 
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicExternal,           //its logical volume
                    "External",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //     
  // Периферийный объем (LAB+PPO)
  //        
  G4double pRMin_pe = 0*mm; //это диаметры
  G4double pRMax_pe = 1871*mm;
  G4double pDz_pe = 1301*mm; 
  G4double pSPhi_pe = 0*deg;
  G4double pDPhi_pe = 360*deg;
  G4Material* peripheral_mat = fScint_LAB;

  G4Tubs* solidPeripheral =    
    new G4Tubs("Peripheral", 
        0.5*pRMin_pe, 0.5*pRMax_pe , 0.5*pDz_pe, pSPhi_pe, pDPhi_pe);
                      
  G4LogicalVolume* logicPeripheral =                         
    new G4LogicalVolume(solidPeripheral,     //its solid
                        peripheral_mat,      //its material
                        "Peripheral");       //its name
               
  G4PVPlacement *physPeripheral = 
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at position
                    logicPeripheral,         //its logical volume
                    "Peripheral",            //its name
                    logicExternal,           //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //     
  // Внутренний бак детектора
  // 
  G4double pRMin_in = 0*mm; //это диаметры
  G4double pRMax_in = 1250*mm;
  G4double pDz_in = 1301*mm;
  G4double pSPhi_in = 0*deg;
  G4double pDPhi_in = 360*deg;

  G4Tubs* solidInternal =    
    new G4Tubs("Internal",                   //its name
        0.5*pRMin_in, 0.5*(pRMax_in + 2*thick), 0.5*pDz_in, pSPhi_in, pDPhi_in); //its size
      
  G4LogicalVolume* logicInternal =                         
    new G4LogicalVolume(solidInternal,       //its solid
                        external_mat,        //its material
                        "Internal");         //its name
               
  G4PVPlacement *physInternal = 
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicInternal,           //its logical volume
                    "Internal",              //its name
                    logicPeripheral,         //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //     
  // Мишень (LAB+PPO+Gd)
  //         
  G4double pRMin_ta = 0*mm; //это диаметры
  G4double pRMax_ta = 1250*mm;
  G4double pDz_ta = 833*mm; 
  G4double pSPhi_ta = 0*deg;
  G4double pDPhi_ta = 360*deg;
  G4Material* target_mat = fScint_LAB_Gd;
  G4ThreeVector target_pos = G4ThreeVector(0, 0, -232*mm);

  G4Tubs* solidTarget =    
    new G4Tubs("Target", 
        0.5*pRMin_ta, 0.5*pRMax_ta, 0.5*pDz_ta, pSPhi_ta, pDPhi_ta);
                      
  G4LogicalVolume* logicTarget =                         
    new G4LogicalVolume(solidTarget,         //its solid
                        target_mat,          //its material
                        "Target");           //its name
               
  G4PVPlacement *physTarget =
  new G4PVPlacement(0,                       //no rotation
                    target_pos,              //at position
                    logicTarget,             //its logical volume
                    "Target",                //its name
                    logicInternal,           //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //     
  // Калибровочный канал
  //


  //     
  // Акриловая мембрана (Акрил)
  //
  //G4Material* membrane_mat = fAcrylic;
  G4Material* membrane_mat = fGlass;
  G4ThreeVector membrane_pos = G4ThreeVector(0, 0, 187.5*mm);

  // Trapezoid shape       
  G4double pRMin_ac = 0*mm; //это диаметры
  G4double pRMax_ac = 1250*mm;
  G4double pDz_ac = 10*mm;
  G4double pSPhi_ac = 0*deg;
  G4double pDPhi_ac = 360*deg;

  G4Tubs* solidMembrane =    
    new G4Tubs ("Membrane", 
        0.5*pRMin_ac, 0.5*pRMax_ac, 0.5*pDz_ac, pSPhi_ac, pDPhi_ac);
                
  G4LogicalVolume* logicMembrane =                         
    new G4LogicalVolume(solidMembrane,       //its solid
                        membrane_mat,        //its material
                        "Membrane");         //its name
               
  G4PVPlacement *physMembrane =
  new G4PVPlacement(0,                       //no rotation
                    membrane_pos,            //at position
                    logicMembrane,           //its logical volume
                    "Membrane",              //its name
                    logicInternal,           //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //     
  // Гамма-кетчер (LAB)
  //
  G4Material* catcher_mat = fLAB;
  G4ThreeVector catcher_pos = G4ThreeVector(0, 0, 419.5*mm);

  // Trapezoid shape       
  G4double pRMin_ca = 0*mm; //это диаметры
  G4double pRMax_ca = 1250*mm;
  G4double pDz_ca = 458*mm;
  G4double pSPhi_ca = 0*deg;
  G4double pDPhi_ca = 360*deg;

  G4Tubs* solidCatcher =    
    new G4Tubs ("Catcher", 
        0.5*pRMin_ca, 0.5*pRMax_ca, 0.5*pDz_ca, pSPhi_ca, pDPhi_ca);
                
  G4LogicalVolume* logicCatcher =                         
    new G4LogicalVolume(solidCatcher,        //its solid
                        catcher_mat,         //its material
                        "Catcher");          //its name
               
  G4PVPlacement *physCatcher =
  new G4PVPlacement(0,                       //no rotation
                    catcher_pos,             //at position
                    logicCatcher,            //its logical volume
                    "Catcher",               //its name
                    logicInternal,           //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //
  //ФЭУ
  //
  G4double Diam_min_in = 1250*mm; // диаметр(внутренний) внутреннего объема
  G4double Diam_min_ex = 1871*mm; // диаметр(внутренний) внешнего объема
  //G4double z_A = 833*mm; // высота объема А
  //G4double z_B = 1301*mm; // высота объема В (высота детектора)
  //G4double thick = 2*mm; // толщина стенок
  G4double Diam_pmt = 200*mm; // диаметр ФЭУ
  G4double thick_pmt = 1*mm; // толщина ФЭУ
  //G4double R1_pmt = (Diam_min_in/8) + 14*mm; // Радиус 1-го круга ФЭУ
  //G4double R2_pmt = (Diam_min_in/4) + (Diam_min_in/8) + 14*mm; // Радиус 2-го круга ФЭУ
  //G4double R3_pmt = (Diam_min_ex/2 - Diam_min_in/2)/2 + Diam_min_in + thick; // Радиус 3-го круга ФЭУ

  G4double R1_pmt = Diam_min_in/8;
  G4double R2_pmt = Diam_min_in/4 + R1_pmt;
  G4double R3_pmt = (Diam_min_ex/2 - Diam_min_in/2)/2 + Diam_min_in/2;

  // Весь объем реплики внутри для 4-х 
  G4ThreeVector pmt_vol_pos = G4ThreeVector(0, 0, 227.5*mm);

  G4Tubs* pmt_vol1_tube = 
    new G4Tubs("pmt_vol1",
        0, Diam_min_in/4, thick_pmt/2, 0, 360*deg);
  G4LogicalVolume* pmt_vol1_log =
    new G4LogicalVolume(pmt_vol1_tube,
                        fN2, 
                        "pmt_vol1");
  G4PVPlacement *pmt_vol1_phys =
  new G4PVPlacement(0,
                    pmt_vol_pos,
                    pmt_vol1_log,
                    "pmt_vol1",
                    logicCatcher,
                    false,
                    0,
                    checkOverlaps);
  /////////////////////////////////посчитать радиусы
    
  // Часть, которая копируется 
  G4Tubs* pmt_tube =
    new G4Tubs("pmt",
        0, Diam_pmt/2, thick_pmt/2, 0, 360*deg);
  pmt_log =
    new G4LogicalVolume(pmt_tube,
                        fGlass,
                        "pmt");
  //G4double x1 = (i-7)*10.*mm;
  //G4double xy = (i-7)*10.*mm;
  //G4double z1 = (i-7)*10.*mm;
    for (G4int i=1;i<5;i++)
    {
        new G4PVPlacement(0,G4ThreeVector(R1_pmt*cos(M_PI*0.5*i*rad),R1_pmt*sin(M_PI*0.5*i*rad),0),
                          pmt_log,
                          "PMT",pmt_vol1_log,
                          false,i,checkOverlaps);
    }
  
  /*//1 круг
  new G4PVReplica("pmt",
                  pmt_log,
                  pmt_vol1_log,
                  kPhi,
                  4,
                  M_PI*0.5*rad); 
   */
  // внутри для 12-и 
  G4Tubs* pmt_vol2_tube =
    new G4Tubs("pmt_vol2",
        Diam_min_in/4, Diam_min_in/2, thick_pmt/2, 0, 360*deg);
  G4LogicalVolume* pmt_vol2_log =
    new G4LogicalVolume(pmt_vol2_tube,
                        fVacuum,
                        "pmt_vol2");
  G4PVPlacement *pmt_vol2_phys =
  new G4PVPlacement(0,
                    pmt_vol_pos,
                    pmt_vol2_log,
                    "pmt_vol2",
                    logicCatcher,
                    false,
                    0,
                    checkOverlaps);

  for (G4int i=1;i<13;i++)
    {
        new G4PVPlacement(0,G4ThreeVector(R2_pmt*cos(M_PI/6*i*rad),R2_pmt*sin(M_PI/6*i*rad),0),
                          pmt_log,
                          "PMT",pmt_vol2_log,
                          false,i,checkOverlaps);
    }

  /*//2 круг
  new G4PVReplica("pmt",
                  pmt_log,
                  pmt_vol2_log,
                  kPhi,
                  12,
                  M_PI/6*rad,
                  0);
   */
  //снаружи 
  G4ThreeVector pmt_vol3_pos = G4ThreeVector(0, 0, 649*mm);

  G4Tubs* pmt_vol3_tube =
    new G4Tubs("pmt_vol3", Diam_min_in/2, Diam_min_ex/2, thick_pmt/2, 0, 360*deg);
  G4LogicalVolume* pmt_vol3_log =
    new G4LogicalVolume(pmt_vol3_tube,
                        fVacuum,
                        "pmt_vol3");
  G4PVPlacement *pmt_vol3_phys =
  new G4PVPlacement(0,
                    pmt_vol3_pos,
                    pmt_vol3_log,
                    "pmt_vol3",
                    logicPeripheral,
                    false,
                    0,
                    checkOverlaps);

  for (G4int i=1;i<13;i++)
    {
        new G4PVPlacement(0,G4ThreeVector(R3_pmt*cos(M_PI/6*i*rad),R3_pmt*sin(M_PI/6*i*rad),0),
                          pmt_log,
                          "PMT",pmt_vol3_log,
                          false,i,checkOverlaps);
    }
    
/*
  //3 круг
  new G4PVReplica("pmt",
                  pmt_log,
                  pmt_vol3_log,
                  kPhi,
                  12,
                  M_PI/6*rad,
                  0);
  */
  //Surface properties
 /*   G4OpticalSurface* scintTank = new G4OpticalSurface("ScintTank");
 
    new G4LogicalBorderSurface("ScintTank", physExternal,
                               physPeripheral,
                               scintTank);
    new G4LogicalBorderSurface("ScintTank", physInternal,
                               physPeripheral,
                               scintTank);
    new G4LogicalBorderSurface("ScintTank", physInternal,
                               physTarget,
                               scintTank);
    new G4LogicalBorderSurface("ScintTank", physInternal,
                               physMembrane,
                               scintTank);
    new G4LogicalBorderSurface("ScintTank", physInternal,
                               physCatcher,
                               scintTank);
    new G4LogicalBorderSurface("ScintTank", physCatcher,
                               pmt_vol1_phys,
                               scintTank);
    new G4LogicalBorderSurface("ScintTank", physCatcher,
                               pmt_vol2_phys,
                               scintTank);
    new G4LogicalBorderSurface("ScintTank", physPeripheral,
                               pmt_vol3_phys,
                               scintTank);
 
    scintTank->SetType(dielectric_metal);
    scintTank->SetFinish(polished);
    scintTank->SetModel(glisur);

    const G4int numSurf = 2;

    G4double pp[numSurf] = {2.0*eV, 3.5*eV};
    G4double reflectivity[numSurf] = {0.4, 0.4};
    G4double efficiency[numSurf] = {0.0, 0.0};
    
    G4MaterialPropertiesTable* scintTankProperty 
      = new G4MaterialPropertiesTable();

    scintTankProperty->AddProperty("REFLECTIVITY",pp,reflectivity,numSurf);
    scintTankProperty->AddProperty("EFFICIENCY",pp,efficiency,numSurf);
    scintTank->SetMaterialPropertiesTable(scintTankProperty);
    */
    //**Sphere surface properties
    const G4int numSurf = 2;

    G4double pp[numSurf] = {2.0*eV, 3.5*eV};
    G4double reflectivity[numSurf] = {0.9, 0.9};
    G4double efficiency[numSurf] = {0.0, 0.0};

    G4MaterialPropertiesTable* scintTankProperty = new G4MaterialPropertiesTable();
    scintTankProperty->AddProperty("REFLECTIVITY", pp, reflectivity, numSurf);
    scintTankProperty->AddProperty("EFFICIENCY", pp, efficiency, numSurf);
    G4OpticalSurface* OpTankSurface =
      new G4OpticalSurface("SphereSurface",unified,polished,dielectric_metal);
    OpTankSurface->SetMaterialPropertiesTable(scintTankProperty);
    
    new G4LogicalSkinSurface("tank_surface_ex",logicExternal,OpTankSurface);
    new G4LogicalSkinSurface("tank_surface_in",logicInternal,OpTankSurface);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // visualization attributes ------------------------------------------------
  /*logicWorld
  logicExternal
  logicPeripheral
  logicInternal
  logicTarget
  logicMembrane
  logicCatcher
  pmt_log*/
    
    G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    visAttributes->SetVisibility(false);
    logicWorld->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.5,0.5,0.5));   // LightGray
    logicExternal->SetVisAttributes(visAttributes);
    logicInternal->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    //visAttributes->SetVisibility(false);
    //logicPeripheral->SetVisAttributes(visAttributes);
    logicMembrane->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.8888,0.0,0.0));
    logicTarget->SetVisAttributes(visAttributes);
    //fHodoscope2Logical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    logicPeripheral->SetVisAttributes(visAttributes);
    //chamber2Logical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    /*visAttributes = new G4VisAttributes(G4Colour(0.0,0.8888,0.0));
    visAttributes->SetVisibility(false);
    fWirePlane1Logical->SetVisAttributes(visAttributes);
    fWirePlane2Logical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    */
    visAttributes = new G4VisAttributes(G4Colour(0.8888,0.0,0.8888));
    //visAttributes->SetVisibility(false);
    logicCatcher->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.0));
    pmt_log->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    /*visAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));
    hadCalorimeterLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    visAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));
    visAttributes->SetVisibility(false);
    HadCalColumnLogical->SetVisAttributes(visAttributes);
    HadCalCellLogical->SetVisAttributes(visAttributes);
    HadCalLayerLogical->SetVisAttributes(visAttributes);
    fHadCalScintiLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    */              
  // Set Shape2 as scoring volume
  //
  //fScoringVolume = logicShape2;

    // sensitive detectors -----------------------------------------------------
  /*  G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    iDREAMPMTSD* PMTs = new iDREAMPMTSD("PMT");
    SDman->AddNewDetector(PMTs);
    pmt_log->SetSensitiveDetector(PMTs);
*/
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void iDREAMDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "iDREAMPMTSD";
  iDREAMPMTSD* PMTSDSD = new iDREAMPMTSD(trackerChamberSDname);
  SetSensitiveDetector( pmt_log,  PMTSDSD );
}