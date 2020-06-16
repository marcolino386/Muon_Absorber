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
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "B1SD.hh"
#include "G4VSensitiveDetector.hh"
#include "G4EquationOfMotion.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4UniformMagField.hh"
#include "G4UniformElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4DormandPrince745.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Cons.hh"
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume1(0),
  fScoringVolume2(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
 

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;



  //
  // World
  //
 
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  G4Box* solidWorld =
    new G4Box("World",                       //its name
      5*m, 5*m, 10*m);     //its size

    logicWorld =
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



//some importante initial values

G4double initial_angle = 2.;
  G4double final_angle = 10.;
 
 # define PI 3.14159265


G4double pos_after_detec = 3.*m;
G4double z_0 = 90.0*cm;
G4double carbon_pDz = 225.0*cm/2;
G4double concrete_pDz = 126.*cm/2;
G4double dzFaFlange = 2.*cm;
 G4double dzCarbonConeS = 108.3*cm/2;
//G4double mag_position = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + 2*tungs3_dz + 2*tungs4_pDz + 2*dzFaWTail1 + 2*dzFaWTail2 + 2*dzFaWTailR + 2*dzFaWTailB + tail_z + 5*m)/2;
G4double mag_position = 5*m;
//_________


  G4double dzFaWPlateF = 2.00*cm/2;
  G4double rInFaQPlateF = 20.50*cm;
  G4double rOuFaQPlateF = 40.05*cm;
  // 1st Central Part 24 deg
  G4double dzFaWPlateC1 = 7.95*cm/2;
  G4double rInFaQPlateC1 = 16.35*cm;
  G4double rOuFaQPlateC1 = rOuFaQPlateF + (dzFaWPlateF * tan(24.*PI/180.00));
  // 2nd Central Part 5 deg
  G4double dzFaWPlateC2 = 1.05*cm/2;
  G4double rInFaQPlateC2 = rInFaQPlateC1 + (dzFaWPlateC1 * tan(initial_angle*PI/180.00));
  G4double rOuFaQPlateC2 = rOuFaQPlateC1 + (dzFaWPlateC1 * tan(24.*PI/180.00));
  G4double rInFaQPlateC3 = 17.94*cm;
  G4double rOuFaQPlateC3 = 44.49*cm;
  // Rear Flange
  G4double dzFaWPlateR = 1.00*cm/2;
  G4double rInFaQPlateR = 21.00*cm;
  G4double rOuFaQPlateR = 42.55*cm;
  // Lenth of Plate - Rear Flange
  G4double dzFaWPlate = dzFaWPlateF + dzFaWPlateC1 + dzFaWPlateC2;

//-------------------------------






//Area with magnetic field

G4Box* solidMag =
    new G4Box("World",                       //its name
      3*m, 3*m, mag_position);

    logicMag =
    new G4LogicalVolume(solidMag,          //its solid
                        world_mat,           //its material
                        "MagField_box"); 


       new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0, mag_position),       //at (0,0,0)
                      logicMag,            //its logical volume
                      "MagField",               //its name
                      logicWorld,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking




 
 //-------- ABSORBER --------------
 //materials

  //Materiais from CERN model

  std::vector<G4String> symbol = {"H", "He", "Li","Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba","La","Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu" ,"Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Ti", "Pb", "Bi"};




 //Carbon Material 

 G4Material* kMedCSh = new G4Material("ABSO_C_C2", 6.,  12.01*g/mole,  1.75*g/cm3);
//matmgr.Material("ABSO", 46, "CARBON2$", 12.01, 6., 1.75, 24.4, 49.9);
 G4cout << kMedCSh << G4endl;

//LEAD - Pb

G4Material* kMedPb = new G4Material("ABSO_PB_C0", 82., 207.19*g/mole,15.9994*g/cm3);
G4cout << kMedPb<< G4endl;
//matmgr.Material("ABSO", 53, "LEAD2$", 207.19, 82., 11.35, .56, 18.5);
 //Air 

  std::vector<G4double> aAir = {12.0107, 14.0067, 15.9994, 39.948};
  std::vector<G4double> zAir = {6., 7., 8., 18.};
  std::vector<G4double> wAir = {0.000124, 0.755267, 0.231781, 0.012827};
  

// Polyethylene
  //
  std::vector<G4double> apoly = {12.01, 1.};
  std::vector<G4double> zpoly = {6., 1.};
  std::vector<G4double> wpoly = {.33, .67};

// Concrete
  //
  std::vector<G4double> aconc = {1.*g/mole, 12.01*g/mole, 15.994*g/mole, 22.99*g/mole, 24.305*g/mole, 26.98*g/mole, 28.086*g/mole, 39.1*g/mole, 40.08*g/mole, 55.85*g/mole};
  std::vector<G4double> zconc = {1., 6., 8., 11., 12., 13., 14., 19., 20., 26.};
  std::vector<G4double> wconc = {.01, .001, .529107, .016, .002, .033872, .337021, .013, .044, .014};
//matmgr.Mixture("ABSO", 57, "CONCRETE2$", aconc, zconc, 2.35, 10, wconc);

std::vector<G4Element *>conc;
conc.reserve(aconc.size());

G4Material* kMedConcSh = new G4Material("CONCRETE CC2", 2.35*g/cm3, aconc.size());


for (std::size_t i = 0; i < aconc.size(); i++) {
   G4String symb = symbol[int(zconc[i]) - 1];
   conc[i] = new G4Element(("conc" + std::to_string(i)),symb,zconc[i], aconc[i]);
   kMedConcSh->AddElement(conc[i], wconc[i]);
   //G4cout << steel[i] << G4endl;

}

G4cout << kMedConcSh << G4endl;

// Steel
  //
  std::vector<G4double> asteel = {55.847*g/mole, 51.9961*g/mole, 58.6934*g/mole, 28.0855*g/mole};
  std::vector<G4double> zsteel = {26., 24., 28., 14.};
  std::vector<G4double> wsteel = {.715, .18, .1, .005};
//
  std::vector<G4Element *>steel;
  steel.reserve(asteel.size());

  G4Material* kMedSteel = new G4Material("STAINLESS STEEL0", 7.88*g/cm3, asteel.size());

  for (std::size_t i = 0; i < asteel.size(); i++) {
   G4String symb = symbol[int(zsteel[i]) - 1];
   steel[i] = new G4Element(("Steel" + std::to_string(i)),symb,zsteel[i], asteel[i]);
   kMedSteel->AddElement(steel[i], wsteel[i]);
   //G4cout << steel[i] << G4endl;

}
 
 G4cout << kMedSteel << G4endl;
 
 // Ni-Cu-W alloy
  //
  std::vector<G4double> aniwcu = {58.6934*g/mole, 183.84*g/mole, 63.546*g/mole};
  std::vector<G4double> zniwcu = {28., 74., 29.};
  std::vector<G4double> wniwcu = {0.015, 0.95, 0.035};
  
  std::vector<G4Element *>Ni_Cu_W;
  Ni_Cu_W.reserve(aniwcu.size());
  
  G4Material* kMedNiW = new G4Material("ABSO_Ni/W0", 18.78*g/cm3, aniwcu.size());
  
  for (std::size_t i = 0; i < aniwcu.size(); i++) {
   G4String symb = symbol[int(zniwcu[i]) - 1];
   G4cout << symb << G4endl;
   Ni_Cu_W[i] = new G4Element(("ABSO_Ni/W0" + std::to_string(i)),symb,zniwcu[i], aniwcu[i]);
   kMedNiW->AddElement(Ni_Cu_W[i], wniwcu[i]);
   //G4cout << steel[i] << G4endl;

}


  //
  // Poly Concrete
  //                      H     Li     F       C      Al     Si      Ca      Pb     O
  std::vector<G4double> aPolyCc = {1.*g/mole, 6.941*g/mole, 18.998*g/mole, 12.01*g/mole, 26.98*g/mole, 28.086*g/mole, 40.078*g/mole, 207.2*g/mole, 15.999*g/mole};
  std::vector<G4double> zPolyCc = {1., 3., 9., 6., 13., 14., 20., 82., 8.};
  std::vector<G4double> wPolyCc = {4.9, 1.2, 1.3, 1.1, 0.15, 0.02, 0.06, 0.7, 1.1};
  G4double wtot = 0;

   for (G4int i = 0; i < 9; i++) {
    wtot += wPolyCc[i];
  }
   for (G4int i = 0; i < 9; i++) {
    wPolyCc[i] /= wtot;
  }

 // matmgr.Mixture("ABSO", 58, "POLYETHYLEN2$", apoly, zpoly, .95, 2, wpoly);
  
  std::vector<G4Element* >CH2;
  CH2.reserve(aPolyCc.size());

  G4Material* KmedCH2Sh = new G4Material("POLYETHYLEN2$", .95*g/mole, aPolyCc.size() );
  
  for (std::size_t i = 0; i < aPolyCc.size(); i++) {
   G4String symb = symbol[int(zPolyCc[i]) - 1];
   CH2[i] = new G4Element(("Poly" + std::to_string(i)),symb,zPolyCc[i], aPolyCc[i]);
   KmedCH2Sh->AddElement(CH2[i], wPolyCc[i]);
   //G4cout << steel[i] << G4endl;

}
   
/*
 d = 8.02*g/cm3;


  //initial values from the aborber paper

 

  //cone trunks



  //carbon cone
 
*/
  //front 
  
  // Inner radius at the front
  G4double carbon1_pRmin1 = 4.5*cm;
  // Outer radius at the front
  G4double carbon1_pRmax1 = (z_0 + dzFaFlange)*tan(final_angle*PI/180.00);
  // Inner radius at start of inner opening cone
  G4double carbon1_pRmin2 = 7.0*cm;
  // Outer radius at start of inner opening cone
  G4double carbon1_pRmax2 = (z_0 + dzFaFlange + 2*dzCarbonConeS)*tan(final_angle*PI/180.00);
  G4double carbon_pSphi = 0.*deg;
  G4double carbon_pDphi = 360.*deg;
  
  G4Cons* carbon1_cons = new G4Cons("carbon1_cons", carbon1_pRmin1, carbon1_pRmax1, 
		                   carbon1_pRmin1, carbon1_pRmax2, dzCarbonConeS,
				    carbon_pSphi, carbon_pDphi);

  G4LogicalVolume* carbon1_Lvolume = new G4LogicalVolume(carbon1_cons, kMedCSh, "carbon1_logical");


 G4double carbon1_z = (z_0  + dzCarbonConeS ) -  mag_position ;


// aaa


 G4VisAttributes* blue = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
 carbon1_Lvolume->SetVisAttributes(blue);

 Logical_volumes.push_back(carbon1_Lvolume);

// inner opening cone

G4double carbon_Dz_2 = (2*carbon_pDz - 2*dzCarbonConeS)/2;
G4double carbon_pRmin1 = carbon1_pRmin2;
G4double carbon_pRmax1 = carbon1_pRmax2;
G4double carbon_pRmin2 = 11.0*cm;
G4double carbon_pRmax2 = (z_0 + 2*carbon_pDz)*tan(final_angle*PI/180.00);

G4double carbon_z = (z_0  + 2*dzCarbonConeS + carbon_Dz_2 ) - mag_position;

G4Cons* carbon_cons = new G4Cons("carbon_cons", carbon_pRmin1, carbon_pRmax1, 
		                   carbon_pRmin2, carbon_pRmax2, carbon_Dz_2,
				    carbon_pSphi, carbon_pDphi);



G4LogicalVolume* carbon_Lvolume = new G4LogicalVolume(carbon_cons, kMedCSh, "carbon_logical");

carbon_Lvolume->SetVisAttributes(blue);
Logical_volumes.push_back(carbon_Lvolume);

 //concrete cone
 ////////////////AAAARRRUMARRRRR/////////////////
 

  G4double concrete_pRmin1 = 11.*cm;  // inner radius at the begining of cone /
  G4double concrete_pRmax1 = carbon_pRmax2; // outer radius at the begining of cone



  G4double concrete_pRmin2 = concrete_pRmin1  + (2*concrete_pDz)*tan(initial_angle*PI/180.00); // inner radius at the end of cone
  G4double concrete_pRmax2 = concrete_pRmax1 + (2*concrete_pDz)*tan(final_angle*PI/180.00); // outer radius at the end of cone
   G4double concrete_pSphi = 0.*deg; //starting angle
  G4double concrete_pDphi = 360.*deg; // eng angle
  
 G4Cons* concrete_cons = new G4Cons("concrete_cons", concrete_pRmin1, concrete_pRmax1, concrete_pRmin2, concrete_pRmax2, concrete_pDz, concrete_pSphi, concrete_pDphi);

 G4LogicalVolume* concrete_Lvolume = new G4LogicalVolume(concrete_cons, kMedConcSh, "concrete_logical");

G4double concrete_z = (z_0 + 2*carbon_pDz + concrete_pDz ) - mag_position;


 


G4VisAttributes* red = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
concrete_Lvolume->SetVisAttributes(red);

 Logical_volumes.push_back(concrete_Lvolume);

 G4Colour  black   (0.0, 0.0, 0.0) ;  // black



 // Pos 4+5
  ///////////////////////////////////
  //    FA W Plate A+B             //
  //    Drawing ALIP2A__0043       //
  ///////////////////////////////////
  // Front Flange
  //TGeoPcon* shFaWPlateA = new TGeoPcon(0., 360., 7);
  //z = 0.;
  // Front Flange
  
  
  G4Cons* plateAB1_cons = new G4Cons("plateAB1_cons", rInFaQPlateF, rOuFaQPlateF, 
		                   rInFaQPlateF, rOuFaQPlateC1, dzFaWPlateF,
				    carbon_pSphi, carbon_pDphi);

  G4LogicalVolume* plateAB1_Lvolume = new G4LogicalVolume(plateAB1_cons, kMedNiW, "plateAB1_logical");

  Logical_volumes.push_back(plateAB1_Lvolume);

  G4double plateAB1_z = (z_0 + dzFaWPlateF) - mag_position;  

  // 24 deg cone

  G4Cons* plateAB2_cons = new G4Cons("plateAB2_cons", rInFaQPlateC1, rOuFaQPlateC1, 
		                   rInFaQPlateC2, rOuFaQPlateC2, dzFaWPlateC1,
				    carbon_pSphi, carbon_pDphi);

  G4LogicalVolume* plateAB2_Lvolume = new G4LogicalVolume(plateAB2_cons, kMedNiW, "plateAB2_logical");

  Logical_volumes.push_back(plateAB2_Lvolume);

   G4double plateAB2_z = (z_0 + 2*dzFaWPlateF + dzFaWPlateC1) - mag_position; 


  
  // 5 deg cone
  G4Cons* plateAB3_cons = new G4Cons("plateAB3_cons", rInFaQPlateC2, rOuFaQPlateC2, 
		                   rInFaQPlateC3, rOuFaQPlateC3, dzFaWPlateC2,
				    carbon_pSphi, carbon_pDphi);

  G4LogicalVolume* plateAB3_Lvolume = new G4LogicalVolume(plateAB3_cons, kMedNiW, "plateAB3_logical");
  
  Logical_volumes.push_back(plateAB3_Lvolume);

  G4double plateAB3_z = (z_0 + 2*dzFaWPlateF + 2*dzFaWPlateC1 + dzFaWPlateC2) - mag_position;
  
 
  // Rear Flange
  G4Cons* plateAB4_cons = new G4Cons("plateAB4_cons", rInFaQPlateR, rOuFaQPlateR, 
		                   rInFaQPlateR, rOuFaQPlateR, dzFaWPlateR,
				    carbon_pSphi, carbon_pDphi);

  G4LogicalVolume* plateAB4_Lvolume = new G4LogicalVolume(plateAB4_cons, kMedNiW, "plateAB4_logical");
  
  Logical_volumes.push_back(plateAB4_Lvolume);

  G4double plateAB4_z = (z_0 + 2*dzFaWPlateF + 2*dzFaWPlateC1 + 2*dzFaWPlateC2 + dzFaWPlateR) - mag_position;
  
  
   
  
   plateAB1_Lvolume->SetVisAttributes(gray);
  plateAB2_Lvolume->SetVisAttributes(gray);
  plateAB3_Lvolume->SetVisAttributes(gray);
  plateAB4_Lvolume->SetVisAttributes(gray);


// Pos 1
  ///////////////////////////////////
  //    FA Steel Envelope          //
  //    Drawing ALIP2A__0036       //
  ///////////////////////////////////
  // Thickness of the envelope
  G4double dSteelEnvelope = 1.5*cm/2;
  // Front cover
  //
  // Length
  G4double dzSteelEnvelopeFC = 4.00*cm/2;
  // Inner Radius
  G4double rInSteelEnvelopeFC1 = 35.90*cm / 2.;
  G4double rInSteelEnvelopeFC2 = rInSteelEnvelopeFC1 + (dzSteelEnvelopeFC * tan(final_angle*PI/180.00));
  // Outer Radius
  G4double rOuSteelEnvelopeFC1 = 88.97*cm / 2.;
  G4double rOuSteelEnvelopeFC2 = rOuSteelEnvelopeFC1 + (dzSteelEnvelopeFC * tan(5.*PI/180.00));
  //
  // 5 deg cone
  G4double dzSteelEnvelopeC5 = 168.9*cm/2;
  G4double rInSteelEnvelopeC5 = rOuSteelEnvelopeFC2 - (dSteelEnvelope / cos(5. * PI/180.00));
  G4double rOuSteelEnvelopeC5 = rOuSteelEnvelopeFC2;
  // 10 deg cone
  G4double dzSteelEnvelopeC10 = (227.1*cm - 4.*cm)/2;
  G4double rInSteelEnvelopeC10 = 116.22*cm / 2.;
  G4double rOuSteelEnvelopeC10 = rInSteelEnvelopeC10 + (dSteelEnvelope / cos(10 *PI/180.00));
  // Rear ring
  G4double dzSteelEnvelopeR = 4.*cm/2;
  G4double rInSteelEnvelopeR2 = 196.3*cm / 2.;
  G4double rOuSteelEnvelopeR2 = 212.0*cm / 2.;
  G4double rInSteelEnvelopeR1 = rInSteelEnvelopeR2 - (dzSteelEnvelopeR * tan(final_angle*PI/180.00));
  G4double rOuSteelEnvelopeR1 = rInSteelEnvelopeR1 + (dSteelEnvelope / cos(final_angle*PI/180.00));
  // Front insert

  G4double dzSteelEnvelopeFI = 1.*cm/2;
  G4double rInSteelEnvelopeFI = 42.0*cm / 2.;
  G4double rOuSteelEnvelopeFI = (85.0*cm / 2.) + 0.06*cm;

  //TGeoPcon* shFaSteelEnvelopeC = new TGeoPcon(0., 360., 7);
  
  // Front cover
  G4VisAttributes* color = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));

  G4Cons* sEnvelope1_cons = new G4Cons("Senvelope1_cons", rInSteelEnvelopeFC1, rOuSteelEnvelopeFC1, 
		                   rInSteelEnvelopeFC2, rOuSteelEnvelopeFC2, dzSteelEnvelopeFC,
				    carbon_pSphi, carbon_pDphi);

  G4LogicalVolume* sEnvelope1_Lvolume = new G4LogicalVolume(sEnvelope1_cons, kMedSteel, "sEnvelope1_logical");

  Logical_volumes.push_back(sEnvelope1_Lvolume);

  G4double sEnvelope1_z = (z_0 + dzSteelEnvelopeFC + 2*dzFaWPlate ) - mag_position;
  
  // 5 deg cone
  
  G4Cons* sEnvelope2_cons = new G4Cons("Senvelope2_cons", rInSteelEnvelopeC5, rOuSteelEnvelopeC5, 
		                   rInSteelEnvelopeC10, rOuSteelEnvelopeC10, dzSteelEnvelopeC5,
				    carbon_pSphi, carbon_pDphi);

  G4LogicalVolume* sEnvelope2_Lvolume = new G4LogicalVolume(sEnvelope2_cons, kMedSteel, "sEnvelope2_logical");

  Logical_volumes.push_back(sEnvelope2_Lvolume);

   G4double sEnvelope2_z = (z_0 + 2*dzSteelEnvelopeFC + dzSteelEnvelopeC5 + 2*dzFaWPlate ) - mag_position;
 
  // 10 deg cone
  
  G4Cons* sEnvelope3_cons = new G4Cons("Senvelope3_cons", rInSteelEnvelopeC10, rOuSteelEnvelopeC10, 
		                   rInSteelEnvelopeR1, rOuSteelEnvelopeR1, dzSteelEnvelopeC10,
				    carbon_pSphi, carbon_pDphi);

  G4LogicalVolume* sEnvelope3_Lvolume = new G4LogicalVolume(sEnvelope3_cons, kMedSteel, "sEnvelope3_logical");

  Logical_volumes.push_back(sEnvelope3_Lvolume);

  G4double sEnvelope3_z = (z_0 + 2*dzSteelEnvelopeFC + 2*dzSteelEnvelopeC5 + dzSteelEnvelopeC10+ 2*dzFaWPlate) - mag_position;

  // Rear Ring
  

  G4Cons* sEnvelope4_cons = new G4Cons("Senvelope4_cons", rInSteelEnvelopeR1, rOuSteelEnvelopeR2, 
		                   rInSteelEnvelopeR2, rOuSteelEnvelopeR2, dzSteelEnvelopeR,
				    carbon_pSphi, carbon_pDphi);

  G4LogicalVolume* sEnvelope4_Lvolume = new G4LogicalVolume(sEnvelope4_cons, kMedSteel, "sEnvelope4_logical");

  Logical_volumes.push_back(sEnvelope4_Lvolume);
 
  G4double sEnvelope4_z = (z_0 + 2*dzSteelEnvelopeFC + 2*dzSteelEnvelopeC5 + 2*dzSteelEnvelopeC10 + dzSteelEnvelopeR+ 2*dzFaWPlate) - mag_position;
    
  // Insert

  G4Tubs* sEnvelope_tub = new G4Tubs("sEnvelope_tubs", rInSteelEnvelopeFI, rOuSteelEnvelopeFI, dzSteelEnvelopeFI, 0.*deg,360.*deg);
  
 G4LogicalVolume* sEnvelopeT_Lvolume = new G4LogicalVolume(sEnvelope_tub, kMedSteel, "sEnvelopeT_logical");//FALTA MATERIAL
 
 Logical_volumes.push_back(sEnvelopeT_Lvolume);
 
 G4double sEnvelopeT_z = (z_0 + 2*dzSteelEnvelopeFC + 2*dzSteelEnvelopeC5 + 2*dzSteelEnvelopeC10 + 2*dzSteelEnvelopeR + 2*dzFaWPlate) - mag_position;




sEnvelope1_Lvolume->SetVisAttributes(color);
sEnvelope2_Lvolume->SetVisAttributes(color);
sEnvelope3_Lvolume->SetVisAttributes(color);
sEnvelope4_Lvolume->SetVisAttributes(color);
sEnvelopeT_Lvolume->SetVisAttributes(color);


  
//  TGeoVolume* voFaSteelEnvelope = new TGeoVolume("AFaSteelEnvelope", shFaSteelEnvelope, kMedSteel);



 // Pos 2
  ///////////////////////////////////
  //    FA End Plate               //
  //    Drawing ALIP2A__0037       //
  ///////////////////////////////////
  //
  //
  //
  //    Outer dimensions dx, dy, dz
  G4double dxEndPlate = 220.0*cm/2;
  G4double dyEndPlate = 220.0*cm/2;
  G4double dzEndPlate = 6.0*cm/2;
  //    Inner radius
  G4double rInEndPlate = 52.5*cm / 2.;
  //    Insert
  G4double rInEndPlateI = 175.3*cm / 2.;
  G4double rOuEndPlateI = 212.2*cm / 2.;
  G4double dzEndPlateI = 2.0*cm/2;


  G4Box* endplate_box =
    new G4Box("endplate_box",                       //its name
      dxEndPlate, dyEndPlate, dzEndPlate);

  G4LogicalVolume* endplate_Lvolume = new G4LogicalVolume(endplate_box, kMedSteel, "sEnvelope1_logical");//MISSING MATERIAL!!!!!!!!!

  Logical_volumes.push_back(endplate_Lvolume);

  G4double endplate_z = (z_0 + 2*dzSteelEnvelopeFC + 2*dzSteelEnvelopeC5 + 2*dzSteelEnvelopeC10 + 2*dzSteelEnvelopeR + 2*dzFaWPlate + dzEndPlate) - mag_position;

   endplate_Lvolume->SetVisAttributes(gray);
   
  G4Tubs* endplate_tub1 = new G4Tubs("endplate1_tub", 0., rInEndPlate, (dzEndPlate + 0.1*cm), 0.*deg,360.*deg);
  
  G4LogicalVolume* endplate2_Lvolume = new G4LogicalVolume(endplate_tub1, kMedSteel, "endplate_tub1_logical");

  Logical_volumes.push_back(endplate2_Lvolume);

   endplate2_Lvolume->SetVisAttributes(G4VisAttributes::Invisible);

  G4Tubs* endplate_tub2 = new G4Tubs("endplate2_tub", rInEndPlateI, rOuEndPlateI, (dzEndPlateI + 0.1), 0.*deg,360.*deg);
  
  G4LogicalVolume* endplate3_Lvolume = new G4LogicalVolume(endplate_tub2, kMedSteel, "endplate_tub2_logical");

  Logical_volumes.push_back(endplate3_Lvolume);

   
  

  /*
  
  TGeoTube* endPlate2 = new TGeoTube(0., rInEndPlate, (dzEndPlate + 0.1) / 2.);
  endPlate2->SetName("endPlate2");
  TGeoTube* endPlate3 = new TGeoTube(rInEndPlateI, rOuEndPlateI, (dzEndPlateI + 0.1) / 2.);
  endPlate3->SetName("endPlate3");

  TGeoTranslation* tPlate = new TGeoTranslation("tPlate", 0., 0., -dzEndPlateI - 0.05);
  tPlate->RegisterYourself();

  TGeoCompositeShape* shFaEndPlate = new TGeoCompositeShape("shFaEndPlate", "endPlate1-(endPlate2+endPlate3:tPlate)");
  TGeoVolume* voFaEndPlate = new TGeoVolume("AFaEndPlate", shFaEndPlate, kMedSteel);
*/


//
  // Inner Tungsten Shield
  // Part 1  99.8 cm
  // Part 2 143.5 cm
  // Part 3  25.0 cm
  // Part 4  31.0 cm
  // ====================
  //        299.3 cm - 0.6 overlap between Part 1 and Part 2
  //        298.7 cm
  // Starting position 499.0 - 298.7 = 200.3
  // Within C cone:    200.3 -  92.0 = 108.3 = end of straight section of the Graphite Cone
  //

  // Pos 6
  ///////////////////////////////////
  //    FA Tungsten Tube Part 1    //
  //    Drawing ALIP2A__0045       //
  ///////////////////////////////////


//first part or part1

G4double tungs1_pDz_1 = 98.8*cm /2;
G4double tungs1_pDz_2 = 1.*cm/2;
G4double tungs1_pDz = tungs1_pDz_1 + tungs1_pDz_2;


G4double tungs1_pRmin1_1 =9.1*cm / 2.; 
G4double tungs1_pRmax1_1 = 13.8*cm / 2.;

G4double tungs1_pRmin2_1 = tungs1_pRmin1_1; // inner radius at the end of cone
G4double tungs1_pRmax2_1 = 20.7*cm / 2.; // outer radius at the end of cone


G4Cons* tungs1_cons_1 = new G4Cons("Tungsten_shield_part1_1", tungs1_pRmin1_1, tungs1_pRmax1_1, 
		                   tungs1_pRmin2_1, tungs1_pRmax2_1, tungs1_pDz_1,
				    carbon_pSphi, carbon_pDphi);

G4LogicalVolume* tungs1_Lvolume_1 = new G4LogicalVolume(tungs1_cons_1, kMedNiW, "tungs1_logical_1");

Logical_volumes.push_back(tungs1_Lvolume_1);


G4double tungs1_z_1 = (z_0 + 2*dzCarbonConeS + tungs1_pDz_1) - mag_position;

//final

G4double tungs1_pRmin1= tungs1_pRmin1_1;
G4double tungs1_pRmax1 = 15.0*cm / 2.;
G4double tungs1_pRmin2 = tungs1_pRmin1;
G4double tungs1_pRmax2 = tungs1_pRmax1;

G4double tungs1_z = (z_0 +  2*dzCarbonConeS + 2*tungs1_pDz_1 + tungs1_pDz_2 + 2*dzFaWPlate) - mag_position;

G4Cons* tungs1_cons = new G4Cons("Tungsten_shield_part1", tungs1_pRmin1, tungs1_pRmax1, 
		                   tungs1_pRmin2, tungs1_pRmax2, tungs1_pDz_2,
				    carbon_pSphi, carbon_pDphi);


G4LogicalVolume* tungs1_Lvolume = new G4LogicalVolume(tungs1_cons, kMedNiW, "tungs1_logical");

Logical_volumes.push_back(tungs1_Lvolume);


G4VisAttributes* aa = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
tungs1_Lvolume_1->SetVisAttributes(aa);
tungs1_Lvolume->SetVisAttributes(aa);


 ///////////////////////////////////
  //    FA Tungsten Tube Part 2   //
  //    Drawing ALIP2A__0045       //
  ///////////////////////////////////


// 1

G4double tungs2_pDz_1 = 142.9*cm/2;
G4double tungs2_pDz_2 = 0.6*cm/2;

G4double tungs2_pDz = tungs2_pDz_1 + tungs2_pDz_2;


G4double tungs2_pRmin1_1 = tungs1_pRmin1; 
G4double tungs2_pRmax1_1 = 20.7*cm/2;

G4double tungs2_pRmin2_1 = tungs1_pRmin1;
G4double tungs2_pRmax2_1 =  20.7*cm/2;

 // inner radius at the end of cone



 G4Cons* tungs2_cons_1 = new G4Cons("Tungsten_shield_part1_1", tungs2_pRmin1_1, tungs2_pRmax1_1, 
		                   tungs2_pRmin2_1, tungs2_pRmax2_1, tungs2_pDz_1,
				    carbon_pSphi, carbon_pDphi);

 G4LogicalVolume* tungs2_Lvolume_1 = new G4LogicalVolume(tungs2_cons_1, kMedNiW, "tungs2_logical_1");

 Logical_volumes.push_back(tungs2_Lvolume_1);


G4double tungs2_z_1 = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + tungs2_pDz_1 ) - mag_position;

//2


G4double tungs2_pRmax2 = 30.72*cm / 2. - 0.05*cm; // outer radius at the end of cone
G4double tungs2_pRmin2 = 12.58*cm / 2.;
G4double tungs2_pRmin1 = 9.10*cm / 2.;
G4double tungs2_pRmax1 = tungs2_pRmax1_1;

G4double tungs2_z = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz_1 + tungs2_pDz_2 + 2*dzFaWPlate) - mag_position;

G4Cons* tungs2_cons = new G4Cons("Tungsten_shield_part2_2", tungs2_pRmin1, tungs2_pRmax1, 
		                   tungs2_pRmin2, tungs2_pRmax2, tungs2_pDz_2,
				    carbon_pSphi, carbon_pDphi);

G4LogicalVolume* tungs2_Lvolume = new G4LogicalVolume(tungs2_cons, kMedNiW, "tungs2_logical_2");

Logical_volumes.push_back(tungs2_Lvolume);


tungs2_Lvolume->SetVisAttributes(aa);

tungs2_Lvolume_1->SetVisAttributes(aa);


// Pos 3
  ///////////////////////////////////
  //    FA Flange                  //
  //    Drawing ALIP2A__0038       //
  ///////////////////////////////////
  // Width of the Flange	
  // Outer radius
  G4double rOuFaFlange = 41.0*cm / 2.;
  // 1st section
  G4double dzFaFlange1 = 0.8*cm/2;
  G4double rInFaFlange1 = 33.4*cm / 2.;
  // 2nd section
  G4double dzFaFlange2 = 1.2*cm/2;
  G4double rInFaFlange2 = 36.4*cm / 2.;

  /*

  TGeoPcon* shFaFlange = new TGeoPcon(0., 360., 4);
  z = 0;
  shFaFlange->DefineSection(0, z, rInFaFlange1, rOuFaFlange);
  z += dzFaFlange1;
  shFaFlange->DefineSection(1, z, rInFaFlange1, rOuFaFlange);
  shFaFlange->DefineSection(2, z, rInFaFlange2, rOuFaFlange);
  z += dzFaFlange2;
  shFaFlange->DefineSection(3, z, rInFaFlange2, rOuFaFlange);

  TGeoVolume* voFaFlange = new TGeoVolume("AFaFlange", shFaFlange, kMedSteel);
  */
  G4double flange1_z = (z_0 + dzFaFlange1) - mag_position;  

  G4Cons* flange1_cons = new G4Cons("flange1", rInFaFlange1, rOuFaFlange, 
		                   rInFaFlange1, rOuFaFlange, dzFaFlange1,
				    carbon_pSphi, carbon_pDphi);


  G4LogicalVolume* flange1_Lvolume = new G4LogicalVolume(flange1_cons, kMedSteel, "flange_1");

  Logical_volumes.push_back(flange1_Lvolume);


 G4VisAttributes* flangeCone = new G4VisAttributes(aa);
  flange1_Lvolume->SetVisAttributes(flangeCone);

 G4double flange_z = (z_0 + 2*dzFaFlange1 + dzFaFlange2) - mag_position; 
  
 G4Cons* flange_cons = new G4Cons("flange", rInFaFlange2, rOuFaFlange, 
		                   rInFaFlange2, rOuFaFlange, dzFaFlange2,
				    carbon_pSphi, carbon_pDphi);

G4LogicalVolume* flange_Lvolume = new G4LogicalVolume(flange_cons, kMedSteel, "flange"); 

 Logical_volumes.push_back(flange_Lvolume);

flange_Lvolume->SetVisAttributes(flangeCone);

// Pos 8
  ///////////////////////////////////
  //    FA Tungsten Tube Part 3    //
  //    Drawing ALIP2A__0047       //
  ///////////////////////////////////


G4double tungs3_dz = 25.0*cm/2;
G4double tungs3_pRmin1 = 12.59*cm / 2.;
G4double tungs3_pRmin2 = 13.23*cm / 2.;
G4double tungs3_pRmax1 = 30.60*cm / 2.;
G4double tungs3_pRmax2 = 32.35*cm / 2.;

G4double tungs3_z = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + tungs3_dz ) - mag_position;

G4Cons* tungs3_cons = new G4Cons("Tungsten_shield_part3", tungs3_pRmin1, tungs3_pRmax1, 
		                   tungs3_pRmin2, tungs3_pRmax2, tungs3_dz,
				    carbon_pSphi, carbon_pDphi);





G4LogicalVolume* tungs3_Lvolume = new G4LogicalVolume(tungs3_cons, kMedNiW, "tungs3_logical");

Logical_volumes.push_back(tungs3_Lvolume);

tungs3_Lvolume->SetVisAttributes(aa);


///////////////////////////////////
  //    FA Tungsten Tube Part 4    //
  //    Drawing ALIP2A__0048       //
  ///////////////////////////////////
 

 G4double tungs4_pDz = 31.0*cm/2;
 G4double tungs4_pRmin1 = 13.23*cm / 2.;
 G4double tungs4_pRmin2 = 13.98*cm / 2.;
 G4double tungs4_pRmax1 = 48.80*cm / 2.;
 G4double tungs4_pRmax2 = 52.05*cm/ 2.;

 G4double tungs4_z = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + 2*tungs3_dz + tungs4_pDz ) - mag_position;

G4Cons* tungs4_cons = new G4Cons("Tungsten_shield_part4", tungs4_pRmin1, tungs4_pRmax1, 
		                   tungs4_pRmin2, tungs4_pRmax2, tungs4_pDz,
				    carbon_pSphi, carbon_pDphi);

G4LogicalVolume* tungs4_Lvolume = new G4LogicalVolume(tungs4_cons, kMedNiW, "tungs4_logical");

Logical_volumes.push_back(tungs4_Lvolume);

tungs4_Lvolume->SetVisAttributes(aa);


// Pos 12
  ///////////////////////////////////
  //    FA Lead Cone               //
  //    Drawing ALIP2A__0077       //
  ///////////////////////////////////

// 5 deg cone
  G4double dzFaPbCone5 = 168.9*cm/2;
  G4double rInFaPbCone5 = 37.35*cm / 2.;
  G4double rOuFaPbCone5 = 85.66*cm / 2.;



  // 10 deg cone
  G4double dzFaPbCone10 = 25.9*cm/2;
  G4double  rInFaPbCone10 = rInFaPbCone5 + 2*(dzFaPbCone5*tan(final_angle*PI/180.00));
  G4double  rOuFaPbCone10 = 115.2*cm / 2.;
  // end
  G4double  rInFaPbConeE = 106.05*cm / 2.;
  G4double  rOuFaPbConeE = 124.35*cm / 2.;
  // Total length
  //G4double dzFaPbCone = dzFaPbCone5 + dzFaPbCone10;


  G4Cons* lead_cons_1 = new G4Cons("lead_cons", rInFaPbCone5, rOuFaPbCone5, 
		                   rInFaPbCone10, rOuFaPbCone10, dzFaPbCone5,
				    carbon_pSphi, carbon_pDphi);

   G4double lead_z_1 = (z_0 + 2*dzFaWPlate + dzFaPbCone5 + 2*dzSteelEnvelopeFC) - mag_position;
 
  G4LogicalVolume* lead_Lvolume_1 = new G4LogicalVolume(lead_cons_1, kMedPb, "lead_logical1");

  Logical_volumes.push_back(lead_Lvolume_1);


  G4Cons* lead_cons = new G4Cons("lead_cons", rInFaPbCone10, rOuFaPbCone10, 
		                   rInFaPbConeE,  rOuFaPbConeE, dzFaPbCone10,
				    carbon_pSphi, carbon_pDphi);

  G4LogicalVolume* lead_Lvolume = new G4LogicalVolume(lead_cons, kMedPb, "lead_logical");

   Logical_volumes.push_back(lead_Lvolume);

  G4double lead_z = (z_0 + 2*dzFaWPlate + 2*dzFaPbCone5 + dzFaPbCone10 + 2*dzSteelEnvelopeFC) - mag_position;

   G4Colour  white   (1.0, 1.0, 1.0);
  
  G4VisAttributes* leadCone = new G4VisAttributes(white);
  lead_Lvolume_1->SetVisAttributes(leadCone);
  lead_Lvolume->SetVisAttributes(leadCone);

 //minor part 


 


// Pos 15
  ///////////////////////////////////
  //    FA Steel Plate 250 mm      //
  //    Drawing ALIP2A__00xx       //
  ///////////////////////////////////

G4double eps = 0.001*cm;
G4double SteelCone25_pDz = 25*cm/2;
G4double SteelCone25_Rmin1 = concrete_pRmin2;
G4double SteelCone25_Rmax1 = concrete_pRmax2;
G4double SteelCone25_Rmin2 = SteelCone25_Rmin1 + (2*SteelCone25_pDz)*tan(initial_angle*PI/180.00);
G4double SteelCone25_Rmax2 = SteelCone25_Rmax1 + (2*SteelCone25_pDz)*tan(final_angle*PI/180.00);


G4Cons* Steel25_cons = new G4Cons("Steel25_cons", (SteelCone25_Rmin1 + eps), (SteelCone25_Rmax1 - eps), 
		                   (SteelCone25_Rmin2 + eps), (SteelCone25_Rmax2 - eps), SteelCone25_pDz,
				    carbon_pSphi, carbon_pDphi);


  
G4LogicalVolume* Steel25_Lvolume = new G4LogicalVolume(Steel25_cons, kMedSteel, "Steel25_logical");

Logical_volumes.push_back(Steel25_Lvolume);

G4double Steel25_z = (z_0 + 2*carbon_pDz + 2*concrete_pDz + SteelCone25_pDz )  - mag_position; 


G4Colour  green   (0.0, 1.0, 0.0);
G4VisAttributes* SteelCone = new G4VisAttributes(green);
Steel25_Lvolume->SetVisAttributes(SteelCone);


  ///////////////////////////////////
  //    FA Steel Plate 310 mm      //
  //    Drawing ALIP2A__00xx       //
  ///////////////////////////////////

G4double SteelCone31_pDz = 31*cm/2;
G4double SteelCone31_Rmin1 = 48.80*cm / 2.;
G4double SteelCone31_Rmax1 = SteelCone25_Rmax2;
G4double SteelCone31_Rmin2 = 52.05*cm / 2.;
G4double SteelCone31_Rmax2 = SteelCone31_Rmax1 + (2*SteelCone31_pDz)*tan(final_angle*PI/180.00);


G4Cons* Steel31_cons = new G4Cons("Steel31_cons", (SteelCone31_Rmin1 + eps), (SteelCone31_Rmax1 - eps), 
		                  ( SteelCone31_Rmin2 + eps), (SteelCone31_Rmax2 - eps), SteelCone31_pDz,
				    carbon_pSphi, carbon_pDphi);


G4LogicalVolume* Steel31_Lvolume = new G4LogicalVolume(Steel31_cons, kMedSteel, "Steel31_logical");

Logical_volumes.push_back(Steel31_Lvolume);

G4double Steel31_z = (z_0 + 2*carbon_pDz + 2*concrete_pDz + 2*SteelCone25_pDz + SteelCone31_pDz )  - mag_position; 

Steel31_Lvolume->SetVisAttributes(SteelCone);


// Pos 14
  ///////////////////////////////////
  //    FA Polyethylene Parts      //
  //    Drawing ALIP2A__0034       //
  ///////////////////////////////////
  G4double dzFaCH2Cone = 201.*cm/2;
  G4double rInFaCH2Cone1 = 106.0*cm / 2.;
  G4double rInFaCH2Cone2 = 176.9*cm / 2.;
  G4double dFaCH2Cone = 7.5*cm / cos(10.*PI/180.00);

  
  G4double poly_z = (z_0 + 2*dzFaWPlate + 2*dzFaPbCone5 + 2*dzFaPbCone10 + 2*dzSteelEnvelopeFC + dzFaCH2Cone ) - mag_position;

 
  G4Cons* polyethylene_cons = new G4Cons("polyethylene_cons", rInFaCH2Cone1, (rInFaCH2Cone1 + dFaCH2Cone), 
		                   rInFaCH2Cone2, (rInFaCH2Cone2 + dFaCH2Cone), dzFaCH2Cone,
				    carbon_pSphi, carbon_pDphi);


  G4LogicalVolume* polyethylene_Lvolume = new G4LogicalVolume(polyethylene_cons, KmedCH2Sh, "polyethylene_logical");

 Logical_volumes.push_back(polyethylene_Lvolume);

  G4Colour brown (0.7, 0.4, 0.1);
G4VisAttributes*copperVisAttributes = new G4VisAttributes(brown);
 polyethylene_Lvolume->SetVisAttributes(copperVisAttributes);



///////////////////////////////////
//    FA Tungsten Tail           //
//    Drawing ALIP2A__0049       //
//    Drawing ALIP2A__0111       //
///////////////////////////////////
//
//    The tail as built is shorter than in drawing ALIP2A__0049. 
//    The CDD data base has to be updated !
//
//    Inner radius at the entrance of the flange
     G4double rInFaWTail1  = 13.98*cm/2.;
//    Outer radius at the entrance of the flange
      G4double rOuFaWTail1  = 52.00*cm/2.;
//    Outer radius at the end of the section inside the FA
      G4double rOuFaWTail2  = 35.27*cm/2.;
//    Length of the Flange section inside the FA
      G4double dzFaWTail1   =  6.00*cm/2;
//    Length of the Flange section ouside the FA
      G4double dzFaWTail2   = 12.70*cm/2;
//    Inner radius at the end of the section inside the FA 
      G4double rInFaWTail2  = rInFaWTail1 +  (dzFaWTail1 * tan(0.71*PI/180.00));
//    Inner radius at the end of the flange
      G4double rInFaWTail3  = rInFaWTail2 +  (dzFaWTail2 * tan(0.71*PI/180.00));
//    Outer radius at the end of the flange
      G4double rOuFaWTail3  = rOuFaWTail2 +  (dzFaWTail2 * tan(2.*PI/180.00));
//    Outer radius of the recess for station 1
      G4double rOuFaWTailR  = 30.8*cm/2.;
//    Length of the recess
      G4double dzFaWTailR   = 36.00*cm/2;
//    Inner radiues at the end of the recess      
      G4double rInFaWTail4  =  rInFaWTail3 +  (dzFaWTailR * tan(0.71*PI/180.00));
//    Outer radius at the end of the recess      
      G4double rOuFaWTail4  =  rOuFaWTail3 +  (dzFaWTailR * tan(2.0*PI/180.00));
//    Inner radius of the straight section
      G4double rInFaWTailS  = 22.30*cm/2.;
//    Length of the bulge
      G4double dzFaWTailB   = 13.0*cm/2;
//    Outer radius at the end of the bulge
      G4double rOuFaWTailB  =  rOuFaWTail4 +  (dzFaWTailB * tan(2.0*PI/180.00));
//    Outer radius at the end of the tail 
      G4double rOuFaWTailE  = 31.6*cm/2.;
//    Total length of the tail
      const G4double dzFaWTail    = 70.7*cm/2;

     // TGeoPcon* shFaWTail = new TGeoPcon(0., 360., 10);
      //z    = 0.;
//    Flange section inside FA

     G4Cons* wtail1_cons = new G4Cons("wtail1_cons", rInFaWTail1, rOuFaWTail1, 
		                   rInFaWTail2, rOuFaWTail2, dzFaWTail1,
				    carbon_pSphi, carbon_pDphi);


     G4LogicalVolume* wtail1_Lvolume = new G4LogicalVolume(wtail1_cons, kMedNiW, "wtail1_logical");

     Logical_volumes.push_back(wtail1_Lvolume);

     G4double wtail1_z = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + 2*tungs3_dz + 2*tungs4_pDz + dzFaWTail1) - mag_position;

//    Flange section outside FA
   
    G4Cons* wtail2_cons = new G4Cons("wtail2_cons", rInFaWTail2, rOuFaWTail2, 
		                   rInFaWTail3, rOuFaWTail3, dzFaWTail2,
				    carbon_pSphi, carbon_pDphi);


     G4LogicalVolume* wtail2_Lvolume = new G4LogicalVolume(wtail2_cons, kMedNiW, "wtail2_logical");
 
     Logical_volumes.push_back(wtail2_Lvolume);
  
    G4double wtail2_z = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + 2*tungs3_dz + 2*tungs4_pDz + 2*dzFaWTail1 + dzFaWTail2) - mag_position;
    

      
//    Recess Station 1
     G4Cons* wtail3_cons = new G4Cons("wtail3_cons", rInFaWTail3, rOuFaWTailR, 
		                   rInFaWTail4, rOuFaWTailR, dzFaWTailR,
				    carbon_pSphi, carbon_pDphi);


     G4LogicalVolume* wtail3_Lvolume = new G4LogicalVolume(wtail3_cons, kMedNiW, "wtail3_logical");

    Logical_volumes.push_back(wtail3_Lvolume);
  
    G4double wtail3_z = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + 2*tungs3_dz + 2*tungs4_pDz + 2*dzFaWTail1 + 2*dzFaWTail2 + dzFaWTailR) - mag_position;
    

    
     
      
//    Bulge
   

  G4Cons* wtail4_cons = new G4Cons("wtail4_cons", rInFaWTailS, rOuFaWTail4, 
		                   rInFaWTailS, rOuFaWTailB, dzFaWTailB,
				    carbon_pSphi, carbon_pDphi);


   G4LogicalVolume* wtail4_Lvolume = new G4LogicalVolume(wtail4_cons, kMedNiW, "wtail4_logical");
  
   Logical_volumes.push_back(wtail4_Lvolume);
  
   G4double wtail4_z = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + 2*tungs3_dz + 2*tungs4_pDz + 2*dzFaWTail1 + 2*dzFaWTail2 + 2*dzFaWTailR + dzFaWTailB) - mag_position;



      
//    End
    //  z =  dzFaWTail;

      G4double tail_z = dzFaWTail - dzFaWTailB - dzFaWTailR - dzFaWTail2 - dzFaWTail1;

    	
   
       G4Cons* wtail5_cons = new G4Cons("wtail5_cons", rInFaWTailS, rOuFaWTailE, 
		                   rInFaWTailS, rOuFaWTailE, tail_z,
				    carbon_pSphi, carbon_pDphi);


   G4LogicalVolume* wtail5_Lvolume = new G4LogicalVolume(wtail5_cons, kMedNiW, "wtail5_logical");

  Logical_volumes.push_back(wtail5_Lvolume);
  
   G4double wtail5_z = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + 2*tungs3_dz + 2*tungs4_pDz + 2*dzFaWTail1 + 2*dzFaWTail2 + 2*dzFaWTailR + 2*dzFaWTailB + tail_z) - mag_position;


//TGeoVolume* voFaWTail = new TGeoVolume("YFaWTail", shFaWTail, kMedNiW);




//--------------STRUCTURE SECTION ----------------------------------
 
build_abs = true;

if (build_abs) {




new G4PVPlacement(0,
		 G4ThreeVector(0, 0, concrete_z),
		 concrete_Lvolume,
		 "concrete_cone",
		 logicMag,
		 false,
		 1);


new G4PVPlacement(0,
		    G4ThreeVector(0,0,carbon1_z),
		    carbon1_Lvolume,
		    "carbon1_cone",
		    logicMag,
		    false,
		    1
		    );




new G4PVPlacement(0,
		    G4ThreeVector(0,0,carbon_z),
		    carbon_Lvolume,
		    "carbon_cone",
		    logicMag,
		    false,
		    1
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,tungs1_z_1),
		    tungs1_Lvolume_1,
		    "tungs1_cone_1",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,tungs1_z),
		    tungs1_Lvolume,
		    "tungs1_cone",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,tungs2_z_1),
		    tungs2_Lvolume_1,
		    "tungs2_cone_1",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,tungs2_z),
		    tungs2_Lvolume,
		    "tungs2_cone",
		    logicMag,
		    false,
		    1
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,tungs3_z),
		    tungs3_Lvolume,
		    "tungs3_cone",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,tungs4_z),
		    tungs4_Lvolume,
		    "tungs4_cone",
		    logicMag,
		    false,
		    1
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,Steel25_z),
		    Steel25_Lvolume,
		    "Steel25_cone",
		    logicMag,
		    false,
		    1
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,Steel31_z),
		    Steel31_Lvolume,
		    "Steel31_cone",
		    logicMag,
		    false,
		    1
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,poly_z),
		    polyethylene_Lvolume,
		    "polyethylene_cone",
		    logicMag,
		    false,
		    1);



new G4PVPlacement(0,
		    G4ThreeVector(0,0,lead_z_1),
		    lead_Lvolume_1,
		    "lead_1",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,lead_z),
		    lead_Lvolume,
		    "lead",
		    logicMag,
		    false,
		    1
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,flange1_z),
		    flange1_Lvolume,
		    "flange1",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,flange_z),
		    flange_Lvolume,
		    "flange",
		    logicMag,
		    false,
		    1
		    );




new G4PVPlacement(0,
		    G4ThreeVector(0,0,sEnvelope1_z),
		    sEnvelope1_Lvolume,
		    "sEnvelope1",
		    logicMag,
		    false,
		    1
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,sEnvelope2_z),
		    sEnvelope2_Lvolume,
		    "sEnvelope2",
		    logicMag,
		    false,
		    1
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,sEnvelope3_z),
		    sEnvelope3_Lvolume,
		    "sEnvelope3",
		    logicMag,
		    false,
		    1
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,sEnvelope4_z),
		    sEnvelope4_Lvolume,
		    "sEnvelope4",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,sEnvelopeT_z),
		    sEnvelopeT_Lvolume,
		    "sEnvelope1T",
		    logicMag,
		    false,
		    1
		    );




new G4PVPlacement(0,
		    G4ThreeVector(0,0,plateAB1_z),
		    plateAB1_Lvolume,
		    "plateAB1",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,plateAB2_z),
		    plateAB2_Lvolume,
		    "plateAB2",
		    logicMag,
		    false,
		    1
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,plateAB3_z),
		    plateAB3_Lvolume,
		    "plateAB3",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,plateAB4_z),
		    plateAB4_Lvolume,
		    "plateAB4",
		    logicMag,
		    false,
		    1
		    );







new G4PVPlacement(0,
		    G4ThreeVector(0,0,endplate_z),
		    endplate_Lvolume,
		    "endplate",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,0),
		    endplate2_Lvolume,
		    "endplate",
		    endplate_Lvolume,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,0),
		    endplate3_Lvolume,
		    "endplate",
		    endplate_Lvolume,
		    false,
		    1
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,wtail1_z),
		    wtail1_Lvolume,
		    "wtail1",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,wtail2_z),
		    wtail2_Lvolume,
		    "wtail2",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,wtail3_z),
		    wtail3_Lvolume,
		    "wtail3",
		    logicMag,
		    false,
		    1
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,wtail4_z),
		    wtail4_Lvolume,
		    "wtail4",
		    logicMag,
		    false,
		    1
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,wtail5_z),
		    wtail5_Lvolume,
		    "wtail5",
		    logicMag,
		    false,
		    1
		    );






}








//detector 1
pos_after_detec = 0.01*m;
G4double detec_length = 0.5*cm;
G4double initial_radius = 0;
G4double final_radius = (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec)*tan(final_angle*PI/180.00);

G4Tubs* detec_tub = new G4Tubs("detec_tubs", initial_radius, final_radius, detec_length, 0.*deg,360.*deg);

G4LogicalVolume* detec_volume = new G4LogicalVolume(detec_tub, world_mat, "detec");

G4double detec_z = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + 2*tungs3_dz + 2*tungs4_pDz + 2*dzFaWTail1 + 2*dzFaWTail2 + 2*dzFaWTailR + 2*dzFaWTailB + 2*tail_z + pos_after_detec + detec_length) - mag_position;


 new G4PVPlacement(0,
		 G4ThreeVector(0,0,detec_z),
		 detec_volume,
		 "detector",
		 logicMag,
		 false,
		 0,
		false
		 );




//detector 2

G4Tubs* detec_tub2 = new G4Tubs("detec_tubs", initial_radius, final_radius, detec_length, 0.*deg,360.*deg);
G4LogicalVolume* detec_volume2 = new G4LogicalVolume(detec_tub2, world_mat, "detec2");
G4double detec_z2 = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + 2*tungs3_dz + 2*tungs4_pDz + 2*dzFaWTail1 + 2*dzFaWTail2 + 2*dzFaWTailR + 2*dzFaWTailB + 2*tail_z + pos_after_detec + 1*m + 3*detec_length) - mag_position;

new G4PVPlacement(0,
		 G4ThreeVector(0,0,detec_z2),
		 detec_volume2,
		 "detector2",
		 logicMag,
		 false,
		 0,
		false
		 );

num_detec=2;



 carbon_Lvolume->SetVisAttributes(color );

  fScoringVolume1 = carbon_Lvolume;
  fScoringVolume2 = concrete_Lvolume;



//sensitive detector1

auto sdman = G4SDManager::GetSDMpointer();
G4String SDname1 = "SD1";
auto sensitive1 = new B1SD(SDname1);
sdman->AddNewDetector(sensitive1);
detec_volume->SetSensitiveDetector(sensitive1);



G4String SDname2 = "SD2";
auto sensitive2 = new B1SD(SDname2);
sdman->AddNewDetector(sensitive2);
detec_volume2->SetSensitiveDetector(sensitive2);



  return physWorld;
}

void B1DetectorConstruction::ConstructSDandField() {
  //creating uniform magnetic field
  
    G4MagneticField* magField = new G4UniformMagField(G4ThreeVector(0.,0.,-.5*tesla));
    G4FieldManager* FieldMgr  =new G4FieldManager(magField);
    logicMag->SetFieldManager(FieldMgr, true);
    G4cout << "Applying magnetic field" << G4endl;
  

  


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
