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
#include "G4SubtractionSolid.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4AssemblyVolume.hh"
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





//Area with magnetic field
G4double mag_position = 4*m;
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

//Magnesium
//matmgr.Material("ABSO", 7, "MAGNESIUM$", 24.31, 12., 1.74, 25.3, 46.0);


//G4Material* kMedMg = new G4Material("Magnesium", 12.,24.31*g/mole,1.74*g/mole);
G4Material* kMedMg = nist->FindOrBuildMaterial("G4_Mg");



//some importante initial values

 # define PI 3.14159265
 
  G4double kDegRad = PI/180.00;

  // Mimimum angle of the tracking region
  G4double angle02 = tan(2.*kDegRad);
  // Opening angle of W rear plug
  G4double angle03 = tan(3.*kDegRad);
  //
  G4double angle05 = tan(5.*kDegRad);
  // Maximum angle of the tracking region
  G4double angle10 = tan(10.*kDegRad);  
  // Opening angle of the FA snout
  G4double angle24 = tan(24.*kDegRad);
  // Opneing angle of the inner cone
  G4double angle71 = tan(0.697*kDegRad);

  // Starting position in z   const Float_t zFa = 90.0;
  G4double zFa = 90.0*cm;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  G4double pos_after_detec = 3.*m;
  G4double angle0 = 0.*deg;
  G4double angle360 = 360.*deg;
  G4Transform3D Tr;

  G4ThreeVector Ta;

  G4RotationMatrix Ra;

  G4double z;

 


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 ///////////////////////////////////
  //                               //
  //        Front Absorber         //
  //        Drawing ALIP2A__0106   //
  //                               //
  //                               //
  ///////////////////////////////////
  //
  // Pos  1 Steel Envelope
  // Pos  2 End Plate
  // Pos  3 Flange (wrong arrow in the drawing)
  // Pos  4 W Plate A
  // Pos  5 W Plate B
  // Pos  6 Tungsten Tube Part 1
  // Pos  7 Tungsten Tube Part 2
  // Pos  8 Tungsten Tube Part 3
  // Pos  9 Tungsten Tube Part 4
  // Pos 10 Tungsten Tail
  // Pos 11 Graphite Cone
  // Pos 12 Pb       Cone
  // Pos 13 Concrete Cone
  // Pos 14 Polyethylene Parts
  // Pos 15 Steel Plate 25 cm
  // Pos 16 Steel Plate 31 cm
  // Pos 17 Magnesium Ring
  // Pos 18 Composite Ring
  //
  //
  /*
  // Mimimum angle of the tracking region
  const Float_t angle02 = TMath::Tan(2. * kDegRad);
  // Maximum angle of the tracking region
  const Float_t angle10 = TMath::Tan(10. * kDegRad);
  // Opening angle of W rear plug
  const Float_t angle03 = TMath::Tan(3. * kDegRad);
  //
  const Float_t angle05 = TMath::Tan(5. * kDegRad);
  // Opening angle of the FA snout
  const Float_t angle24 = TMath::Tan(24. * kDegRad);
  // Opneing angle of the inner cone
  const Float_t angle71 = TMath::Tan(0.697 * kDegRad);
  // Starting position in z
  const Float_t zFa = 90.0;
  */



//Construction of solids
   
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
  G4double rInSteelEnvelopeFC2 = rInSteelEnvelopeFC1 + (2*dzSteelEnvelopeFC * angle10);
  // Outer Radius
  G4double rOuSteelEnvelopeFC1 = 88.97*cm / 2.;
  G4double rOuSteelEnvelopeFC2 = rOuSteelEnvelopeFC1 + (2*dzSteelEnvelopeFC * angle05);
  //
  // 5 deg cone
  G4double dzSteelEnvelopeC5 = 168.9*cm/2;
  G4double rInSteelEnvelopeC5 = rOuSteelEnvelopeFC2 - (2*dSteelEnvelope / cos(5. * kDegRad));
  G4double rOuSteelEnvelopeC5 = rOuSteelEnvelopeFC2;
  // 10 deg cone
  G4double dzSteelEnvelopeC10 = (227.1*cm - 4.*cm)/2;
  G4double rInSteelEnvelopeC10 = 116.22*cm / 2.;
  G4double rOuSteelEnvelopeC10 = rInSteelEnvelopeC10 + (2*dSteelEnvelope / cos(10 * kDegRad));
  // Rear ring
  G4double dzSteelEnvelopeR = 4.*cm/2;
  G4double rInSteelEnvelopeR2 = 196.3*cm / 2.;
  G4double rOuSteelEnvelopeR2 = 212.0*cm / 2.;
  G4double rInSteelEnvelopeR1 = rInSteelEnvelopeR2 - (2*dzSteelEnvelopeR * angle10);
  G4double rOuSteelEnvelopeR1 = rInSteelEnvelopeR1 + (2*dSteelEnvelope / cos(10 *kDegRad));
  // Front insert

  G4double dzSteelEnvelopeFI = 1.*cm/2;
  G4double rInSteelEnvelopeFI = 42.0*cm / 2.;
  G4double rOuSteelEnvelopeFI = (85.0*cm / 2.) + 0.06*cm;

  
  

  

  // Front cover
  

  G4Cons* shFaSteelEnvelopeC1_cone = new G4Cons("shFaSteelEnvelopeC1_cone", rInSteelEnvelopeFC1, rOuSteelEnvelopeFC1,  rInSteelEnvelopeFC2, rOuSteelEnvelopeFC2, dzSteelEnvelopeFC, angle0, angle360);
 

  G4LogicalVolume* shFaSteelEnvelopeC1 = new G4LogicalVolume(shFaSteelEnvelopeC1_cone, kMedSteel, "shFaSteelEnvelopeC1");  

  Logical_volumes.push_back(shFaSteelEnvelopeC1);

  // 5 deg cone
  
  G4Cons* shFaSteelEnvelopeC2_cone = new G4Cons("shFaSteelEnvelopeC2_cone", rInSteelEnvelopeC5, rOuSteelEnvelopeC5, 
		                   rInSteelEnvelopeC10, rOuSteelEnvelopeC10, dzSteelEnvelopeC5,
				    angle0, angle360);

  G4LogicalVolume* shFaSteelEnvelopeC2 = new G4LogicalVolume(shFaSteelEnvelopeC2_cone, kMedSteel, "shFaSteelEnvelopeC2");

  Logical_volumes.push_back(shFaSteelEnvelopeC2);

  // 10 deg cone
  
  G4Cons* shFaSteelEnvelopeC3_cone = new G4Cons("shFaSteelEnvelopeC3_cone", rInSteelEnvelopeC10, rOuSteelEnvelopeC10, 
		                   rInSteelEnvelopeR1, rOuSteelEnvelopeR1, dzSteelEnvelopeC10,
				    angle0, angle360);

  G4LogicalVolume* shFaSteelEnvelopeC3 = new G4LogicalVolume(shFaSteelEnvelopeC3_cone, kMedSteel, "shFaSteelEnvelopeC3");

  Logical_volumes.push_back(shFaSteelEnvelopeC3);

 

  // Rear Ring
  

  G4Cons* shFaSteelEnvelopeC4_cone = new G4Cons("shFaSteelEnvelopeC4_cone", rInSteelEnvelopeR1, rOuSteelEnvelopeR2, 
		                   rInSteelEnvelopeR2, rOuSteelEnvelopeR2, dzSteelEnvelopeR,
				    angle0, angle360);

  G4LogicalVolume* shFaSteelEnvelopeC4 = new G4LogicalVolume(shFaSteelEnvelopeC4_cone, kMedSteel, "shFaSteelEnvelopeC4");

  Logical_volumes.push_back(shFaSteelEnvelopeC4);
 
  
  // Insert


  G4Tubs* shFaSteelEnvelopeT_tub = new G4Tubs("shFaSteelEnvelopeT_tub", rInSteelEnvelopeFI, rOuSteelEnvelopeFI, dzSteelEnvelopeFI, 0.*deg,360.*deg);
  
 G4LogicalVolume* shFaSteelEnvelopeT = new G4LogicalVolume(shFaSteelEnvelopeT_tub, kMedSteel, "shFaSteelEnvelopeT");//FALTA MATERIAL
 
 Logical_volumes.push_back(shFaSteelEnvelopeT);

  
 G4AssemblyVolume* voFaSteelEnvelope = new G4AssemblyVolume();

  z = dzSteelEnvelopeFC;
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaSteelEnvelope->AddPlacedVolume(shFaSteelEnvelopeC1, Tr);


  z += dzSteelEnvelopeFC + dzSteelEnvelopeC5;
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaSteelEnvelope->AddPlacedVolume(shFaSteelEnvelopeC2, Tr);

  z +=  dzSteelEnvelopeC5 + dzSteelEnvelopeC10;
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaSteelEnvelope->AddPlacedVolume(shFaSteelEnvelopeC3, Tr);

  z += dzSteelEnvelopeC10 + dzSteelEnvelopeR;
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaSteelEnvelope->AddPlacedVolume(shFaSteelEnvelopeC4, Tr);



//setting color to the solids

G4VisAttributes* color = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));

shFaSteelEnvelopeC1->SetVisAttributes(color);
shFaSteelEnvelopeC2->SetVisAttributes(color);
shFaSteelEnvelopeC3->SetVisAttributes(color);
shFaSteelEnvelopeC4->SetVisAttributes(color);
shFaSteelEnvelopeT->SetVisAttributes(color);


  



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


  G4Box* endPlate1 =
    new G4Box("endplate_box",                       //its name
      dxEndPlate, dyEndPlate, dzEndPlate);


   
  G4Tubs* endPlate2 = new G4Tubs("endPlate2", 0., rInEndPlate, (dzEndPlate+ 0.1*cm), 0.*deg,360.*deg);
  

  G4Tubs* endPlate3 = new G4Tubs("endPlate3", rInEndPlateI, rOuEndPlateI, (dzEndPlateI + 0.1*cm), 0.*deg,360.*deg);
  


  

  G4VSolid* hole_square0 = new G4SubtractionSolid("Box+Cylinder", endPlate1, endPlate2);  

  
  G4VSolid* shFaEndPlate_solid = new G4SubtractionSolid("shFaEndPlate", hole_square0, endPlate3, 0,G4ThreeVector(0,0, (-2*dzEndPlateI -0.05*cm)));  


  G4LogicalVolume* voFaEndPlate = new G4LogicalVolume(shFaEndPlate_solid, kMedSteel, "shFaEndPlate");
   
 
 //setting color to the solid
 voFaEndPlate->SetVisAttributes(gray);



// Pos 3
  ///////////////////////////////////
  //    FA Flange                  //
  //    Drawing ALIP2A__0038       //
  ///////////////////////////////////

  // Width of the Flange
  G4double dzFaFlange = 2.*cm/2;
 	
  // Outer radius
  G4double rOuFaFlange = 41.0*cm / 2.;
  // 1st section
  G4double dzFaFlange1 = 0.8*cm/2;
  G4double rInFaFlange1 = 33.4*cm / 2.;
  // 2nd section
  G4double dzFaFlange2 = 1.2*cm/2;
  G4double rInFaFlange2 = 36.4*cm / 2.;



  G4Cons* shFaFlange1_cone = new G4Cons("shFaFlange1_cone", rInFaFlange1, rOuFaFlange, 
		                   rInFaFlange1, rOuFaFlange, dzFaFlange1,
				    angle0, angle360);


  G4LogicalVolume* shFaFlange1 = new G4LogicalVolume(shFaFlange1_cone, kMedSteel, "shFaFlange1");

  Logical_volumes.push_back(shFaFlange1);

  
 G4Cons* shFaFlange2_cone = new G4Cons("shFaFlange2_cone", rInFaFlange2, rOuFaFlange, 
		                   rInFaFlange2, rOuFaFlange, dzFaFlange2,
				   angle0, angle360);

G4LogicalVolume* shFaFlange2 = new G4LogicalVolume(shFaFlange2_cone, kMedSteel, "flange"); 

 Logical_volumes.push_back(shFaFlange2);

 
G4AssemblyVolume* voFaFlange = new G4AssemblyVolume();

z = dzFaFlange1;
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaFlange->AddPlacedVolume(shFaFlange1, Tr);

z += dzFaFlange1 + dzFaFlange2;
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaFlange->AddPlacedVolume(shFaFlange2, Tr);



//setting color to solids
shFaFlange1->SetVisAttributes(gray);
shFaFlange2->SetVisAttributes(gray);


 // Pos 4+5
  ///////////////////////////////////
  //    FA W Plate A+B             //
  //    Drawing ALIP2A__0043       //
  ///////////////////////////////////
  // Front Flange
  G4double dzFaWPlateF = 2.00*cm/2.;
  G4double rInFaQPlateF = 20.50*cm;
  G4double rOuFaQPlateF = 40.05*cm;
  // 1st Central Part 24 deg
  G4double dzFaWPlateC1 = 7.95*cm/2;
  G4double rInFaQPlateC1 = 16.35*cm;
  G4double rOuFaQPlateC1 = rOuFaQPlateF + (2*dzFaWPlateF * angle24);
  // 2nd Central Part 5 deg
  G4double dzFaWPlateC2 = 1.05*cm/2;
  G4double rInFaQPlateC2 = rInFaQPlateC1 + (2*dzFaWPlateC1 * angle10);
  G4double rOuFaQPlateC2 = rOuFaQPlateC1 + (2*dzFaWPlateC1 * angle24);
  G4double rInFaQPlateC3 = 17.94*cm;
  G4double rOuFaQPlateC3 = 44.49*cm;
  // Rear Flange
  G4double dzFaWPlateR = 1.00*cm/2;
  G4double rInFaQPlateR = 21.00*cm;
  G4double rOuFaQPlateR = 42.55*cm;
  // Lenth of Plate - Rear Flange
  G4double dzFaWPlate = dzFaWPlateF + dzFaWPlateC1 + dzFaWPlateC2;
  



  G4Cons* shFaWPlateA1_cone = new G4Cons("shFaWPlateA1_cone", rInFaQPlateF, rOuFaQPlateF, 
		                   rInFaQPlateF, rOuFaQPlateC1, dzFaWPlateF,
				   angle0, angle360);

  G4LogicalVolume* shFaWPlateA1 = new G4LogicalVolume(shFaWPlateA1_cone, kMedNiW, "shFaWPlateA1");

  Logical_volumes.push_back(shFaWPlateA1);


  // 24 deg cone

  G4Cons* shFaWPlateA2_cone = new G4Cons("shFaWPlateA2_cone", rInFaQPlateC1, rOuFaQPlateC1, 
		                   rInFaQPlateC2, rOuFaQPlateC2, dzFaWPlateC1,
				   angle0, angle360);

  G4LogicalVolume* shFaWPlateA2 = new G4LogicalVolume(shFaWPlateA2_cone, kMedNiW, "shFaWPlateA2");

  Logical_volumes.push_back(shFaWPlateA2);

 

  
  // 5 deg cone
  G4Cons* shFaWPlateA3_cone = new G4Cons("shFaWPlateA3_cone", rInFaQPlateC2, rOuFaQPlateC2, 
		                   rInFaQPlateC3, rOuFaQPlateC3, dzFaWPlateC2,
				    angle0, angle360);

  G4LogicalVolume* shFaWPlateA3 = new G4LogicalVolume(shFaWPlateA3_cone, kMedNiW, "shFaWPlateA3");
  
  Logical_volumes.push_back(shFaWPlateA3);

  
  // Rear Flange
  G4Cons* shFaWPlateA4_cone = new G4Cons("shFaWPlateA4_cone", rInFaQPlateR, rOuFaQPlateR, 
		                   rInFaQPlateR, rOuFaQPlateR, dzFaWPlateR,
				    angle0, angle360);

  G4LogicalVolume* shFaWPlateA4 = new G4LogicalVolume(shFaWPlateA4_cone, kMedNiW, "shFaWPlateA4");
  
  Logical_volumes.push_back(shFaWPlateA4);




  G4AssemblyVolume* voFaWPlateteA = new G4AssemblyVolume();

  z = dzFaWPlateF;
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaWPlateteA->AddPlacedVolume(shFaWPlateA1, Tr);

  z+= dzFaWPlateF + dzFaWPlateC1; 
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaWPlateteA->AddPlacedVolume(shFaWPlateA2, Tr);

  z+=  dzFaWPlateC1 + dzFaWPlateC2; 
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaWPlateteA->AddPlacedVolume(shFaWPlateA3, Tr);

  z+=  dzFaWPlateC2 + dzFaWPlateR; 
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaWPlateteA->AddPlacedVolume(shFaWPlateA4, Tr);

 
  
 //setting color  
  
 shFaWPlateA1->SetVisAttributes(gray);
 shFaWPlateA2->SetVisAttributes(gray);
 shFaWPlateA3->SetVisAttributes(gray);
 shFaWPlateA4->SetVisAttributes(gray);


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

//
  // Inner radius
  G4double rInFaWTube1C1 = 9.1*cm / 2.;
  // Central part
  G4double dzFaWTube1C = 98.8*cm/2;
  G4double rOuFaWTube1C1 = 13.8*cm / 2.;
  G4double rOuFaWTube1C2 = 20.7*cm / 2.;
  // Rear Flange
  G4double dzFaWTube1R = 1.0*cm/2;
  G4double rOuFaWTube1R = 15.0*cm / 2.;
  // Total length
  G4double dzFaWTube1 = dzFaWTube1C + dzFaWTube1R;


// Central Part

G4Cons* shFaWTube1_1cone = new G4Cons("shFaWTube1_1cone", rInFaWTube1C1, rOuFaWTube1C1, 
		                   rInFaWTube1C1, rOuFaWTube1C2, dzFaWTube1C,
				   angle0, angle360);

G4LogicalVolume* shFaWTube1_1 = new G4LogicalVolume(shFaWTube1_1cone, kMedNiW, "shFaWTube1_1");

Logical_volumes.push_back(shFaWTube1_1);



// Rear Flange


G4Cons* shFaWTube1_2cone = new G4Cons("shFaWTube1_2cone", rInFaWTube1C1, rOuFaWTube1R, 
		                   rInFaWTube1C1, rOuFaWTube1R, dzFaWTube1R,
				    angle0, angle360);


G4LogicalVolume* shFaWTube1_2 = new G4LogicalVolume(shFaWTube1_2cone, kMedNiW, "shFaWTube1_2");

Logical_volumes.push_back(shFaWTube1_2);


//setting color

G4VisAttributes* aa = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
shFaWTube1_1->SetVisAttributes(aa);
shFaWTube1_2->SetVisAttributes(aa);


// Pos 7
///////////////////////////////////
  //    FA Tungsten Tube Part 2   //
  //    Drawing ALIP2A__0045       //
  ///////////////////////////////////

//

  // Central part
  G4double dzFaWTube2C = 142.9*cm/2;
  G4double rInFaWTube2C1 = 9.10*cm / 2.;
  G4double rInFaWTube2C2 = 12.58*cm / 2.;
  G4double rOuFaWTube2C1 = 20.70*cm / 2.;
  G4double rOuFaWTube2C2 = (30.72*cm / 2.) - 0.05*cm;
  // Front Flange
  G4double dzFaWTube2F = 0.6*cm/2;
  G4double rInFaWTube2F = 15.4*cm / 2.;
  // Total length
  G4double dzFaWTube2 = dzFaWTube2C + dzFaWTube2F;


// Front Flange


 G4Cons* shFaWTube2_1cone = new G4Cons("shFaWTube2_1cone", rInFaWTube2F, rOuFaWTube2C1, 
		                   rInFaWTube2F, rOuFaWTube2C1, dzFaWTube2F,
				    angle0, angle360);

 G4LogicalVolume* shFaWTube2_1 = new G4LogicalVolume(shFaWTube2_1cone, kMedNiW, "shFaWTube2_1");

 Logical_volumes.push_back(shFaWTube2_1);


// Central part


G4Cons* shFaWTube2_2cone = new G4Cons("shFaWTube2_2cone", rInFaWTube2C1, rOuFaWTube2C1, 
		                   rInFaWTube2C2, rOuFaWTube2C2, dzFaWTube2C,
				    angle0, angle360);

G4LogicalVolume* shFaWTube2_2 = new G4LogicalVolume(shFaWTube2_2cone, kMedNiW, "shFaWTube2_2");

Logical_volumes.push_back(shFaWTube2_2);


//setting color

shFaWTube2_2->SetVisAttributes(aa);

shFaWTube2_1->SetVisAttributes(aa);


// Pos 8
  ///////////////////////////////////
  //    FA Tungsten Tube Part 3    //
  //    Drawing ALIP2A__0047       //
  ///////////////////////////////////


G4double dzFaWTube3 = 25.0*cm/2;
G4double rInFaWTube3C1 = 12.59*cm / 2.;
G4double rInFaWTube3C2 = 13.23*cm / 2.;
G4double rOuFaWTube3C1 = 30.60*cm / 2.;
G4double rOuFaWTube3C2 = 32.35*cm / 2.;


G4Cons* voFaWTube3_cone = new G4Cons("voFaWTube3_cone", rInFaWTube3C1, rOuFaWTube3C1, 
		                   rInFaWTube3C2, rOuFaWTube3C2, dzFaWTube3,
				   angle0, angle360);


G4LogicalVolume* voFaWTube3 = new G4LogicalVolume(voFaWTube3_cone, kMedNiW, "voFaWTube3");

Logical_volumes.push_back(voFaWTube3);

//setting color
voFaWTube3->SetVisAttributes(aa);

// Pos 9	
///////////////////////////////////
  //    FA Tungsten Tube Part 4    //
  //    Drawing ALIP2A__0048       //
  ///////////////////////////////////
 

 G4double dzFaWTube4 = 31.0*cm/2;
 G4double rInFaWTube4C1 = 13.23*cm / 2.;
 G4double rInFaWTube4C2 = 13.98*cm / 2.;
 G4double rOuFaWTube4C1 = 48.80*cm / 2.;
 G4double rOuFaWTube4C2 = 52.05*cm/ 2.;


G4Cons* voFaWTube4_cone = new G4Cons("voFaWTube4_cone", rInFaWTube4C1, rOuFaWTube4C1, 
		                   rInFaWTube4C2, rOuFaWTube4C2, dzFaWTube4,
				    angle0, angle360);

G4LogicalVolume* voFaWTube4 = new G4LogicalVolume(voFaWTube4_cone, kMedNiW, "voFaWTube4");

Logical_volumes.push_back(voFaWTube4);

voFaWTube4->SetVisAttributes(aa);


//pos 10
//
// This section is partially in Shield

G4double dzFaWTube5 = 6.0*cm/2;
G4double kAngle0071 = tan(0.71 * kDegRad);
G4double rInFaWTube5C1 = rInFaWTube4C2;
G4double rInFaWTube5C2 = rInFaWTube4C2 + (2*dzFaWTube4 * kAngle0071);
G4double rOuFaWTube5C1 = rOuFaWTube4C2;



G4Cons* voFaWTube5_cone = new G4Cons("voFaWTube5_cone", rInFaWTube5C1, rOuFaWTube5C1, 
		                   rInFaWTube5C2, rOuFaWTube5C1, dzFaWTube5,
				   angle0, angle360);


G4LogicalVolume* voFaWTube5 = new G4LogicalVolume(voFaWTube5_cone, kMedNiW, "voFaWTube5");



  //
  // Pos 11
  ///////////////////////////////////
  //    FA Graphite Cone           //
  //    Drawing ALIP2_0002         //
  ///////////////////////////////////
  //
  // Total length
 G4double dzFaGraphiteCone = 225.0*cm/2;
  // Straight section = start of the 2deg inner cone
  G4double dzFaGraphiteConeS = 108.3*cm/2;
  // Inner radius at the front
 G4double rInFaGraphiteCone1 = 4.5*cm;
  // Outer radius at the front
 G4double rOuFaGraphiteCone1 = (zFa + 2*dzFaFlange) * angle10;
  // Inner radius at start of inner opening cone
 G4double rInFaGraphiteCone2 = 7.0*cm;
  // Outer radius at start of inner opening cone
 G4double rOuFaGraphiteCone2 = (zFa + 2*dzFaFlange + 2*dzFaGraphiteConeS) * angle10;
  // Inner radius the rear
 G4double rInFaGraphiteCone3 = 11.0*cm;
  // Ouer radius at the rear
 G4double rOuFaGraphiteCone3 = (zFa + 2*dzFaFlange + 2*dzFaGraphiteCone) * angle10;
  
  // Straight section
  
  G4Cons* shFaGraphiteCone1_cone = new G4Cons("shFaGraphiteCone1_cone", rInFaGraphiteCone1,  rOuFaGraphiteCone1, 
		                   rInFaGraphiteCone1, rOuFaGraphiteCone2, dzFaGraphiteConeS,
				   angle0, angle360);

  G4LogicalVolume* shFaGraphiteCone1 = new G4LogicalVolume(shFaGraphiteCone1_cone, kMedCSh, "shFaGraphiteCone1");

 
  Logical_volumes.push_back(shFaGraphiteCone1);


// 2 deg opening cone

G4Cons* shFaGraphiteCone2_cone  = new G4Cons("shFaGraphiteCone2_cone ", rInFaGraphiteCone2, rOuFaGraphiteCone2, 
		                   rInFaGraphiteCone3, rOuFaGraphiteCone3, (dzFaGraphiteCone - dzFaGraphiteConeS),
				   angle0, angle360);



G4LogicalVolume* shFaGraphiteCone2 = new G4LogicalVolume(shFaGraphiteCone2_cone, kMedCSh, "shFaGraphiteCone2");


Logical_volumes.push_back(shFaGraphiteCone2);



 //setting color

 G4VisAttributes* blue = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
 shFaGraphiteCone1->SetVisAttributes(blue);
 shFaGraphiteCone2->SetVisAttributes(blue);


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
  G4double  rInFaPbCone10 = rInFaPbCone5 + (2*dzFaPbCone5*angle10);
  G4double  rOuFaPbCone10 = 115.2*cm / 2.;
// end
  G4double  rInFaPbConeE = 106.05*cm / 2.;
  G4double  rOuFaPbConeE = 124.35*cm / 2.;
//Total length
  G4double dzFaPbCone = dzFaPbCone5 + dzFaPbCone10;

	
 // 5 deg cone
  G4Cons* shFaPbCone1_cone = new G4Cons("shFaPbCone1_cones", rInFaPbCone5, rOuFaPbCone5, 
		                   rInFaPbCone10, rOuFaPbCone10, dzFaPbCone5,
				    angle0, angle360);
 
  G4LogicalVolume* shFaPbCone1 = new G4LogicalVolume(shFaPbCone1_cone, kMedPb, "shFaPbCone1");

  Logical_volumes.push_back(shFaPbCone1);
  
 // 10 deg cone

  G4Cons* shFaPbCone2_cone = new G4Cons("shFaPbCone2_cone", rInFaPbCone10, rOuFaPbCone10, 
		                   rInFaPbConeE,  rOuFaPbConeE, dzFaPbCone10,
				    angle0, angle360);

  G4LogicalVolume* shFaPbCone2 = new G4LogicalVolume(shFaPbCone2_cone, kMedPb, "shFaPbCone2");

   Logical_volumes.push_back(shFaPbCone2);


  
  G4AssemblyVolume* voFaPbCone = new G4AssemblyVolume();

  z = dzFaPbCone5;
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaPbCone->AddPlacedVolume(shFaPbCone1, Tr);

  z += dzFaPbCone5 + dzFaPbCone10;
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaPbCone->AddPlacedVolume(shFaPbCone2, Tr);
  


  //setting color

   G4Colour  white   (1.0, 1.0, 1.0);
  
  G4VisAttributes* leadCone = new G4VisAttributes(white);
  shFaPbCone1->SetVisAttributes(leadCone);
  shFaPbCone2->SetVisAttributes(leadCone);



 // Pos 13
  ///////////////////////////////////
  //    FA Concrete Cone           //
  //    Drawing ALIP2A__00xx       //
  ///////////////////////////////////
  G4double dzFaConcreteCone = 126.*cm/2;
  G4double rOuFaConcreteCone1 = rOuFaGraphiteCone3;
  G4double rInFaConcreteCone1 = 11.*cm;
  G4double rOuFaConcreteCone2 = rOuFaConcreteCone1 + 2*dzFaConcreteCone * angle10;
  G4double rInFaConcreteCone2 = rInFaConcreteCone1 + 2*dzFaConcreteCone * angle02;
 
  
 G4Cons* voFaConcreteCone_cone = new G4Cons("voFaConcreteCone_cone", rInFaConcreteCone1, rOuFaConcreteCone1, rInFaConcreteCone2, rOuFaConcreteCone2, dzFaConcreteCone, angle0, angle360);

 G4LogicalVolume* voFaConcreteCone = new G4LogicalVolume(voFaConcreteCone_cone, kMedConcSh, "voFaConcreteCone");

  Logical_volumes.push_back(voFaConcreteCone);

G4VisAttributes* red = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
voFaConcreteCone->SetVisAttributes(red);


// Pos 14
  ///////////////////////////////////
  //    FA Polyethylene Parts      //
  //    Drawing ALIP2A__0034       //
  ///////////////////////////////////
  G4double dzFaCH2Cone = 201.*cm/2;
  G4double rInFaCH2Cone1 = 106.0*cm / 2.;
  G4double rInFaCH2Cone2 = 176.9*cm / 2.;
  G4double dFaCH2Cone = 7.5*cm / cos(10.*kDegRad);
 
  G4Cons* voFaCH2Cone_cone = new G4Cons("voFaCH2Cone_cone", rInFaCH2Cone1, (rInFaCH2Cone1 + dFaCH2Cone), 
		                   rInFaCH2Cone2, (rInFaCH2Cone2 + dFaCH2Cone), dzFaCH2Cone,
				    angle0, angle360);


  G4LogicalVolume* voFaCH2Cone = new G4LogicalVolume(voFaCH2Cone_cone, KmedCH2Sh, "voFaCH2Cone");

 Logical_volumes.push_back(voFaCH2Cone);

//setting color
  G4Colour brown (0.7, 0.4, 0.1);
 G4VisAttributes*copperVisAttributes = new G4VisAttributes(brown);
 voFaCH2Cone->SetVisAttributes(copperVisAttributes);



// Pos 15
  ///////////////////////////////////
  //    FA Steel Plate 250 mm      //
  //    Drawing ALIP2A__00xx       //
  ///////////////////////////////////

G4double dzFaSteelCone25 = 25.*cm/2;
G4double eps = 0.001*cm;
G4double rInFaSteelCone25A = rInFaConcreteCone2;
G4double rOuFaSteelCone25A = rOuFaConcreteCone2;
G4double rInFaSteelCone25B = rInFaSteelCone25A + 2*dzFaSteelCone25 * angle02;
G4double rOuFaSteelCone25B = rOuFaSteelCone25A + 2*dzFaSteelCone25 * angle10;


G4Cons* voFaSteelCone25_cone = new G4Cons("voFaSteelCone25_cone", (rInFaSteelCone25A + eps), (rOuFaSteelCone25A - eps), 
		                   (rInFaSteelCone25B + eps), (rOuFaSteelCone25B - eps), dzFaSteelCone25,
				    angle0, angle360);


  
G4LogicalVolume* voFaSteelCone25 = new G4LogicalVolume(voFaSteelCone25_cone, kMedSteel, "voFaSteelCone25");

Logical_volumes.push_back(voFaSteelCone25);

//setting color

G4Colour  green   (0.0, 1.0, 0.0);
G4VisAttributes* SteelCone = new G4VisAttributes(green);
voFaSteelCone25->SetVisAttributes(SteelCone);




  // Pos 16
  ///////////////////////////////////
  //    FA Steel Plate 310 mm      //
  //    Drawing ALIP2A__00xx       //
  ///////////////////////////////////

G4double dzFaSteelCone31 = 31*cm/2;
G4double rInFaSteelCone31A = rOuFaWTube4C1;
G4double rOuFaSteelCone31A = rOuFaSteelCone25B;
G4double rInFaSteelCone31B = rOuFaWTube4C2;
G4double rOuFaSteelCone31B = rOuFaSteelCone31A + 2*dzFaSteelCone31 * angle10;


G4Cons* voFaSteelCone31_cone = new G4Cons("voFaSteelCone31_cone", (rInFaSteelCone31A + eps), (rOuFaSteelCone31A - eps), 
		                    (rInFaSteelCone31B + eps), (rOuFaSteelCone31B - eps), dzFaSteelCone31,
				    angle0, angle360);


G4LogicalVolume* voFaSteelCone31 = new G4LogicalVolume(voFaSteelCone31_cone, kMedSteel, "voFaSteelCone31");

Logical_volumes.push_back(voFaSteelCone31);

//setting color
voFaSteelCone31->SetVisAttributes(SteelCone);


 ///////////////////////////////////
  //    FA Composite Ring          //
  //    Drawing ALIP2A__0126       //
  ///////////////////////////////////
  // 1st section
  G4double dzFaCompRing1 = 0.8*cm/2;
  G4double rInFaCompRing1 = 11.0*cm / 2.;
  G4double rOuFaCompRing1 = 32.4*cm / 2.;
  // 2nd section
  G4double dzFaCompRing2 = 1.2*cm/2;
  G4double rInFaCompRing2 = 14.0*cm / 2.;
  G4double rOuFaCompRing2 = 35.3*cm / 2.;



  
  G4Cons* shFaCompRing1_cone = new G4Cons("shFaCompRing1_cone", rInFaCompRing1,  rOuFaCompRing1,rInFaCompRing1, rOuFaCompRing1, dzFaCompRing1, angle0, angle360);

  G4LogicalVolume* shFaCompRing1 = new G4LogicalVolume(shFaCompRing1_cone, kMedCSh, "shFaCompRing1");
  

  
  G4Cons* shFaCompRing2_cone = new G4Cons("shFaCompRing2_cone", rInFaCompRing2,  rOuFaCompRing2,rInFaCompRing2, rOuFaCompRing2, dzFaCompRing2, angle0, angle360);

  G4LogicalVolume* shFaCompRing2 = new G4LogicalVolume(shFaCompRing2_cone, kMedCSh, "shFaCompRing2");


  G4AssemblyVolume* voFaCompRing = new G4AssemblyVolume();

  // 1st section
  z = dzFaCompRing1;
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaCompRing->AddPlacedVolume(shFaCompRing1, Tr);
  
  // 2nd section
  z += dzFaCompRing1 + dzFaCompRing2;
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaCompRing->AddPlacedVolume(shFaCompRing2, Tr);
  

 //setting color
shFaCompRing1->SetVisAttributes(SteelCone);
shFaCompRing2->SetVisAttributes(SteelCone); 


///////////////////////////////////
  //    FA Magnesium Ring          //
  //    Drawing ALIP2A__0127       //
  ///////////////////////////////////
  //
  // The inner radii
  // section 1+3
  G4double dzFaMgRingO = 0.7*cm/2;
  G4double rInFaMgRingO = 3.0*cm;
  // section 2
  G4double dzFaMgRingI = 0.6*cm/2;
  G4double rInFaMgRingI = 3.5*cm;




  G4Cons* shFaMgRing1_cone = new G4Cons("shFaMgRing1_cone", rInFaMgRingO,  rInFaCompRing1,rInFaMgRingO, rInFaCompRing1, dzFaMgRingO, angle0, angle360);

  G4LogicalVolume* shFaMgRing1 = new G4LogicalVolume(shFaMgRing1_cone, kMedMg, "shFaMgRing1");


  G4Cons* shFaMgRing2_cone = new G4Cons("shFaMgRing2_cone", rInFaMgRingI,  rInFaCompRing1,rInFaMgRingI, rInFaCompRing1, dzFaMgRingI / 2., angle0, angle360);

  G4LogicalVolume* shFaMgRing2 = new G4LogicalVolume(shFaMgRing2_cone, kMedMg, "shFaMgRing2");


  G4Cons* shFaMgRing3_cone = new G4Cons("shFaMgRing2_cone", rInFaMgRingI,  rInFaCompRing2,rInFaMgRingI, rInFaCompRing2, dzFaMgRingI / 2., angle0, angle360);

  G4LogicalVolume* shFaMgRing3 = new G4LogicalVolume(shFaMgRing3_cone, kMedMg, "shFaMgRing3");

 
  G4Cons* shFaMgRing4_cone = new G4Cons("shFaMgRing4_cone", rInFaMgRingO, rInFaCompRing2,rInFaMgRingO, rInFaCompRing2, dzFaMgRingO, angle0, angle360);

  G4LogicalVolume* shFaMgRing4 = new G4LogicalVolume(shFaMgRing4_cone, kMedMg, "shFaMgRing4");



  G4AssemblyVolume* voFaMgRing = new G4AssemblyVolume();

  // 1st section
  z = dzFaMgRingO;
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaMgRing->AddPlacedVolume(shFaMgRing1, Tr);
  
  // 2nd section
  z += dzFaMgRingO + (dzFaMgRingI / 2.);
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaMgRing->AddPlacedVolume(shFaMgRing2, Tr);
  
  // 3rd section

  z += (dzFaMgRingI / 2.) + (dzFaMgRingI / 2.);
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaMgRing->AddPlacedVolume(shFaMgRing3, Tr);  
  
  // 4th section

  z += (dzFaMgRingI / 2.) + dzFaMgRingO;
  Ta.setZ(z);
  Tr = G4Transform3D(Ra,Ta);
  voFaMgRing->AddPlacedVolume(shFaMgRing4, Tr);  

  

//setting color
shFaMgRing1->SetVisAttributes(copperVisAttributes);
shFaMgRing2->SetVisAttributes(copperVisAttributes);
shFaMgRing3->SetVisAttributes(copperVisAttributes);
shFaMgRing4->SetVisAttributes(copperVisAttributes);


///---------MAKING THE ASSEMBLY VOLUMES----------------------///


// voFaGraphiteCone, voFaConcreteCone, voFaSteelCone25, voFaSteelCone31


G4AssemblyVolume* voFaAccM = new G4AssemblyVolume();

z = dzFaGraphiteConeS;
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaAccM->AddPlacedVolume(shFaGraphiteCone1, Tr); 

z = z + dzFaGraphiteConeS + (dzFaGraphiteCone - dzFaGraphiteConeS);
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaAccM->AddPlacedVolume(shFaGraphiteCone2, Tr);

z = z + (dzFaGraphiteCone - dzFaGraphiteConeS) + dzFaConcreteCone;
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaAccM->AddPlacedVolume(voFaConcreteCone, Tr);
 
z = z+ dzFaConcreteCone + dzFaSteelCone25;
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaAccM->AddPlacedVolume(voFaSteelCone25, Tr);

z = z + dzFaSteelCone25 + dzFaSteelCone31;
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaAccM->AddPlacedVolume(voFaSteelCone31, Tr);


//inner shield

G4AssemblyVolume* voFaInnerShield = new G4AssemblyVolume();

z = dzFaWTube1C;
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaInnerShield->AddPlacedVolume(shFaWTube1_1, Tr);

z += dzFaWTube1C + dzFaWTube1R;
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaInnerShield->AddPlacedVolume(shFaWTube1_2, Tr);

z = 2*dzFaWTube1 + dzFaWTube2F - 0.6*cm;
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaInnerShield->AddPlacedVolume(shFaWTube2_1, Tr);

z+= dzFaWTube2F + dzFaWTube2C; 
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaInnerShield->AddPlacedVolume(shFaWTube2_2, Tr);

z+= dzFaWTube2C + dzFaWTube3; 
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaInnerShield->AddPlacedVolume(voFaWTube3, Tr);

z+= dzFaWTube3 + dzFaWTube4; 
Ta.setZ(z);
Tr = G4Transform3D(Ra,Ta);
voFaInnerShield->AddPlacedVolume(voFaWTube4, Tr);


//Position of volumes in world

G4double voFaAccM_z = zFa + 2*dzFaFlange - mag_position; 

G4double voFaInnerShield_z = zFa + 2*dzFaGraphiteConeS + 2*dzFaFlange - mag_position;


G4double voFaWTube5_z = zFa + 412.*cm - mag_position;





G4double voFaWPlateteA_z = zFa - mag_position;


G4double voFaSteelEnvelope_z = zFa + 2*dzFaWPlate - mag_position;


//G4double shFaSteelEnvelopeT = ... (Onde...)


G4double voFaPbCone_z = voFaSteelEnvelope_z + 2*dzSteelEnvelopeFC; 


G4double voFaCH2Cone_z = voFaPbCone_z + 2*dzFaPbCone + dzFaCH2Cone ;

G4double voFaFlange_z = zFa - mag_position;

G4double dzFa = (dzFaFlange + dzFaGraphiteCone + dzFaConcreteCone + dzFaSteelCone25 + dzFaSteelCone31)*2;

G4double voFaEndPlate_z = zFa + (dzFa + dzEndPlate) - mag_position;

G4double voFaCompRing_z = zFa - mag_position;

G4double voFaMgRing_z = zFa - mag_position;


//--------------STRUCTURE SECTION ----------------------------------
 
build_abs = true;
//build_abs = false;

build_magnetic = true;
//build_magnetic = false;

if (build_abs) {

Ta.setZ(voFaAccM_z);
Tr = G4Transform3D(Ra,Ta);
voFaAccM->MakeImprint(logicMag, Tr);


Ta.setZ(voFaInnerShield_z);
Tr = G4Transform3D(Ra,Ta);
voFaInnerShield->MakeImprint(logicMag, Tr);

Ta.setZ(voFaSteelEnvelope_z);
Tr = G4Transform3D(Ra,Ta);
voFaSteelEnvelope->MakeImprint(logicMag, Tr);

Ta.setZ(voFaWPlateteA_z);
Tr = G4Transform3D(Ra,Ta);
voFaWPlateteA->MakeImprint(logicMag, Tr);

Ta.setZ(voFaFlange_z);
Tr = G4Transform3D(Ra,Ta);
voFaFlange->MakeImprint(logicMag, Tr);

Ta.setZ(voFaPbCone_z);
Tr = G4Transform3D(Ra,Ta);
voFaPbCone->MakeImprint(logicMag, Tr);


Ta.setZ(voFaCompRing_z);
Tr = G4Transform3D(Ra,Ta);
voFaCompRing->MakeImprint(logicMag, Tr);


Ta.setZ(voFaMgRing_z);
Tr = G4Transform3D(Ra,Ta);
voFaMgRing->MakeImprint(logicMag, Tr);


new G4PVPlacement(0,
		 G4ThreeVector(0, 0,voFaCH2Cone_z),
		 voFaCH2Cone,
		 "shFaGraphiteCone1",
		 logicMag,
		 false,
		 1, checkOverlaps);



new G4PVPlacement(0,
		 G4ThreeVector(0, 0,voFaEndPlate_z),
		 voFaEndPlate,
		 "voFaEndPlate",
		 logicMag,
		 false,
		 1, checkOverlaps);



new G4PVPlacement(0,
		 G4ThreeVector(0, 0,voFaWTube5_z),
		 voFaWTube5,
		 "voFaWTube5",
		 logicMag,
		 false,
		 1, checkOverlaps);






/*
new G4PVPlacement(0,
		 G4ThreeVector(0, 0,shFaGraphiteCone1_z),
		 shFaGraphiteCone1,
		 "shFaGraphiteCone1",
		 logicMag,
		 false,
		 1, checkOverlaps);

new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaGraphiteCone2_z),
		    shFaGraphiteCone2,
		    "shFaGraphiteCone2",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );




new G4PVPlacement(0,
		    G4ThreeVector(0,0,voFaConcreteCone_z),
		    voFaConcreteCone,
		    "voFaConcreteCone",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,voFaSteelCone25_z),
		    voFaSteelCone25,
		    "voFaSteelCone25_z",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,voFaSteelCone31_z),
		    voFaSteelCone31,
		    "voFaSteelCone31",
		    logicMag,
		    false,
		    1
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaWTube1_1_z),
		    shFaWTube1_1,
		    "shFaWTube1_1",
		    logicMag,
		    false,
		    1,checkOverlaps
		 );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaWTube1_2_z),
		     shFaWTube1_2,
		    " shFaWTube1_2",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );




new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaWTube2_1_z),
		    shFaWTube2_1,
		    "shFaWTube2_1",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaWTube2_2_z),
		    shFaWTube2_2,
		    "shFaWTube2_2",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,voFaWTube3_z),
		    voFaWTube3,
		    "voFaWTube3",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,voFaWTube4_z),
		    voFaWTube4,
		    "voFaWTube4",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaWPlateA1_z),
		    shFaWPlateA1,
		    "shFaWPlateA1",
		    logicMag,
		    false,
		    1, checkOverlaps);



new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaWPlateA2_z),			
		    shFaWPlateA2,
		    "shFaWPlateA2",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaWPlateA3_z),
		    shFaWPlateA3,
		    "shFaWPlateA3",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaWPlateA4_z),
		    shFaWPlateA4,
		    "shFaWPlateA4",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaSteelEnvelopeC1_z),
		    shFaSteelEnvelopeC1,
		    "shFaSteelEnvelopeC1",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );




new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaSteelEnvelopeC2_z),
		    shFaSteelEnvelopeC2,
		    "shFaSteelEnvelopeC2",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaSteelEnvelopeC3_z),
		    shFaSteelEnvelopeC3,
		    "shFaSteelEnvelopeC3",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaSteelEnvelopeC4_z),
		    shFaSteelEnvelopeC4,
		    "shFaSteelEnvelopeC4",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaPbCone1_z),
		    shFaPbCone1,
		    "shFaPbCone1",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,shFaPbCone2_z),
		    shFaPbCone2,
		    "shFaPbCone2",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );



/*
new G4PVPlacement(0,
		    G4ThreeVector(0,0,plateAB1_z),
		    plateAB1_Lvolume,
		    "plateAB1",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,plateAB2_z),
		    plateAB2_Lvolume,
		    "plateAB2",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,plateAB3_z),
		    plateAB3_Lvolume,
		    "plateAB3",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );




new G4PVPlacement(0,
		    G4ThreeVector(0,0,plateAB4_z),
		    plateAB4_Lvolume,
		    "plateAB4",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );







new G4PVPlacement(0,
		    G4ThreeVector(0,0,endplate_z),
		    endplate_Lvolume,
		    "endplate1",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,0),
		    endplate2_Lvolume,
		    "endplate2",
		    endplate_Lvolume,
		    false,
		    1,checkOverlaps
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,0),
		    endplate3_Lvolume,
		    "endplate3",
		    endplate_Lvolume,
		    false,
		    1,checkOverlaps
		    );


new G4PVPlacement(0,
		    G4ThreeVector(0,0,wtail1_z),
		    wtail1_Lvolume,
		    "wtail1",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,wtail2_z),
		    wtail2_Lvolume,
		    "wtail2",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,wtail3_z),
		    wtail3_Lvolume,
		    "wtail3",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );



new G4PVPlacement(0,
		    G4ThreeVector(0,0,wtail4_z),
		    wtail4_Lvolume,
		    "wtail4",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );

new G4PVPlacement(0,
		    G4ThreeVector(0,0,wtail5_z),
		    wtail5_Lvolume,
		    "wtail5",
		    logicMag,
		    false,
		    1,checkOverlaps
		    );

	


*/

}








//detector 1
pos_after_detec = 0.01*m;
G4double detec_length = 0.5*cm;
G4double initial_radius = 5*cm;
G4double final_radius = (5.979*m)*tan(9*PI/180.00);

G4Tubs* detec_tub = new G4Tubs("detec_tubs", initial_radius, final_radius, detec_length, 0.*deg,360.*deg);

G4LogicalVolume* detec_volume = new G4LogicalVolume(detec_tub, world_mat, "detec");

G4double detec_z = zFa + (dzFa + 2*dzEndPlate) + 1*cm - mag_position;


 new G4PVPlacement(0,
		 G4ThreeVector(0,0,detec_z),
		 detec_volume,
		 "detector",
		 logicMag,
		 false,
		 0,
		checkOverlaps
		 );




//detector 2

G4Tubs* detec_tub2 = new G4Tubs("detec_tubs", initial_radius, 3*final_radius, detec_length, 0.*deg,360.*deg);
G4LogicalVolume* detec_volume2 = new G4LogicalVolume(detec_tub2, world_mat, "detec2");
//G4double detec_z2 = (z_0  + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + 2*tungs3_dz + 2*tungs4_pDz + 2*dzFaWTail1 + 2*dzFaWTail2 + 2*dzFaWTailR + 2*dzFaWTailB + 2*tail_z + pos_after_detec + 1*m + 3*detec_length) - mag_position;

/*

new G4PVPlacement(0,
		 G4ThreeVector(0,0,detec_z2),
		 detec_volume2,
		 "detector2",
		 logicMag,
		 false,
		 0,
		false
		 );

*/

num_detec=2;



  //fScoringVolume1 = carbon_Lvolume;
  //fScoringVolume2 = concrete_Lvolume;



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
    if (build_magnetic) {
        logicMag->SetFieldManager(FieldMgr, true);
        G4cout << "Applying magnetic field" << G4endl;
     } else {
        G4cout << "magnetic field disabled" << G4endl;
    }

  


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......*dzFaWPlate
