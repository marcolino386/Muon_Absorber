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
  G4double env_sizeXY = 10*m, env_sizeZ = 10*m;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;



  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  G4Box* solidWorld =
    new G4Box("World",                       //its name
      5*m, 5*m, 9*m);     //its size

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


G4double pos_after_detec = 3.*m;
G4double z_0 = 0.9*m;
G4double carbon_pDz = 1.125*m;
G4double concrete_pDz = 0.76*m;
//Area with magnetic field

G4Box* solidMag =
    new G4Box("World",                       //its name
      3*m, 3*m, (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec +0.25*m)/2);

    logicMag =
    new G4LogicalVolume(solidMag,          //its solid
                        world_mat,           //its material
                        "MagField_box"); 

G4VPhysicalVolume* physMag =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0, (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec +0.25*m)/2),       //at (0,0,0)
                      logicMag,            //its logical volume
                      "MagField",               //its name
                      logicWorld,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


 //-------- ABSORBER --------------
 //materials
 
  //carbon

  G4Material* carbon = nist->FindOrBuildMaterial("G4_C");
  
  //concrete
  
  G4Material* concrete = nist->FindOrBuildMaterial("G4_CONCRETE");

  //initial values from the aborber paper

  G4double initial_angle = 2.;
  G4double final_angle = 10.;
  z_0 = 0.9*m;
# define PI 3.14159265

  //cone trunks

  //carbon cone
  carbon_pDz = 1.125*m;
  G4double carbon_pRmin1 = z_0*tan(initial_angle*PI/180.00);
  G4double carbon_pRmax1 = z_0*tan(final_angle*PI/180.00);
  G4double carbon_pRmin2 = (z_0 + 2*carbon_pDz)*tan(initial_angle*PI/180.00);
  G4double carbon_pRmax2 = (z_0 + 2*carbon_pDz)*tan(final_angle*PI/180.00);
  G4double carbon_pSphi = 0.*deg;
  G4double carbon_pDphi = 360.*deg;

  G4Cons* carbon_cons = new G4Cons("carbon_cons", carbon_pRmin1, carbon_pRmax1, 
		                   carbon_pRmin2, carbon_pRmax2, carbon_pDz,
				    carbon_pSphi, carbon_pDphi);

  G4LogicalVolume* carbon_Lvolume = new G4LogicalVolume(carbon_cons, carbon, "carbon_logical");


 G4double carbon_z = z_0 + carbon_pDz - (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec +0.25*m)/2 ;

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,carbon_z),
		    carbon_Lvolume,
		    "carbon_cone",
		    logicMag,
		    false,
		    1
		    );


 G4VisAttributes* blue = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
 carbon_Lvolume->SetVisAttributes(blue);


 //concrete cone
 
 concrete_pDz = 0.76*m; // half lenght of cone
  G4double concrete_pRmin1 = carbon_pRmin2;  // inner radius at the begining of cone /
  G4double concrete_pRmax1 = carbon_pRmax2; // outer radius at the begining of cone
  G4double concrete_pRmin2 = (z_0 + 2*carbon_pDz + 2*concrete_pDz)*tan(initial_angle*PI/180.00); // inner radius at the end of cone
  G4double concrete_pRmax2 = (z_0 + 2*carbon_pDz + 2*concrete_pDz)*tan(final_angle*PI/180.00); // outer radius at the end of cone
    G4double concrete_pSphi = 0.*deg; //starting angle
  G4double concrete_pDphi = 360.*deg; // eng angle
  
 G4Cons* concrete_cons = new G4Cons("concrete_cons", concrete_pRmin1, concrete_pRmax1, concrete_pRmin2, concrete_pRmax2, concrete_pDz, concrete_pSphi, concrete_pDphi);

 G4LogicalVolume* concrete_Lvolume = new G4LogicalVolume(concrete_cons, concrete, "concrete_logical");

G4double concrete_z = (z_0 + 2*carbon_pDz + concrete_pDz) - (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec +0.25*m)/2;

 new G4PVPlacement(0,
		 G4ThreeVector(0, 0, concrete_z),
		 concrete_Lvolume,
		 "concrete_cone",
		 logicMag,
		 false,
		 1);

G4VisAttributes* red = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
concrete_Lvolume->SetVisAttributes(red);
  

//detector
pos_after_detec = 3.*m;
G4double detec_length = 0.5*cm;
G4double initial_radius = (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec)*tan(initial_angle*PI/180.00);
G4double final_radius = (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec)*tan(final_angle*PI/180.00);

G4Tubs* detec_tub = new G4Tubs("detec_tubs", initial_radius, final_radius, detec_length, 0.*deg,360.*deg);

G4LogicalVolume* detec_volume = new G4LogicalVolume(detec_tub, world_mat, "detec");

G4double detec_z = (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec) - (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec +0.25*m)/2;

 new G4PVPlacement(0,
		 G4ThreeVector(0,0,detec_z),
		 detec_volume,
		 "detector",
		 logicMag,
		 false,
		 0,
		false
		 );


G4VisAttributes* color = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
 carbon_Lvolume->SetVisAttributes(color );

  fScoringVolume1 = carbon_Lvolume;
  fScoringVolume2 = concrete_Lvolume;

//sensitive detector

auto sdman = G4SDManager::GetSDMpointer();
G4String SDname = "SD";
auto sensitive = new B1SD(SDname);
sdman->AddNewDetector(sensitive);
detec_volume->SetSensitiveDetector(sensitive);

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
