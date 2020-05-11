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



//some importante initial values

G4double pos_after_detec = 3.*m;
G4double z_0 = 90.0*cm;
G4double carbon_pDz = 225.0*cm/2;
G4double concrete_pDz = 126.*cm/2;
G4double dzFaFlange = 2.*cm;
 G4double dzCarbonConeS = 108.3*cm/2;
G4double mag_position = (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec +0.25*m + dzFaFlange )/2;
//Area with magnetic field

G4Box* solidMag =
    new G4Box("World",                       //its name
      3*m, 3*m, mag_position);

    logicMag =
    new G4LogicalVolume(solidMag,          //its solid
                        world_mat,           //its material
                        "MagField_box"); 

G4VPhysicalVolume* physMag =
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

  G4double A;  // atomic mass
     G4double Z;  // atomic number
    G4double d;  // density
   
  






  //carbon

  G4Material* carbon = nist->FindOrBuildMaterial("G4_C");
  
  //concrete
  
  G4Material* concrete = nist->FindOrBuildMaterial("G4_CONCRETE");
 
  G4Material* elSi = nist->FindOrBuildMaterial("G4_Si");
 
  G4Material* elMn = nist->FindOrBuildMaterial("G4_Mn");
  
  G4Material* elCr = nist->FindOrBuildMaterial("G4_Cr");

  G4Material* elNi = nist->FindOrBuildMaterial("G4_Ni");
  
  G4Material* elFe = nist->FindOrBuildMaterial("G4_Fe");
  
  G4Material* elW = nist->FindOrBuildMaterial("G4_W");

// Stainless steel (Medical Physics, Vol 25, No 10, Oct 1998)



 d = 8.02*g/cm3;

   G4Material* matsteel = new G4Material("Stainless steel",d,5);
   matsteel->AddMaterial(elMn, 0.02);
  matsteel->AddMaterial(elSi, 0.01);
   matsteel->AddMaterial(elCr, 0.19);
   matsteel->AddMaterial(elNi, 0.10);
   matsteel->AddMaterial(elFe, 0.68);

  //initial values from the aborber paper

  G4double initial_angle = 2.;
  G4double final_angle = 10.;
 
# define PI 3.14159265

  //cone trunks

/*

  //carbon cone
 // Inner radius at the front
  Float_t rInFaGraphiteCone1 = 4.5;
  // Outer radius at the front
  Float_t rOuFaGraphiteCone1 = (zFa + dzFaFlange) * angle10;
  // Inner radius at start of inner opening cone
  Float_t rInFaGraphiteCone2 = 7.0;
  // Outer radius at start of inner opening cone
  Float_t rOuFaGraphiteCone2 = (zFa + dzFaFlange + dzFaGraphiteConeS) * angle10;
  // Inner radius the rear
  Float_t rInFaGraphiteCone3 = 11.0;
  // Ouer radius at the rear
  Float_t rOuFaGraphiteCone3 = (zFa + dzFaFlange + dzFaGraphiteCone) * angle10;
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
		                   carbon1_pRmin2, carbon1_pRmax2, dzCarbonConeS,
				    carbon_pSphi, carbon_pDphi);

  G4LogicalVolume* carbon1_Lvolume = new G4LogicalVolume(carbon1_cons, carbon, "carbon1_logical");


 G4double carbon1_z = (z_0 + dzFaFlange + dzCarbonConeS) -  mag_position ;


 


 G4VisAttributes* blue = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
 carbon1_Lvolume->SetVisAttributes(blue);

// inner opening cone

G4double carbon_Dz_2 = (2*carbon_pDz - 2*dzCarbonConeS)/2;
G4double carbon_pRmin1 = carbon1_pRmin2;
G4double carbon_pRmax1 = carbon1_pRmax2;
G4double carbon_pRmin2 = 11.0*cm;
G4double carbon_pRmax2 = (z_0 + dzFaFlange + 2*carbon_pDz)*tan(final_angle*PI/180.00);

G4double carbon_z = (z_0 + dzFaFlange + 2*dzCarbonConeS + carbon_Dz_2) - mag_position;

G4Cons* carbon_cons = new G4Cons("carbon_cons", carbon_pRmin1, carbon_pRmax1, 
		                   carbon_pRmin2, carbon_pRmax2, carbon_Dz_2,
				    carbon_pSphi, carbon_pDphi);



G4LogicalVolume* carbon_Lvolume = new G4LogicalVolume(carbon_cons, carbon, "carbon_logical");

carbon_Lvolume->SetVisAttributes(blue);


 //concrete cone
 

  G4double concrete_pRmin1 = 11.*cm;  // inner radius at the begining of cone /
  G4double concrete_pRmax1 = carbon_pRmax2; // outer radius at the begining of cone



  G4double concrete_pRmin2 = concrete_pRmin1  + (2*concrete_pDz)*tan(initial_angle*PI/180.00); // inner radius at the end of cone
  G4double concrete_pRmax2 = concrete_pRmax1 + (2*concrete_pDz)*tan(final_angle*PI/180.00); // outer radius at the end of cone
    G4double concrete_pSphi = 0.*deg; //starting angle
  G4double concrete_pDphi = 360.*deg; // eng angle
  
 G4Cons* concrete_cons = new G4Cons("concrete_cons", concrete_pRmin1, concrete_pRmax1, concrete_pRmin2, concrete_pRmax2, concrete_pDz, concrete_pSphi, concrete_pDphi);

 G4LogicalVolume* concrete_Lvolume = new G4LogicalVolume(concrete_cons, concrete, "concrete_logical");

G4double concrete_z = (z_0 + 2*carbon_pDz + concrete_pDz + dzFaFlange) - mag_position;


 


G4VisAttributes* red = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0));
concrete_Lvolume->SetVisAttributes(red);


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

G4LogicalVolume* tungs1_Lvolume_1 = new G4LogicalVolume(tungs1_cons_1, elW, "tungs1_logical_1");


G4double tungs1_z_1 = (z_0 + dzFaFlange + 2*dzCarbonConeS + tungs1_pDz_1) - mag_position;

//final

G4double tungs1_pRmin1= tungs1_pRmin1_1;
G4double tungs1_pRmax1 = 15.0*cm / 2.;
G4double tungs1_pRmin2 = tungs1_pRmin1;
G4double tungs1_pRmax2 = tungs1_pRmax1;

G4double tungs1_z = (z_0 + dzFaFlange + 2*dzCarbonConeS + 2*tungs1_pDz_1 + tungs1_pDz_2 ) - mag_position;

G4Cons* tungs1_cons = new G4Cons("Tungsten_shield_part1", tungs1_pRmin1, tungs1_pRmax1, 
		                   tungs1_pRmin2, tungs1_pRmax2, tungs1_pDz_2,
				    carbon_pSphi, carbon_pDphi);


G4LogicalVolume* tungs1_Lvolume = new G4LogicalVolume(tungs1_cons, elW, "tungs1_logical");


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

 G4LogicalVolume* tungs2_Lvolume_1 = new G4LogicalVolume(tungs2_cons_1, elW, "tungs2_logical_1");


G4double tungs2_z_1 = (z_0 + dzFaFlange + 2*dzCarbonConeS + 2*tungs1_pDz + tungs2_pDz_1) - mag_position;

//2


G4double tungs2_pRmax2 = 30.72*cm / 2. - 0.05*cm; // outer radius at the end of cone
G4double tungs2_pRmin2 = 12.58*cm / 2.;
G4double tungs2_pRmin1 = 9.10*cm / 2.;
G4double tungs2_pRmax1 = tungs2_pRmax1_1;

G4double tungs2_z = (z_0 + dzFaFlange + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz_1 + tungs2_pDz_2) - mag_position;

G4Cons* tungs2_cons = new G4Cons("Tungsten_shield_part2_2", tungs2_pRmin1, tungs2_pRmax1, 
		                   tungs2_pRmin2, tungs2_pRmax2, tungs2_pDz_2,
				    carbon_pSphi, carbon_pDphi);

G4LogicalVolume* tungs2_Lvolume = new G4LogicalVolume(tungs2_cons, elW, "tungs2_logical_2");

tungs2_Lvolume->SetVisAttributes(aa);

tungs2_Lvolume_1->SetVisAttributes(aa);

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

G4double tungs3_z = (z_0 + dzFaFlange + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + tungs3_dz) - mag_position;

G4Cons* tungs3_cons = new G4Cons("Tungsten_shield_part3", tungs3_pRmin1, tungs3_pRmax1, 
		                   tungs3_pRmin2, tungs3_pRmax2, tungs3_dz,
				    carbon_pSphi, carbon_pDphi);





G4LogicalVolume* tungs3_Lvolume = new G4LogicalVolume(tungs3_cons, elW, "tungs3_logical");

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

 G4double tungs4_z = (z_0 + dzFaFlange + 2*dzCarbonConeS + 2*tungs1_pDz + 2*tungs2_pDz + 2*tungs3_dz + tungs4_pDz ) - mag_position;

G4Cons* tungs4_cons = new G4Cons("Tungsten_shield_part4", tungs4_pRmin1, tungs4_pRmax1, 
		                   tungs4_pRmin2, tungs4_pRmax2, tungs4_pDz,
				    carbon_pSphi, carbon_pDphi);

G4LogicalVolume* tungs4_Lvolume = new G4LogicalVolume(tungs4_cons, elW, "tungs4_logical");

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

/*

  // 10 deg cone
  Float_t dzFaPbCone10 = 25.9;
  Float_t rInFaPbCone10 = rInFaPbCone5 + dzFaPbCone5 * angle10;
  Float_t rOuFaPbCone10 = 115.2 / 2.;
  // end
  Float_t rInFaPbConeE = 106.05 / 2.;
  Float_t rOuFaPbConeE = 124.35 / 2.;
  // Total length
  Float_t dzFaPbCone = dzFaPbCone5 + dzFaPbCone10;

*/





// Pos 15
  ///////////////////////////////////
  //    FA Steel Plate 250 mm      //
  //    Drawing ALIP2A__00xx       //
  ///////////////////////////////////


G4double SteelCone25_pDz = 25*cm/2;
G4double SteelCone25_Rmin1 = concrete_pRmin2;
G4double SteelCone25_Rmax1 = concrete_pRmax2;
G4double SteelCone25_Rmin2 = SteelCone25_Rmin1 + (2*SteelCone25_pDz)*tan(initial_angle*PI/180.00);
G4double SteelCone25_Rmax2 = SteelCone25_Rmax1 + (2*SteelCone25_pDz)*tan(final_angle*PI/180.00);


G4Cons* Steel25_cons = new G4Cons("Steel25_cons", SteelCone25_Rmin1, SteelCone25_Rmax1, 
		                   SteelCone25_Rmin2, SteelCone25_Rmax2, SteelCone25_pDz,
				    carbon_pSphi, carbon_pDphi);


  
G4LogicalVolume* Steel25_Lvolume = new G4LogicalVolume(Steel25_cons, matsteel, "Steel25_logical");

G4double Steel25_z = (z_0 + 2*carbon_pDz + 2*concrete_pDz + SteelCone25_pDz + dzFaFlange)  - mag_position; 


G4VisAttributes* SteelCone = new G4VisAttributes(G4Colour(1.0, 6.0, 1.0));
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


G4Cons* Steel31_cons = new G4Cons("Steel31_cons", SteelCone31_Rmin1, SteelCone31_Rmax1, 
		                   SteelCone31_Rmin2, SteelCone31_Rmax2, SteelCone31_pDz,
				    carbon_pSphi, carbon_pDphi);


G4LogicalVolume* Steel31_Lvolume = new G4LogicalVolume(Steel31_cons, matsteel, "Steel31_logical");

G4double Steel31_z = (z_0 + 2*carbon_pDz + 2*concrete_pDz + 2*SteelCone25_pDz + SteelCone31_pDz)  - (z_0 + 2*carbon_pDz + 2*concrete_pDz + pos_after_detec +0.25*m)/2; 

Steel31_Lvolume->SetVisAttributes(SteelCone);





//--------------STRUCTURE SECTION ----------------------------------
 

G4bool build_abs = true;

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


}









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
