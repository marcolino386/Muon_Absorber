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
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1RunAction.hh"
#include "B1Hits.hh"
#include "B1DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep_mu_plus(0.),
  fEdep_mu_minus(0.)

{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep_mu_plus = 0.;
  fEdep_mu_minus = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event* event)
{   
  // accumulate statistics in run action

  if(fEdep_mu_plus != 0.) {
     fRunAction->AddEdep1(fEdep_mu_plus);
     fRunAction->AddMu_plus();
} else if (fEdep_mu_minus != 0.) {
     fRunAction->AddEdep2(fEdep_mu_minus);
     fRunAction->AddMu_minus();

}
  

 static int CHCID = -1;
if (CHCID < 0) {
    CHCID = G4SDManager::GetSDMpointer()->GetCollectionID("SD");//Descobrir isso
  }



 // pega as collections ID's
 
 G4SDManager * SDman = G4SDManager::GetSDMpointer();
 G4HCofThisEvent* HCE = event->GetHCofThisEvent();
 
 B1HitsCollection* HitsCol = 0;


  if(HCE) {
    HitsCol = (B1HitsCollection*)(HCE->GetHC(CHCID));
  }

 
    int n_hit = HitsCol->entries();
     //G4cout << "My detector has " << n_hit << "hits" << G4endl;
     B1Hits* hit = new B1Hits;
    

     double n_mu_p = 0.0;
     double n_mu_m = 0.0;
     G4double total_energy_mu_p = 0;
     G4double total_energy_mu_m = 0;

    const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
 

   const G4ParticleGun* particleGun = generatorAction->GetParticleGun();


   G4double particleEnergy = particleGun->GetParticleEnergy();
     
     std::ofstream mu_p_pos("position_mu_p.dat",std::ios_base::app);
 
  std::ofstream mu_m_pos("position_mu_m.dat",std::ios_base::app);

     for(int i1 = 0; i1 < n_hit; i1++) {
      B1Hits* hit = (*HitsCol)[i1];
      const G4String name = hit->getParticleInTarget();
      //G4cout << name << G4endl;
      if (name == "mu+" || name=="mu-") {
	 G4double energy = hit->getParticleEnergy();
         G4ThreeVector position = hit->getParticlePos();
	if (name == "mu+") {
	n_mu_p++;
        mu_p_pos << position.x()/(m) << "  " << position.y()/(m) << "\n";
 	total_energy_mu_p += energy;
	} else if (name=="mu-") {
        total_energy_mu_m += energy;
        mu_m_pos << position.x()/(m) << "  " << position.y()/(m) << "\n";
	n_mu_m++;

	}
    }
      
}

 mu_p_pos.close();
mu_m_pos.close();


 if(n_mu_p != 0.0) {
       G4double value = total_energy_mu_p/n_mu_p;
       //G4cout << value/(GeV) << G4endl;
       fRunAction->AddE_mup(value);
    } else {fRunAction->AddE_mup(total_energy_mu_p);}
 
 if(n_mu_m != 0.0) {
       G4double value = total_energy_mu_m/n_mu_m;
       fRunAction->AddE_mum(value);
    } else {fRunAction->AddE_mum(total_energy_mu_m);}

} 


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
