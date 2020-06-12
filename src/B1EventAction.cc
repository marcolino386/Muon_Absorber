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
     //fRunAction->AddMu_plus();
} else if (fEdep_mu_minus != 0.) {
     fRunAction->AddEdep2(fEdep_mu_minus);
    // fRunAction->AddMu_minus();

}
  
    

     double n_mu_p = 0.0;
     double n_mu_m = 0.0;
     G4double total_energy_mu_p = 0;
     G4double total_energy_mu_m = 0;

    const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
 

   const G4ParticleGun* particleGun = generatorAction->GetParticleGun();


   G4double particleEnergy = particleGun->GetParticleEnergy();

  G4double y = particleGun->GetParticleMomentumDirection().y();
  G4double z = particleGun->GetParticleMomentumDirection().z();
  
   #define PI 3.14159265
    G4double angle = atan(y/z)*180/PI;
  
     const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
   
     G4bool sim_struct = detectorConstruction->get_sim_state();

  
 // pega as collections ID's
  
 
   G4SDManager * SDman = G4SDManager::GetSDMpointer();
   G4HCofThisEvent* HCE = event->GetHCofThisEvent();
 

   std::vector<G4int> col;
   std::vector<B1HitsCollection*> HitsCol;

   
      
     G4int num = 2;
     col.reserve(num + 1);
     HitsCol.reserve(num + 1);

    


 for(G4int i=0; i < num; i++) { 
   
   std::ofstream mu_p_pos("data_mu_plus" + std::to_string(i + 1) +  "/Energy" + std::to_string(particleEnergy/GeV) + "_Angle_" + std::to_string(angle) + ".dat",std::ios_base::app);
   std::ofstream mu_m_pos("data_mu_minus" + std::to_string(i + 1) +  "/Energy" + std::to_string(particleEnergy/GeV) + "_Angle_" + std::to_string(angle) + ".dat",std::ios_base::app);

   col[i] = SDman->GetCollectionID("SD" + std::to_string(i + 1));

   if(HCE) {
     HitsCol[i] = (B1HitsCollection*)(HCE->GetHC(col[i]));
   }

      
     int n_hit = HitsCol[i]->entries();
     G4cout << n_hit << G4endl;
     
     
     //G4cout << "My detector has " << n_hit << "hits" << G4endl;
     B1Hits* hit = new B1Hits;

 for(int i1 = 0; i1 < n_hit; i1++) {
      B1Hits* hit = (*HitsCol[i])[i1];
      //G4cout << "aaa" << G4endl;
      const G4String name = hit->getParticleInTarget();
      
      //G4cout << name << G4endl;
      if (name == "mu+" || name=="mu-") {
	 G4double energy = hit->getParticleEnergy();
         G4ThreeVector position = hit->getParticlePos();
         G4ThreeVector momentum = hit->getParticleMomentum();
	if (name == "mu+") {
	
          if (sim_struct) {
		//G4cout << "heyy" << G4endl;
         
               //store position
      		fRunAction->add_number_of_event(i);
                G4int numb_of_event = fRunAction->get_n_event(i);
      		mu_p_pos << numb_of_event << " " << position.x()/(m) << "  " << position.y()/(m) << " "  << energy/GeV << " " << momentum.x()/GeV << " " << momentum.y()/(GeV)<< " " <<  momentum.z()/(GeV) << "\n";
      		
                total_energy_mu_p += energy;
                fRunAction->AddMu_plus();
                n_mu_p++;

         } else {
		//store position
      		fRunAction->add_number_of_event(i);
                G4int numb_of_event = fRunAction->get_n_event(i);
                if (numb_of_event < 2){
		    mu_p_pos << 0 << " " << position.x()/(m) << "  " << position.y()/(m) << " " << particleEnergy/GeV << " " << momentum.x()/GeV << " " << momentum.y()/(GeV)<< " " <<  momentum.z()/(GeV) << "\n";
		}
      		
      		//mu_p_pos0.close();

             }
 	
	} else if (name=="mu-") {
        
        
         if (sim_struct) {

		//store position
      		//std::ofstream mu_m_pos("data_mu_minus/Energy" + std::to_string(particleEnergy/GeV) + "_" + std::to_string(angle) + ".dat",std::ios_base::app);
                fRunAction->add_number_of_event(i);
                G4int numb_of_event = fRunAction->get_n_event(i);
      		mu_m_pos << numb_of_event << " " << position.x()/(m) << "  " << position.y()/(m) << " " << energy/GeV << " " << momentum.x()/GeV << " " << momentum.y()/(GeV)<< " " <<  momentum.z()/(GeV) << "\n";
		//mu_m_pos.close();
 		total_energy_mu_m += energy;
                fRunAction->AddMu_minus();
                 n_mu_m++;

          } else {

                //store position
     		fRunAction->add_number_of_event(i);
                G4int numb_of_event = fRunAction->get_n_event(i);
                if(numb_of_event < 2) {
		   mu_m_pos << 0 << " " <<position.x()/(m) << "  " << position.y()/(m) << " " << particleEnergy/GeV << " " << momentum.x()/GeV << " " << momentum.y()/(GeV)<< " " <<  momentum.z()/(GeV) << "\n";
		}
      		
      		//mu_m_pos.close();

      }
	

	}
    }
      
}


mu_m_pos.close();
mu_p_pos.close();

}
   



 if(n_mu_p != 0.0 && sim_struct == true) {
       G4double value = total_energy_mu_p/n_mu_p;
       //G4cout << value/(GeV) << G4endl;
       fRunAction->AddE_mup(value);
    } else {fRunAction->AddE_mup(total_energy_mu_p);}
 
 if(n_mu_m != 0.0 && sim_struct == true) {
       G4double value = total_energy_mu_m/n_mu_m;
       fRunAction->AddE_mum(value);
    } else {fRunAction->AddE_mum(total_energy_mu_m);}

} 


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
