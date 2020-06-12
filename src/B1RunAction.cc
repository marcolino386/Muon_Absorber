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
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
// #include "B1Run.hh"
#include <stdio.h>      /* printf */
#include <math.h> 
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction(),
  fEdep1(0.),
  fEdep2(0.),
  fE_mup(0.),
  fE_mum(0.),
  n_of_mu_plus(0.),
  n_of_mu_minus(0.)
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep1);
  accumulableManager->RegisterAccumulable(fEdep2);
  accumulableManager->RegisterAccumulable(fE_mup);
  accumulableManager->RegisterAccumulable(fE_mum);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
	
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();
  n_of_mu_plus = 0.;
  n_of_mu_minus = 0.;
  

  const B1DetectorConstruction* detectorConstruction
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

   
   
   for (G4int i=0; i<detectorConstruction->getNumDetec();i++) {

       num_event_detec.push_back(0);
}

  

  //Get angle to write in file
   

    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep1.GetValue() + fEdep2.GetValue();
  G4double edep_mu_plus = fEdep1.GetValue();
   G4double edep_mu_minus = fEdep2.GetValue();
  //G4double edep2 = fEdep2.GetValue();
  G4double mu_p_energy = fE_mup.GetValue();
  G4double mu_m_energy = fE_mum.GetValue();
  //G4double rms = edep2 - edep*edep/nofEvents;
 // if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const B1DetectorConstruction* detectorConstruction
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  //G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/nofEvents;
  G4double mu_plus_edep_av;
  G4double mu_minus_edep_av;
  G4double mu_p_en = mu_p_energy/nofEvents;
  G4double mu_m_en = mu_m_energy/nofEvents;
  //G4double rmsDose = rms/mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    //G4ParticleMomentum direction = particleGun->GetParticleMomentumDirection();
    G4double y = particleGun->GetParticleMomentumDirection().y();
    G4double z = particleGun->GetParticleMomentumDirection().z();
    
    runCondition += G4BestUnit(particleEnergy,"Energy");
   
   #define PI 3.14159265
    G4double angle = atan(y/z)*180/PI;
   


  std::ofstream mu_p_energy("energy_mu_plus_data/" + std::to_string(particleEnergy/GeV) + ".dat",std::ios_base::app);
 
  std::ofstream mu_m_energy("energy_mu_minus_data/" + std::to_string(particleEnergy/GeV) + ".dat",std::ios_base::app);


  if(n_of_mu_plus != 0.) {
     mu_plus_edep_av =  edep_mu_plus/n_of_mu_plus;
   //  mu_p_energy << mu_plus_edep_av/GeV << " " << angle << "\n";

} else {
  mu_plus_edep_av = 0;

} if (n_of_mu_minus !=0.) {
   mu_minus_edep_av =  edep_mu_minus/n_of_mu_minus;
  // mu_m_energy << mu_minus_edep_av/GeV << " " << angle << "\n";

} else {
 mu_minus_edep_av = 0;

}
        
  

  mu_p_energy.close();
  mu_m_energy.close();
  }
 
  //write on files
 
  
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Average Cumulated energy in volumes : " 
     << G4BestUnit(dose,"Energy") << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << "Average energy deposited by mu+ : " 
     << G4BestUnit(mu_plus_edep_av,"Energy") << G4endl
     << G4endl
     << "Average energy deposited by mu- : " 
     << G4BestUnit(mu_minus_edep_av,"Energy") << G4endl
     << G4endl
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::add_number_of_event(G4int detec_id)
{
  num_event_detec[detec_id] += 1;
 // G4cout << num_event_detec[detec_id]  << G4endl;
  
}

void B1RunAction::AddEdep1(G4double edep)
{
  fEdep1  += edep;
  
}

void B1RunAction::AddEdep2(G4double edep)
{
  fEdep2  += edep;
  
}
void B1RunAction::AddE_mup(G4double edep)
{
  fE_mup  += edep;
  
}
void B1RunAction::AddE_mum(G4double edep)
{
  fE_mum  += edep;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

