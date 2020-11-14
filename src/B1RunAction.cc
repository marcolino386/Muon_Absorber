
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
#include "G4Timer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction(),
  n_of_mu_plus(0.),
  n_of_mu_minus(0.),
  nEvent(0)
{ 
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
 G4cout << "BEGIN RUN" << G4endl;
	
  n_of_mu_plus = 0.;
  n_of_mu_minus = 0.;
  

  const B1DetectorConstruction* detectorConstruction
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

   
   
   for (G4int i=0; i<detectorConstruction->getNumDetec();i++) {

       num_event_detec.push_back(0);
}

  
  
  timer = new G4Timer();
  timer->Start();
     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
   
  timer->Stop();
 
 
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
  
  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance  

  const B1DetectorConstruction* detectorConstruction
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());



  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  
  G4double angle_a;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    G4double x = particleGun->GetParticleMomentumDirection().x();
    G4double y = particleGun->GetParticleMomentumDirection().y();
    G4double z = particleGun->GetParticleMomentumDirection().z();
    
    runCondition += G4BestUnit(particleEnergy,"Energy");
   
   #define PI 3.14159265
    angle_a = acos((z/sqrt(x*x+y*y+z*z)))*180/PI;
   
   G4cout << nEvent << G4endl;
        

  }
 
 

  num_event_detec.clear();
  
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
     << "Time of Run: " << timer->GetRealElapsed()
     << G4endl
     << "Angle: " << angle_a << " graus"
     << G4endl;
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::add_number_of_event(G4int detec_id)
{
  num_event_detec[detec_id] += 1;
 // G4cout << num_event_detec[detec_id]  << G4endl;
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

