
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
  fRunAction(runAction)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event*)
{    
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::EndOfEventAction(const G4Event* event)
{   
  
    fRunAction->add_event();

     double n_mu_p = 0.0;
     double n_mu_m = 0.0;
    

    const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
 

   const G4ParticleGun* particleGun = generatorAction->GetParticleGun();

   
   G4double part_en = (particleGun->GetParticleEnergy())/GeV;
   
  
   

  G4double x = particleGun->GetParticleMomentumDirection().x();
  G4double y = particleGun->GetParticleMomentumDirection().y();
  G4double z = particleGun->GetParticleMomentumDirection().z();

  G4double P = particleGun->GetParticleMomentum()/GeV;
  

   #define PI 3.14159265
    
   // G4double angle_a = atan(y/z)*180/PI;
    
   G4double angle = acos((z/sqrt(x*x+y*y+z*z)));
   
   G4double angle_a = acos((z/sqrt(x*x+y*y+z*z)))*180/PI;
    
   G4double Pz = P*cos(angle);
  
     const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
   
     G4bool sim_struct = detectorConstruction->get_sim_state();
 
   G4SDManager * SDman = G4SDManager::GetSDMpointer();
   G4HCofThisEvent* HCE = event->GetHCofThisEvent();
 

   std::vector<G4int> col;
   std::vector<B1HitsCollection*> HitsCol;

   
      
     G4int num = 4;
     col.reserve(num + 1);
     HitsCol.reserve(num + 1);



   for(G4int i=0; i < 2; i++) {
     
       for(G4int j=1; j <= 2; j++){

       std::stringstream energy_k;
       std::stringstream angle;
       std::stringstream momentum_z;
       angle << std::setprecision(4) << angle_a;
       energy_k << std::setprecision(4) << part_en;
       momentum_z << std::setprecision(4) << Pz;
         
       std::ofstream mu_p_pos("plus" + std::to_string(i + 1) + std::to_string(j) +  "/Energy" + momentum_z.str() + "_Angle_" + angle.str() + ".dat",std::ios_base::app);
   std::ofstream mu_m_pos("minus" + std::to_string(i + 1) +  "/Energy" + momentum_z.str() + "_Angle_" + angle.str() + ".dat",std::ios_base::app);
   
       col[i] = SDman->GetCollectionID("SD" + std::to_string(i + 1) + std::to_string(j));

       if(HCE) {
          HitsCol[i] = (B1HitsCollection*)(HCE->GetHC(col[i]));
          }

      
          int n_hit = HitsCol[i]->entries();
           //G4cout << n_hit << G4endl;
     
     
          //G4cout << "My detector has " << n_hit << "hits" << G4endl;
          B1Hits* hit = new B1Hits;

          for(int i1 = 0; i1 < n_hit; i1++) {
           B1Hits* hit = (*HitsCol[i])[i1];
           G4double energy = hit->getParticleEnergy();
           G4ThreeVector position = hit->getParticlePos();
           G4ThreeVector momentum = hit->getParticleMomentum();
           const G4String name = hit->getParticleInTarget();
           if (name == "mu+") {
              mu_p_pos << position.x()/(m) << "  " << position.y()/(m) << " "  << energy/GeV << " " << momentum.x()/GeV << " " << momentum.y()/(GeV)<< " " <<  momentum.z()/(GeV) << " " << energy_k.str() << " " << angle.str() <<"\n";
             } else if (name == "mu-") {
               mu_m_pos << position.x()/(m) << "  " << position.y()/(m) << " "  << energy/GeV << " " << momentum.x()/GeV << " " << momentum.y()/(GeV)<< " " <<  momentum.z()/(GeV) << " " << energy_k.str() << " " << angle.str() <<"\n";
         
             }
            
           
           } 



  }
}

  /*  

 for(G4int i=0; i < num; i++) { 
   std::stringstream energy_k;
   std::stringstream angle;
   std::stringstream momentum_z;
   angle << std::setprecision(4) << angle_a;
   energy_k << std::setprecision(4) << part_en;
   momentum_z << std::setprecision(4) << Pz;


  
   std::ofstream mu_p_pos("plus" + std::to_string(i + 1) +  "/Energy" + momentum_z.str() + "_Angle_" + angle.str() + ".dat",std::ios_base::app);
   std::ofstream mu_m_pos("minus" + std::to_string(i + 1) +  "/Energy" + momentum_z.str() + "_Angle_" + angle.str() + ".dat",std::ios_base::app);
   

   col[i] = SDman->GetCollectionID("SD" + std::to_string(i + 1));

   if(HCE) {
     HitsCol[i] = (B1HitsCollection*)(HCE->GetHC(col[i]));
   }

      
     int n_hit = HitsCol[i]->entries();
     //G4cout << n_hit << G4endl;
     
     
     //G4cout << "My detector has " << n_hit << "hits" << G4endl;
     B1Hits* hit = new B1Hits;

 for(int i1 = 0; i1 < n_hit; i1++) {
      B1Hits* hit = (*HitsCol[i])[i1];
      
      const G4String name = hit->getParticleInTarget();
      
      //G4cout << name << G4endl;
      if (name == "mu+" || name=="mu-") {
	 G4double energy = hit->getParticleEnergy();
         G4ThreeVector position = hit->getParticlePos();
         G4ThreeVector momentum = hit->getParticleMomentum();
	if (name == "mu+") {
	
          if (sim_struct) {
         
               //store position
      		fRunAction->add_number_of_event(i);
                G4int numb_of_event = fRunAction->get_n_event(i);
      		mu_p_pos << numb_of_event << " " << position.x()/(m) << "  " << position.y()/(m) << " "  << energy/GeV << " " << momentum.x()/GeV << " " << momentum.y()/(GeV)<< " " <<  momentum.z()/(GeV) << " " << energy_k.str() << " " << angle.str() <<"\n";
      		           
               

         } else {
		//store position
      		fRunAction->add_number_of_event(i);
                G4int numb_of_event = fRunAction->get_n_event(i);
                if (numb_of_event < 2){
		    mu_p_pos << 0 << " " << position.x()/(m) << "  " << position.y()/(m) << " " << energy/GeV << " " << momentum.x()/GeV << " " << momentum.y()/(GeV)<< " " <<  momentum.z()/(GeV)<< " " << energy_k.str() << " " << angle.str() <<"\n";
		}

	   }

          }  else if (name=="mu-") {
        
        
         if (sim_struct) {

		//store position
      	
                fRunAction->add_number_of_event(i);
                G4int numb_of_event = fRunAction->get_n_event(i);
      		mu_m_pos << numb_of_event << " " << position.x()/(m) << "  " << position.y()/(m) << " " << energy/GeV << " " << momentum.x()/GeV << " " << momentum.y()/(GeV)<< " " <<  momentum.z()/(GeV) << " " << energy_k.str() << " " << angle.str() <<"\n";
		

          } else {

                //store position
     		fRunAction->add_number_of_event(i);
                G4int numb_of_event = fRunAction->get_n_event(i);
                if(numb_of_event < 2) {
		   mu_m_pos << 0 << " " <<position.x()/(m) << "  " << position.y()/(m) << " " << energy/GeV << " " << momentum.x()/GeV << " " << momentum.y()/(GeV)<< " " <<  momentum.z()/(GeV) << " " << energy_k.str() << " " << angle.str() <<"\n";
		}

      }
	

	}
    }
    
}


mu_m_pos.close();
mu_p_pos.close();

}


*/

}
   


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
