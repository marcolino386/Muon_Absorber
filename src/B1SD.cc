#include "B1SD.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"


/// \file B1SD.cc
/// \brief Implementation of the B1SD class.


B1SD::B1SD(G4String SDname): G4VSensitiveDetector(SDname),
  hitCollection(nullptr), HCID(-1) {
  //cria a hit collection
  G4cout << "Criando Hit Collection com nome: " << SDname <<G4endl;
  collectionName.insert(SDname);
  track_id = 0;
  sdname = SDname;
}

B1SD::~B1SD() {

}

G4bool B1SD::ProcessHits(G4Step* step, G4TouchableHistory* ROhist) {
  G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();
  //pega o nome da partícula
  const G4String particle_name = step->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
  //pega o Track relacionado à ela
  G4int track = step->GetTrack()->GetTrackID();
//pega energia da particula
 const G4double particle_energy = step->GetPreStepPoint()->GetTotalEnergy();
//pega posição da particula
const G4ThreeVector particle_position = step->GetPreStepPoint()->GetPosition();
//pega momentum da particula
const G4ThreeVector particle_momentum = step->GetPreStepPoint()->GetMomentum();

const G4StepPoint* thePreVL = step->GetPreStepPoint();
const G4VPhysicalVolume* thePostVL = step->GetPostStepPoint()->GetPhysicalVolume();

	
//Checa se a partícula é repetida
if (thePreVL->GetStepStatus() == fGeomBoundary && track == 1) {
      
  
      B1Hits* hit = new B1Hits();
      hit->set_partdef(particle_name);
      hit->set_energy(particle_energy);
      hit->set_position(particle_position);
      hit->set_momentum(particle_momentum);
      //coloca o valor na hitCollecion
      hitCollection->insert(hit);
      track_id = track;
      return true;

}

}

void B1SD::Initialize(G4HCofThisEvent* HCE) {
  hitCollection = new B1HitsCollection(SensitiveDetectorName,collectionName[0]);
  if (HCID < 0) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, hitCollection);
}

void B1SD::EndOfEvent(G4HCofThisEvent* HCE) {
  track_id = 0;
}

