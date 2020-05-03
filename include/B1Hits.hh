#include "G4VHit.hh"
#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"

class B1Hits: public G4VHit {

public:
  B1Hits();
  ~B1Hits();

  void print();
  void set_partdef (const G4String particle_name);
  void set_energy (const G4double particle_en);
  void set_position (const G4ThreeVector particle_pos);
  const G4String getParticleInTarget() {
    return fParticleInTarget;
  }
  const G4double getParticleEnergy() {
    return fparticle_energy;
  }
  const G4ThreeVector getParticlePos() {
    return fparticle_position;;
  }
  
  
private:

 G4String fParticleInTarget;
 G4double fparticle_energy;
 G4ThreeVector fparticle_position;

};

typedef G4THitsCollection<B1Hits> B1HitsCollection;

