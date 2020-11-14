
/// \file B1DetectorConstruction.hh
/// \brief Definition of the B1DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include 	<vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class B1SD;
class G4VSensitiveDetector;


using namespace std;

///brief Class that build the absorber Geometry. It is a structure composed of Carbon, Concrete, Lead, Magnesium and Polyethilene. 

class B1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    B1DetectorConstruction();
    virtual ~B1DetectorConstruction();
 
    virtual void ConstructSDandField();
    /// Add magnetic field to a specific logical volume.
  
    virtual G4VPhysicalVolume* Construct();
    
    G4LogicalVolume* GetScoringVolume1() const { return fScoringVolume1; }
    G4LogicalVolume* GetScoringVolume2() const { return fScoringVolume2; }

    G4bool get_sim_state() const {return build_abs;}
    /** Returns the boolean value to know if the absorber should be constructed or not */
 
    std::vector<G4LogicalVolume* > GetVolumes() const {return Logical_volumes;}
    /** Put all Logical volume in a vector */    

    G4int getNumDetec() const {return num_detec;}
   /** Get the number of detectors */

  protected:
    G4LogicalVolume*  fScoringVolume1;
    G4bool build_abs;
    G4int num_detec;
    G4LogicalVolume* fScoringVolume2;
    std::vector<G4LogicalVolume *> Logical_volumes;
    G4LogicalVolume* logicMag;
    G4LogicalVolume* logicWorld;
    G4bool build_magnetic;


};
    

#endif
