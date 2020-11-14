
/// \file B1SteppingAction.hh
/// \brief Definition of the B1SteppingAction class

#ifndef B1SteppingAction_h
#define B1SteppingAction_h 

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include 	<vector>
#include "G4TrackStatus.hh"

class B1EventAction;

class G4LogicalVolume;

/// Stepping action class


using namespace std;	
///The stepping Action is used primary to kill the secondary particles to improve the time of run.  
class B1SteppingAction : public G4UserSteppingAction
{
  public:
    B1SteppingAction(B1EventAction* eventAction);
    virtual ~B1SteppingAction();

    /// method from the base class
    virtual void UserSteppingAction(const G4Step*);
    

  private:
    B1EventAction*  fEventAction;
    std::vector<G4LogicalVolume* > Log_volumes;
    G4LogicalVolume* fScoringVolume1;
    G4LogicalVolume* fScoringVolume2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
