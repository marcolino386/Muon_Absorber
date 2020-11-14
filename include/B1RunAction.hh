
/// \file B1RunAction.hh
/// \brief Definition of the B1RunAction class

#ifndef B1RunAction_h
#define B1RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"
#include <G4Timer.hh>
#include <vector>

class G4Run;




using namespace std;
/// Run action class
class B1RunAction : public G4UserRunAction
{
  public:
    B1RunAction();
    virtual ~B1RunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void add_number_of_event(G4int detec_id);
    /// add nEvent value in vector
    void add_event(){nEvent++;}
    /// increase 1 in nEvent variable
    G4int get_n_event(G4int detec_id) {
        ///get the number of the event to store in the output file
        return num_event_detec[detec_id];	
       // G4cout << num << G4endl;
       }
    G4Timer* timer;
  private:
    
    G4int nEvent;
    G4double n_of_mu_plus;
    G4double  n_of_mu_minus;
    std::vector<G4int >num_event_detec;
};

#endif

