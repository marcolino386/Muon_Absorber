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
/// \file B1RunAction.hh
/// \brief Definition of the B1RunAction class

#ifndef B1RunAction_h
#define B1RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"
#include <vector>

class G4Run;

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.

using namespace std;

class B1RunAction : public G4UserRunAction
{
  public:
    B1RunAction();
    virtual ~B1RunAction();

    // virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    void AddEdep1 (G4double edep); 
    void AddEdep2 (G4double edep);
    void AddE_mum (G4double edep) ; 
    void AddE_mup (G4double edep) ;
    void AddMu_plus() {n_of_mu_plus++;}
    void AddMu_minus() {n_of_mu_minus++;}
    void add_number_of_event(G4int detec_id);
    void add_event(){nEvent++;}
    G4int get_n_event(G4int detec_id) {
        return num_event_detec[detec_id];	
       // G4cout << num << G4endl;
       }
  private:
    G4Accumulable<G4double> fEdep1;
    G4Accumulable<G4double> fEdep2;
    G4Accumulable<G4double> fE_mum;
    G4Accumulable<G4double> fE_mup;
    G4int nEvent;
    G4double n_of_mu_plus;
    G4double  n_of_mu_minus;
    std::vector<G4int >num_event_detec;
};

#endif

