\mainpage
Being constructed.
## Brief Introduction

The Muon Absorber project is a Geant4 Simulation which builds the ALICE's Absorber. The main goal of this application is to study the influence of the Absorber on the muon trajectory using a sensitive detector. This project was built based on Geant4's exampleB1.

## Geometry

The geometry is define in the B1DetectorConstruction class. The Absorber is a structure made of Carbon, concrete, lead, magnesium and polyethilene, written based on ALICE's gitlab project in this link: <https://github.com/AliceO2Group/AliceO2/blob/dev/Detectors/Passive/src/Absorber.cxx>

We have also included a Magnetic field of 0.5T, similar to the ALICE's environment.

## Physics list

The Physics List used is the FTFP_BERT, recommended for high energy simulations.

## Primary Generator
   The primary generator is defined in the B1PrimaryGeneratorAction class. The default energy is a mu+ of 1 Gev.

## Sensitive Detector

The sensitive Detector is defined in the B1SD class. It gets the values of momentum, energy, angle and position of the primary muon. It does not store the values of any secondary particle.

##Stepping Action

The Stepping Action is defined in the B1SteppingAction class. It is used to eliminate the secondary particles generated to improve time of run.

##Event Action

The Event Action is defined in the B1EventAction class. In this class the data stored in the hit collection is written in the output files.

##Output file

The application writes the data in 2 .dat file, storing the data in collunms. The index 0 indicates that those are the initial values. The file with "position" is the one that contains the data about positions and the file with "momentum" contains the data about energy, angle and momentum. The relation between collunm and variable is in the "README.txt"



================================================================

Authors: Derós, MAO and Pereira, LG 

Project: <https://github.com/marcolino386/Muon_Absorber>

High Energy Physics Simulations (HEPSIM) - Instituto de Física UFRGS

<https://www.ufrgs.br/hepsim/> 


