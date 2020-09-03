#ifndef MUON_GENERATOR_H
#define MUON_GENERATOR_H

#include <G4VPrimaryGenerator.hh>

#include "G4Box.hh"

class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;


// class BaseGeometry;

class MuonGenerator: public G4VPrimaryGenerator
{
public:
/// Constructor
MuonGenerator();
/// Destructor
~MuonGenerator();

/// This method is invoked at the beginning of the event. It sets
/// a primary vertex (that is, a particle in a given position and time)
/// in the event.
void GeneratePrimaryVertex(G4Event*);

private:

/// Generate a random kinetic energy with flat probability in
//  the interval [energy_min, energy_max].
G4double RandomEnergy() const;
G4String MuonCharge() const;
G4double GetPhi() const;
G4double GetTheta() const;

private:
G4GenericMessenger* msg_;

G4ParticleDefinition* particle_definition_;

G4double energy_min_; ///< Minimum kinetic energy
G4double energy_max_; ///< Maximum kinetic energy

// G4String region_;

// const BaseGeometry* geom_; ///< Pointer to the detector geometry

G4double momentum_X_;
G4double momentum_Y_;
G4double momentum_Z_;

G4double position_X_;
G4double position_Y_;
G4double position_Z_;

double detector_length_x_;
double detector_length_y_;
double detector_length_z_;

G4Box* detector_solid_vol_;

};

#endif
