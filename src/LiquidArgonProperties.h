// ----------------------------------------------------------------------------
// nexus | LiquidArgonProperties.cc
//
// This class collects the relevant physical properties of liquid xenon.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef XENON_LIQUID_PROPERTIES_H
#define XENON_LIQUID_PROPERTIES_H 1

#include <globals.hh>
#include <vector>

class G4MaterialPropertiesTable;


class LiquidArgonProperties
{
public:
/// Constructor
LiquidArgonProperties();
/// Destructor
~LiquidArgonProperties();

// G4double Density();

/// Return the refractive index of xenon gas for a given photon energy
G4double RefractiveIndex(G4double energy);

G4double Scintillation(G4double energy);
void Scintillation(std::vector<G4double>& energy, std::vector<G4double>& intensity);

private:
// G4double density_;

};

#endif