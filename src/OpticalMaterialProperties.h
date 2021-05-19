// ----------------------------------------------------------------------------
// nexus | OpticalMaterialProperties.h
//
// Optical properties of relevant materials.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef OPTICAL_MATERIAL_PROPERTIES_H
#define OPTICAL_MATERIAL_PROPERTIES_H 1

#include <globals.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

class G4MaterialPropertiesTable;

using namespace CLHEP;

// This is a stateless class where all methods are static functions.
class OpticalMaterialProperties
{
public:

static G4MaterialPropertiesTable* OpticalLAr();


static G4MaterialPropertiesTable* FakeGrid(G4double transparency=.9,
                                            G4double thickness=1.*mm,
                                            G4int sc_yield=25510/MeV);

static G4MaterialPropertiesTable* TPB();

static G4MaterialPropertiesTable* PTFE();


private:

static constexpr G4double optPhotMinE_ =  0.2  * eV;
static constexpr G4double optPhotMaxE_ = 11.5  * eV;
static constexpr G4double noAbsLength_ = 1.e8  * m;

// Constant that allows to convert nm to eV:
// nm_to_eV_ / wavelength (nm) = energy (eV)
static constexpr G4double nm_to_eV_ = h_Planck * c_light * 1.e6;


private:
// Constructor (hidden)
OpticalMaterialProperties();

// Destructor (hidden)
~OpticalMaterialProperties();
};

#endif