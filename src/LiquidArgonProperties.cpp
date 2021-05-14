// ----------------------------------------------------------------------------
// nexus | LiquidArgonProperties.cc
//
// This class collects the relevant physical properties of liquid xenon.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "LiquidArgonProperties.h"

// #include <G4AnalyticalPolSolver.hh>
#include <G4MaterialPropertiesTable.hh>

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"
// #include <stdexcept>
// #include <cassert>

using CLHEP::h_Planck;
using CLHEP::c_light;

using namespace CLHEP;

LiquidArgonProperties::LiquidArgonProperties()
{
}



LiquidArgonProperties::~LiquidArgonProperties()
{
}



G4double LiquidArgonProperties::RefractiveIndex(G4double energy)
{
    // Formula for the refractive index taken from
    // https://arxiv.org/pdf/2002.09346.pdf
    G4double a0  = 0.335;
    G4double aUV = 0.099;
    G4double aIR = 0.008;
    G4double lambdaUV = 106.6 * nm;
    G4double lambdaIR = 908.3 * nm;
   
    G4double wavelength = h_Planck * c_light / energy;
    G4double x = a0 + aUV * pow (wavelength, 2) / (pow (wavelength, 2) - pow (lambdaUV, 2)) + aIR * pow (wavelength, 2) / (pow (wavelength, 2) - pow (lambdaIR, 2));
    G4double rindex_element = sqrt (1 + 3 * x / (3 - x));

    return rindex_element;
}



G4double LiquidArgonProperties::Scintillation(G4double energy)
{
    // Gaussian centered at 128 with a FWHM of 10
    // values taken from the LBL page
    // https://lar.bnl.gov/properties/
    G4double Wavelength_peak = 128.*nm;
    G4double Wavelength_FWHM = 10.*nm;
    G4double Wavelength_sigma = Wavelength_FWHM/2.35;

    G4double Energy_peak = (h_Planck*c_light/Wavelength_peak);
    G4double Energy_sigma = (h_Planck*c_light*Wavelength_sigma/pow(Wavelength_peak,2));

    G4double intensity =
        exp(-pow(Energy_peak/eV-energy/eV,2)/(2*pow(Energy_sigma/eV, 2)))/(Energy_sigma/eV*sqrt(pi*2.));

    return intensity;
}



void LiquidArgonProperties::Scintillation
(std::vector<G4double>& energy, std::vector<G4double>& intensity)
{
    for (unsigned i=0; i<energy.size(); i++)
        intensity.push_back(Scintillation(energy[i]));
}
