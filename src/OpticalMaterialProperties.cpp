// ----------------------------------------------------------------------------
// nexus | OpticalMaterialProperties.cc
//
// Optical properties of relevant materials.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "OpticalMaterialProperties.h"
#include "LiquidArgonProperties.h"

#include <G4MaterialPropertiesTable.hh>

#include <assert.h>

using namespace CLHEP;

G4MaterialPropertiesTable* OpticalMaterialProperties::OpticalLAr()
{
  LiquidArgonProperties LAr_prop;
  G4MaterialPropertiesTable* LAr_mpt = new G4MaterialPropertiesTable();

  const G4int ri_entries = 200;
  G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

  std::vector<G4double> ri_energy;
  for (int i=0; i<ri_entries; i++) 
  {
    ri_energy.push_back(optPhotMinE_ + i * eWidth);
  }

  std::vector<G4double> ri_index;
  for (G4int i=0; i<ri_entries; i++) 
  {
    ri_index.push_back(LAr_prop.RefractiveIndex(ri_energy[i]));
  }

  assert(ri_energy.size() == ri_index.size());
  LAr_mpt->AddProperty("RINDEX", ri_energy.data(), ri_index.data(), ri_energy.size());

  // Sampling from 128 +- 25 (103 nm --> 153 nm) [12.04 eV --> 8.104 eV]
  const G4int sc_entries = 500;
  eWidth = (12.04*eV - 8.104*eV) / sc_entries;
  
  std::vector<G4double> sc_energy;
  for (int j=0; j<sc_entries; j++)
  {
    sc_energy.push_back(8.104*eV + j * eWidth);
  }

  std::vector<G4double> intensity;
  LAr_prop.Scintillation(sc_energy, intensity);

  assert(sc_energy.size() == intensity.size());
  LAr_mpt->AddProperty("FASTCOMPONENT", sc_energy.data(), intensity.data(), sc_energy.size());
  LAr_mpt->AddProperty("SLOWCOMPONENT", sc_energy.data(), intensity.data(), sc_energy.size());

  // LAr_mpt->AddConstProperty("SCINTILLATIONYIELD", 51282.1/MeV);
  LAr_mpt->AddConstProperty("SCINTILLATIONYIELD", 2./MeV);
  LAr_mpt->AddConstProperty("RESOLUTIONSCALE", 1);
  LAr_mpt->AddConstProperty("RAYLEIGH", 95.*cm);
  LAr_mpt->AddConstProperty("FASTTIMECONSTANT", 6.*ns);
  LAr_mpt->AddConstProperty("SLOWTIMECONSTANT", 1590.*ns);

  // LAr_mpt->AddConstProperty("IONIZATIONENERGY", 42372.9/MeV);
  // LAr_mpt->AddConstProperty("FANOFACTOR", 0.107);
  // for (int i=0; i<100; ++i){G4cout << "MAT" <<"\t"<< LAr_mpt->GetConstProperty("IONIZATIONENERGY") <<  G4endl;} 

  // Yeild ratio from here
  // https://journals.aps.org/prd/pdf/10.1103/PhysRevD.103.043001
  LAr_mpt->AddConstProperty("YIELDRATIO", 0.84);

  std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
  std::vector<G4double> abs_length = {noAbsLength_, noAbsLength_};

  assert(abs_energy.size() == abs_length.size());
  LAr_mpt->AddProperty("ABSLENGTH", abs_energy.data(), abs_length.data(), abs_energy.size());

  return LAr_mpt;
}




/// Fake Grid ///
G4MaterialPropertiesTable* OpticalMaterialProperties::FakeGrid(G4double transparency,
                                                               G4double thickness,
                                                               G4int    sc_yield)
{
  G4MaterialPropertiesTable* mpt      = new G4MaterialPropertiesTable();

  // PROPERTIES FROM Liquid Argon
  G4MaterialPropertiesTable* LAr_pt = OpticalLAr();

  mpt->AddProperty("RINDEX",        LAr_pt->GetProperty("RINDEX"));
  mpt->AddProperty("FASTCOMPONENT", LAr_pt->GetProperty("FASTCOMPONENT"));
  mpt->AddProperty("SLOWCOMPONENT", LAr_pt->GetProperty("SLOWCOMPONENT"));

  mpt->AddConstProperty("SCINTILLATIONYIELD", LAr_pt->GetConstProperty("SCINTILLATIONYIELD"));
  mpt->AddConstProperty("RESOLUTIONSCALE",    LAr_pt->GetConstProperty("RESOLUTIONSCALE"));
  mpt->AddConstProperty("FASTTIMECONSTANT",   LAr_pt->GetConstProperty("FASTTIMECONSTANT"));
  mpt->AddConstProperty("SLOWTIMECONSTANT",   LAr_pt->GetConstProperty("SLOWTIMECONSTANT"));
  mpt->AddConstProperty("YIELDRATIO",         LAr_pt->GetConstProperty("YIELDRATIO"));
  mpt->AddConstProperty("ATTACHMENT",         LAr_pt->GetConstProperty("ATTACHMENT"));

  // ABSORPTION LENGTH
  G4double abs_length   = -thickness/log(transparency);
  G4double abs_energy[] = {optPhotMinE_, optPhotMaxE_};
  G4double absLength[]  = {abs_length, abs_length};
  mpt->AddProperty("ABSLENGTH", abs_energy, absLength, 2);

  // PHOTOELECTRIC REEMISSION
  // https://aip.scitation.org/doi/10.1063/1.1708797
  G4double stainless_wf = 4.3 * eV; // work function
  mpt->AddConstProperty("WORK_FUNCTION", stainless_wf);
  
  return mpt;
}




/// PTFE (== TEFLON) ///
G4MaterialPropertiesTable* OpticalMaterialProperties::PTFE()
{
  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

  // REFLECTIVITY
  const G4int REFL_NUMENTRIES = 7;
  G4double ENERGIES[REFL_NUMENTRIES] = {
    optPhotMinE_,  2.8 * eV,  3.5 * eV,  4. * eV,
    6. * eV,       7.2 * eV,  optPhotMaxE_
  };
  G4double REFLECTIVITY[REFL_NUMENTRIES] = {
    .98,  .98,  .98,  .98,
    .72,  .72,  .72
  };
  mpt->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY, REFL_NUMENTRIES);

  // REFLEXION BEHAVIOR
  const G4int NUMENTRIES = 2;
  G4double ENERGIES_2[NUMENTRIES]    = {optPhotMinE_, optPhotMaxE_};
  // Specular reflection about the normal to a microfacet.
  // Such a vector is chosen according to a gaussian distribution with
  // sigma = SigmaAlhpa (in rad) and centered in the average normal.
  G4double specularlobe[NUMENTRIES]  = {0., 0.};
  // specular reflection about the average normal
  G4double specularspike[NUMENTRIES] = {0., 0.};
  // 180 degrees reflection.
  G4double backscatter[NUMENTRIES]   = {0., 0.};
  // 1 - the sum of these three last parameters is the percentage of Lambertian reflection

  mpt->AddProperty("SPECULARLOBECONSTANT", ENERGIES_2, specularlobe,  NUMENTRIES);
  mpt->AddProperty("SPECULARSPIKECONSTANT",ENERGIES_2, specularspike, NUMENTRIES);
  mpt->AddProperty("BACKSCATTERCONSTANT",  ENERGIES_2, backscatter,   NUMENTRIES);

  // REFRACTIVE INDEX
  G4double rIndex[] = {1.41, 1.41};
  mpt->AddProperty("RINDEX", ENERGIES_2, rIndex, NUMENTRIES);

  return mpt;
}





/// TPB (tetraphenyl butadiene) ///
G4MaterialPropertiesTable* OpticalMaterialProperties::TPB()
{
  // Data from https://doi.org/10.1140/epjc/s10052-018-5807-z
  G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

  // REFRACTIVE INDEX
  const G4int rIndex_numEntries = 2;
  G4double rIndex_energies[rIndex_numEntries] = {optPhotMinE_, optPhotMaxE_};
  G4double TPB_rIndex[rIndex_numEntries]      = {1.67    , 1.67};
  mpt->AddProperty("RINDEX", rIndex_energies,
                   TPB_rIndex, rIndex_numEntries);

  // ABSORPTION LENGTH
  // Assuming no absorption except WLS
  G4double abs_energy[] = {optPhotMinE_, optPhotMaxE_};
  G4double absLength[]  = {noAbsLength_, noAbsLength_};
  mpt->AddProperty("ABSLENGTH", abs_energy, absLength, 2);

  // WLS ABSORPTION LENGTH (Version NoSecWLS)
  // The NoSecWLS is forced by setting the WLS_absLength to noAbsLength_
  // for wavelengths higher than 380 nm where the WLS emission spectrum starts.
  G4double WLS_abs_energy[] = {
    optPhotMinE_,
    h_Planck * c_light / (380. * nm),  h_Planck * c_light / (370. * nm),
    h_Planck * c_light / (360. * nm),  h_Planck * c_light / (330. * nm),
    h_Planck * c_light / (320. * nm),  h_Planck * c_light / (310. * nm),
    h_Planck * c_light / (300. * nm),  h_Planck * c_light / (270. * nm),
    h_Planck * c_light / (250. * nm),  h_Planck * c_light / (230. * nm),
    h_Planck * c_light / (210. * nm),  h_Planck * c_light / (190. * nm),
    h_Planck * c_light / (170. * nm),  h_Planck * c_light / (150. * nm),
    h_Planck * c_light / (100. * nm),  optPhotMaxE_
  };
  const G4int WLS_abs_entries = sizeof(WLS_abs_energy) / sizeof(G4double);

  G4double WLS_absLength[] = {
    noAbsLength_,
    noAbsLength_,   50. * nm,     // 380 , 370 nm
     30. * nm,      30. * nm,     // 360 , 330 nm
     50. * nm,      80. * nm,     // 320 , 310 nm
    100. * nm,     100. * nm,     // 300 , 270 nm
    400. * nm,     400. * nm,     // 250 , 230 nm
    350. * nm,     250. * nm,     // 210 , 190 nm
    350. * nm,     400. * nm,     // 170 , 150 nm
    400. * nm,     noAbsLength_   // 100 nm
  };

  //for (int i=0; i<WLS_abs_entries; i++)
  //  G4cout << "* TPB WLS absLength:  " << std::setw(8) << WLS_abs_energy[i] / eV
  //         << " eV  ==  " << std::setw(8) << (h_Planck * c_light / WLS_abs_energy[i]) / nm
  //         << " nm  ->  " << std::setw(6) << WLS_absLength[i] / nm << " nm" << G4endl;
  assert(sizeof(WLS_absLength) == sizeof(WLS_abs_energy));
  mpt->AddProperty("WLSABSLENGTH", WLS_abs_energy,
                   WLS_absLength,  WLS_abs_entries);

  // WLS EMISSION SPECTRUM
  // Implemented with formula (7), with parameter values in table (3)
  // Sampling from ~380 nm to 600 nm <--> from 2.06 to 3.26 eV
  const G4int WLS_emi_entries = 120;
  G4double WLS_emi_energy[WLS_emi_entries];
  for (int i=0; i<WLS_emi_entries; i++)
    WLS_emi_energy[i] = 2.06 * eV + 0.01 * i * eV;

  G4double WLS_emiSpectrum[WLS_emi_entries];
  G4double A      = 0.782;
  G4double alpha  = 3.7e-2;
  G4double sigma1 = 15.43;
  G4double mu1    = 418.10;
  G4double sigma2 = 9.72;
  G4double mu2    = 411.2;

  for (int i=0; i<WLS_emi_entries; i++) {
    G4double wl = (h_Planck * c_light / WLS_emi_energy[i]) / nm;
    WLS_emiSpectrum[i] = A * (alpha/2.) * exp((alpha/2.) *
                         (2*mu1 + alpha*pow(sigma1,2) - 2*wl)) *
                         erfc((mu1 + alpha*pow(sigma1,2) - wl) / (sqrt(2)*sigma1)) +
                         (1-A) * (1 / sqrt(2*pow(sigma2,2)*3.1416)) *
                         exp((-pow(wl-mu2,2)) / (2*pow(sigma2,2)));
    // G4cout << "* TPB WLSemi:  " << std::setw(4)
    //        << wl << " nm -> " << WLS_emiSpectrum[i] << G4endl;
  };
  assert(sizeof(WLS_emiSpectrum) == sizeof(WLS_emi_energy));
  mpt->AddProperty("WLSCOMPONENT", WLS_emi_energy,
                   WLS_emiSpectrum, WLS_emi_entries);

  // WLS Delay
  mpt->AddConstProperty("WLSTIMECONSTANT", 1.2 * ns);

  // WLS Quantum Efficiency
  // According to the paper, the QE of TPB depends on the incident wavelength.
  // As Geant4 doesn't allow this possibility, it is set to the value corresponding
  // to Xe scintillation spectrum peak.
  mpt->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.65);

  return mpt;
}

