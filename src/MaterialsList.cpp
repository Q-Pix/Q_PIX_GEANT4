// ----------------------------------------------------------------------------
// nexus | MaterialsList.cc
//
// Definition of materials of common use.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "MaterialsList.h"

#include <G4Material.hh>
#include <G4Element.hh>
#include <G4NistManager.hh>

#include <CLHEP/Units/SystemOfUnits.h>



G4Material* MaterialsList::OLAr()
{
  G4String name = "OLAr";

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    mat = new G4Material(name, 1.3973*CLHEP::g/CLHEP::cm3, 1,
			 kStateLiquid, 87*CLHEP::kelvin, 1*CLHEP::bar);

    G4Element* Ar = nist->FindOrBuildElement("Ar");

    mat->AddElement(Ar,1);
  }

  return mat;
}

G4Material* MaterialsList::Steel()
{
  // Composition ranges correspond to stainless steel grade 304L

  G4String name = "Steel";

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {

    mat = new G4Material(name, 8000*kg/m3, 4);

    G4NistManager* nist = G4NistManager::Instance();

    G4Element* Fe = nist->FindOrBuildElement("Fe");
    mat->AddElement(Fe, 0.66);

    G4Element* Cr = nist->FindOrBuildElement("Cr");
    mat->AddElement(Cr, 0.20);

    G4Element* Mn = nist->FindOrBuildElement("Mn");
    mat->AddElement(Mn, 0.02);

    G4Element* Ni = nist->FindOrBuildElement("Ni");
    mat->AddElement(Ni, 0.12);
  }

  return mat;
}



G4Material* MaterialsList::Steel316Ti()
{
  G4String name = "Steel316Ti";

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (!mat) {
    mat = new G4Material(name, 8000*kg/m3, 6);

    G4NistManager* nist = G4NistManager::Instance();

    G4Element* Fe = nist->FindOrBuildElement("Fe");
    mat->AddElement(Fe, 0.636);

    G4Element* Cr = nist->FindOrBuildElement("Cr");
    mat->AddElement(Cr, 0.18);

    G4Element* Ni = nist->FindOrBuildElement("Ni");
    mat->AddElement(Ni, 0.14);

    G4Element* Mo = nist->FindOrBuildElement("Mo");
    mat->AddElement(Mo, 0.03);

    G4Element* Si = nist->FindOrBuildElement("Si");
    mat->AddElement(Si, 0.007);

    G4Element* Ti = nist->FindOrBuildElement("Ti");
    mat->AddElement(Ti, 0.007);
  }

  return mat;
}



G4Material* MaterialsList::Epoxy()
{
  // Definition taken from the Geant4 advanced example "Composite Calorimeter"
  // (Geant4/examples/advanced/composite_calorimeter/dataglobal/material.cms).

  G4String name = "Epoxy";

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");
    G4Element* O = nist->FindOrBuildElement("O");

    mat = new G4Material(name, 1.3*g/cm3, 3);
    mat->AddElement(H, 44);
    mat->AddElement(C, 15);
    mat->AddElement(O, 7);
  }

  return mat;
}



G4Material* MaterialsList::Kovar()
{
  G4String name = "Kovar";

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* Fe = nist->FindOrBuildElement("Fe");
    G4Element* Ni = nist->FindOrBuildElement("Ni");
    G4Element* Co = nist->FindOrBuildElement("Co");

    mat = new G4Material(name, 8.35*g/cm3, 3);
    mat->AddElement(Fe, 54);
    mat->AddElement(Ni, 29);
    mat->AddElement(Co, 17);
  }

  return mat;
}



G4Material* MaterialsList::PEEK()
{
  G4String name = "PEEK";

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");
    G4Element* O = nist->FindOrBuildElement("O");

    mat = new G4Material(name, 1.3*g/cm3, 3, kStateSolid);
    mat->AddElement(H, 12);
    mat->AddElement(C, 18);
    mat->AddElement(O,  3);
  }

  return mat;
}



G4Material* MaterialsList::Sapphire()
{
  G4Material* mat =
    G4NistManager::Instance()->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");

  return mat;
}



G4Material* MaterialsList::FusedSilica()
{
  G4Material* mat =
    G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  return mat;
}



G4Material* MaterialsList::PS() // polystyrene
{
  G4String name = "PS";

  // Check whether material exists already in the materials table
  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");
    //      G4Element* O = nist->FindOrBuildElement("O");

    //      mat = new G4Material(name, 1320*kg/m3, 3, kStateSolid);
    mat = new G4Material(name, 1.05*g/cm3, 2, kStateSolid);
    mat->AddElement(H, 8);
    mat->AddElement(C, 8);
    //      mat->AddElement(O,  3);
  }

  return mat;
}



G4Material* MaterialsList::TPB()
{
  G4String name = "TPB"; // Tetraphenyl butadiene

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");

    mat = new G4Material(name, 1*g/cm3, 2, kStateSolid);
    mat->AddElement(H, 22);
    mat->AddElement(C, 28);
  }

  return mat;
}


G4Material* MaterialsList::ITO()
{
  // Density data from  https://doi.org/10.1116/1.1371326
  // Composition data from: https://en.wikipedia.org/wiki/Indium_tin_oxide

  G4String name = "ITO";

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* In = nist->FindOrBuildElement("In");
    G4Element* O  = nist->FindOrBuildElement("O");
    G4Element* Sn = nist->FindOrBuildElement("Sn");

    mat = new G4Material(name, 6.8*g/cm3, 3, kStateSolid);
    mat->AddElement(In, 74. * perCent);
    mat->AddElement(O , 18. * perCent);
    mat->AddElement(Sn,  8. * perCent);
  }
  return mat;
}


G4Material* MaterialsList::PVT()
{
  G4String name = "PVT"; //

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");

    mat = new G4Material(name, 1*g/cm3, 2, kStateSolid);
    mat->AddElement(H, 4);
    mat->AddElement(C, 2);
  }

  return mat;

}


G4Material* MaterialsList::Kevlar()
{
  G4String name = "Kevlar"; //

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");
    G4Element* N = nist->FindOrBuildElement("N");
    G4Element* O = nist->FindOrBuildElement("O");

    mat = new G4Material(name, 1.44*g/cm3, 4, kStateSolid);
    mat->AddElement(H, 10);
    mat->AddElement(C, 14);
    mat->AddElement(O, 2);
    mat->AddElement(N, 2);
  }

  return mat;

}


G4Material* MaterialsList::HDPE()
{
  G4String name = "HDPE";

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");

    mat = new G4Material(name, .954*g/cm3, 2, kStateSolid);
    mat->AddElement(H, 4);
    mat->AddElement(C, 2);
  }

  return mat;

}


G4Material* MaterialsList::OpticalSilicone()
{
  // Silicone resin with a methyl group
  // (https://en.wikipedia.org/wiki/Silicone_resin)

  G4String name = "OpticalSilicone";

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H  = nist->FindOrBuildElement("H");
    G4Element* C  = nist->FindOrBuildElement("C");
    G4Element* O  = nist->FindOrBuildElement("O");
    G4Element* Si = nist->FindOrBuildElement("Si");

    mat = new G4Material(name, 1.05*g/cm3, 4, kStateSolid);
    mat->AddElement(H,  3);
    mat->AddElement(C,  1);
    mat->AddElement(Si, 1);
    mat->AddElement(O,  1);
  }

  return mat;
}


G4Material* MaterialsList::SeF6(G4double pressure, G4double temperature)
{
  // Composition ranges correspond to Selenium Hexafluoride

  G4String name = "SeF6";
  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {

    G4NistManager* nist = G4NistManager::Instance();
    G4double density = 7.887*kg/m3;

    if (pressure/bar > 9.5 && pressure/bar < 10.5)
      density *= 10.;
    else
      G4cout  << "[MaterialsList] Pressure " << pressure/bar
              << " bar not recognized for SeF6! ... Assuming 1bar. " << G4endl;

    mat = new G4Material(name, density, 2, kStateGas, temperature, pressure);
    G4Element* Se = nist->FindOrBuildElement("Se");
    mat->AddElement(Se, 1);
    G4Element* F = nist->FindOrBuildElement("F");
    mat->AddElement(F, 6);
  }

  return mat;
}


G4Material* MaterialsList::FR4()
{
  // FR-4 is a composite material widely used for printed circuits boards.
  // It consists of woven fiberglass cloth with an epoxy resin binder that is
  // flame resistant. Typical proportions are 60% fused silica and 40% epoxy.

  G4String name = "FR4";
  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    mat = new G4Material(name, 1.850*g/cm3, 2, kStateSolid);
    mat->AddMaterial(MaterialsList::FusedSilica(), 0.60);
    mat->AddMaterial(MaterialsList::Epoxy(),       0.40);
  }

  return mat;
}


G4Material* MaterialsList::Limestone()
{
  // Limestone modelled as pure calcium carbonate
  G4String name = "Limestone";

  G4Material* mat = G4Material::GetMaterial(name, false);
  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* Ca = nist->FindOrBuildElement("Ca");
    G4Element* C  = nist->FindOrBuildElement("C");
    G4Element* O  = nist->FindOrBuildElement("O");
    mat = new G4Material(name, 2.711*g/cm3, 3, kStateSolid);
    mat->AddElement(Ca, 1);
    mat->AddElement(C , 1);
    mat->AddElement(O , 1);
  }
  return mat;
}


G4Material* MaterialsList::CopyMaterial(G4Material* original, const G4String& newname)
{
  G4Material* newmat = G4Material::GetMaterial(newname, false);

  if (newmat == 0) {

    G4double density     = original->GetDensity();
    G4double temperature = original->GetTemperature();
    G4double pressure    = original->GetPressure();
    G4State  state       = original->GetState();
    G4int    n_elem      = original->GetNumberOfElements();

    if (n_elem == 1) {
      G4double z = original->GetZ();
      G4double a = original->GetA();
      newmat = new G4Material(newname, z, a, density, state, temperature, pressure);
    }
    else {
      const G4double* fractions = original->GetFractionVector();
      newmat = new G4Material(newname, density, n_elem, state, temperature, pressure);
      for (G4int i = 0; i < n_elem; ++i)
        newmat->AddElement(new G4Element(original->GetElement(i)->GetName(),
                                         original->GetElement(i)->GetSymbol(),
                                         original->GetElement(i)->GetZ(),
                                         original->GetElement(i)->GetA()),
                                         fractions[i]);
      }
  }

  return newmat;
}


G4Material* MaterialsList::FakeDielectric(G4Material* model_mat, G4String name)
{
  return CopyMaterial(model_mat, name);
}



G4Material* MaterialsList::TPH()
{
  G4String name = "TPH"; // p - terphenyl

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");

    mat = new G4Material(name, 1.24*g/cm3, 2, kStateSolid);
    mat->AddElement(H, 14);
    mat->AddElement(C, 18);
  }

  return mat;
}


// WLS EJ-280
G4Material* MaterialsList::EJ280()
{
  G4String name = "EJ280"; //

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    // The base is Polyvinyltoluene
    // Linear formula: [CH2CH(C6H4CH3)]n
    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");

    mat = new G4Material(name, 1.023*g/cm3, 2, kStateSolid);
    mat->AddElement(H, 10);
    mat->AddElement(C, 9);
  }

  return mat;
}


// Kuraray Y-11
G4Material* MaterialsList::Y11()
{
  G4String name = "Y11"; //

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    // The base is Polystyrene
    // Linear formula: (C8H8)n
    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");

    mat = new G4Material(name, 1.05* g / cm3, 2, kStateSolid);
    mat->AddElement(H, 8);
    mat->AddElement(C, 8);
  }

  return mat;
}


// Pethylene (cladding material)
G4Material* MaterialsList::Pethylene()
{
  G4String name = "Pethylene"; //

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");

    mat = new G4Material(name, 1.200 * g/cm3, 2, kStateSolid);
    mat->AddElement(H, 4);
    mat->AddElement(C, 2);
  }

  return mat;
}


// FPethylene (fluorinated polyethylene - cladding material)
G4Material* MaterialsList::FPethylene()
{
  G4String name = "FPethylene"; //

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");

    mat = new G4Material(name, 1.400 * g/cm3, 2, kStateSolid);
    mat->AddElement(H, 4);
    mat->AddElement(C, 2);
  }

  return mat;
}



// PMMA == PolyMethylmethacrylate (cladding)
G4Material* MaterialsList::PMMA()
{
  G4String name = "PMMA"; //

  G4Material* mat = G4Material::GetMaterial(name, false);

  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");
    G4Element* O = nist->FindOrBuildElement("O");

    mat = new G4Material(name, 1.190 * g/cm3, 3, kStateSolid);
    mat->AddElement(H, 8);
    mat->AddElement(C, 5);
    mat->AddElement(O, 2);
  }

  return mat;
}


G4Material* MaterialsList::PEDOT()
{
  // Poly(3,4-ethylenedioxythiophene) is a conductive polymer.
  // Linear Formula: (C6H4O2S)n
  // Density of dry coating with poly(styrenesulfonate): 1.011 g/cm3

  G4String name = "PEDOT";
  G4Material* material = G4Material::GetMaterial(name, false);

  if (!material) {
    G4NistManager* nist = G4NistManager::Instance();

    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");
    G4Element* O = nist->FindOrBuildElement("O");
    G4Element* S = nist->FindOrBuildElement("S");

    material = new G4Material(name, 1.011*g/cm3, 4, kStateSolid);
    material->AddElement(H,4);
    material->AddElement(C,6);
    material->AddElement(O,2);
    material->AddElement(S,1);
  }

  return material;
}