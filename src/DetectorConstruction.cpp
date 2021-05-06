// -----------------------------------------------------------------------------
//  G4_QPIX | DetectorConstruction.cpp
//
//  Definition of detector geometry and materials.
//   * Author: Justo Martin-Albo
//   * Creation date: 14 Aug 2019
// -----------------------------------------------------------------------------

#include "DetectorConstruction.h"
#include "TrackingSD.h"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"


DetectorConstruction::DetectorConstruction(): G4VUserDetectorConstruction()
{}

DetectorConstruction::~DetectorConstruction()
{}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // WORLD /////////////////////////////////////////////////

  G4double world_size = 20.*cm;
  G4Material* world_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

  G4Box* world_solid_vol =
    new G4Box("world.solid", world_size/2., world_size/2., world_size/2.);

  G4LogicalVolume* world_logic_vol =
    new G4LogicalVolume(world_solid_vol, world_mat, "world.logical");
  world_logic_vol->SetVisAttributes(G4VisAttributes::Invisible);

  G4VPhysicalVolume* world_phys_vol =
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
                      world_logic_vol, "world.physical", 0, false, 0, true);

  // DETECTOR //////////////////////////////////////////////
  // rough LAr volume in the UTAH 
  // accounts for both light tubes and mesh rings...
  G4double detector_width   = 2.3*m;
  G4double detector_height  = 6.0*m;
  G4double detector_length  = 3.6*m;
  G4double r_min    = 0.0*mm;
  G4double r_max    = 46.0*mm;
  G4double z_half   = 125.0*mm;
  G4double phi_min  = 0.0;
  G4double phi_max  = 2*M_PI;

  G4Material* detector_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");

  G4Tubs* detector_solid_vol =
    // new G4Box("detector.solid", detector_width/2., detector_height/2., detector_length/2.);
    new G4Tubs("detector.solid", r_min, r_max, z_half, phi_min, phi_max);

  G4LogicalVolume* detector_logic_vol =
    new G4LogicalVolume(detector_solid_vol, detector_mat, "detector.logical");

  // G4ThreeVector offset(detector_width/2., detector_height/2., detector_length/2.);

  // new G4PVPlacement(0, 0,
  //                   detector_logic_vol, "detector.physical", world_logic_vol, false, 0, true);

  //////////////////////////////////////////////////////////



  return world_phys_vol;
}

void DetectorConstruction::ConstructSDandField()
{
  // SENSITIVE DETECTOR ////////////////////////////////////

  TrackingSD* tracking_sd = new TrackingSD("/G4QPIX/TRACKING", "TrackingHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(tracking_sd);

  G4LogicalVolume* detector_logic_vol =
    G4LogicalVolumeStore::GetInstance()->GetVolume("detector.logical");

  SetSensitiveDetector(detector_logic_vol, tracking_sd);

  //////////////////////////////////////////////////////////
}
