// -----------------------------------------------------------------------------
//  G4_QPIX | DetectorConstruction.cpp
//
//  Definition of detector geometry and materials.
//   * Author: Justo Martin-Albo
//   * Creation date: 14 Aug 2019
// -----------------------------------------------------------------------------

#include "DetectorConstruction.h"
#include "TrackingSD.h"
#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "GenericPhotosensor.h"


#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GenericMessenger.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

DetectorConstruction::DetectorConstruction(): G4VUserDetectorConstruction(), Detector_Geometry_("APA")
{
  messenger_ = new G4GenericMessenger(this, "/Inputs/", "Set the geometry.");
  messenger_->DeclareProperty("Geometry", Detector_Geometry_,  "which geometry?");
}

DetectorConstruction::~DetectorConstruction()
{
  delete messenger_;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  // WORLD /////////////////////////////////////////////////
  G4double world_size = 55.*cm;
  G4Material* world_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");

  G4Box* world_solid_vol =
    new G4Box("world.solid", world_size/2., world_size/2., world_size/2.);

  G4LogicalVolume* world_logic_vol =
    new G4LogicalVolume(world_solid_vol, world_mat, "world.logical");
  world_logic_vol->SetVisAttributes(G4VisAttributes::Invisible);

  G4VPhysicalVolume* world_phys_vol =
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
                      world_logic_vol, "world.physical", 0, false, 0, true);


  // DETECTOR //////////////////////////////////////////////
  G4double r_min    = 0.0*cm;
  G4double r_max    = 4.6*cm;
  G4double z_half   = 12.5*cm;
  G4double phi_min  = 0.0;
  G4double phi_max  = 2*M_PI;

  G4Material* detector_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
  detector_mat->SetMaterialPropertiesTable(OpticalMaterialProperties::LAr());

  G4Tubs* detector_solid_vol =
    new G4Tubs("detector.solid", r_min, r_max, z_half, phi_min, phi_max);

  G4LogicalVolume* detector_logic_vol =
    new G4LogicalVolume(detector_solid_vol, detector_mat, "detector.logical");

  G4ThreeVector offset(0, 0, z_half);

  new G4PVPlacement(0, offset,
                    detector_logic_vol, "detector.physical", world_logic_vol, false, 0, true);
  //////////////////////////////////////////////////////////



  // LIGHT TUBE ////////////////////////////////////////////
  G4double teflon_thickn_ = 2*cm;
  G4double tpb_thickn_    = 1*micrometer;
  G4double r_inner  = r_max; 
  G4double r_outter = r_inner + teflon_thickn_;

  // THE TEFLON //
  G4Material* teflon_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");

  G4Tubs* teflon_drift_solid =
    new G4Tubs("LIGHT_TUBE", r_inner, r_outter, z_half, phi_min, phi_max);

  G4LogicalVolume* teflon_drift_logic =
    new G4LogicalVolume(teflon_drift_solid, teflon_, "LIGHT_TUBE");

  new G4PVPlacement(0, offset,
                    teflon_drift_logic, "LIGHT_TUBE", world_logic_vol, false, 0, false);

  // THE TPB ON TEFLON//
  G4Material* tpb_ = MaterialsList::TPB();
  tpb_->SetMaterialPropertiesTable(OpticalMaterialProperties::TPB());

  G4Tubs* tpb_drift_solid =
    new G4Tubs("TPB_TUBE", r_inner, (r_inner + tpb_thickn_), z_half, phi_min, phi_max);

  G4LogicalVolume* tpb_drift_logic =
    new G4LogicalVolume(tpb_drift_solid, tpb_, "TPB_TUBE");
  G4VPhysicalVolume* tpb_drift_phys =
    new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), tpb_drift_logic,
                      "TPB_TUBE", teflon_drift_logic, false, 0, false);

  new G4PVPlacement(0, offset,
                    tpb_drift_logic, "TPB_TUBE", world_logic_vol, false, 0, false);

  /// Optical surface on teflon ///
  G4OpticalSurface* refl_Surf =
    new G4OpticalSurface("refl_Surf", unified, ground, dielectric_metal, .01);
  refl_Surf->SetMaterialPropertiesTable(OpticalMaterialProperties::PTFE());
  new G4LogicalSkinSurface("refl_teflon_surf", teflon_drift_logic, refl_Surf);
  // new G4LogicalSkinSurface("refl_teflon_surf", teflon_buffer_logic, refl_Surf);

  /// Optical surface between LAr and TPB to model roughness ///
  G4OpticalSurface* lar_tpb_teflon_surf =
    new G4OpticalSurface("gas_tpb_teflon_surf", glisur, ground,
                         dielectric_dielectric, .01);

  new G4LogicalBorderSurface("gas_tpb_teflon_surf", tpb_drift_phys, world_phys_vol,
                             lar_tpb_teflon_surf);
  new G4LogicalBorderSurface("gas_tpb_teflon_surf", world_phys_vol, tpb_drift_phys,
                             lar_tpb_teflon_surf);
                             
                             
                            
  // PHOTOMULTIPLIER ///////////////////////////////////////////////
  GenericPhotosensor gen_pmt('pmt', 10*cm, 5*mm, 1*mm);
  gen_pmt.SetWithWLSCoating(true);
  gen_pmt.SetSensorDepth(1);
  // gen_pmt.SetOpticalProperties(gen_pmt.GetMyPhotOptSurf());
  gen_pmt.Construct();


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
