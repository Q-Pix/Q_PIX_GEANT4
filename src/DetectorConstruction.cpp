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
// #include "GenericPhotosensor.h"
#include "PmtSD.h"


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

DetectorConstruction::DetectorConstruction(): G4VUserDetectorConstruction()//, Detector_Geometry_("APA")
{
  // messenger_ = new G4GenericMessenger(this, "/Inputs/", "Set the geometry.");
  // messenger_->DeclareProperty("Geometry", Detector_Geometry_,  "which geometry?");
}

DetectorConstruction::~DetectorConstruction()
{
  // delete messenger_;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  // WORLD /////////////////////////////////////////////////
  G4double world_size = 55.*cm;
  // G4Material* world_mat_hold = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
  // G4Material* world_mat = MaterialsList::CopyMaterial( world_mat_hold, 'null_LAr');

  // G4Material* world_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4Material* world_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_lAr");
  // world_mat->SetMaterialPropertiesTable(OpticalMaterialProperties::OpticalLAr());

  // G4Material* world_mat = MaterialsList::OLAr();
  // world_mat->SetMaterialPropertiesTable(OpticalMaterialProperties::OpticalLAr());
  // G4Tubs* world_solid_vol =
  //   new G4Tubs("world.solid", 0.*cm, 8*cm, 20*cm, 0., 2*M_PI);

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

  // G4Material* detector_mat = MaterialsList::CopyMaterial( world_mat, 'OpticalLAr');

  detector_mat->SetMaterialPropertiesTable(OpticalMaterialProperties::OpticalLAr());

  // G4MaterialPropertiesTable* mpt = detector_mat->GetMaterialPropertiesTable();
  // for (int i=0; i<100; ++i){G4cout << "MAT" <<"\t"<< mpt->GetConstProperty("IONIZATIONENERGY") <<  G4endl;} 

  G4Tubs* detector_solid_vol =
    new G4Tubs("detector.solid", r_min, r_max, z_half, phi_min, phi_max);

  G4LogicalVolume* detector_logic_vol =
    new G4LogicalVolume(detector_solid_vol, detector_mat, "detector.logical");
  // G4LogicalVolume* detector_logic_vol =
  //   new G4LogicalVolume(detector_solid_vol, world_mat, "detector.logical");

  // G4ThreeVector offset(0, 0, z_half);
  G4ThreeVector offset(0, 0, 0.);

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

  /// Optical surface between LAr and TPB to model roughness ///
  G4OpticalSurface* lar_tpb_teflon_surf =
    new G4OpticalSurface("gas_tpb_teflon_surf", glisur, ground,
                         dielectric_dielectric, .01);

  new G4LogicalBorderSurface("gas_tpb_teflon_surf", tpb_drift_phys, world_phys_vol,
                             lar_tpb_teflon_surf);
  new G4LogicalBorderSurface("gas_tpb_teflon_surf", world_phys_vol, tpb_drift_phys,
                             lar_tpb_teflon_surf);
                             
                             
                            
  // PHOTOMULTIPLIER ///////////////////////////////////////////////

  // GenericPhotosensor gen_pmt('pmt', 10*cm, 5*mm, 1*mm);
  // gen_pmt.SetWithWLSCoating(true);
  // gen_pmt.SetSensorDepth(1);
  // gen_pmt.Construct();


  // PHOTOSENSITIVE AREA /////////////////////////////////////////////

  // Sensitive /////
  G4Material* sensitive_mat_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

  G4Box* sensarea_solid_vol =
    new G4Box("SENSAREA", 10.*cm, 10.*cm, 5.*mm);

  G4LogicalVolume* sensarea_logic_vol =
    new G4LogicalVolume(sensarea_solid_vol, sensitive_mat_, "SENSAREA");

  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -1.*z_half), sensarea_logic_vol, "SENSAREA", world_logic_vol, false, 0, false);

  // making it photon sensitive
  G4MaterialPropertiesTable* sensitive_mpt_ = GetMyPhotOptSurf();
  G4OpticalSurface* sensitive_opsurf =
    new G4OpticalSurface("SENSAREA_optSurf", unified, polished, dielectric_metal);
  sensitive_opsurf->SetMaterialPropertiesTable(sensitive_mpt_);
  new G4LogicalSkinSurface("SENSAREA_optSurf", sensarea_logic_vol, sensitive_opsurf);


  return world_phys_vol;
}

void DetectorConstruction::ConstructSDandField()
{
  // SENSITIVE DETECTOR ////////////////////////////////////
  // PmtSD* pmtsd = new PmtSD("/G4QPIX/PHOTON", "PmtHitsCollection");
  // G4SDManager::GetSDMpointer()->AddNewDetector(pmtsd);

  // G4LogicalVolume* optical_logic_vol =
  //   G4LogicalVolumeStore::GetInstance()->GetVolume("SENSAREA");

  // SetSensitiveDetector(optical_logic_vol, pmtsd);




  // TrackingSD* tracking_sd = new TrackingSD("/G4QPIX/TRACKING", "TrackingHitsCollection");
  // G4SDManager::GetSDMpointer()->AddNewDetector(tracking_sd);

  // G4LogicalVolume* detector_logic_vol =
  //   G4LogicalVolumeStore::GetInstance()->GetVolume("detector.logical");

  // SetSensitiveDetector(detector_logic_vol, tracking_sd);

  //////////////////////////////////////////////////////////
}




// G4OpticalSurface* DetectorConstruction::GetMyPhotOptSurf()
G4MaterialPropertiesTable* DetectorConstruction::GetMyPhotOptSurf()
{
  const G4int entries = 57;
  G4double ENERGIES[entries] =
    {	2.06640321682 *eV, 2.10142700016 *eV,  2.13765850016 *eV,  2.17516128086 *eV,
  2.21400344659 *eV, 2.25425805471 *eV,  2.29600357425 *eV,  2.3393243964 *eV,
  2.38431140402 *eV, 2.43106260802 *eV,  2.47968386018 *eV,  2.53028965325 *eV,
  2.58300402103 *eV, 2.63796155339 *eV,  2.69530854368 *eV,  2.75520428909 *eV,
  2.81782256839 *eV, 2.8833533258 *eV,   2.95200459546 *eV,  3.02400470754 *eV,
  3.09960482523 *eV, 3.17908187203 *eV,  3.2627419213 *eV,   3.35092413538 *eV,
  3.44400536137 *eV, 3.54240551455 *eV,  3.64659391204 *eV,  3.75709675786 *eV,
  3.87450603154 *eV, 3.99949009707 *eV,  4.13280643364 *eV,  4.27531700032 *eV,
  4.42800689319 *eV, 4.59200714849 *eV,  4.76862280805 *eV,  4.95936772037 *eV,
  5.16600804205 *eV, 5.39061708736 *eV,  5.63564513678 *eV,  5.90400919092 *eV,
  6.19920965046 *eV, 6.358163744063561 *eV, 6.525483842591549 *eV, 6.7018482707697 *eV,
  6.888010722735524 *eV, 7.084811029099397 *eV, 7.293187824072908 *eV, 7.314701652462505 *eV,
  7.336342781611801 *eV, 7.358112344761985 *eV, 7.380011488645204 *eV, 7.402041373685937 *eV,
  7.424203174205955 *eV, 7.446498078633 *eV, 7.4689272897132195 *eV, 7.491492024727459 *eV,
  7.51419351571 *eV};

  G4double EFFICIENCY[entries] =
    { 0.0530,	0.0625, 0.0720, 0.0850,
  0.1050, 0.1190, 0.1335, 0.1550,
  0.1770, 0.1970, 0.2100, 0.2200,
  0.2300,	0.2430, 0.2580, 0.2770,
  0.2920,	0.3050, 0.3150, 0.3270,
  0.3320, 0.3400, 0.3480, 0.3500,
  0.3530,	0.3600, 0.3680, 0.3650,
  0.3640, 0.3640, 0.3560, 0.3420,

  0.3280, 0.3180, 0.3050, 0.2980,
  0.2920,	0.2900, 0.2920, 0.2945,
  0.3100,	0.3280, 0.3560, 0.3880,
  0.3920,	0.3900, 0.4040, 0.3930,
  0.3700,	0.3500, 0.3300, 0.3150,
  0.2950,	0.2750, 0.2550, 0.2450,
  0.2400 };

  G4double REFLECTIVITY[entries] =
    { 0., 0., 0., 0., 0.,
    0., 0., 0., 0., 0.,
    0., 0., 0., 0., 0.,
    0., 0., 0., 0., 0.,
    0., 0., 0., 0., 0.,
    0., 0., 0., 0., 0. };

  G4MaterialPropertiesTable* phcath_mpt = new G4MaterialPropertiesTable();
  phcath_mpt->AddProperty("EFFICIENCY", ENERGIES, EFFICIENCY, entries);
  phcath_mpt->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY, entries);
  return phcath_mpt;
  // G4OpticalSurface* opt_surf =
  //   new G4OpticalSurface("PHOTOCATHODE", unified, polished, dielectric_metal);
  // opt_surf->SetMaterialPropertiesTable(phcath_mpt);

  // return opt_surf;
}