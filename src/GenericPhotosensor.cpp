// -----------------------------------------------------------------------------
//  nexus | GenericPhotosensor.cc
//
//  Geometry of a configurable box-shaped photosensor.
//
//  The NEXT Collaboration
// -----------------------------------------------------------------------------

#include "GenericPhotosensor.h"

#include "MaterialsList.h"
#include "PmtSD.h"
#include "OpticalMaterialProperties.h"
#include "Visibilities.h"

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4SDManager.hh>
#include <G4GenericMessenger.hh>



GenericPhotosensor::GenericPhotosensor(G4String name,
                                       G4double width,
                                       G4double height,
                                       G4double thickness):
  BaseGeometry        (),
  name_               (name),
  width_              (width),           // Width of the Sensitive Area
  height_             (height),          // Height of the Sensitive Area
  thickness_          (thickness),       // Thickness of the whole sensor
  window_thickness_   (0.2 * mm),        // Window thickness    (similar to Sensl SiPMs)
  sensarea_thickness_ (0.2 * mm),        // Sensitive thickness (similar to Sensl SiPMs)
  wls_thickness_      (1. * micrometer), // WLS thickness = 1 micron by definition)
  with_wls_coating_   (false),
  window_rindex_      (nullptr),
  sensitive_mpt_      (nullptr),
  sensor_depth_       (-1),
  mother_depth_       (0),
  naming_order_       (0),
  time_binning_       (1.0 * us),
  visibility_         (false)
{
}


GenericPhotosensor::GenericPhotosensor(G4String name, G4double size):
  GenericPhotosensor(name, size, size)
{
}


GenericPhotosensor::~GenericPhotosensor()
{
}


void GenericPhotosensor::ComputeDimensions()
{
  // Reduced size for components inside the case except the WLS
  reduced_width_  = width_  - 1. * micrometer;
  reduced_height_ = height_ - 1. * micrometer;

  if (!with_wls_coating_) wls_thickness_ = 0.;

  // Check that components (window + sensitive + wls) fits into the case
  if ((window_thickness_ + sensarea_thickness_ + wls_thickness_) > thickness_) {
    G4Exception("[GenericPhotosensor]", "ComputeDimensions()", FatalException,
                ("Sensor size too small. Required thickness >= " +
                 std::to_string(window_thickness_ + sensarea_thickness_ + wls_thickness_) + 
                 " mm").data());
  }
}


void GenericPhotosensor::DefineMaterials()
{
  // Case /////
  case_mat_ = MaterialsList::FR4();

  // Window /////
  // The optical properties of a given material, is common for the whole geometry so,
  // making a copy of the window material restricts its use to this photosensor
  // and prevents interferences with other possible uses.
  window_mat_ =
    MaterialsList::CopyMaterial(MaterialsList::OpticalSilicone(),
                                name_ + "_WINDOW_MATERIAL");
  G4MaterialPropertiesTable* window_optProp = new G4MaterialPropertiesTable();

  // In the default behavior of this class, the refractive index of WLS and window
  // are matched to avoid reflection losses in the interfaces. If the user sets
  // explicitly a refractive index for the window, it won't be used, raising a warning.
  if (with_wls_coating_) {
    if (window_rindex_)
      G4Exception("[GenericPhotosensor]", "DefineMaterials()", JustWarning,
                  "Window rindex set, but NOT USED. Using TPB rindex.");

    window_optProp->AddProperty("RINDEX",
                                OpticalMaterialProperties::TPB()->GetProperty("RINDEX"));
    window_mat_->SetMaterialPropertiesTable(window_optProp);
  }

  // If the sensor has NOT WLS coating the window must have the rindex
  // set by the user. If it is not set, an exception raises.
  else {
    if (!window_rindex_)
      G4Exception("[GenericPhotosensor]", "DefineMaterials()", FatalException,
                  "Window rindex must be set before constructing");

    window_optProp->AddProperty("RINDEX", window_rindex_);
    window_mat_->SetMaterialPropertiesTable(window_optProp);
  }

  // Sensitive /////
  sensitive_mat_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

  // WLS coating /////
  wls_mat_ = MaterialsList::TPB();
  wls_mat_->SetMaterialPropertiesTable(OpticalMaterialProperties::TPB());
}


void GenericPhotosensor::Construct()
{
  ComputeDimensions();

  DefineMaterials();

  // PHOTOSENSOR CASE //////////////////////////////////
  G4String name = name_ + "_CASE";

  G4Box* case_solid_vol =
    new G4Box(name, width_/2., height_/2., thickness_/2.);

  G4LogicalVolume* case_logic_vol =
    new G4LogicalVolume(case_solid_vol, case_mat_, name);

  BaseGeometry::SetLogicalVolume(case_logic_vol);


  // OPTICAL WINDOW ////////////////////////////////////////
  name = name_ + "_WINDOW";

  G4Box* window_solid_vol =
    new G4Box(name, reduced_width_/2., reduced_height_/2., window_thickness_/2.);

  G4LogicalVolume* window_logic_vol =
    new G4LogicalVolume(window_solid_vol, window_mat_, name);

  G4double window_zpos = thickness_/2. - wls_thickness_ - window_thickness_/2.;

  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., window_zpos), window_logic_vol,
                    name, case_logic_vol, false, 0, false);


  // PHOTOSENSITIVE AREA /////////////////////////////////////////////
  name = name_ + "_SENSAREA";

  G4Box* sensarea_solid_vol =
    new G4Box(name, reduced_width_/2., reduced_height_/2., sensarea_thickness_/2.);

  G4LogicalVolume* sensarea_logic_vol =
    new G4LogicalVolume(sensarea_solid_vol, sensitive_mat_, name);

  G4double sensarea_zpos = thickness_/2. - wls_thickness_ - window_thickness_ -
                           sensarea_thickness_/2.;

  new G4PVPlacement(nullptr, G4ThreeVector(0., 0., sensarea_zpos), sensarea_logic_vol,
                    name, case_logic_vol, false, 0, false);


  // WLS_COATING /////////////////////////////////////////////
  G4LogicalVolume* wls_logic_vol;
  if (with_wls_coating_) {
    name = name_ + "_WLS";

    G4Box* wls_solid_vol =
      new G4Box(name, width_/2., height_/2., wls_thickness_/2.);

    wls_logic_vol = new G4LogicalVolume(wls_solid_vol, wls_mat_, name);

    G4double wls_zpos = thickness_/2. - wls_thickness_/2.;

    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., wls_zpos), wls_logic_vol,
                      name, case_logic_vol, false, 0, false);

    G4OpticalSurface* wls_optSurf = new G4OpticalSurface(name + "_OPSURF",
                                                         glisur, ground,
                                                         dielectric_dielectric, .01);

    new G4LogicalSkinSurface(name + "_OPSURF", wls_logic_vol, wls_optSurf);

    wls_logic_vol->SetVisAttributes(G4VisAttributes::Invisible);
  }


  // VISIBILITIES /////////////////////////////////////////////
  if (visibility_) {
    window_logic_vol  ->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes red = Red();
    red.SetForceSolid(true);
    sensarea_logic_vol->SetVisAttributes(red);
    case_logic_vol->SetVisAttributes(G4VisAttributes::Invisible);
  }
  else {
    window_logic_vol  ->SetVisAttributes(G4VisAttributes::Invisible);
    sensarea_logic_vol->SetVisAttributes(G4VisAttributes::Invisible);
    case_logic_vol->SetVisAttributes(G4VisAttributes::Invisible);
  }


  // SENSOR OPTICAL PROPERTIES ////////////////////////////////////////
  G4MaterialPropertiesTable* sensitive_mpt_ = GetMyPhotOptSurf();
  if (!sensitive_mpt_)
    G4Exception("[GenericPhotosensor]", "Construct()", FatalException,
                "Sensor Optical Properties must be set before constructing");
  G4OpticalSurface* sensitive_opsurf =
    new G4OpticalSurface(name + "_optSurf", unified, polished, dielectric_metal);
  sensitive_opsurf->SetMaterialPropertiesTable(sensitive_mpt_);
  new G4LogicalSkinSurface(name + "_optSurf", sensarea_logic_vol, sensitive_opsurf);


  // SENSITIVE DETECTOR //////////////////////////////////////////////
  G4String sdname = "/GENERIC_PHOTOSENSOR/" + name_;
  G4SDManager* sdmgr = G4SDManager::GetSDMpointer();

  if (!sdmgr->FindSensitiveDetector(sdname, false)) {
    PmtSD* sensdet = new PmtSD(sdname, "PmtHitsCollection");
    if (sensor_depth_ == -1)
      G4Exception("[GenericPhotosensor]", "Construct()", FatalException,
                  "Sensor Depth must be set before constructing");

    if ((naming_order_ > 0) & (mother_depth_ == 0))
      G4Exception("[GenericPhotosensor]", "Construct()", FatalException,
                  "naming_order set without setting mother_depth");

    sensdet->SetDetectorVolumeDepth(sensor_depth_);
    sensdet->SetMotherVolumeDepth  (mother_depth_);
    sensdet->SetDetectorNamingOrder(naming_order_);
    sensdet->SetTimeBinning        (time_binning_);

    G4SDManager::GetSDMpointer()->AddNewDetector(sensdet);
    sensarea_logic_vol->SetSensitiveDetector(sensdet);
  }
}


// G4OpticalSurface* GenericPhotosensor::GetMyPhotOptSurf()
G4MaterialPropertiesTable* GenericPhotosensor::GetMyPhotOptSurf()
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