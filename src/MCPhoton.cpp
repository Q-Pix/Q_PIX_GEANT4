// -----------------------------------------------------------------------------
//  MCPhoton.cpp
//
//  Class definition of MCPhoton
//   * Author: Everybody is an author!
//   * Creation date: 7 August 2020
// -----------------------------------------------------------------------------

#include "MCPhoton.h"

// GEANT4 includes
#include "G4VProcess.hh"

//-----------------------------------------------------------------------------
MCPhoton::MCPhoton()
{}

//-----------------------------------------------------------------------------
MCPhoton::~MCPhoton()
{}

//-----------------------------------------------------------------------------
// void MCPhoton::AddPhotonPoint(const PhotonPoint & point)
// {
//     Photon_.push_back(point);
// }

//-----------------------------------------------------------------------------
void MCPhoton::AddPhotonHit(const PhotonHit & hit)
{
    photon_hits_.push_back(hit);
}

//-----------------------------------------------------------------------------
void MCPhoton::AddPhotonHit(const G4Step * step)
{
    // initialize hit
    PhotonHit hit;
    hit.start_ = { step->GetPreStepPoint()->GetPosition().x()/CLHEP::cm,
                   step->GetPreStepPoint()->GetPosition().y()/CLHEP::cm,
                   step->GetPreStepPoint()->GetPosition().z()/CLHEP::cm, };
    hit.end_   = { step->GetPostStepPoint()->GetPosition().x()/CLHEP::cm,
                   step->GetPostStepPoint()->GetPosition().y()/CLHEP::cm,
                   step->GetPostStepPoint()->GetPosition().z()/CLHEP::cm, };
    // hit.energy_deposit_ = step->GetTotalEnergyDeposit() / CLHEP::MeV;
    hit.start_time_ = step->GetPreStepPoint()->GetGlobalTime() / CLHEP::ns;
    hit.end_time_   = step->GetPostStepPoint()->GetGlobalTime() / CLHEP::ns;
    // hit.track_id_   = step->GetTrack()->GetTrackID();
    // hit.pdg_code_   = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
    // hit.length_     = step->GetStepLength() / CLHEP::cm;
    // hit.process_    = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

    // add hit
    this->AddPhotonHit(hit);
}
