// ----------------------------------------------------------------------------
// nexus | PmtSD.cc
//
// This class is the sensitive detector that allows for the registration
// of the charge detected by a photosensor.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "PmtSD.h"
#include "MCParticle.h"
#include "MCTruthManager.h"


#include <G4OpticalPhoton.hh>
#include <G4SDManager.hh>
#include <G4ProcessManager.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4RunManager.hh>
#include <G4RunManager.hh>


PmtSD::PmtSD(G4String sdname, const G4String& hc_name):
G4VSensitiveDetector(sdname),
naming_order_(0), sensor_depth_(0), mother_depth_(0),
boundary_(0)
{
    // Register the name of the collection of hits
    collectionName.insert(GetCollectionUniqueName());
    collectionName.insert(hc_name);
}



PmtSD::~PmtSD()
{
}



G4String PmtSD::GetCollectionUniqueName()
{
    return "PmtHitsCollection";
}



// void PmtSD::Initialize(G4HCofThisEvent* HCE)
// {
//     // for (int i=0; i<100; ++i){G4cout << "WHOOOO" << G4endl;}
//     // Create a new collection of PMT hits
//     HC_ = new PmtHitsCollection(this->GetName(), this->GetCollectionName(0));

//     G4int HCID = G4SDManager::GetSDMpointer()->
//         GetCollectionID(this->GetName()+"/"+this->GetCollectionName(0));

//     HCE->AddHitsCollection(HCID, HC_);
// }



G4bool PmtSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    // for (int i=0; i<10; ++i){G4cout << "FKN WORK" << G4endl;} 
    // G4cout << "FKN WORK" << G4endl;

    // Check whether the track is an optical photon
    G4ParticleDefinition* pdef = step->GetTrack()->GetDefinition();
    // G4cout << pdef << G4endl;
    if (pdef != G4OpticalPhoton::Definition()) return false;
    // G4cout << "I is photon" << G4endl;

    // Retrieve the pointer to the optical boundary process, if it has
    // not been done yet (i.e., if the pointer is not defined)
    if (!boundary_) 
    {
        // Get the list of processes defined for the optical photon
        // and loop through it to find the optical boundary process.
        G4ProcessVector* pv = pdef->GetProcessManager()->GetProcessList();
        for (size_t i=0; i<pv->size(); i++) 
        {
            if ((*pv)[i]->GetProcessName() == "OpBoundary") 
            {
                boundary_ = (G4OpBoundaryProcess*) (*pv)[i];
                break;
            }       
        }
    }


    // Check if the photon has reached a geometry boundary
    if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) 
    {

        // Check whether the photon has been detected in the boundary
        if (boundary_->GetStatus() == Detection) 
        {
            const G4VTouchable* touchable =
                step->GetPostStepPoint()->GetTouchable();

            // G4int pmt_id = FindPmtID(touchable);

            // PmtHit* hit = 0;
            // for (size_t i=0; i<HC_->entries(); i++) 
            // {
            //     if ((*HC_)[i]->GetPmtID() == pmt_id) 
            //     {
            //         hit = (*HC_)[i];
            //         break;
            //     }
            // }

            // // If no hit associated to this sensor exists already,
            // // create it and set main properties
            // if (!hit) 
            // {
            //     hit = new PmtHit();
            //     hit->SetPmtID(pmt_id);
            //     hit->SetBinSize(timebinning_);
            //     hit->SetPosition(touchable->GetTranslation());
            //     HC_->insert(hit);
            // }

            //---------------------------------------------------------------------------
            // begin add hit to MCParticle
            //---------------------------------------------------------------------------

            // get MC truth manager
            MCTruthManager * mc_truth_manager = MCTruthManager::Instance();

            // get MC particle
            MCParticle * particle = mc_truth_manager->GetMCParticle(step->GetTrack()->GetTrackID());

            // add hit to MC particle
            particle->AddPmtHit(step);



            // G4double time = step->GetPreStepPoint()->GetGlobalTime()/CLHEP::ns;
            // G4double MyX  = step->GetPreStepPoint()->GetPosition().x()/CLHEP::cm;
            // G4double MyY  = step->GetPreStepPoint()->GetPosition().y()/CLHEP::cm;
            // G4double MyZ  = step->GetPreStepPoint()->GetPosition().z()/CLHEP::cm;
            // G4cout << "HIT TIME, " << time <<"\t"<< MyX <<"\t"<< MyY <<"\t"<< MyZ << G4endl;

            // time = step->GetPostStepPoint()->GetGlobalTime()/CLHEP::ns;
            // MyX  = step->GetPostStepPoint()->GetPosition().x()/CLHEP::cm;
            // MyY  = step->GetPostStepPoint()->GetPosition().y()/CLHEP::cm;
            // MyZ  = step->GetPostStepPoint()->GetPosition().z()/CLHEP::cm;
            // G4cout << "HIT TIME, " << time <<"\t"<< MyX <<"\t"<< MyY <<"\t"<< MyZ << G4endl;

            // hit->Fill(time);
        }
    }

    return true;
}



G4int PmtSD::FindPmtID(const G4VTouchable* touchable)
{
    G4int pmtid = touchable->GetCopyNumber(sensor_depth_);
    if (naming_order_ != 0) {
        G4int motherid = touchable->GetCopyNumber(mother_depth_);
        pmtid = naming_order_ * motherid + pmtid;
    }
    return pmtid;
}


void PmtSD::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{
//  int HCID = G4SDManager::GetSDMpointer()->
//    GetCollectionID(this->GetCollectionName(0));
//  // }
// HCE->AddHitsCollection(HCID, HC_);


}

