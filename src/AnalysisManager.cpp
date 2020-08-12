// -----------------------------------------------------------------------------
//  AnalysisManager.cpp
//
//  Class definition of the analysis manager
//   * Author: Everybody is an author!
//   * Creation date: 4 August 2020
// -----------------------------------------------------------------------------

#include "AnalysisManager.h"

AnalysisManager * AnalysisManager::instance_ = 0;

//-----------------------------------------------------------------------------
AnalysisManager::AnalysisManager()
{
#ifdef G4ANALYSIS_USE
#endif
}

//-----------------------------------------------------------------------------
AnalysisManager::~AnalysisManager()
{
#ifdef G4ANALYSIS_USE
#endif
}

//-----------------------------------------------------------------------------
AnalysisManager * AnalysisManager::Instance()
{
    if (instance_ == 0) instance_ = new AnalysisManager();
    return instance_;
}

//-----------------------------------------------------------------------------
void AnalysisManager::Book(std::string const file_path)
{
    // ROOT output file
    tfile_ = new TFile(file_path.data(), "recreate", "q_pix_geant");

    // GEANT4 event tree
    ttree_ = new TTree("g4_event_tree", "GEANT4 event tree");

    ttree_->Branch("run",   &run_,   "run/I");
    ttree_->Branch("event", &event_, "event/I");

    ttree_->Branch("energy_deposit", &energy_deposit_, "energy_deposit/D");

    ttree_->Branch("particle_track_id",        &particle_track_id_);
    ttree_->Branch("particle_parent_track_id", &particle_parent_track_id_);
    ttree_->Branch("particle_pdg_code",        &particle_pdg_code_);
    ttree_->Branch("particle_mass",            &particle_mass_);
    ttree_->Branch("particle_charge",          &particle_charge_);
    ttree_->Branch("particle_process_key",     &particle_process_key_);
    ttree_->Branch("particle_total_occupancy", &particle_total_occupancy_);

    ttree_->Branch("particle_initial_x", &particle_initial_x_);
    ttree_->Branch("particle_initial_y", &particle_initial_y_);
    ttree_->Branch("particle_initial_z", &particle_initial_z_);
    ttree_->Branch("particle_initial_t", &particle_initial_t_);

    ttree_->Branch("particle_initial_px", &particle_initial_px_);
    ttree_->Branch("particle_initial_py", &particle_initial_py_);
    ttree_->Branch("particle_initial_pz", &particle_initial_pz_);
    ttree_->Branch("particle_initial_energy", &particle_initial_energy_);

    ttree_->Branch("hit_track_id",       &hit_track_id_);

    ttree_->Branch("hit_start_x",        &hit_start_x_);
    ttree_->Branch("hit_start_y",        &hit_start_y_);
    ttree_->Branch("hit_start_z",        &hit_start_z_);
    ttree_->Branch("hit_start_t",        &hit_start_t_);

    ttree_->Branch("hit_end_x",          &hit_end_x_);
    ttree_->Branch("hit_end_y",          &hit_end_y_);
    ttree_->Branch("hit_end_z",          &hit_end_z_);
    ttree_->Branch("hit_end_t",          &hit_end_t_);

    ttree_->Branch("hit_energy_deposit", &hit_energy_deposit_);
    ttree_->Branch("hit_length",         &hit_length_);
    ttree_->Branch("hit_process_key",    &hit_process_key_);

    // MARLEY event tree
    marley_event_tree_ = new TTree("marley_event_tree", "MARLEY event tree");

    marley_event_tree_->Branch("run", &run_, "run/I");
    marley_event_tree_->Branch("event", &event_, "event/I");

    marley_event_tree_->Branch("initial_particle_px", &marley_initial_particle_px_);
    marley_event_tree_->Branch("initial_particle_py", &marley_initial_particle_py_);
    marley_event_tree_->Branch("initial_particle_pz", &marley_initial_particle_pz_);
    marley_event_tree_->Branch("initial_particle_energy", &marley_initial_particle_energy_);
    marley_event_tree_->Branch("initial_particle_pdg_code", &marley_initial_particle_pdg_code_);
    marley_event_tree_->Branch("initial_particle_mass", &marley_initial_particle_mass_);
    marley_event_tree_->Branch("initial_particle_charge", &marley_initial_particle_charge_);

    marley_event_tree_->Branch("final_particle_px", &marley_final_particle_px_);
    marley_event_tree_->Branch("final_particle_py", &marley_final_particle_py_);
    marley_event_tree_->Branch("final_particle_pz", &marley_final_particle_pz_);
    marley_event_tree_->Branch("final_particle_energy", &marley_final_particle_energy_);
    marley_event_tree_->Branch("final_particle_pdg_code", &marley_final_particle_pdg_code_);
    marley_event_tree_->Branch("final_particle_mass", &marley_final_particle_mass_);
    marley_event_tree_->Branch("final_particle_charge", &marley_final_particle_charge_);
}

//-----------------------------------------------------------------------------
void AnalysisManager::Save()
{
    // write TTree objects to file and close file
    ttree_->Write();
    marley_event_tree_->Write();
    tfile_->Close();
}

//-----------------------------------------------------------------------------
void AnalysisManager::EventReset()
{
    // reset event variables after filling TTree objects per event
    event_ = -1;
    energy_deposit_ = 0;

    particle_track_id_.clear();
    particle_parent_track_id_.clear();
    particle_pdg_code_.clear();
    particle_mass_.clear();
    particle_charge_.clear();
    particle_process_key_.clear();
    particle_total_occupancy_.clear();

    particle_initial_x_.clear();
    particle_initial_y_.clear();
    particle_initial_z_.clear();
    particle_initial_t_.clear();

    particle_initial_px_.clear();
    particle_initial_py_.clear();
    particle_initial_pz_.clear();
    particle_initial_energy_.clear();

    hit_track_id_.clear();
    hit_start_x_.clear();
    hit_start_y_.clear();
    hit_start_z_.clear();
    hit_start_t_.clear();
    hit_end_x_.clear();
    hit_end_y_.clear();
    hit_end_z_.clear();
    hit_end_t_.clear();
    hit_energy_deposit_.clear();
    hit_length_.clear();
    hit_process_key_.clear();

    marley_initial_particle_px_.clear();
    marley_initial_particle_py_.clear();
    marley_initial_particle_pz_.clear();
    marley_initial_particle_energy_.clear();
    marley_initial_particle_pdg_code_.clear();
    marley_initial_particle_mass_.clear();
    marley_initial_particle_charge_.clear();

    marley_final_particle_px_.clear();
    marley_final_particle_py_.clear();
    marley_final_particle_pz_.clear();
    marley_final_particle_energy_.clear();
    marley_final_particle_pdg_code_.clear();
    marley_final_particle_mass_.clear();
    marley_final_particle_charge_.clear();
}

//-----------------------------------------------------------------------------
void AnalysisManager::EventFill()
{
    // fill TTree objects per event
    ttree_->Fill();
    marley_event_tree_->Fill();
}

//-----------------------------------------------------------------------------
void AnalysisManager::SetRun(int const value)
{
    run_ = value;
}

//-----------------------------------------------------------------------------
void AnalysisManager::SetEvent(int const value)
{
    event_ = value;
}

//-----------------------------------------------------------------------------
void AnalysisManager::AddMCParticle(MCParticle const * particle)
{
    particle_track_id_.push_back(particle->TrackID());
    particle_parent_track_id_.push_back(particle->ParentTrackID());
    particle_pdg_code_.push_back(particle->PDGCode());
    particle_mass_.push_back(particle->Mass());
    particle_charge_.push_back(particle->Charge());
    particle_process_key_.push_back(this->ProcessToKey(particle->Process()));
    particle_total_occupancy_.push_back(particle->TotalOccupancy());

    // inline TLorentzVector InitialPosition() const { return initial_position_; }
    // inline TLorentzVector InitialMomentum() const { return initial_momentum_; }

    particle_initial_x_.push_back(particle->InitialPosition().X());
    particle_initial_y_.push_back(particle->InitialPosition().Y());
    particle_initial_z_.push_back(particle->InitialPosition().Z());
    particle_initial_t_.push_back(particle->InitialPosition().T());

    particle_initial_px_.push_back(particle->InitialMomentum().X());
    particle_initial_py_.push_back(particle->InitialMomentum().Y());
    particle_initial_pz_.push_back(particle->InitialMomentum().Z());
    particle_initial_energy_.push_back(particle->InitialMomentum().E());

    std::vector< TrajectoryHit > const hits = particle->Hits();

    for (auto const & hit : hits)
    {
        energy_deposit_ += hit.Energy();

        hit_track_id_.push_back(hit.TrackID());

        hit_start_x_.push_back(hit.StartPoint().X());
        hit_start_y_.push_back(hit.StartPoint().Y());
        hit_start_z_.push_back(hit.StartPoint().Z());
        hit_start_t_.push_back(hit.StartTime());

        hit_end_x_.push_back(hit.EndPoint().X());
        hit_end_y_.push_back(hit.EndPoint().Y());
        hit_end_z_.push_back(hit.EndPoint().Z());
        hit_end_t_.push_back(hit.EndTime());

        hit_length_.push_back(hit.Length());
        hit_energy_deposit_.push_back(hit.Energy());

        hit_process_key_.push_back(this->ProcessToKey(hit.Process()));
    }
}

//-----------------------------------------------------------------------------
void AnalysisManager::AddInitialMARLEYParticle(marley::Particle const & particle)
{
    marley_initial_particle_px_.push_back(particle.px());
    marley_initial_particle_py_.push_back(particle.py());
    marley_initial_particle_pz_.push_back(particle.pz());
    marley_initial_particle_energy_.push_back(particle.total_energy());
    marley_initial_particle_pdg_code_.push_back(particle.pdg_code());
    marley_initial_particle_mass_.push_back(particle.mass());
    marley_initial_particle_charge_.push_back(particle.charge());
}

//-----------------------------------------------------------------------------
void AnalysisManager::AddFinalMARLEYParticle(marley::Particle const & particle)
{
    marley_final_particle_px_.push_back(particle.px());
    marley_final_particle_py_.push_back(particle.py());
    marley_final_particle_pz_.push_back(particle.pz());
    marley_final_particle_energy_.push_back(particle.total_energy());
    marley_final_particle_pdg_code_.push_back(particle.pdg_code());
    marley_final_particle_mass_.push_back(particle.mass());
    marley_final_particle_charge_.push_back(particle.charge());
}

//-----------------------------------------------------------------------------
int AnalysisManager::ProcessToKey(std::string const & process)
{
    int key = 0;

    if      (process.compare("eIoni")                == 0) key =  1;
    else if (process.compare("msc")                  == 0) key =  2;
    else if (process.compare("compt")                == 0) key =  3;
    else if (process.compare("phot")                 == 0) key =  4;
    else if (process.compare("eBrem")                == 0) key =  5;
    else if (process.compare("ionIoni")              == 0) key =  6;
    else if (process.compare("hIoni")                == 0) key =  7;
    else if (process.compare("RadioactiveDecayBase") == 0) key =  8;
    else if (process.compare("CoulombScat")          == 0) key =  9;
    else if (process.compare("Rayl")                 == 0) key = 10;
    else if (process.compare("Transportation")       == 0) key = 11;
    else if (process.compare("annihil")              == 0) key = 12;
    else if (process.compare("conv")                 == 0) key = 13;
    else if (process.compare("hadElastic")           == 0) key = 14;
    else if (process.compare("nCapture")             == 0) key = 15;
    else if (process.compare("neutronInelastic")     == 0) key = 16;
    else if (process.compare("photonNuclear")        == 0) key = 17;

    return key;
}
