#include "PrimaryGeneration.h"

#include <G4VPrimaryGenerator.hh>
#include <G4Event.hh>


#include "MuonGenerator.h"

PrimaryGeneration::PrimaryGeneration():
  G4VUserPrimaryGeneratorAction(), generator_(0)
{
  msg_ = new G4GenericMessenger(this, "/Generator/");
  msg_->DeclareProperty("RegisterGenerator", name_, "");

}



PrimaryGeneration::~PrimaryGeneration()
{
  delete msg_;
}



void PrimaryGeneration::GeneratePrimaries(G4Event* event)
{
  G4VPrimaryGenerator* generator_ = 0;
  if      (name_ == "MUON")            generator_ = new MuonGenerator();

  if (!generator_)
    G4Exception("[PrimaryGeneration]", "GeneratePrimaries()",
                FatalException, "Generator not set!");

  generator_->GeneratePrimaryVertex(event);
}






// // -----------------------------------------------------------------------------
// //  G4Basic | PrimaryGeneration.cpp
// //
// //  Class for the definition of the primary generation action.
// //   * Author: Justo Martin-Albo
// //   * Creation date: 14 Aug 2019
// // -----------------------------------------------------------------------------

// #include "PrimaryGeneration.h"

// // Q-Pix includes
// #include "MARLEYManager.h"
// #include "MCTruthManager.h"

// // MARLEY includes
// #include "marley/Event.hh"
// #include "marley/Particle.hh"

// // GEANT4 includes
// #include "G4PhysicalConstants.hh"
// #include "G4LogicalVolumeStore.hh"

// #include "G4ParticleDefinition.hh"
// #include "G4SystemOfUnits.hh"
// #include "G4IonTable.hh"
// #include "G4PrimaryParticle.hh"
// #include "G4PrimaryVertex.hh"
// #include "G4Event.hh"
// #include "G4Electron.hh"
// #include "G4MuonPlus.hh"
// #include "G4Proton.hh"


// #include "G4GenericMessenger.hh"
// #include "Randomize.hh"


// #include "G4ParticleGun.hh"
// #include "G4ParticleTable.hh"
// #include "globals.hh"



// // C++ includes
// #include <stdlib.h>
// #include <math.h>


// PrimaryGeneration::PrimaryGeneration():
//   G4VUserPrimaryGeneratorAction()
// {
//   msg_ = new G4GenericMessenger(this, "/Inputs/", "Control commands of the ion primary generator.");

//   msg_->DeclareProperty("Event_Type", Event_Type_,  "which particle?");

//   msg_->DeclareProperty("Particle_Type", Particle_Type_,  "which particle?");

//   //msg_->DeclareProperty("Particle_energy", Particle_Energy_,  "Energy of the particle.");

//   // get dictionary of particles
//   particle_table_ = G4ParticleTable::GetParticleTable();
// }


// PrimaryGeneration::~PrimaryGeneration()
// {
//   delete msg_;
// }


// void PrimaryGeneration::GeneratePrimaries(G4Event* event)
// {
//   // get detector dimensions
//   if (!detector_solid_vol_)
//   {
//     G4LogicalVolume* detector_logic_vol
//       = G4LogicalVolumeStore::GetInstance()->GetVolume("detector.logical");
//     if (detector_logic_vol) detector_solid_vol_ = dynamic_cast<G4Box*>(detector_logic_vol->GetSolid());
//   }
//   if (detector_solid_vol_)
//   {
//     detector_length_x_ = detector_solid_vol_->GetXHalfLength() * 2.;
//     detector_length_y_ = detector_solid_vol_->GetYHalfLength() * 2.;
//     detector_length_z_ = detector_solid_vol_->GetZHalfLength() * 2.;
//     // G4cout << "det. dim.: " << detector_length_x_ << " m × "
//     //                         << detector_length_y_ << " m × "
//     //                         << detector_length_z_ << " m"
//     //        << G4endl;
//   }

//   MCTruthManager * mc_truth_manager = MCTruthManager::Instance();

//   if ( Event_Type_ = "Particle");
//   {
//     G4int nParticles = 1;
//     particleGun = new G4ParticleGun(nParticles);
//     G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
//     G4ParticleDefinition* particle = particleTable->FindParticle(Particle_Type_);

//     particleGun->SetParticleDefinition(particle);
//     particleGun->SetParticleEnergy(50*keV);
//     particleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
//     // std::cout<<"HadronPrimaryGeneratorAction: setting position "<<std::endl;
//     // G4double zVertex = -5.*CLHEP::mm;
//     particleGun->SetParticlePosition(G4ThreeVector(detector_length_x_/2 , detector_length_y_/2 , detector_length_z_/2));

//     particleGun->GeneratePrimaryVertex(event);

//     for (int i; i<20; i++){G4cout<<"AUSTIN_HERE"<<G4endl;}
    

//   }





//   // if (Particle_Type_ ==  "Ar39")
//   // {
//   //   G4ParticleDefinition* pdef = G4IonTable::GetIonTable()->GetIon(18, 39, 0.); // Ar39
//   //   if (!pdef)G4Exception("SetParticleDefinition()", "[IonGun]",FatalException, " can not create ion ");

//   //   G4PrimaryParticle* particle = new G4PrimaryParticle(pdef);
//   //   particle->SetMomentumDirection(G4ThreeVector(0.,1.,0.));
//   //   particle->SetKineticEnergy(1.*eV); // just an ion sitting

//   //   double Ran_X = G4UniformRand() * detector_length_x_/2.;
//   //   double Ran_Y = G4UniformRand() * detector_length_y_/2.;
//   //   double Ran_Z = G4UniformRand() * detector_length_z_/2.;
//   //   if (G4UniformRand()>0.5){ Ran_X *= -1; }
//   //   if (G4UniformRand()>0.5){ Ran_Y *= -1; }
//   //   if (G4UniformRand()>0.5){ Ran_Z *= -1; }
//   //   G4PrimaryVertex* vertex = new G4PrimaryVertex(G4ThreeVector(Ran_X,Ran_Y,Ran_Z), 0.);
//   //   vertex->SetPrimary(particle);
//   //   event->AddPrimaryVertex(vertex);
//   // }
//   // else if (Particle_Type_ ==  "Muon")
//   // {
//   //   double Ran_X = G4UniformRand() * detector_length_x_/2.;
//   //   double Ran_Y = G4UniformRand() * detector_length_y_/2.;
//   //   double Ran_Z = G4UniformRand() * detector_length_z_/2.;
//   //   if (G4UniformRand()>0.5){ Ran_X *= -1; }
//   //   if (G4UniformRand()>0.5){ Ran_Y *= -1; }
//   //   if (G4UniformRand()>0.5){ Ran_Z *= -1; }
//   //   G4ParticleDefinition* pdef = G4MuonPlus::Definition();

//   //   G4PrimaryParticle* particle = new G4PrimaryParticle(pdef);
//   //   particle->SetMomentumDirection(G4ThreeVector(0.,1.,0.));
//   //   particle->SetKineticEnergy(10.*MeV); 


//   //   G4PrimaryVertex* vertex = new G4PrimaryVertex(G4ThreeVector(0, -detector_length_y_/2., Ran_Z), 0.);
//   //   vertex->SetPrimary(particle);
//   //   event->AddPrimaryVertex(vertex);



//   // }
//   // else if (Particle_Type_ ==  "Proton")
//   // {

//   //   double Ran_X = G4UniformRand() * detector_length_x_/2.;
//   //   double Ran_Y = G4UniformRand() * detector_length_y_/2.;
//   //   double Ran_Z = G4UniformRand() * detector_length_z_/2.;
//   //   if (G4UniformRand()>0.5){ Ran_X *= -1; }
//   //   if (G4UniformRand()>0.5){ Ran_Y *= -1; }
//   //   if (G4UniformRand()>0.5){ Ran_Z *= -1; }
//   //   double NorXY = sqrt(Ran_X*Ran_X+Ran_Y*Ran_Y) ;

//   //   G4ParticleDefinition* pdef = G4Proton::Definition();

//   //   G4PrimaryParticle* particle = new G4PrimaryParticle(pdef);
//   //   particle->SetMomentumDirection(G4ThreeVector(Ran_X/NorXY , Ran_Y/NorXY, 0.));
//   //   particle->SetKineticEnergy(100.*MeV); // just an ion sitting

//   //   G4PrimaryVertex* vertex = new G4PrimaryVertex(G4ThreeVector(Ran_X,Ran_Y,Ran_Z), 0.);
//   //   vertex->SetPrimary(particle);
//   //   event->AddPrimaryVertex(vertex);
//   // }
//   // else if (Particle_Type_ ==  "Electron")
//   // {
//   //   G4ParticleDefinition* pdef = G4Electron::Definition();

//   //   G4PrimaryParticle* particle = new G4PrimaryParticle(pdef);
//   //   particle->SetMomentumDirection(G4ThreeVector(0.,1.,0.));
//   //   particle->SetKineticEnergy(3.*MeV); 


//   //   G4PrimaryVertex* vertex = new G4PrimaryVertex(G4ThreeVector(0, 0 ,0.5), 0.);
//   //   vertex->SetPrimary(particle);
//   //   event->AddPrimaryVertex(vertex);
//   // }

//   // else if (Particle_Type_ ==  "MARLEY")
//   // {
//   //   // get MARLEY manager and generator
//   //   MARLEYManager * marley_manager = MARLEYManager::Instance();
//   //   marley::Generator & marley_generator = marley_manager->Generator();

//   //   // G4PrimaryVertex* vertex = new G4PrimaryVertex(G4ThreeVector(0., 0., 0.), 0.);

//   //   G4ThreeVector offset(detector_length_x_/2.,
//   //                        detector_length_y_/2.,
//   //                        detector_length_z_/2.);

//   //   G4PrimaryVertex* vertex = new G4PrimaryVertex(offset, 0.);

//   //   // Generate a new MARLEY event using the owned marley::Generator object
//   //   marley::Event ev = marley_generator.create_event();

//   //   // // print MARLEY event information
//   //   // ev.print_human_readable(G4cout);

//   //   // Loop over each of the initial particles in the MARLEY event
//   //   for (const auto& ip : ev.get_initial_particles())
//   //   {
//   //     // add initial MARLEY particle to the MC truth manager
//   //     mc_truth_manager->AddInitialMARLEYParticle(*ip);
//   //   }

//   //   // Loop over each of the final particles in the MARLEY event
//   //   for (const auto& fp : ev.get_final_particles())
//   //   {
//   //     // add final MARLEY particle to the MC truth manager
//   //     mc_truth_manager->AddFinalMARLEYParticle(*fp);

//   //     // Convert each one from a marley::Particle into a G4PrimaryParticle.
//   //     // Do this by first setting the PDG code and the 4-momentum components.

//   //     // get dictionary of particles if necessary
//   //     if (particle_table_ == 0)
//   //     {
//   //       particle_table_ = G4ParticleTable::GetParticleTable();
//   //     }

//   //     // initialize particle definition
//   //     G4ParticleDefinition* pdef;

//   //     // get PDG code of marley::Particle
//   //     int const pdg_code = fp->pdg_code();

//   //     if (pdg_code == 0)
//   //     {
//   //       pdef = particle_table_->FindParticle("opticalphoton");
//   //     }
//   //     else
//   //     {
//   //       pdef = particle_table_->FindParticle(pdg_code);
//   //     }

//   //     // if the particle is a nucleus
//   //     if (pdg_code > 1000000000)
//   //     {
//   //       if (!pdef)
//   //       {
//   //         int const Z = (pdg_code % 10000000) / 10000; // atomic number
//   //         int const A = (pdg_code % 10000) / 10; // mass number
//   //         pdef = particle_table_->GetIonTable()->GetIon(Z, A, 0.);
//   //       }
//   //     } // if the particle is a nucleus

//   //     if (pdef == 0)
//   //     {
//   //       std::string message = "\nLine "
//   //                           + std::to_string(__LINE__)
//   //                           + " of file "
//   //                           + __FILE__
//   //                           + "\n\nUnknown PDG code: "
//   //                           + std::to_string(pdg_code)
//   //                           + "\n";
//   //       G4Exception("PrimaryGeneration::GeneratePrimaries", "Error",
//   //                   FatalException, message.c_str());
//   //     }

//   //     G4PrimaryParticle* particle = new G4PrimaryParticle(pdef,
//   //                                                         fp->px(),
//   //                                                         fp->py(),
//   //                                                         fp->pz());

//   //     // Also set the charge of the G4PrimaryParticle appropriately
//   //     particle->SetCharge( fp->charge() );

//   //     // particle->SetPolarization(); ??

//   //     // Add the fully-initialized G4PrimaryParticle to the primary vertex
//   //     vertex->SetPrimary( particle );
//   //   }

//   //   // The primary vertex has been fully populated with all final-state particles
//   //   // from the MARLEY event. Add it to the G4Event object so that Geant4 can
//   //   // begin tracking the particles through the simulated geometry.
//   //   event->AddPrimaryVertex( vertex );

//   // }

//   // else
//   // {
//   //   exit (EXIT_FAILURE);
//   //   //G4Exception(FatalException, " Pick a defined particle... ");
//   // }


// }
