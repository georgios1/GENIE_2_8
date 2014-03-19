//____________________________________________________________________________
/*
 Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory - October 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 07, 2009 - CA
   Removed call to AddTargetNucleusRemnant(). This simulation step is now
   performed further upstream in the processing chain.  
 @ Mar 03, 2009 - CA
   Moved into the new QEL package from its previous location (EVGModules)
 @ Mar 19, 2013 GC
   Add the coherent length 

*/
//____________________________________________________________________________

#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "QEL/QELHadronicSystemGenerator.h"
#include "Algorithm/AlgFactory.h"
#include "Algorithm/AlgConfigPool.h"

#include "Utils/PhysUtils.h"
#include "Utils/KineUtils.h"
#include "Interaction/Kinematics.h"

using namespace genie;
using namespace genie::utils;
//___________________________________________________________________________
QELHadronicSystemGenerator::QELHadronicSystemGenerator() :
HadronicSystemGenerator("genie::QELHadronicSystemGenerator")
{

}
//___________________________________________________________________________
QELHadronicSystemGenerator::QELHadronicSystemGenerator(string config) :
HadronicSystemGenerator("genie::QELHadronicSystemGenerator", config)
{

}
//___________________________________________________________________________
QELHadronicSystemGenerator::~QELHadronicSystemGenerator()
{

}
//___________________________________________________________________________
void QELHadronicSystemGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state hadronic system

  // Add the recoil baryon 
  // (p or n - Lambda_c+,Sigma_c+,Sigma_c++ in charm/QEL)
  // Its 4-momentum is computed by requiring the energy + momentum to be
  // conserved.
  this->AddRecoilBaryon(evrec);

  //-- Simulate the coh length
  this->SimulateCohLength(evrec);
}
//___________________________________________________________________________
void QELHadronicSystemGenerator::AddRecoilBaryon(GHepRecord * evrec) const
{
  //-- Determine the pdg & status code of the recoil baryon
  Interaction * interaction = evrec->Summary();
  const XclsTag & xcls = interaction->ExclTag();
  int pdgc = 0;
  if(xcls.IsCharmEvent()) { pdgc = xcls.CharmHadronPdg();           }
  else                    { pdgc = interaction->RecoilNucleonPdg(); }
  assert(pdgc!=0);

  //-- Determine the status code
  const Target & tgt = interaction->InitState().Tgt();
  GHepStatus_t ist = (tgt.IsNucleus()) ?
                          kIStHadronInTheNucleus : kIStStableFinalState;

  //-- Get the vtx position
  GHepParticle * neutrino  = evrec->Probe();
  const TLorentzVector & vtx = *(neutrino->X4());

  //-- Get nucleon 4-momentum (in the LAB frame) & position
  TLorentzVector p4 = this->Hadronic4pLAB(evrec);

  //-- Get mother position
  int mom = evrec->HitNucleonPosition();

  //-- Add the final state recoil baryon at the EventRecord
  LOG("QELHadronicVtx", pINFO) 
      << "Adding recoil baryon [pdgc = " << pdgc << "]";

  GHepParticle p(pdgc, ist, mom,-1,-1,-1, p4, vtx);
  double w = xcls.IsCharmEvent() ? 
                0. : evrec->Particle(mom)->RemovalEnergy();
  p.SetRemovalEnergy(w);
  evrec->AddParticle(p);
}
//___________________________________________________________________________
void QELHadronicSystemGenerator::SimulateCohLength(GHepRecord * evrec) const
{

  // Controls the coh length simulation
  if(!fruncohlength)
    return;

  // Print out some useful info
  LOG("QELHadronicVtx", pDEBUG)
    << "Simulating coh length for the QES hadronic system";

  
  // Not found the nuclear target? - Do not simulate the formation zone
  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("QELHadronicVtx", pDEBUG)
      << "No nuclear target was found - No need to simulate coh length";
    return;
  }

  /*
  // QEL events only
  Interaction * interaction = evrec->Summary();
  const ProcessInfo & proc = interaction->ProcInfo();
  if(!proc.IsQuasiElastic()) {
    LOG("QELHadronicVtx", pINFO) << "Not a QEL event - The QELHadronicVtx exits";
    return;
  }
  */

  // Compute the nuclear radius & how far away a particle is being tracked by
  // the intranuclear hadron transport
  assert(nucltgt && pdg::IsIon(nucltgt->Pdg()));
  double A = nucltgt->A();
  double R = fR0 * TMath::Power(A, 1./3.);
  R *= TMath::Max(fNR,1.); // particle is tracked much further outside the nuclear boundary as the density is non-zero

  // Decay very short living particles
  this->PreHadronTransportDecays(evrec);

  TObjArrayIter piter(evrec);
  GHepParticle *recoil = 0;

  while( (recoil = (GHepParticle *) piter.Next()) ){
    // Hardcoded!!! 
    if(recoil->Status() != 14){continue;}
    
    if (!recoil) {
      LOG("QELHadronicVtx", pDEBUG)
	<< "No nuclear target daughter was found - No need to simulate coh length";
      return;
    }
    
    // Calculation
    const TLorentzVector & p4 = *(recoil->P4());
    TVector3 p3hadr = recoil->P4()->Vect(); // (px,py,pz)
    
    // Get the momentum transfer
    GHepParticle * neutrino = evrec->Probe();
    assert(neutrino);
    GHepParticle * fsl = evrec->FinalStatePrimaryLepton();
    assert(fsl);
    
    const TLorentzVector & k1 = *(neutrino->P4());                     // v 4-p (k1)
    const TLorentzVector & k2 = *(fsl->P4());                          // l 4-p (k2)
    
    TLorentzVector q  = k1-k2;                     // q=k1-k2, 4-p transfer

    // Clculate the coh length, t= E/|p*q| - arXiv:1202.4197
    double fz = phys::QELCohLength(p4, q);

    //-- Apply the coh length step
    double step = fz;
    TVector3 dr = recoil->P4()->Vect().Unit();            // unit vector along its direction
    // double c  = kLightSpeed / (units::fm/units::ns); // c in fm/nsec
    dr.SetMag(step);                                 // spatial step size
    // double dt = step/c;                              // temporal step:
    double dt = 0;
    TLorentzVector dx4(dr,dt);                       // 4-vector step
    TLorentzVector x4new = *(recoil->X4()) + dx4;         // new position
    
    //-- If the coh length was large enough that the particle is now outside
    //   the nucleus make sure that it is not placed further away from the
    //   (max distance particles tracked by intranuclear cascade) + ~2 fm
    double epsilon = 2; // fm
    double r       = x4new.Vect().Mag(); // fm
    double rmax    = R+epsilon;
    if(r > rmax) {
      LOG("QELHadronicVtx", pINFO)
	<< "Particle was stepped too far away (r = " << r << " fm)";
      LOG("QELHadronicVtx", pINFO)
	<< "Placing it ~2 fm away from the furthermost position tracked "
	<< "by intranuclear cascades (r' = " << rmax << " fm)";
      
      double scale = rmax/r;
      x4new *= scale;
    }

    // Apply a scale parameter - 1.0 by default
    x4new *= fscaleparameter; 
    
    recoil->SetPosition(x4new);
  }
}
//___________________________________________________________________________
void QELHadronicSystemGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELHadronicSystemGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void QELHadronicSystemGenerator::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  //-- Get parameters controlling the nuclear sizes
  // fm
  fR0  = fConfig->GetDoubleDef ("R0", gc->GetDouble("NUCL-R0"));
  fNR  = fConfig->GetDoubleDef ("NR", gc->GetDouble("NUCL-NR"));

  fruncohlength   = fConfig->GetBoolDef("RunCohLength",    false);
  fscaleparameter = fConfig->GetDoubleDef ("ScaleParameter", 1.0);
}
//____________________________________________________________________________
