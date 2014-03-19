//____________________________________________________________________________
/*!

\class    genie::QELHadronicSystemGenerator

\brief    Generates the final state hadronic system in v QEL interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 03, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _QEL_HADRONIC_SYSTEM_GENERATOR_H_
#define _QEL_HADRONIC_SYSTEM_GENERATOR_H_

#include "EVGModules/HadronicSystemGenerator.h"

namespace genie {

class QELHadronicSystemGenerator : public HadronicSystemGenerator {

public :
  QELHadronicSystemGenerator();
  QELHadronicSystemGenerator(string config);
 ~QELHadronicSystemGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event_rec) const;

  // Simulate the formation zone
  void SimulateCohLength(GHepRecord * evrec) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void AddRecoilBaryon    (GHepRecord * event_rec) const;

  double fR0;             ///< param controling nuclear size
  double fNR;             ///< how far beyond the nuclear boundary does the particle tracker goes?
  bool fruncohlength;     ///< run the coh length simulation 
  double fscaleparameter; ///< scale parameter 

  void LoadConfig (void);
};

}      // genie namespace
#endif // _QEL_HADRONIC_SYSTEM_GENERATOR_H_
