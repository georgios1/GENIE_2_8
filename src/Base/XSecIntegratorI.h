//____________________________________________________________________________
/*!

\class    genie::XSecIntegratorI

\brief    Cross Section Integrator Interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  May 03, 2004

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _XSEC_INTEGRATOR_I_H_
#define _XSEC_INTEGRATOR_I_H_

#include "Algorithm/Algorithm.h"
#include "Base/XSecAlgorithmI.h"
#include "Interaction/Interaction.h"

namespace genie {

class IntegratorI;
class XSecIntegratorI : public Algorithm {

public:
  virtual ~XSecIntegratorI();

  virtual double Integrate(const XSecAlgorithmI * model, 
                           const Interaction * interaction 
                       /*, const KPhaseSpaceCut * cut=0*/) const= 0;
protected:
  XSecIntegratorI();
  XSecIntegratorI(string name);
  XSecIntegratorI(string name, string config);

  const IntegratorI * fIntegrator; ///< GENIE numerical integrator 

  string fGSLIntgType; ///< name of GSL numerical integrator
  double fGSLRelTol;   ///< required relative tolerance (error)
};

}       // genie namespace
#endif  // _XSEC_INTEGRATOR_I_H_
