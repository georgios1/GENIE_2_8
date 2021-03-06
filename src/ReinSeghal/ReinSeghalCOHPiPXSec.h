//____________________________________________________________________________
/*!

\class    genie::ReinSeghalCOHPiPXSec

\brief    Computes the double differential cross section for CC & NC coherent
          pion production according to the \b Rein-Seghal model.
          v(vbar)A->v(vbar)Api0, vA->l-Api+, vbarA->l+Api-

          The t-dependence of the triple differential cross (d^3xsec/dxdydt)
          is analytically integrated out.

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      D.Rein and L.M.Seghal, Coherent pi0 production in neutrino
          reactions, Nucl.Phys.B223:29-144 (1983)

          D.Rein and L.M.Sehgal, PCAC and the Deficit of Forward Muons in pi+ 
          Production by Neutrinos, hep-ph/0606185

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  March 11, 2005

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _REIN_SEGHAL_COHPI_PXSEC_H_
#define _REIN_SEGHAL_COHPI_PXSEC_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class ReinSeghalCOHPiPXSec : public XSecAlgorithmI {

public:
  ReinSeghalCOHPiPXSec();
  ReinSeghalCOHPiPXSec(string config);
  virtual ~ReinSeghalCOHPiPXSec();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void LoadConfig(void);

  //-- private data members loaded from config Registry or set to defaults
  double fMa;      ///< axial mass
  double fReIm;    ///< Re/Im {forward pion scattering amplitude}
  double fRo;      ///< nuclear size scale parameter
  bool   fModPCAC; ///< use modified PCAC (including f/s lepton mass)

  const XSecIntegratorI * fXSecIntegrator;
};

}       // genie namespace

#endif  // _REIN_SEGHAL_COHPI_PXSEC_H_
