#ifndef __DECAYMODEL_HH__
#define __DECAYMODEL_HH__

#include <cfit/pdfmodel.hh>
#include <cfit/variable.hh>
#include <cfit/amplitude.hh>
#include <cfit/phasespace.hh>

class DecayModel : public PdfModel
{
protected:
  Variable _mSq12;
  Variable _mSq13;
  Variable _mSq23;

  Amplitude  _amp;
  PhaseSpace _ps;
};

#endif
