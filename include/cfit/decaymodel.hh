#ifndef __DECAYMODEL_HH__
#define __DECAYMODEL_HH__

#include <cfit/pdfmodel.hh>
#include <cfit/variable.hh>
#include <cfit/amplitude.hh>
#include <cfit/phasespace.hh>

class DecayModel : public PdfModel
{
protected:
  Amplitude  _amp;
  PhaseSpace _ps;
public:
  DecayModel( const Variable&   mSq12,
              const Variable&   mSq13,
              const Variable&   mSq23,
              const Amplitude&  amp  ,
              const PhaseSpace& ps    );

  double mSq12() const { return getVar( 0 ).value(); }
  double mSq13() const { return getVar( 1 ).value(); }
  double mSq23() const { return getVar( 2 ).value(); }
};

#endif
