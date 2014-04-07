#ifndef __DECAYMODEL_HH__
#define __DECAYMODEL_HH__

#include <cfit/pdfmodel.hh>
#include <cfit/variable.hh>
#include <cfit/amplitude.hh>
#include <cfit/phasespace.hh>


class FunctionMinimum;

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

  virtual DecayModel* copy() const = 0;

  virtual void setPars( const std::vector< double >& pars              ) throw( PdfException );
  virtual void setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException );
  virtual void setPars( const FunctionMinimum&                    pars ) throw( PdfException );

  const std::string mSq12name() const { return getVar( 0 ).name(); }
  const std::string mSq13name() const { return getVar( 1 ).name(); }
  const std::string mSq23name() const { return getVar( 2 ).name(); }

  const double mSq12() const { return getVar( 0 ).value(); }
  const double mSq13() const { return getVar( 1 ).value(); }
  const double mSq23() const { return getVar( 2 ).value(); }
};

#endif
