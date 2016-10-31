#ifndef __DECAYMODEL_HH__
#define __DECAYMODEL_HH__

#include <cfit/pdfmodel.hh>
#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/amplitude.hh>
#include <cfit/phasespace.hh>

class FunctionMinimum;

class DecayModel : public PdfModel
{
protected:
  Amplitude  _amp;
  PhaseSpace _ps;

  // One or more functions to define the efficiency.
  std::vector< Function > _funcs;

public:
  DecayModel( const Variable&   mSq12,
              const Variable&   mSq13,
              const Variable&   mSq23,
              const Amplitude&  amp  ,
              const PhaseSpace& ps    );

  virtual DecayModel* copy() const = 0;

  void setPars( const std::vector< double >&              pars ) throw( PdfException );
  void setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException );
  void setPars( const FunctionMinimum&                    pars ) throw( PdfException );

  const std::string mSq12name() const { return getVar( 0 ).name(); }
  const std::string mSq13name() const { return getVar( 1 ).name(); }
  const std::string mSq23name() const { return getVar( 2 ).name(); }
};

#endif
