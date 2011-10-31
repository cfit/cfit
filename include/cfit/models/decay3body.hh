#ifndef __DECAYMIXING3BODY_HH__
#define __DECAYMIXING3BODY_HH__

#include <vector>

#include <cfit/decaymodel.hh>
#include <cfit/variable.hh>
#include <cfit/amplitude.hh>
#include <cfit/phasespace.hh>

class Decay3Body : public DecayModel
{
private:
  double _norm;

public:
  Decay3Body( const Variable&   mSq12,
	      const Variable&   mSq13,
	      const Variable&   mSq23,
	      const Amplitude&  amp  ,
	      const PhaseSpace& ps     );

  void setPars( const std::map< std::string, Parameter >& pars );

  void cache();
  double evaluate(                                   ) const throw( PdfException );
  double evaluate( const std::vector< double >& vars ) const throw( PdfException );
};

#endif
