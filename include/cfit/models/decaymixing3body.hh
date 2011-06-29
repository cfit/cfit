#ifndef __DECAYMIXING3BODY_HH__
#define __DECAYMIXING3BODY_HH__

#include <vector>

#include <cfit/decaymodel.hh>
#include <cfit/variable.hh>
#include <cfit/amplitude.hh>
#include <cfit/phasespace.hh>

class DecayMixing3Body : public DecayModel
{
private:
  const Variable  _t;
  const Parameter _tau;
  const Parameter _x;
  const Parameter _y;

  std::complex< double > e1( const double& gt ) const;
  std::complex< double > e2( const double& gt ) const;
public:
  DecayMixing3Body( const Variable&   mSq12,
		    const Variable&   mSq13,
		    const Variable&   mSq23,
		    const Variable&   t    ,
		    const Amplitude&  amp  ,
		    const Parameter&  tau  ,
		    const Parameter&  x    ,
		    const Parameter&  y    ,
		    const PhaseSpace& ps     );

  double evaluate(                                   ) const throw( PdfException );
  double evaluate( const std::vector< double >& vars ) const throw( PdfException );
};

#endif
