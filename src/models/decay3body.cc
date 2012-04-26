
#include <cfit/models/decay3body.hh>

Decay3Body::Decay3Body( const Variable&   mSq12,
			const Variable&   mSq13,
			const Variable&   mSq23,
			const Amplitude&  amp  ,
			const PhaseSpace& ps     )
  : DecayModel( mSq12, mSq13, mSq23, amp, ps ), _norm( 1. )
{
  // Do calculations common to all values of variables
  //    (usually compute norm).
  cache();
}


void Decay3Body::cache()
{
  // Compute the value of _norm.
  _norm = 0.0;

  // Define the properties of the integration method.
  const int    nBins = 400;
  const double min   = pow( _ps.m1()      + _ps.m2(), 2 );
  const double max   = pow( _ps.mMother() - _ps.m3(), 2 );
  const double step  = ( max - min ) / double( nBins );

  const double mSqSum = _ps.mSqMother() + _ps.mSq1() + _ps.mSq2() + _ps.mSq3();

  // Define the variables at each bin.
  double mSq12;
  double mSq13;
  double mSq23;

  // Compute the integral on the grid.
  for ( int binX = 0; binX < nBins; ++binX )
    for ( int binY = 0; binY < nBins; ++binY )
    {
      mSq12 = min + step * ( binX + .5 );
      mSq13 = min + step * ( binY + .5 );
      mSq23 = mSqSum - mSq12 - mSq13;

      // Proceed only if the point lies inside the kinematically allowed Dalitz region.
      // std::norm returns the squared modulus of the complex number, not its norm.
      if ( _ps.contains( mSq12, mSq13, mSq23 ) )
        _norm += std::norm( _amp.evaluate( _ps, mSq12, mSq13, mSq23 ) );
    }

  _norm *= pow( step, 2 );

  return;
}


double Decay3Body::evaluate() const throw( PdfException )
{
  // Phase space amplitude of the decay of the particle.
  std::complex< double > amp = _amp.evaluate( _ps, mSq12(), mSq13(), mSq23() );

  // std::norm returns the squared modulus of the complex number, not its norm.
  return std::norm( amp ) / _norm;
}


double Decay3Body::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  throw PdfException( "Do not use evaluate( vars ) in decay3body" );
}

