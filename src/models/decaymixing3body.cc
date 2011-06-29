
#include <cfit/models/decaymixing3body.hh>

DecayMixing3Body::DecayMixing3Body( const Variable&   mSq12,
				    const Variable&   mSq13,
				    const Variable&   mSq23,
				    const Variable&   t    ,
				    const Amplitude&  amp  ,
				    const Parameter&  tau  ,
				    const Parameter&  x    ,
				    const Parameter&  y    ,
				    const PhaseSpace& ps     )
  : _t( t ), _tau( tau ), _x( x ), _y( y )
{
  _mSq12 = mSq12;
  _mSq13 = mSq13;
  _mSq23 = mSq23;
  _amp   = amp;
  _ps    = ps;
}

double DecayMixing3Body::evaluate() const throw( PdfException )
{
  // Phase space amplitudes of the decay of the particle and that of the antiparticle.
  std::complex< double > ampDirect = _amp.evaluate( _ps, _mSq12.value(), _mSq13.value(), _mSq23.value() );
  std::complex< double > ampConjug = _amp.evaluate( _ps, _mSq13.value(), _mSq12.value(), _mSq23.value() );

  // Assume there's no direct CP violation.
  std::complex< double > barAOverA = ampConjug / ampDirect;

  // TEMPORARY.
  double _qp = 1;
  double _norm = 1;

  // CP violation in the interference. _qp implements CP violation in the mixing.
  std::complex< double > chi = _qp * barAOverA;

  // Gamma t. Notice that Gamma is not 1 / tau, since it also depends on x and y,
  //    as well as on integrals of products of amplitudes over the Dalitz plot.
  double gt = _t.value() / _tau.value();

  // Compute time dependent amplitude.
  std::complex< double > amp = .5 * ampDirect * ( ( 1. + chi ) * e1( gt ) + ( 1. - chi ) * e2( gt ) );

  // std::norm returns the squared modulus of the complex number, not its norm.
  return std::norm( amp ) / _norm;
}


double DecayMixing3Body::evaluate( const std::vector< double >& vars ) const throw( PdfException )
{
  throw PdfException( "Do not use evaluate( vars ) in decaymixing3body" );
}


// < f | H | D^0 (t) > = 1/2 * [ ( 1 + \chi_f ) * A_f * e_1(gt) + ( 1 - \chi_f ) * A_f * e_2(gt) ]
std::complex< double > DecayMixing3Body::e1( const double& gt ) const
{
  double real = exp ( - ( 1. + _y.value() ) * gt / 2. ) * cos(   _x.value() * gt / 2. );
  double imag = exp ( - ( 1. + _y.value() ) * gt / 2. ) * sin( - _x.value() * gt / 2. );

  return std::complex< double >( real, imag );
}

std::complex< double > DecayMixing3Body::e2( const double& gt ) const
{
  double real = exp ( - ( 1. - _y.value() ) * gt / 2. ) * cos( _x.value() * gt / 2. );
  double imag = exp ( - ( 1. - _y.value() ) * gt / 2. ) * sin( _x.value() * gt / 2. );

  return std::complex< double >( real, imag );
}
