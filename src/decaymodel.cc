
#include <cfit/decaymodel.hh>

DecayModel::DecayModel( const Variable&   mSq12,
                        const Variable&   mSq13,
                        const Variable&   mSq23,
                        const Amplitude&  amp  ,
                        const PhaseSpace& ps    )
  : _amp( amp ), _ps( ps )
{
  const std::map< std::string, Parameter >& pars = _amp.getParameters();
  typedef std::map< const std::string, Parameter >::const_iterator pIter;
  for ( pIter par = pars.begin(); par != pars.end(); ++par )
    push( par->second );

  push( mSq12 );
  push( mSq13 );
  push( mSq23 );
}

