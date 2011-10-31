
#include <cfit/decaymodel.hh>

DecayModel::DecayModel( const Variable&   mSq12,
                        const Variable&   mSq13,
                        const Variable&   mSq23,
                        const Amplitude&  amp  ,
                        const PhaseSpace& ps    )
  : _amp( amp ), _ps( ps )
{
  push( mSq12 );
  push( mSq13 );
  push( mSq23 );
}

