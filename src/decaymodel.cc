
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


void DecayModel::setPars( const std::vector< double >& pars ) throw( PdfException )
{
  if ( _parMap.size() != pars.size() )
    throw PdfException( "Number of arguments passed does not match number of required arguments." );

  typedef std::map< std::string, Parameter >::iterator pIter;
  int index = 0;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars[ index++ ] );

  _amp.setPars( _parMap );
}

void DecayModel::setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException )
{
  typedef std::map< std::string, Parameter >::iterator pIter;
  for ( pIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars.find( par->first )->second.value() );

  _amp.setPars( pars );
}
