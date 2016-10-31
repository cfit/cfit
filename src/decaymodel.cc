
#include <cfit/decaymodel.hh>
#include <cfit/function.hh>

#include <Minuit/FunctionMinimum.h>
#include <Minuit/MnUserParameters.h>
#include <Minuit/MinuitParameter.h>

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

  push( amp );
}


// Need to overwrite setters defined in PdfModel, since function parameters may need to be set.
// Set the parameters to those given as argument.
// They must be sorted alphabetically, since it's how MnUserParameters are passed
//    in the minimize function of the minimizers. It must be so, because pushing
//    two parameters with the same name would create confusion otherwise.
void DecayModel::setPars( const std::vector< double >& pars ) throw( PdfException )
{
  setParMap( pars );

  _amp.setPars( _parMap );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setPars( _parMap );

  setParExpr();
}


// Need to overwrite setters defined in PdfModel, since function parameters may need to be set.
void DecayModel::setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException )
{
  setParMap( pars );

  _amp.setPars( pars );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setPars( _parMap );

  setParExpr();
}


// Need to overwrite setters defined in PdfModel, since function parameters may need to be set.
void DecayModel::setPars( const FunctionMinimum& min ) throw( PdfException )
{
  setParMap( min );

  _amp.setPars( _parMap );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setPars( _parMap );

  setParExpr();
}

