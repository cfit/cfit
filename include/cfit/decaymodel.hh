#ifndef __DECAYMODEL_HH__
#define __DECAYMODEL_HH__

#include <cfit/pdfmodel.hh>
#include <cfit/variable.hh>
#include <cfit/parameter.hh>
#include <cfit/amplitude.hh>
#include <cfit/binnedamplitude.hh>
#include <cfit/phasespace.hh>


#include <cfit/function.hh>

class FunctionMinimum;

template< class AmplitudeClass >
class DecayModel : public PdfModel
{
protected:
  AmplitudeClass _amp;
  PhaseSpace     _ps;

  // One or more functions to define the efficiency.
  std::vector< Function > _funcs;

public:
  DecayModel< AmplitudeClass >( const Variable&       mSq12,
                                const Variable&       mSq13,
                                const Variable&       mSq23,
                                const AmplitudeClass& amp  ,
                                const PhaseSpace&     ps    )
    : _amp( amp ), _ps( ps )
  {
    push( mSq12 );
    push( mSq13 );
    push( mSq23 );

    push( amp );
  }

  virtual DecayModel< AmplitudeClass >* copy() const = 0;

  void setPars( const std::vector< double >&              pars ) throw( PdfException );
  void setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException );
  void setPars( const FunctionMinimum&                    pars ) throw( PdfException );

  const std::string mSq12name() const { return getVar( 0 ).name(); }
  const std::string mSq13name() const { return getVar( 1 ).name(); }
  const std::string mSq23name() const { return getVar( 2 ).name(); }

  const double evaluateFuncs( const double& mSq12, const double& mSq13, const double& mSq23 ) const;
  const double evaluateFuncs( const double& mSq12, const double& mSq13                      ) const;
};


// Need to overwrite setters defined in PdfModel, since function parameters may need to be set.
// Set the parameters to those given as argument.
// They must be sorted alphabetically, since it's how MnUserParameters are passed
//    in the minimize function of the minimizers. It must be so, because pushing
//    two parameters with the same name would create confusion otherwise.
template < class AmplitudeClass >
inline
void DecayModel< AmplitudeClass >::setPars( const std::vector< double >& pars ) throw( PdfException )
{
  setParMap( pars );

  _amp.setPars( _parMap );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setPars( _parMap );

  setParExpr();
}


// Need to overwrite setters defined in PdfModel, since function parameters may need to be set.
template < class AmplitudeClass >
inline
void DecayModel< AmplitudeClass >::setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException )
{
  setParMap( pars );

  _amp.setPars( pars );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setPars( _parMap );

  setParExpr();
}


// Need to overwrite setters defined in PdfModel, since function parameters may need to be set.
template < class AmplitudeClass >
inline
void DecayModel< AmplitudeClass >::setPars( const FunctionMinimum& min ) throw( PdfException )
{
  setParMap( min );

  _amp.setPars( _parMap );

  typedef std::vector< Function >::iterator fIter;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
    func->setPars( _parMap );

  setParExpr();
}


template < class AmplitudeClass >
inline
const double DecayModel< AmplitudeClass >::evaluateFuncs( const double& mSq12, const double& mSq13, const double& mSq23 ) const
{
  double value = 1.0;

  const std::string& name12 = mSq12name(); // mSq12
  const std::string& name13 = mSq13name(); // mSq13
  const std::string& name23 = mSq23name(); // mSq23

  typedef std::vector< Function >::const_iterator fIter;

  std::map< std::string, double > varMap;
  for ( fIter func = _funcs.begin(); func != _funcs.end(); ++func )
  {
    if ( func->dependsOn( name12 ) ) varMap[ name12 ] = mSq12;
    if ( func->dependsOn( name13 ) ) varMap[ name13 ] = mSq13;
    if ( func->dependsOn( name23 ) ) varMap[ name23 ] = mSq23;

    value *= func->evaluate( varMap );
  }

  // Always return a non-negative value. Default to zero.
  return std::max( value, 0.0 );
}


template < class AmplitudeClass >
inline
const double DecayModel< AmplitudeClass >::evaluateFuncs( const double& mSq12, const double& mSq13 ) const
{
  const double& mSq23 = _ps.mSqSum() - mSq12 - mSq13;

  return evaluateFuncs( mSq12, mSq13, mSq23 );
};

#endif
