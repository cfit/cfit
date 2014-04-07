#ifndef __DECAY3BODY_HH__
#define __DECAY3BODY_HH__

#include <vector>

#include <cfit/decaymodel.hh>
#include <cfit/variable.hh>
#include <cfit/amplitude.hh>
#include <cfit/phasespace.hh>
#include <cfit/function.hh>

#include <Minuit/FunctionMinimum.h>


class Decay3Body : public DecayModel
{
private:
  double _norm;

  std::vector< Function > _funcs;

  const double evaluateFuncs() const;
  const double evaluateFuncs( const double& mSq12, const double& mSq13, const double& mSq23 ) const;

  // Auxiliary function to compute the center of a bin.
  static const double binCenter( const unsigned& bin, const unsigned& nbins, const double& min, const double& max )
  {
    return ( max - min ) / double( nbins ) * ( bin + 0.5 ) + min;
  }

public:
  Decay3Body( const Variable&   mSq12,
	      const Variable&   mSq13,
	      const Variable&   mSq23,
	      const Amplitude&  amp  ,
	      const PhaseSpace& ps     );

  Decay3Body* copy() const;

  // Need to define own setters, since function variables and parameters may need to be set, too.
  void setVars( const std::vector< double >&              vars ) throw( PdfException );
  void setVars( const std::map< std::string, Variable >&  vars ) throw( PdfException );
  void setVars( const std::map< std::string, double >&    vars ) throw( PdfException );
  void setPars( const std::vector< double >&              pars ) throw( PdfException );
  void setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException );
  void setPars( const FunctionMinimum&                    min  ) throw( PdfException );

  void cache();
  const double evaluate(                                                               ) const throw( PdfException );
  const double evaluate( const double& mSq12, const double& mSq13, const double& mSq23 ) const throw( PdfException );
  const double evaluate( const double& mSq12, const double& mSq13                      ) const throw( PdfException );
  const double evaluate( const std::vector< double >& vars                             ) const throw( PdfException );

  const double project ( const std::string& varName, const double& value ) const throw( PdfException );

  const std::map< std::string, double > generate() const throw( PdfException );

  friend const Decay3Body  operator* (       Decay3Body left, const Function&  right );
  friend const Decay3Body  operator* ( const Function&  left,       Decay3Body right );
  const        Decay3Body& operator*=( const Function& right ) throw( PdfException );
};

#endif
