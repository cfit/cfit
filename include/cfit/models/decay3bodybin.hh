#ifndef __DECAY3BODYBIN_HH__
#define __DECAY3BODYBIN_HH__

#include <vector>
#include <map>

#include <cfit/decaymodel.hh>
#include <cfit/variable.hh>
#include <cfit/amplitude.hh>
#include <cfit/phasespace.hh>
#include <cfit/function.hh>
#include <cfit/binnedamplitude.hh>

#include <Minuit/FunctionMinimum.h>

class Dataset;

class Decay3BodyBin : public DecayModel< BinnedAmplitude >
{
private:
  Binning       _binning;

  bool          _hasKappa;
  ParameterExpr _kappa;

  CoefExpr      _phi;

  // Constants to speed up norm calculation.
  double _nDir;
  double _nXed;
  double _norm;

  bool   _fixedAmp;

  // // Maximum value of the pdf.
  // double _maxPdf;

  // Index of the cached bins.
  unsigned _binIndex;

  std::vector< unsigned > _binCache;

  // const double evaluateUnnorm( const int& bin ) const throw( PdfException );
  const double evaluateUnnorm( const double& mSq12, const double& mSq13 ) const throw( PdfException );

  const std::map< unsigned, std::vector< double > > cacheReal( const Dataset& data );

  void setParExpr();

  // One or more functions to define the efficiency.
  std::vector< Function > _funcs;

public:
  Decay3BodyBin( const Variable&        mSq12         ,
                 const Variable&        mSq13         ,
                 const Variable&        mSq23         ,
                 const BinnedAmplitude& amp           ,
                 const Binning&         binning       ,
                 const CoefExpr&        z             ,
                 const PhaseSpace&      ps            ,
                 bool                   docache = true );

  Decay3BodyBin( const Variable&        mSq12         ,
                 const Variable&        mSq13         ,
                 const Variable&        mSq23         ,
                 const BinnedAmplitude& amp           ,
                 const Binning&         binning       ,
                 const CoefExpr&        z             ,
                 const Parameter&       kappa         ,
                 const PhaseSpace&      ps            ,
                 bool                   docache = true );

  Decay3BodyBin( const Variable&        mSq12         ,
                 const Variable&        mSq13         ,
                 const Variable&        mSq23         ,
                 const BinnedAmplitude& amp           ,
                 const Binning&         binning       ,
                 const CoefExpr&        z             ,
                 const ParameterExpr&   kappa         ,
                 const PhaseSpace&      ps            ,
                 bool                   docache = true );

  Decay3BodyBin* copy() const;

  // Getters.
  const std::complex< double > phi()   const { return             _phi  .evaluate();       }
  const std::complex< double > z()     const { return             std::tanh( phi() );      }
  const double                 kappa() const { return _hasKappa ? _kappa.evaluate() : 1.0; }

  // Norm components getters.
  const double& nDir() const { return _nDir; }
  const double& nXed() const { return _nXed; }

  // Norm components setters.
  void setNormComponents( const double& nDir, const double& nXed )
  {
    _fixedAmp = true; //_amp.isFixed();

    if ( _fixedAmp )
    {
      _nDir = nDir;
      _nXed = nXed;
    }
  }

  // Caching functions.
  void cacheNormComponents();
  void cache();

  const double evaluate( const double& mSq12, const double& mSq13, const double& ) const throw( PdfException );
  const double evaluate( const double& mSq12, const double& mSq13                ) const throw( PdfException );
  const double evaluate( const std::vector< double >& vars ) const throw( PdfException );

  const double evaluate( const std::vector< double >&                 vars  ,
                         const std::vector< double >&                 cacheR,
                         const std::vector< std::complex< double > >& cacheC ) const throw( PdfException );

  friend const Decay3BodyBin  operator* (       Decay3BodyBin left, const Function&     right );
  friend const Decay3BodyBin  operator* ( const Function&     left,       Decay3BodyBin right );
  const        Decay3BodyBin& operator*=( const Function&     right                           ) throw( PdfException );
};

#endif
