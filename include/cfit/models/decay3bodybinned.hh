#ifndef __DECAY3BODYBINNED_HH__
#define __DECAY3BODYBINNED_HH__

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

class Decay3BodyBinned : public PdfModel
{
private:
  BinnedAmplitude _amp;

  bool          _hasKappa;
  ParameterExpr _kappa;

  CoefExpr  _z;

  // Constants to speed up norm calculation.
  double _nDir;
  double _nXed;
  double _norm;

  bool   _fixedAmp;

  // // Maximum value of the pdf.
  // double _maxPdf;

  // Indices of the cached direct and conjugated amplitudes.
  bool     _cacheAmps;
  unsigned _ampDirCache;
  unsigned _ampCnjCache;

  // const double evaluateFuncs() const;
  // const double evaluateFuncs( const int& bin ) const;

  const double evaluateUnnorm( const int& bin ) const throw( PdfException );

  // const std::map< unsigned, std::vector< std::complex< double > > > cacheComplex( const Dataset& data );

  void setParExpr();

  // One or more functions to define the efficiency.
  std::vector< Function > _funcs;

public:
  Decay3BodyBinned( const Variable&        bin           ,
                    const BinnedAmplitude& amp           ,
                    const CoefExpr&        z             ,
                    bool                   docache = true );

  Decay3BodyBinned( const Variable&        bin           ,
                    const BinnedAmplitude& amp           ,
                    const CoefExpr&        z             ,
                    const Parameter&       kappa         ,
                    bool                   docache = true );

  Decay3BodyBinned( const Variable&        bin           ,
                    const BinnedAmplitude& amp           ,
                    const CoefExpr&        z             ,
                    const ParameterExpr&   kappa         ,
                    bool                   docache = true );

  Decay3BodyBinned* copy() const;

  // Getters.
  const std::complex< double > z()     const { return             _z    .evaluate(); }
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

  // const double evaluate(                                   ) const throw( PdfException );
  const double evaluate( const int&                   bin  ) const throw( PdfException );
  const double evaluate( const std::vector< double >& vars ) const throw( PdfException );

  const double evaluate( const std::vector< double >&                 vars  ,
                         const std::vector< double >&                 cacheR,
                         const std::vector< std::complex< double > >& cacheC ) const throw( PdfException );

  // const double project ( const std::string& varName, const double& value ) const throw( PdfException );

  // const std::map< std::string, double > generate() const throw( PdfException );

  // friend const Decay3BodyBinned  operator* (       Decay3BodyBinned left, const Function&        right );
  // friend const Decay3BodyBinned  operator* ( const Function&        left,       Decay3BodyBinned right );
  // const        Decay3BodyBinned& operator*=( const Function&        right                              ) throw( PdfException );
};

#endif
