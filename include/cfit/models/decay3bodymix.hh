#ifndef __DECAY3BODYMIX_HH__
#define __DECAY3BODYMIX_HH__

#include <vector>
#include <map>

#include <cfit/pdfmodel.hh>
#include <cfit/variable.hh>
#include <cfit/amplitude.hh>
#include <cfit/phasespace.hh>
#include <cfit/function.hh>

#include <Minuit/FunctionMinimum.h>


class Dataset;

class Decay3BodyMix : public PdfModel
{
private:
  // Names of the squared invariant mass and life time variables.
  std::string _mSq12;
  std::string _mSq13;
  std::string _mSq23;
  std::string _t;

  ParameterExpr _width;

  Amplitude _amp;

  CoefExpr  _z;

  PhaseSpace _ps;

  // Constants to speed up norm calculation.
  double                 _nDir;
  double                 _nCnj;
  std::complex< double > _nXed;
  double                 _norm;

  bool                   _fixedAmp;

  // Maximum value of the pdf.
  double _maxPdf;

  // Indices of the cached direct and conjugated amplitudes.
  bool     _cacheAmps;
  unsigned _ampDirCache;
  unsigned _ampCnjCache;


  std::vector< Function > _funcs;

  // const double evaluateFuncs() const;
  const double evaluateFuncs( const double& mSq12, const double& mSq13, const double& mSq23 ) const;
  const double evaluateFuncs( const double& mSq12, const double& mSq13                      ) const;

  // Auxiliary function to compute the center of a bin.
  static const double binCenter( const unsigned& bin, const unsigned& nbins, const double& min, const double& max )
  {
    return ( max - min ) / double( nbins ) * ( bin + 0.5 ) + min;
  }

  const double evaluateUnnorm( const double& mSq12, const double& mSq13, const double& mSq23, const double& t ) const;

  void cacheNormComponents();

  const std::map< unsigned, std::vector< std::complex< double > > > cacheComplex( const Dataset& data );

  const double                 psip( const double& t ) const;
  const double                 psim( const double& t ) const;
  const std::complex< double > psii( const double& t ) const;

public:
  Decay3BodyMix( const Variable&      mSq12         ,
                 const Variable&      mSq13         ,
                 const Variable&      mSq23         ,
                 const Variable&      t             ,
                 const ParameterExpr& width         ,
                 const Amplitude&     amp           ,
                 const CoefExpr&      z             ,
                 const PhaseSpace&    ps            ,
                 bool                 docache = true  );

  Decay3BodyMix* copy() const;

  // Getters.
  const double                 gamma() const { return _width.evaluate();          }
  const double                 tau()   const { return 1.0 / gamma();              }
  const double                 x()     const { return std::real( _z.evaluate() ); }
  const double                 y()     const { return std::imag( _z.evaluate() ); }
  const std::complex< double > z()     const { return _z.evaluate();              }

  // Getters for the norm components.
  const double&                 nDir() const { return _nDir; }
  const double&                 nCnj() const { return _nCnj; }
  const std::complex< double >& nXed() const { return _nXed; }

  // Need to define own setters, since function variables and parameters may need to be set, too.
  void setPars( const std::vector< double >&              pars ) throw( PdfException );
  void setPars( const std::map< std::string, Parameter >& pars ) throw( PdfException );
  void setPars( const FunctionMinimum&                    min  ) throw( PdfException );

  // Norm components setters.
  void setNormComponents( const double& nDir, const double& nCnj, const std::complex< double >& nXed )
  {
    _fixedAmp = _amp.isFixed();

    if ( _fixedAmp )
    {
      _nDir = nDir;
      _nCnj = nCnj;
      _nXed = nXed;
    }

    cache();
  }

  void setNormComponents( const double& nDir, const std::complex< double >& nXed )
  {
    _fixedAmp = _amp.isFixed();

    if ( _fixedAmp )
    {
      _nDir = nDir;
      _nCnj = nDir;
      _nXed = nXed;
    }

    cache();
  }


  void cache();
  const double evaluate( const double& mSq12, const double& mSq13, const double& mSq23, const double& t ) const;
  const double evaluate( const double& mSq12, const double& mSq13                     , const double& t ) const;

  const double evaluate( const std::vector< double >&                 vars    ) const throw( PdfException );
  const double evaluate( const std::vector< double >&                 vars  ,
                         const std::vector< double >&                 cacheR,
                         const std::vector< std::complex< double > >& cacheC  ) const throw( PdfException );

  const double project ( const std::string& varName, const double& value ) const throw( PdfException );

  void setMaxPdf( const double& max ) { _maxPdf = max; }
  const std::map< std::string, double > generate() const throw( PdfException );

  friend const Decay3BodyMix  operator* (       Decay3BodyMix left, const Function&     right );
  friend const Decay3BodyMix  operator* ( const Function&     left,       Decay3BodyMix right );
  const        Decay3BodyMix& operator*=( const Function& right );
};

#endif
