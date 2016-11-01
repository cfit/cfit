#ifndef __MINIMIZER_HH__
#define __MINIMIZER_HH__

#include <vector>

#include <Minuit/FCNBase.h>
#include <Minuit/FunctionMinimum.h>

#include <cfit/variable.hh>
#include <cfit/dataset.hh>
#include <cfit/pdfbase.hh>


class Minimizer : public FCNBase
{
private:
  void cache();

protected:
  PdfBase*       _pdf;
  const Dataset& _data;

  // Minimizer variation to produce uncertaities at a given number of sigmas.
  //    Notice that, if the user wants n-sigma uncertainties, up = n^2.
  double _up;

  bool   _verbose;

  // Maps of cached expressions.
  std::map< unsigned, std::vector< double >                 > _cacheR;
  std::map< unsigned, std::vector< std::complex< double > > > _cacheC;

public:
  Minimizer( const PdfBase& pdf, const Dataset& data )
    : _pdf    ( pdf.copy() ),
      _data   ( data       ),
      _up     ( -1.0       ),
      _verbose( false      )
  {
    cache();
  }

  // Copy constructor.
  Minimizer( const Minimizer& minimizer )
    : _pdf    ( minimizer._pdf->copy() ),
      _data   ( minimizer._data        ),
      _up     ( minimizer._up          ),
      _verbose( minimizer._verbose     ),
      _cacheR ( minimizer._cacheR      ),
      _cacheC ( minimizer._cacheC      )
    {}

  virtual Minimizer* copy() const = 0;

  virtual ~Minimizer()
  {
    delete _pdf;
  }

  const PdfBase& pdf()  const { return *_pdf; }
  const Dataset& data() const { return _data; }

  // Should be const double, but minuit declares these functions as double.
  double up() const throw( MinimizerException );
  double operator()( const std::vector<double>& par ) const throw( PdfException ) = 0;

  // Setters.
  void setUp  ( const double& up         ) { _up      = up;  }
  void verbose( const bool&   val = true ) { _verbose = val; }

  FunctionMinimum minimize() const;
};

#endif
