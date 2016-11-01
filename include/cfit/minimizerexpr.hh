#ifndef __MINIMIZEREXPR_HH__
#define __MINIMIZEREXPR_HH__

#include <string>
#include <vector>

#include <Minuit/FCNBase.h>
#include <Minuit/FunctionMinimum.h>

#include <cfit/minimizer.hh>

class MinimizerExpr : public FCNBase
{
private:
  double _up;
  bool   _verbose;

  std::vector< const Minimizer* > _minimizers;
  std::map< std::string, Parameter > _parMap;

  void clear()
  {
    for ( std::vector< const Minimizer* >::iterator mmzr = _minimizers.begin(); mmzr != _minimizers.end(); ++mmzr )
      delete *mmzr;

    _parMap.clear();
  }

public:
  MinimizerExpr()
    : _up( -1.0 ), _verbose( false )
    {}

  ~MinimizerExpr()
  {
    clear();
  }

  // Getters.
  double up() const throw( MinimizerException );

  // Setters.
  void setUp  ( const double& up  ) { _up      = up;  }
  void verbose( const bool&   val ) { _verbose = val; }

  double operator()( const std::vector< double >& par ) const throw( PdfException );

  FunctionMinimum minimize() const;

  // Assignment operators. No need to define operator=( const MinimizerExpr& ) because
  //    the default one is just fine.
  MinimizerExpr& operator= ( const Minimizer&     right );
  MinimizerExpr& operator+=( const Minimizer&     right );
  MinimizerExpr& operator+=( const MinimizerExpr& right );

  // Addition operator.
  friend MinimizerExpr operator+( const Minimizer& left, const Minimizer& right );
};

// Need to declare the operator here, or otherwise there would be no declaration.
MinimizerExpr operator+( const Minimizer& left, const Minimizer& right );

#endif
