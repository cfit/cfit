
#include <stack>
#include <algorithm>

#include <cfit/phasespace.hh>
#include <cfit/amplitude.hh>
#include <cfit/coef.hh>
#include <cfit/resonance.hh>
#include <cfit/operation.hh>

void Amplitude::append( const double& ctnt )
{
  _ctnts.push_back( std::complex< double >( ctnt ) );
  _expression += "c"; // c = constant.
}

void Amplitude::append( const std::complex< double >& ctnt )
{
  _ctnts.push_back( ctnt );
  _expression += "c"; // c = constant.
}

void Amplitude::append( const Coef& coef )
{
  _coefs.push_back( coef );

  _parMap[ coef.real().name() ] = coef.real();
  _parMap[ coef.imag().name() ] = coef.imag();

  _expression += "k"; // k = coefficient.
}

void Amplitude::append( const CoefExpr& expr )
{
  _ctnts.insert( _ctnts.end(), expr._ctnts.begin(), expr._ctnts.end() );
  _parms.insert( _parms.end(), expr._parms.begin(), expr._parms.end() );
  _coefs.insert( _coefs.end(), expr._coefs.begin(), expr._coefs.end() );
  _opers.insert( _opers.end(), expr._opers.begin(), expr._opers.end() );

  typedef std::vector< Parameter >::const_iterator pIter;
  for ( pIter par = expr._parms.begin(); par != expr._parms.end(); ++par )
    _parMap[ par->name() ] = *par;

  typedef std::vector< Coef >::const_iterator cIter;
  for ( cIter coef = expr._coefs.begin(); coef != expr._coefs.end(); ++coef )
  {
    _parMap[ coef->real().name() ] = coef->real();
    _parMap[ coef->imag().name() ] = coef->imag();
  }

  _expression += expr._expression;
}

void Amplitude::append( const Resonance& reso )
{
  _resos.push_back( reso.copy() );

  _parMap.insert( reso._parMap.begin(), reso._parMap.end() );

  _expression += "r"; // r = resonance.
}

void Amplitude::append( const Fvector& fvec )
{
  _fvecs.push_back( fvec );

  _parMap.insert( fvec._parMap.begin(), fvec._parMap.end() );

  _expression += "F"; // F = element of F vector.
}

void Amplitude::append( const Amplitude& ampl )
{
  _ctnts.insert( _ctnts.end(), ampl._ctnts.begin(), ampl._ctnts.end() );
  _parms.insert( _parms.end(), ampl._parms.begin(), ampl._parms.end() );
  _coefs.insert( _coefs.end(), ampl._coefs.begin(), ampl._coefs.end() );
  _opers.insert( _opers.end(), ampl._opers.begin(), ampl._opers.end() );

  std::transform( ampl._resos.begin(), ampl._resos.end(),
                  std::back_inserter( _resos ), std::mem_fun( &Resonance::copy ) );
  _fvecs.insert( _fvecs.end(), ampl._fvecs.begin(), ampl._fvecs.end() );

  _parMap.insert( ampl._parMap.begin(), ampl._parMap.end() );

  _expression += ampl._expression;
}

void Amplitude::append( const Operation::Op& oper )
{
  _opers.push_back( oper );
  _expression += "b"; // b = binary operation.
}


void Amplitude::setPars( const std::map< std::string, Parameter >& pars )
{
  typedef std::map< std::string, Parameter >::iterator mIter;
  for ( mIter par = _parMap.begin(); par != _parMap.end(); ++par )
    par->second.setValue( pars.find( par->first )->second.value() );

  // Set the values of the parameters.
  typedef std::vector< Parameter >::iterator pIter;
  for ( pIter par = _parms.begin(); par != _parms.end(); ++par )
    par->setValue( pars.find( par->name() )->second.value() );

  // Set the values of the coefficients.
  typedef std::vector< Coef >::iterator cIter;
  for ( cIter coef = _coefs.begin(); coef != _coefs.end(); ++coef )
    coef->setValue( pars.find( coef->real().name() )->second.value(),
                    pars.find( coef->imag().name() )->second.value() );

  // Propagate the values to the list of resonances.
  typedef std::vector< Resonance* >::iterator rIter;
  for ( rIter reso = _resos.begin(); reso != _resos.end(); ++reso )
    (*reso)->setPars( pars );

  // Propagate the values to the list of fvector components.
  typedef std::vector< Fvector >::iterator fIter;
  for ( fIter fvec = _fvecs.begin(); fvec != _fvecs.end(); ++fvec )
    fvec->setPars( pars );
}



// Evaluate the amplitude at the given point, with the current values of its parameters.
std::complex< double > Amplitude::evaluate( const PhaseSpace& ps,
                                            const double&     mSq12,
                                            const double&     mSq13,
                                            const double&     mSq23 ) const throw( PdfException )
{
  if ( ! ps.contains( mSq12, mSq13, mSq23 ) )
    return 0.;

  std::stack< std::complex< double > > values;

  std::complex< double > x;
  std::complex< double > y;

  std::vector< std::complex< double > >::const_iterator ctt = _ctnts.begin();
  std::vector< Parameter              >::const_iterator par = _parms.begin();
  std::vector< Coef                   >::const_iterator coe = _coefs.begin();
  std::vector< Resonance*             >::const_iterator res = _resos.begin();
  std::vector< Fvector                >::const_iterator fvc = _fvecs.begin();
  std::vector< Operation::Op          >::const_iterator ops = _opers.begin();

  // Parsing loop.
  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ++ch )
    if ( *ch == 'c' )
      values.push( *ctt++ );
    else if ( *ch == 'p' )
      values.push( std::complex< double >( par++->value(), 0. ) );
    else if ( *ch == 'k' )
      values.push( coe++->value() );
    else if ( *ch == 'r' )
      values.push( (*res++)->evaluate( ps, mSq12, mSq13, mSq23 ) );
    else if ( *ch == 'F' )
      values.push( fvc++->evaluate( ps, mSq12, mSq13, mSq23 ) );
    else
    {
      if ( *ch == 'b' ) // Binary operation with complex numbers.
      {
        if ( values.size() < 2 )
          throw PdfException( "Parse error: not enough values in the stack." );
        y = values.top();
        values.pop();
        x = values.top();
        values.pop();

        values.push( Operation::operate( x, y, *ops++ ) );
      }
      else if ( *ch == 'u' ) // Unary operation with complex numbers.
      {
        if ( values.empty() )
          throw PdfException( "Parse error: not enough values in the stack." );
        x = values.top();
        values.pop();

        values.push( Operation::operate( x, *ops++ ) );
      }
      else
        throw PdfException( std::string( "Parse error: unknown operation " ) + *ch + "." );
    }

  if ( values.size() != 1 )
    throw PdfException( "Parse error: too many values have been supplied." );

  return values.top();
}



// Assignment operations.
const Amplitude& Amplitude::operator=( const Resonance& right )
{
  append( right );
  return *this;
}

const Amplitude& Amplitude::operator=( const Fvector& right )
{
  append( right );
  return *this;
}

const Amplitude& Amplitude::operator=( const Amplitude& right )
{
  append( right );
  return *this;
}

const Amplitude& Amplitude::operator+=( const double& right )
{
  if ( _expression.empty() )
    append( right );
  else
  {
    append( right           );
    append( Operation::plus );
  }

  return *this;
}

const Amplitude& Amplitude::operator-=( const double& right )
{
  if ( _expression.empty() )
    throw PdfException( "Attempting to subtract a constant from an empty amplitude." );

  append( right            );
  append( Operation::minus );
  return *this;
}

const Amplitude& Amplitude::operator*=( const double& right )
{
  if ( _expression.empty() )
    throw PdfException( "Attempting to multiply an empty amplitude by a constant." );

  append( right           );
  append( Operation::mult );
  return *this;
}

const Amplitude& Amplitude::operator/=( const double& right )
{
  if ( _expression.empty() )
    throw PdfException( "Attempting to divide an empty amplitude by a constant." );

  append( right          );
  append( Operation::div );
  return *this;
}


const Amplitude& Amplitude::operator+=( const std::complex< double >& right )
{
  if ( _expression.empty() )
    append( right );
  else
  {
    append( right           );
    append( Operation::plus );
  }

  return *this;
}

const Amplitude& Amplitude::operator-=( const std::complex< double >& right )
{
  if ( _expression.empty() )
    throw PdfException( "Attempting to subtract a constant from an empty amplitude." );

  append( right            );
  append( Operation::minus );
  return *this;
}

const Amplitude& Amplitude::operator*=( const std::complex< double >& right )
{
  if ( _expression.empty() )
    throw PdfException( "Attempting to multiply an empty amplitude by a constant." );

  append( right           );
  append( Operation::mult );
  return *this;
}

const Amplitude& Amplitude::operator/=( const std::complex< double >& right )
{
  if ( _expression.empty() )
    throw PdfException( "Attempting to divide an empty amplitude by a constant." );

  append( right          );
  append( Operation::div );
  return *this;
}


const Amplitude& Amplitude::operator+=( const Coef& right )
{
  if ( _expression.empty() )
    append( right );
  else
  {
    append( right           );
    append( Operation::plus );
  }

  return *this;
}

const Amplitude& Amplitude::operator-=( const Coef& right )
{
  append( right            );
  append( Operation::minus );
  return *this;
}

const Amplitude& Amplitude::operator*=( const Coef& right )
{
  append( right           );
  append( Operation::mult );
  return *this;
}

const Amplitude& Amplitude::operator/=( const Coef& right )
{
  append( right          );
  append( Operation::div );
  return *this;
}



const Amplitude& Amplitude::operator+=( const CoefExpr& right )
{
  if ( _expression.empty() )
    append( right );
  else
  {
    append( right           );
    append( Operation::plus );
  }

  return *this;
}

const Amplitude& Amplitude::operator-=( const CoefExpr& right )
{
  append( right            );
  append( Operation::minus );
  return *this;
}

const Amplitude& Amplitude::operator*=( const CoefExpr& right )
{
  append( right           );
  append( Operation::mult );
  return *this;
}

const Amplitude& Amplitude::operator/=( const CoefExpr& right )
{
  append( right          );
  append( Operation::div );
  return *this;
}



const Amplitude& Amplitude::operator+=( const Resonance& right )
{
  if ( _expression.empty() )
    append( right );
  else
  {
    append( right           );
    append( Operation::plus );
  }

  return *this;
}

const Amplitude& Amplitude::operator-=( const Resonance& right )
{
  append( right            );
  append( Operation::minus );
  return *this;
}


const Amplitude& Amplitude::operator+=( const Fvector& right )
{
  if ( _expression.empty() )
    append( right );
  else
  {
    append( right           );
    append( Operation::plus );
  }

  return *this;
}

const Amplitude& Amplitude::operator-=( const Fvector& right )
{
  append( right            );
  append( Operation::minus );
  return *this;
}


const Amplitude& Amplitude::operator+=( const Amplitude& right )
{
  if ( _expression.empty() )
    append( right );
  else
  {
    append( right           );
    append( Operation::plus );
  }

  return *this;
}

const Amplitude& Amplitude::operator-=( const Amplitude& right )
{
  append( right            );
  append( Operation::minus );
  return *this;
}



const Amplitude operator+( const Resonance& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Resonance& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator+( const Fvector& left, const Fvector& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Fvector& left, const Fvector& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator+( const Amplitude& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Amplitude& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::minus );
}



// Operations of constants and resonances.
const Amplitude operator+( const double& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const double& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const double& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator+( const Resonance& left, const double& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Resonance& left, const double& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const Resonance& left, const double& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator/( const Resonance& left, const double& right )
{
  return Amplitude( left, right, Operation::div );
}


// Operations with constants and amplitudes.
const Amplitude operator+( const double& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const double& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const double& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator+( const Amplitude& left, const double& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Amplitude& left, const double& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const Amplitude& left, const double& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator/( const Amplitude& left, const double& right )
{
  return Amplitude( left, right, Operation::div );
}



// Operations of coefficients and resonances.
const Amplitude operator+( const Coef& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Coef& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const Coef& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator+( const Resonance& left, const Coef& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Resonance& left, const Coef& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const Resonance& left, const Coef& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator/( const Resonance& left, const Coef& right )
{
  return Amplitude( left, right, Operation::div );
}



// Operations with coefficients and amplitudes.
const Amplitude operator+( const Coef& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Coef& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const Coef& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator+( const Amplitude& left, const Coef& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Amplitude& left, const Coef& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const Amplitude& left, const Coef& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator/( const Amplitude& left, const Coef& right )
{
  return Amplitude( left, right, Operation::div );
}



// Operations of coefficient expressions and resonances.
const Amplitude operator+( const CoefExpr& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const CoefExpr& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const CoefExpr& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator+( const Resonance& left, const CoefExpr& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Resonance& left, const CoefExpr& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const Resonance& left, const CoefExpr& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator/( const Resonance& left, const CoefExpr& right )
{
  return Amplitude( left, right, Operation::div );
}



// Operations with coefficient expressions and amplitudes.
const Amplitude operator+( const CoefExpr& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const CoefExpr& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const CoefExpr& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator+( const Amplitude& left, const CoefExpr& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Amplitude& left, const CoefExpr& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const Amplitude& left, const CoefExpr& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator/( const Amplitude& left, const CoefExpr& right )
{
  return Amplitude( left, right, Operation::div );
}



// Operations with resonances and amplitudes.
const Amplitude operator+( const Resonance& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Resonance& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator+( const Amplitude& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Amplitude& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::minus );
}



const Amplitude operator+( const Fvector& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Fvector& left, const Amplitude& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator+( const Amplitude& left, const Fvector& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Amplitude& left, const Fvector& right )
{
  return Amplitude( left, right, Operation::minus );
}

