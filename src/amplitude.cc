
#include <stack>

#include <cfit/amplitude.hh>
#include <cfit/coef.hh>
#include <cfit/resonance.hh>

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
  _expression += "k"; // k = coefficient.
}

void Amplitude::append( const Resonance& reso )
{
  _resos.push_back( reso.copy() );
  _expression += "r"; // r = resonance.
}

void Amplitude::append( const Amplitude& ampl )
{
  _ctnts.insert( _ctnts.end(), ampl._ctnts.begin(), ampl._ctnts.end() );
  _coefs.insert( _coefs.end(), ampl._coefs.begin(), ampl._coefs.end() );
  _opers.insert( _opers.end(), ampl._opers.begin(), ampl._opers.end() );

  std::transform( ampl._resos.begin(), ampl._resos.end(),
                  std::back_inserter( _resos ), std::mem_fun( &Resonance::copy ) );

  _expression += ampl._expression;
}

void Amplitude::append( const Operation::Op& oper )
{
  _opers.push_back( oper );
  _expression += "b"; // b = binary operation.
}


void Amplitude::setPars( const std::map< std::string, Parameter >& pars )
{
  // Set the values of the coefficients.
  typedef std::vector< Coef >::iterator cIter;
  for ( cIter coef = _coefs.begin(); coef != _coefs.end(); ++coef )
    coef->setValue( pars.find( coef->real().name() )->second.value(),
                    pars.find( coef->imag().name() )->second.value() );

  // Propagate the values to the list of resonances.
  typedef std::vector< Resonance* >::iterator rIter;
  for ( rIter reso = _resos.begin(); reso != _resos.end(); ++reso )
    (*reso)->setPars( pars );
}


// Binary operations with complex numbers.
std::complex< double > Amplitude::operate( const std::complex< double >& x,
                                           const std::complex< double >& y,
                                           const Operation::Op&          oper ) const throw( PdfException )
{
  if ( oper == Operation::plus )
    return x + y;
  if ( oper == Operation::minus )
    return x - y;
  if ( oper == Operation::mult )
    return x * y;
  if ( oper == Operation::div )
    return x / y;

  throw PdfException( std::string( "Parse error: unknown binary operation " ) + Operation::tostring( oper ) + "." );
}


// Evaluate the amplitude at the given point, with the current values of its parameters.
std::complex< double > Amplitude::evaluate( const PhaseSpace& ps,
                                            const double&     mSq12,
                                            const double&     mSq13,
                                            const double&     mSq23 ) const throw( PdfException )
{
  std::stack< std::complex< double > > values;

  std::complex< double > x;
  std::complex< double > y;

  std::vector< std::complex< double > >::const_iterator ctt = _ctnts.begin();
  std::vector< Coef                   >::const_iterator coe = _coefs.begin();
  std::vector< Resonance*             >::const_iterator res = _resos.begin();
  std::vector< Operation::Op          >::const_iterator ops = _opers.begin();

  // Parsing loop.
  typedef std::string::const_iterator eIter;
  for ( eIter ch = _expression.begin(); ch != _expression.end(); ch++ )
    if ( *ch == 'c' )
      values.push( *ctt++ );
    else if ( *ch == 'k' )
      values.push( coe++->value() );
    else if ( *ch == 'r' )
      values.push( (*res++)->evaluate( ps, mSq12, mSq13, mSq23 ) );
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

        values.push( operate( x, y, *ops++ ) );
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

const Amplitude& Amplitude::operator=( const Amplitude& right )
{
  append( right );
  return *this;
}

const Amplitude& Amplitude::operator+=( const double& right )
{
  append( right           );
  append( Operation::plus );
  return *this;
}

const Amplitude& Amplitude::operator-=( const double& right )
{
  append( right            );
  append( Operation::minus );
  return *this;
}

const Amplitude& Amplitude::operator*=( const double& right )
{
  append( right           );
  append( Operation::mult );
  return *this;
}

const Amplitude& Amplitude::operator/=( const double& right )
{
  append( right          );
  append( Operation::div );
  return *this;
}

const Amplitude& Amplitude::operator+=( const Coef& right )
{
  append( right           );
  append( Operation::plus );
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


const Amplitude& Amplitude::operator+=( const Resonance& right )
{
  append( right           );
  append( Operation::plus );
  return *this;
}

const Amplitude& Amplitude::operator-=( const Resonance& right )
{
  append( right            );
  append( Operation::minus );
  return *this;
}

const Amplitude& Amplitude::operator+=( const Amplitude& right )
{
  append( right           );
  append( Operation::plus );
  return *this;
}

const Amplitude& Amplitude::operator-=( const Amplitude& right )
{
  append( right            );
  append( Operation::minus );
  return *this;
}




// Operations of objects with themselves.
const Amplitude operator+( const Coef& left, const Coef& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Coef& left, const Coef& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const Coef& left, const Coef& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator/( const Coef& left, const Coef& right )
{
  return Amplitude( left, right, Operation::div );
}

const Amplitude operator+( const Resonance& left, const Resonance& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Resonance& left, const Resonance& right )
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


// Operations of constants and coefficients.
const Amplitude operator+( const double& left, const Coef& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const double& left, const Coef& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const double& left, const Coef& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator/( const double& left, const Coef& right )
{
  return Amplitude( left, right, Operation::div );
}

const Amplitude operator+( const Coef& left, const double& right )
{
  return Amplitude( left, right, Operation::plus );
}

const Amplitude operator-( const Coef& left, const double& right )
{
  return Amplitude( left, right, Operation::minus );
}

const Amplitude operator*( const Coef& left, const double& right )
{
  return Amplitude( left, right, Operation::mult );
}

const Amplitude operator/( const Coef& left, const double& right )
{
  return Amplitude( left, right, Operation::div );
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

